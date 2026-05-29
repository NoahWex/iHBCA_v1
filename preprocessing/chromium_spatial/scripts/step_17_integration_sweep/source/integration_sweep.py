#!/usr/bin/env python3
"""
step_17_integration_sweep / integration_sweep.py

Train one integration config (scVI / scVI_austin / Harmony / PCA / CONCORD) on the
FIVE-WAY-filtered atlas-quality cell set produced by step 16 (or a per-
compartment subset). Writes an integrated.h5ad under the canonical sweep
output tree.

The sweep enumerates 9 configs x {full, epi, str, imm} targets via
wrapper.py + the SLURM array in run/run_step_17_sweep.sh. This script
receives ONE (config, target) pair and produces one integrated.h5ad.

Canonical inputs (resolved by wrapper.py via paths.yaml):
  --vf-dir          post-FIVE-WAY per-sample VFs from step 16
                    (CFG_PREPROCESSING_STEP16_VF_DIR)
  --cell-list       either central_cell_status.csv (full target) or one of
                    the compartment cell lists under step 14 outputs
  --sample-list     raw data manifest (CFG_RAW_DATA_SAMPLE_MANIFEST)
  --labels-csv      optional, used downstream by scIB scoring; not read here
                    beyond being forwarded

Outputs:
  {output_dir}/integrated.h5ad
  {output_dir}/integration_metadata.csv
  {output_dir}/latent_embedding.csv   (if --output-embedding)

"""

# ----------------------------------------------------------------------------
# Numba cache patch — disables cache kwargs for numba.njit/vectorize because
# the scGPT container has a read-only /home and numba tries to cache under
# site-packages. Must run BEFORE any import of numba-using modules.
# ----------------------------------------------------------------------------
try:
    import numba
    _o = numba.njit
    def _nc(*a, **k):  # noqa: E306
        k.pop("cache", None)
        return _o(*a, **k)
    numba.njit = _nc
    _ov = numba.vectorize
    def _ncv(*a, **k):  # noqa: E306
        k.pop("cache", None)
        return _ov(*a, **k)
    numba.vectorize = _ncv
except ImportError:
    pass

import argparse
import gc
import glob
import os
import sys
import time

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc


# ============================================================================
# Validation preamble helpers
# ============================================================================

def _die(msg, hint=None):
    """Fail loud with an optional diagnosis hint."""
    sys.stderr.write(f"FATAL: {msg}\n")
    if hint:
        sys.stderr.write(f"  hint: {hint}\n")
    sys.exit(2)


def _assert_file(path, label):
    if not path:
        _die(f"{label} path is empty")
    if not os.path.exists(path):
        _die(f"{label} not found: {path}",
             hint="verify CFG_* paths and that upstream step outputs were written")
    if os.path.getsize(path) == 0:
        _die(f"{label} is empty (0 bytes): {path}")


def _assert_dir(path, label):
    if not path:
        _die(f"{label} path is empty")
    if not os.path.isdir(path):
        _die(f"{label} not a directory: {path}")


def _load_cell_list(cell_list_path):
    """Read a passing-cell CSV into a set of cell IDs.

    Accepts any of: `cell_id`, `obs_name`, or first column.
    """
    cell_df = pd.read_csv(cell_list_path)
    if "cell_id" in cell_df.columns:
        col = "cell_id"
    elif "obs_name" in cell_df.columns:
        col = "obs_name"
    else:
        col = cell_df.columns[0]
    ids = set(cell_df[col].astype(str))
    if not ids:
        _die(f"cell list {cell_list_path} produced 0 IDs "
             f"(col={col}, rows={len(cell_df)})")
    return ids, col


def _collect_vfs(vf_dir, vf_subdir):
    """Glob per-sample VF files under vf_dir/*/vf_subdir/variable_features.txt.

    Also accepts a flat directory of *_vfs_filtered.txt files (step 16 layout)
    when vf_subdir is empty or 'filtered_vfs'. Returns union VF set.

    Fail-loud addition: distinguishes "dir exists but no VF files found" from
    "dir missing" and prints both patterns tried so the reviewer can debug.
    """
    _assert_dir(vf_dir, "--vf-dir")
    patterns_tried = []

    # Pattern A: classic pipeline_comparison_20260401 layout
    # {vf_dir}/{sample}/{vf_subdir}/variable_features.txt
    pattern_a = os.path.join(vf_dir, f"*/{vf_subdir}/variable_features.txt")
    patterns_tried.append(pattern_a)
    vf_files = glob.glob(pattern_a)

    # Pattern B: canonical step 16 layout
    # {vf_dir}/{sample}_vfs_filtered.txt
    if not vf_files:
        pattern_b = os.path.join(vf_dir, "*_vfs_filtered.txt")
        patterns_tried.append(pattern_b)
        vf_files = glob.glob(pattern_b)

    # Pattern C: caller already passed a glob
    if not vf_files:
        patterns_tried.append(vf_dir)
        vf_files = glob.glob(vf_dir)

    if not vf_files:
        _die(f"No VF files found under {vf_dir}",
             hint=f"tried patterns: {patterns_tried}. For canonical step 17, "
                  f"--vf-dir should point at CFG_PREPROCESSING_STEP16_VF_DIR "
                  f"and --vf-subdir should be 'filtered_vfs' or empty.")

    all_vfs = set()
    for vf_file in vf_files:
        with open(vf_file) as fh:
            vfs = [ln.strip() for ln in fh if ln.strip()]
        all_vfs.update(vfs)

    if not all_vfs:
        _die(f"VF files found but all empty: {vf_files[:3]}...")

    print(f"Union VFs: {len(all_vfs)} from {len(vf_files)} files "
          f"(subdir hint: {vf_subdir!r})")
    return all_vfs


# ============================================================================
# Data loading
# ============================================================================

def load_per_sample_data(sample_list_path, cell_list_path, vf_dir,
                         vf_subdir="filtered_vfs"):
    """Load per-sample raw counts, subset to passing cells, merge.

    Modernizations:
      - cell list loader extracted to _load_cell_list (used in two places)
      - VF collection extracted to _collect_vfs with multi-pattern fallback
      - fail-loud on 0-cell intersections with an explicit hint
    """
    _assert_file(sample_list_path, "--sample-list")
    _assert_file(cell_list_path, "--cell-list")

    samples = pd.read_csv(sample_list_path, sep="\t", header=None,
                          names=["sample_id", "h5_path"])
    print(f"Samples: {len(samples)}")

    passing, passing_col = _load_cell_list(cell_list_path)
    print(f"Passing cells: {len(passing)} (col={passing_col})")

    all_vfs = _collect_vfs(vf_dir, vf_subdir)

    # Load per-sample, subset, merge
    adatas = []
    samples_missing = []
    samples_zero_passing = []
    for _, row in samples.iterrows():
        h5_path = row["h5_path"]
        sample_id = row["sample_id"]

        if not os.path.exists(h5_path):
            samples_missing.append(sample_id)
            continue

        adata = sc.read_10x_h5(h5_path)
        adata.var_names_make_unique()
        adata.obs_names = [f"{sample_id}_{bc}" for bc in adata.obs_names]
        adata.obs["sample_id"] = sample_id
        adata.obs["patient_id"] = sample_id.split("_")[0]

        keep = [n for n in adata.obs_names if n in passing]
        if not keep:
            samples_zero_passing.append(sample_id)
            continue

        adata = adata[keep].copy()
        print(f"  {sample_id}: {len(keep)} cells")
        adatas.append(adata)

    if not adatas:
        _die("No samples contributed cells to the merged object",
             hint=f"missing H5: {samples_missing[:5]}; "
                  f"zero-passing: {samples_zero_passing[:5]}; "
                  f"cell list size: {len(passing)}. "
                  f"Check that sample IDs in cell list match barcode prefixes.")

    print(f"Merging {len(adatas)} samples "
          f"({len(samples_missing)} missing, {len(samples_zero_passing)} zero-passing)...")
    merged = ad.concat(adatas, join="outer", fill_value=0)
    merged.var_names_make_unique()
    print(f"  Merged: {merged.n_obs} cells x {merged.n_vars} genes")

    shared_vfs = [g for g in all_vfs if g in merged.var_names]
    if not shared_vfs:
        _die("0 VFs in merged gene set",
             hint="VF file gene symbols do not intersect the merged H5 var_names. "
                  "Check that VF files come from the same gene space "
                  "(same 10x reference) as the raw H5 files.")
    print(f"  VFs in merged: {len(shared_vfs)}")

    del adatas
    gc.collect()
    return merged, shared_vfs


def load_merged_input(merged_input, cell_list_path, vf_dir, vf_subdir):
    """Load pre-assembled H5AD and apply cell-list filter + VF intersection.

    """
    _assert_file(merged_input, "--merged-input")
    _assert_file(cell_list_path, "--cell-list")

    print(f"Reading pre-assembled counts: {merged_input}")
    adata = sc.read_h5ad(merged_input)
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    keep_cells, col = _load_cell_list(cell_list_path)
    n_before = adata.n_obs
    mask = adata.obs_names.isin(keep_cells)
    adata = adata[mask].copy()
    print(f"  Cell-list filter: {n_before} -> {adata.n_obs} "
          f"(list={len(keep_cells)} col={col})")
    if adata.n_obs == 0:
        _die("cell-list filter left 0 cells",
             hint="cell list IDs do not match merged H5AD obs_names. "
                  "Inspect head(adata.obs_names) vs head(cell_list) formats.")

    all_vfs = _collect_vfs(vf_dir, vf_subdir)
    vf_genes = [g for g in all_vfs if g in adata.var_names]
    if not vf_genes:
        _die("0 VFs intersect merged input var_names")
    print(f"  VFs in merged: {len(vf_genes)}")
    return adata, vf_genes


# ============================================================================
# Integration methods
# ============================================================================

def run_scvi(adata, vf_genes, n_latent, batch_key,
             max_epochs=300, early_stopping_patience=15, batch_size=128,
             n_layers=3, dropout_rate=0.1):
    """Run scVI integration, writing latent representation to obsm['X_emb'].


    NOTE on scVI GPU stochasticity:
      scVI training is GPU-stochastic — two runs with the same hyperparameters
      and seed will produce slightly different latents because cuDNN kernels
      are non-deterministic by default. The winner scVI model checkpoint is
      therefore treated as CANONICAL per user directive (see L3_spec.md:14-16):
      it is captured once and NOT regenerated, because bit-exact reproduction
      from (counts + config) is impossible. Once a winner is selected by
      scib_aggregate.py + svd_diagnostic.py, freeze the model directory
      (model/model.pt) and point all downstream steps at it.
    """
    adata_vf = adata[:, vf_genes].copy()

    import scvi as scvi_module
    scvi_module.model.SCVI.setup_anndata(adata_vf, batch_key=batch_key)
    model = scvi_module.model.SCVI(
        adata_vf,
        n_latent=n_latent,
        n_layers=n_layers,
        dropout_rate=dropout_rate,
    )
    print(f"  Training scVI (n_latent={n_latent}, max_epochs={max_epochs}, "
          f"n_layers={n_layers}, dropout={dropout_rate}, batch_size={batch_size})")
    model.train(
        max_epochs=max_epochs,
        batch_size=batch_size,
        early_stopping=True,
        early_stopping_patience=early_stopping_patience,
    )
    latent = model.get_latent_representation()
    adata.obsm["X_emb"] = latent
    print(f"  scVI latent: {latent.shape}")

    del model, adata_vf
    gc.collect()
    return adata


def run_harmony(adata, vf_genes, batch_key, n_comps=50, theta=2.0):
    """Harmony-on-PCA via harmonypy (bypasses sc.external wrapper).

    Kept verbatim including the Z_corr shape fix for harmonypy version drift.
    Extended to accept n_comps and theta for sweep grid coverage.

    Args:
        n_comps: PCA components fed into Harmony (analogous to n_latent in scVI).
        theta: diversity clustering penalty (harmonypy default 2.0). Higher values
               push Harmony to mix batches more aggressively; lower values preserve
               more local batch structure. Passed as a scalar — harmonypy broadcasts
               it across all batch variables in vars_use.
    """
    adata_work = adata[:, vf_genes].copy()
    sc.pp.normalize_total(adata_work, target_sum=1e4)
    sc.pp.log1p(adata_work)
    sc.pp.scale(adata_work, max_value=10)
    sc.tl.pca(adata_work, n_comps=n_comps)

    pca = adata_work.obsm["X_pca"]
    import harmonypy as hm
    ho = hm.run_harmony(pca, adata_work.obs, batch_key, theta=theta)
    z = np.asarray(ho.Z_corr)
    # harmonypy <0.1 returns (n_pcs, n_cells); >=0.2 returns (n_cells, n_pcs)
    if z.shape[0] == pca.shape[1] and z.shape[1] == pca.shape[0]:
        z = z.T

    adata.obsm["X_emb"] = z
    print(f"  Harmony embedding: {adata.obsm['X_emb'].shape}  (n_comps={n_comps}, theta={theta})")
    del adata_work
    gc.collect()
    return adata


def run_concord(adata, vf_genes, batch_key, latent_dim=50):
    """CONCORD contrastive learning integration, writing latent to obsm['X_emb'].

    Reference: Nature Biotechnology 2025 (s41587-025-02950-z), Gartner-Lab/Concord.
    Package: concord-sc (pip install concord-sc), installed in scgpt PYTHONUSERBASE.

    CONCORD trains a single-hidden-layer self-supervised contrastive network.
    domain_key corrects batch effects via dataset-aware sampling.

    Neighbor graph notes:
      CONCORD documentation recommends cosine distance for downstream UMAP.
      For scIB benchmark comparability, Euclidean (default) is used here so all
      methods are scored on the same metric. The cosine recommendation applies to
      visualization, not to the scIB neighbor-graph-based metrics.

    GPU note: auto-detects CUDA via torch.cuda.is_available(). Falls back to CPU
    if no GPU is available (slower, ~5-10x).
    """
    import torch
    import concord as ccd

    adata_vf = adata[:, vf_genes].copy()
    sc.pp.normalize_total(adata_vf, target_sum=1e4)
    sc.pp.log1p(adata_vf)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"  CONCORD device: {device}")

    con = ccd.Concord(
        adata_vf,
        save_dir=None,       # suppress model/log files during sweep
        domain_key=batch_key,
        latent_dim=latent_dim,
        normalize_total=False,   # already normalized above
        log1p=False,             # already log-transformed above
        device=device,
    )
    con.fit_transform(output_key="X_concord", save_model=False)
    adata.obsm["X_emb"] = np.array(adata_vf.obsm["X_concord"])
    print(f"  CONCORD latent: {adata.obsm['X_emb'].shape}")

    del adata_vf, con
    gc.collect()
    return adata


def run_pca_baseline(adata, vf_genes):
    """Unintegrated PCA baseline. Required as an SVD anchor config.

    """
    adata_work = adata[:, vf_genes].copy()
    sc.pp.normalize_total(adata_work, target_sum=1e4)
    sc.pp.log1p(adata_work)
    sc.pp.scale(adata_work, max_value=10)
    sc.tl.pca(adata_work, n_comps=50)
    adata.obsm["X_emb"] = adata_work.obsm["X_pca"]
    adata.obsm["X_pca"] = adata_work.obsm["X_pca"]
    print(f"  PCA baseline embedding: {adata.obsm['X_emb'].shape}")
    del adata_work
    gc.collect()
    return adata


# ============================================================================
# Main
# ============================================================================

def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--method", required=True, choices=["scvi", "harmony", "pca", "concord"])
    parser.add_argument("--config-label", required=True,
                        help="Config identifier (e.g. scvi_n50, scvi_austin_n100). "
                             "Used as the output subdirectory name.")
    parser.add_argument("--n-latent", type=int, default=50)
    parser.add_argument("--vf-dir", required=True,
                        help="Directory containing per-sample VF files. "
                             "For canonical step 17, point at "
                             "CFG_PREPROCESSING_STEP16_VF_DIR.")
    parser.add_argument("--vf-subdir", default="filtered_vfs",
                        help="VF subdirectory name (default 'filtered_vfs'). "
                             "For canonical step 17 the filtered_vfs/ directory "
                             "from step 16 is already flat, so the empty string "
                             "or 'filtered_vfs' both work — _collect_vfs falls "
                             "back through multiple patterns.")
    parser.add_argument("--cell-list", required=True,
                        help="CSV with passing cell IDs (cell_id column). "
                             "For full target: central_cell_status.csv. "
                             "For compartment: compartment cell list.")
    parser.add_argument("--sample-list", required=True,
                        help="TSV: sample_id<tab>h5_path. "
                             "CFG_RAW_DATA_SAMPLE_MANIFEST.")
    parser.add_argument("--batch-key", default="patient_id")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory (per-config). wrapper.py builds "
                             "this under CFG_CANONICAL_SWEEP_{FULL_OBJECT,COMPARTMENTS}/.")
    parser.add_argument("--max-epochs", type=int, default=300)
    parser.add_argument("--n-layers", type=int, default=3)
    parser.add_argument("--dropout-rate", type=float, default=0.1)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--early-stopping-patience", type=int, default=15)
    parser.add_argument("--theta", type=float, default=2.0,
                        help="Harmony diversity clustering penalty (default 2.0). "
                             "Ignored for non-harmony methods.")
    parser.add_argument("--leiden-scoring-resolution", type=float, default=1.0,
                        help="Single Leiden resolution computed during the "
                             "sweep stage for scIB NMI/ARI scoring. The full "
                             "multi-resolution set is computed only on the "
                             "winner by postprocess_winner.py.")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--test", action="store_true",
                        help="Subsample to 5000 cells, 10 epochs, for smoke testing.")
    parser.add_argument("--merged-input", default=None,
                        help="Pre-assembled H5AD. Skips per-sample loading.")
    parser.add_argument("--output-embedding", action="store_true",
                        help="Write latent_embedding.csv alongside integrated.h5ad")
    parser.add_argument("--dry-run", action="store_true",
                        help="Validate inputs, print plan, exit 0 before training")
    return parser.parse_args()


def _dry_run_report(args):
    print("=" * 70)
    print("=== DRY RUN MODE ===")
    print(f"Step:    step_17_integration_sweep")
    print(f"Config:  {args.config_label}")
    print(f"Method:  {args.method}  n_latent={args.n_latent}")
    print(f"Output:  {args.output_dir}")
    print("Inputs:")
    print(f"  --sample-list  : {args.sample_list}")
    print(f"  --cell-list    : {args.cell_list}")
    print(f"  --vf-dir       : {args.vf_dir}  (subdir={args.vf_subdir})")
    if args.merged_input:
        print(f"  --merged-input : {args.merged_input}")
    # Validate files
    _assert_file(args.sample_list, "--sample-list")
    _assert_file(args.cell_list, "--cell-list")
    _assert_dir(args.vf_dir, "--vf-dir")
    if args.merged_input:
        _assert_file(args.merged_input, "--merged-input")
    print("VALIDATION PASSED")


def main():
    args = _parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print(f"step_17 integration sweep: {args.config_label}")
    print(f"  method={args.method}  n_latent={args.n_latent}")
    print(f"  batch_key={args.batch_key}")
    print(f"  output_dir={args.output_dir}")
    print("=" * 70)

    if args.dry_run:
        _dry_run_report(args)
        return

    # ------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------
    if args.merged_input:
        adata, vf_genes = load_merged_input(
            args.merged_input, args.cell_list, args.vf_dir, args.vf_subdir)
    else:
        adata, vf_genes = load_per_sample_data(
            args.sample_list, args.cell_list, args.vf_dir,
            vf_subdir=args.vf_subdir)

    if args.test:
        n_test = min(5000, adata.n_obs)
        np.random.seed(42)
        idx = np.random.choice(adata.n_obs, n_test, replace=False)
        adata = adata[idx].copy()
        args.max_epochs = 10
        print(f"  TEST: {adata.n_obs} cells, {args.max_epochs} epochs")

    # ------------------------------------------------------------
    # Run integration
    # ------------------------------------------------------------
    t0 = time.time()
    if args.method == "scvi":
        adata = run_scvi(
            adata, vf_genes, args.n_latent, args.batch_key,
            max_epochs=args.max_epochs,
            early_stopping_patience=args.early_stopping_patience,
            batch_size=args.batch_size,
            n_layers=args.n_layers,
            dropout_rate=args.dropout_rate,
        )
    elif args.method == "harmony":
        adata = run_harmony(adata, vf_genes, args.batch_key,
                            n_comps=args.n_latent, theta=args.theta)
    elif args.method == "concord":
        adata = run_concord(adata, vf_genes, args.batch_key, latent_dim=args.n_latent)
    elif args.method == "pca":
        adata = run_pca_baseline(adata, vf_genes)
    elapsed = time.time() - t0
    print(f"  Integration time: {elapsed/60:.1f} min")

    # ------------------------------------------------------------
    # Neighbors + UMAP + sweep-stage Leiden (single resolution for scoring)
    # The full multi-resolution Leiden set is computed only on the winner
    # by postprocess_winner.py per the locked Leiden decision.
    # ------------------------------------------------------------
    print("  Neighbors + UMAP + sweep-stage Leiden...")
    sc.pp.neighbors(adata, use_rep="X_emb", n_neighbors=args.n_neighbors)
    sc.tl.umap(adata)

    res = args.leiden_scoring_resolution
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
    n_clust = adata.obs[f"leiden_{res}"].nunique()
    print(f"    leiden res {res}: {n_clust} clusters")

    # ------------------------------------------------------------
    # Write integrated H5AD
    # ------------------------------------------------------------
    h5ad_out = os.path.join(args.output_dir, "integrated.h5ad")
    obs_clean = adata.obs.copy()
    for col in obs_clean.columns:
        if hasattr(obs_clean[col], "cat"):
            obs_clean[col] = obs_clean[col].astype(str).replace("nan", "")
        elif obs_clean[col].dtype == object or obs_clean[col].isna().any():
            obs_clean[col] = obs_clean[col].fillna("").astype(str)
    out = ad.AnnData(
        X=adata.X,
        obs=obs_clean,
        var=adata.var[[]].copy(),
        obsm={k: np.array(adata.obsm[k]) for k in adata.obsm.keys()},
    )
    out.write_h5ad(h5ad_out)
    del out
    print(f"  H5AD: {h5ad_out}")

    # ------------------------------------------------------------
    # Latent embedding CSV (optional)
    # ------------------------------------------------------------
    if args.output_embedding and "X_emb" in adata.obsm:
        emb = adata.obsm["X_emb"]
        emb_cols = [f"{args.method}_{i+1}" for i in range(emb.shape[1])]
        emb_df = pd.DataFrame(emb, index=adata.obs_names, columns=emb_cols)
        emb_df.index.name = "cell_id"
        emb_path = os.path.join(args.output_dir, "latent_embedding.csv")
        emb_df.to_csv(emb_path)
        print(f"  Embedding: {emb_path} ({emb.shape[0]} x {emb.shape[1]})")

    # ------------------------------------------------------------
    # Integration metadata sidecar
    # ------------------------------------------------------------
    meta_out = os.path.join(args.output_dir, "integration_metadata.csv")
    pd.DataFrame({
        "config_label": args.config_label,
        "method": args.method,
        "n_latent": args.n_latent,
        "n_cells": adata.n_obs,
        "n_vfs": len(vf_genes),
        "n_batches": adata.obs[args.batch_key].nunique(),
        "batch_key": args.batch_key,
        "elapsed_min": round(elapsed / 60, 1),
    }, index=[0]).to_csv(meta_out, index=False)
    print(f"  Metadata: {meta_out}")

    print(f"\n{'='*70}")
    print(f"Done: {args.config_label}")


if __name__ == "__main__":
    main()
