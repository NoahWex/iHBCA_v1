#!/usr/bin/env python3
"""
Per-compartment scVI integration on the iHBCA v1 atlas (2.12M cells).

Trains a scVI model with negative binomial likelihood on a single compartment
(Immune / Epithelial / Stromal). Reads counts from a sparse NPZ matrix +
gene_data.csv + cell metadata CSV — avoids loading a multi-GB h5ad. Saves the
trained model (with adata.h5ad embedded for downstream scANVI), 50D latent
representation, UMAP coordinates, and Leiden clusterings at 10 resolutions.

Canonical configuration (driven by the SLURM wrapper run_scvi_full.sh):
    n_latent=50, n_layers=3, n_hvg=4000, max_epochs=300 (early-stopping
    patience=20), seed=42, batch_key=dataset.

Usage:
    python 01_compartment_scvi_full.py \
        --counts-npz <counts.npz> \
        --gene-data <gene_data.csv> \
        --metadata <cell_metadata.csv> \
        --compartment Immune \
        --n-latent 50 \
        --output-dir <out_dir>
"""

import argparse
import gc
import os
import resource
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp


def log_peak_rss(label):
    """Print peak RSS so far (Linux: KB; macOS: bytes)."""
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"  [RSS] {label}: {peak_kb / 1024 / 1024:.1f} GB peak", flush=True)

# Numba cache patch
try:
    import numba
    _o_njit = numba.njit
    def _njit_nc(*a, **k): k.pop("cache", None); return _o_njit(*a, **k)
    numba.njit = _njit_nc
    _o_vec = numba.vectorize
    def _vec_nc(*a, **k): k.pop("cache", None); return _o_vec(*a, **k)
    numba.vectorize = _vec_nc
except ImportError:
    pass

import scanpy as sc
import scvi

# Canonical scVI parameters
SCVI_PARAMS = {
    "n_latent": 50,
    "n_layers": 3,
    "max_epochs": 300,
    "batch_size": 128,
    "early_stopping": True,
    "early_stopping_patience": 20,
}

RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]


def load_npz_counts(npz_path, gene_data_path):
    """Load NPZ as raw integer counts (no normalization — scVI handles that)."""
    print(f"Loading NPZ: {npz_path}")
    t0 = time.time()
    counts = sp.load_npz(npz_path)
    if not sp.isspmatrix_csr(counts):
        counts = counts.tocsr()

    gene_df = pd.read_csv(gene_data_path, index_col=0)

    n_cells, n_genes = counts.shape
    if n_genes != len(gene_df) and n_cells == len(gene_df):
        counts = counts.T.tocsr()
        n_cells, n_genes = counts.shape

    assert n_genes == len(gene_df), f"Gene mismatch: {n_genes} vs {len(gene_df)}"

    # Convert to float32 for scVI but keep integer values
    if counts.dtype != np.float32:
        counts = counts.astype(np.float32)

    print(f"  {n_cells} cells x {n_genes} genes, {time.time()-t0:.1f}s")
    return counts, gene_df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-npz", required=True)
    parser.add_argument("--gene-data", required=True)
    parser.add_argument("--metadata", required=True,
                        help="cell_metadata_aligned.csv with level0_annotation, dataset, native labels")
    parser.add_argument("--compartment", required=True, choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--n-latent", type=int, default=50)
    parser.add_argument("--n-hvg", type=int, default=4000,
                        help="HVG count (default: 4000 for expanded gene space; "
                             "was 2000 for original 15K-gene inner join)")
    parser.add_argument("--seed", type=int, default=42,
                        help="Random seed for reproducibility (numpy + torch + scvi)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--test", action="store_true", help="Subset to 5000 cells, 10 epochs")
    args = parser.parse_args()

    # --- Seed everything for reproducibility ---
    np.random.seed(args.seed)
    try:
        import torch
        torch.manual_seed(args.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(args.seed)
    except ImportError:
        pass
    try:
        scvi.settings.seed = args.seed
    except AttributeError:
        pass

    SCVI_PARAMS["n_latent"] = args.n_latent
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 70)
    print(f"Full-Size Per-Compartment scVI: {args.compartment}, n_latent={args.n_latent}")
    print("=" * 70)

    # --- Load counts ---
    counts, gene_df = load_npz_counts(args.counts_npz, args.gene_data)

    # --- Load metadata ---
    print("Loading metadata...")
    # Probe available columns so we tolerate metadata variants where the
    # per-study native_* label columns may or may not be present.
    meta_header = pd.read_csv(args.metadata, nrows=0).columns.tolist()
    expected_native = ["native_kumar", "native_nee", "native_gray",
                       "native_murrow", "native_twigger", "native_pal", "native_reed"]
    native_cols = [c for c in expected_native if c in meta_header]
    missing_native = [c for c in expected_native if c not in meta_header]
    if missing_native:
        print(f"  NOTE: native_* columns missing from metadata: {missing_native}")
        print(f"  Proceeding without native labels in compartment h5ad obs.")
    required_cols = ["numeric_id", "dataset", "level0_annotation"]
    missing_req = [c for c in required_cols if c not in meta_header]
    if missing_req:
        raise ValueError(f"Required columns missing from metadata: {missing_req}")
    meta = pd.read_csv(args.metadata,
                       usecols=required_cols + native_cols,
                       low_memory=False)
    print(f"  Metadata: {len(meta)} rows, native_cols={len(native_cols)}/7")

    # --- Subset to compartment ---
    comp_mask = meta["level0_annotation"] == args.compartment
    comp_ids = meta.loc[comp_mask, "numeric_id"].values
    print(f"  {args.compartment}: {len(comp_ids)} cells")

    # Subset counts matrix by row index (numeric_id IS the row index)
    comp_counts = counts[comp_ids]
    comp_meta = meta.loc[comp_mask].reset_index(drop=True)

    # Free full matrix
    del counts
    gc.collect()
    print(f"  Freed full counts matrix")

    # Build AnnData
    adata = sc.AnnData(X=comp_counts, var=pd.DataFrame(index=gene_df.index))
    adata.var["symbol"] = gene_df["symbol"].values

    # Add metadata
    adata.obs["dataset"] = comp_meta["dataset"].values
    for col in native_cols:
        adata.obs[col] = comp_meta[col].values
    adata.obs_names = [str(i) for i in comp_ids]

    del comp_counts, comp_meta
    gc.collect()

    print(f"  AnnData: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  Studies: {adata.obs['dataset'].nunique()}")

    if args.test:
        n_test = min(5000, adata.n_obs)
        np.random.seed(42)
        idx = np.random.choice(adata.n_obs, n_test, replace=False)
        adata = adata[idx].copy()
        SCVI_PARAMS["max_epochs"] = 10
        print(f"  TEST: subsetted to {adata.n_obs}, max_epochs=10")

    # --- Filter genes ---
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  Genes after filtering: {adata.n_vars}")

    log_peak_rss("after AnnData construction")

    # --- HVG selection + scVI ---
    n_studies = adata.obs["dataset"].nunique()
    if n_studies < 2:
        print("  Only 1 study — PCA + Leiden only")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg)
        sc.tl.pca(adata, n_comps=30)
        embed_key = "X_pca"
    else:
        params = SCVI_PARAMS.copy()

        # Save raw counts, then normalize for HVG selection
        raw_X = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=args.n_hvg, batch_key="dataset")
        hvg = adata.var_names[adata.var["highly_variable"]]
        print(f"  HVGs: {len(hvg)} (requested {args.n_hvg})")
        log_peak_rss("after HVG selection")

        # Build raw-counts HVG AnnData for scVI
        hvg_idx = adata.var["highly_variable"].values
        adata_raw = sc.AnnData(
            X=raw_X[:, hvg_idx],
            obs=adata.obs.copy(),
            var=adata.var.loc[hvg].copy(),
        )
        del raw_X

        # Train scVI
        scvi.model.SCVI.setup_anndata(adata_raw, batch_key="dataset")
        model = scvi.model.SCVI(
            adata_raw,
            n_latent=params["n_latent"],
            n_layers=params["n_layers"],
        )
        print(f"  Training scVI (n_latent={params['n_latent']}, "
              f"n_layers={params['n_layers']}, max_epochs={params['max_epochs']})...")
        model.train(
            max_epochs=params["max_epochs"],
            batch_size=params["batch_size"],
            early_stopping=params["early_stopping"],
            early_stopping_patience=params["early_stopping_patience"],
        )

        latent = model.get_latent_representation()
        # The output h5ad is anchored on the 4000-HVG matrix that scVI was
        # actually trained on (adata_raw), not the 27K-gene pre-HVG matrix.
        # This keeps the scVI and scANVI output h5ads in the same gene space
        # so scib.metrics.pcr_comparison is comparable across methods.
        adata_raw.obsm["X_scVI"] = latent
        embed_key = "X_scVI"
        print(f"  scVI done. Latent: {latent.shape}")
        log_peak_rss("after scVI training")

        # Save the trained scVI model so scANVI can reuse it via SCANVI.from_scvi_model()
        scvi_model_dir = os.path.join(args.output_dir, "scvi_model")
        model.save(scvi_model_dir, overwrite=True, save_anndata=True)
        print(f"  scVI model saved: {scvi_model_dir}")

        # Swap: discard the 27K-gene wrapper, promote adata_raw to the
        # canonical `adata` for downstream neighbors/UMAP/leiden/save.
        del adata, model
        adata = adata_raw
        del adata_raw
        gc.collect()

    # --- UMAP + Leiden ---
    print("  Neighbors + UMAP...")
    sc.pp.neighbors(adata, use_rep=embed_key, n_neighbors=30)
    sc.tl.umap(adata)

    print("  Leiden clustering...")
    for res in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
        n_cl = adata.obs[f"leiden_{res}"].nunique()
        print(f"    res {res}: {n_cl} clusters")

    # --- Save ---
    comp_short = {"Immune": "imm", "Epithelial": "epi", "Stromal": "str"}[args.compartment]

    # H5AD
    h5ad_out = os.path.join(args.output_dir, f"{comp_short}_scvi.h5ad")
    obs_clean = adata.obs.copy()
    for col in obs_clean.columns:
        if hasattr(obs_clean[col], "cat"):
            obs_clean[col] = obs_clean[col].astype(str).replace("nan", "")
        elif obs_clean[col].dtype == object or obs_clean[col].isna().any():
            obs_clean[col] = obs_clean[col].fillna("").astype(str)
    out_adata = sc.AnnData(
        X=adata.X,
        obs=obs_clean,
        var=adata.var[["symbol"]].copy(),
        obsm={k: np.array(adata.obsm[k]) for k in adata.obsm.keys()},
    )
    out_adata.write_h5ad(h5ad_out)
    del out_adata
    print(f"  H5AD saved: {h5ad_out}")
    print(f"    obsm: {list(adata.obsm.keys())}")

    # Leiden CSV
    leiden_cols = [f"leiden_{r}" for r in RESOLUTIONS]
    leiden_df = adata.obs[leiden_cols].copy()
    leiden_df.insert(0, "cell_id", adata.obs_names)
    leiden_df.to_csv(os.path.join(args.output_dir, f"{comp_short}_leiden.csv"), index=False)

    # UMAP CSV
    umap_df = pd.DataFrame({
        "cell_id": adata.obs_names,
        "UMAP_1": adata.obsm["X_umap"][:, 0],
        "UMAP_2": adata.obsm["X_umap"][:, 1],
    })
    umap_df.to_csv(os.path.join(args.output_dir, f"{comp_short}_scvi_umap.csv"), index=False)

    print(f"\n{'='*70}")
    print("Done.")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
