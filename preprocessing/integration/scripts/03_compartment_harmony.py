#!/usr/bin/env python3
"""
Per-compartment Harmony integration — comparative scib benchmark only.

Reads the same per-compartment NPZ + metadata inputs as scVI/scANVI and runs:
  1. Normalize + log1p
  2. HVG selection
  3. PCA (n_pcs=50)
  4. Harmony correction on PCA (batch_key=dataset, max_iter_harmony=20)
  5. Neighbors + UMAP + Leiden on the Harmony embedding

Included as a comparative integration baseline in the scib benchmark suite
(see 05_scib_benchmark.py and the 05*/06* sidecars). Not part of the
canonical pipeline — scANVI at n_latent=50 is the published winner.

Output mirrors scVI: H5AD with X_pca + X_pca_harmony + X_umap + leiden_*.

Usage:
    python 03_compartment_harmony.py \
        --counts-npz <{comp}_counts.npz> \
        --gene-data <gene_data.csv> \
        --metadata <{comp}_metadata.csv> \
        --compartment Immune \
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

# Numba cache patch (same as scVI)
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
import scanpy.external as sce


RESOLUTIONS = [0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0]


def log_peak_rss(label):
    peak_kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(f"  [RSS] {label}: {peak_kb / 1024 / 1024:.1f} GB peak", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-npz", required=True)
    parser.add_argument("--gene-data", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--compartment", required=True,
                        choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--n-hvg", type=int, default=4000)
    parser.add_argument("--n-pcs", type=int, default=50,
                        help="Number of PCs for PCA + Harmony (default: 50)")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)
    os.makedirs(args.output_dir, exist_ok=True)
    np.random.seed(args.seed)
    sc.settings.seed = args.seed

    print("=" * 70, flush=True)
    print(f"Per-Compartment Harmony: {args.compartment}, n_pcs={args.n_pcs}",
          flush=True)
    print("=" * 70, flush=True)

    # --- Load counts (same pattern as scVI) ---
    print(f"Loading NPZ: {args.counts_npz}", flush=True)
    t0 = time.time()
    counts = sp.load_npz(args.counts_npz)
    if not sp.isspmatrix_csr(counts):
        counts = counts.tocsr()
    gene_df = pd.read_csv(args.gene_data, index_col=0)
    n_cells, n_genes = counts.shape
    if n_genes != len(gene_df) and n_cells == len(gene_df):
        counts = counts.T.tocsr()
        n_cells, n_genes = counts.shape
    if counts.dtype != np.float32:
        counts = counts.astype(np.float32)
    print(f"  {n_cells} cells x {n_genes} genes, {time.time() - t0:.1f}s", flush=True)

    # --- Load metadata ---
    print("Loading metadata...", flush=True)
    meta = pd.read_csv(
        args.metadata,
        usecols=["numeric_id", "dataset", "level0_annotation"],
        low_memory=False,
    )
    print(f"  Metadata: {len(meta)} rows", flush=True)

    # --- Subset to compartment ---
    comp_mask = meta["level0_annotation"] == args.compartment
    comp_ids = meta.loc[comp_mask, "numeric_id"].values
    print(f"  {args.compartment}: {len(comp_ids)} cells", flush=True)
    comp_counts = counts[comp_ids]
    comp_meta = meta.loc[comp_mask].reset_index(drop=True)
    del counts
    gc.collect()

    # --- Build AnnData ---
    adata = sc.AnnData(X=comp_counts, var=pd.DataFrame(index=gene_df.index))
    adata.var["symbol"] = gene_df["symbol"].values
    adata.obs["dataset"] = comp_meta["dataset"].values
    adata.obs_names = [str(i) for i in comp_ids]
    del comp_counts, comp_meta
    gc.collect()
    print(f"  AnnData: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)
    print(f"  Studies: {adata.obs['dataset'].nunique()}", flush=True)
    log_peak_rss("after AnnData construction")

    # --- Filter + normalize + HVG (same recipe as scVI for fair comparison) ---
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"  Genes after filtering: {adata.n_vars}", flush=True)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=args.n_hvg, batch_key="dataset"
    )
    hvg = adata.var_names[adata.var["highly_variable"]]
    print(f"  HVGs: {len(hvg)} (requested {args.n_hvg})", flush=True)
    log_peak_rss("after HVG selection")

    # Subset to HVGs for PCA (matches scVI which trains on HVG)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_hvg, max_value=10)
    log_peak_rss("after scale")

    # --- PCA ---
    print(f"\nPCA (n_pcs={args.n_pcs})...", flush=True)
    t0 = time.time()
    sc.tl.pca(adata_hvg, n_comps=args.n_pcs, random_state=args.seed)
    print(f"  PCA in {time.time() - t0:.1f}s", flush=True)
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
    log_peak_rss("after PCA")

    # --- Harmony ---
    print(f"\nHarmony correction...", flush=True)
    t0 = time.time()
    sce.pp.harmony_integrate(
        adata,
        key="dataset",
        basis="X_pca",
        adjusted_basis="X_pca_harmony",
        max_iter_harmony=20,
    )
    print(f"  Harmony in {(time.time() - t0) / 60:.1f} min", flush=True)
    log_peak_rss("after Harmony")

    del adata_hvg
    gc.collect()

    # --- Neighbors + UMAP + Leiden on harmony embedding ---
    print("\nNeighbors + UMAP on X_pca_harmony...", flush=True)
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=30)
    sc.tl.umap(adata)

    print("Leiden clustering...", flush=True)
    for res in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
        n_cl = adata.obs[f"leiden_{res}"].nunique()
        print(f"    res {res}: {n_cl} clusters", flush=True)

    # --- Save outputs (same naming convention as scVI) ---
    comp_short = {"Immune": "imm", "Epithelial": "epi", "Stromal": "str"}[
        args.compartment
    ]
    h5ad_out = os.path.join(args.output_dir, f"{comp_short}_harmony.h5ad")
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
    print(f"  H5AD saved: {h5ad_out}", flush=True)
    print(f"    obsm: {list(adata.obsm.keys())}", flush=True)

    # CSVs
    leiden_cols = [f"leiden_{r}" for r in RESOLUTIONS]
    leiden_df = adata.obs[leiden_cols].copy()
    leiden_df.insert(0, "cell_id", adata.obs_names)
    leiden_df.to_csv(
        os.path.join(args.output_dir, f"{comp_short}_leiden.csv"), index=False
    )
    pd.DataFrame({
        "cell_id": adata.obs_names,
        "UMAP_1": adata.obsm["X_umap"][:, 0],
        "UMAP_2": adata.obsm["X_umap"][:, 1],
    }).to_csv(
        os.path.join(args.output_dir, f"{comp_short}_harmony_umap.csv"), index=False
    )

    print("\nDone.", flush=True)


if __name__ == "__main__":
    main()
