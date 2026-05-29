#!/usr/bin/env python3
"""
PCR sidecar — recompute scib PCR on a consistent 4000-HVG baseline.

PCR (principal component regression) compares variance-explained-by-batch
between the input matrix and the integrated embedding. To produce comparable
values across methods, all configurations must use the same input gene space.
This script anchors PCR_before on the 4000-HVG training matrix that scVI
actually saw — saved by SCVI.save(..., save_anndata=True) inside the model
directory — and attaches the relevant embedding (X_scVI from the scVI output
or X_scANVI from the scANVI output) before calling pcr_comparison.

No retraining. No scVI / scANVI re-inference. In-memory reassembly + one PCR
call + one sidecar CSV write.

Usage:
    python 05d_scib_pcr_recompute.py \\
        --method scvi --compartment Immune --n-latent 50 \\
        --scvi-training-adata <scvi_model/adata.h5ad> \\
        --embedding-h5ad <{comp}_{method}.h5ad> \\
        --embedding-key <X_scVI|X_scANVI> \\
        --batch-key dataset \\
        --output-csv <out_pcr.csv>
"""

import argparse
import os
import sys
import time

# Numba cache patch — required before scanpy/scib import.
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

import anndata as ad  # type: ignore[import-not-found]  # container-only
import numpy as np
import pandas as pd
import scanpy as sc  # type: ignore[import-not-found]  # container-only
import scib  # type: ignore[import-not-found]  # container-only


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--method", required=True, choices=["scvi", "scanvi"])
    parser.add_argument("--compartment", required=True,
                        choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--n-latent", type=int, required=True)
    parser.add_argument("--scvi-training-adata", required=True,
                        help="Path to scvi_model/adata.h5ad (the 4000-HVG "
                             "training adata, source of truth for the gene "
                             "space scVI was trained on)")
    parser.add_argument("--embedding-h5ad", required=True,
                        help="Path to the output h5ad containing the "
                             "integration embedding (scvi output or scanvi "
                             "output)")
    parser.add_argument("--embedding-key", required=True,
                        help="obsm key (X_scVI or X_scANVI)")
    parser.add_argument("--batch-key", default="dataset")
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    print("=" * 70, flush=True)
    print(f"PCR recompute: {args.method}_n{args.n_latent} ({args.compartment})",
          flush=True)
    print("=" * 70, flush=True)
    t_total = time.time()

    # --- Load the 4000-HVG training adata (source of truth for gene space) ---
    print(f"Loading scvi training adata: {args.scvi_training_adata}", flush=True)
    t0 = time.time()
    adata = ad.read_h5ad(args.scvi_training_adata)
    print(f"  Shape: {adata.n_obs:,} x {adata.n_vars:,} in "
          f"{time.time() - t0:.1f}s", flush=True)
    print(f"  X dtype: {adata.X.dtype}, nnz: {adata.X.nnz:,}", flush=True)
    if args.batch_key not in adata.obs.columns:
        raise ValueError(
            f"Batch key '{args.batch_key}' not in training adata. "
            f"Available: {list(adata.obs.columns)}"
        )
    print(f"  Batch: {adata.obs[args.batch_key].nunique()} unique", flush=True)

    # --- Normalize + log1p the training adata (scib PCR expects normalized X) ---
    # scvi training adata is saved with RAW counts (SCVI.setup_anndata convention).
    # PCR_before runs sc.tl.pca on adata.X, which is meaningful only on
    # normalized data. Match the preprocessing scVI used during HVG selection.
    print("Normalizing + log1p (for PCR_before baseline)...", flush=True)
    t0 = time.time()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    print(f"  Done in {time.time() - t0:.1f}s", flush=True)

    # --- Load the embedding from the output h5ad and attach to training adata ---
    print(f"Loading embedding from: {args.embedding_h5ad}", flush=True)
    t0 = time.time()
    emb_adata = ad.read_h5ad(args.embedding_h5ad, backed='r')
    if args.embedding_key not in emb_adata.obsm:
        raise ValueError(
            f"Embedding key '{args.embedding_key}' not in obsm. "
            f"Available: {list(emb_adata.obsm.keys())}"
        )
    if emb_adata.n_obs != adata.n_obs:
        raise ValueError(
            f"Cell count mismatch: training adata has {adata.n_obs:,}, "
            f"embedding adata has {emb_adata.n_obs:,}"
        )
    # Align cells by obs_names to be safe
    if list(emb_adata.obs_names) != list(adata.obs_names):
        print("  obs_names differ — realigning by obs_names", flush=True)
        order = [list(emb_adata.obs_names).index(n) for n in adata.obs_names]
        latent = emb_adata.obsm[args.embedding_key][:][order]
    else:
        latent = emb_adata.obsm[args.embedding_key][:]
    emb_adata.file.close()
    adata.obsm[args.embedding_key] = np.asarray(latent)
    print(f"  Attached {args.embedding_key}: shape {latent.shape} in "
          f"{time.time() - t0:.1f}s", flush=True)

    # --- Run pcr_comparison ---
    # scib.metrics.pcr_comparison:
    #   pcr_before = scib.metrics.pcr(adata_pre, ..., recompute_pca=True) on adata.X
    #   pcr_after  = scib.metrics.pcr(adata_post, ..., embed=embed) on obsm[embed]
    #   scale=True: (pcr_before - pcr_after) / pcr_before, clipped to 0 if negative
    print("Running scib.metrics.pcr_comparison...", flush=True)
    t0 = time.time()
    pcr_val = scib.metrics.pcr_comparison(
        adata_pre=adata,
        adata_post=adata,
        covariate=args.batch_key,
        embed=args.embedding_key,
        n_comps=50,
        scale=True,
        verbose=False,
    )
    elapsed_min = (time.time() - t0) / 60
    print(f"  PCR (scaled) = {pcr_val:.6f} in {elapsed_min:.1f} min",
          flush=True)

    # --- Also compute the unscaled raw values for diagnosis ---
    print("Running unscaled pcr_comparison for diagnostic context...", flush=True)
    pcr_unscaled = scib.metrics.pcr_comparison(
        adata_pre=adata,
        adata_post=adata,
        covariate=args.batch_key,
        embed=args.embedding_key,
        n_comps=50,
        scale=False,
        verbose=False,
    )
    print(f"  PCR (unscaled, before-after) = {pcr_unscaled:.6f}", flush=True)

    # --- Write sidecar CSV ---
    out_row = {
        "method": f"{args.method}_n{args.n_latent}",
        "compartment": args.compartment,
        "embedding_key": args.embedding_key,
        "n_cells": int(adata.n_obs),
        "n_genes_baseline": int(adata.n_vars),
        "pcr_scaled": float(pcr_val),
        "pcr_unscaled": float(pcr_unscaled),
        "pcr_elapsed_min": round(elapsed_min, 2),
        "scib_version": getattr(scib, "__version__", "unknown"),
    }
    pd.DataFrame([out_row]).to_csv(args.output_csv, index=False)
    print(f"\nSaved: {args.output_csv}", flush=True)
    print(f"Total: {(time.time() - t_total) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
