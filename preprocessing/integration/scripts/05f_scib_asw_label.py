#!/usr/bin/env python3
"""
ASW_label sidecar — bio silhouette for one integrated h5ad.

Computes scib.metrics.silhouette(adata, label_key, embed=embedding_key) as a
dedicated SLURM task. ASW is O(N^2) in cell count and is the slowest metric
in the scIB suite at million-cell scale (≈3.5h on the 991K-cell epithelial
compartment), so it is split into its own task to allow per-compartment
walltime budgeting and failure isolation.

Usage:
    python 05f_scib_asw_label.py \\
        --h5ad <{comp}_{method}.h5ad> \\
        --embedding-key <X_scVI|X_scANVI|...> \\
        --batch-key dataset \\
        --label-key level1_annotation \\
        --metadata <{comp}_metadata.csv> \\
        --method-name <method_label> \\
        --compartment Immune \\
        --output-csv <out_asw_label.csv>
"""

import argparse
import os
import sys
import time

# Numba cache patch — required before any scib/scanpy import in the
# read-only container.
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

import anndata as ad  # type: ignore[import-not-found]
import pandas as pd
import scib  # type: ignore[import-not-found]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--embedding-key", required=True)
    parser.add_argument("--batch-key", default="dataset")
    parser.add_argument("--label-key", default="level1_annotation")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--method-name", required=True)
    parser.add_argument("--compartment", required=True,
                        choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    print("=" * 70, flush=True)
    print(f"ASW_label sidecar: {args.method_name} ({args.compartment})", flush=True)
    print(f"scib version: {scib.__version__}", flush=True)
    print("=" * 70, flush=True)
    t_total = time.time()

    # --- Load h5ad ---
    print(f"Loading {args.h5ad}...", flush=True)
    t0 = time.time()
    adata = ad.read_h5ad(args.h5ad)
    print(f"  Shape: {adata.n_obs:,} x {adata.n_vars:,} in {time.time() - t0:.1f}s",
          flush=True)

    if args.embedding_key not in adata.obsm:
        raise ValueError(
            f"Embedding key '{args.embedding_key}' not in obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )

    # --- Label join ---
    if args.label_key not in adata.obs.columns:
        print(f"  Label key '{args.label_key}' missing from obs — "
              f"joining from {args.metadata}", flush=True)
        ext_meta = pd.read_csv(
            args.metadata,
            usecols=["cell_id", "numeric_id", args.label_key],
            low_memory=False,
        )
        ext_meta["numeric_id"] = ext_meta["numeric_id"].astype(str)
        ext_meta["cell_id"] = ext_meta["cell_id"].astype(str)
        obs_names_str = adata.obs_names.astype(str)
        for join_col in ("numeric_id", "cell_id"):
            lookup = ext_meta.set_index(join_col)[args.label_key]
            joined = obs_names_str.map(lookup)
            n_matched = int(joined.notna().sum())
            print(f"  Trying join on '{join_col}': "
                  f"{n_matched:,} / {len(obs_names_str):,} matched", flush=True)
            if n_matched > 0.95 * len(obs_names_str):
                adata.obs[args.label_key] = joined.values
                break
        else:
            raise ValueError("Label fallback failed: no valid join column")
        n_null = int(adata.obs[args.label_key].isna().sum())
        if n_null > 0:
            print(f"  Dropping {n_null:,} null-label cells", flush=True)
            adata = adata[~adata.obs[args.label_key].isna()].copy()

    print(f"  Labels: {adata.obs[args.label_key].nunique()} unique", flush=True)

    # --- Compute ASW_label ---
    print(f"\nComputing ASW_label (scib.metrics.silhouette)...", flush=True)
    t0 = time.time()
    asw_label = float(scib.metrics.silhouette(
        adata, label_key=args.label_key, embed=args.embedding_key))
    elapsed_min = (time.time() - t0) / 60
    print(f"  ASW_label = {asw_label:.6f} ({elapsed_min:.1f} min)", flush=True)

    # --- Write sidecar CSV ---
    out_row = {
        "method": args.method_name,
        "compartment": args.compartment,
        "embedding_key": args.embedding_key,
        "n_cells": int(adata.n_obs),
        "ASW_label": asw_label,
        "asw_label_elapsed_min": round(elapsed_min, 2),
        "scib_version": scib.__version__,
    }
    pd.DataFrame([out_row]).to_csv(args.output_csv, index=False)
    print(f"\nSaved: {args.output_csv}", flush=True)
    print(f"Total: {(time.time() - t_total) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
