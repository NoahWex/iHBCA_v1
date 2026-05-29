#!/usr/bin/env python3
"""
graph_connectivity sidecar — kNN-graph connectivity for one integrated h5ad.

Computes scib.metrics.graph_connectivity(adata, label_key). Requires a
neighbors graph (obsp['connectivities']) computed on the integration
embedding; this script builds that graph via sc.pp.neighbors before the
metric call.

Usage:
    python 05h_scib_graph_conn.py \\
        --h5ad <{comp}_{method}.h5ad> \\
        --embedding-key <X_scVI|X_scANVI|...> \\
        --batch-key dataset \\
        --label-key level1_annotation \\
        --metadata <{comp}_metadata.csv> \\
        --method-name <method_label> \\
        --compartment Immune \\
        --output-csv <out_graph_conn.csv>
"""

import argparse
import os
import sys
import time

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
import scanpy as sc  # type: ignore[import-not-found]
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
    parser.add_argument("--n-neighbors", type=int, default=15,
                        help="k for sc.pp.neighbors (default 15, matches "
                             "scib.metrics.metrics() internal default)")
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    print("=" * 70, flush=True)
    print(f"graph_conn sidecar: {args.method_name} ({args.compartment})",
          flush=True)
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

    # scib.metrics.graph_connectivity calls adata.obs[label_key].cat.categories,
    # so the label column must be pandas Categorical. After a fallback join
    # from CSV it lands as object dtype — cast it.
    if not isinstance(adata.obs[args.label_key].dtype, pd.CategoricalDtype):
        adata.obs[args.label_key] = adata.obs[args.label_key].astype("category")

    print(f"  Labels: {adata.obs[args.label_key].nunique()} unique", flush=True)

    # --- Build neighbors graph on the integrated embedding ---
    # graph_connectivity reads adata.obsp['connectivities'], which is
    # populated by sc.pp.neighbors. n_neighbors=15 matches the internal
    # default used by scib.metrics.metrics().
    print(f"\nBuilding scanpy neighbors (k={args.n_neighbors}) on "
          f"{args.embedding_key}...", flush=True)
    t0 = time.time()
    sc.pp.neighbors(adata, use_rep=args.embedding_key,
                    n_neighbors=args.n_neighbors)
    print(f"  neighbors built in {(time.time() - t0) / 60:.1f} min", flush=True)

    # --- Compute graph_conn ---
    print(f"\nComputing graph_connectivity...", flush=True)
    t0 = time.time()
    graph_conn = float(scib.metrics.graph_connectivity(
        adata, label_key=args.label_key))
    elapsed_min = (time.time() - t0) / 60
    print(f"  graph_conn = {graph_conn:.6f} ({elapsed_min:.1f} min)", flush=True)

    # --- Write sidecar CSV ---
    out_row = {
        "method": args.method_name,
        "compartment": args.compartment,
        "embedding_key": args.embedding_key,
        "n_cells": int(adata.n_obs),
        "graph_conn": graph_conn,
        "graph_conn_elapsed_min": round(elapsed_min, 2),
        "n_neighbors": args.n_neighbors,
        "scib_version": scib.__version__,
    }
    pd.DataFrame([out_row]).to_csv(args.output_csv, index=False)
    print(f"\nSaved: {args.output_csv}", flush=True)
    print(f"Total: {(time.time() - t_total) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
