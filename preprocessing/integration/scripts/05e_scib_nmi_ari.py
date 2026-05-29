#!/usr/bin/env python3
"""
NMI / ARI sidecar over pre-computed Leiden columns.

Scans the 10 leiden_* resolution columns saved by 01_compartment_scvi_full.py /
02_compartment_scanvi.py and computes NMI and ARI against the L1 label column
for each. Reports the best NMI and best ARI with their source resolutions.

Substitutes for scib's internal opt_louvain by sweeping the pre-computed
Leiden grid [0.1..5.0] directly — produces comparable values without a
second clustering pass.

Usage:
    python 05e_scib_nmi_ari.py \\
        --h5ad <{comp}_{method}.h5ad> \\
        --embedding-key <X_scVI|X_scANVI|...> \\
        --batch-key dataset \\
        --label-key level1_annotation \\
        --metadata <{comp}_metadata.csv> \\
        --method-name <method_label> \\
        --compartment Immune \\
        --output-csv <out_nmi_ari.csv>
"""

import argparse
import os
import re
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
    parser.add_argument("--metadata", required=True,
                        help="cell_metadata_enriched.csv for label join fallback")
    parser.add_argument("--method-name", required=True)
    parser.add_argument("--compartment", required=True,
                        choices=["Immune", "Epithelial", "Stromal"])
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
    os.makedirs(os.path.dirname(args.output_csv), exist_ok=True)

    print("=" * 70, flush=True)
    print(f"NMI/ARI sidecar: {args.method_name} ({args.compartment})", flush=True)
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

    # --- Label join (same pattern as 05b kBET sidecar) ---
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

    # --- Find leiden columns ---
    leiden_cols = {}
    for col in adata.obs.columns:
        m = re.match(r"^leiden_([\d.]+)$", col)
        if m:
            leiden_cols[float(m.group(1))] = col
    if not leiden_cols:
        raise ValueError(
            f"No leiden_* columns found in obs. "
            f"Available: {list(adata.obs.columns)}"
        )
    print(f"  Found {len(leiden_cols)} leiden resolutions: "
          f"{sorted(leiden_cols.keys())}", flush=True)

    # --- NMI/ARI scan ---
    # scib.metrics.nmi uses average_method='arithmetic' by default (matches
    # the internal call in scib.metrics.metrics()). ARI is sklearn's
    # adjusted_rand_score via scib.metrics.ari.
    print(f"\nScanning NMI/ARI across {len(leiden_cols)} resolutions...",
          flush=True)
    t0 = time.time()
    best_nmi, best_ari = -1.0, -1.0
    best_nmi_res, best_ari_res = None, None
    for res in sorted(leiden_cols.keys()):
        col = leiden_cols[res]
        adata.obs[col] = adata.obs[col].astype(str).astype("category")
        nmi_val = scib.metrics.nmi(adata, col, args.label_key)
        ari_val = scib.metrics.ari(adata, col, args.label_key)
        n_cl = adata.obs[col].nunique()
        print(f"  res={res:5.1f}  n_clusters={n_cl:4d}  "
              f"NMI={nmi_val:.6f}  ARI={ari_val:.6f}", flush=True)
        if nmi_val > best_nmi:
            best_nmi, best_nmi_res = nmi_val, res
        if ari_val > best_ari:
            best_ari, best_ari_res = ari_val, res
    elapsed_min = (time.time() - t0) / 60
    print(f"\n  Best NMI: {best_nmi:.6f} (res={best_nmi_res})", flush=True)
    print(f"  Best ARI: {best_ari:.6f} (res={best_ari_res})", flush=True)
    print(f"  Elapsed: {elapsed_min:.1f} min", flush=True)

    # --- Write sidecar CSV ---
    out_row = {
        "method": args.method_name,
        "compartment": args.compartment,
        "embedding_key": args.embedding_key,
        "n_cells": int(adata.n_obs),
        "NMI": best_nmi,
        "ARI": best_ari,
        "NMI_best_res": best_nmi_res,
        "ARI_best_res": best_ari_res,
        "nmi_ari_elapsed_min": round(elapsed_min, 2),
        "scib_version": scib.__version__,
    }
    pd.DataFrame([out_row]).to_csv(args.output_csv, index=False)
    print(f"\nSaved: {args.output_csv}", flush=True)
    print(f"Total: {(time.time() - t_total) / 60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
