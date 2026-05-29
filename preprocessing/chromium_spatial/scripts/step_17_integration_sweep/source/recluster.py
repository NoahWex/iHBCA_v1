#!/usr/bin/env python3
"""
step_17_integration_sweep / recluster.py

Re-cluster an existing integration H5AD without re-running scVI. Loads
the integrated.h5ad, runs sc.pp.neighbors + sc.tl.leiden at user-specified
resolutions, and writes a flat sidecar CSV (cell_id + new <prefix>_r<res>
columns) to a patches directory next to the H5AD. The H5AD is never mutated.

This script is a side utility for Phase B — it exists so that downstream
Rmd reports can pick up additional Leiden resolutions on the winner
integration without rebuilding the full postprocess pipeline.

provenance: Spatial_HBCA_preprocessing/scripts/06c_recluster_integration.py
            (moved from root scripts/ to step_17 source/ per numbering_plan.md)

Modernizations over the source:
  - none required; the source already had --dry-run, validation, and clean
    argument handling. Kept as-is with a relocation header.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True,
                   help="Input integration h5ad (read-only)")
    p.add_argument("--patches-dir", required=True,
                   help="Output directory for the sidecar CSV "
                        "(typically {winner}/patches/)")
    p.add_argument("--patch-id", required=True,
                   help="Identifier for this patch; used as filename and "
                        "embedded in column names")
    p.add_argument("--use-rep", default="X_emb",
                   help="obsm key for neighbor graph (default: X_emb)")
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--resolutions", required=True,
                   help="Comma-separated list of leiden resolutions, "
                        "e.g. '0.4,0.6,1.5,2.0'")
    p.add_argument("--col-prefix", default=None,
                   help="Column name prefix; defaults to --patch-id")
    p.add_argument("--random-seed", type=int, default=0,
                   help="leiden random_state (default: 0)")
    p.add_argument("--dry-run", action="store_true",
                   help="Validate inputs, print plan, exit before clustering")
    args = p.parse_args()

    resolutions = [float(r) for r in args.resolutions.split(",") if r.strip()]
    col_prefix = args.col_prefix or args.patch_id

    print("=== step_17 recluster ===")
    print(f"Input h5ad:  {args.h5ad}")
    print(f"Patches dir: {args.patches_dir}")
    print(f"Patch id:    {args.patch_id}")
    print(f"use_rep:     {args.use_rep}")
    print(f"n_neighbors: {args.n_neighbors}")
    print(f"resolutions: {resolutions}")
    print(f"col_prefix:  {col_prefix}")
    print(f"random_seed: {args.random_seed}")

    if not os.path.exists(args.h5ad):
        sys.exit(f"FATAL: input h5ad not found: {args.h5ad}")

    if args.dry_run:
        print("VALIDATION PASSED (dry run)")
        sys.exit(0)

    out = Path(args.patches_dir)
    out.mkdir(parents=True, exist_ok=True)
    out_csv = out / f"{args.patch_id}.csv"

    print("Reading h5ad...")
    adata = ad.read_h5ad(args.h5ad)
    print(f"  shape: {adata.shape}")
    if args.use_rep not in adata.obsm.keys():
        sys.exit(f"FATAL: '{args.use_rep}' not in obsm ({list(adata.obsm.keys())})")

    print(f"Building neighbor graph (use_rep={args.use_rep}, "
          f"n_neighbors={args.n_neighbors})...")
    sc.pp.neighbors(adata, use_rep=args.use_rep, n_neighbors=args.n_neighbors)

    cluster_cols = {}
    for res in resolutions:
        col = f"{col_prefix}_r{res}"
        print(f"  leiden res={res} -> {col}")
        sc.tl.leiden(adata, resolution=res, key_added=col,
                     random_state=args.random_seed)
        n_clust = adata.obs[col].nunique()
        cluster_cols[col] = adata.obs[col].astype(str).values
        print(f"    clusters: {n_clust}")

    df = pd.DataFrame(cluster_cols)
    df.insert(0, "cell_id", list(adata.obs_names))
    df.to_csv(out_csv, index=False)
    print(f"Wrote {out_csv} ({df.shape[0]} rows x {df.shape[1] - 1} cluster cols)")

    prov = {
        "patch_id": args.patch_id,
        "source_h5ad": args.h5ad,
        "use_rep": args.use_rep,
        "n_neighbors": args.n_neighbors,
        "resolutions": resolutions,
        "col_prefix": col_prefix,
        "random_seed": args.random_seed,
        "n_cells": int(adata.n_obs),
        "cluster_cols": list(cluster_cols.keys()),
    }
    with open(out / f"{args.patch_id}.json", "w") as f:
        json.dump(prov, f, indent=2)
    print(f"Provenance: {out / (args.patch_id + '.json')}")
    print("Done.")


if __name__ == "__main__":
    main()
