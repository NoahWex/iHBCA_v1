#!/usr/bin/env python3
# Step 20b — Promote full-target scvi_n100 to integration_intermediate/full/scvi_n100/
#
# Extracts latent + obs + umap from the step 17 sweep integrated.h5ad (canonical
# 258K cells, 5-way filtered). Writes the minimal schema consumed by Milo build.
# Counts matrix is NOT written here — can be added later via h5ad extraction if
# needed for marker work.
#
# Pattern: minimal subset of 20a_concat_full_object.py; avoids the concat path
# because sweep full target has 258,317 cells vs compartment-sum 256,858.

import argparse
import json
import os
import shutil
import sys
from datetime import date

import anndata as ad
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sweep-h5ad", required=True,
                        help="sweep/full_object/scvi_n100/integrated.h5ad")
    parser.add_argument("--winner-dir", required=True,
                        help="17_ScviIntegration/winner/full_object/")
    parser.add_argument("--out-dir", required=True,
                        help="integration_intermediate/full/scvi_n100/")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print(f"Loading h5ad: {args.sweep_h5ad}")
    adata = ad.read_h5ad(args.sweep_h5ad, backed="r")
    n_cells, n_genes = adata.shape
    print(f"  {n_cells:,} cells x {n_genes:,} genes")
    print(f"  obsm keys: {list(adata.obsm.keys())}")

    latent_key = "X_emb"
    if latent_key not in adata.obsm:
        sys.exit(f"ERROR: {latent_key} not in obsm")
    latent = adata.obsm[latent_key]
    print(f"  latent shape: {latent.shape}")

    # latent.csv: cell_id index + dims
    latent_df = pd.DataFrame(
        latent,
        index=adata.obs_names,
        columns=[f"emb_{i+1}" for i in range(latent.shape[1])],
    )
    latent_df.index.name = "cell_id"
    latent_path = os.path.join(args.out_dir, "latent.csv")
    latent_df.to_csv(latent_path)
    print(f"Wrote: {latent_path}")

    # obs.csv: use metadata_with_umap.csv from winner (has obs + leiden_multi + UMAP)
    # Add position column derived from sample_id (strip "{patient_id}_" prefix).
    # Required by downstream Milo summary + Fig 4/5 anatomic DA contrasts.
    meta_src = os.path.join(args.winner_dir, "metadata_with_umap.csv")
    obs_path = os.path.join(args.out_dir, "obs.csv")
    obs = pd.read_csv(meta_src)
    obs["position"] = [
        sid[len(pid) + 1:] if sid.startswith(pid + "_") else sid
        for sid, pid in zip(obs["sample_id"].astype(str), obs["patient_id"].astype(str))
    ]
    obs.to_csv(obs_path, index=False)
    print(f"Wrote: {obs_path}  (added position column)")

    # umap.csv: copy from winner
    umap_src = os.path.join(args.winner_dir, "umap.csv")
    umap_path = os.path.join(args.out_dir, "umap.csv")
    shutil.copy(umap_src, umap_path)
    print(f"Wrote: {umap_path}")

    # manifest.json
    obs_df = pd.read_csv(obs_path)
    leiden_cols = [c for c in obs_df.columns if c.startswith("leiden_")]
    leiden_res = sorted(float(c.replace("leiden_", "")) for c in leiden_cols)

    manifest = {
        "source": "step 17 full-target sweep winner (MANUAL_OVERRIDE=scvi_n100)",
        "date": str(date.today()),
        "config_label": "scvi_n100",
        "n_cells": int(n_cells),
        "n_genes": int(n_genes),
        "latent_key": latent_key,
        "latent_dims": int(latent.shape[1]),
        "umap_key": "X_umap",
        "obs_columns": list(obs_df.columns),
        "leiden_resolutions": leiden_res,
        "cell_filter": "canonical 5-way (step_01_umi, step_02_mad, step_03_doublet, step_10b_manifold, step_15)",
        "counts_status": "not yet extracted — can be regenerated from sweep_h5ad if needed",
    }

    manifest_path = os.path.join(args.out_dir, "manifest.json")
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote: {manifest_path}")

    print(f"\nDone. Promoted to: {args.out_dir}")


if __name__ == "__main__":
    main()
