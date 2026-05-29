#!/usr/bin/env python3
"""
Export a compartment integration h5ad to a flat intermediate directory
that R can consume directly (no reticulate).

Outputs (under --output-dir):
    counts.mtx.gz   # cells x genes sparse, Matrix::readMM after gunzip
    cells.tsv       # cell IDs (one per line, no header)
    genes.tsv       # gene names (one per line, no header)
    obs.csv         # full obs frame, cell_id as first column
    latent.csv      # latent embedding (cell_id + latent_1..latent_N)
    umap.csv        # umap embedding (cell_id + UMAP_1, UMAP_2)
    manifest.json   # source path, dimensions, key metadata for downstream sanity

The R Rmd consumes this directory and never touches Python.
"""
import argparse
import gzip
import json
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse


LATENT_CANDIDATES = ("X_emb", "X_scVI", "X_scvi", "X_latent", "X_integrated")
UMAP_CANDIDATES = ("X_umap", "X_UMAP")


def pick_obsm(adata, candidates):
    for k in candidates:
        if k in adata.obsm.keys():
            return k
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True, help="Input integration h5ad")
    p.add_argument("--output-dir", required=True, help="Intermediate output dir")
    p.add_argument("--compartment", required=True, help="Compartment label")
    p.add_argument("--config-label", required=True, help="Integration config label")
    p.add_argument("--dry-run", action="store_true",
                   help="Validate inputs and exit before reading h5ad")
    args = p.parse_args()

    print(f"=== 08a_export_h5ad_intermediate ===")
    print(f"Input:       {args.h5ad}")
    print(f"Output dir:  {args.output_dir}")
    print(f"Compartment: {args.compartment}")
    print(f"Config:      {args.config_label}")

    if not os.path.exists(args.h5ad):
        sys.exit(f"FATAL: input h5ad not found: {args.h5ad}")

    if args.dry_run:
        print("VALIDATION PASSED (dry run)")
        sys.exit(0)

    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Reading h5ad...")
    adata = ad.read_h5ad(args.h5ad)
    print(f"  shape: {adata.shape} (cells x genes)")
    print(f"  layers: {list(adata.layers.keys())}")
    print(f"  obsm: {list(adata.obsm.keys())}")
    print(f"  obs columns: {list(adata.obs.columns)}")

    # Counts: prefer layers['counts'] if present, else .X (verified raw in v2 sweep)
    if "counts" in adata.layers.keys():
        print("Using layers['counts'] as raw counts")
        X = adata.layers["counts"]
    else:
        print("Using .X as counts")
        X = adata.X

    if not scipy.sparse.issparse(X):
        X = scipy.sparse.csr_matrix(X)

    # Sanity: counts must be integer-valued. Fail loudly if normalized.
    sample = X.data[:1000] if X.nnz > 1000 else X.data
    if not np.allclose(sample, np.round(sample)):
        sys.exit("FATAL: counts contain non-integer values - .X appears normalized")
    print(f"Counts integer-check: PASS (nnz={X.nnz})")

    # Pick latent + umap
    latent_key = pick_obsm(adata, LATENT_CANDIDATES)
    umap_key = pick_obsm(adata, UMAP_CANDIDATES)
    if latent_key is None:
        sys.exit(f"FATAL: no latent obsm found among {LATENT_CANDIDATES}")
    if umap_key is None:
        sys.exit(f"FATAL: no umap obsm found among {UMAP_CANDIDATES}")
    print(f"Latent: {latent_key} ({adata.obsm[latent_key].shape})")
    print(f"UMAP:   {umap_key} ({adata.obsm[umap_key].shape})")

    # Write counts MTX (gzipped)
    mtx_path = out / "counts.mtx"
    print(f"Writing {mtx_path}.gz ...")
    scipy.io.mmwrite(str(mtx_path), X.tocoo())
    with open(mtx_path, "rb") as f_in, gzip.open(str(mtx_path) + ".gz", "wb") as f_out:
        f_out.writelines(f_in)
    os.remove(mtx_path)

    # Cells, genes
    cell_ids = list(adata.obs_names)
    gene_ids = list(adata.var_names)
    (out / "cells.tsv").write_text("\n".join(cell_ids) + "\n")
    (out / "genes.tsv").write_text("\n".join(gene_ids) + "\n")
    print(f"  cells.tsv: {len(cell_ids)}")
    print(f"  genes.tsv: {len(gene_ids)}")

    # obs (with cell_id as first column)
    obs_out = adata.obs.copy()
    obs_out.insert(0, "cell_id", cell_ids)
    obs_out.to_csv(out / "obs.csv", index=False)
    print(f"  obs.csv: {obs_out.shape}")

    # Latent
    latent = np.asarray(adata.obsm[latent_key])
    latent_df = pd.DataFrame(
        latent,
        columns=[f"latent_{i+1}" for i in range(latent.shape[1])],
    )
    latent_df.insert(0, "cell_id", cell_ids)
    latent_df.to_csv(out / "latent.csv", index=False)
    print(f"  latent.csv: {latent_df.shape}")

    # UMAP
    umap = np.asarray(adata.obsm[umap_key])
    umap_df = pd.DataFrame(
        umap,
        columns=[f"UMAP_{i+1}" for i in range(umap.shape[1])],
    )
    umap_df.insert(0, "cell_id", cell_ids)
    umap_df.to_csv(out / "umap.csv", index=False)
    print(f"  umap.csv: {umap_df.shape}")

    # Manifest
    manifest = {
        "source_h5ad": str(args.h5ad),
        "compartment": args.compartment,
        "config_label": args.config_label,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "latent_key": latent_key,
        "umap_key": umap_key,
        "obs_columns": list(adata.obs.columns),
    }
    with open(out / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"  manifest.json written")
    print("Done.")


if __name__ == "__main__":
    main()
