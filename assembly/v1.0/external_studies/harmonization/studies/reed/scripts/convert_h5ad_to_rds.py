#!/usr/bin/env python3
"""
convert_h5ad_to_rds.py - Phase 1: Extract H5AD components
==========================================================
Extracts count matrix, metadata, and embeddings from CELLxGENE H5AD
into intermediate files for R-side Seurat construction.

Outputs (in --output-dir):
  counts.mtx.gz    - sparse count matrix (Market Matrix format)
  barcodes.tsv.gz  - cell barcodes
  features.tsv.gz  - gene names
  metadata.csv     - cell-level metadata from .obs
  pca.csv          - PCA embedding if present

Usage:
  python convert_h5ad_to_rds.py --input reed.h5ad --output-dir /tmp/reed_extract
"""

# The H5AD walks into a bar. The bartender says "We don't serve your type here."
# The H5AD replies: "That's fine, I'm converting."

import argparse
import os
import sys
import gzip

import scanpy as sc
import scipy.io
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Extract H5AD to intermediate files")
    parser.add_argument("--input", required=True, help="Path to H5AD file")
    parser.add_argument("--output-dir", required=True, help="Output directory for intermediates")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"FATAL: Input file not found: {args.input}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("PHASE 1: H5AD EXTRACTION")
    print("=" * 60)

    # Load H5AD
    print(f"\n[1] Loading H5AD: {args.input}")
    adata = sc.read_h5ad(args.input)
    print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")
    print(f"  Obs columns: {list(adata.obs.columns)[:10]}...")
    print(f"  Obsm keys: {list(adata.obsm.keys())}")

    # Extract raw counts (prefer .raw if available, else .X)
    print("\n[2] Extracting count matrix...")
    if adata.raw is not None:
        print("  Using adata.raw.X")
        counts = adata.raw.X
        genes = list(adata.raw.var_names)
    else:
        print("  Using adata.X (no .raw layer)")
        counts = adata.X
        genes = list(adata.var_names)

    from scipy.sparse import issparse
    if not issparse(counts):
        from scipy.sparse import csr_matrix
        counts = csr_matrix(counts)
        print("  Converted dense to sparse")

    # Transpose: scanpy is cells x genes, writeMM expects genes x cells for R
    counts_t = counts.T.tocsc()
    print(f"  Matrix: {counts_t.shape[0]} genes x {counts_t.shape[1]} cells")

    # Write to /tmp first (CRSP stale file handle workaround)
    tmp_mtx = os.path.join(args.output_dir, "counts.mtx")
    scipy.io.mmwrite(tmp_mtx, counts_t)
    print(f"  Wrote: counts.mtx")

    # Gzip
    with open(tmp_mtx, 'rb') as f_in:
        with gzip.open(tmp_mtx + ".gz", 'wb') as f_out:
            f_out.writelines(f_in)
    os.remove(tmp_mtx)
    print(f"  Compressed: counts.mtx.gz")

    # Barcodes
    print("\n[3] Writing barcodes...")
    barcodes = list(adata.obs_names)
    bc_path = os.path.join(args.output_dir, "barcodes.tsv.gz")
    with gzip.open(bc_path, 'wt') as f:
        f.write("\n".join(barcodes) + "\n")
    print(f"  Wrote: {len(barcodes)} barcodes")

    # Features
    print("\n[4] Writing features...")
    feat_path = os.path.join(args.output_dir, "features.tsv.gz")
    with gzip.open(feat_path, 'wt') as f:
        for g in genes:
            f.write(f"{g}\t{g}\tGene Expression\n")
    print(f"  Wrote: {len(genes)} features")

    # Metadata
    print("\n[5] Writing metadata...")
    meta_path = os.path.join(args.output_dir, "metadata.csv")
    adata.obs.to_csv(meta_path)
    print(f"  Wrote: metadata.csv ({len(adata.obs.columns)} columns)")

    # PCA embedding
    pca_keys = [k for k in adata.obsm.keys() if "pca" in k.lower() or k == "X_pca"]
    if pca_keys:
        print(f"\n[6] Writing PCA embedding ({pca_keys[0]})...")
        pca = adata.obsm[pca_keys[0]]
        pca_df = pd.DataFrame(pca, index=adata.obs_names,
                              columns=[f"PC_{i+1}" for i in range(pca.shape[1])])
        pca_path = os.path.join(args.output_dir, "pca.csv")
        pca_df.to_csv(pca_path)
        print(f"  Wrote: {pca.shape[0]} x {pca.shape[1]} PCA embedding")
    else:
        print("\n[6] No PCA embedding found in obsm")

    print("\n" + "=" * 60)
    print("PHASE 1 COMPLETE")
    print(f"  Cells: {adata.shape[0]}")
    print(f"  Genes: {len(genes)}")
    print(f"  Output: {args.output_dir}")
    print("=" * 60)


if __name__ == "__main__":
    main()
