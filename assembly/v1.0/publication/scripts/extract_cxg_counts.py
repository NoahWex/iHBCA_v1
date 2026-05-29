#!/usr/bin/env python3
"""
extract_cxg_counts.py - Extract raw integer counts from CxG h5ad files
=======================================================================
Reads the original CELLxGENE h5ad for a study and writes intermediates
(counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz) with integer counts.

Fixes A8: gray/twigger Seurat RDS extraction produced log-normalized floats
instead of raw counts. CxG h5ads have verified integer counts in raw.X.

Usage:
  python extract_cxg_counts.py \
      --cxg-h5ad /path/to/gray.h5ad \
      --output-dir /path/to/intermediates/gray \
      [--dry-run]

"""

import argparse
import gzip
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import scipy.io
import scipy.sparse as sp


def inspect_h5ad(path):
    """Inspect h5ad structure and report count matrix location + dtype."""
    print(f"\n--- Inspecting: {path} ---")

    with h5py.File(path, "r") as f:
        # Check X
        if "X" in f:
            x_obj = f["X"]
            if isinstance(x_obj, h5py.Group) and "data" in x_obj:
                dtype = x_obj["data"].dtype
                vals = x_obj["data"][:100]
                is_int = np.all(vals == vals.astype(int))
                print(f"  X: dtype={dtype}, all_int={is_int}, sample_range=[{vals.min()}, {vals.max()}]")
            elif isinstance(x_obj, h5py.Dataset):
                print(f"  X: dense dataset, dtype={x_obj.dtype}")

        # Check raw.X
        if "raw" in f and "X" in f["raw"]:
            rx = f["raw"]["X"]
            if isinstance(rx, h5py.Group) and "data" in rx:
                dtype = rx["data"].dtype
                vals = rx["data"][:100]
                is_int = np.all(vals == vals.astype(int))
                print(f"  raw.X: dtype={dtype}, all_int={is_int}, sample_range=[{vals.min()}, {vals.max()}]")

        # Check var structure
        if "var" in f:
            var_cols = list(f["var"].keys())
            print(f"  var columns: {var_cols}")
            if "_index" in f["var"]:
                idx = f["var"]["_index"][:3]
                print(f"  var index sample: {[x.decode() if isinstance(x, bytes) else x for x in idx]}")
        if "raw" in f and "var" in f["raw"]:
            raw_var_cols = list(f["raw"]["var"].keys())
            print(f"  raw.var columns: {raw_var_cols}")

        # Check obs
        if "obs" in f:
            if "_index" in f["obs"]:
                idx = f["obs"]["_index"][:3]
                print(f"  obs index sample: {[x.decode() if isinstance(x, bytes) else x for x in idx]}")


def extract_counts(cxg_path, output_dir, dry_run=False):
    """Extract raw integer counts from CxG h5ad to intermediates format."""

    print(f"\n{'=' * 70}")
    print(f"EXTRACTING RAW COUNTS FROM CxG H5AD")
    print(f"  Source: {cxg_path}")
    print(f"  Output: {output_dir}")
    print(f"{'=' * 70}")

    # First inspect
    inspect_h5ad(cxg_path)

    # Determine count matrix source: prefer raw.X (unfiltered genes), fall back to X
    print(f"\n[1/4] Loading count matrix...")
    adata = ad.read_h5ad(cxg_path, backed="r")

    # Determine which matrix has integer counts
    use_raw = False
    if adata.raw is not None:
        # Check raw.X dtype via h5py (faster than loading full matrix)
        with h5py.File(cxg_path, "r") as f:
            if "raw" in f and "X" in f["raw"]:
                rx = f["raw"]["X"]
                if isinstance(rx, h5py.Group) and "data" in rx:
                    raw_dtype = rx["data"].dtype
                    raw_vals = rx["data"][:1000]
                    if np.issubdtype(raw_dtype, np.integer) or np.all(raw_vals == raw_vals.astype(int)):
                        use_raw = True
                        print(f"  Using raw.X (dtype={raw_dtype}, integer counts)")

    if not use_raw:
        with h5py.File(cxg_path, "r") as f:
            x_obj = f["X"]
            if isinstance(x_obj, h5py.Group) and "data" in x_obj:
                x_dtype = x_obj["data"].dtype
                x_vals = x_obj["data"][:1000]
                if np.issubdtype(x_dtype, np.integer) or np.all(x_vals == x_vals.astype(int)):
                    print(f"  Using X (dtype={x_dtype}, integer counts)")
                else:
                    print(f"  ERROR: Neither X nor raw.X contains integer counts!")
                    print(f"    X dtype={x_dtype}, sample values: {x_vals[:5]}")
                    sys.exit(1)

    adata.file.close()

    # Now load fully (not backed) to get the matrix
    print(f"  Loading full h5ad into memory...")
    adata = ad.read_h5ad(cxg_path)

    if use_raw:
        mat = adata.raw.X
        var_df = adata.raw.var
    else:
        mat = adata.X
        var_df = adata.var

    # Ensure sparse CSR
    if not sp.issparse(mat):
        mat = sp.csr_matrix(mat)
    elif not sp.isspmatrix_csr(mat):
        mat = mat.tocsr()

    # Ensure integer dtype
    if not np.issubdtype(mat.dtype, np.integer):
        if np.all(mat.data == mat.data.astype(int)):
            print(f"  Converting float → int (values are integer-valued)")
            mat.data = mat.data.astype(np.int64)
        else:
            print(f"  ERROR: Matrix values are not integers!")
            print(f"    dtype={mat.dtype}, sample: {mat.data[:10]}")
            sys.exit(1)

    n_cells, n_genes = mat.shape
    print(f"  Shape: {n_cells:,} cells x {n_genes:,} genes")
    print(f"  Non-zero: {mat.nnz:,}")
    print(f"  Dtype: {mat.dtype}")
    print(f"  Value range: [{mat.data.min()}, {mat.data.max()}]")

    # Get gene names (symbols preferred for compatibility with pipeline)
    print(f"\n[2/4] Extracting gene names...")
    # CxG h5ads typically have var.index = Ensembl IDs, with feature_name column for symbols
    gene_names = None
    for col in ["feature_name", "gene_symbol", "gene_name"]:
        if col in var_df.columns:
            gene_names = var_df[col].values.tolist()
            print(f"  Using var['{col}'] as gene names")
            break

    if gene_names is None:
        # Fall back to var index
        gene_names = var_df.index.tolist()
        print(f"  Using var.index as gene names")

    print(f"  Gene name sample: {gene_names[:3]}")
    print(f"  Total genes: {len(gene_names)}")

    # Get cell barcodes
    print(f"\n[3/4] Extracting cell barcodes...")
    cell_ids = adata.obs.index.tolist()
    print(f"  Cell ID sample: {cell_ids[:3]}")
    print(f"  Total cells: {len(cell_ids)}")

    if dry_run:
        print(f"\n  [DRY RUN] Would write to {output_dir}")
        return

    # Write intermediates
    print(f"\n[4/4] Writing intermediates...")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write to /tmp first (CRSP stale file handle avoidance)
    import tempfile
    import shutil
    tmpdir = Path(tempfile.mkdtemp(prefix="cxg_extract_"))

    # counts.mtx.gz — genes x cells (transpose for R convention)
    mat_t = mat.T.tocoo()  # genes x cells, COO for MatrixMarket
    tmp_mtx = tmpdir / "counts.mtx"
    scipy.io.mmwrite(str(tmp_mtx), mat_t, field="integer")
    print(f"  Wrote counts.mtx ({tmp_mtx.stat().st_size:,} bytes)")

    import subprocess
    subprocess.run(["gzip", "-f", str(tmp_mtx)], check=True)
    tmp_mtx_gz = tmpdir / "counts.mtx.gz"
    dst_mtx = output_dir / "counts.mtx.gz"
    shutil.copy(str(tmp_mtx_gz), str(dst_mtx))
    print(f"  Compressed + copied: {dst_mtx} ({dst_mtx.stat().st_size:,} bytes)")

    # features.tsv.gz
    tmp_feat = tmpdir / "features.tsv.gz"
    with gzip.open(str(tmp_feat), "wt") as f:
        for name in gene_names:
            f.write(f"{name}\n")
    shutil.copy(str(tmp_feat), str(output_dir / "features.tsv.gz"))
    print(f"  Wrote features.tsv.gz ({len(gene_names)} genes)")

    # barcodes.tsv.gz
    tmp_bc = tmpdir / "barcodes.tsv.gz"
    with gzip.open(str(tmp_bc), "wt") as f:
        for cid in cell_ids:
            f.write(f"{cid}\n")
    shutil.copy(str(tmp_bc), str(output_dir / "barcodes.tsv.gz"))
    print(f"  Wrote barcodes.tsv.gz ({len(cell_ids)} cells)")

    # Cleanup tmpdir
    shutil.rmtree(str(tmpdir))

    # Verification
    print(f"\n--- Verification ---")
    verify_mat = scipy.io.mmread(str(dst_mtx))
    verify_mat = sp.csr_matrix(verify_mat.T)
    print(f"  Re-read shape: {verify_mat.shape}")
    print(f"  Re-read dtype: {verify_mat.dtype}")
    print(f"  Re-read nnz: {verify_mat.nnz:,}")
    print(f"  Values are integers: {np.all(verify_mat.data == verify_mat.data.astype(int))}")

    if verify_mat.shape != (n_cells, n_genes):
        print(f"  ERROR: Shape mismatch! Expected ({n_cells}, {n_genes}), got {verify_mat.shape}")
        sys.exit(1)
    if verify_mat.nnz != mat.nnz:
        print(f"  ERROR: NNZ mismatch! Expected {mat.nnz:,}, got {verify_mat.nnz:,}")
        sys.exit(1)

    print(f"\n{'=' * 70}")
    print(f"EXTRACTION COMPLETE")
    print(f"  Cells: {n_cells:,}")
    print(f"  Genes: {n_genes:,}")
    print(f"  Output: {output_dir}")
    print(f"{'=' * 70}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract raw integer counts from CxG h5ad to intermediates format"
    )
    parser.add_argument("--cxg-h5ad", required=True, help="Path to original CxG h5ad")
    parser.add_argument("--output-dir", required=True, help="Output intermediates directory")
    parser.add_argument("--dry-run", action="store_true", help="Inspect only, don't write")
    args = parser.parse_args()

    extract_counts(args.cxg_h5ad, args.output_dir, args.dry_run)


if __name__ == "__main__":
    main()
