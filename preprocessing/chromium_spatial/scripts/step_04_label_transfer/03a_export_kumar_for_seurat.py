#!/usr/bin/env python3
"""
Export Kumar reference from NPZ counts + metadata for Seurat FindTransferAnchors.

Subsets the 2.12M-cell NPZ to dataset=="kumar" (714K cells), writes MTX + CSV.
Much faster than loading the 23GB h5ad.

Usage:
    python 03a_export_kumar_for_seurat.py \
        --counts-npz /path/to/preintegration_HBCA_inner_ENSEMBL_simple_counts_matrix.npz \
        --gene-data /path/to/gene_data.csv \
        --metadata /path/to/cell_metadata_aligned.csv \
        --output-dir outputs/kumar_reference/
"""

try:
    import numba
    _o = numba.njit
    def _nc(*a, **k): k.pop("cache", None); return _o(*a, **k)
    numba.njit = _nc
    _ov = numba.vectorize
    def _ncv(*a, **k): k.pop("cache", None); return _ov(*a, **k)
    numba.vectorize = _ncv
except ImportError:
    pass

import argparse
import gzip
import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts-npz", required=True)
    parser.add_argument("--gene-data", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--label-col", default="native_kumar")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Load NPZ
    print(f"Loading NPZ: {args.counts_npz}")
    t0 = time.time()
    data = np.load(args.counts_npz, allow_pickle=True)
    counts = sp.csr_matrix((data["data"], data["indices"], data["indptr"]),
                            shape=tuple(data["shape"]))
    print(f"  {counts.shape[0]} cells x {counts.shape[1]} genes, {time.time()-t0:.1f}s")

    # Load gene data
    gene_df = pd.read_csv(args.gene_data, index_col=0)
    print(f"  Genes: {len(gene_df)}")

    # Load metadata — subset to Kumar
    print("Loading metadata...")
    meta = pd.read_csv(args.metadata,
                       usecols=["numeric_id", "dataset", args.label_col],
                       low_memory=False)
    kumar_mask = meta["dataset"] == "kumar"
    kumar_ids = meta.loc[kumar_mask, "numeric_id"].values
    kumar_labels = meta.loc[kumar_mask, args.label_col].values
    print(f"  Kumar cells: {len(kumar_ids)}")

    # Subset counts
    print("Subsetting counts to Kumar...")
    kumar_counts = counts[kumar_ids]
    print(f"  Kumar counts: {kumar_counts.shape}")

    # Filter out cells with no label
    has_label = pd.notna(kumar_labels) & (kumar_labels != "") & (kumar_labels.astype(str) != "nan")
    kumar_counts = kumar_counts[has_label]
    kumar_labels = kumar_labels[has_label]
    kumar_ids = kumar_ids[has_label]
    print(f"  With labels: {len(kumar_ids)} ({np.unique(kumar_labels).shape[0]} unique types)")

    # Write MTX (genes x cells for Seurat)
    print("Writing MTX...")
    t0 = time.time()
    kumar_counts_t = kumar_counts.T.tocsc()
    mtx_path = os.path.join(args.output_dir, "kumar_sc_counts.mtx")
    mmwrite(mtx_path, kumar_counts_t)
    with open(mtx_path, "rb") as f_in:
        with gzip.open(mtx_path + ".gz", "wb") as f_out:
            f_out.writelines(f_in)
    os.remove(mtx_path)
    print(f"  MTX: {kumar_counts_t.shape[0]} genes x {kumar_counts_t.shape[1]} cells, {time.time()-t0:.1f}s")

    # Write barcodes
    bc_path = os.path.join(args.output_dir, "kumar_sc_barcodes.csv")
    pd.DataFrame({"barcode": [str(i) for i in kumar_ids]}).to_csv(bc_path, index=False)

    # Write features
    feat_path = os.path.join(args.output_dir, "kumar_sc_features.csv")
    pd.DataFrame({"gene": gene_df.index.values}).to_csv(feat_path, index=False)

    # Write metadata
    meta_path = os.path.join(args.output_dir, "kumar_sc_metadata.csv")
    pd.DataFrame({
        "barcode": [str(i) for i in kumar_ids],
        "celltype": kumar_labels,
    }).to_csv(meta_path, index=False)

    print(f"\nSaved to {args.output_dir}:")
    print(f"  {len(kumar_ids)} cells, {np.unique(kumar_labels).shape[0]} types")
    print("Done.")


if __name__ == "__main__":
    main()
