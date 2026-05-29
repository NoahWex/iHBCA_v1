"""03a_pseudobulk.py — aggregate counts to (cluster × library) pseudobulk per compartment.

FLEX is single-study (Spatial HBCA), so there is no study covariate. FLEX
`sample_id` is the library identifier; `patient_id` provides the dupCor block
downstream. Hierarchy is patient → library:
    cluster        = one-vs-rest target (the limma `condition` factor)
    library_id     = design covariate in `~ library_id + condition`
    patient_id     = dupCor block (biological replicate) when n_patients >= 3

Outputs:
    {Compartment}_pseudobulk_counts.csv  rows=sample_id (cluster__library_id), cols=genes
    {Compartment}_pseudobulk_meta.csv    sample_id, cluster, library_id, patient_id, n_cells
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import gzip
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("pseudobulk")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compartment", required=True, choices=["Epithelial", "Immune", "Stromal"])
    p.add_argument("--bundle-dir", required=True, type=Path)
    p.add_argument("--clusters-csv", required=True, type=Path)
    p.add_argument("--resolution", required=True,
                   help="Resolution suffix; clusters CSV must have column leiden_<res>")
    p.add_argument("--min-cells", type=int, default=10,
                   help="Min cells per (cluster × library) sample (default 10)")
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def load_mtx_gz(mtx_path: Path) -> sp.csr_matrix:
    log.info("Loading mtx.gz: %s", mtx_path)
    with gzip.open(mtx_path, "rb") as fh:
        m = sio.mmread(fh)
    if not sp.issparse(m):
        m = sp.csr_matrix(m)
    return m.tocsr()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    cells_tsv = args.bundle_dir / "cells.tsv"
    genes_tsv = args.bundle_dir / "genes.tsv"
    mtx = args.bundle_dir / "counts.mtx.gz"
    for f in (cells_tsv, genes_tsv, mtx):
        if not f.exists():
            sys.exit(f"ERROR: missing input: {f}")

    cells = pd.read_csv(cells_tsv, sep="\t", header=None, names=["cell_id"])
    genes = pd.read_csv(genes_tsv, sep="\t", header=None, names=["gene"])
    log.info("Cells: %d, Genes: %d", len(cells), len(genes))

    X = load_mtx_gz(mtx)
    if X.shape == (len(cells), len(genes)):
        log.info("  mtx oriented cells × genes")
    elif X.shape == (len(genes), len(cells)):
        log.info("  mtx oriented genes × cells; transposing to cells × genes")
        X = X.T.tocsr()
    else:
        sys.exit(f"ERROR: mtx shape {X.shape} does not match cells={len(cells)} genes={len(genes)}")

    clusters = pd.read_csv(args.clusters_csv)
    log.info("Clusters CSV: %d rows", len(clusters))
    leiden_col = f"leiden_{args.resolution}"
    if leiden_col not in clusters.columns:
        sys.exit(f"ERROR: clusters CSV missing column {leiden_col}; have {list(clusters.columns)}")
    clusters = clusters.rename(columns={leiden_col: "cluster"})
    if not {"cell_id", "cluster", "sample_id", "patient_id"}.issubset(clusters.columns):
        sys.exit(f"ERROR: clusters CSV missing required columns; have {list(clusters.columns)}")

    # Align cells to cluster CSV
    common = pd.Index(clusters["cell_id"]).intersection(pd.Index(cells["cell_id"]))
    log.info("Common cells: %d", len(common))
    if len(common) == 0:
        sys.exit("ERROR: no cells in common between bundle and clusters CSV")

    cell_to_idx = {c: i for i, c in enumerate(cells["cell_id"])}
    keep_idx = np.array([cell_to_idx[c] for c in common], dtype=np.int64)
    X = X[keep_idx]
    cluster_idx = clusters.set_index("cell_id").loc[common]

    cluster_per_cell = cluster_idx["cluster"].astype(str).values
    library_per_cell = cluster_idx["sample_id"].astype(str).values
    patient_per_cell = cluster_idx["patient_id"].astype(str).values

    log.info("Aggregating pseudobulk by (cluster × library_id)...")
    df_meta = pd.DataFrame({
        "cluster": cluster_per_cell,
        "library_id": library_per_cell,
        "patient_id": patient_per_cell,
    })
    groups = df_meta.groupby(["cluster", "library_id"], observed=True).indices

    samples = []
    count_rows = {}
    for key, idx in groups.items():
        clust, lib = key  # type: ignore[misc]
        if len(idx) < args.min_cells:
            continue
        sample_id = f"{clust}__{lib}"
        counts = np.asarray(X[idx].sum(axis=0)).flatten()
        count_rows[sample_id] = counts
        patients = np.unique(patient_per_cell[idx])
        if len(patients) > 1:
            log.warning("  library_id %s has multiple patient_id %s — using first",
                        lib, patients.tolist())
        samples.append({
            "sample_id": sample_id,
            "cluster": clust,
            "library_id": lib,
            "patient_id": patients[0],
            "n_cells": len(idx),
        })

    if not samples:
        sys.exit(f"ERROR: no samples passed min_cells={args.min_cells}")

    counts_df = pd.DataFrame(count_rows, index=genes["gene"].tolist()).T
    meta_df = pd.DataFrame(samples)
    log.info("Pseudobulk: %d samples × %d genes", len(counts_df), counts_df.shape[1])
    log.info("  clusters: %d, libraries: %d, patients: %d",
             meta_df["cluster"].nunique(),
             meta_df["library_id"].nunique(),
             meta_df["patient_id"].nunique())

    tag = f"{args.compartment}_leiden_{args.resolution}"
    counts_path = args.out_dir / f"{tag}_pseudobulk_counts.csv"
    meta_path = args.out_dir / f"{tag}_pseudobulk_meta.csv"
    counts_df.to_csv(counts_path)
    meta_df.to_csv(meta_path, index=False)
    log.info("Wrote %s", counts_path)
    log.info("Wrote %s", meta_path)


if __name__ == "__main__":
    main()
