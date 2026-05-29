#!/usr/bin/env python3
"""Step 03a — Per-position, per-boundary count matrices.

For one Xenium sample directory, reads transcripts.parquet + cells.parquet and
splits transcripts by `overlaps_nucleus` to build three count matrices:

  whole        = all QV-passing transcripts assigned to a cell
  nuclear      = overlaps_nucleus == 1  (inside nuclear polygon)
  cytoplasmic  = overlaps_nucleus == 0  (expansion zone)

Per-boundary QC metrics are computed; a soft `qc_pass` flag applies the plan's
nuclear thresholds (nFeature_Nuclear > 5 AND transcripts_per_gene_Nuclear > 2).

Outputs (to --out-dir/{xenium_id}/):
  counts_whole.mtx.gz        cells × genes, int32
  counts_nuclear.mtx.gz
  counts_cytoplasmic.mtx.gz
  genes.tsv                   one gene per line (shared across matrices)
  cells.tsv                   cell_id per row, prefixed with {xenium_id}__
  obs.csv                     per-cell metadata + QC metrics + qc_pass

(reused logic for overlaps_nucleus split and QC metric definitions).
"""
import argparse
import gzip
import os
import shutil
import warnings

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def _read_parquet_or_csv(base, stem):
    """Read {stem}.parquet if present else {stem}.csv.gz. Decode bytes cols."""
    pq = os.path.join(base, f"{stem}.parquet")
    csv = os.path.join(base, f"{stem}.csv.gz")
    if os.path.isfile(pq):
        try:
            df = pd.read_parquet(pq)
        except Exception as e:
            if os.path.isfile(csv):
                warnings.warn(f"corrupt {pq}, fallback csv: {e}")
                df = pd.read_csv(csv)
            else:
                raise
    elif os.path.isfile(csv):
        df = pd.read_csv(csv)
    else:
        raise FileNotFoundError(f"no {stem}.parquet or .csv.gz in {base}")
    for col in df.columns:
        if df[col].dtype == object and len(df) and isinstance(df[col].iloc[0], bytes):
            df[col] = df[col].str.decode("utf-8")
    return df


def _build_matrix(tx, cell_to_idx, gene_to_idx):
    """transcripts df → sparse (cells × genes) int32."""
    mask = tx["feature_name"].isin(gene_to_idx) & tx["cell_id"].isin(cell_to_idx)
    t = tx[mask]
    rows = t["cell_id"].map(cell_to_idx).to_numpy()
    cols = t["feature_name"].map(gene_to_idx).to_numpy()
    data = np.ones(len(rows), dtype=np.int32)
    return sp.coo_matrix(
        (data, (rows, cols)),
        shape=(len(cell_to_idx), len(gene_to_idx)),
        dtype=np.int32,
    ).tocsr()


def _write_mtx_gz(X, path):
    tmp = f"/tmp/_counts_{os.getpid()}.mtx"
    scipy.io.mmwrite(tmp, X, field="integer")
    with open(tmp, "rb") as src, gzip.open(path, "wb") as dst:
        shutil.copyfileobj(src, dst)
    os.remove(tmp)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-path", required=True, help="Xenium sample directory")
    p.add_argument("--xenium-id", required=True, help="e.g. Pat1_P1_xenium_1")
    p.add_argument("--patient-id", required=True)
    p.add_argument("--out-dir", required=True, help="parent; subdir per xenium_id")
    p.add_argument("--qv-cutoff", type=int, default=20)
    args = p.parse_args()

    out = os.path.join(args.out_dir, args.xenium_id)
    os.makedirs(out, exist_ok=True)

    print(f"[{args.xenium_id}] loading cells + transcripts from {args.sample_path}")
    cells = _read_parquet_or_csv(args.sample_path, "cells")
    tx = _read_parquet_or_csv(args.sample_path, "transcripts")

    # QV + control filter
    n0 = len(tx)
    tx = tx[
        (tx["qv"] >= args.qv_cutoff)
        & (~tx["feature_name"].str.startswith("BLANK_"))
        & (~tx["feature_name"].str.startswith("NegControl"))
    ]
    print(f"  transcripts after QV{args.qv_cutoff} + control filter: "
          f"{len(tx)}/{n0} ({100*len(tx)/n0:.1f}%)")

    # Gene + cell indexes
    genes = np.sort(tx["feature_name"].unique())
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    cell_ids = cells["cell_id"].astype(str).to_numpy()
    cell_to_idx = {c: i for i, c in enumerate(cell_ids)}
    n_cells, n_genes = len(cell_ids), len(genes)
    print(f"  cells={n_cells}, genes={n_genes}")

    # Split by overlaps_nucleus
    tx_nuc = tx[tx["overlaps_nucleus"] == 1]
    tx_cyt = tx[tx["overlaps_nucleus"] == 0]
    print(f"  nuclear={len(tx_nuc)}, cytoplasmic={len(tx_cyt)}")

    X_whole = _build_matrix(tx, cell_to_idx, gene_to_idx)
    X_nuc = _build_matrix(tx_nuc, cell_to_idx, gene_to_idx)
    X_cyt = _build_matrix(tx_cyt, cell_to_idx, gene_to_idx)

    # Per-boundary QC metrics
    def _nc(X):
        return np.asarray(X.sum(axis=1)).ravel().astype(np.int64)

    def _nf(X):
        return np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int64)

    obs = pd.DataFrame({
        "cell_id": [f"{args.xenium_id}__{c}" for c in cell_ids],
        "raw_cell_id": cell_ids,
        "xenium_id": args.xenium_id,
        "patient_id": args.patient_id,
    })
    for nm, X in [("Whole", X_whole), ("Nuclear", X_nuc), ("Cytoplasmic", X_cyt)]:
        nc = _nc(X)
        nf = _nf(X)
        obs[f"nCount_{nm}"] = nc
        obs[f"nFeature_{nm}"] = nf
        obs[f"tx_per_gene_{nm}"] = np.where(nf > 0, nc / np.maximum(nf, 1), 0.0)

    # Fractions
    nc_w = obs["nCount_Whole"].to_numpy().clip(min=1)
    obs["nuclear_fraction"] = obs["nCount_Nuclear"] / nc_w
    obs["cytoplasmic_fraction"] = obs["nCount_Cytoplasmic"] / nc_w

    # Centroids + cell metadata from cells.parquet
    for col in ("x_centroid", "y_centroid", "cell_area", "nucleus_area"):
        if col in cells.columns:
            obs[col] = cells[col].to_numpy()
    if "nucleus_area" in obs.columns and "cell_area" in obs.columns:
        obs["nucleus_ratio"] = obs["nucleus_area"] / obs["cell_area"].clip(lower=1)

    # Plan QC: nuclear gate. Soft flag — cells retained, downstream decides.
    obs["qc_pass_nuclear"] = (
        (obs["nFeature_Nuclear"] > 5) & (obs["tx_per_gene_Nuclear"] > 2)
    )
    obs["qc_pass_whole"] = (
        (obs["nFeature_Whole"] > 5) & (obs["tx_per_gene_Whole"] > 2)
    )

    # Write MTX triplet
    print(f"  writing MTX bundle → {out}")
    _write_mtx_gz(X_whole, os.path.join(out, "counts_whole.mtx.gz"))
    _write_mtx_gz(X_nuc, os.path.join(out, "counts_nuclear.mtx.gz"))
    _write_mtx_gz(X_cyt, os.path.join(out, "counts_cytoplasmic.mtx.gz"))

    with open(os.path.join(out, "genes.tsv"), "w") as f:
        for g in genes:
            f.write(g + "\n")
    with open(os.path.join(out, "cells.tsv"), "w") as f:
        for c in obs["cell_id"]:
            f.write(c + "\n")
    obs.to_csv(os.path.join(out, "obs.csv"), index=False)

    # Summary
    n_pass_nuc = int(obs["qc_pass_nuclear"].sum())
    n_pass_wh = int(obs["qc_pass_whole"].sum())
    print(f"\n[{args.xenium_id}] summary")
    print(f"  nuclear QC pass: {n_pass_nuc}/{n_cells} ({100*n_pass_nuc/n_cells:.1f}%)")
    print(f"  whole   QC pass: {n_pass_wh}/{n_cells} ({100*n_pass_wh/n_cells:.1f}%)")
    print(f"  median nuclear_fraction: {obs['nuclear_fraction'].median():.3f}")
    print("Done.")


if __name__ == "__main__":
    main()
