#!/usr/bin/env python3
"""Re-ingest UCI604 per-position count bundles.

Legacy Exp 03 produced empty bundles for UCI604 because its raw cell_id is
int32 while other patients use string hashes. The legacy script cast
cells["cell_id"] to str but left transcripts["cell_id"] as int, so the
`tx.cell_id.isin(cell_to_idx)` mask returned all False → zero transcripts
mapped.

This script is the legacy 03a_build_boundary_counts.py with one fix:
both sides of the match are coerced to str. Outputs land in
`01_qc/outputs/xenium_counts_v2/{xenium_id}/` to avoid mutating `_legacy/`.
"""
import argparse, gzip, os, shutil, warnings
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def _read_parquet_or_csv(base, stem):
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


def process_one(sample_path, xenium_id, patient_id, out_root, qv_cutoff=20):
    out = os.path.join(out_root, xenium_id)
    os.makedirs(out, exist_ok=True)

    print(f"[{xenium_id}] loading from {sample_path}", flush=True)
    cells = _read_parquet_or_csv(sample_path, "cells")
    tx = _read_parquet_or_csv(sample_path, "transcripts")

    # Dtype fix: both cell_id sides as str.
    cells["cell_id"] = cells["cell_id"].astype(str)
    tx["cell_id"] = tx["cell_id"].astype(str)

    n0 = len(tx)
    tx = tx[
        (tx["qv"] >= qv_cutoff)
        & (~tx["feature_name"].str.startswith("BLANK_"))
        & (~tx["feature_name"].str.startswith("NegControl"))
    ]
    print(f"  transcripts after QV{qv_cutoff}+ctrl: {len(tx)}/{n0} "
          f"({100*len(tx)/n0:.1f}%)", flush=True)

    genes = np.sort(tx["feature_name"].unique())
    gene_to_idx = {g: i for i, g in enumerate(genes)}
    cell_ids = cells["cell_id"].to_numpy()
    cell_to_idx = {c: i for i, c in enumerate(cell_ids)}
    n_cells, n_genes = len(cell_ids), len(genes)
    print(f"  cells={n_cells}, genes={n_genes}", flush=True)

    tx_nuc = tx[tx["overlaps_nucleus"] == 1]
    tx_cyt = tx[tx["overlaps_nucleus"] == 0]
    print(f"  nuclear={len(tx_nuc)}, cytoplasmic={len(tx_cyt)}", flush=True)

    X_whole = _build_matrix(tx, cell_to_idx, gene_to_idx)
    X_nuc = _build_matrix(tx_nuc, cell_to_idx, gene_to_idx)
    X_cyt = _build_matrix(tx_cyt, cell_to_idx, gene_to_idx)

    def _nc(X): return np.asarray(X.sum(axis=1)).ravel().astype(np.int64)
    def _nf(X): return np.asarray((X > 0).sum(axis=1)).ravel().astype(np.int64)

    obs = pd.DataFrame({
        "cell_id": [f"{xenium_id}__{c}" for c in cell_ids],
        "raw_cell_id": cell_ids,
        "xenium_id": xenium_id,
        "patient_id": patient_id,
    })
    for nm, X in [("Whole", X_whole), ("Nuclear", X_nuc), ("Cytoplasmic", X_cyt)]:
        nc = _nc(X); nf = _nf(X)
        obs[f"nCount_{nm}"] = nc
        obs[f"nFeature_{nm}"] = nf
        obs[f"tx_per_gene_{nm}"] = np.where(nf > 0, nc / np.maximum(nf, 1), 0.0)

    nc_w = obs["nCount_Whole"].to_numpy().clip(min=1)
    obs["nuclear_fraction"] = obs["nCount_Nuclear"] / nc_w
    obs["cytoplasmic_fraction"] = obs["nCount_Cytoplasmic"] / nc_w

    for col in ("x_centroid", "y_centroid", "cell_area", "nucleus_area"):
        if col in cells.columns:
            obs[col] = cells[col].to_numpy()
    if "nucleus_area" in obs.columns and "cell_area" in obs.columns:
        obs["nucleus_ratio"] = obs["nucleus_area"] / obs["cell_area"].clip(lower=1)

    obs["qc_pass_nuclear"] = (obs["nFeature_Nuclear"] > 5) & (obs["tx_per_gene_Nuclear"] > 2)
    obs["qc_pass_whole"]   = (obs["nFeature_Whole"]   > 5) & (obs["tx_per_gene_Whole"]   > 1)

    print(f"  writing MTX bundle → {out}", flush=True)
    _write_mtx_gz(X_whole, os.path.join(out, "counts_whole.mtx.gz"))
    _write_mtx_gz(X_nuc,   os.path.join(out, "counts_nuclear.mtx.gz"))
    _write_mtx_gz(X_cyt,   os.path.join(out, "counts_cytoplasmic.mtx.gz"))

    with open(os.path.join(out, "genes.tsv"), "w") as f:
        for g in genes: f.write(g + "\n")
    with open(os.path.join(out, "cells.tsv"), "w") as f:
        for c in obs["cell_id"]: f.write(c + "\n")
    obs.to_csv(os.path.join(out, "obs.csv"), index=False)

    print(f"[{xenium_id}] pass_whole={int(obs['qc_pass_whole'].sum())}/{n_cells} "
          f"({100*obs['qc_pass_whole'].mean():.1f}%)", flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True)
    p.add_argument("--patient-filter", default="UCI604",
                   help="only process rows where patient_id == this value")
    p.add_argument("--out-root", required=True)
    p.add_argument("--qv-cutoff", type=int, default=20)
    args = p.parse_args()

    os.makedirs(args.out_root, exist_ok=True)
    man = pd.read_csv(args.manifest, sep="\t")
    sub = man[man.iloc[:, 1] == args.patient_filter]
    print(f"manifest: {len(sub)} rows matching patient_id={args.patient_filter}",
          flush=True)

    for _, row in sub.iterrows():
        xid = row.iloc[0]; pid = row.iloc[1]; txp = row.iloc[2]
        sample_path = os.path.dirname(txp)
        try:
            process_one(sample_path, xid, pid, args.out_root, args.qv_cutoff)
        except Exception as e:
            print(f"[{xid}] FAILED: {e}", flush=True)


if __name__ == "__main__":
    main()
