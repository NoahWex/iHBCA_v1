#!/usr/bin/env python3
"""Leave-one-patient-out validation of the cell-level density-weighted NMF.

For each of the 4 patients, refit a fresh NMF (k=9) on the cells from the other
3 patients, then project the held-out patient's cells through the refit H
matrix and assign argmax programs. Compare to the assignments produced by the
full-cohort H. Two summaries:

  1. H-matrix concordance: per full-cohort program, the best-matching held-out
     program by cosine similarity. Median across folds reflects vocabulary
     reproducibility under patient holdout.
  2. Cell-level retention: among held-out cells, the fraction whose argmax
     under the refit-and-mapped H equals their full-cohort argmax.

Outputs (data/validation/s38_holdout_patient/):
  - s38_holdout_per_cell.parquet      sampled per-cell agreement record
                                       (up to 5000 cells per fold)
  - s38_holdout_H_concordance.csv     per full-program x fold best-match cosine
  - s38_holdout_summary.csv           per fold: train/test counts, recon err,
                                       median H cosine, retention frac
  - s38_audit.json                    n_folds, patient list, aggregate stats
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from sklearn.decomposition import NMF


NONEPI_21 = [
    "Fb", "Fb_Activated", "Fb_SFRP4", "Adipo",
    "Plas", "T-NK", "cDC1", "CD4T", "Neu", "Mast", "Treg",
    "pDC", "Mac", "cDC2", "cDC", "CD8T", "B",
    "EC", "PV", "LEC", "Vas-cap",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True)
    ap.add_argument("--config", default=None)
    ap.add_argument("--cells-parquet", default=None)
    ap.add_argument("--basis", default=None)
    ap.add_argument("--out-dir", default=None)
    ap.add_argument("--k", type=int, default=9)
    ap.add_argument("--n-init", type=int, default=5)
    return ap.parse_args()


def build_weighted_input(df):
    feat_cols = [f"nonepi_neigh_frac_{l}" for l in NONEPI_21]
    X = np.clip(df[feat_cols].to_numpy(dtype=np.float64), 0, None)
    d = df["local_density"].to_numpy(dtype=np.float64)
    d_log = np.log10(d + 1e-10)
    d_rank = np.argsort(np.argsort(d_log)) / max(len(d_log) - 1, 1)
    d_rank = np.clip(d_rank, 1e-3, 1.0)
    ec = df["epi_context"].to_numpy(dtype=np.float64)
    nonepi_share = np.clip(1.0 - ec, 1e-3, 1.0)
    sqw = np.sqrt(d_rank * nonepi_share).astype(np.float32)
    return X, X * sqw[:, None], sqw


def project_batched(X, H):
    HHT = H @ H.T
    inv = np.linalg.inv(HHT + 1e-8 * np.eye(H.shape[0]))
    return np.clip(X @ H.T @ inv, 0, None).astype(np.float32)


def main() -> int:
    args = parse_args()
    sys.path.insert(0, os.path.join(args.project_root, "config"))
    from load_paths import load_paths, resolve  # noqa: E402

    paths = load_paths(args.project_root, args.config)
    cells_parquet = (args.cells_parquet
                     or os.path.join(resolve(paths, "data.inputs"),
                                     "joint_cells.parquet"))
    basis_path = (args.basis
                  or os.path.join(resolve(paths, "data.nmf"), "basis.csv"))
    out_dir = (args.out_dir
               or os.path.join(resolve(paths, "data.validation"),
                               "s38_holdout_patient"))
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    cells = pd.read_parquet(cells_parquet)
    patients = sorted(cells["patient_id"].unique())

    basis = pd.read_csv(basis_path)
    H_full = basis[NONEPI_21].to_numpy(dtype=np.float64)
    H_full_norm = H_full / (H_full.sum(axis=1, keepdims=True) + 1e-12)
    k = H_full.shape[0]

    X_all, _, _ = build_weighted_input(cells)
    U_real = project_batched(X_all, H_full)
    U_real_norm = U_real / (U_real.sum(axis=1, keepdims=True) + 1e-12)
    real_argmax = np.argmax(U_real_norm, axis=1)

    fold_rows = []
    H_concordance_rows = []
    per_cell_rows = []

    for held_out in patients:
        train_mask = cells["patient_id"] != held_out
        test_mask = cells["patient_id"] == held_out
        cells_train = cells[train_mask]
        cells_test = cells[test_mask]
        n_train, n_test = len(cells_train), len(cells_test)
        if n_train < args.k * 2 or n_test == 0:
            continue

        _, X_train_w, _ = build_weighted_input(cells_train)
        best_err = np.inf
        best_H = None
        rng = np.random.RandomState(42)
        for _ in range(args.n_init):
            s = int(rng.randint(0, 2**31 - 1))
            try:
                nmf = NMF(n_components=args.k, init="nndsvda",
                          max_iter=400, random_state=s, tol=1e-5)
                nmf.fit(X_train_w)
                if nmf.reconstruction_err_ < best_err:
                    best_err = float(nmf.reconstruction_err_)
                    best_H = nmf.components_
            except Exception:
                continue
        if best_H is None:
            continue
        H_train = best_H
        H_train_norm = H_train / (H_train.sum(axis=1, keepdims=True) + 1e-12)

        Hf_n = H_full_norm / (np.linalg.norm(H_full_norm, axis=1,
                                             keepdims=True) + 1e-12)
        Ht_n = H_train_norm / (np.linalg.norm(H_train_norm, axis=1,
                                              keepdims=True) + 1e-12)
        sim = Ht_n @ Hf_n.T

        for p in range(k):
            row = sim[:, p]
            best_train_idx = int(row.argmax())
            H_concordance_rows.append({
                "held_out": held_out,
                "full_program": f"P{p}",
                "top_l1p5_full": NONEPI_21[int(np.argmax(H_full[p]))],
                "best_train_program": f"P{best_train_idx}",
                "best_train_top_l1p5":
                    NONEPI_21[int(np.argmax(H_train[best_train_idx]))],
                "best_cos": float(row.max()),
            })

        X_test, _, _ = build_weighted_input(cells_test)
        U_test = project_batched(X_test, H_train)
        U_test_norm = U_test / (U_test.sum(axis=1, keepdims=True) + 1e-12)
        test_argmax_train = np.argmax(U_test_norm, axis=1)
        train_to_full = np.zeros(k, dtype=np.int32)
        for q in range(k):
            train_to_full[q] = int(np.argmax(sim[q]))
        test_argmax_mapped = train_to_full[test_argmax_train]
        test_real = real_argmax[test_mask.values]
        n_match = int((test_argmax_mapped == test_real).sum())

        fold_rows.append({
            "held_out": held_out,
            "n_train_cells": n_train,
            "n_test_cells": n_test,
            "H_recon_err": best_err,
            "frac_retained": n_match / max(n_test, 1),
            "median_H_cos": float(np.median(
                [r["best_cos"] for r in H_concordance_rows
                 if r["held_out"] == held_out])),
            "n_programs_match_at_0.95": int(np.sum(
                [r["best_cos"] >= 0.95 for r in H_concordance_rows
                 if r["held_out"] == held_out])),
        })

        idx_sample = np.random.RandomState(0).choice(
            len(test_real), size=min(5000, len(test_real)), replace=False)
        for ix in idx_sample:
            per_cell_rows.append({
                "held_out": held_out,
                "real_program": f"P{int(test_real[ix])}",
                "holdout_program": f"P{int(test_argmax_mapped[ix])}",
                "agree": bool(test_argmax_mapped[ix] == test_real[ix]),
            })

    fold_df = pd.DataFrame(fold_rows)
    fold_df.to_csv(os.path.join(out_dir, "s38_holdout_summary.csv"), index=False)
    Hcd = pd.DataFrame(H_concordance_rows)
    Hcd.to_csv(os.path.join(out_dir, "s38_holdout_H_concordance.csv"), index=False)
    pd.DataFrame(per_cell_rows).to_parquet(
        os.path.join(out_dir, "s38_holdout_per_cell.parquet"), index=False)

    audit = {
        "n_folds": int(len(fold_rows)),
        "patients_tested": list(fold_df["held_out"]),
        "median_retention_frac": float(fold_df["frac_retained"].median()),
        "median_H_cosine_per_program": float(Hcd["best_cos"].median()),
        "n_perfect_program_match_at_0.95": int((Hcd["best_cos"] >= 0.95).sum()),
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(out_dir, "s38_audit.json"), "w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    print(json.dumps(audit, default=str), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
