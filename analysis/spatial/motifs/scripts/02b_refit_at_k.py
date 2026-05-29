#!/usr/bin/env python3
"""Refit NMF at alternate K's and Hungarian-match to the canonical K=9 basis.

Companion to 01_cell_nmf.py. Loads the persisted NMF input substrate
(`nmf_input_features.parquet` — written by 01_cell_nmf.py), refits NMF at
each K in --k-list with the same hyperparameters as canonical (n_init=15,
init='nndsvda', max_iter=400, tol=1e-5), and reports per-K Hungarian-matched
cosine similarity to the canonical basis.csv (K=9).

Outputs (in --out-dir):
  - basis_K{N}.csv (one per K, same schema as canonical basis.csv)
  - rank_neighborhood_cosine.csv (long-format: k_refit, refit_motif,
    canonical_motif, cosine, matched (bool))
  - rank_neighborhood_audit.json

Interpretation:
  - K < canonical: matched programs at high cosine indicate "stable" canonical
    motifs; unmatched canonical motifs indicate motifs that were absorbed
    into a shared lower-K program. Examine the row in the lower-K basis
    whose Hungarian-match has the lowest cosine to identify what collapsed.
  - K > canonical: matched canonical motifs at high cosine indicate stable
    composition; unmatched refit-K programs are candidate "splits" of an
    existing canonical motif. Examine cosine of each canonical motif to
    multiple refit-K programs to identify the split.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import NMF


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--features", required=True,
                    help="nmf_input_features.parquet from 01_cell_nmf.py")
    ap.add_argument("--canonical-basis", required=True,
                    help="basis.csv from 01_cell_nmf.py (canonical K=9)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--k-list", default="6,11",
                    help="Comma-separated list of K values to refit "
                    "(default: 6,11)")
    ap.add_argument("--n-init", type=int, default=15)
    ap.add_argument("--random-state", type=int, default=0)
    return ap.parse_args()


def cosine_matrix(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Pairwise cosine similarity between rows of A and rows of B.

    A: (m, d), B: (n, d). Returns (m, n) matrix.
    """
    A_norm = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    B_norm = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)
    return A_norm @ B_norm.T


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    print(f"[load] features", flush=True)
    feat = pd.read_parquet(args.features)
    xw_cols = [c for c in feat.columns if c.startswith("Xw_")]
    if not xw_cols:
        sys.exit("ERROR: features parquet missing Xw_* columns")
    X_w = feat[xw_cols].to_numpy(dtype=np.float64)
    print(f"  X_w shape: {X_w.shape}", flush=True)

    print(f"[load] canonical basis", flush=True)
    canonical = pd.read_csv(args.canonical_basis)
    feat_cols = [c for c in canonical.columns if c != "program"]
    H_canonical = canonical[feat_cols].to_numpy(dtype=np.float64)
    canon_programs = canonical["program"].tolist()
    print(f"  canonical basis: {H_canonical.shape}, programs={canon_programs}",
          flush=True)

    if len(feat_cols) != len(xw_cols):
        sys.exit(f"ERROR: feature dim mismatch: features has {len(xw_cols)} "
                 f"cols, canonical basis has {len(feat_cols)}")

    k_list = [int(k.strip()) for k in args.k_list.split(",")]
    rng = np.random.RandomState(args.random_state)
    long_rows = []
    audit = {"k_list": k_list, "n_init": args.n_init,
             "n_cells": int(X_w.shape[0]), "k_results": {}}

    for k in k_list:
        print(f"\n=== K = {k} ===", flush=True)
        best_err = np.inf
        best_model = None
        for i in range(args.n_init):
            s = int(rng.randint(0, 2**31 - 1))
            try:
                nmf = NMF(n_components=k, init="nndsvda",
                          max_iter=400, random_state=s, tol=1e-5)
                nmf.fit(X_w)
                err = nmf.reconstruction_err_
            except Exception as e:
                print(f"  init {i} failed: {e}", flush=True)
                continue
            if err < best_err:
                best_err = err
                best_model = nmf
        if best_model is None:
            print(f"  K={k}: all inits failed", flush=True)
            continue
        H_k = best_model.components_
        H_k_norm = H_k / (H_k.sum(axis=1, keepdims=True) + 1e-12)
        print(f"  recon_err={best_err:.4f}; basis shape {H_k_norm.shape}",
              flush=True)

        basis_df = pd.DataFrame(
            H_k_norm, columns=feat_cols
        ).assign(program=[f"P{i}" for i in range(k)])
        basis_df = basis_df[["program"] + feat_cols]
        basis_path = os.path.join(args.out_dir, f"basis_K{k}.csv")
        basis_df.to_csv(basis_path, index=False)
        print(f"  [write] {basis_path}", flush=True)

        # Pairwise cosines: rows = K-refit programs, cols = canonical K=9
        cos_mat = cosine_matrix(H_k_norm, H_canonical)
        # Hungarian match: max cosine ↔ negate for assignment minimizer
        n_pair = min(k, len(canon_programs))
        row_ind, col_ind = linear_sum_assignment(-cos_mat)
        matched_pairs = set(zip(row_ind, col_ind))
        for ri in range(k):
            for ci in range(len(canon_programs)):
                long_rows.append({
                    "k_refit": k,
                    "refit_motif": f"P{ri}",
                    "canonical_motif": canon_programs[ci],
                    "cosine": float(cos_mat[ri, ci]),
                    "matched": (ri, ci) in matched_pairs,
                })

        matched_cos = float(np.mean([cos_mat[ri, ci] for ri, ci
                                       in zip(row_ind, col_ind)]))
        unmatched_refit = [f"P{ri}" for ri in range(k) if ri not in row_ind]
        unmatched_canon = [canon_programs[ci] for ci in range(len(canon_programs))
                            if ci not in col_ind]
        audit["k_results"][k] = {
            "recon_err":            float(best_err),
            "mean_matched_cosine":  matched_cos,
            "n_matched_pairs":      int(n_pair),
            "unmatched_refit":      unmatched_refit,
            "unmatched_canonical":  unmatched_canon,
        }
        print(f"  K={k} mean_matched_cosine={matched_cos:.3f}; "
              f"unmatched_refit={unmatched_refit}; "
              f"unmatched_canonical={unmatched_canon}", flush=True)

    long_df = pd.DataFrame(long_rows)
    long_path = os.path.join(args.out_dir, "rank_neighborhood_cosine.csv")
    long_df.to_csv(long_path, index=False)
    print(f"\n[write] {long_path} ({len(long_df):,} rows)", flush=True)

    audit["elapsed_sec"] = float(time.time() - t0)
    audit_path = os.path.join(args.out_dir, "rank_neighborhood_audit.json")
    with open(audit_path, "w") as f:
        json.dump(audit, f, indent=2)
    print(f"[write] {audit_path}", flush=True)
    print(f"[done] {audit['elapsed_sec']:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
