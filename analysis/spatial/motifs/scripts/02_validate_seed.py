#!/usr/bin/env python3
"""Seed stability validation for the cell-level density-weighted NMF.

Refits the same NMF (k=9 by default) on the canonical input matrix with N
independent random seeds. For each pair of seeds, computes the row-wise cosine
similarity between H matrices. Stable programs match across seeds at
cosine >= 0.95 with the reference basis; deterministic factorisations match
at cosine = 1.0.

Outputs (data/validation/s36_seed_stability/):
  - s36_pairwise_program_cosine.csv  full pairwise (n_seeds * k) x (n_seeds * k)
                                     cosine matrix
  - s36_program_stability_score.csv  per program in seed 0: best-match cosine
                                     against every other seed
  - s36_vs_canonical_basis.csv       per program in the canonical basis: best
                                     match cosine across all seeds
  - s36_seed_stability_audit.json    n_seeds, k, per-seed recon error,
                                     summary statistics

Verifies the reproducibility of the canonical 9-program decomposition without
requiring re-fitting on every downstream consumer's machine.
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
    ap.add_argument("--cells-parquet", default=None,
                    help="Default: ${project_root}/data/inputs/joint_cells.parquet")
    ap.add_argument("--reference-basis", default=None,
                    help="Default: ${project_root}/data/nmf/basis.csv")
    ap.add_argument("--out-dir", default=None,
                    help="Default: ${project_root}/data/validation/s36_seed_stability")
    ap.add_argument("--k", type=int, default=9)
    ap.add_argument("--n-seeds", type=int, default=20)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    sys.path.insert(0, os.path.join(args.project_root, "config"))
    from load_paths import load_paths, resolve  # noqa: E402

    paths = load_paths(args.project_root, args.config)
    cells_parquet = (args.cells_parquet
                     or os.path.join(resolve(paths, "data.inputs"),
                                     "joint_cells.parquet"))
    reference_basis = (args.reference_basis
                       or os.path.join(resolve(paths, "data.nmf"), "basis.csv"))
    out_dir = (args.out_dir
               or os.path.join(resolve(paths, "data.validation"),
                               "s36_seed_stability"))
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(cells_parquet):
        raise FileNotFoundError(cells_parquet)
    if not os.path.exists(reference_basis):
        raise FileNotFoundError(reference_basis)

    t0 = time.time()
    cells = pd.read_parquet(cells_parquet)
    feat_cols = [f"nonepi_neigh_frac_{l}" for l in NONEPI_21]
    missing = [c for c in feat_cols if c not in cells.columns]
    if missing:
        raise RuntimeError(f"Missing feature cols: {missing}")
    X = np.clip(cells[feat_cols].to_numpy(dtype=np.float64), 0, None)
    d = cells["local_density"].to_numpy(dtype=np.float64)
    d_log = np.log10(d + 1e-10)
    d_rank = np.argsort(np.argsort(d_log)) / max(len(d_log) - 1, 1)
    d_rank = np.clip(d_rank, 1e-3, 1.0)
    ec = cells["epi_context"].to_numpy(dtype=np.float64)
    nonepi_share = np.clip(1.0 - ec, 1e-3, 1.0)
    sqw = np.sqrt(d_rank * nonepi_share).astype(np.float32)
    X_w = X * sqw[:, None]

    rng = np.random.RandomState(0)
    seeds = rng.randint(1, 2**31 - 1, size=args.n_seeds)

    Hs = []
    diag = []
    for i, s in enumerate(seeds):
        try:
            nmf = NMF(n_components=args.k, init="nndsvda",
                      max_iter=400, random_state=int(s), tol=1e-5)
            nmf.fit(X_w)
            H = nmf.components_
            H_norm = H / (H.sum(axis=1, keepdims=True) + 1e-12)
            Hs.append(H_norm)
            diag.append({"seed_idx": i, "seed": int(s),
                         "recon_err": float(nmf.reconstruction_err_)})
        except Exception as e:
            diag.append({"seed_idx": i, "seed": int(s), "error": str(e)})

    if not Hs:
        raise SystemExit("All NMF inits failed")

    n_seeds = len(Hs)
    H_stack = np.vstack(Hs)
    norms = np.linalg.norm(H_stack, axis=1, keepdims=True) + 1e-12
    H_n = H_stack / norms
    sim_full = H_n @ H_n.T
    pd.DataFrame(
        sim_full,
        index=[f"s{i}_P{p}" for i in range(n_seeds) for p in range(args.k)],
        columns=[f"s{i}_P{p}" for i in range(n_seeds) for p in range(args.k)],
    ).to_csv(os.path.join(out_dir, "s36_pairwise_program_cosine.csv"))

    # Per-program stability vs seed 0
    ref = Hs[0]
    ref_norm = ref / (np.linalg.norm(ref, axis=1, keepdims=True) + 1e-12)
    program_scores = []
    for p in range(args.k):
        ref_vec = ref_norm[p]
        best_matches = []
        for i in range(1, n_seeds):
            Hi = Hs[i]
            Hi_n = Hi / (np.linalg.norm(Hi, axis=1, keepdims=True) + 1e-12)
            best_matches.append(float((Hi_n @ ref_vec).max()))
        if best_matches:
            program_scores.append({
                "program": f"P{p}",
                "ref_top_l1p5": NONEPI_21[int(np.argmax(ref[p]))],
                "n_compared_seeds": len(best_matches),
                "median_best_cos": float(np.median(best_matches)),
                "min_best_cos": float(np.min(best_matches)),
                "max_best_cos": float(np.max(best_matches)),
                "frac_above_0.95": float(np.mean(
                    [c >= 0.95 for c in best_matches])),
                "frac_above_0.90": float(np.mean(
                    [c >= 0.90 for c in best_matches])),
            })
    pd.DataFrame(program_scores).to_csv(
        os.path.join(out_dir, "s36_program_stability_score.csv"), index=False)

    # vs canonical basis
    canonical_compare = []
    canonical_basis_df = pd.read_csv(reference_basis)
    canonical_H = canonical_basis_df[NONEPI_21].to_numpy(dtype=np.float64)
    canonical_n = canonical_H / (
        np.linalg.norm(canonical_H, axis=1, keepdims=True) + 1e-12)
    for p in range(canonical_basis_df.shape[0]):
        matches = []
        for i in range(n_seeds):
            Hi = Hs[i]
            Hi_n = Hi / (np.linalg.norm(Hi, axis=1, keepdims=True) + 1e-12)
            matches.append(float((Hi_n @ canonical_n[p]).max()))
        canonical_compare.append({
            "canonical_program": canonical_basis_df["program"].iloc[p],
            "top_l1p5": NONEPI_21[int(np.argmax(canonical_H[p]))],
            "median_best_cos_across_seeds": float(np.median(matches)),
            "min_best_cos": float(np.min(matches)),
            "frac_seeds_match_at_0.95": float(np.mean(
                [c >= 0.95 for c in matches])),
            "frac_seeds_match_at_0.90": float(np.mean(
                [c >= 0.90 for c in matches])),
        })
    pd.DataFrame(canonical_compare).to_csv(
        os.path.join(out_dir, "s36_vs_canonical_basis.csv"), index=False)

    audit = {
        "n_seeds": int(n_seeds),
        "k": int(args.k),
        "seed_diagnostics": diag,
        "stability_summary": program_scores,
        "vs_canonical_basis": canonical_compare,
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(out_dir, "s36_seed_stability_audit.json"), "w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    print(json.dumps({"n_seeds": n_seeds, "k": int(args.k),
                      "min_match_above_0.95": float(np.min(
                          [r["frac_seeds_match_at_0.95"] for r in canonical_compare])),
                      "elapsed_sec": round(time.time() - t0, 1)}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
