#!/usr/bin/env python3
"""S1 — Per-patient NMF reproducibility.

For each of the 4 Xenium patients (Pat1, Pat2, UCI604, UCI220228), refit
the cell-level density-weighted KNN-composition NMF at k=9 on that patient's
cells alone. Hungarian-match the recovered basis to the canonical cohort
basis (data/nmf/basis.csv). Report per-motif cosine similarity.

If the same 9 motifs emerge across all 4 patients, the cohort-level motif
decomposition reflects shared biological structure, not patient-specific
artifacts.

Outputs:
  - per_patient_basis_cosine.csv   (4 x 9 matrix)
  - per_patient_match_audit.json   (per-patient match indices + recon_err)
  - S1_per_patient_nmf.pdf          (heatmap panel)
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.decomposition import NMF


NONEPI_21 = [
    "Fb", "Fb_Activated", "Fb_SFRP4", "Adipo",
    "Plas", "T-NK", "cDC1", "CD4T", "Neu", "Mast", "Treg",
    "pDC", "Mac", "cDC2", "cDC", "CD8T", "B",
    "EC", "PV", "LEC", "Vas-cap",
]


def fit_nmf_best(X_w, k, n_init, base_seed):
    best_err = np.inf
    best_H = None
    rng = np.random.RandomState(base_seed)
    for _ in range(n_init):
        s = int(rng.randint(0, 2**31 - 1))
        try:
            nmf = NMF(n_components=k, init="nndsvda", max_iter=400,
                      random_state=s, tol=1e-5)
            nmf.fit(X_w)
            if nmf.reconstruction_err_ < best_err:
                best_err = nmf.reconstruction_err_
                best_H = nmf.components_
        except Exception:
            continue
    return best_H, best_err


def cosine_matrix(A, B):
    """Row-wise cosine similarity between basis matrices A (kA x d) and B (kB x d)."""
    An = A / (np.linalg.norm(A, axis=1, keepdims=True) + 1e-12)
    Bn = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)
    return An @ Bn.T  # kA x kB


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joint-cells", required=True,
                    help="joint_cells.parquet")
    ap.add_argument("--canonical-basis", required=True,
                    help="data/nmf/basis.csv (9 x 21 with `program` col)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--k", type=int, default=9)
    ap.add_argument("--n-init", type=int, default=15)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("[load] joint_cells.parquet", flush=True)
    jc = pd.read_parquet(args.joint_cells)
    feat_cols = [f"nonepi_neigh_frac_{lab}" for lab in NONEPI_21]
    print(f"[load] {len(jc):,} cells", flush=True)

    print("[load] canonical basis", flush=True)
    basis_df = pd.read_csv(args.canonical_basis)
    canonical_H = basis_df[NONEPI_21].to_numpy()
    print(f"       shape {canonical_H.shape}", flush=True)

    patients = ["Pat1", "Pat2", "UCI604", "UCI220228"]
    cosine_results = []  # rows = patient, cols = canonical motif Mi
    audit = []

    for pid in patients:
        sub = jc[jc["patient_id"] == pid].copy()
        if len(sub) < 1000:
            print(f"[skip] {pid}: only {len(sub)} cells", flush=True)
            continue
        X = sub[feat_cols].to_numpy()
        w = np.sqrt(np.maximum(
            sub["local_density"].to_numpy() * (1 - sub["epi_context"].to_numpy()),
            0)).astype(np.float32)
        X_w = X * w[:, None]
        # Drop zero rows (no neighbors) so NMF doesn't choke
        keep = X_w.sum(axis=1) > 0
        X_w = X_w[keep]
        print(f"[fit] {pid}: {len(X_w):,} non-zero cells × {X_w.shape[1]}-D",
              flush=True)

        H, err = fit_nmf_best(X_w, args.k, args.n_init, args.seed)
        print(f"      recon_err={err:.4f}", flush=True)

        # Hungarian-match canonical_H rows → patient H rows (maximize cosine)
        C = cosine_matrix(canonical_H, H)  # 9 x 9
        row_ind, col_ind = linear_sum_assignment(-C)
        matched_cos = [float(C[i, j]) for i, j in zip(row_ind, col_ind)]
        cosine_results.append({"patient": pid, **{
            f"M{i}": matched_cos[i] for i in range(args.k)}})
        audit.append({"patient": pid, "n_cells": int(len(X_w)),
                      "recon_err": float(err),
                      "match_indices": [int(j) for j in col_ind],
                      "matched_cosines": matched_cos,
                      "median_cosine": float(np.median(matched_cos)),
                      "min_cosine": float(np.min(matched_cos))})

    # Write CSV + audit
    cos_df = pd.DataFrame(cosine_results)
    cos_df.to_csv(out / "per_patient_basis_cosine.csv", index=False)
    with open(out / "per_patient_match_audit.json", "w") as f:
        json.dump(audit, f, indent=2)
    print(f"[write] {out / 'per_patient_basis_cosine.csv'}", flush=True)

    # Render panel
    mat = cos_df[[f"M{i}" for i in range(args.k)]].to_numpy()
    fig, ax = plt.subplots(figsize=(3.0, 1.6), dpi=300)
    im = ax.imshow(mat, cmap="viridis", vmin=0, vmax=1, aspect="auto")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=4, color="white" if v < 0.5 else "black")
    ax.set_xticks(range(args.k))
    ax.set_xticklabels([f"M{i}" for i in range(args.k)], fontsize=6)
    ax.set_yticks(range(len(cos_df)))
    ax.set_yticklabels(cos_df["patient"].tolist(), fontsize=6)
    ax.tick_params(length=1.0, pad=1.0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("basis cosine", fontsize=6)
    cbar.ax.tick_params(labelsize=5, length=1.0, pad=1.0)
    for spine in ax.spines.values():
        spine.set_linewidth(0.3)

    pdf = out / "S1_per_patient_nmf.pdf"
    fig.savefig(pdf, bbox_inches="tight", dpi=300)
    plt.close(fig)
    print(f"[write] {pdf}", flush=True)


if __name__ == "__main__":
    main()
