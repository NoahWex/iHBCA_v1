#!/usr/bin/env python3
"""Spatial null for the cell-level density-weighted NMF.

For each Xenium SF instance, shuffle cell xy positions WITHIN the instance.
Cell identities (L1.5 labels) are preserved; only their spatial arrangement is
permuted. Re-compute the K=10 nearest-neighbor non-epithelial composition under
the shuffled geometry, then project the shuffled feature vectors through the
fixed canonical H matrix (closed-form non-negative least squares) to obtain
shuffled program loadings.

Per program, retention fraction = (cells whose argmax program is unchanged
after shuffle) / (cells with that argmax in the real data). Programs that
depend on real spatial organisation should drop substantially below the
observed argmax distribution; programs that are mostly a statement about cell-
type frequency should retain.

Outputs (data/validation/s37_xy_shuffle/):
  - s37_per_cell_real_vs_perm.parquet  per-cell real argmax vs shuffled argmax
                                       + max loadings
  - s37_program_stability.csv          per-program retention fraction and
                                       shuffled-argmax distribution
  - s37_audit.json                     cohort concordance + per-program summary
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


EPI_TYPES = {"LASP-basal", "LASP", "LHS", "BMYO-myo"}
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
    ap.add_argument("--k-nn", type=int, default=10)
    ap.add_argument("--n-perms", type=int, default=1)
    ap.add_argument("--random-state", type=int, default=42)
    return ap.parse_args()


def compute_features_from_shuffled(df_sample, k_nn, rng, nonepi_to_idx):
    n = len(df_sample)
    feat = np.zeros((n, len(NONEPI_21)), dtype=np.float32)
    labels = df_sample["l1p5_short"].to_numpy()
    inst_arr = df_sample["instance_id"].to_numpy()
    xy_orig = df_sample[["x", "y"]].to_numpy(dtype=np.float32)

    xy_shuf = xy_orig.copy()
    for inst_id in np.unique(inst_arr):
        idx = np.where(inst_arr == inst_id)[0]
        if len(idx) < 2:
            continue
        perm = rng.permutation(len(idx))
        xy_shuf[idx] = xy_orig[idx[perm]]

    for inst_id in np.unique(inst_arr):
        idx = np.where(inst_arr == inst_id)[0]
        if len(idx) < 2:
            continue
        sub_xy = xy_shuf[idx]
        tree = cKDTree(sub_xy)
        k_eff = min(k_nn + 1, len(idx))
        _, nn_idx = tree.query(sub_xy, k=k_eff)
        nn_idx = nn_idx[:, 1:]
        for row_local in range(len(idx)):
            row_global = idx[row_local]
            neigh_labels = labels[idx[nn_idx[row_local]]]
            for nl in neigh_labels:
                if nl in nonepi_to_idx:
                    feat[row_global, nonepi_to_idx[nl]] += 1

    row_sum = feat.sum(axis=1, keepdims=True)
    nz = row_sum.flatten() > 0
    feat_frac = np.zeros_like(feat)
    feat_frac[nz] = feat[nz] / row_sum[nz]
    return feat_frac


def project_batched(X, H):
    """Closed-form non-negative least squares projection onto H rows."""
    HHT = H @ H.T
    inv = np.linalg.inv(HHT + 1e-8 * np.eye(H.shape[0]))
    U = X @ H.T @ inv
    return np.clip(U, 0, None).astype(np.float32)


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
                               "s37_xy_shuffle"))
    os.makedirs(out_dir, exist_ok=True)

    t0 = time.time()
    cells = pd.read_parquet(cells_parquet)
    feat_cols = [f"nonepi_neigh_frac_{l}" for l in NONEPI_21]
    X_real = np.clip(cells[feat_cols].to_numpy(dtype=np.float64), 0, None)

    basis = pd.read_csv(basis_path)
    H = basis[NONEPI_21].to_numpy(dtype=np.float64)
    k = H.shape[0]

    U_real = project_batched(X_real, H)
    U_real_norm = U_real / (U_real.sum(axis=1, keepdims=True) + 1e-12)
    real_argmax = np.argmax(U_real_norm, axis=1)

    rng = np.random.RandomState(args.random_state)
    nonepi_to_idx = {l: i for i, l in enumerate(NONEPI_21)}

    cells_with_meta = cells.copy()
    cells_with_meta["real_program_idx"] = real_argmax

    perm_argmax_all = []
    perm_loadings_all = []
    for _ in range(args.n_perms):
        all_feat = np.zeros_like(X_real, dtype=np.float64)
        for sid in cells["sample_id"].unique():
            df_s = cells_with_meta[cells_with_meta["sample_id"] == sid]
            local_idx = df_s.index.values
            all_feat[local_idx] = compute_features_from_shuffled(
                df_s, args.k_nn, rng, nonepi_to_idx)
        U_perm = project_batched(np.clip(all_feat, 0, None), H)
        U_perm_norm = U_perm / (U_perm.sum(axis=1, keepdims=True) + 1e-12)
        perm_argmax_all.append(np.argmax(U_perm_norm, axis=1))
        perm_loadings_all.append(U_perm_norm)

    perm_argmax = perm_argmax_all[0]
    perm_loadings = perm_loadings_all[0]

    out_df = pd.DataFrame({
        "sample_id": cells["sample_id"].values,
        "patient_id": cells["patient_id"].values,
        "instance_id": cells["instance_id"].values,
        "l1p5_short": cells["l1p5_short"].values,
        "real_program": [f"P{i}" for i in real_argmax],
        "perm_program": [f"P{i}" for i in perm_argmax],
        "retained": real_argmax == perm_argmax,
        "real_max_load": U_real_norm.max(axis=1),
        "perm_max_load": perm_loadings.max(axis=1),
    })
    out_df.to_parquet(
        os.path.join(out_dir, "s37_per_cell_real_vs_perm.parquet"), index=False)

    stab_rows = []
    for p in range(k):
        m = real_argmax == p
        n_real = int(m.sum())
        n_retain = int((perm_argmax[m] == p).sum())
        conf = (pd.Series(perm_argmax[m])
                .value_counts(normalize=True).round(4).to_dict())
        stab_rows.append({
            "real_program": f"P{p}",
            "top_l1p5": NONEPI_21[int(np.argmax(H[p]))],
            "n_cells_real": n_real,
            "n_retain_under_perm": n_retain,
            "retention_frac": n_retain / max(n_real, 1),
            "perm_argmax_distribution": json.dumps(
                {f"P{int(k_)}": float(v) for k_, v in conf.items()}),
        })
    pd.DataFrame(stab_rows).to_csv(
        os.path.join(out_dir, "s37_program_stability.csv"), index=False)

    n_total = len(cells)
    cohort_concord = float((real_argmax == perm_argmax).sum()) / n_total

    audit = {
        "n_cells": int(n_total),
        "n_perms": int(args.n_perms),
        "cohort_concordance": cohort_concord,
        "per_program_retention": [
            {k_: v for k_, v in r.items() if k_ != "perm_argmax_distribution"}
            for r in stab_rows
        ],
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(out_dir, "s37_audit.json"), "w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    print(json.dumps({"cohort_concordance": cohort_concord,
                      "n_cells": int(n_total),
                      "elapsed_sec": round(time.time() - t0, 1)}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
