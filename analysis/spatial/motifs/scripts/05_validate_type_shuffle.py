#!/usr/bin/env python3
"""Cell-type null for the cell-level density-weighted NMF.

Complementary to the xy-shuffle null (script 03). For each Xenium SF instance,
the `l1p5_short` labels of the cells inside that instance are permuted while
their xy positions are held fixed. Each cell's KNN composition is then
recomputed under the shuffled identity assignment and projected through the
canonical H matrix to obtain a shuffled argmax program.

Interpretation: xy-shuffle (script 03) destroys spatial structure but preserves
local cell-type frequencies, so a program whose definition is essentially "this
cell type tends to be here" will retain. Type-shuffle preserves spatial
structure but destroys cell-type identity, so all programs are expected to
collapse, with retention near 1/k for cells in instances of moderate size. The
two nulls together demonstrate that programs depend jointly on spatial
arrangement and on cell-type identity, which is the requirement for a
microanatomic motif.

Outputs (data/validation/s05_type_shuffle/):
  - type_shuffle_per_cell.parquet     per-cell real argmax vs type-shuffled argmax
  - type_shuffle_program_stability.csv per-program retention fraction
  - type_shuffle_audit.json           cohort concordance + per-program summary
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


def compute_features_from_type_shuffled(df_sample, k_nn, rng, nonepi_to_idx):
    """Within each instance, permute labels (keep xy fixed); recompute KNN
    composition under the shuffled labels.
    """
    n = len(df_sample)
    feat = np.zeros((n, len(NONEPI_21)), dtype=np.float32)
    labels_orig = df_sample["l1p5_short"].to_numpy()
    inst_arr = df_sample["instance_id"].to_numpy()
    xy = df_sample[["x", "y"]].to_numpy(dtype=np.float32)

    labels_shuf = labels_orig.copy()
    for inst_id in np.unique(inst_arr):
        idx = np.where(inst_arr == inst_id)[0]
        if len(idx) < 2:
            continue
        perm = rng.permutation(len(idx))
        labels_shuf[idx] = labels_orig[idx[perm]]

    for inst_id in np.unique(inst_arr):
        idx = np.where(inst_arr == inst_id)[0]
        if len(idx) < 2:
            continue
        sub_xy = xy[idx]
        tree = cKDTree(sub_xy)
        k_eff = min(k_nn + 1, len(idx))
        _, nn_idx = tree.query(sub_xy, k=k_eff)
        nn_idx = nn_idx[:, 1:]
        for row_local in range(len(idx)):
            row_global = idx[row_local]
            neigh_labels = labels_shuf[idx[nn_idx[row_local]]]
            for nl in neigh_labels:
                if nl in nonepi_to_idx:
                    feat[row_global, nonepi_to_idx[nl]] += 1

    row_sum = feat.sum(axis=1, keepdims=True)
    nz = row_sum.flatten() > 0
    feat_frac = np.zeros_like(feat)
    feat_frac[nz] = feat[nz] / row_sum[nz]
    return feat_frac


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
                               "s05_type_shuffle"))
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

    perm_argmax_all = []
    perm_loadings_all = []
    for _ in range(args.n_perms):
        all_feat = np.zeros_like(X_real, dtype=np.float64)
        for sid in cells["sample_id"].unique():
            df_s = cells[cells["sample_id"] == sid]
            local_idx = df_s.index.values
            all_feat[local_idx] = compute_features_from_type_shuffled(
                df_s, args.k_nn, rng, nonepi_to_idx)
        U_perm = project_batched(np.clip(all_feat, 0, None), H)
        U_perm_norm = U_perm / (U_perm.sum(axis=1, keepdims=True) + 1e-12)
        perm_argmax_all.append(np.argmax(U_perm_norm, axis=1))
        perm_loadings_all.append(U_perm_norm)

    perm_argmax = perm_argmax_all[0]
    perm_loadings = perm_loadings_all[0]

    pd.DataFrame({
        "sample_id": cells["sample_id"].values,
        "patient_id": cells["patient_id"].values,
        "instance_id": cells["instance_id"].values,
        "l1p5_short": cells["l1p5_short"].values,
        "real_program": [f"P{i}" for i in real_argmax],
        "type_shuf_program": [f"P{i}" for i in perm_argmax],
        "retained": real_argmax == perm_argmax,
        "real_max_load": U_real_norm.max(axis=1),
        "type_shuf_max_load": perm_loadings.max(axis=1),
    }).to_parquet(
        os.path.join(out_dir, "type_shuffle_per_cell.parquet"), index=False)

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
            "n_retain_under_type_shuffle": n_retain,
            "retention_frac": n_retain / max(n_real, 1),
            "type_shuf_argmax_distribution": json.dumps(
                {f"P{int(k_)}": float(v) for k_, v in conf.items()}),
        })
    pd.DataFrame(stab_rows).to_csv(
        os.path.join(out_dir, "type_shuffle_program_stability.csv"), index=False)

    n_total = len(cells)
    cohort_concord = float((real_argmax == perm_argmax).sum()) / n_total

    audit = {
        "null_type": "l1p5_within_instance_permutation",
        "n_cells": int(n_total),
        "n_perms": int(args.n_perms),
        "cohort_concordance": cohort_concord,
        "per_program_retention": [
            {kk: vv for kk, vv in r.items()
             if kk != "type_shuf_argmax_distribution"}
            for r in stab_rows
        ],
        "elapsed_sec": round(time.time() - t0, 1),
    }
    with open(os.path.join(out_dir, "type_shuffle_audit.json"), "w") as fh:
        json.dump(audit, fh, indent=2, default=str)

    print(json.dumps({"cohort_concordance": cohort_concord,
                      "n_cells": int(n_total),
                      "elapsed_sec": round(time.time() - t0, 1)}), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
