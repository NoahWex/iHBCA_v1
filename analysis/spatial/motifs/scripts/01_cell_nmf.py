#!/usr/bin/env python3
"""Cell-level density-weighted NMF on non-epithelial L1.5 neighborhoods.

For every Xenium cell within a SpaceFlow (SF) structure instance:
  1. Identify K=10 spatial nearest neighbors among all cells inside the same instance.
  2. Build a 21-dimensional feature vector recording the fraction of those
     neighbors carrying each non-epithelial L1.5 type. Epithelial neighbors are
     summarised separately as a scalar `epi_context` so the factor model is not
     dominated by the epi axis.
  3. Drop cells whose feature vector is all-zero (i.e. K nearest neighbors are
     entirely epithelial). The remaining cells form the NMF input.
  4. Row-weight each observation by sqrt(density-percentile-rank * (1 - epi_context))
     so factor learning is driven by cells in dense, genuinely non-epi pockets.

A k-sweep over k in [k_min, k_max] is run; k is picked by maximum second-difference
of reconstruction error (knee selection). For the canonical run the procedure
selects k=9 programs. The basis matrix H is row-normalised so each program is a
distribution over 21 non-epi L1.5 types.

Inputs (resolved via config/paths.yaml):
  - joint_l1p5 CSV (cell-level L1.5 labels, Xenium subset)
  - manifest TSV (Xenium sample roster)
  - annotated h5ad directory (per-sample Xenium objects with xy + structure ids)

Outputs (outputs/):
  - basis.csv                         k x 21 H matrix (program x L1.5 distributions)
  - diagnostics.csv                   per-k recon error + second-difference
  - joint_l1p5_xenium_cells_v2.parquet per-cell features + W loadings + dominant program
  - program_epi_context_stats.csv     per-program orthogonal epi-context summary
  - nmf_input_features.parquet        persisted X_w + sqw for rank-stability refits

Notes for reviewers:
  - The canonical run is deterministic for a fixed --random-state and a fixed
    n-init. The substrate at outputs/basis.csv is the output of this script
    with k=9 and the defaults below.
  - Validation companions: 02 (seed stability across 20 random inits), 03
    (within-instance xy-shuffle null), 04 (leave-one-patient-out hold-out),
    05 (within-instance cell-type-label shuffle null). Rank stability across
    neighboring K is assessed by 02b against the persisted features.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from sklearn.decomposition import NMF


EPI_TYPES = ["LASP-basal", "LASP", "LHS", "BMYO-myo"]
NONEPI_21 = [
    "Fb", "Fb_Activated", "Fb_SFRP4", "Adipo",
    "Plas", "T-NK", "cDC1", "CD4T", "Neu", "Mast", "Treg",
    "pDC", "Mac", "cDC2", "cDC", "CD8T", "B",
    "EC", "PV", "LEC", "Vas-cap",
]
ARTIFACT_LABELS = {"Stromal_art", "Immune_art", "BMYO_art",
                   "Epithelial_art", "Mast_art", "NA", "nan"}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", required=True,
                    help="Pipeline root (contains config/paths.yaml)")
    ap.add_argument("--config", default=None,
                    help="Path to paths.yaml (default: ${project_root}/config/paths.yaml)")
    ap.add_argument("--joint-l1p5", required=True,
                    help="CSV with cell_id, platform, library_id, l1p5_short, is_uoq")
    ap.add_argument("--annotated-dir", required=True,
                    help="Directory with {sample}_annotated.h5ad files (per-Xenium)")
    ap.add_argument("--manifest", required=True,
                    help="TSV with xenium_id + patient_id columns")
    ap.add_argument("--out-dir", default=None,
                    help="Output dir (default: ${project_root}/data/nmf)")
    ap.add_argument("--k-nn", type=int, default=10)
    ap.add_argument("--k-min", type=int, default=4)
    ap.add_argument("--k-max", type=int, default=12)
    ap.add_argument("--n-init", type=int, default=15)
    ap.add_argument("--random-state", type=int, default=42)
    return ap.parse_args()


def per_sample_features(sample_id, h5ad_path, js_sample, k_nn):
    """Build per-cell features for one Xenium sample.

    Returns a DataFrame keyed on cell_id with columns: xy, instance_id,
    container, epi_context, local_density, has_nonepi_neighbor, and 21
    nonepi_neigh_frac_* columns summing to 1 (or 0 if no non-epi neighbors).
    """
    if not os.path.exists(h5ad_path):
        return pd.DataFrame()
    adata = ad.read_h5ad(h5ad_path)
    if "spatial" not in adata.obsm:
        return pd.DataFrame()
    xy = np.asarray(adata.obsm["spatial"]).astype(np.float32)
    sid_obs = adata.obs.index.astype(str).values

    js = js_sample.copy()
    js["barcode_tail"] = js["cell_id"].str.split("__").str[-1]

    obs_df = pd.DataFrame({
        "cell_id_full": sid_obs,
        "barcode_tail": pd.Series(sid_obs).str.split("__").str[-1].values,
        "x": xy[:, 0], "y": xy[:, 1],
        "instance_id": adata.obs["structure_label_discrete"].astype(str).values,
        "container": adata.obs["structure_label"].astype(str).values,
    })
    merged = obs_df.merge(js[["barcode_tail", "l1p5_short"]],
                          on="barcode_tail", how="left")
    merged = merged[(merged["container"] != "Epidermal")
                    & merged["l1p5_short"].notna()].copy().reset_index(drop=True)
    if len(merged) < k_nn + 1:
        return pd.DataFrame()

    epi_set = set(EPI_TYPES)
    nonepi_to_idx = {l: i for i, l in enumerate(NONEPI_21)}
    labels = merged["l1p5_short"].to_numpy()
    xyc = merged[["x", "y"]].to_numpy(dtype=np.float32)
    inst_arr = merged["instance_id"].to_numpy()

    n = len(merged)
    feat = np.zeros((n, len(NONEPI_21)), dtype=np.float32)
    epi_context = np.zeros(n, dtype=np.float32)
    local_density = np.zeros(n, dtype=np.float32)

    for inst_id in np.unique(inst_arr):
        idx = np.where(inst_arr == inst_id)[0]
        if len(idx) < 2:
            continue
        sub_xy = xyc[idx]
        sub_tree = cKDTree(sub_xy)
        k_eff = min(k_nn + 1, len(idx))
        dists, nn_idx = sub_tree.query(sub_xy, k=k_eff)
        nn_idx = nn_idx[:, 1:]
        nn_dists = dists[:, 1:]
        for row_local in range(len(idx)):
            row_global = idx[row_local]
            neigh_labels = labels[idx[nn_idx[row_local]]]
            k_real = len(neigh_labels)
            n_epi = 0
            for nl in neigh_labels:
                if nl in epi_set:
                    n_epi += 1
                elif nl in nonepi_to_idx:
                    feat[row_global, nonepi_to_idx[nl]] += 1
            epi_context[row_global] = n_epi / k_real if k_real > 0 else 0
            r_k = nn_dists[row_local, -1]
            if r_k > 0:
                local_density[row_global] = k_real / (np.pi * r_k * r_k)

    row_sum = feat.sum(axis=1, keepdims=True)
    nonzero = row_sum > 0
    feat_frac = np.zeros_like(feat)
    feat_frac[nonzero[:, 0]] = feat[nonzero[:, 0]] / row_sum[nonzero[:, 0]]

    merged["sample_id"] = sample_id
    merged["epi_context"] = epi_context
    merged["local_density"] = local_density
    merged["nonepi_neigh_count"] = feat.sum(axis=1)
    merged["has_nonepi_neighbor"] = (feat.sum(axis=1) > 0)
    for i, l in enumerate(NONEPI_21):
        merged[f"nonepi_neigh_frac_{l}"] = np.round(feat_frac[:, i], 5)
    merged["is_artifact_cell"] = merged["l1p5_short"].isin(ARTIFACT_LABELS)
    return merged.rename(columns={"cell_id_full": "cell_id"}).drop(columns=["barcode_tail"])


def main() -> int:
    args = parse_args()
    sys.path.insert(0, os.path.join(args.project_root, "config"))
    from load_paths import load_paths, resolve  # noqa: E402

    paths = load_paths(args.project_root, args.config)
    out_dir = args.out_dir or resolve(paths, "data.nmf")
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(args.joint_l1p5):
        raise FileNotFoundError(f"joint_l1p5 not found: {args.joint_l1p5}")
    if not os.path.exists(args.manifest):
        raise FileNotFoundError(f"manifest not found: {args.manifest}")

    joint = pd.read_csv(args.joint_l1p5,
                        usecols=["cell_id", "platform", "library_id",
                                 "l1p5_short", "is_uoq"])
    joint = joint[joint["platform"] == "xenium"].copy()
    print(f"[joint] {len(joint):,} xenium cells", flush=True)

    manifest = pd.read_csv(args.manifest, sep="\t")
    pieces = []
    t0 = time.time()
    for _, row in manifest.iterrows():
        sid = row["xenium_id"]
        h5 = os.path.join(args.annotated_dir, f"{sid}_annotated.h5ad")
        js = joint[joint["library_id"] == sid][["cell_id", "l1p5_short"]]
        df = per_sample_features(sid, h5, js, args.k_nn)
        if len(df) > 0:
            pieces.append(df)
            if len(pieces) % 10 == 0:
                print(f"  ... {len(pieces)} samples done at "
                      f"{time.time()-t0:.1f}s", flush=True)
    cells = pd.concat(pieces, ignore_index=True)
    print(f"[cohort] {len(cells):,} cells in {time.time()-t0:.1f}s", flush=True)

    cells = cells.merge(
        manifest[["xenium_id", "patient_id"]].rename(columns={"xenium_id": "sample_id"}),
        on="sample_id", how="left")
    sample_uoq = (joint[["library_id", "is_uoq"]]
                  .drop_duplicates("library_id")
                  .rename(columns={"library_id": "sample_id"}))
    cells = cells.merge(sample_uoq, on="sample_id", how="left")

    cells_real = cells[(~cells["is_artifact_cell"])
                       & cells["has_nonepi_neighbor"]].copy().reset_index(drop=True)
    n_dropped = int((cells["has_nonepi_neighbor"] == False).sum())
    print(f"[NMF input] {len(cells_real):,} cells "
          f"(dropped {n_dropped:,} with all-epi K-nn neighborhoods)", flush=True)

    feat_cols = [f"nonepi_neigh_frac_{l}" for l in NONEPI_21]
    X = np.clip(cells_real[feat_cols].to_numpy(dtype=np.float64), 0, None)

    # Row weight: density percentile-rank x non-epi share, sqrt-transformed
    d = cells_real["local_density"].to_numpy(dtype=np.float64)
    d_log = np.log10(d + 1e-10)
    d_rank = np.argsort(np.argsort(d_log)) / max(len(d_log) - 1, 1)
    d_rank = np.clip(d_rank, 1e-3, 1.0)
    ec = cells_real["epi_context"].to_numpy(dtype=np.float64)
    nonepi_share = np.clip(1.0 - ec, 1e-3, 1.0)
    w = d_rank * nonepi_share
    sqw = np.sqrt(w).astype(np.float32)
    X_w = X * sqw[:, None]

    # Persist NMF input substrate so downstream alternate-K refits (e.g.
    # 02b_refit_at_k.py) can refit without rebuilding the KNN/density
    # substrate from scratch.
    pd.DataFrame(
        X_w, columns=[f"Xw_{l}" for l in NONEPI_21]
    ).assign(
        cell_id=cells_real["cell_id"].to_numpy(),
        sqw=sqw,
    ).to_parquet(
        os.path.join(out_dir, "nmf_input_features.parquet"), index=False)

    rng = np.random.RandomState(args.random_state)
    diag = []
    best_models = {}
    for k in range(args.k_min, args.k_max + 1):
        best_err = np.inf
        best_model = None
        for _ in range(args.n_init):
            s = int(rng.randint(0, 2**31 - 1))
            try:
                nmf = NMF(n_components=k, init="nndsvda", max_iter=400,
                          random_state=s, tol=1e-5)
                nmf.fit(X_w)
                err = nmf.reconstruction_err_
            except Exception:
                continue
            if err < best_err:
                best_err = err
                best_model = nmf
        diag.append({"k": k, "recon_err": float(best_err)})
        best_models[k] = best_model
        print(f"  k={k:2d} recon_err={best_err:.4f}", flush=True)

    errs = np.array([d_["recon_err"] for d_ in diag])
    ks = np.array([d_["k"] for d_ in diag])
    d2 = np.full(len(ks), np.nan)
    if len(ks) >= 3:
        d2[1:-1] = errs[2:] - 2 * errs[1:-1] + errs[:-2]
    diag_df = pd.DataFrame({"k": ks, "recon_err": errs, "second_diff": d2})
    diag_df.to_csv(os.path.join(out_dir, "diagnostics.csv"), index=False)

    k_chosen = (int(ks[int(np.nanargmax(d2))])
                if np.any(np.isfinite(d2)) else int(ks[np.argmin(errs)]))
    print(f"[picked k={k_chosen}]", flush=True)

    model = best_models[k_chosen]
    H = model.components_
    H_norm = H / (H.sum(axis=1, keepdims=True) + 1e-12)
    basis = pd.DataFrame(H_norm, columns=NONEPI_21)
    basis.insert(0, "program", [f"P{i}" for i in range(k_chosen)])
    basis.to_csv(os.path.join(out_dir, "basis.csv"), index=False)

    # Exit gate: basis shape and W row count match expectations
    assert basis.shape == (k_chosen, 22), f"basis shape unexpected: {basis.shape}"

    U_w = model.transform(X_w)
    U = U_w / sqw[:, None]
    U_norm = U / (U.sum(axis=1, keepdims=True) + 1e-12)
    cells_real["program"] = [f"P{i}" for i in np.argmax(U_norm, axis=1)]
    for i in range(k_chosen):
        cells_real[f"loading_P{i}"] = np.round(U_norm[:, i], 5)
    cells_real.to_parquet(
        os.path.join(out_dir, "joint_l1p5_xenium_cells_v2.parquet"),
        index=False)

    ec_rows = []
    for prog in cells_real["program"].unique():
        sub = cells_real[cells_real["program"] == prog]
        ec_rows.append({
            "program": prog,
            "n_cells": int(len(sub)),
            "mean_epi_context": round(float(sub["epi_context"].mean()), 4),
            "median_epi_context": round(float(sub["epi_context"].median()), 4),
            "frac_with_zero_epi": round(float((sub["epi_context"] == 0).mean()), 4),
            "frac_with_high_epi": round(float((sub["epi_context"] >= 0.5).mean()), 4),
            "mean_local_density": round(float(sub["local_density"].mean()), 6),
        })
    pd.DataFrame(ec_rows).sort_values("mean_epi_context").to_csv(
        os.path.join(out_dir, "program_epi_context_stats.csv"), index=False)

    print(json.dumps({
        "k_chosen": int(k_chosen),
        "basis_shape": list(basis.shape),
        "W_rows": int(len(cells_real)),
        "n_dropped_no_nonepi": n_dropped,
    }), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
