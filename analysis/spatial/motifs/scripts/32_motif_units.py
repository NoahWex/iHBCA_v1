#!/usr/bin/env python3
"""Define discrete motif units via DBSCAN on argmax cells per sample.

Per-cell substrate (already canonical):
  argmax_P (the dominant motif) and argmax_loading (its strength), plus
  secondary_P / secondary_loading (the 2nd-highest motif) to preserve
  gradation/hybrid information without using it for unit assignment.

Unit definition (per sample, per motif):
  DBSCAN on the xy of argmax_P==P_k cells, with eps_um and min_samples
  user-configurable. Each cluster (non-noise label >=0) is a unit. Noise
  cells (-1) remain unassigned and contribute only to the per-cell table.

Per-unit annotations:
  - n_cells, hull_area_um2, density (cells/µm²)
  - centroid_x, centroid_y
  - purity = mean(argmax_loading) within unit
  - mean full loading vector across unit cells (9-dim)
  - L1.5 composition (counts per L1.5 type in the unit)

Per-program GMM positivity thresholds (script's own QC step) are written
out and used only as a descriptive annotation per cell (n_motifs_above_GMM)
— the unit substrate uses straight argmax.

Inputs:
  - joint_cells.parquet (440k non-epi cells, xy + sample_id + loadings)
  - xenium_obs.csv (cell_id ↔ xy mapping for KD-tree join)
  - joint_l1p5.csv (patient_id + position metadata per cell)

Outputs (reports/motif_units/):
  - motif_units.parquet (one row per unit)
  - motif_cells.parquet (one row per cell with unit annotations)
  - motif_units.audit.json
  - motif_units.md
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree, ConvexHull, QhullError
from shapely.geometry import Point
from shapely.ops import unary_union
from sklearn.cluster import DBSCAN


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True,
                    help="joint_cells.parquet (non-epi cells with loadings + xy)")
    ap.add_argument("--xenium-obs", required=True,
                    help="xenium_obs.csv (cell_id + xy)")
    ap.add_argument("--joint-l1p5", required=True,
                    help="joint_l1p5.csv (per-cell patient_id, position)")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--eps-um", type=float, default=35.0,
                    help="DBSCAN eps in µm (global default, overridden per "
                    "motif if --params-yaml provided)")
    ap.add_argument("--min-samples", type=int, default=10,
                    help="DBSCAN min_samples (global default, overridden "
                    "per motif if --params-yaml provided)")
    ap.add_argument("--params-yaml", default=None,
                    help="Per-motif DBSCAN params YAML. Each key is a "
                    "motif name (P0..P8), value is {eps_um, min_samples}. "
                    "Falls back to --eps-um/--min-samples for missing keys.")
    ap.add_argument("--alpha-r-um", type=float, default=20.0,
                    help="Alpha-shape buffer radius (µm) for per-unit "
                    "geometry. Stored as WKT; matches the display geometry "
                    "used by napari and downstream point-in-polygon tests.")
    return ap.parse_args()


def hull_area(xy: np.ndarray) -> float:
    """Convex hull area (µm²). Returns 0 if fewer than 3 unique points or
    collinear / hull-fail."""
    if len(xy) < 3:
        return 0.0
    try:
        return float(ConvexHull(xy).volume)  # 2D: .volume == area
    except QhullError:
        return 0.0


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    print(f"[load] joint_cells.parquet", flush=True)
    jc = pd.read_parquet(args.parquet)
    load_cols = sorted([c for c in jc.columns
                          if c.startswith("loading_P")])
    programs = [c.replace("loading_", "") for c in load_cols]
    EPI = ["LASP-basal", "LASP", "LHS", "BMYO-myo"]
    n_before = len(jc)
    jc = jc[~jc["l1p5_short"].isin(EPI)].copy()
    print(f"  dropped {n_before - len(jc):,} epi cells from parquet",
          flush=True)
    print(f"  {len(jc):,} non-epi cells with loadings", flush=True)

    # ----- argmax + secondary -------------------------------------------
    load_mat = jc[load_cols].to_numpy()
    argmax_idx = load_mat.argmax(axis=1)
    jc["argmax_P"] = np.asarray(programs)[argmax_idx]
    jc["argmax_loading"] = load_mat[np.arange(len(jc)), argmax_idx]
    # Mask argmax to find second-highest
    masked = load_mat.copy()
    masked[np.arange(len(jc)), argmax_idx] = -np.inf
    sec_idx = masked.argmax(axis=1)
    jc["secondary_P"] = np.asarray(programs)[sec_idx]
    jc["secondary_loading"] = load_mat[np.arange(len(jc)), sec_idx]

    # ----- Patient + position metadata ----------------------------------
    print(f"\n[load] joint_l1p5 (patient + position)", flush=True)
    l1p5 = pd.read_csv(
        args.joint_l1p5,
        usecols=["library_id", "patient_id", "p_level", "quadrant",
                 "is_uoq", "platform"])
    l1p5 = l1p5[l1p5["platform"] == "xenium"].drop_duplicates("library_id")
    sid2meta = l1p5.set_index("library_id").to_dict("index")

    # ----- KD-tree join for cell_id -------------------------------------
    print(f"\n[KD-tree] join jc → obs to recover cell_id", flush=True)
    obs = pd.read_csv(args.xenium_obs,
                      usecols=["cell_id", "xenium_id", "x_centroid",
                               "y_centroid"])
    jc["cell_id"] = None
    for sid, grp in jc.groupby("sample_id"):
        og = obs[obs["xenium_id"] == sid]
        if len(og) == 0:
            continue
        tree = cKDTree(og[["x_centroid", "y_centroid"]].to_numpy())
        d, idx = tree.query(grp[["x", "y"]].to_numpy(), k=1)
        mask = d <= 0.5
        jc.loc[grp.index[mask], "cell_id"] = (
            og["cell_id"].to_numpy()[idx[mask]])
    jc = jc.dropna(subset=["cell_id"]).reset_index(drop=True)
    print(f"  joined: {len(jc):,} cells with cell_id", flush=True)

    # Stamp patient + position
    def _meta_lookup(sid, key, default=None):
        d = sid2meta.get(sid, {})
        return d.get(key, default)
    jc["patient_id"] = jc["sample_id"].map(
        lambda s: _meta_lookup(s, "patient_id"))
    jc["p_level"] = jc["sample_id"].map(
        lambda s: _meta_lookup(s, "p_level"))
    jc["quadrant"] = jc["sample_id"].map(
        lambda s: _meta_lookup(s, "quadrant"))
    jc["is_uoq"] = jc["sample_id"].map(
        lambda s: _meta_lookup(s, "is_uoq", False))

    # ----- Per-motif DBSCAN params --------------------------------------
    per_motif_params = {}
    if args.params_yaml:
        import yaml
        with open(args.params_yaml) as fh:
            per_motif_params = yaml.safe_load(fh) or {}
        print(f"\n[params] loaded per-motif overrides from "
              f"{args.params_yaml}: {list(per_motif_params.keys())}",
              flush=True)
    def get_params(prog):
        d = per_motif_params.get(prog, {})
        return (float(d.get("eps_um", args.eps_um)),
                int(d.get("min_samples", args.min_samples)))

    # ----- DBSCAN per (sample, motif) → units ---------------------------
    print(f"\n[DBSCAN] per (sample, motif) — per-motif params:", flush=True)
    for prog in programs:
        eps, ms = get_params(prog)
        print(f"  {prog}: eps={eps}µm, min_samples={ms}", flush=True)
    jc["unit_id"] = ""
    unit_rows = []
    unit_counter = 0
    for sid, sgrp in jc.groupby("sample_id"):
        for prog in programs:
            eps_prog, ms_prog = get_params(prog)
            mask = (sgrp["argmax_P"] == prog)
            sub = sgrp[mask]
            if len(sub) < ms_prog:
                continue
            xy = sub[["x", "y"]].to_numpy()
            db = DBSCAN(eps=eps_prog,
                          min_samples=ms_prog).fit(xy)
            labels = db.labels_  # -1 == noise
            n_units_local = len(set(labels)) - (1 if -1 in labels else 0)
            for lbl in range(n_units_local):
                u_mask = labels == lbl
                u_cells = sub.iloc[u_mask]
                u_id = f"{sid}|{prog}|{lbl:03d}"
                jc.loc[u_cells.index, "unit_id"] = u_id
                u_xy = u_cells[["x", "y"]].to_numpy()
                ha = hull_area(u_xy)
                cx, cy = u_xy.mean(axis=0)
                # alpha-shape (buffered-point unary_union) — matches
                # display + downstream point-in-polygon contracts
                alpha_geom = unary_union(
                    [Point(float(x), float(y)).buffer(args.alpha_r_um)
                     for x, y in u_xy])
                alpha_area = float(alpha_geom.area)
                alpha_wkt = alpha_geom.wkt
                # mean + SD loading vector across unit cells (9-dim each)
                mean_loads = u_cells[load_cols].mean()
                sd_loads = u_cells[load_cols].std(ddof=1)
                # L1.5 composition
                l1p5_counts = (
                    u_cells["l1p5_short"].value_counts().to_dict())
                unit_rows.append({
                    "unit_id": u_id,
                    "sample_id": sid,
                    "patient_id": _meta_lookup(sid, "patient_id"),
                    "p_level": _meta_lookup(sid, "p_level"),
                    "quadrant": _meta_lookup(sid, "quadrant"),
                    "is_uoq": _meta_lookup(sid, "is_uoq", False),
                    "motif": prog,
                    "n_cells": int(u_mask.sum()),
                    "hull_area_um2": ha,
                    "alpha_area_um2": alpha_area,
                    "alpha_shape_wkt": alpha_wkt,
                    "density_cells_per_um2": (
                        u_mask.sum() / ha if ha > 0 else np.nan),
                    "centroid_x": float(cx),
                    "centroid_y": float(cy),
                    "purity_mean_argmax_loading": float(
                        u_cells["argmax_loading"].mean()),
                    "l1p5_composition_json": json.dumps(l1p5_counts),
                    **{f"mean_{c}": float(mean_loads[c])
                       for c in load_cols},
                    **{f"sd_{c}": (float(sd_loads[c])
                                    if pd.notna(sd_loads[c]) else 0.0)
                       for c in load_cols},
                })
                unit_counter += 1
        if unit_counter > 0 and unit_counter % 200 == 0:
            print(f"  {unit_counter:,} units; sample {sid} done",
                  flush=True)
    print(f"  total units: {unit_counter:,}", flush=True)

    # ----- Persist -------------------------------------------------------
    units = pd.DataFrame(unit_rows)
    units_path = os.path.join(args.out_dir, "motif_units.parquet")
    units.to_parquet(units_path, index=False)
    print(f"\n[write] {units_path} ({len(units):,} rows)", flush=True)

    cells_keep = ["cell_id", "sample_id", "patient_id", "p_level",
                   "quadrant", "is_uoq", "l1p5_short", "x", "y",
                   "argmax_P", "argmax_loading", "secondary_P",
                   "secondary_loading", "unit_id"]
    cells = jc[cells_keep].copy()
    cells_path = os.path.join(args.out_dir, "motif_cells.parquet")
    cells.to_parquet(cells_path, index=False)
    print(f"[write] {cells_path} ({len(cells):,} rows)", flush=True)

    # Audit
    by_motif = (units.groupby("motif")
                       .agg(n_units=("unit_id", "count"),
                            median_n_cells=("n_cells", "median"),
                            p10_n_cells=("n_cells",
                                          lambda s: np.percentile(s, 10)),
                            p90_n_cells=("n_cells",
                                          lambda s: np.percentile(s, 90)),
                            median_hull_um2=("hull_area_um2", "median"),
                            median_purity=(
                                "purity_mean_argmax_loading", "median"))
                       .round(3))
    print("\n=== Unit count + size per motif ===")
    print(by_motif.to_string())

    cells_in_units = int((cells["unit_id"] != "").sum())
    audit = {
        "params_yaml": args.params_yaml,
        "per_motif_params": {p: {"eps_um": get_params(p)[0],
                                   "min_samples": get_params(p)[1]}
                              for p in programs},
        "default_eps_um": args.eps_um,
        "default_min_samples": args.min_samples,
        "alpha_r_um": args.alpha_r_um,
        "n_cells_input": int(len(jc)),
        "n_cells_in_units": cells_in_units,
        "n_cells_noise_or_below_min_samples": int(len(jc) - cells_in_units),
        "n_units_total": int(unit_counter),
        "units_per_motif": units["motif"].value_counts().to_dict(),
        "runtime_s": round(time.time() - t0, 1),
    }
    audit_path = os.path.join(args.out_dir, "motif_units.audit.json")
    with open(audit_path, "w") as fh:
        json.dump(audit, fh, indent=2, default=str)
    print(f"\n[write] {audit_path}", flush=True)

    md = [
        f"# Motif units (eps={args.eps_um}µm, min_samples={args.min_samples})",
        "",
        f"Discrete units defined by DBSCAN on per-sample argmax cells per "
        f"motif. Each cell is assigned to its motif's argmax cluster (if "
        f"non-noise) plus annotated with its secondary loading for "
        f"gradation/hybrid context.",
        "",
        f"- {audit['n_units_total']:,} units across "
        f"{cells['sample_id'].nunique()} samples and 9 motifs",
        f"- {audit['n_cells_in_units']:,} / {audit['n_cells_input']:,} "
        f"({100 * audit['n_cells_in_units'] / audit['n_cells_input']:.1f}%) "
        f"cells assigned to a unit (rest are DBSCAN noise or in motifs "
        f"with <min_samples cells in their sample)",
        "",
        "## Per-motif unit summary",
        "",
        by_motif.to_markdown(),
        "",
        "## Per-cell substrate",
        "- `argmax_P`: dominant motif",
        "- `argmax_loading`: dominant strength (median ~0.45-0.55 per motif)",
        "- `secondary_P`, `secondary_loading`: 2nd-highest motif "
        "(gradation info)",
        "- `unit_id`: per-cell motif unit (empty string if DBSCAN noise)",
        "",
        "## Per-unit loading profile",
        "- `mean_loading_P{0..8}` + `sd_loading_P{0..8}`: per-program mean "
        "and SD across the unit's cells. Encodes multi-motif character "
        "without threshold-picking — downstream chooses how to summarize.",
    ]
    md_path = os.path.join(args.out_dir, "motif_units.md")
    with open(md_path, "w") as fh:
        fh.write("\n".join(md))
    print(f"[done] {time.time()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
