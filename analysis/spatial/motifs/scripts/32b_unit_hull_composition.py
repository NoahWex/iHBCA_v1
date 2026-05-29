#!/usr/bin/env python3
"""Per-unit full-L1.5 hull composition via point-in-polygon over alpha hulls.

Companion to 32_motif_units.py. The canonical motif_units.parquet defines
`l1p5_composition_json` as the L1.5 counts of cluster cells (those with
argmax_P == motif) — a descriptor of motif identity. This script adds a
separate descriptor of the unit's full spatial footprint by point-in-
polygon-testing every Xenium cell (any cell type, any motif assignment,
including epi) against each unit's stored alpha_shape_wkt.

The DIFFERENCE between the two descriptors is informative:
  - cluster_composition (l1p5_composition_json in motif_units.parquet) =
    cells assigned to THIS motif by argmax
  - hull_composition (this script) = cells in the unit's spatial polygon
    REGARDLESS of motif assignment, including epi cells
  - context_composition = hull_composition - cluster_composition (cells
    that share the unit's footprint but are NOT argmax-assigned to it).
    Captures the immediate cellular surroundings of the motif's
    signal carriers.

Inputs:
  - motif_units.parquet (canonical; contains alpha_shape_wkt + sample_id +
    unit_id + alpha_area_um2)
  - xenium_obs.csv (cell_id, xenium_id == sample_id, x_centroid, y_centroid)
  - joint_l1p5.csv (cell_id, platform, l1p5_short, is_artifact)

Outputs (same out-dir as 32_motif_units.py by default):
  - unit_hull_composition.parquet (one row per unit; columns: unit_id,
    sample_id, motif, alpha_area_um2, hull_n_total, hull_n_epi,
    hull_n_nonepi, hull_density_per_um2, hull_composition_json)
  - unit_hull_composition.audit.json
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from shapely import wkt as shapely_wkt
from shapely.geometry import Point
from shapely.prepared import prep


EPI_TYPES = {"LASP-basal", "LASP", "LHS", "BMYO-myo"}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--motif-units", required=True)
    ap.add_argument("--xenium-obs", required=True)
    ap.add_argument("--joint-l1p5", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--no-artifact-filter", action="store_true",
                    help="Include is_artifact==True cells (default: drop)")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    t0 = time.time()

    print(f"[load] motif_units.parquet", flush=True)
    units = pd.read_parquet(args.motif_units)
    if "alpha_shape_wkt" not in units.columns:
        sys.exit("ERROR: motif_units.parquet missing alpha_shape_wkt")
    print(f"  {len(units):,} units", flush=True)

    print(f"[load] joint_l1p5.csv", flush=True)
    l1p5 = pd.read_csv(args.joint_l1p5,
                        usecols=["cell_id", "platform", "l1p5_short",
                                 "is_artifact"])
    if not args.no_artifact_filter:
        l1p5 = l1p5[~l1p5["is_artifact"].fillna(False)]
    xen = l1p5[l1p5["platform"] == "xenium"].copy()
    print(f"  {len(xen):,} xenium cells (all L1.5)", flush=True)

    print(f"[load] xenium_obs.csv", flush=True)
    obs = pd.read_csv(args.xenium_obs,
                       usecols=["cell_id", "xenium_id",
                                "x_centroid", "y_centroid"])
    xen = xen.merge(obs, on="cell_id", how="inner")
    print(f"  {len(xen):,} cells with xy", flush=True)

    # Track full L1.5 vocabulary observed (drives column order in output)
    all_l1p5 = sorted(xen["l1p5_short"].dropna().unique())
    print(f"  L1.5 vocabulary ({len(all_l1p5)} types): {all_l1p5}",
          flush=True)

    rows = []
    samples = sorted(units["sample_id"].dropna().unique())
    for s_idx, sid in enumerate(samples, start=1):
        u_sid = units[units["sample_id"] == sid]
        x_sid = xen[xen["xenium_id"] == sid]
        if len(u_sid) == 0 or len(x_sid) == 0:
            continue
        x_xy = x_sid[["x_centroid", "y_centroid"]].to_numpy()
        x_pts = [Point(float(x), float(y)) for x, y in x_xy]
        x_l1p5 = x_sid["l1p5_short"].to_numpy()

        for _, ur in u_sid.iterrows():
            try:
                geom = shapely_wkt.loads(str(ur["alpha_shape_wkt"]))
            except Exception as e:
                print(f"  WARN: unit {ur['unit_id']} WKT parse failed: {e}",
                      flush=True)
                continue
            prepared = prep(geom)
            inside = np.array([prepared.contains(p) for p in x_pts])
            if inside.sum() == 0:
                comp = {k: 0 for k in all_l1p5}
            else:
                vc = pd.Series(x_l1p5[inside]).value_counts().to_dict()
                comp = {k: int(vc.get(k, 0)) for k in all_l1p5}
            n_total = int(sum(comp.values()))
            n_epi = int(sum(v for k, v in comp.items() if k in EPI_TYPES))
            rows.append({
                "unit_id": ur["unit_id"],
                "sample_id": sid,
                "motif": ur["motif"],
                "alpha_area_um2": float(ur["alpha_area_um2"]),
                "hull_n_total":  n_total,
                "hull_n_epi":    n_epi,
                "hull_n_nonepi": n_total - n_epi,
                "hull_density_per_um2": (
                    n_total / float(ur["alpha_area_um2"])
                    if ur["alpha_area_um2"] > 0 else np.nan),
                "hull_composition_json": json.dumps(comp),
            })
        if s_idx % 5 == 0 or s_idx == len(samples):
            print(f"  [{s_idx:>2}/{len(samples)}] {sid}: "
                  f"{len(u_sid)} units done", flush=True)

    out = pd.DataFrame(rows)
    out_path = os.path.join(args.out_dir, "unit_hull_composition.parquet")
    out.to_parquet(out_path, index=False)
    print(f"\n[write] {out_path} ({len(out):,} rows)", flush=True)

    audit = {
        "n_units_input":         int(len(units)),
        "n_units_with_hull_hits": int((out["hull_n_total"] > 0).sum()),
        "n_units_empty":         int((out["hull_n_total"] == 0).sum()),
        "n_samples":             int(len(samples)),
        "mean_cells_per_unit":   float(out["hull_n_total"].mean()),
        "median_cells_per_unit": float(out["hull_n_total"].median()),
        "mean_frac_epi":         float(
            (out["hull_n_epi"] / out["hull_n_total"].clip(lower=1)).mean()),
        "l1p5_vocabulary":       list(all_l1p5),
        "elapsed_sec":           float(time.time() - t0),
    }
    audit_path = os.path.join(args.out_dir, "unit_hull_composition.audit.json")
    with open(audit_path, "w") as f:
        json.dump(audit, f, indent=2)
    print(f"[write] {audit_path}", flush=True)
    print(f"[done] {audit['elapsed_sec']:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
