#!/usr/bin/env python3
"""Whole-sample inventory render — L1.5 cell-type + motif argmax per sample.

For each Xenium sample, render two scatter panels:
  - L1.5 cell type (categorical color by l1p5_short)
  - motif argmax (categorical color by program)

Aspect ratio preserved per-sample (samples have variable shape/extent).
Output: per-sample PDF pair `{sample_id}__inventory_{l1p5,argmax}.pdf`.

These are the inventory tier of the FOV composition (whole-tissue scale,
not zoomed). The score gallery uses the zoomed render_from_saved_views.py
output instead.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Reuse the renderer's palettes for consistency with 3d/3e/score gallery
sys.path.insert(0, str(Path(__file__).parent))
from render_from_saved_views import L1P5_PALETTE, MOTIF_PAL, ARTIFACT


def render_sample(jc_s, sample_id, out_dir, long_axis_in=4.0,
                   point_size=0.10, dpi=600):
    """Render L1.5 + argmax panels for one sample at whole-tissue scale.

    Aspect ratio preserved by computing tile dims from bbox of cell positions.
    """
    x = jc_s["x"].to_numpy()
    y = jc_s["y"].to_numpy()
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    if x_range > y_range:
        w = long_axis_in
        h = long_axis_in * (y_range / x_range)
    else:
        h = long_axis_in
        w = long_axis_in * (x_range / y_range)

    # --- L1.5 panel ---
    fig, ax = plt.subplots(figsize=(w, h), dpi=dpi, facecolor="black")
    ax.set_facecolor("black")
    colors_l1p5 = jc_s["l1p5_short"].map(L1P5_PALETTE).fillna(ARTIFACT)
    ax.scatter(x, y, c=colors_l1p5, s=point_size, linewidth=0, rasterized=True)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.axis("off")
    fig.savefig(out_dir / f"{sample_id}__inventory_l1p5.pdf",
                 bbox_inches="tight", pad_inches=0.02, dpi=dpi,
                 facecolor="black")
    plt.close(fig)

    # --- argmax panel ---
    fig, ax = plt.subplots(figsize=(w, h), dpi=dpi, facecolor="black")
    ax.set_facecolor("black")
    def program_to_color(p):
        if isinstance(p, str) and (p.startswith("P") or p.startswith("M")):
            try:
                idx = int(p[1:])
                return MOTIF_PAL[f"M{idx}"]
            except (ValueError, KeyError):
                return ARTIFACT
        return ARTIFACT
    colors_argmax = jc_s["program"].map(program_to_color)
    ax.scatter(x, y, c=colors_argmax, s=point_size, linewidth=0,
                rasterized=True)
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.axis("off")
    fig.savefig(out_dir / f"{sample_id}__inventory_argmax.pdf",
                 bbox_inches="tight", pad_inches=0.02, dpi=dpi,
                 facecolor="black")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--joint-cells", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--samples", default=None,
                    help="comma-separated sample_ids; default = all")
    ap.add_argument("--long-axis-in", type=float, default=4.0)
    ap.add_argument("--point-size", type=float, default=0.10)
    ap.add_argument("--dpi", type=int, default=600)
    args = ap.parse_args()

    print("[load] joint_cells", flush=True)
    jc = pd.read_parquet(args.joint_cells)
    print(f"       {len(jc):,} cells, {jc['sample_id'].nunique()} samples",
          flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    samples = (args.samples.split(",") if args.samples
                else sorted(jc["sample_id"].dropna().unique()))
    for sid in samples:
        sub = jc[jc["sample_id"] == sid]
        if len(sub) < 100:
            print(f"[skip] {sid}: only {len(sub)} cells", flush=True)
            continue
        print(f"[render] {sid}: {len(sub):,} cells", flush=True)
        render_sample(sub, sid, out_dir, args.long_axis_in,
                       args.point_size, args.dpi)


if __name__ == "__main__":
    main()
