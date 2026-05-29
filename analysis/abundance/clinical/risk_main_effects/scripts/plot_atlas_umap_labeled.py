#!/usr/bin/env python3
"""
plot_atlas_umap_labeled.py - UMAP scatter with text labels at group medians.

Color cells by:
  --color-by ng_viable   (limma-viable NhoodGroups for a contrast; needs --contrast)
  --color-by L2          (parent_L2_joint per labels_full)

Labels placed at median(UMAP_1), median(UMAP_2) per group. Cells outside the
viable set rendered as small grey backdrop.

Args:
  --inquiry-dir
  --contrast (only used for color-by=ng_viable)
  --color-by {ng_viable,L2}
  --atlas-nhoods-dir (default <inquiry>/outputs/atlas_nhoods)
  --labels-full
  --out-pdf
"""
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path
from collections import Counter

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.load_paths import load_inquiry  # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--contrast", default=None)
    p.add_argument("--color-by", choices=["ng_viable", "L2"], required=True)
    p.add_argument("--atlas-nhoods-dir", type=Path, default=None)
    p.add_argument("--labels-full", default=None,
                   help="Override path; default resolves via paths.yaml inputs.labels")
    p.add_argument("--out-pdf", required=True)
    p.add_argument("--max-points", type=int, default=400_000,
                   help="downsample to this many cells for plotting (memory)")
    p.add_argument("--point-size", type=float, default=0.5)
    p.add_argument("--label-fontsize", type=float, default=5.5)
    return p.parse_args()


def make_palette(n: int) -> list[str]:
    """Distinguishable categorical palette via HSV stride."""
    np.random.seed(7)
    hues = (np.linspace(0, 1, n + 1)[:-1] + np.random.uniform(0, 0.04, n)) % 1.0
    perm = np.random.permutation(n)  # spread similar IDs
    hsv = np.zeros((n, 3))
    hsv[:, 0] = hues[perm]
    hsv[:, 1] = 0.65
    hsv[:, 2] = 0.85
    rgb = matplotlib.colors.hsv_to_rgb(hsv)
    return [matplotlib.colors.rgb2hex(c) for c in rgb]


def main():
    a = parse_args()
    inq = a.inquiry_dir
    if not a.labels_full:
        a.labels_full = load_inquiry(inq)["paths"]["inputs"]["labels"]
    nhoods_dir = a.atlas_nhoods_dir if a.atlas_nhoods_dir else (inq / "outputs/atlas_nhoods")

    # ---- atlas cells + UMAP ----
    cells_atlas = pd.read_csv(nhoods_dir / "cells.tsv", header=None)[0].values
    umap_df = pd.read_csv(nhoods_dir / "atlas_umap.tsv", sep="\t")
    print(f"atlas cells: {len(cells_atlas)}; umap_df: {len(umap_df)}")
    # Index UMAP by cell_id
    umap_df = umap_df.set_index("cell_id").loc[cells_atlas].reset_index()
    coords = umap_df[["UMAP_1", "UMAP_2"]].values

    # ---- Pre-load labels_full (used by both modes for L2 fallback / direct color) ----
    labels = pd.read_csv(a.labels_full,
                          usecols=["cell_id", "compartment", "label", "is_artifact"])
    labels["is_artifact"] = labels["is_artifact"].astype(str).str.lower().isin(
        {"true", "1"})
    labels["L2_joint"] = (labels["compartment"].astype(str) + "::"
                            + labels["label"].astype(str))
    labels.loc[labels.is_artifact, "L2_joint"] = "__artifact__"
    cell_to_l2 = dict(zip(labels.cell_id.values, labels.L2_joint.values))
    l2_per_cell = np.array([cell_to_l2.get(c, "__none__") for c in cells_atlas],
                              dtype=object)

    # ---- per-cell group label ----
    if a.color_by == "ng_viable":
        if not a.contrast:
            raise SystemExit("--contrast required for --color-by ng_viable")
        # Per-NG metadata: parent_L2 + size_rank_within_L2
        sum_path = inq / "outputs/stageF1_nhoodgroups" / a.contrast / "nhood_groups_summary.csv"
        ng_summary = pd.read_csv(sum_path, low_memory=False)
        ng_summary["NhoodGroup"] = pd.to_numeric(ng_summary["NhoodGroup"], errors="coerce")
        ng_summary = ng_summary.dropna(subset=["NhoodGroup"])
        ng_summary["NhoodGroup"] = ng_summary["NhoodGroup"].astype(int).astype(str)
        ng_summary = ng_summary[ng_summary.n_nhoods_in_group >= 10]
        # Compute per-L2 rank by descending NG size
        ng_summary = ng_summary.sort_values(
            ["parent_L2_joint", "n_nhoods_in_group"],
            ascending=[True, False])
        ng_summary["rank_in_l2"] = (ng_summary.groupby("parent_L2_joint")
                                                 .cumcount() + 1)
        # Build label: "<parent_label>_<rank>" (compartment-stripped)
        ng_summary["ng_display"] = (ng_summary.parent_label.astype(str) + "_"
                                       + ng_summary.rank_in_l2.astype(str))
        l2s_with_viable_ngs = set(ng_summary.parent_L2_joint.unique())
        print(f"Viable NGs: {len(ng_summary)} across {len(l2s_with_viable_ngs)} L2s")

        # Per-Nhood -> label map
        ng_csv = inq / "outputs/stageF1_nhoodgroups" / a.contrast / "nhood_groups.csv"
        ng_map = pd.read_csv(ng_csv, low_memory=False)
        ng_map["NhoodGroup"] = pd.to_numeric(ng_map["NhoodGroup"], errors="coerce")
        ng_map = ng_map.dropna(subset=["NhoodGroup"])
        ng_map["NhoodGroup"] = ng_map["NhoodGroup"].astype(int).astype(str)
        ng_map = ng_map.merge(
            ng_summary[["NhoodGroup", "parent_L2_joint", "ng_display"]],
            on="NhoodGroup", how="inner")
        nhood_to_label = dict(zip(ng_map.Nhood.astype(int),
                                     ng_map.ng_display.astype(str)))

        print("Loading nhoods.mtx ...")
        nh = sio.mmread(str(nhoods_dir / "nhoods.mtx")).tocsr()
        viable_nhood_ids = sorted([n for n in nhood_to_label.keys()
                                       if 0 < n <= nh.shape[1]])
        col_idx = [n - 1 for n in viable_nhood_ids]
        nlabel_per_col = np.array([nhood_to_label[n] for n in viable_nhood_ids])
        nh_sub = nh[:, col_idx]
        groups = np.array(["__none__"] * nh.shape[0], dtype=object)
        nh_lil = nh_sub.tolil()
        for i in range(nh_sub.shape[0]):
            cols = nh_lil.rows[i]
            if not cols:
                continue
            ngs = nlabel_per_col[cols]
            c = Counter(ngs)
            groups[i] = c.most_common(1)[0][0]

        # Fallback: cells with no dominant viable NG that belong to L2s WITHOUT
        # any viable NG -> show L2 cluster as fallback (parent_label only).
        # Cells in L2s that DO have viable NGs but landed outside -> grey.
        no_ng_mask = (groups == "__none__")
        for i in np.where(no_ng_mask)[0]:
            l2 = l2_per_cell[i]
            if l2 in {"__none__", "__artifact__"}:
                continue
            if l2 not in l2s_with_viable_ngs:
                # Strip compartment prefix; mark with __ prefix to signal L2-fallback
                groups[i] = l2.split("::")[1] + "_L2"
        title = f"NhoodGroup (viable) - {a.contrast}"
    else:
        groups = l2_per_cell.copy()
        title = "L2 cell-type label"

    # ---- Downsample for plotting ----
    rng = np.random.default_rng(42)
    if len(coords) > a.max_points:
        keep_ix = rng.choice(len(coords), a.max_points, replace=False)
    else:
        keep_ix = np.arange(len(coords))
    sub_coords = coords[keep_ix]
    sub_groups = groups[keep_ix]

    # ---- Group medians (for label placement; computed on FULL data, not sub) ----
    unique_groups = [g for g in sorted(set(groups)) if g not in {"__none__", "__artifact__"}]
    print(f"Unique groups: {len(unique_groups)}")
    group_medians = {}
    for g in unique_groups:
        mask = (groups == g)
        if mask.sum() == 0:
            continue
        group_medians[g] = (float(np.median(coords[mask, 0])),
                              float(np.median(coords[mask, 1])))

    # ---- Render ----
    pal = make_palette(len(unique_groups))
    group_colors = dict(zip(unique_groups, pal))
    # Assign integer IDs 1..N for on-plot labeling + side legend
    group_to_num = {g: i + 1 for i, g in enumerate(unique_groups)}

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.2, 1.0], wspace=0.04,
                            left=0.04, right=0.99, top=0.97, bottom=0.04)
    ax = fig.add_subplot(gs[0, 0])
    ax_leg = fig.add_subplot(gs[0, 1])

    # backdrop: __none__ + __artifact__
    backdrop_mask = np.isin(sub_groups, ["__none__", "__artifact__"])
    if backdrop_mask.any():
        ax.scatter(sub_coords[backdrop_mask, 0], sub_coords[backdrop_mask, 1],
                     s=a.point_size, c="#dddddd", alpha=0.4, linewidths=0,
                     rasterized=True)
    # foreground per-group
    for g in unique_groups:
        m = (sub_groups == g)
        if not m.any():
            continue
        ax.scatter(sub_coords[m, 0], sub_coords[m, 1],
                     s=a.point_size, c=group_colors[g], alpha=0.7, linewidths=0,
                     rasterized=True)

    # NUMBERS at median coords (white halo)
    for g, (gx, gy) in group_medians.items():
        ax.text(gx, gy, str(group_to_num[g]),
                  fontsize=a.label_fontsize + 2, color="black", fontweight="bold",
                  ha="center", va="center",
                  path_effects=[matplotlib.patheffects.Stroke(linewidth=1.6, foreground="white"),
                                  matplotlib.patheffects.Normal()])

    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp_name in ("top", "right"):
        ax.spines[sp_name].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.set_xlabel("UMAP_1", fontsize=7)
    ax.set_ylabel("UMAP_2", fontsize=7)

    # ---- Side legend: numbered list with color chips ----
    ax_leg.set_xlim(0, 1); ax_leg.set_ylim(0, 1)
    ax_leg.set_xticks([]); ax_leg.set_yticks([])
    for sp_name in ("top", "right", "bottom", "left"):
        ax_leg.spines[sp_name].set_visible(False)
    n_groups = len(unique_groups)
    # Two columns if many groups; else single column
    n_cols = 2 if n_groups > 30 else 1
    per_col = int(np.ceil(n_groups / n_cols))
    for i, g in enumerate(unique_groups):
        col = i // per_col
        row = i % per_col
        x = 0.02 + col * 0.50
        y = 0.99 - row / max(1, per_col - 1) * 0.97
        # Number
        ax_leg.text(x, y, f"{group_to_num[g]}", fontsize=6, fontweight="bold",
                      ha="left", va="center", transform=ax_leg.transAxes,
                      color="#222")
        # Color chip
        ax_leg.add_patch(plt.Rectangle((x + 0.04, y - 0.005), 0.022, 0.012,
                                         facecolor=group_colors[g],
                                         transform=ax_leg.transAxes,
                                         edgecolor="#666", linewidth=0.2))
        # Label (strip compartment prefix; strip _L2 fallback marker)
        short = g.replace("epi::", "").replace("imm::", "").replace("str::", "")
        if short.endswith("_L2"):
            short = short[:-3] + " (L2)"
        ax_leg.text(x + 0.075, y, short, fontsize=5.5, ha="left", va="center",
                      transform=ax_leg.transAxes, color="#333")

    fig.savefig(a.out_pdf, bbox_inches="tight", dpi=300)
    print(f"Wrote {a.out_pdf}")


if __name__ == "__main__":
    import matplotlib.patheffects
    main()
