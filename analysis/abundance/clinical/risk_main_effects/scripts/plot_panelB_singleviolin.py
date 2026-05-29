#!/usr/bin/env python3
"""plot_panelB_singleviolin.py - single-violin per-L2 forest for main-effect contrasts.

Single-side half-violin (hollow grey) per L2 + NhoodGroup forest dots
(median + IQR, rank-colored, sized by n_nhoods). For risk_main_effects
contrasts (BR1/HRS/BR2 vs AR, HRS_no_cancerhx) where only one coef is
tested per contrast (no parity x risk interaction).

Mirrors plot_panelB_split_forest.R aesthetic, drops the split.

Args:
  --inquiry-dir
  --contrast
  --labels-full (default: publication labels_full.csv)
  --out-pdf
"""
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.load_paths import load_inquiry  # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

RANK_PALETTE = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
                "#CC6677", "#AA4499", "#882255", "#999933", "#661100"]

COMPARTMENT_ORDER = ["epi", "imm", "str"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--contrast", required=True)
    p.add_argument("--labels-full", default=None,
                   help="Override path; default resolves via paths.yaml inputs.labels")
    p.add_argument("--out-pdf", required=True)
    p.add_argument("--width-in", type=float, default=12.0)
    p.add_argument("--height-in", type=float, default=4.5)
    p.add_argument("--violin-width", type=float, default=0.40,
                   help="half-violin max width (in L2-grid units)")
    p.add_argument("--max-dot-area", type=float, default=70.0)
    return p.parse_args()


def load_data(inq, contrast, labels_full_path):
    """Load da_results + nhood_groups; filter artifacts via labels_full."""
    da = pd.read_csv(inq / "outputs/stageD_da_results" / contrast / "da_results.csv",
                       low_memory=False)
    ng_map = pd.read_csv(inq / "outputs/stageF1_nhoodgroups" / contrast / "nhood_groups.csv",
                           low_memory=False)
    # Build L2_joint if missing (risk inquiry da_results lacks it; pxr has it)
    if "L2_joint" not in da.columns:
        da["L2_joint"] = (da["compartment"].astype(str) + "::"
                            + da["label"].astype(str))
    # Keep canonical fields
    da = da[["Nhood", "logFC", "SpatialFDR", "compartment", "label", "L2_joint",
              "is_artifact"]].copy()
    da["is_artifact"] = da["is_artifact"].astype(str).str.lower().isin({"true", "1"})

    # Artifact L2 filter (drop L2s with >=50% is_artifact in labels_full)
    labels = pd.read_csv(labels_full_path,
                          usecols=["compartment", "label", "is_artifact"])
    labels["is_artifact"] = labels["is_artifact"].astype(str).str.lower().isin(
        {"true", "1"})
    artifact_pct = (labels.groupby(["compartment", "label"])["is_artifact"]
                          .mean().reset_index())
    artifact_l2s = set(artifact_pct[artifact_pct.is_artifact >= 0.50]
                        .apply(lambda r: f"{r['compartment']}::{r['label']}", axis=1))
    da = da[~da.L2_joint.isin(artifact_l2s)].copy()

    # NG mapping (drop NA-NhoodGroup rows; cast to int-string)
    ng_map["NhoodGroup"] = pd.to_numeric(ng_map["NhoodGroup"], errors="coerce")
    ng_map = ng_map.dropna(subset=["NhoodGroup"])
    ng_map["NhoodGroup"] = ng_map["NhoodGroup"].astype(int).astype(str)

    return da, ng_map


def per_l2_aggregate(da, ng_map):
    """Per-L2: collect nhood logFCs (for violin); per-NG: aggregate (for forest dots)."""
    # Filter da to non-artifact rows already done
    # Per-L2 nhood logFCs
    by_l2 = (da.groupby(["compartment", "L2_joint"])
                 .agg(logfcs=("logFC", list),
                      n_nhoods=("Nhood", "count"))
                 .reset_index())

    # Per-NG aggregates: median, q25, q75, n_nhoods, parent_L2
    # ng_map already carries per-nhood logFC (annotated by F.1); no merge needed.
    per_ng = (ng_map.groupby(["parent_compartment", "parent_L2_joint",
                                  "NhoodGroup"], as_index=False)
                       .agg(med=("logFC", "median"),
                            q25=("logFC", lambda x: np.percentile(x, 25)),
                            q75=("logFC", lambda x: np.percentile(x, 75)),
                            n=("logFC", "count")))
    # Viability: n_nhoods >= 10
    per_ng = per_ng[per_ng.n >= 10]
    # Rank within parent_L2 by descending n
    per_ng = per_ng.sort_values(["parent_L2_joint", "n"],
                                  ascending=[True, False]).reset_index(drop=True)
    per_ng["rank_in_l2"] = per_ng.groupby("parent_L2_joint").cumcount() + 1
    return by_l2, per_ng


def order_l2s_per_compartment(by_l2):
    """Order L2s by median nhood logFC desc within compartment."""
    by_l2 = by_l2.copy()
    by_l2["med_lfc"] = by_l2["logfcs"].apply(lambda v: float(np.median(v)) if len(v) else 0.0)
    out = []
    for comp in COMPARTMENT_ORDER:
        sub = by_l2[by_l2.compartment == comp].sort_values("med_lfc", ascending=False)
        out.append(sub)
    return pd.concat(out, ignore_index=True)


def main():
    a = parse_args()
    if not a.labels_full:
        a.labels_full = load_inquiry(a.inquiry_dir)["paths"]["inputs"]["labels"]
    da, ng_map = load_data(a.inquiry_dir, a.contrast, a.labels_full)
    by_l2, per_ng = per_l2_aggregate(da, ng_map)
    ordered = order_l2s_per_compartment(by_l2)
    print(f"L2s after artifact filter: {len(ordered)}, "
          f"viable NGs: {len(per_ng)}")

    # Globally unique x positions across all L2s
    ordered["x_pos"] = np.arange(len(ordered))
    l2_to_x = dict(zip(ordered.L2_joint, ordered.x_pos))

    # Compartment facet ranges
    comp_ranges = {}
    for comp in COMPARTMENT_ORDER:
        sub = ordered[ordered.compartment == comp]
        if len(sub):
            comp_ranges[comp] = (int(sub.x_pos.min()), int(sub.x_pos.max()))

    # ---- Figure with one row, faceted by compartment (free_x) ----
    fig = plt.figure(figsize=(a.width_in, a.height_in))
    width_ratios = [comp_ranges[c][1] - comp_ranges[c][0] + 1
                       for c in COMPARTMENT_ORDER if c in comp_ranges]
    gs = gridspec.GridSpec(1, len(width_ratios), figure=fig,
                             width_ratios=width_ratios,
                             wspace=0.04,
                             left=0.05, right=0.99, top=0.92, bottom=0.20)

    rank_to_color = {i + 1: c for i, c in enumerate(RANK_PALETTE)}

    # Y-axis range (shared across facets)
    all_lfcs = np.concatenate([np.asarray(v) for v in by_l2["logfcs"]])
    y_min = np.percentile(all_lfcs, 0.5)
    y_max = np.percentile(all_lfcs, 99.5)
    y_range = y_max - y_min

    # Plot each compartment facet
    facet_idx = 0
    for comp in COMPARTMENT_ORDER:
        if comp not in comp_ranges:
            continue
        ax = fig.add_subplot(gs[0, facet_idx])
        facet_idx += 1
        x0, x1 = comp_ranges[comp]
        sub = ordered[ordered.compartment == comp]
        ax.set_xlim(x0 - 0.5, x1 + 0.5)
        ax.set_ylim(y_min - 0.05 * y_range, y_max + 0.05 * y_range)
        ax.axhline(0, color="grey", linestyle="--", linewidth=0.4, alpha=0.6)

        # Per-L2 half violin + forest dots
        for _, row in sub.iterrows():
            x = row.x_pos
            lfcs = np.asarray(row.logfcs)
            if len(lfcs) >= 5 and lfcs.std() > 1e-6:
                kde = gaussian_kde(lfcs)
                y_grid = np.linspace(y_min, y_max, 200)
                d = kde(y_grid)
                # Normalize to violin_width
                d_max = d.max() if d.max() > 0 else 1.0
                d_scaled = d / d_max * a.violin_width
                # Full bilateral violin: mirror around x
                xs = np.r_[(x - d_scaled), (x + d_scaled)[::-1]]
                ys = np.r_[y_grid, y_grid[::-1]]
                ax.add_patch(Polygon(np.c_[xs, ys], closed=True,
                                       facecolor="none",
                                       edgecolor="#555555",
                                       linewidth=0.35))

            # Forest dots: per NG belonging to this L2
            ngs = per_ng[per_ng.parent_L2_joint == row.L2_joint]
            for _, ng in ngs.iterrows():
                rk = int(ng.rank_in_l2)
                col = rank_to_color.get(rk, "#888")
                size = max(8, min(a.max_dot_area, np.sqrt(ng.n) * 4.5))
                # IQR vertical line (q25 - q75), forest dot at median
                ax.plot([x, x], [ng.q25, ng.q75], color=col, linewidth=0.7,
                          alpha=0.85, solid_capstyle="butt")
                ax.scatter([x], [ng.med], s=size, color=col,
                            edgecolor="white", linewidth=0.4, zorder=5)

        # X-axis: L2 names rotated 90°
        ax.set_xticks(sub.x_pos.values)
        ax.set_xticklabels([l.split("::")[1] for l in sub.L2_joint.values],
                            fontsize=5.5, rotation=90, ha="center")
        ax.tick_params(axis="x", which="both", length=0, pad=1.5)
        ax.tick_params(axis="y", labelsize=6, length=2)

        # Compartment title (top)
        ax.text(0.5, 0.99, comp, transform=ax.transAxes,
                  fontsize=8, fontweight="bold", ha="center", va="top",
                  color="#222")
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.spines["left"].set_linewidth(0.4)
        ax.spines["bottom"].set_linewidth(0.4)

        if facet_idx == 1:
            ax.set_ylabel(f"per-nhood logFC ({a.contrast})", fontsize=7)

    # No big title, no on-plot caption (Illustrator-clean)
    fig.savefig(a.out_pdf, bbox_inches="tight", dpi=300)
    print(f"Wrote {a.out_pdf}")


if __name__ == "__main__":
    main()
