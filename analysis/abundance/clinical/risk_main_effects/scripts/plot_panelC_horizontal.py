#!/usr/bin/env python3
"""
plot_panelC_horizontal.py - Panel C horizontal dotplot, multi-L2 stacked.

Sharp publication-style: NGs as rows (left to right: rank chip + direction
chip + label), genes as cols. One strip per L2, stacked vertically with
shared legend on the right.

  size  = fraction of cells in NG with count > 0
  color = z(mean log1p CPM) per gene across NGs in that L2

Args:
  --inquiry-dir
  --L2-list  comma-separated L2s
  --contrast (default parity_x_HR_BRCA1)
  --genes-per-l2 (default 25 - top by score per L2)
  --out-pdf
"""
from __future__ import annotations
import argparse
import gzip
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

RANK_PALETTE = ["#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
                "#CC6677", "#AA4499", "#882255", "#999933", "#661100"]
DIV_CMAP = LinearSegmentedColormap.from_list(
    "div", ["#1f77b4", "#f7f7f7", "#d62728"])
EXPR_CMAP = LinearSegmentedColormap.from_list(
    "expr", ["#2166AC", "white", "#B2182B"])


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--L2-list", required=True,
                   help="comma-separated L2 names")
    p.add_argument("--contrast", default="parity_x_HR_BRCA1")
    p.add_argument("--genes-per-ng", type=int, default=5,
                   help="genes per source-NG within an L2 strip (balanced)")
    p.add_argument("--out-pdf", required=True)
    p.add_argument("--max-dot-area", type=float, default=70.0)
    p.add_argument("--per-ng-height-in", type=float, default=0.30,
                   help="absolute inches per NG row in the dot grid")
    p.add_argument("--strip-overhead-in", type=float, default=0.55,
                   help="inches reserved for L2 title + gene labels per strip")
    p.add_argument("--legend-height-in", type=float, default=0.55)
    p.add_argument("--width-in", type=float, default=10.0,
                   help="Figure width (Nature double-col=7.2; wider for breathing room)")
    p.add_argument("--highlight-genes", default="",
                   help="comma-separated symbols to spike-in/highlight (always shown, "
                         "labels rendered in red bold)")
    return p.parse_args()


def load_features(inq, L2, contrast):
    L2_safe = L2.replace("::", "__").replace("-", "_")
    fd = inq / "outputs/panelC_features" / f"{L2_safe}_{contrast}"
    cells = pd.read_csv(fd / "cells.tsv", sep="\t")
    genes = pd.read_csv(fd / "genes.tsv", sep="\t")
    ng_summary = pd.read_csv(fd / "ng_summary.csv")
    with gzip.open(fd / "counts.mtx.gz", "rb") as f:
        m = sio.mmread(f)
    X = sp.csr_matrix(m).astype(float)
    return dict(L2=L2, cells=cells, genes=genes, ng_summary=ng_summary, X=X)


def select_genes(d, genes_per_ng, total_target, highlight_set=None):
    """Per-L2 gene selection.

    Stage 1: balanced top-K per source-NG (diagonal block pattern).
    Stage 2: pad to total_target by next-highest-score across NGs.
    Stage 3: spike-in highlight genes (rescue if not already selected).
    """
    g = d["genes"].copy()
    g["source_rank"] = g["source_rank"].astype(int)
    g["score"] = g["score"].astype(float)
    balanced = (g.sort_values("score", ascending=False)
                  .groupby("source_rank", sort=True)
                  .head(genes_per_ng))
    remaining = g.loc[~g.index.isin(balanced.index)].sort_values(
        "score", ascending=False)
    need = max(0, total_target - len(balanced))
    extra = remaining.head(need)
    out = pd.concat([balanced, extra], ignore_index=True)

    # Spike-in: if any highlight gene isn't already in `out`, rescue from full
    # genes table. Use the gene's source_NG (whatever the extractor recorded)
    # so it lives in the right block.
    if highlight_set:
        in_panel = set(out["symbol"])
        missing = [s for s in highlight_set if s not in in_panel]
        if missing:
            rescue = g[g["symbol"].isin(missing)]
            if len(rescue):
                out = pd.concat([out, rescue], ignore_index=True)
                print(f"  spike-in rescued: {sorted(set(rescue['symbol']))}")

    out = out.sort_values(["source_rank", "score"],
                            ascending=[True, False]).reset_index(drop=True)
    out = out.drop_duplicates("symbol", keep="first").reset_index(drop=True)
    return out


def aggregate(d, picked_genes):
    ng_summary = d["ng_summary"].sort_values("rank").reset_index(drop=True)
    ng_summary["NhoodGroup"] = ng_summary["NhoodGroup"].astype(str)
    g = picked_genes

    cells = d["cells"].copy()
    cells["NG_dominant"] = cells["NG_dominant"].astype(str)

    X = d["X"]
    Xsub = X[:, g["row_idx"].astype(int).values].toarray()
    libsize = X.sum(axis=1).A1
    libsize[libsize == 0] = 1.0
    Xn = np.log1p(Xsub / libsize[:, None] * 1e4)

    n_ngs = len(ng_summary)
    n_genes = len(g)
    pct = np.zeros((n_ngs, n_genes))
    mean_expr = np.zeros((n_ngs, n_genes))
    for i, ng in enumerate(ng_summary["NhoodGroup"]):
        mask = (cells["NG_dominant"] == ng).values
        if mask.sum() == 0:
            continue
        pct[i, :] = (Xsub[mask, :] > 0).mean(axis=0)
        mean_expr[i, :] = Xn[mask, :].mean(axis=0)

    mu = mean_expr.mean(axis=0, keepdims=True)
    sd = mean_expr.std(axis=0, keepdims=True) + 1e-9
    z = np.clip((mean_expr - mu) / sd, -2.5, 2.5)

    return dict(
        L2=d["L2"], pct=pct, z=z,
        gene_labels=g["symbol"].tolist(),
        source_rank_per_gene=g["source_rank"].tolist(),
        ng_order=ng_summary["NhoodGroup"].tolist(),
        ng_rank=ng_summary["rank"].astype(int).tolist(),
        ng_med_lfc=ng_summary["group_med_lfc"].astype(float).tolist(),
        ng_n_cells=ng_summary["n_cells"].astype(int).tolist())


def draw_strip(fig, gs_cell, agg, max_dot_area):
    """One L2 strip: header, gene labels (top), NG rows (left chips + dots)."""
    # Two-column inner: left=labels (NG name + chips), right=dotgrid+gene-labels
    inner = gs_cell.subgridspec(2, 2,
                                  height_ratios=[0.28, 1.0],
                                  width_ratios=[0.18, 1.0],
                                  hspace=0.0, wspace=0.012)
    ax_title = fig.add_subplot(inner[0, :])
    ax_lbl = fig.add_subplot(inner[1, 0])
    ax_dots = fig.add_subplot(inner[1, 1])

    n_ngs = len(agg["ng_order"])
    n_genes = len(agg["gene_labels"])

    # ---- Title row: L2 + gene labels (top) ----
    ax_title.set_xlim(0, 1); ax_title.set_ylim(0, 1)
    ax_title.set_xticks([]); ax_title.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax_title.spines[sp].set_visible(False)
    L2_short = agg["L2"]
    ax_title.text(0.0, 0.05, L2_short, fontsize=8, fontweight="bold",
                    ha="left", va="bottom", color="#222")

    # Gene labels above the dot grid (right side of title row)
    # Use a sibling axis aligned with dots
    ax_genes = ax_dots.twiny()
    ax_genes.set_xlim(-0.5, n_genes - 0.5)
    ax_genes.set_xticks(np.arange(n_genes))
    ax_genes.set_xticklabels(agg["gene_labels"],
                                fontsize=4.6, rotation=90, ha="center", va="bottom")
    ax_genes.tick_params(axis="x", which="both", length=0, pad=1, top=True,
                            labeltop=True, bottom=False, labelbottom=False)
    # (highlight styling intentionally disabled - spike-in still rescues genes
    # into the panel via select_genes; visual call-out left to coordinator/PI)
    for sp in ("top", "right", "bottom", "left"):
        ax_genes.spines[sp].set_visible(False)

    # ---- Left labels: rank in rank-color + (logfc) ----
    rank_to_color = {i + 1: c for i, c in enumerate(RANK_PALETTE)}
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(-0.5, n_ngs - 0.5)
    ax_lbl.invert_yaxis()
    ax_lbl.set_xticks([]); ax_lbl.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax_lbl.spines[sp].set_visible(False)
    for i, (rk, m) in enumerate(zip(agg["ng_rank"], agg["ng_med_lfc"])):
        ax_lbl.text(0.95, i,
                      f"{rk} ({m:+.2f})",
                      fontsize=6, fontweight="bold",
                      ha="right", va="center",
                      color="#222")

    # ---- Dot grid ----
    xs, ys = np.meshgrid(np.arange(n_genes), np.arange(n_ngs), indexing="xy")
    sizes = agg["pct"] * max_dot_area
    sizes = np.where(sizes < 1.5, 0.0, sizes)
    sc = ax_dots.scatter(xs.ravel(), ys.ravel(),
                            s=sizes.ravel(),
                            c=agg["z"].ravel(), cmap=EXPR_CMAP, vmin=-2.5, vmax=2.5,
                            edgecolors="#666", linewidths=0.18)

    ax_dots.set_xlim(-0.5, n_genes - 0.5)
    ax_dots.set_ylim(n_ngs - 0.5, -0.5)
    ax_dots.set_xticks([])
    ax_dots.set_yticks([])
    ax_dots.tick_params(axis="both", which="both", length=0)
    for sp in ("top", "right"):
        ax_dots.spines[sp].set_visible(False)
    ax_dots.spines["left"].set_linewidth(0.3)
    ax_dots.spines["bottom"].set_linewidth(0.3)

    # Gene block dividers (between source-NG groups)
    src_rank = np.array(agg["source_rank_per_gene"])
    breaks = np.where(np.diff(src_rank) != 0)[0] + 0.5
    for b in breaks:
        ax_dots.axvline(b, color="#bbb", linewidth=0.4, alpha=0.7,
                          linestyle=(0, (1, 2)))
    return sc


def main():
    a = parse_args()
    L2_list = [s.strip() for s in a.L2_list.split(",")]
    highlight_set = set(s.strip() for s in a.highlight_genes.split(",") if s.strip())
    # Two-pass: load features, find target gene count, then aggregate uniform
    raw = [load_features(a.inquiry_dir, L2, a.contrast) for L2 in L2_list]
    balanced_sizes = []
    for d in raw:
        g = d["genes"].copy()
        g["source_rank"] = g["source_rank"].astype(int)
        balanced = (g.groupby("source_rank").head(a.genes_per_ng))
        balanced_sizes.append(len(balanced))
    target_n = max(balanced_sizes)
    print(f"per-L2 balanced sizes: {balanced_sizes}; target n cols: {target_n}")
    picks = [select_genes(d, a.genes_per_ng, target_n, highlight_set) for d in raw]
    aggs = [aggregate(d, p) for d, p in zip(raw, picks)]
    # Carry highlight set into rendering
    for ag in aggs:
        ag["highlight_set"] = highlight_set

    n_strips = len(L2_list)
    # Absolute per-strip heights (LHS-major gets ~3in for 10 NGs).
    n_ngs_per = [len(agg["ng_order"]) for agg in aggs]
    strip_heights_in = [n * a.per_ng_height_in + a.strip_overhead_in
                          for n in n_ngs_per]
    total_strip_h = sum(strip_heights_in)
    fig_h = total_strip_h + 0.4
    fig = plt.figure(figsize=(a.width_in, fig_h))
    print(f"per-strip heights (in): {[round(h,2) for h in strip_heights_in]}, fig_h={fig_h:.2f}")

    # Strips on the left; right side reserved for stacked legends (no bottom row)
    outer = gridspec.GridSpec(n_strips, 2, figure=fig,
                                height_ratios=strip_heights_in,
                                width_ratios=[1.0, 0.13],
                                hspace=0.0, wspace=0.04,
                                left=0.005, right=0.99,
                                top=0.97, bottom=0.03)
    last_sc = None
    for i, agg in enumerate(aggs):
        sc = draw_strip(fig, outer[i, 0], agg, a.max_dot_area)
        last_sc = sc

    # ------ Stacked legends on the right side ------
    # Top: colorbar (z-CPM)
    cbar_ax = fig.add_axes((0.91, 0.62, 0.014, 0.18))
    cbar = fig.colorbar(last_sc, cax=cbar_ax)
    cbar.set_label("z(mean\nlog1p CPM)", fontsize=6, rotation=0,
                    labelpad=14, loc="center")
    cbar.ax.tick_params(labelsize=5, length=2, pad=1.5)
    cbar.outline.set_linewidth(0.3)

    # Below colorbar: size legend (vertical stack of dots)
    size_ax = fig.add_axes((0.90, 0.20, 0.08, 0.32))
    size_ax.set_xlim(0, 1); size_ax.set_ylim(0, 1)
    size_ax.set_xticks([]); size_ax.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        size_ax.spines[sp].set_visible(False)

    legend_dot_area = a.max_dot_area * 2.8
    size_ax.text(0.5, 0.96, "% cells\nexpressing", fontsize=6, ha="center",
                   va="top", transform=size_ax.transAxes)
    pcts = [1.00, 0.80, 0.50, 0.25, 0.10]
    y_top, y_bot = 0.78, 0.05
    for i, p in enumerate(pcts):
        y = y_top - i * (y_top - y_bot) / (len(pcts) - 1)
        size_ax.scatter([0.30], [y], s=p * legend_dot_area, color="#888",
                          edgecolor="#666", linewidth=0.25,
                          transform=size_ax.transAxes)
        size_ax.text(0.62, y, f"{int(p*100)}%", fontsize=6,
                       ha="left", va="center", transform=size_ax.transAxes)

    fig.savefig(a.out_pdf, bbox_inches="tight", dpi=300)
    print(f"Wrote {a.out_pdf}")


if __name__ == "__main__":
    main()
