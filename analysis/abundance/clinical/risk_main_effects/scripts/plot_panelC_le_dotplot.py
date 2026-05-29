#!/usr/bin/env python3
"""
plot_panelC_le_dotplot.py - Stage-2 dotplot from leading-edge features.

Reads outputs/panelC_features/<L2_safe>_<contrast>/. Per-NG aggregation:
  size  = fraction of cells in NG with count > 0     (pct expressing)
  color = z-score of mean(log1p(CPM)) across NGs     (relative expression)

Cols = NhoodGroups (ordered by within-L2 size rank).
Rows = genes (ordered by source-NG rank, then score).
"""
from __future__ import annotations

import argparse
import gzip
import logging
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("panelC_le_dotplot")

RANK_PALETTE = [
    "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
    "#CC6677", "#AA4499", "#882255", "#999933", "#661100",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--features-dir", required=True, type=Path)
    p.add_argument("--out-pdf", default=None)
    p.add_argument("--max-dot-area", type=float, default=120.0)
    return p.parse_args()


def main():
    a = parse_args()
    fd = a.features_dir
    cells = pd.read_csv(fd / "cells.tsv", sep="\t")
    genes = pd.read_csv(fd / "genes.tsv", sep="\t")
    ng_summary = pd.read_csv(fd / "ng_summary.csv")
    log.info("cells=%d, genes=%d, NGs=%d", len(cells), len(genes), len(ng_summary))

    cells["NG_dominant"] = cells["NG_dominant"].astype(str)
    ng_summary["NhoodGroup"] = ng_summary["NhoodGroup"].astype(str)
    ng_summary = ng_summary.sort_values("rank").reset_index(drop=True)

    # Load counts + log1p(CPM)
    log.info("Loading counts.mtx.gz")
    with gzip.open(fd / "counts.mtx.gz", "rb") as f:
        m = sio.mmread(f)
    X = sp.csr_matrix(m).astype(float)
    libsize = np.asarray(X.sum(axis=1)).ravel()
    libsize[libsize == 0] = 1.0
    # Per-cell normalize then log1p - keep sparse-friendly: dense for the small subset
    Xd = X.toarray()
    Xn = np.log1p(Xd / libsize[:, None] * 1e4)  # cells x genes
    log.info("normalized matrix: %s", Xn.shape)

    # Aggregate per NG: pct expr (count > 0), mean log1p-CPM
    ng_order = ng_summary["NhoodGroup"].tolist()
    pct = np.zeros((len(ng_order), Xn.shape[1]))
    mean_expr = np.zeros((len(ng_order), Xn.shape[1]))
    for i, ng in enumerate(ng_order):
        rows = (cells["NG_dominant"] == ng).values
        n = rows.sum()
        if n == 0:
            continue
        pct[i, :] = (Xd[rows, :] > 0).mean(axis=0)
        mean_expr[i, :] = Xn[rows, :].mean(axis=0)

    # Z-score mean_expr per gene across NGs, clip
    mu = mean_expr.mean(axis=0, keepdims=True)
    sd = mean_expr.std(axis=0, keepdims=True) + 1e-9
    z = (mean_expr - mu) / sd
    z = np.clip(z, -2.5, 2.5)

    # Row order = source-NG rank then score
    genes["source_rank"] = genes["source_rank"].astype(int)
    genes["score"] = genes["score"].astype(float)
    row_order = genes.sort_values(["source_rank", "score"],
                                    ascending=[True, False]).index.values
    z = z[:, row_order]
    pct = pct[:, row_order]
    genes = genes.iloc[row_order].reset_index(drop=True)

    # ---- Render ----
    n_genes = len(genes)
    n_ngs = len(ng_order)
    fig_h = max(4.5, 0.18 * n_genes + 1.8)
    fig_w = max(4.5, 0.5 * n_ngs + 3.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.5, n_genes],
                           hspace=0.04)
    ax_dir = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=ax_dir)

    # Direction strip: per-NG group_med_lfc, red(+)/blue(-) diverging
    med_lfcs = ng_summary["group_med_lfc"].values.astype(float)
    lfc_lim = max(abs(med_lfcs).max(), 0.5)
    div_cmap = LinearSegmentedColormap.from_list(
        "div", ["#1f77b4", "#f7f7f7", "#d62728"])
    ax_dir.imshow(med_lfcs.reshape(1, -1), aspect="auto",
                   cmap=div_cmap, vmin=-lfc_lim, vmax=lfc_lim,
                   interpolation="nearest")
    for i, v in enumerate(med_lfcs):
        ax_dir.text(i, 0, f"{v:+.2f}", ha="center", va="center",
                     fontsize=6, fontweight="bold",
                     color="white" if abs(v) > lfc_lim * 0.5 else "#222")
    ax_dir.set_yticks([0])
    ax_dir.set_yticklabels(["med logFC"], fontsize=6)
    ax_dir.set_xticks([])
    for sp_name in ("top", "right", "left", "bottom"):
        ax_dir.spines[sp_name].set_visible(False)

    # Build grid coords (cols=NGs, rows=genes)
    xs, ys = np.meshgrid(np.arange(n_ngs), np.arange(n_genes), indexing="xy")
    sizes = pct.T * a.max_dot_area  # transpose -> (genes, NGs); pct is (NGs, genes)
    sizes = np.where(sizes < 1.0, 0.0, sizes)  # hide near-zero dots

    cmap = LinearSegmentedColormap.from_list("rb", ["#2166AC", "white", "#B2182B"])
    sc = ax.scatter(xs.ravel(), ys.ravel(),
                     s=sizes.ravel(),
                     c=z.T.ravel(), cmap=cmap, vmin=-2.5, vmax=2.5,
                     edgecolors="#666666", linewidths=0.2)

    # Axes
    ax.set_xticks(np.arange(n_ngs))
    # Use NG name_short ("BMYO-basal_1") and color tick by rank
    short_labels = []
    rank_to_color = {i + 1: c for i, c in enumerate(RANK_PALETTE)}
    for _, r in ng_summary.iterrows():
        nf = str(r["name_full"])
        nf_parts = nf.split("_")
        rank_part = nf_parts[-2] if len(nf_parts) >= 2 else str(r["rank"])
        short_labels.append(f"{rank_part}")
    ax.set_xticklabels(short_labels, fontsize=6, rotation=0)
    for tick, rk in zip(ax.get_xticklabels(), ng_summary["rank"]):
        tick.set_color(rank_to_color.get(int(rk), "#444"))
        tick.set_fontweight("bold")

    ax.set_yticks(np.arange(n_genes))
    ax.set_yticklabels(genes["symbol"].tolist(), fontsize=5.5)

    # Row dividers between source-NG groups
    row_rank = genes["source_rank"].values
    row_breaks = np.where(np.diff(row_rank) != 0)[0] + 0.5
    for b in row_breaks:
        ax.axhline(b, color="black", linewidth=0.4, alpha=0.25)

    ax.set_xlim(-0.5, n_ngs - 0.5)
    ax.set_ylim(n_genes - 0.5, -0.5)
    ax.set_xlabel("NhoodGroup (rank within L2)", fontsize=7)
    ax.tick_params(axis="both", which="both", length=0)
    for sp_name in ("top", "right"):
        ax.spines[sp_name].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.spines["bottom"].set_linewidth(0.4)

    # Color bar
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("z(mean log1p CPM)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    # Size legend (separate axes)
    fig.subplots_adjust(right=0.82)
    legend_pcts = [0.1, 0.25, 0.5, 0.75, 1.0]
    handles = [plt.scatter([], [], s=p * a.max_dot_area, color="#888888",
                            edgecolor="#666666", linewidth=0.2,
                            label=f"{int(p*100)}%") for p in legend_pcts]
    ax.legend(handles=handles, title="% cells\nexpressing",
               bbox_to_anchor=(1.18, 0.5), loc="center left",
               frameon=False, fontsize=6, title_fontsize=6, labelspacing=0.8)

    out = Path(a.out_pdf) if a.out_pdf else (
        fd / f"panelC_le_dotplot_{fd.name}.pdf")
    fig.savefig(out, bbox_inches="tight")
    log.info("Wrote %s", out)


if __name__ == "__main__":
    main()
