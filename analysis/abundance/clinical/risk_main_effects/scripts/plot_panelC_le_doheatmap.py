#!/usr/bin/env python3
"""
plot_panelC_le_doheatmap.py - cell-grain DoHeatmap (Stage-2 plot).

Reads the feature artifact at outputs/panelC_features/<L2_safe>_<contrast>/
written by extract_panelC_features.py. No counts.npz / milo loading here -
plot iterations are seconds, not minutes.

Per-cell log1p(CPM) heatmap with cells as columns grouped by NhoodGroup,
genes as rows ordered by source-NG rank. Top annotation strip = NG identity
colored by within-L2 size rank.

Args:
  --features-dir   path to outputs/panelC_features/<L2_safe>_<contrast>/
  --max-cells-per-ng   downsample cap per NG (default 200)
  --out-pdf        output PDF path (default: same dir / panelC_<safe>_<contrast>_doheatmap.pdf)
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
log = logging.getLogger("panelC_le_doheatmap")


RANK_PALETTE = [
    "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
    "#CC6677", "#AA4499", "#882255", "#999933", "#661100",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--features-dir", required=True, type=Path)
    p.add_argument("--max-cells-per-ng", type=int, default=200)
    p.add_argument("--out-pdf", default=None)
    return p.parse_args()


def load_counts(features_dir: Path) -> sp.csr_matrix:
    mtx_path = features_dir / "counts.mtx.gz"
    log.info("Loading %s", mtx_path)
    with gzip.open(mtx_path, "rb") as f:
        m = sio.mmread(f)
    return sp.csr_matrix(m)


def main():
    a = parse_args()
    fd = a.features_dir
    if not fd.exists():
        raise SystemExit(f"features dir not found: {fd}")

    cells = pd.read_csv(fd / "cells.tsv", sep="\t")
    genes = pd.read_csv(fd / "genes.tsv", sep="\t")
    ng_summary = pd.read_csv(fd / "ng_summary.csv")
    log.info("cells=%d, genes=%d, NGs=%d", len(cells), len(genes), len(ng_summary))

    # NhoodGroup -> rank lookup
    ng_summary["NhoodGroup"] = ng_summary["NhoodGroup"].astype(str)
    ng_rank = dict(zip(ng_summary.NhoodGroup, ng_summary["rank"]))
    ng_name_full = dict(zip(ng_summary.NhoodGroup, ng_summary.name_full))
    ng_med_lfc = dict(zip(ng_summary.NhoodGroup, ng_summary.group_med_lfc))

    cells["NG_dominant"] = cells["NG_dominant"].astype(str)
    cells["NG_rank"] = cells["NG_dominant"].map(ng_rank)

    # Downsample per NG
    rng = np.random.default_rng(42)
    keep_rows = []
    for ng in ng_summary.NhoodGroup:
        ix = cells.index[cells.NG_dominant == ng].to_numpy()
        if len(ix) > a.max_cells_per_ng:
            ix = rng.choice(ix, a.max_cells_per_ng, replace=False)
        keep_rows.extend(ix.tolist())
    keep_rows = np.array(keep_rows, dtype=int)
    cells_sub = cells.iloc[keep_rows].reset_index(drop=True)
    log.info("cells after downsample: %d", len(cells_sub))

    # Order columns: by NG rank, stable within NG
    col_order = np.argsort(cells_sub["NG_rank"].values, kind="stable")
    cells_sub = cells_sub.iloc[col_order].reset_index(drop=True)
    keep_rows = keep_rows[col_order]

    # Load counts (cells x genes)
    X = load_counts(fd)
    log.info("counts: %s, nnz=%d", X.shape, X.nnz)
    sub = X[keep_rows, :].tocsr()

    # log1p(CPM-like) per-cell normalization
    libsize = np.asarray(sub.sum(axis=1)).ravel().astype(float)
    libsize[libsize == 0] = 1.0
    sub_norm = np.log1p((sub.toarray().T / libsize * 1e4).T)

    # z-score per gene; clip to [-3,3]
    z = (sub_norm - sub_norm.mean(axis=0)) / (sub_norm.std(axis=0) + 1e-9)
    z = np.clip(z.T, -3, 3)  # rows = genes, cols = cells

    # Order rows by source-NG rank (block-grouped genes), then by gene_id
    genes["source_rank"] = genes["source_rank"].astype(int)
    row_order = np.argsort(genes["source_rank"].values, kind="stable")
    z = z[row_order, :]
    genes_sub = genes.iloc[row_order].reset_index(drop=True)

    # ---- Render ----
    fig_h = max(4.5, 0.13 * z.shape[0] + 1.5)
    fig_w = max(7.5, 0.012 * z.shape[1] + 4.0)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(3, 2, height_ratios=[0.5, 0.5, 20],
                          width_ratios=[40, 1],
                          hspace=0.05, wspace=0.05)
    ax_dir = fig.add_subplot(gs[0, 0])
    ax_ann = fig.add_subplot(gs[1, 0], sharex=ax_dir)
    ax_hm = fig.add_subplot(gs[2, 0], sharex=ax_dir)
    ax_cb = fig.add_subplot(gs[2, 1])

    # Direction strip: per-cell column = parent NG's group_med_lfc, red(+)/blue(-)
    rank_to_color = {i + 1: c for i, c in enumerate(RANK_PALETTE)}
    col_rank = cells_sub["NG_rank"].values
    col_lfc = cells_sub["NG_dominant"].map(ng_med_lfc).values.astype(float)
    lfc_lim = max(np.nanmax(np.abs(col_lfc)), 0.5)
    div_cmap = LinearSegmentedColormap.from_list(
        "div", ["#1f77b4", "#f7f7f7", "#d62728"])
    ax_dir.imshow(col_lfc.reshape(1, -1), aspect="auto",
                   cmap=div_cmap, vmin=-lfc_lim, vmax=lfc_lim,
                   interpolation="nearest")
    ax_dir.set_yticks([0]); ax_dir.set_xticks([])
    ax_dir.set_yticklabels(["med\nlogFC"], fontsize=6)

    # Rank annotation strip (column-side)
    ann_colors = np.array([rank_to_color.get(int(r), "#888888") for r in col_rank])
    rgb = [matplotlib.colors.to_rgb(c) for c in ann_colors]
    ax_ann.imshow(np.array(rgb).reshape(1, -1, 3), aspect="auto",
                   interpolation="nearest")
    ax_ann.set_yticks([0]); ax_ann.set_xticks([])
    ax_ann.set_yticklabels(["NG"], fontsize=6)

    # Heatmap
    cmap = LinearSegmentedColormap.from_list("rb", ["#2166AC", "white", "#B2182B"])
    im = ax_hm.imshow(z, aspect="auto", cmap=cmap, vmin=-3, vmax=3,
                       interpolation="nearest")
    ax_hm.set_yticks(np.arange(z.shape[0]))
    ax_hm.set_yticklabels(genes_sub["symbol"].tolist(), fontsize=5.5)
    ax_hm.set_xticks([])
    ax_hm.set_xlabel(f"cells (n={z.shape[1]}), grouped by NhoodGroup",
                      fontsize=7)

    # NG block dividers (vertical lines)
    breaks = np.where(np.diff(col_rank.astype(int)) != 0)[0] + 0.5
    for b in breaks:
        ax_hm.axvline(b, color="white", linewidth=0.8)
        ax_ann.axvline(b, color="white", linewidth=0.8)

    # Row block dividers
    row_rank = genes_sub["source_rank"].values.astype(int)
    row_breaks = np.where(np.diff(row_rank) != 0)[0] + 0.5
    for b in row_breaks:
        ax_hm.axhline(b, color="black", linewidth=0.4, alpha=0.3)

    # NG block labels at top of annotation strip
    boundaries = np.r_[0, breaks, z.shape[1]]
    for i in range(len(boundaries) - 1):
        mid = (boundaries[i] + boundaries[i + 1]) / 2
        rk = int(col_rank[int(boundaries[i])])
        ng_id_at_block = cells_sub.iloc[int(boundaries[i])].NG_dominant
        nf = ng_name_full.get(ng_id_at_block, f"NG_{rk}")
        nf_short = nf.split("_")[-2] if nf.count("_") >= 2 else f"NG_{rk}"
        # name_full like "BMYO-basal_1_BR1" -> show "1"
        med = ng_med_lfc.get(ng_id_at_block, np.nan)
        label = f"{nf_short}\n({med:+.2f})" if not pd.isna(med) else nf_short
        ax_ann.text(mid, -0.6, label, ha="center", va="bottom",
                     fontsize=6, fontweight="bold",
                     color=rank_to_color.get(rk, "#444"))

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("z(log1p CPM)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    out = Path(a.out_pdf) if a.out_pdf else (
        fd / f"panelC_le_doheatmap_{fd.name}.pdf")
    fig.savefig(out, bbox_inches="tight")
    log.info("Wrote %s", out)


if __name__ == "__main__":
    main()
