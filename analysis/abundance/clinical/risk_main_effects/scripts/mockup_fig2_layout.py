#!/usr/bin/env python3
"""
mockup_fig2_layout.py - Fig 2 layout mockup at Nature double-column.

Renders all 5 panels into one PDF showing proposed composition:
  A  cohort alluvial (placeholder box)
  B  per-L2 split-forest (placeholder box)
  C  horizontal multi-L2 dotplot (REAL data from feature artifacts)
  D  pathway NES heatmap (REAL data from pathways.csv sidecars)
  E  caption / size+color legends

Args:
  --inquiry-dir
  --L2-list  comma-separated L2 names (default: 3 selected for main fig)
  --contrast (default parity_x_HR_BRCA1)
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
from matplotlib.patches import Rectangle, FancyBboxPatch
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import scipy.io as sio

RANK_PALETTE = [
    "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
    "#CC6677", "#AA4499", "#882255", "#999933", "#661100",
]

DIV_CMAP = LinearSegmentedColormap.from_list(
    "div", ["#1f77b4", "#f7f7f7", "#d62728"])
EXPR_CMAP = LinearSegmentedColormap.from_list(
    "expr", ["#2166AC", "white", "#B2182B"])


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--L2-list", default="epi::BMYO-basal,epi::LHS-major,str::Fibro-major")
    p.add_argument("--contrast", default="parity_x_HR_BRCA1")
    p.add_argument("--genes-per-l2", type=int, default=20,
                   help="genes per L2 strip in horizontal dotplot")
    p.add_argument("--out-pdf", required=True)
    return p.parse_args()


def load_features(inq, L2, contrast):
    L2_safe = L2.replace("::", "__").replace("-", "_")
    fd = inq / "outputs/panelC_features" / f"{L2_safe}_{contrast}"
    cells = pd.read_csv(fd / "cells.tsv", sep="\t")
    genes = pd.read_csv(fd / "genes.tsv", sep="\t")
    ng_summary = pd.read_csv(fd / "ng_summary.csv")
    pathways = pd.read_csv(fd / "pathways.csv") if (fd / "pathways.csv").exists() else None
    with gzip.open(fd / "counts.mtx.gz", "rb") as f:
        m = sio.mmread(f)
    import scipy.sparse as sp
    X = sp.csr_matrix(m).astype(float)
    return dict(L2=L2, L2_safe=L2_safe, fd=fd,
                 cells=cells, genes=genes, ng_summary=ng_summary,
                 pathways=pathways, X=X)


def aggregate_dotplot_data(d, top_n_genes):
    """Returns pct (NGs x genes), z (NGs x genes), gene_labels, ng_order, ng_med_lfc, ng_rank."""
    ng_summary = d["ng_summary"].sort_values("rank").reset_index(drop=True)
    ng_summary["NhoodGroup"] = ng_summary["NhoodGroup"].astype(str)
    ng_order = ng_summary["NhoodGroup"].tolist()
    ng_rank = ng_summary["rank"].astype(int).tolist()
    ng_med_lfc = ng_summary["group_med_lfc"].astype(float).tolist()

    # Pick top N genes by score; keep block order by source_NG rank
    g = d["genes"].copy()
    g["source_rank"] = g["source_rank"].astype(int)
    g["score"] = g["score"].astype(float)
    g = g.sort_values(["source_rank", "score"], ascending=[True, False]).head(top_n_genes)
    gene_idx = g["row_idx"].astype(int).values  # row indices in counts.mtx

    cells = d["cells"].copy()
    cells["NG_dominant"] = cells["NG_dominant"].astype(str)

    X = d["X"]
    # Subset to top-N genes (cols of X)
    Xsub = X[:, gene_idx].toarray()  # cells x genes
    libsize = X.sum(axis=1).A1
    libsize[libsize == 0] = 1.0
    Xn = np.log1p(Xsub / libsize[:, None] * 1e4)

    pct = np.zeros((len(ng_order), len(gene_idx)))
    mean_expr = np.zeros((len(ng_order), len(gene_idx)))
    for i, ng in enumerate(ng_order):
        mask = (cells["NG_dominant"] == ng).values
        n = mask.sum()
        if n == 0:
            continue
        pct[i, :] = (Xsub[mask, :] > 0).mean(axis=0)
        mean_expr[i, :] = Xn[mask, :].mean(axis=0)

    mu = mean_expr.mean(axis=0, keepdims=True)
    sd = mean_expr.std(axis=0, keepdims=True) + 1e-9
    z = np.clip((mean_expr - mu) / sd, -2.5, 2.5)

    return dict(pct=pct, z=z, gene_labels=g["symbol"].tolist(),
                 ng_order=ng_order, ng_rank=ng_rank, ng_med_lfc=ng_med_lfc,
                 source_rank_per_gene=g["source_rank"].astype(int).tolist())


def draw_placeholder(ax, label, subtitle, color="#f0f0f0"):
    """Labeled rectangle placeholder for unrendered panels."""
    ax.add_patch(Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                             facecolor=color, edgecolor="#888", linewidth=0.5,
                             linestyle="--"))
    ax.text(0.5, 0.55, label, transform=ax.transAxes,
              ha="center", va="center", fontsize=10, fontweight="bold",
              color="#444")
    ax.text(0.5, 0.40, subtitle, transform=ax.transAxes,
              ha="center", va="center", fontsize=7, fontstyle="italic",
              color="#666")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax.spines[sp].set_visible(False)


def draw_horizontal_dotplot_strip(fig, gs_cell, agg, L2_label, max_dot_area=70.0,
                                    show_gene_labels=False):
    """Draws a horizontal strip: NGs as rows, genes as cols.

    Layout within the cell:
      [direction bar | rank strip | dotplot grid | NG name labels]
    """
    # Subdivide the cell vertically
    inner = gs_cell.subgridspec(3, 2,
                                  height_ratios=[0.5, 0.5, 5],
                                  width_ratios=[0.18, 1.0],
                                  hspace=0.05, wspace=0.02)
    ax_lbl = fig.add_subplot(inner[2, 0])  # NG name labels (left)
    ax_dir = fig.add_subplot(inner[0, 1])  # direction bar (top)
    ax_rk = fig.add_subplot(inner[1, 1])   # rank strip
    ax = fig.add_subplot(inner[2, 1])      # dotplot grid

    n_ngs = len(agg["ng_order"])
    n_genes = len(agg["gene_labels"])

    # Direction bar (1 row, n_ngs cols)
    med = np.array(agg["ng_med_lfc"])
    lfc_lim = max(np.abs(med).max(), 0.5)
    ax_dir.imshow(np.full((1, n_genes), np.nan), aspect="auto")
    # Actually paint per-NG bars: but they're in rows. So show a SIDE bar (transpose) - LEFT side
    # Wait: NG is rows here. So direction bar should be the LEFT column.
    # Let me redo: direction & rank strips go ON THE LEFT (column annotations for rows).
    ax_dir.set_visible(False); ax_rk.set_visible(False)

    # Replace with left-side strips: use ax_lbl for combined name + direction
    ax_lbl.clear()
    # Build per-NG label panel: small colored box + name
    rank_to_color = {i + 1: c for i, c in enumerate(RANK_PALETTE)}
    for i, (ng, rk, m) in enumerate(zip(agg["ng_order"], agg["ng_rank"],
                                          agg["ng_med_lfc"])):
        # rank color box
        ax_lbl.add_patch(Rectangle((0.0, n_ngs - i - 1), 0.15, 1,
                                     facecolor=rank_to_color.get(int(rk), "#888"),
                                     edgecolor="white", linewidth=0.5))
        # direction box
        col = DIV_CMAP((m + lfc_lim) / (2 * lfc_lim))
        ax_lbl.add_patch(Rectangle((0.18, n_ngs - i - 1), 0.15, 1,
                                     facecolor=col, edgecolor="white", linewidth=0.5))
        # text
        ax_lbl.text(0.40, n_ngs - i - 0.5, f"{rk}",
                     ha="left", va="center", fontsize=6, fontweight="bold",
                     color=rank_to_color.get(int(rk), "#222"))
        ax_lbl.text(0.55, n_ngs - i - 0.5, f"{m:+.1f}",
                     ha="left", va="center", fontsize=5.5, color="#444")
    ax_lbl.set_xlim(0, 1); ax_lbl.set_ylim(0, n_ngs)
    ax_lbl.set_xticks([]); ax_lbl.set_yticks([])
    for sp in ("top", "right", "bottom", "left"):
        ax_lbl.spines[sp].set_visible(False)
    ax_lbl.text(-0.05, n_ngs + 0.3, L2_label, fontsize=7, fontweight="bold",
                 ha="left", va="bottom", transform=ax_lbl.transData)

    # Dotplot grid
    xs, ys = np.meshgrid(np.arange(n_genes), np.arange(n_ngs)[::-1], indexing="xy")
    sizes = agg["pct"] * max_dot_area
    sizes = np.where(sizes < 1.0, 0.0, sizes)
    sc = ax.scatter(xs.ravel(), ys.ravel(),
                     s=sizes.ravel(),
                     c=agg["z"].ravel(), cmap=EXPR_CMAP, vmin=-2.5, vmax=2.5,
                     edgecolors="#666666", linewidths=0.15)

    ax.set_xlim(-0.5, n_genes - 0.5)
    ax.set_ylim(-0.5, n_ngs - 0.5)

    # Source-NG block dividers (vertical lines)
    src_rank = np.array(agg["source_rank_per_gene"])
    breaks = np.where(np.diff(src_rank) != 0)[0] + 0.5
    for b in breaks:
        ax.axvline(b, color="#aaaaaa", linewidth=0.4, alpha=0.5)

    if show_gene_labels:
        ax.set_xticks(np.arange(n_genes))
        ax.set_xticklabels(agg["gene_labels"], fontsize=4.5, rotation=90, ha="center")
    else:
        ax.set_xticks([])
    ax.set_yticks([])
    ax.tick_params(axis="both", which="both", length=0, pad=1)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.3)
    ax.spines["bottom"].set_linewidth(0.3)
    return sc, ax


def draw_pathway_heatmap(ax, datasets, top_paths_per_l2=8):
    """Pathway NES heatmap: rows = pathways, cols = (L2, NG).
    Picks top pathways per L2 by max(|NES|) across that L2's NGs.
    """
    # Aggregate pathway records across L2s
    all_p = []
    for d in datasets:
        if d["pathways"] is None:
            continue
        p = d["pathways"].copy()
        p["L2"] = d["L2"]
        all_p.append(p)
    if not all_p:
        ax.text(0.5, 0.5, "no pathway records",
                  transform=ax.transAxes, ha="center", va="center")
        return
    df = pd.concat(all_p, ignore_index=True)

    # Pick top-N pathways per L2 by max |NES| within L2
    top_paths = (df.assign(absNES=df.NES.abs())
                  .sort_values("absNES", ascending=False)
                  .groupby("L2", sort=False)
                  .head(top_paths_per_l2))
    pathways_to_show = top_paths.pathway.unique().tolist()

    # Build matrix: rows = pathways, cols = (L2, NG_rank) ordered
    # Cols in (L2, rank) order
    col_keys = []
    for d in datasets:
        ng_summary = d["ng_summary"].sort_values("rank")
        for _, r in ng_summary.iterrows():
            col_keys.append((d["L2"], int(r["rank"]), r["NhoodGroup"]))
    n_cols = len(col_keys)
    n_rows = len(pathways_to_show)
    mat = np.full((n_rows, n_cols), np.nan)

    for ci, (L2, rank, ng) in enumerate(col_keys):
        sub = df[(df.L2 == L2) & (df.NhoodGroup.astype(str) == str(ng))]
        for _, row in sub.iterrows():
            if row["pathway"] in pathways_to_show:
                ri = pathways_to_show.index(row["pathway"])
                mat[ri, ci] = row["NES"]

    nes_lim = max(np.nanmax(np.abs(mat)), 0.5) if not np.all(np.isnan(mat)) else 1.5
    im = ax.imshow(mat, aspect="auto", cmap=DIV_CMAP, vmin=-nes_lim, vmax=nes_lim,
                    interpolation="nearest")
    ax.set_xticks(np.arange(n_cols))
    col_labels = [f"{L2.split('::')[1].split('-')[0]}.{rank}"
                    for L2, rank, _ in col_keys]
    ax.set_xticklabels(col_labels, fontsize=4.5, rotation=90)
    ax.set_yticks(np.arange(n_rows))
    short_paths = [p.replace("HALLMARK_", "H_").replace("REACTOME_", "R_")[:30]
                    for p in pathways_to_show]
    ax.set_yticklabels(short_paths, fontsize=4.5)
    ax.tick_params(axis="both", which="both", length=0, pad=1)

    # L2 group dividers
    boundaries = []
    cur = None
    for ci, (L2, _, _) in enumerate(col_keys):
        if L2 != cur:
            if ci > 0:
                boundaries.append(ci - 0.5)
            cur = L2
    for b in boundaries:
        ax.axvline(b, color="white", linewidth=1.5)

    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.3)
    ax.spines["bottom"].set_linewidth(0.3)
    return im


def main():
    a = parse_args()
    L2_list = [s.strip() for s in a.L2_list.split(",")]
    print("L2s:", L2_list)
    datasets = [load_features(a.inquiry_dir, L2, a.contrast) for L2 in L2_list]
    aggs = [aggregate_dotplot_data(d, a.genes_per_l2) for d in datasets]

    # ---- Figure: Nature double-column = 183mm = 7.2in wide; ~9.5in tall page ----
    fig = plt.figure(figsize=(7.2, 9.5))
    outer = gridspec.GridSpec(5, 1, figure=fig,
                                height_ratios=[1.5, 2.2, 3.6, 2.0, 0.4],
                                hspace=0.35, left=0.06, right=0.97,
                                top=0.97, bottom=0.03)

    # ----- Panel A: cohort alluvial (placeholder box) -----
    ax_A = fig.add_subplot(outer[0])
    draw_placeholder(ax_A,
                       "Panel A   cohort design alluvial",
                       "study -> facs -> genotype -> ch -> risk_class -> parity\n"
                       "+ comparison arcs: parity_in_AR <- paired -> parity_x_HR_BRCA1",
                       color="#fff8e8")
    ax_A.text(0.01, 1.04, "a", transform=ax_A.transAxes,
                fontsize=9, fontweight="bold", va="bottom")

    # ----- Panel B: split-violin per-L2 forest (placeholder) -----
    ax_B = fig.add_subplot(outer[1])
    draw_placeholder(ax_B,
                       "Panel B   per-L2 split-violin + NhoodGroup forest dots",
                       "L2 = x (rotated, 90deg);  logFC = y;  faceted by compartment\n"
                       "L=parity_in_AR (parity main effect)   |   R=parity_x_HR_BRCA1 (BRCA1 mod)\n"
                       "forest dots = NG median +- IQR, rank-colored, sized by n_nhoods",
                       color="#eaf3fb")
    ax_B.text(0.01, 1.04, "b", transform=ax_B.transAxes,
                fontsize=9, fontweight="bold", va="bottom")

    # ----- Panel C: horizontal multi-L2 dotplot (REAL) -----
    panelC_outer = outer[2].subgridspec(len(L2_list) + 1, 1,
                                          height_ratios=[1.0] * len(L2_list) + [0.25],
                                          hspace=0.25)
    last_sc = None
    for i, (L2, agg) in enumerate(zip(L2_list, aggs)):
        is_last = (i == len(L2_list) - 1)
        sc, ax = draw_horizontal_dotplot_strip(fig, panelC_outer[i], agg,
                                                  L2_label=L2,
                                                  show_gene_labels=is_last)
        last_sc = sc
        if i == 0:
            ax.text(-0.02, 1.32, "c", transform=ax.transAxes,
                      fontsize=9, fontweight="bold", va="bottom", ha="right")
            ax.text(0.5, 1.20,
                      "Panel C   per-NhoodGroup leading-edge program (interaction story)",
                      transform=ax.transAxes, fontsize=8, ha="center",
                      fontweight="bold")

    # ----- Panel D: pathway NES heatmap (REAL) -----
    ax_D = fig.add_subplot(outer[3])
    im_D = draw_pathway_heatmap(ax_D, datasets, top_paths_per_l2=8)
    ax_D.text(-0.02, 1.04, "d", transform=ax_D.transAxes,
                fontsize=9, fontweight="bold", va="bottom", ha="right")
    ax_D.text(0.5, 1.04,
                "Panel D   top-NES pathways per NG (Hallmark + Reactome)",
                transform=ax_D.transAxes, fontsize=8, ha="center", fontweight="bold")

    # ----- Panel E: legends (real, shared across C+D) -----
    ax_E = fig.add_subplot(outer[4])
    ax_E.axis("off")
    # color legend (color = z-mean log1p CPM, same for C and D after rescale)
    cax_color = fig.add_axes([0.10, 0.04, 0.20, 0.012])
    cbar = fig.colorbar(last_sc, cax=cax_color, orientation="horizontal")
    cbar.set_label("z(mean log1p CPM)  /  NES", fontsize=6)
    cbar.ax.tick_params(labelsize=5)

    # size legend (% cells expressing)
    ax_E.text(0.42, 0.55, "% cells expressing:", fontsize=6,
                ha="left", va="center", transform=ax_E.transAxes)
    pcts = [0.10, 0.25, 0.50, 0.80]
    for i, p in enumerate(pcts):
        ax_E.scatter(0.55 + i * 0.06, 0.55, s=p * 70.0,
                       color="#888", edgecolor="#666", linewidth=0.15,
                       transform=ax_E.transAxes)
        ax_E.text(0.55 + i * 0.06, 0.20, f"{int(p*100)}%", fontsize=5,
                    ha="center", transform=ax_E.transAxes)

    # NG direction key
    ax_E.text(0.82, 0.55, "NG direction:", fontsize=6,
                ha="left", va="center", transform=ax_E.transAxes)
    for i, (lbl, c) in enumerate([("med-", "#1f77b4"), ("0", "#f7f7f7"), ("med+", "#d62728")]):
        ax_E.add_patch(Rectangle((0.92 + i * 0.02, 0.50), 0.018, 0.10,
                                    transform=ax_E.transAxes,
                                    facecolor=c, edgecolor="#888", linewidth=0.3))

    fig.savefig(a.out_pdf, bbox_inches="tight", dpi=300)
    print(f"Wrote {a.out_pdf}")


if __name__ == "__main__":
    main()
