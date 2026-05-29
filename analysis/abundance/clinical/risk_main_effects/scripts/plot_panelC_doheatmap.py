#!/usr/bin/env python3
"""
plot_panelC_doheatmap.py - cell × gene heatmap (DoHeatmap-style) per L2.

Cell-grain version of Panel C. Cells (columns) grouped by their dominant
NhoodGroup; genes (rows) are top markers per NhoodGroup from F.3 limma-voom.
Top column-annotation strip = NG identity colored by within-L2 size rank.

Reads pre-extracted nhoods sparse matrix (from export_milo_nhoods.R) so we
don't reload milo per L2. Counts come from per-compartment iHBCAv1 assembly
counts.npz (read directly).

Args:
  --L2 e.g. "epi::BMYO-basal"
  --contrast e.g. "parity_x_HR_BRCA1"
  --inquiry-dir
  --top-n-per-ng (default 12)
  --max-cells-per-ng (default 200)
  --counts-root (default to iHBCAv1 assembly outputs/components)
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from collections import Counter

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")
import sys as _sys
from pathlib import Path as _Path
_sys.path.insert(0, str(_Path(__file__).resolve().parent))
from lib.load_paths import load_inquiry  # noqa: E402

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%H:%M:%S")
log = logging.getLogger("panelC_doheatmap")


RANK_PALETTE = [
    "#332288", "#117733", "#44AA99", "#88CCEE", "#DDCC77",
    "#CC6677", "#AA4499", "#882255", "#999933", "#661100",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--L2", required=True)
    p.add_argument("--contrast", default="parity_x_HR_BRCA1")
    p.add_argument("--inquiry-dir", required=True, type=Path)
    p.add_argument("--top-n-per-ng", type=int, default=12)
    p.add_argument("--max-cells-per-ng", type=int, default=200)
    p.add_argument("--counts-root", type=Path,
                   default=None,
                   help="Override; default resolves via paths.yaml inputs.components_root")
    return p.parse_args()


def main():
    a = parse_args()
    inq = Path(a.inquiry_dir)
    L2_compartment, L2_label = a.L2.split("::", 1)
    L2_safe = a.L2.replace("::", "__").replace("-", "_")
    comp_short = {"epi": "epi", "imm": "imm", "str": "str"}.get(L2_compartment)
    if comp_short is None:
        raise SystemExit(f"unknown compartment: {L2_compartment}")

    # ---- Load lookup (viable NGs in target L2 + contrast) ----
    lookup = pd.read_csv(inq / "outputs/nhoodgroup_renaming/lookup.csv")
    lookup = lookup[(lookup.viable) &
                    (lookup.contrast == a.contrast) &
                    (lookup.parent_L2_joint == a.L2)].copy()
    lookup["NhoodGroup"] = lookup["NhoodGroup"].astype(str)
    lookup = lookup.sort_values("size_rank_within_L2").reset_index(drop=True)
    if lookup.empty:
        raise SystemExit(f"No viable NGs for {a.L2} in {a.contrast}")
    log.info("Viable NGs: %d", len(lookup))
    log.info("\n%s", lookup[["NhoodGroup", "n_nhoods_in_group", "group_med_lfc",
                              "size_rank_within_L2"]].to_string(index=False))

    # ---- Load top markers per NG ----
    marker_dir = inq / "outputs/stageF3_markers" / a.contrast
    top_markers = []
    for _, r in lookup.iterrows():
        ng = r.NhoodGroup
        rk = r.size_rank_within_L2
        fname = f"{L2_compartment}__{L2_label.replace('-', '_')}_{ng}_vs_parent.csv"
        fp = marker_dir / fname
        if not fp.exists():
            log.warning("missing marker file: %s", fp.name); continue
        m = pd.read_csv(fp)
        m = m[(m.is_dissoc_stress != True) &
              (m.symbol.notna()) & (m.symbol != "") &
              (m.FDR < 0.05) & (m.logFC.abs() > np.log2(1.25))]
        m = m.sort_values("logFC", key=lambda s: -s.abs()).head(a.top_n_per_ng)
        m["rank"] = rk
        m["ng"] = ng
        top_markers.append(m)
    markers_df = pd.concat(top_markers, ignore_index=True)
    log.info("Marker rows: %d, unique symbols: %d",
             len(markers_df), markers_df.symbol.nunique())

    # ---- Load atlas cells + nhoods sparse ----
    nhoods_dir = inq / "outputs/atlas_nhoods"
    cells_atlas = pd.read_csv(nhoods_dir / "cells.tsv", header=None)[0].values
    log.info("atlas cells: %d", len(cells_atlas))
    log.info("Loading atlas_nhoods.mtx ...")
    nh_mat = sio.mmread(str(nhoods_dir / "nhoods.mtx")).tocsr()
    log.info("atlas_nhoods sparse: %s, nnz=%d", nh_mat.shape, nh_mat.nnz)

    # nhood_groups.csv: Nhood -> NhoodGroup map (only entries present here)
    ng_csv = inq / "outputs/stageF1_nhoodgroups" / a.contrast / "nhood_groups.csv"
    ng_map_full = pd.read_csv(ng_csv)
    ng_map_full["NhoodGroup"] = ng_map_full["NhoodGroup"].astype(str)
    viable_set = set(lookup.NhoodGroup)
    ng_map_full = ng_map_full[ng_map_full.NhoodGroup.isin(viable_set)]
    log.info("nhoods in viable NGs: %d", len(ng_map_full))
    nhood_to_ng = dict(zip(ng_map_full["Nhood"].astype(int),
                            ng_map_full["NhoodGroup"]))

    # ---- Load counts.npz for compartment ----
    counts_npz = a.counts_root / f"{comp_short}_counts.npz"
    log.info("Loading %s ...", counts_npz)
    loaded = np.load(counts_npz, allow_pickle=True)
    if "data" in loaded.files:
        X = sp.csr_matrix(
            (loaded["data"], loaded["indices"], loaded["indptr"]),
            shape=tuple(loaded["shape"]),
        )
    else:
        raise SystemExit(f"counts.npz unexpected keys: {loaded.files}")
    log.info("counts: %s (cells x genes)", X.shape)

    # ---- Companion metadata: cell barcodes + gene ENSG IDs ----
    # By convention from 01_compute_markers.py the cell row order matches a
    # `<comp>_cells.csv` (or barcoded leiden csv). gene_data.csv supplies ENSG.
    # Per-compartment cell metadata in the assembly outputs:
    #   {comp}_metadata.csv: cell_id (or barcode) row order matches counts.npz rows
    cand_paths = [
        a.counts_root / f"{comp_short}_metadata.csv",
        a.counts_root / f"{comp_short}_metadata_enriched.csv",
        a.counts_root / f"{comp_short}_cells.csv",
    ]
    cells_comp = None
    for cp in cand_paths:
        if cp.exists():
            cells_comp = pd.read_csv(cp)
            log.info("compartment metadata: %s", cp)
            break
    if cells_comp is None:
        raise SystemExit(f"No compartment metadata file: tried {cand_paths}")
    if "cell_id" not in cells_comp.columns:
        # try barcode column
        for cand in ("barcode", "cell", "obs_names", "Cell"):
            if cand in cells_comp.columns:
                cells_comp = cells_comp.rename(columns={cand: "cell_id"}); break
    cells_comp_ids = cells_comp["cell_id"].values if "cell_id" in cells_comp.columns \
                      else cells_comp.iloc[:, 0].values
    log.info("compartment cells (counts row order): %d", len(cells_comp_ids))
    if X.shape[0] != len(cells_comp_ids):
        log.warning("counts rows (%d) != cells_comp (%d)",
                    X.shape[0], len(cells_comp_ids))

    gene_data = pd.read_csv(a.counts_root / "gene_data.csv", index_col=0)
    if "symbol" not in gene_data.columns:
        for cand in ("gene_symbol", "Symbol", "symbol_"):
            if cand in gene_data.columns:
                gene_data = gene_data.rename(columns={cand: "symbol"}); break
    log.info("gene_data: %d entries", len(gene_data))

    # ---- Filter cells to target L2 (via labels_full.csv) ----
    cfg = load_inquiry(a.inquiry_dir)["paths"]["inputs"]
    if a.counts_root is None:
        a.counts_root = Path(cfg["components_root"])
    labels_path = Path(cfg["labels"])
    labels = pd.read_csv(labels_path,
                          usecols=["cell_id", "compartment", "label", "is_artifact"])
    labels.is_artifact = labels.is_artifact.astype(str).str.lower().isin({"true", "1"})
    l2_cells = set(labels[(labels.compartment == L2_compartment) &
                          (labels.label == L2_label) &
                          (~labels.is_artifact)].cell_id)
    log.info("Cells in %s per labels_full: %d", a.L2, len(l2_cells))

    # Cell row indices in counts.npz that are in target L2
    in_l2_mask = np.isin(cells_comp_ids, list(l2_cells))
    l2_idx_in_counts = np.where(in_l2_mask)[0]
    l2_cells_ordered = cells_comp_ids[l2_idx_in_counts]
    log.info("L2 cells in counts: %d", len(l2_idx_in_counts))
    if len(l2_idx_in_counts) == 0:
        raise SystemExit("No L2 cells found in counts barcodes")

    # ---- Cell -> dominant NhoodGroup via atlas nhoods ----
    # Match L2 cell barcodes to atlas (milo) cell index
    atlas_cell_to_idx = {c: i for i, c in enumerate(cells_atlas)}
    l2_atlas_idx = np.array([atlas_cell_to_idx.get(c, -1) for c in l2_cells_ordered])
    in_atlas_mask = l2_atlas_idx >= 0
    if not in_atlas_mask.all():
        log.warning("Cells not in atlas nhoods: %d", (~in_atlas_mask).sum())
    keep_l2 = np.where(in_atlas_mask)[0]
    l2_atlas_idx = l2_atlas_idx[in_atlas_mask]

    # nhoods columns we care about (in viable NGs)
    nh_keep_cols = sorted([n - 1 for n in nhood_to_ng.keys()  # 1-indexed -> 0
                           if (n - 1) < nh_mat.shape[1]])
    log.info("nhood columns of interest: %d", len(nh_keep_cols))
    nh_sub = nh_mat[l2_atlas_idx, :][:, nh_keep_cols]
    log.info("L2 cells x viable nhoods sparse: %s", nh_sub.shape)

    # Map column index back to NhoodGroup
    ng_per_col = np.array([nhood_to_ng[c + 1] for c in nh_keep_cols])

    # Dominant NG per cell
    dominant_ng = np.array(["__none__"] * nh_sub.shape[0], dtype=object)
    nh_sub_lil = nh_sub.tolil()
    for i in range(nh_sub.shape[0]):
        cols = nh_sub_lil.rows[i]
        if not cols: continue
        ngs = ng_per_col[cols]
        c = Counter(ngs); dominant_ng[i] = c.most_common(1)[0][0]
    has_ng = dominant_ng != "__none__"
    log.info("L2 cells with viable NG assignment: %d / %d",
             has_ng.sum(), len(has_ng))
    keep_l2 = keep_l2[has_ng]
    dominant_ng = dominant_ng[has_ng]
    final_l2_idx = l2_idx_in_counts[keep_l2]
    final_cell_ids = l2_cells_ordered[keep_l2]
    log.info("cells per NG:")
    for ng, n in Counter(dominant_ng).items():
        log.info("  %s : %d", ng, n)

    # Downsample per NG for plot legibility
    rng = np.random.default_rng(42)
    keep_idx = []
    keep_ng = []
    keep_cells = []
    for ng in lookup.NhoodGroup:
        ix = np.where(dominant_ng == ng)[0]
        if len(ix) > a.max_cells_per_ng:
            ix = rng.choice(ix, a.max_cells_per_ng, replace=False)
        keep_idx.extend(final_l2_idx[ix].tolist())
        keep_ng.extend([ng] * len(ix))
        keep_cells.extend(final_cell_ids[ix].tolist())
    keep_idx = np.array(keep_idx)
    keep_ng = np.array(keep_ng)
    log.info("Heatmap cells (downsampled): %d", len(keep_idx))

    # ---- Marker-gene subset (ENSG -> column index) ----
    sym2ens = dict(zip(markers_df.symbol, markers_df.gene_id))
    ens_genes = [sym2ens.get(s) for s in markers_df.symbol]
    # gene_data index is ENSG; row id in counts is gene_data row order
    ens_index = list(gene_data.index)
    ens_to_col = {e: i for i, e in enumerate(ens_index)}
    gene_cols = []
    display_syms = []
    for ens, sym in zip(ens_genes, markers_df.symbol):
        if ens in ens_to_col:
            gene_cols.append(ens_to_col[ens])
            display_syms.append(sym)
    gene_cols = np.array(gene_cols)
    log.info("Genes available: %d / %d", len(gene_cols), len(markers_df))

    # Order rows by source-NG rank (preserve group blocks)
    gene_to_rank = (markers_df.groupby("symbol", sort=False).rank.first()
                     .to_dict())
    row_rank = np.array([gene_to_rank.get(s, max(lookup.size_rank_within_L2) + 1)
                          for s in display_syms])
    sort_rows = np.argsort(row_rank, kind="stable")
    gene_cols = gene_cols[sort_rows]
    display_syms = [display_syms[i] for i in sort_rows]
    row_rank = row_rank[sort_rows]
    # Deduplicate genes (keep first occurrence)
    seen = set(); dedup = []
    for i, s in enumerate(display_syms):
        if s not in seen:
            seen.add(s); dedup.append(i)
    gene_cols = gene_cols[dedup]
    display_syms = [display_syms[i] for i in dedup]
    row_rank = row_rank[dedup]

    # ---- Slice counts (cells x marker-genes) ----
    sub = X[keep_idx, :][:, gene_cols].toarray().astype(float)
    log.info("expr submatrix: %s", sub.shape)

    # log1p of CPM-like normalization (per-cell library size)
    libsize = np.asarray(X[keep_idx, :].sum(axis=1)).ravel()
    libsize[libsize == 0] = 1
    sub_norm = np.log1p((sub.T / libsize * 1e4).T)

    # z-score per gene
    z = (sub_norm - sub_norm.mean(axis=0)) / (sub_norm.std(axis=0) + 1e-9)
    z = np.clip(z.T, -3, 3)  # rows = genes, cols = cells

    # Order columns by NG rank then within-NG randomly (stable)
    ng_rank_lookup = dict(zip(lookup.NhoodGroup,
                                lookup.size_rank_within_L2))
    col_rank = np.array([ng_rank_lookup.get(g, 99) for g in keep_ng])
    col_order = np.argsort(col_rank, kind="stable")
    z = z[:, col_order]
    keep_ng = keep_ng[col_order]
    col_rank = col_rank[col_order]

    # ---- Render ----
    fig_h = max(4.0, 0.13 * z.shape[0] + 1.5)
    fig_w = max(7.0, 0.012 * z.shape[1] + 3.5)
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.5, 20], width_ratios=[40, 1],
                          hspace=0.05, wspace=0.05)
    ax_ann = fig.add_subplot(gs[0, 0])
    ax_hm  = fig.add_subplot(gs[1, 0], sharex=ax_ann)
    ax_cb  = fig.add_subplot(gs[1, 1])

    # NG annotation strip
    rank_to_color = {i + 1: c for i, c in enumerate(RANK_PALETTE)}
    ann_colors = np.array([rank_to_color.get(int(r), "#888888")
                            for r in col_rank])
    # Draw as vertical lines / image
    rgb = [matplotlib.colors.to_rgb(c) for c in ann_colors]
    ax_ann.imshow(np.array(rgb).reshape(1, -1, 3), aspect="auto")
    ax_ann.set_yticks([]); ax_ann.set_xticks([])
    ax_ann.set_ylabel("NG", rotation=0, ha="right", va="center", fontsize=7)

    # Heatmap
    cmap = LinearSegmentedColormap.from_list("rb", ["#2166AC", "white", "#B2182B"])
    im = ax_hm.imshow(z, aspect="auto", cmap=cmap, vmin=-3, vmax=3,
                       interpolation="nearest")
    ax_hm.set_yticks(np.arange(z.shape[0]))
    ax_hm.set_yticklabels(display_syms, fontsize=5.5)
    ax_hm.set_xticks([])
    ax_hm.set_xlabel(f"cells in {a.L2}, grouped by NhoodGroup (n={z.shape[1]})",
                     fontsize=7)

    # NG block dividers
    breaks = np.where(np.diff(col_rank) != 0)[0] + 0.5
    for b in breaks:
        ax_hm.axvline(b, color="white", linewidth=0.8)
        ax_ann.axvline(b, color="white", linewidth=0.8)

    # NG block labels at top of annotation
    boundaries = np.r_[0, breaks, z.shape[1]]
    for i in range(len(boundaries) - 1):
        mid = (boundaries[i] + boundaries[i + 1]) / 2
        rk = int(col_rank[int(boundaries[i])])
        ax_ann.text(mid, -0.6, f"NG_{rk}", ha="center", va="bottom",
                     fontsize=6.5, fontweight="bold",
                     color=rank_to_color.get(rk, "#444"))

    cb = fig.colorbar(im, cax=ax_cb)
    cb.set_label("z(log1p CPM)", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    fig.suptitle("", fontsize=8)
    out_dir = inq / "reports/figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / f"panelC_{L2_safe}_{a.contrast}_doheatmap.pdf"
    fig.savefig(out_pdf, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)


if __name__ == "__main__":
    main()
