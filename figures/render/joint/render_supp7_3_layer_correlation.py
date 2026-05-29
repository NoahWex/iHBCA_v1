"""
render_supp7_3_layer_correlation.py — Fig 3 Supp 7-3.

Cross-platform expression correlation: per-cell-type scatter of mean
log-norm marker-gene expression in Xenium nuclear-segmented counts (x)
vs FLEX scvi_n100 per-compartment counts (y). Pearson r + Spearman p
quantify cross-platform agreement at the cell-type level.

Cytoplasmic-segmented Xenium counts are intentionally excluded — the
Xenium annotation cascade is built on nuclear counts, so nuclear is the
canonical reference for cross-platform concordance. Nuclear-vs-cyto
agreement is a separate within-Xenium-pipeline question.

Inputs:
    paths.xenium_nuclear_bundle        xenium_cells.tsv + xenium_genes.tsv
                                        + xenium_*counts.mtx.gz
                                        (+ xenium_obs.csv ignored here)
    paths.flex_compartment_bundle_root  per-compartment FLEX scvi_n100
                                        bundles {Immune,Epithelial,
                                        Stromal}/scvi_n100/{counts.mtx.gz,
                                        genes.tsv, cells.tsv}
    paths.xenium_compartment_dir        all_compartments.csv (cell_id,
                                        compartment, label_short,
                                        is_artifact, platform, ...)

Outputs:
    supp7_3_layer_correlation.pdf          per-cell-type scatter grid
    supp7_3_layer_correlation_summary.pdf  all cell types overlaid
    supp7_3_layer_correlation_stats.csv    per-cell-type r, rho, p, n

Pattern source:
    Spatial_HBCA_Xenium/pipeline/scripts/10_annotation_panels/
    06_render_layer_correlation.py — MARKERS / FULL_NAMES / LABEL_COMP /
    LABEL_ORDER ported verbatim; FLEX/Nuc loaders ported with minor
    refactors. Inline COMP_COLOR replaced with get_palette('compartment'),
    inline figsize replaced with get_dimensions, cyto comparison dropped,
    Spearman p added to per-cell-type stats.

Framework conformance: Option B per CP_supp7_python_framework_gap.
"""

from __future__ import annotations

import argparse
import gzip
import io
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import scipy.stats

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
    get_dimensions,
    get_matplotlib_theme,
    get_palette,
    load_aesthetics,
)
from load_paths import load_paths, resolve_path  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp7_3")


COMP_DIRS = {"Immune": "Immune", "Epithelial": "Epithelial", "Stromal": "Stromal"}

# Map l15_canonical (column in fig3_dotplot_features.csv) -> label_short
# (column in all_compartments.csv). The promoted features file uses the
# fine-grained L1.5 canonical names; the data uses the abbreviated label_short.
# BMYO-myo and BMYO-basal both roll up to a single label_short (BMYO and BMNC
# respectively) — the data annotation cascade does this rollup.
L15_CANONICAL_TO_LABEL_SHORT = {
    "BMYO-myo":              "BMYO",
    "BMYO-basal":            "BMNC",
    "LASP":                  "LASP",
    "LASP-basal":            "LASP-basal",
    "LHS":                   "LHS",
    "Fibroblast":            "Fb",
    "Fibroblast_activated":  "Fb_Activated",
    "Fibroblast_SFRP4":      "Fb_SFRP4",
    "Endothelial":           "EC",
    "Lymphatic Endothelial": "LEC",
    "Pericyte":              "PV",
    "Adipocyte":             "Adipo",
    "Macrophage":            "Mac",
    "cDC1":                  "cDC1",
    "cDC2":                  "cDC2",
    "pDC":                   "pDC",
    "Mast cell":             "Mast",
    "Neutrophil":            "Neu",
    "CD4 T cell":            "CD4T",
    "CD4 Treg":              "Treg",
    "CD8 T cell":            "CD8T",
    "T-NK":                  "T-NK",
    "B cell":                "B",
    "Plasma cell":           "Plas",
}


def load_canonical_markers(features_csv: Path) -> tuple[dict, list, dict, dict, list]:
    """Load `fig3_dotplot_features.csv` and build:

        markers       : dict label_short -> [feature, ...]   (canonical set per type)
        label_order   : ordered list of label_short values in the CSV's row order
        full_names    : dict label_short -> human display name (the l15_canonical
                        value, deduplicated to first-seen)
        label_comp    : dict label_short -> compartment (from lbridge_compartment)
        all_markers   : ordered union of every feature in the CSV
    """
    csv = pd.read_csv(features_csv)
    required = {"feature", "lbridge_compartment", "l15_canonical"}
    missing = required - set(csv.columns)
    if missing:
        raise KeyError(f"fig3_dotplot_features.csv missing columns: {missing}")

    markers: dict[str, list[str]] = {}
    label_order: list[str] = []
    full_names: dict[str, str] = {}
    label_comp: dict[str, str] = {}
    all_markers: list[str] = []
    unmapped_l15: set[str] = set()

    for _, row in csv.iterrows():
        l15 = str(row["l15_canonical"]).strip()
        feature = str(row["feature"]).strip()
        comp = str(row["lbridge_compartment"]).strip()
        lbl = L15_CANONICAL_TO_LABEL_SHORT.get(l15)
        if lbl is None:
            unmapped_l15.add(l15)
            continue
        if lbl not in markers:
            markers[lbl] = []
            label_order.append(lbl)
            full_names[lbl] = l15
            label_comp[lbl] = comp
        if feature not in markers[lbl]:
            markers[lbl].append(feature)
        if feature not in all_markers:
            all_markers.append(feature)

    if unmapped_l15:
        log.warning("l15_canonical values without label_short mapping: %s",
                    sorted(unmapped_l15))
    return markers, label_order, full_names, label_comp, all_markers


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 5 cell types for quick sanity render.")
    return p.parse_args()


def load_xenium_means(
    bundle_path: Path, ann: pd.DataFrame,
    all_markers: list, label_order: list,
) -> tuple[dict, list]:
    """Per-cell-type mean log-norm marker expression from Xenium nuclear bundle."""
    genes = pd.read_csv(bundle_path / "xenium_genes.tsv", header=None)[0].tolist()
    cells = pd.read_csv(bundle_path / "xenium_cells.tsv", header=None)[0].tolist()
    gene_idx = {g: i for i, g in enumerate(genes)}
    cell_idx = {c: i for i, c in enumerate(cells)}

    mtx_files = list(bundle_path.glob("xenium_*counts.mtx.gz"))
    if not mtx_files:
        raise FileNotFoundError(f"No xenium_*counts.mtx.gz in {bundle_path}")
    mtx_file = mtx_files[0]
    log.info("loading nuclear MTX (%d cells × %d genes): %s", len(cells), len(genes), mtx_file)
    with gzip.open(mtx_file, "rb") as gz:
        mat = sp.csr_matrix(scipy.io.mmread(io.BytesIO(gz.read())))

    present = [g for g in all_markers if g in gene_idx]
    marker_cols = [gene_idx[g] for g in present]
    log.info("nuclear markers present: %d / %d", len(present), len(all_markers))

    xen_ann = ann[ann["platform"] == "xenium"].copy()
    xen_ann = xen_ann[xen_ann["cell_id"].isin(cell_idx)]
    log.info("nuclear: matched %d xenium cells in annotation", len(xen_ann))

    row_idx = np.array([cell_idx[c] for c in xen_ann["cell_id"]])
    labels = xen_ann["label_short"].values
    row_totals = np.asarray(mat[row_idx].sum(axis=1)).flatten()
    row_totals[row_totals == 0] = 1

    mat_sub = mat[:, marker_cols]
    means = {}
    for lbl in label_order:
        mask = labels == lbl
        if mask.sum() == 0:
            means[lbl] = np.zeros(len(present))
            continue
        sub = mat_sub[row_idx[mask]].toarray().astype(float)
        totals = row_totals[mask, np.newaxis]
        sub_norm = np.log1p(sub / totals * 1e4)
        means[lbl] = sub_norm.mean(axis=0)

    return means, present


def load_flex_means(
    flex_root: Path, ann: pd.DataFrame,
    all_markers: list, label_order: list,
) -> tuple[dict, list]:
    """Per-cell-type mean log-norm marker expression from FLEX per-compartment bundles."""
    flex_ann = ann[ann["platform"] == "flex"].copy()
    means: dict[str, dict[str, float]] = {lbl: {} for lbl in label_order}

    for comp, subdir in COMP_DIRS.items():
        cdir = flex_root / subdir / "scvi_n100"
        if not cdir.exists():
            log.warning("FLEX/%s bundle missing: %s", comp, cdir)
            continue
        genes = pd.read_csv(cdir / "genes.tsv", header=None)[0].tolist()
        cells = pd.read_csv(cdir / "cells.tsv", header=None)[0].tolist()
        gene_idx = {g: i for i, g in enumerate(genes)}
        cell_idx = {c: i for i, c in enumerate(cells)}

        comp_ann = flex_ann[flex_ann["compartment"] == comp].copy()
        comp_ann = comp_ann[comp_ann["cell_id"].isin(cell_idx)]
        if len(comp_ann) == 0:
            log.warning("FLEX/%s: no annotated cells", comp)
            continue

        log.info("loading FLEX/%s MTX (%d cells × %d genes)", comp, len(cells), len(genes))
        with gzip.open(cdir / "counts.mtx.gz", "rb") as gz:
            mat = sp.csr_matrix(scipy.io.mmread(io.BytesIO(gz.read())))

        row_idx = np.array([cell_idx[c] for c in comp_ann["cell_id"]])
        labels = comp_ann["label_short"].values
        row_totals = np.asarray(mat[row_idx].sum(axis=1)).flatten()
        row_totals[row_totals == 0] = 1

        comp_markers = [g for g in all_markers if g in gene_idx]
        comp_cols = [gene_idx[g] for g in comp_markers]
        mat_sub = mat[np.ix_(row_idx, comp_cols)]

        for lbl in label_order:
            mask = labels == lbl
            if mask.sum() == 0:
                continue
            sub = mat_sub[mask].toarray().astype(float)
            totals = row_totals[mask, np.newaxis]
            sub_norm = np.log1p(sub / totals * 1e4)
            for gi, g in enumerate(comp_markers):
                means[lbl][g] = float(sub_norm[:, gi].mean())

    aligned = {lbl: np.array([means[lbl].get(g, 0.0) for g in all_markers])
               for lbl in label_order}
    return aligned, list(all_markers)


def _stats_pair(x: np.ndarray, y: np.ndarray) -> dict:
    """Pearson r + Spearman ρ + Spearman p, restricted to (x>0)|(y>0) points."""
    mask = (x > 0) | (y > 0)
    n = int(mask.sum())
    if n < 4:
        return {"r": float("nan"), "rho": float("nan"), "p": float("nan"), "n": n}
    r, _ = scipy.stats.pearsonr(x[mask], y[mask])
    rho, p = scipy.stats.spearmanr(x[mask], y[mask])
    return {"r": float(r), "rho": float(rho), "p": float(p), "n": n}


def main() -> int:
    args = parse_args()
    config_dir = str(args.project_root / "publication" / "config")
    paths = load_paths(config_dir=config_dir)
    aes = load_aesthetics(config_dir=config_dir)

    # Load canonical markers from the promoted fig3 dotplot features CSV
    features_path = Path(resolve_path("fig3_dotplot_features", paths))
    log.info("loading canonical markers: %s", features_path)
    markers, label_order_full, full_names, label_comp, all_markers = (
        load_canonical_markers(features_path)
    )
    log.info("loaded %d cell types, %d unique features",
             len(label_order_full), len(all_markers))

    ann_path = Path(resolve_path("xenium_compartment_dir", paths)) / "all_compartments.csv"
    log.info("loading annotations: %s", ann_path)
    ann = pd.read_csv(ann_path)
    log.info("  rows: %d | columns: %s", len(ann), list(ann.columns))

    # is_artifact may be string or bool depending on writer; normalise
    if ann["is_artifact"].dtype == object:
        ann = ann[ann["is_artifact"].astype(str).str.lower() == "false"].copy()
    else:
        ann = ann[~ann["is_artifact"].astype(bool)].copy()
    log.info("non-artifact cells: %d", len(ann))

    label_order = label_order_full
    if args.test:
        label_order = label_order_full[:5]
        log.info("--test: limiting to %d cell types", len(label_order))

    nuc_bundle = Path(resolve_path("xenium_nuclear_bundle", paths))
    flex_root = Path(resolve_path("flex_compartment_bundle_root", paths))

    nuc_means, nuc_genes = load_xenium_means(nuc_bundle, ann, all_markers, label_order)
    flex_means, flex_genes = load_flex_means(flex_root, ann, all_markers, label_order)

    common_genes = [g for g in all_markers if g in nuc_genes and g in flex_genes]
    nuc_gi = [nuc_genes.index(g) for g in common_genes]
    flex_gi = [flex_genes.index(g) for g in common_genes]
    log.info("genes in both layers: %d / %d", len(common_genes), len(all_markers))

    # ---- Theme + palette ----
    plt.rcParams.update(get_matplotlib_theme(aes))
    comp_palette = get_palette("compartment", aes)

    # ---- Per-cell-type stats ----
    stats_rows = []
    for lbl in label_order:
        nuc_v = nuc_means[lbl][nuc_gi]
        flex_v = flex_means[lbl][flex_gi]
        s = _stats_pair(nuc_v, flex_v)
        stats_rows.append({
            "label_short": lbl,
            "label_full": full_names.get(lbl, lbl),
            "compartment": label_comp[lbl],
            **s,
        })
    stats_df = pd.DataFrame(stats_rows)
    log.info("computed stats for %d cell types", len(stats_df))

    # ---- Layout: 4 rows × 6 cols per-cell-type grid ----
    # 24 cell types in 4×6 = 24 slots, no empty. Each ~1.2 × 1.2 in =>
    # 7.2 wide × 4.8 tall. Plus summary panel below.
    heat_dims = get_dimensions("heatmap_tile", aes)
    fig_w = float(heat_dims["width"])  # 7.2
    grid_h = 4.8
    summary_h = 3.5
    fig_h = grid_h + summary_h + 0.7  # spacing + bottom legend

    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[grid_h, summary_h + 0.7],
        hspace=0.20,
        left=0.06, right=0.98, top=0.97, bottom=0.04,
    )
    grid = outer[0, 0].subgridspec(4, 6, hspace=0.6, wspace=0.45)
    summary_block = outer[1, 0].subgridspec(
        nrows=2, ncols=1, height_ratios=[summary_h, 0.5], hspace=0.10
    )
    ax_summary = fig.add_subplot(summary_block[0, 0])
    ax_summary_legend = fig.add_subplot(summary_block[1, 0])
    ax_summary_legend.axis("off")
    log.info("figure: %.2fin x %.2fin", fig_w, fig_h)

    # ---- Per-cell-type panels: FLEX vs Nuclear ----
    for idx, lbl in enumerate(label_order):
        r, c = divmod(idx, 6)
        ax = fig.add_subplot(grid[r, c])
        comp = label_comp[lbl]
        color = comp_palette.get(comp, "#666666")
        canonical = set(markers[lbl])
        nuc_v = nuc_means[lbl][nuc_gi]
        flex_v = flex_means[lbl][flex_gi]

        bg_mask = np.array([g not in canonical for g in common_genes])
        ax.scatter(nuc_v[bg_mask], flex_v[bg_mask],
                   s=4, c="lightgrey", linewidths=0, alpha=0.7, zorder=1,
                   rasterized=True)
        ax.scatter(nuc_v[~bg_mask], flex_v[~bg_mask],
                   s=10, c=color, linewidths=0, alpha=0.9, zorder=2,
                   rasterized=True)

        vmax = max(nuc_v.max(), flex_v.max(), 0.1)
        ax.plot([0, vmax], [0, vmax], lw=0.4, color="lightgrey", zorder=0)

        s = stats_df[stats_df["label_short"] == lbl].iloc[0]
        if not np.isnan(s["r"]):
            ax.text(0.03, 0.95,
                    f"r={s['r']:.2f}\np={s['p']:.1e}",
                    transform=ax.transAxes, fontsize=4.5,
                    ha="left", va="top")

        ax.set_xlabel("Nuclear", fontsize=5)
        if c == 0:
            ax.set_ylabel("FLEX", fontsize=5)
        ax.tick_params(labelsize=4, length=2, pad=1)
        ax.set_title(full_names.get(lbl, lbl), fontsize=5.5,
                     fontweight="bold", pad=2)
        for sp_loc in ax.spines.values():
            sp_loc.set_linewidth(0.4)

    # ---- Summary panel: all cell types overlaid, coloured by compartment ----
    all_r = []
    for lbl in label_order:
        comp = label_comp[lbl]
        color = comp_palette.get(comp, "#666666")
        nuc_v = nuc_means[lbl][nuc_gi]
        flex_v = flex_means[lbl][flex_gi]
        ax_summary.scatter(nuc_v, flex_v, s=4, c=color,
                           alpha=0.5, linewidths=0, rasterized=True)
        s = stats_df[stats_df["label_short"] == lbl].iloc[0]
        if not np.isnan(s["r"]):
            all_r.append(s["r"])

    summary_vmax = max(
        max(nuc_means[lbl][nuc_gi].max() for lbl in label_order),
        max(flex_means[lbl][flex_gi].max() for lbl in label_order),
        0.1,
    )
    ax_summary.plot([0, summary_vmax], [0, summary_vmax],
                    lw=0.5, color="lightgrey", zorder=0)
    mean_r = float(np.mean(all_r)) if all_r else 0.0
    ax_summary.text(0.03, 0.95, f"mean r = {mean_r:.2f}",
                    transform=ax_summary.transAxes,
                    fontsize=7, ha="left", va="top")
    ax_summary.set_xlabel("Xenium nuclear (mean log-norm)", fontsize=7)
    ax_summary.set_ylabel("FLEX (mean log-norm)", fontsize=7)
    ax_summary.tick_params(labelsize=6, length=2, pad=1)
    for sp_loc in ax_summary.spines.values():
        sp_loc.set_linewidth(0.5)

    # Compartment legend
    import matplotlib.patches as mpatches  # local import to keep top tight
    leg_handles = [
        mpatches.Patch(color=comp_palette[c], label=c)
        for c in ("Epithelial", "Immune", "Stromal") if c in comp_palette
    ]
    ax_summary_legend.legend(
        handles=leg_handles, loc="center", ncol=len(leg_handles),
        fontsize=7, frameon=False, handlelength=1.2, handletextpad=0.5,
        columnspacing=2.0,
        title="Compartment", title_fontsize=8,
    )

    # ---- Save ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / "supp7_3_layer_correlation.pdf"
    summary_pdf_path = args.out_dir / "supp7_3_layer_correlation_summary.pdf"
    stats_csv_path = args.out_dir / "supp7_3_layer_correlation_stats.csv"

    _raster_dpi = int(aes.get("rendering", {}).get("dpi_umap_raster", 600))
    fig.savefig(pdf_path, bbox_inches="tight", dpi=_raster_dpi)
    plt.close(fig)
    log.info("wrote %s", pdf_path)

    # Summary-only PDF: re-render just the summary panel for Illustrator
    # composition flexibility. Reuses the data already computed.
    fig2 = plt.figure(figsize=(3.5, 3.5))
    ax2 = fig2.add_subplot(111)
    for lbl in label_order:
        comp = label_comp[lbl]
        color = comp_palette.get(comp, "#666666")
        nuc_v = nuc_means[lbl][nuc_gi]
        flex_v = flex_means[lbl][flex_gi]
        ax2.scatter(nuc_v, flex_v, s=5, c=color, alpha=0.5,
                    linewidths=0, rasterized=True)
    ax2.plot([0, summary_vmax], [0, summary_vmax], lw=0.5, color="lightgrey")
    ax2.text(0.03, 0.95, f"mean r = {mean_r:.2f}",
             transform=ax2.transAxes, fontsize=7, ha="left", va="top")
    ax2.set_xlabel("Xenium nuclear (mean log-norm)", fontsize=7)
    ax2.set_ylabel("FLEX (mean log-norm)", fontsize=7)
    ax2.tick_params(labelsize=6, length=2)
    for sp_loc in ax2.spines.values():
        sp_loc.set_linewidth(0.5)
    ax2.legend(handles=leg_handles, loc="upper left", frameon=False,
               fontsize=6, handlelength=1.0, handletextpad=0.4)
    fig2.savefig(summary_pdf_path, bbox_inches="tight", dpi=_raster_dpi)
    plt.close(fig2)
    log.info("wrote %s", summary_pdf_path)

    stats_df.to_csv(stats_csv_path, index=False)
    log.info("wrote %s (%d cell types)", stats_csv_path, len(stats_df))

    return 0


if __name__ == "__main__":
    sys.exit(main())
