"""
Render S7 LASP-override duplex panels — standalone panel a + panel b producers.

Panel a: 6-marker feature plot grid on Epithelial scANVI UMAP (977K cells).
Panel b: Stacked violin (6 markers x 8 res-5.0 sub-clusters, grouped by majority L2).

This module both produces standalone PDFs when invoked directly, and exposes
constants + renderer functions consumed by `render_s7_composite.py` (which
assembles the locked s1.5 composite panel from the same per-panel logic).

Substrate (gzipped CSV, resolved via --substrate-dir):
  s7a_epi_umap_feature_data.csv.gz  (977K Epi cells x 6 markers + UMAP + L2_label)
  s7b_lasp_violin_long.csv.gz       (LASP cells, long-form for violin)

Output (under --out-dir):
  s_ihbca_lasp_feature_plots.{pdf,png,panel.yaml}     (panel a)
  s_ihbca_lasp_subcluster_violins.{pdf,png,panel.yaml} (panel b)
"""
from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle


# ----------------------------------------------------------------------------
# Locate project root + load aesthetics framework at module import time
# ----------------------------------------------------------------------------
# Module-level config load lets sibling modules (composite, panel_b_heatmap)
# consume the lasp_override palette via attribute access on this module
# (e.g. `base.BLOCK_TINT`) without each having to re-load aesthetics.yaml.

def _find_project_root() -> Path:
    here = Path(__file__).resolve().parent
    for p in [here] + list(here.parents):
        if (p / "publication" / "config" / "aesthetics.yaml").exists():
            return p
    raise FileNotFoundError("Could not locate publication/config/aesthetics.yaml")


_PROJECT_ROOT = _find_project_root()
sys.path.insert(0, str(_PROJECT_ROOT / "publication" / "config"))

from load_aesthetics import (  # noqa: E402
    get_matplotlib_theme,
    get_palette,
    load_aesthetics,
)

_CONFIG = load_aesthetics(config_dir=str(_PROJECT_ROOT / "publication" / "config"))
_LASP_PALETTE = get_palette("lasp_override", _CONFIG)

# Panel constants derived from the lasp_override palette in aesthetics.yaml.
# Composite + panel_b_heatmap sibling modules read these as `base.<name>`.
BLOCK_TINT = {k: _LASP_PALETTE[k]["tint"] for k in ("LASP-major", "LASP-basal", "LASP-KIT")}
BLOCK_FACE = {k: _LASP_PALETTE[k]["face"] for k in ("LASP-major", "LASP-basal", "LASP-KIT")}
BACKGROUND_GRAY = _LASP_PALETTE["background_gray"]
ZERO_GRAY = _LASP_PALETTE["zero_gray"]


# ----------------------------------------------------------------------------
# Invariant constants (locked per render plan coordination/plans/s7_lasp_render_20260516.md)
# ----------------------------------------------------------------------------

MARKERS = ["PTN", "KRT5", "KRT14", "KIT", "KRT8", "KRT18"]

# Sub-cluster columns for Panel b violin (subsampled from substrate per Option A:
# 5 representative LASP-major + 2 LASP-basal + 1 LASP-KIT; n is cell count).
PANEL_B_CLUSTERS = [
    (1,  "LASP-major", 23391),
    (2,  "LASP-major", 22318),
    (8,  "LASP-major", 18004),
    (9,  "LASP-major", 17785),
    (13, "LASP-major", 16773),
    (6,  "LASP-basal", 18649),
    (84, "LASP-basal",  3995),
    (62, "LASP-KIT",    6919),
]


def apply_theme() -> None:
    """Apply the Nature-spec matplotlib rcParams from the project aesthetics framework."""
    plt.rcParams.update(get_matplotlib_theme(_CONFIG))


def make_zero_gray_viridis():
    """Custom cmap: gray at value=0, transitioning into viridis for >0.

    Anchors zero to ZERO_GRAY so unexpressed LASP cells are visually distinct
    from dark-purple low-viridis (avoids confusion with background-gray cells).
    """
    base = mpl.colormaps["viridis"]
    n = 256
    colors = [ZERO_GRAY] + [base(i / (n - 1)) for i in range(int(0.05 * n), n)]
    return LinearSegmentedColormap.from_list("zero_gray_viridis", colors, N=n)


# ----------------------------------------------------------------------------
# Panel a — feature plot grid
# ----------------------------------------------------------------------------

def render_panel_a(df: pd.DataFrame, out_pdf: Path, out_png: Path) -> dict:
    """6-marker feature plot grid on Epi UMAP (rasterized scatter, 600 DPI)."""
    cmap = make_zero_gray_viridis()

    is_lasp = df["L2_label"].str.startswith("LASP")
    bg = df.loc[~is_lasp]
    fg = df.loc[is_lasp]

    print(f"Panel a: background n={len(bg):,}, LASP foreground n={len(fg):,}", flush=True)

    fig_w_in = 7.2  # Nature double-column
    fig_h_in = 4.9  # 2 rows x 3 cols, slightly compressed
    fig, axes = plt.subplots(
        2, 3, figsize=(fig_w_in, fig_h_in),
        gridspec_kw={"hspace": 0.12, "wspace": 0.08},
    )
    axes = axes.ravel()

    xmin, xmax = df["UMAP_X"].min(), df["UMAP_X"].max()
    ymin, ymax = df["UMAP_Y"].min(), df["UMAP_Y"].max()
    pad_x = 0.02 * (xmax - xmin)
    pad_y = 0.02 * (ymax - ymin)

    bg_xy = bg[["UMAP_X", "UMAP_Y"]].to_numpy()
    fg_xy = fg[["UMAP_X", "UMAP_Y"]].to_numpy()

    for i, marker in enumerate(MARKERS):
        ax = axes[i]
        col = f"log1p_{marker}"
        vals = fg[col].to_numpy()
        # Per-marker independent scale; vmax = 99th percentile (on LASP only)
        vmax = float(np.quantile(vals, 0.99)) if len(vals) else 1.0
        if vmax <= 0:
            vmax = float(vals.max()) if len(vals) and vals.max() > 0 else 1.0

        # Background layer (single gray)
        ax.scatter(
            bg_xy[:, 0], bg_xy[:, 1],
            s=0.05, c=BACKGROUND_GRAY, alpha=0.15,
            linewidths=0, marker=".",
            rasterized=True,
        )
        # Foreground layer (LASP cells colored by expression).
        # Plot zero-expression cells first (in gray) then non-zero cells colored,
        # so colored points overlay zeros where they coexist.
        zero_mask = vals <= 0
        if zero_mask.any():
            ax.scatter(
                fg_xy[zero_mask, 0], fg_xy[zero_mask, 1],
                s=0.08, c=ZERO_GRAY, alpha=0.6,
                linewidths=0, marker=".",
                rasterized=True,
            )
        nz_mask = ~zero_mask
        if nz_mask.any():
            # Sort by value ascending so high-expression points draw on top
            order = np.argsort(vals[nz_mask])
            xs = fg_xy[nz_mask, 0][order]
            ys = fg_xy[nz_mask, 1][order]
            cs = vals[nz_mask][order]
            sc = ax.scatter(
                xs, ys,
                s=0.10, c=cs, cmap=cmap, vmin=0, vmax=vmax,
                alpha=0.9, linewidths=0, marker=".",
                rasterized=True,
            )
        else:
            sc = None

        ax.set_xlim(xmin - pad_x, xmax + pad_x)
        ax.set_ylim(ymin - pad_y, ymax + pad_y)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Minimal UMAP1/UMAP2 arrow labels (lower-left of each sub-panel)
        ax.annotate("", xy=(0.18, 0.02), xytext=(0.02, 0.02),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=0.4, color="black"))
        ax.annotate("", xy=(0.02, 0.18), xytext=(0.02, 0.02),
                    xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="->", lw=0.4, color="black"))
        ax.text(0.20, 0.02, "UMAP1", transform=ax.transAxes,
                fontsize=5, ha="left", va="bottom")
        ax.text(0.02, 0.20, "UMAP2", transform=ax.transAxes,
                fontsize=5, ha="left", va="bottom", rotation=90)

        # Marker symbol in upper-right (7pt bold italic)
        ax.text(0.98, 0.97, marker, transform=ax.transAxes,
                fontsize=7, fontweight="bold", style="italic",
                ha="right", va="top")

        # Per-marker colorbar (compact, inset to lower-right)
        if sc is not None:
            cax = ax.inset_axes([0.78, 0.10, 0.18, 0.025])
            cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
            cb.set_ticks([0, vmax])
            cb.set_ticklabels(["0", f"{vmax:.1f}"])
            cb.ax.tick_params(labelsize=5, pad=1, length=2, width=0.4)
            cb.outline.set_linewidth(0.4)

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {
        "bg_n": int(len(bg)),
        "fg_n": int(len(fg)),
        "markers": MARKERS,
    }


# ----------------------------------------------------------------------------
# Panel b — stacked violins
# ----------------------------------------------------------------------------

def render_panel_b(df_long: pd.DataFrame, out_pdf: Path, out_png: Path) -> dict:
    """6 markers x 8 sub-clusters stacked violins, column-blocked by majority L2."""
    cluster_ids = [c[0] for c in PANEL_B_CLUSTERS]
    cluster_to_block = {c[0]: c[1] for c in PANEL_B_CLUSTERS}

    # Filter long-form to the 8 selected sub-clusters
    df_long = df_long[df_long["leiden_5.0"].isin(cluster_ids)].copy()
    print(f"Panel b: filtered to {len(df_long):,} (cell × marker) rows across 8 clusters", flush=True)

    # Truncate to 99th percentile per marker
    caps = {}
    for marker in MARKERS:
        sub = df_long.loc[df_long["marker_gene"] == marker, "log1p_expression"]
        cap = float(np.quantile(sub, 0.99))
        caps[marker] = max(cap, 0.1)
    df_long["log1p_expression_trunc"] = df_long.apply(
        lambda r: min(r["log1p_expression"], caps[r["marker_gene"]]), axis=1
    )

    # Build positional index for the 8 columns (in PANEL_B_CLUSTERS order)
    col_pos = {cid: i for i, cid in enumerate(cluster_ids)}

    fig_w_in = 3.5  # Nature single-column
    fig_h_in = 6.8  # tall for 6 rows
    fig, axes = plt.subplots(
        len(MARKERS), 1, figsize=(fig_w_in, fig_h_in), sharex=True,
        gridspec_kw={"hspace": 0.20},
    )

    # Block boundaries (in positional units): LASP-major cols 0-4, basal 5-6, KIT 7
    block_spans = [
        ("LASP-major", -0.5, 4.5, 5),
        ("LASP-basal", 4.5, 6.5, 2),
        ("LASP-KIT",   6.5, 7.5, 1),
    ]

    for ri, marker in enumerate(MARKERS):
        ax = axes[ri]
        data_per_col = []
        positions = []
        for cid in cluster_ids:
            sub = df_long[
                (df_long["leiden_5.0"] == cid) & (df_long["marker_gene"] == marker)
            ]["log1p_expression_trunc"].to_numpy()
            data_per_col.append(sub if len(sub) > 0 else np.array([0.0]))
            positions.append(col_pos[cid])

        # Background tint shading for each block (behind violins)
        ymax = caps[marker] * 1.08
        for block_name, x0, x1, _n in block_spans:
            ax.add_patch(Rectangle(
                (x0, -0.05), x1 - x0, ymax + 0.10,
                facecolor=BLOCK_TINT[block_name], edgecolor="none", zorder=0,
            ))

        # Violin
        parts = ax.violinplot(
            data_per_col, positions=positions,
            widths=0.85, showmeans=False, showmedians=True, showextrema=False,
        )
        # Style violin bodies — color by block (face per block_face in palette)
        for body, cid in zip(parts["bodies"], cluster_ids):
            block = cluster_to_block[cid]
            body.set_facecolor(BLOCK_FACE[block])
            body.set_edgecolor("black")
            body.set_linewidth(0.3)
            body.set_alpha(0.85)
        if "cmedians" in parts:
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(0.6)

        # Vertical separators between blocks
        for x in (4.5, 6.5):
            ax.axvline(x, color="black", linewidth=0.4, zorder=2)

        ax.set_xlim(-0.5, 7.5)
        ax.set_ylim(0, ymax)
        ax.set_yticks([0, caps[marker]])
        ax.set_yticklabels(["0", f"{caps[marker]:.1f}"])
        ax.tick_params(axis="y", labelsize=5, length=2, pad=1, width=0.4)
        ax.tick_params(axis="x", length=0)
        for spine_name in ("top", "right"):
            ax.spines[spine_name].set_visible(False)
        ax.spines["left"].set_linewidth(0.4)
        ax.spines["bottom"].set_linewidth(0.4)

        # Row label: marker symbol (italic, 7pt) outside the left spine
        ax.set_ylabel(marker, fontsize=7, fontstyle="italic", rotation=0,
                      ha="right", va="center", labelpad=8)

        # Top-of-figure block tags (only on the top row)
        if ri == 0:
            for block_name, x0, x1, n in block_spans:
                ax.text(
                    (x0 + x1) / 2, ymax * 1.10,
                    f"{block_name} (n={n})",
                    ha="center", va="bottom",
                    fontsize=6, fontweight="bold",
                    clip_on=False,
                )
            # Extend y-limit so the labels render above the patch
            ax.set_ylim(0, ymax * 1.18)

    # Bottom-row x-tick labels: sub-cluster IDs
    axes[-1].set_xticks(list(range(len(cluster_ids))))
    axes[-1].set_xticklabels([str(c) for c in cluster_ids], fontsize=6)
    axes[-1].set_xlabel("Sub-cluster ID (Epi scANVI leiden res-5.0)",
                        fontsize=6, labelpad=4)

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {
        "n_rows": len(MARKERS),
        "n_cols": len(cluster_ids),
        "caps_99p": caps,
        "cluster_ids": cluster_ids,
    }


# ----------------------------------------------------------------------------
# Panel YAML writers
# ----------------------------------------------------------------------------

PANEL_A_YAML = """panel_id: s_ihbca_lasp_feature_plots
description: |
  LASP-override evidence panel — 6-marker feature plot grid on the Epithelial scANVI
  UMAP (977,541 cells). PTN (override anchor for LASP-basal), KRT5/KRT14 (basal
  keratins, intermediate in LASP-basal), KIT (LASP-KIT anchor), and KRT8/KRT18
  (luminal keratins retained across all LASP types). LASP cells (LASP-major,
  LASP-basal, LASP-KIT) shown in foreground with viridis-on-gray expression
  scale; other Epithelial cells (LHS, BMYO, Lactocyte) shown as faded gray
  background.
source_data:
  - publication/figures/data/fig1/lasp_override/s7a_epi_umap_feature_data.csv.gz
tags: [supplemental, fig1, annotation_method, lasp_override, feature_plot, umap]
prose_token: "[SUPP:s_ihbca_lasp_feature_plots]"
nature_spec: single_panel, double_column, vector_pdf, 600_dpi_rasterized
palette: lasp_override (publication/config/aesthetics.yaml)
parent_figure: Fig 1
"""

PANEL_B_YAML = """panel_id: s_ihbca_lasp_subcluster_violins
description: |
  LASP-override evidence panel — stacked violin plot of log1p expression for 6
  markers (PTN, KRT5, KRT14, KIT, KRT8, KRT18) across 8 representative Epithelial
  scANVI leiden res-5.0 sub-clusters in the LASP region. Sub-clusters grouped by
  majority L2 label: 5 representative LASP-major clusters (largest by cell count),
  2 LASP-basal clusters (all sub-clusters with n>=100 in the LASP-basal majority
  L2), 1 LASP-KIT cluster (sole sub-cluster with n>=100 in LASP-KIT). Violins
  truncated at marker-wise 99th percentile.
source_data:
  - publication/figures/data/fig1/lasp_override/s7b_lasp_violin_long.csv.gz
  - publication/figures/data/fig1/lasp_override/s7b_lasp_subcluster_majority.csv
tags: [supplemental, fig1, annotation_method, lasp_override, violin, leiden_res5]
prose_token: "[SUPP:s_ihbca_lasp_subcluster_violins]"
nature_spec: single_panel, single_column, vector_pdf
palette: lasp_override (publication/config/aesthetics.yaml)
parent_figure: Fig 1
"""


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--substrate-dir", required=True, type=Path,
                   help="Directory containing s7a_epi_umap_feature_data.csv.gz and s7b_lasp_violin_long.csv.gz.")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Destination directory for standalone panel PDFs/PNGs/YAMLs.")
    p.add_argument("--project-root", required=False, type=Path, default=None,
                   help="iHBCA_publication root (auto-detected from script location if omitted).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    apply_theme()

    sub_dir = args.substrate_dir

    # ----- Panel a -----
    print("Loading s7a feature data (977K rows × 11 cols, gzipped)...", flush=True)
    df_a = pd.read_csv(sub_dir / "s7a_epi_umap_feature_data.csv.gz")
    print(f"  shape: {df_a.shape}", flush=True)
    pdf_a = args.out_dir / "s_ihbca_lasp_feature_plots.pdf"
    png_a = args.out_dir / "s_ihbca_lasp_feature_plots.png"
    meta_a = render_panel_a(df_a, pdf_a, png_a)
    print(f"  wrote {pdf_a} and {png_a}", flush=True)
    del df_a
    gc.collect()

    # ----- Panel b -----
    print("Loading s7b violin long (gzipped)...", flush=True)
    df_b = pd.read_csv(sub_dir / "s7b_lasp_violin_long.csv.gz")
    print(f"  shape: {df_b.shape}", flush=True)
    pdf_b = args.out_dir / "s_ihbca_lasp_subcluster_violins.pdf"
    png_b = args.out_dir / "s_ihbca_lasp_subcluster_violins.png"
    meta_b = render_panel_b(df_b, pdf_b, png_b)
    print(f"  wrote {pdf_b} and {png_b}", flush=True)

    # ----- Panel YAML stubs -----
    (args.out_dir / "s_ihbca_lasp_feature_plots.panel.yaml").write_text(PANEL_A_YAML)
    (args.out_dir / "s_ihbca_lasp_subcluster_violins.panel.yaml").write_text(PANEL_B_YAML)

    print("Done.", flush=True)
    print(f"Panel a meta: {meta_a}", flush=True)
    print(f"Panel b meta: {meta_b}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
