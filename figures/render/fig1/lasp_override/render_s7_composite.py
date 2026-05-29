"""
Render s1.5 LASP-override composite supplementary figure.

Composes the two locked S7 panels into a single multi-panel PDF.

  +-----------------------------------------------------------+
  |  a   6-marker feature plot grid on Epi scANVI UMAP        |  ~120mm
  |      (2 rows x 3 cols: PTN/KRT5/KRT14, KIT/KRT8/KRT18)    |
  +-----------------------------------------------------------+
  |  b   res-5.0 stacked violins (6 markers x 8 sub-clusters) |   ~70mm
  +-----------------------------------------------------------+

Total: 183mm (Nature double column) x ~200mm.

Reuses constants, helpers, and per-panel logic from `render_s7.py` — the
plotting functions here adapt that logic to draw into supplied axes
inside a shared gridspec figure (so the composite can re-use the panel-b
heatmap geometry rather than duplicating it).

Substrate (gzipped CSV, resolved via --substrate-dir):
  s7a_epi_umap_feature_data.csv.gz  (977K Epi cells x 6 markers + UMAP + L2)
  s7b_lasp_violin_long.csv.gz       (LASP cells, long-form for violin/heatmap)

Output (under --out-dir):
  s1_5_lasp_basal_override.{pdf,png,panel.yaml}
"""
from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

# Ensure co-located `render_s7` is importable regardless of working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

# Sibling module: constants (MARKERS, BLOCK_TINT, BACKGROUND_GRAY, ZERO_GRAY) +
# helpers (apply_theme, make_zero_gray_viridis). Loaded from aesthetics.yaml
# at base's import time.
import render_s7 as base  # noqa: E402


MM_PER_IN = 25.4


# ----------------------------------------------------------------------------
# Panel renderers that draw INTO supplied axes (no figure creation)
# ----------------------------------------------------------------------------

def render_panel_a_into(fig: plt.Figure, axes: list[plt.Axes], df: pd.DataFrame) -> dict:
    """Render panel a (6-marker feature plots) into a supplied list of 6 axes.

    Mirrors `render_s7.render_panel_a` plotting logic, minus figure creation
    and savefig. The original function is left intact so the standalone script
    keeps working.
    """
    cmap = base.make_zero_gray_viridis()

    is_lasp = df["L2_label"].str.startswith("LASP")
    bg = df.loc[~is_lasp]
    fg = df.loc[is_lasp]

    xmin, xmax = df["UMAP_X"].min(), df["UMAP_X"].max()
    ymin, ymax = df["UMAP_Y"].min(), df["UMAP_Y"].max()
    pad_x = 0.02 * (xmax - xmin)
    pad_y = 0.02 * (ymax - ymin)

    bg_xy = bg[["UMAP_X", "UMAP_Y"]].to_numpy()
    fg_xy = fg[["UMAP_X", "UMAP_Y"]].to_numpy()

    for i, marker in enumerate(base.MARKERS):
        ax = axes[i]
        col = f"log1p_{marker}"
        vals = fg[col].to_numpy()
        vmax = float(np.quantile(vals, 0.99)) if len(vals) else 1.0
        if vmax <= 0:
            vmax = float(vals.max()) if len(vals) and vals.max() > 0 else 1.0

        ax.scatter(
            bg_xy[:, 0], bg_xy[:, 1],
            s=0.05, c=base.BACKGROUND_GRAY, alpha=0.15,
            linewidths=0, marker=".",
            rasterized=True,
        )
        zero_mask = vals <= 0
        if zero_mask.any():
            ax.scatter(
                fg_xy[zero_mask, 0], fg_xy[zero_mask, 1],
                s=0.08, c=base.ZERO_GRAY, alpha=0.6,
                linewidths=0, marker=".",
                rasterized=True,
            )
        nz_mask = ~zero_mask
        sc = None
        if nz_mask.any():
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

        ax.set_xlim(xmin - pad_x, xmax + pad_x)
        ax.set_ylim(ymin - pad_y, ymax + pad_y)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_visible(False)

        # UMAP1/UMAP2 arrows (lower-left of each sub-panel)
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

        # Marker symbol — upper-LEFT in the composite (the upper-right position
        # collides with LASP cells in the upper Epi UMAP region; spec explicitly
        # allows relocation).
        ax.text(0.02, 0.97, marker, transform=ax.transAxes,
                fontsize=7, fontweight="bold", style="italic",
                ha="left", va="top",
                bbox=dict(facecolor="white", edgecolor="none",
                          alpha=0.75, pad=1.0))

        # Per-marker colorbar (compact, inset to lower-right)
        if sc is not None:
            cax = ax.inset_axes([0.78, 0.10, 0.18, 0.025])
            cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
            cb.set_ticks([0, vmax])
            cb.set_ticklabels(["0", f"{vmax:.1f}"])
            cb.ax.tick_params(labelsize=5, pad=1, length=2, width=0.4)
            cb.outline.set_linewidth(0.4)

    return {"bg_n": int(len(bg)), "fg_n": int(len(fg)), "markers": list(base.MARKERS)}


def render_panel_b_into(
    fig: plt.Figure,
    axes,
    df_long: pd.DataFrame,
    label_rotation: float = 45.0,
    label_fontsize: int = 6,
) -> dict:
    """Render panel b as a per-marker row-z-score heatmap of ALL 108 res-5.0
    LASP-region sub-clusters, **primary-blocked by majority_L2** and
    **secondary-ordered by combined LASP-marker signal within each block**.

    Rows = 6 markers (PTN, KRT5, KRT14, KIT, KRT8, KRT18) top-to-bottom.
    Columns = all leiden res-5.0 sub-clusters in the substrate, partitioned
    LEFT->RIGHT into three blocks (LASP-basal, LASP-major, LASP-KIT) — v5
    reordering that flanks the dominant LASP-major mass with the two
    override-target blocks. Within each block, columns are sorted by
    descending Σ(mean log1p) across the 6 markers.

    Cell value = mean log1p expression per (marker, cluster) after a per-marker
    99th-percentile cap on the input. Color scale = per-row z-score on viridis.
    Pale block-tint rectangles sit BEHIND the heatmap cells; thin black
    separators mark block boundaries; one bold block label sits above each
    column range.
    """
    MARKERS = base.MARKERS
    BLOCK_TINT = base.BLOCK_TINT
    BLOCK_ORDER = ("LASP-basal", "LASP-major", "LASP-KIT")

    # Accept either a single Axes (preferred) or a 1-element list/tuple.
    if isinstance(axes, (list, tuple)):
        if len(axes) != 1:
            raise ValueError(
                f"render_panel_b_into expects a single Axes, got {len(axes)}"
            )
        ax = axes[0]
    else:
        ax = axes

    # ---- Per-marker 99th-percentile cap, then mean per (marker, cluster) ----
    caps = {}
    for marker in MARKERS:
        sub = df_long.loc[df_long["marker_gene"] == marker, "log1p_expression"]
        caps[marker] = max(float(np.quantile(sub, 0.99)), 0.1)

    cap_series = df_long["marker_gene"].map(caps).astype(float)
    df_long = df_long.assign(
        expr_capped=np.minimum(
            df_long["log1p_expression"].to_numpy(), cap_series.to_numpy()
        )
    )

    mean_tbl = (
        df_long.groupby(["marker_gene", "leiden_5.0"], observed=True)["expr_capped"]
        .mean()
        .unstack("leiden_5.0")
    )
    mean_tbl = mean_tbl.reindex(MARKERS).fillna(0.0)

    # ---- Resolve cluster -> majority_L2 (each cluster has one majority_L2) ----
    cluster_to_block = (
        df_long.groupby("leiden_5.0")["majority_L2"].first().to_dict()
    )

    # ---- Column ordering: primary block, secondary signal-descending ----
    col_signal = mean_tbl.sum(axis=0)
    ordered_cols: list = []
    block_ranges: dict[str, tuple[int, int]] = {}
    cursor = 0
    for block_name in BLOCK_ORDER:
        in_block = [
            c for c in mean_tbl.columns
            if cluster_to_block.get(c) == block_name
        ]
        in_block_sorted = sorted(
            in_block, key=lambda c: float(col_signal[c]), reverse=True
        )
        ordered_cols.extend(in_block_sorted)
        block_ranges[block_name] = (cursor, len(in_block_sorted))
        cursor += len(in_block_sorted)

    mean_tbl = mean_tbl[ordered_cols]
    n_cols = len(ordered_cols)
    n_rows = len(MARKERS)

    # ---- Row z-score ----
    M = mean_tbl.to_numpy(dtype=float)
    row_mean = M.mean(axis=1, keepdims=True)
    row_std = M.std(axis=1, keepdims=True)
    row_std = np.where(row_std < 1e-9, 1.0, row_std)
    Z = (M - row_mean) / row_std

    zabs = float(np.nanmax(np.abs(Z))) if np.isfinite(Z).any() else 1.0
    zabs = max(zabs, 0.5)
    vmin, vmax = -zabs, zabs

    # ---- Coord system: image extent (-0.5..n_cols-0.5) x (-0.5..n_rows-0.5).
    # Block tints are drawn BEFORE the imshow with zorder=0 so they sit below.
    ax.set_xlim(-0.5, n_cols - 0.5)
    ax.set_ylim(n_rows - 0.5, -0.5)  # origin upper

    for block_name in BLOCK_ORDER:
        start, count = block_ranges[block_name]
        if count == 0:
            continue
        ax.add_patch(Rectangle(
            (start - 0.5, -0.5), count, n_rows,
            facecolor=BLOCK_TINT.get(block_name, "#DDDDDD"),
            edgecolor="none", linewidth=0.0, alpha=0.15, zorder=0,
        ))

    # ---- Heatmap (zorder above tints) ----
    im = ax.imshow(
        Z,
        aspect="auto",
        cmap="viridis",
        vmin=vmin, vmax=vmax,
        interpolation="nearest",
        origin="upper",
        extent=(-0.5, n_cols - 0.5, n_rows - 0.5, -0.5),
        zorder=2,
    )

    # ---- Block separators: thin vertical black lines at boundaries ----
    cum = 0
    for block_name in BLOCK_ORDER[:-1]:
        cum += block_ranges[block_name][1]
        ax.axvline(
            x=cum - 0.5,
            color="black", linewidth=0.6,
            ymin=0.0, ymax=1.0,
            zorder=3, clip_on=True,
        )

    # ---- Block labels above the top row ----
    label_y = -0.85
    for block_name in BLOCK_ORDER:
        start, count = block_ranges[block_name]
        if count == 0:
            continue
        center_x = start + (count - 1) / 2.0
        ax.text(
            center_x, label_y,
            f"{block_name} (n={count})",
            ha="center", va="bottom",
            fontsize=6, fontweight="bold",
            clip_on=False,
        )

    # ---- Row labels (italic markers) ----
    ax.set_yticks(list(range(n_rows)))
    ax.set_yticklabels(MARKERS, fontsize=7, fontstyle="italic")
    ax.tick_params(axis="y", length=0, pad=2)

    # Hide column tick labels (108 IDs unreadable).
    ax.set_xticks([])
    ax.tick_params(axis="x", length=0)

    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.4)
    ax.spines["bottom"].set_linewidth(0.4)

    ax.set_xlabel(
        f"Epi scANVI res-5.0 sub-clusters, blocked by majority L2, "
        f"ordered by combined LASP-marker signal within block (n = {n_cols})",
        fontsize=6, labelpad=3,
    )

    # ---- Compact horizontal colorbar (row-z scale, bottom-right inset) ----
    cax = ax.inset_axes([0.78, -0.28, 0.20, 0.045])
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_ticks([vmin, 0.0, vmax])
    cb.set_ticklabels([f"{vmin:.1f}", "0", f"{vmax:.1f}"])
    cb.ax.tick_params(labelsize=5, pad=1, length=2, width=0.4)
    cb.outline.set_linewidth(0.4)
    cax.set_title("row z-score", fontsize=5, pad=1)

    block_counts = {name: block_ranges[name][1] for name in BLOCK_ORDER}

    return {
        "n_rows": n_rows,
        "n_cols": n_cols,
        "caps_99p": caps,
        "cluster_ids_ordered": [int(c) for c in ordered_cols],
        "block_counts": block_counts,
        "block_ranges": {k: list(v) for k, v in block_ranges.items()},
        "row_zscore_range": [vmin, vmax],
        "label_rotation": label_rotation,
        "label_fontsize": label_fontsize,
        "layout": "single_axes_heatmap_block_then_signal_v5",
    }


# ----------------------------------------------------------------------------
# Composite figure builder
# ----------------------------------------------------------------------------

# Locked composite dimensions (mm). Outer panel heights are gridspec-row mm
# allocations; matplotlib's tight bbox can trim outer whitespace, so we set
# total figure height generously and rely on relative row ratios.
PANEL_A_HEIGHT_MM = 120.0
PANEL_B_HEIGHT_MM = 70.0
COMPOSITE_WIDTH_MM = 183.0
COMPOSITE_HEIGHT_MM = PANEL_A_HEIGHT_MM + PANEL_B_HEIGHT_MM + 10.0  # +10mm inter-row gutter


def build_composite(
    df_a: pd.DataFrame,
    df_b: pd.DataFrame,
    out_pdf: Path,
    out_png: Path,
) -> dict:
    fig_w_in = COMPOSITE_WIDTH_MM / MM_PER_IN
    fig_h_in = COMPOSITE_HEIGHT_MM / MM_PER_IN

    fig = plt.figure(figsize=(fig_w_in, fig_h_in))

    # Outer grid: 2 rows (panel a top, panel b bottom)
    outer = GridSpec(
        2, 1,
        height_ratios=[PANEL_A_HEIGHT_MM, PANEL_B_HEIGHT_MM],
        figure=fig,
        hspace=0.18,
        left=0.05, right=0.985, top=0.985, bottom=0.045,
    )

    # Panel a sub-grid: 2 rows x 3 cols of feature plots
    a_grid = outer[0].subgridspec(2, 3, hspace=0.10, wspace=0.06)
    a_axes = [fig.add_subplot(a_grid[r, c]) for r in range(2) for c in range(3)]

    # Panel b: single Axes covering the 108-cluster heatmap.
    b_grid = outer[1].subgridspec(1, 1)
    b_ax = fig.add_subplot(b_grid[0, 0])

    meta_a = render_panel_a_into(fig, a_axes, df_a)
    meta_b = render_panel_b_into(
        fig, b_ax, df_b,
        label_rotation=45.0,
        label_fontsize=6,
    )

    fig.savefig(out_pdf, dpi=600, bbox_inches="tight", format="pdf")
    fig.savefig(out_png, dpi=300, bbox_inches="tight", format="png")
    plt.close(fig)

    return {
        "width_mm": COMPOSITE_WIDTH_MM,
        "height_mm": COMPOSITE_HEIGHT_MM,
        "panel_a": meta_a,
        "panel_b": meta_b,
    }


# ----------------------------------------------------------------------------
# Panel YAML
# ----------------------------------------------------------------------------

COMPOSITE_PANEL_YAML = """panel_id: s_ihbca_lasp_composite
output_slug: s1_5_lasp_basal_override
description: |
  S7 LASP-override evidence composite — (a) 6-marker feature plot grid on the
  Epithelial scANVI UMAP (977,541 cells; PTN, KRT5, KRT14, KIT, KRT8, KRT18 with
  per-marker independent 99th-percentile vmax and zero-anchored viridis so
  unexpressed LASP cells render light gray distinct from background gray);
  (b) Per-marker row-z-score heatmap of the same 6 markers across ALL 108
  LASP-region Epi scANVI res-5.0 sub-clusters. Columns are primary-blocked by
  majority L2 (LASP-basal | LASP-major | LASP-KIT, left-to-right; v5 ordering
  flanks the dominant LASP-major mass with the two override-target blocks) and
  secondary-ordered by descending Σ(mean log1p) across the 6 markers within
  each block. Pale BLOCK_TINT background rectangles run behind each block at
  low alpha; thin black vertical lines mark block boundaries; one bold label
  per block sits above the heatmap.
source_data:
  - publication/figures/data/fig1/lasp_override/s7a_epi_umap_feature_data.csv.gz
  - publication/figures/data/fig1/lasp_override/s7b_lasp_subcluster_majority.csv
  - publication/figures/data/fig1/lasp_override/s7b_lasp_violin_long.csv.gz
tags: [supplemental, fig1, annotation_method, lasp_override, composite, feature_plot, violin, umap]
prose_token: "[SUPP:s_ihbca_lasp_composite]"
nature_spec: composite, double_column, vector_pdf, 600_dpi_rasterized
palette: lasp_override (publication/config/aesthetics.yaml)
parent_figure: Fig 1
"""


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--substrate-dir", required=True, type=Path,
                   help="Directory containing s7a_epi_umap_feature_data.csv.gz + s7b_lasp_violin_long.csv.gz.")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Destination directory for s1_5_lasp_basal_override.{pdf,png,panel.yaml}.")
    p.add_argument("--project-root", required=False, type=Path, default=None,
                   help="iHBCA_publication root (auto-detected via render_s7 module if omitted).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    base.apply_theme()

    sub_dir = args.substrate_dir

    print("Loading s7a feature data (gzipped)...", flush=True)
    df_a = pd.read_csv(sub_dir / "s7a_epi_umap_feature_data.csv.gz")
    print(f"  shape: {df_a.shape}", flush=True)

    print("Loading s7b violin long (gzipped)...", flush=True)
    df_b = pd.read_csv(sub_dir / "s7b_lasp_violin_long.csv.gz")
    print(f"  shape: {df_b.shape}", flush=True)

    pdf_path = args.out_dir / "s1_5_lasp_basal_override.pdf"
    png_path = args.out_dir / "s1_5_lasp_basal_override.png"
    print("Building composite figure...", flush=True)
    meta = build_composite(df_a, df_b, pdf_path, png_path)
    print(f"  wrote {pdf_path}", flush=True)
    print(f"  wrote {png_path}", flush=True)

    yaml_path = args.out_dir / "s1_5_lasp_basal_override.panel.yaml"
    yaml_path.write_text(COMPOSITE_PANEL_YAML)
    print(f"  wrote {yaml_path}", flush=True)

    del df_a, df_b
    gc.collect()

    print("Done.", flush=True)
    print(f"Composite meta: {meta}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
