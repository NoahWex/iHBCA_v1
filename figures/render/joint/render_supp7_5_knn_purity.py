"""
render_supp7_5_knn_purity.py — Fig 3 Supp 7-5.

Two-panel composite documenting the joint Concord integration kNN purity
benchmark used to select the canonical joint config (joint_nuc_n100).

Panel a: mean kNN purity (k=30, L0.5) across 7 sweep configs (3 source
spaces × 2 latent dims + FLEX-only baseline). Horizontal dashed reference
at FLEX-only baseline. Winner highlighted.

Panel b: per-L0.5-cell-type mean kNN purity for the winner config
(joint_nuc_n100), 13 cell types color-coded by compartment.

Inputs:
    paths.joint_knn_purity            knn_purity_summary.csv (7 × 17)
    paths.joint_knn_purity_per_cell   knn_purity.csv (1.8M rows; per-cell
                                       purity needed to compute per-celltype
                                       means — summary file medians collapse
                                       to 1.0 for every type, mean shows the
                                       0.94–0.99 spread)

Outputs:
    supp7_5_knn_purity_composite.pdf  two-panel composite, 89mm width
    supp7_5_knn_purity_composite.png  raster companion
    supp7_5_knn_purity_data.csv       rendered values (long-format)

Pattern source:
    render_supp7_3_layer_correlation.py — load_paths / load_aesthetics
    scaffolding + matplotlib gridspec composition pattern.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
log = logging.getLogger("supp7_5")


# Config display order (joint sweep first, baseline last)
CONFIG_ORDER = ["nuc_n50", "nuc_n100", "whl_n50", "whl_n100",
                "cyto_n50", "cyto_n100", "flex_only"]

# L0.5 cell type → compartment (for panel b color coding)
L0P5_COMPARTMENT = {
    "Adipocyte": "Stromal",
    "BMYO":      "Epithelial",
    "B_cell":    "Immune",
    "FB":        "Stromal",
    "LASP":      "Epithelial",
    "LE":        "Stromal",
    "LHS":       "Epithelial",
    "Mast":      "Immune",
    "Myeloid":   "Immune",
    "PV":        "Stromal",
    "Plasma":    "Immune",
    "T_cell":    "Immune",
    "VE":        "Stromal",
}

# Display order within each compartment (Epi → Stromal → Immune)
L0P5_ORDER = ["LASP", "LHS", "BMYO",
              "FB", "PV", "Adipocyte", "VE", "LE",
              "T_cell", "B_cell", "Plasma", "Myeloid", "Mast"]

WINNER = "nuc_n100"
BASELINE = "flex_only"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    config_dir = str(args.project_root / "publication" / "config")
    paths = load_paths(config_dir=config_dir)
    aes = load_aesthetics(config_dir=config_dir)

    summary_path = Path(resolve_path("joint_knn_purity", paths))
    per_cell_path = Path(resolve_path("joint_knn_purity_per_cell", paths))

    log.info("loading summary: %s", summary_path)
    summary = pd.read_csv(summary_path).set_index("config")

    log.info("loading per-cell purity (~1.8M rows): %s", per_cell_path)
    per_cell = pd.read_csv(per_cell_path)
    log.info("  rows: %d | configs: %s", len(per_cell),
             sorted(per_cell["config"].unique()))

    # Per-config mean (panel a) — read from promoted summary
    panel_a = summary.loc[CONFIG_ORDER, "mean"].rename("mean_purity")
    baseline_value = float(panel_a.loc[BASELINE])

    # Per-config × per-celltype mean (panel b) — aggregate from per-cell
    log.info("aggregating per-celltype means")
    per_celltype = (
        per_cell.groupby(["config", "l0p5"])["purity"]
        .mean().reset_index().rename(columns={"purity": "mean_purity"})
    )
    panel_b = (
        per_celltype[per_celltype["config"] == WINNER]
        .set_index("l0p5")
        .loc[L0P5_ORDER, "mean_purity"]
    )

    # Theme
    plt.rcParams.update(get_matplotlib_theme(aes))
    comp_palette = get_palette("compartment", aes)
    winner_color = comp_palette.get("Epithelial", "#7570B3")  # accent
    neutral_color = "#666666"

    # Figure: 89mm single column (3.5in wide), two stacked panels
    single_w_in = float(aes["dimensions"]["single_column_mm"]) / 25.4
    panel_a_h = 1.9
    panel_b_h = 2.1
    fig_h = panel_a_h + panel_b_h + 0.6  # + spacing

    fig = plt.figure(figsize=(single_w_in, fig_h))
    gs = fig.add_gridspec(
        nrows=2, ncols=1,
        height_ratios=[panel_a_h, panel_b_h],
        hspace=0.55,
        left=0.16, right=0.97, top=0.97, bottom=0.16,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])

    # ---- Panel a: per-config mean ----
    bar_colors_a = [winner_color if c == WINNER else neutral_color
                    for c in CONFIG_ORDER]
    x_a = np.arange(len(CONFIG_ORDER))
    ax_a.bar(x_a, panel_a.values, color=bar_colors_a,
             edgecolor="none", width=0.7)
    ax_a.axhline(baseline_value, ls="--", lw=0.5, color=neutral_color, zorder=0)
    ax_a.text(len(CONFIG_ORDER) - 0.5, baseline_value, " baseline",
              va="center", ha="left", fontsize=5, color=neutral_color)

    ax_a.set_xticks(x_a)
    ax_a.set_xticklabels(CONFIG_ORDER, rotation=45, ha="right", fontsize=5)
    ax_a.set_ylim(0.97, 1.00)
    ax_a.set_ylabel("Mean kNN purity\n(k=30, L0.5)", fontsize=6)
    ax_a.tick_params(axis="y", labelsize=5, length=2, pad=1)
    ax_a.tick_params(axis="x", length=2, pad=1)
    for sp in ax_a.spines.values():
        sp.set_linewidth(0.4)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    # ---- Panel b: per-celltype mean for winner ----
    bar_colors_b = [comp_palette.get(L0P5_COMPARTMENT[lbl], neutral_color)
                    for lbl in L0P5_ORDER]
    x_b = np.arange(len(L0P5_ORDER))
    ax_b.bar(x_b, panel_b.values, color=bar_colors_b,
             edgecolor="none", width=0.7)

    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels(L0P5_ORDER, rotation=45, ha="right", fontsize=5)
    ax_b.set_ylim(0.93, 1.00)
    ax_b.set_ylabel("Mean kNN purity\n(joint_nuc_n100)", fontsize=6)
    ax_b.tick_params(axis="y", labelsize=5, length=2, pad=1)
    ax_b.tick_params(axis="x", length=2, pad=1)
    for sp in ax_b.spines.values():
        sp.set_linewidth(0.4)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    # Compartment legend (top-right of panel b)
    import matplotlib.patches as mpatches
    leg_handles = [
        mpatches.Patch(color=comp_palette[c], label=c)
        for c in ("Epithelial", "Stromal", "Immune")
    ]
    ax_b.legend(handles=leg_handles, loc="lower right",
                fontsize=5, frameon=False, handlelength=0.8,
                handletextpad=0.3, ncol=3, columnspacing=0.8,
                bbox_to_anchor=(1.0, -0.45))

    # ---- Save ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / "supp7_5_knn_purity_composite.pdf"
    png_path = args.out_dir / "supp7_5_knn_purity_composite.png"
    data_path = args.out_dir / "supp7_5_knn_purity_data.csv"

    raster_dpi = int(aes.get("rendering", {}).get("dpi", 1200))
    fig.savefig(pdf_path, dpi=raster_dpi)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    log.info("wrote %s", pdf_path)
    log.info("wrote %s", png_path)

    # Companion data CSV (long format)
    panel_a_rows = pd.DataFrame({
        "panel": "a", "config": CONFIG_ORDER,
        "cell_type": "all", "mean_purity": panel_a.values,
    })
    panel_b_rows = pd.DataFrame({
        "panel": "b", "config": WINNER,
        "cell_type": L0P5_ORDER, "mean_purity": panel_b.values,
    })
    pd.concat([panel_a_rows, panel_b_rows], ignore_index=True).to_csv(
        data_path, index=False, float_format="%.6f")
    log.info("wrote %s", data_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
