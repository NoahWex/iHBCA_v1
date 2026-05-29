"""
Render s1.3 — per-cell-type study contribution stacked bar (iHBCA L2 × 7 studies).

Substrate (resolved via --substrate-dir):
  s4_study_contribution_long.csv  (278 rows: L2_label, study, pct_within_L2, ...)
  s4_l2_order.csv                 (42 L2 with sort_index, L0_compartment, L1_lineage)

Output (written under --out-dir):
  s1_3_study_contribution.pdf
  s1_3_study_contribution.png
  s1_3_study_contribution.panel.yaml

Render spec: Nature double-column (183 mm); no on-plot title / caption /
panel letter (added in Illustrator at composition).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.transforms import blended_transform_factory


# ---- Constants -------------------------------------------------------------

PANEL_ID = "s_ihbca_study_contribution_stackbar"
OUTPUT_SLUG = "s1_3_study_contribution"
COMPARTMENT_LABEL = {"epi": "Epithelial", "str": "Stromal", "imm": "Immune"}


# ---- CLI -------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--substrate-dir", required=True, type=Path,
                   help="Directory containing s4_study_contribution_long.csv and s4_l2_order.csv.")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Destination directory for PDF + PNG + panel.yaml.")
    p.add_argument("--project-root", required=True, type=Path,
                   help="iHBCA_publication root (for publication/config/load_aesthetics.py).")
    return p.parse_args()


# ---- Main ------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    sys.path.insert(0, str(args.project_root / "publication" / "config"))
    from load_aesthetics import (  # noqa: E402
        get_dimensions,
        get_matplotlib_theme,
        get_palette,
        load_aesthetics,
    )

    config = load_aesthetics(config_dir=str(args.project_root / "publication" / "config"))
    plt.rcParams.update(get_matplotlib_theme(config))

    study_palette = get_palette("study", config)        # dict: study -> hex
    study_order = list(study_palette.keys())            # canonical insertion order

    args.out_dir.mkdir(parents=True, exist_ok=True)
    long_csv = args.substrate_dir / "s4_study_contribution_long.csv"
    order_csv = args.substrate_dir / "s4_l2_order.csv"
    pdf_path = args.out_dir / f"{OUTPUT_SLUG}.pdf"
    png_path = args.out_dir / f"{OUTPUT_SLUG}.png"
    yaml_path = args.out_dir / f"{OUTPUT_SLUG}.panel.yaml"

    # --- Load and prepare data ---------------------------------------------
    long_df = pd.read_csv(long_csv)
    order_df = pd.read_csv(order_csv).sort_values("sort_index").reset_index(drop=True)

    # Pivot to (L2 × study) percentage matrix, in canonical order
    pct = (
        long_df.pivot_table(
            index="L2_label", columns="study", values="pct_within_L2", aggfunc="sum"
        )
        .reindex(order_df["L2_label"])
        .reindex(columns=study_order)
        .fillna(0.0)
    )

    n_l2 = len(order_df)

    # --- Figure geometry ----------------------------------------------------
    dims_dbl = get_dimensions("stacked_bar", config)
    fig_w = dims_dbl["width"]   # 7.2 in (183 mm)
    fig_h = 4.2                 # taller than default 3.0 to fit 45deg L2 labels +
                                # compartment brackets + horizontal legend

    fig = plt.figure(figsize=(fig_w, fig_h))
    # Leave headroom on top for compartment brackets and bottom for rotated labels
    ax = fig.add_axes((0.07, 0.32, 0.78, 0.55))

    # --- Stacked bars -------------------------------------------------------
    x = np.arange(n_l2)
    bottom = np.zeros(n_l2)
    bar_w = 0.86
    for study in study_order:
        heights = pct[study].to_numpy()
        ax.bar(
            x,
            heights,
            bottom=bottom,
            width=bar_w,
            color=study_palette[study],
            label=study,
            linewidth=0,
        )
        bottom += heights

    # --- Y axis -------------------------------------------------------------
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("% within L2")
    # Light horizontal gridlines at each 25%
    ax.yaxis.grid(True, color="#CCCCCC", linewidth=0.4, linestyle="-", zorder=0)
    ax.set_axisbelow(True)

    # --- X axis -------------------------------------------------------------
    ax.set_xlim(-0.6, n_l2 - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels(
        order_df["L2_label"].tolist(),
        rotation=45,
        ha="right",
        rotation_mode="anchor",
        fontsize=5.5,
    )
    ax.tick_params(axis="x", length=2, pad=1)
    ax.tick_params(axis="y", length=2)

    # Spines: keep left + bottom, drop top + right
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)

    # --- Thin separators between L1 lineage blocks --------------------------
    l1_seq = order_df["L1_lineage"].tolist()
    l1_breaks = [i for i in range(1, n_l2) if l1_seq[i] != l1_seq[i - 1]]
    for b in l1_breaks:
        ax.axvline(
            x=b - 0.5,
            color="#444444",
            linewidth=0.35,
            zorder=2,
            ymin=0,
            ymax=1,
        )

    # --- Thick brackets above x-axis marking L0 compartment transitions -----
    # Identify contiguous spans of each compartment in canonical order
    l0_seq = order_df["L0_compartment"].tolist()
    spans: list[tuple[str, int, int]] = []
    start = 0
    for i in range(1, n_l2 + 1):
        if i == n_l2 or l0_seq[i] != l0_seq[start]:
            spans.append((l0_seq[start], start, i - 1))
            start = i

    # Vertical L0-boundary cues on the plot itself (slightly thicker than L1)
    for _, _, end in spans[:-1]:
        ax.axvline(
            x=end + 0.5,
            color="#000000",
            linewidth=0.9,
            zorder=3,
            ymin=0,
            ymax=1,
        )

    # Draw brackets in axes-data coords just above the bars (y in axes-fraction
    # is easier; use a twin transform via ax.transData for x and ax.transAxes for y).
    trans = blended_transform_factory(ax.transData, ax.transAxes)
    bracket_y = 1.05      # just above top spine
    bracket_h = 0.04      # short vertical tick at each end
    text_y = 1.12

    for comp, i0, i1 in spans:
        x0 = i0 - 0.5 + 0.05
        x1 = i1 + 0.5 - 0.05
        # Horizontal bar
        ax.plot(
            [x0, x1],
            [bracket_y, bracket_y],
            color="black",
            linewidth=1.4,
            transform=trans,
            clip_on=False,
            solid_capstyle="butt",
        )
        # End ticks
        for xe in (x0, x1):
            ax.plot(
                [xe, xe],
                [bracket_y - bracket_h, bracket_y],
                color="black",
                linewidth=1.4,
                transform=trans,
                clip_on=False,
                solid_capstyle="butt",
            )
        ax.text(
            (x0 + x1) / 2,
            text_y,
            COMPARTMENT_LABEL[comp],
            ha="center",
            va="bottom",
            fontsize=7,
            fontweight="bold",
            transform=trans,
            clip_on=False,
        )

    # --- Legend (compact, right side) ---------------------------------------
    # Reverse handles so top-of-stack (last drawn) appears at top of legend
    handles, labels = ax.get_legend_handles_labels()
    leg = ax.legend(
        handles[::-1],
        labels[::-1],
        title="Study",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        handlelength=1.2,
        handleheight=1.0,
        labelspacing=0.4,
        borderaxespad=0,
        fontsize=6,
        title_fontsize=7,
    )
    leg.get_title().set_fontweight("bold")

    # --- Save ---------------------------------------------------------------
    # Vector PDF (axes/text editable) — Nature primary deliverable
    fig.savefig(pdf_path, format="pdf", dpi=600)
    # 300 DPI PNG preview
    fig.savefig(png_path, format="png", dpi=300)
    plt.close(fig)

    # --- Per-panel YAML stub ------------------------------------------------
    panel_meta = {
        "panel_id": PANEL_ID,
        "output_slug": OUTPUT_SLUG,
        "description": (
            "Per-cell-type contribution by contributing study, stacked-to-100% per L2 "
            "(42 L2 × 7 studies). Backs Fig 1 §2 prose study-imbalance claims."
        ),
        "source_data": [
            "publication/figures/data/fig1/study_contribution/s4_study_contribution_long.csv",
            "publication/figures/data/fig1/study_contribution/s4_l2_order.csv",
        ],
        "tags": [
            "supplemental",
            "fig1",
            "annotation_method",
            "study_contribution",
        ],
        "prose_token": f"[SUPP:{PANEL_ID}]",
        "nature_spec": "single_panel, double_column, vector_pdf, 600_dpi",
        "palette": "study (publication/config/aesthetics.yaml, 7-color Paul Tol bright)",
    }
    with open(yaml_path, "w", encoding="utf-8") as fh:
        yaml.safe_dump(panel_meta, fh, sort_keys=False)

    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")
    print(f"Wrote {yaml_path}")


if __name__ == "__main__":
    main()
