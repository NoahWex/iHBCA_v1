"""
render_supp7_4_l2s_l1p5_confusion.py — Fig 3 Supp 7-4.

L2S x L1.5 cross-cascade confusion matrix for FLEX cells: each FLEX cell
carries an L2S label (from the FLEX-side annotation cascade in
publication/analysis/annotation/flex/) and an L1.5 label (from the joint
FLEX+Xenium cascade in publication/analysis/annotation/xenium/). Diagonal
mass indicates the two cascades agree on cell-type identity; off-diagonal
mass shows where they split or merge.

Inputs:
    paths.flex_l2s_labels      cell_id, l2s_label per FLEX cell
    paths.xenium_joint_l1p5    cell_id, l1p5_label, compartment, platform per
                               joint cell (filter platform == 'flex')

Outputs:
    supp7_4_l2s_l1p5_confusion.pdf
    supp7_4_l2s_l1p5_confusion_data.csv  (row-normalized proportions)

Pattern source:
    iHBCA_V1/Analysis/stages/V1_Abundance/.dev/paper_supps_20260530/
    1S3_confusion_matrices/scripts/render_1S3_confusion_matrices.py
    diagonal_sort + plot_confusion ported with two changes:
        (a) compartment colour strip uses get_palette('compartment') instead of
            the inline L0_COLORS dict
        (b) figure dimensions come from get_dimensions('heatmap_tile') with
            density-rule height scaling

Framework conformance: Option B per CP_supp7_python_framework_gap.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]  # publication/figures/render/joint -> 3 up
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
log = logging.getLogger("supp7_4")


MIN_JOINED_ROWS = 100_000  # CP threshold per handoff §174 schema-verification gate


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--project-root", required=True, type=Path,
                   help="Path to iHBCA_publication root.")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Output directory for PDF + data CSV.")
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K FLEX cells before crosstab.")
    return p.parse_args()


def diagonal_sort(props: pd.DataFrame) -> pd.DataFrame:
    """Sort rows so the dominant column for each row forms a diagonal.

    Ported verbatim from
    iHBCA_V1/.../1S3_confusion_matrices/scripts/render_1S3_confusion_matrices.py
    (diagonal_sort, lines 82-97).

    For each row, find the column with the highest proportion (argmax). Sort
    rows first by that argmax column index, then within each column group by
    descending row sum (puts more populous / more concentrated labels first).
    """
    col_names = list(props.columns)
    col_index = {name: i for i, name in enumerate(col_names)}
    argmax_col = props.idxmax(axis=1)
    sort_key = argmax_col.map(col_index).astype(int)
    row_sums = props.sum(axis=1)
    order = (
        pd.DataFrame({"col_idx": sort_key, "row_sum": -row_sums})
        .sort_values(by=["col_idx", "row_sum"])
        .index
    )
    return props.loc[order]


def main() -> int:
    args = parse_args()
    config_dir = str(args.project_root / "publication" / "config")
    paths = load_paths(config_dir=config_dir)
    aes = load_aesthetics(config_dir=config_dir)

    # ---- Schema verification gate (per handoff §174) ----
    flex_path = resolve_path("flex_l2s_labels", paths)
    joint_path = resolve_path("xenium_joint_l1p5", paths)

    log.info("loading flex L2S labels: %s", flex_path)
    flex = pd.read_csv(flex_path, usecols=["cell_id", "l2s_label"])
    log.info("  flex rows: %d | L2S unique values: %d", len(flex), flex["l2s_label"].nunique())

    log.info("loading joint L1.5 labels: %s", joint_path)
    # l1p5_short is the compact display form (e.g. 'Fb' for 'Fibroblast');
    # using short labels on the column axis keeps the matrix legible at
    # heatmap_tile width.
    joint = pd.read_csv(
        joint_path,
        usecols=["cell_id", "platform", "compartment", "l1p5_short", "is_artifact"],
    )
    log.info(
        "  joint rows: %d | platforms: %s | L1.5 unique values: %d",
        len(joint),
        sorted(joint["platform"].dropna().unique()),
        joint["l1p5_short"].nunique(),
    )

    # Filter joint to FLEX platform; drop artifact cells (not biology-facing)
    joint_flex = joint[(joint["platform"] == "flex") & (~joint["is_artifact"].astype(bool))].copy()
    log.info("  joint after platform=='flex' & !is_artifact filter: %d", len(joint_flex))

    # Subset for --test
    if args.test:
        n_sub = min(50_000, len(joint_flex))
        joint_flex = joint_flex.sample(n=n_sub, random_state=42)
        log.info("  --test subsample: %d", len(joint_flex))

    merged = flex.merge(joint_flex, on="cell_id", how="inner")
    log.info(
        "  inner-join on cell_id: %d rows (FLEX %d -> joined %d, retention %.1f%%)",
        len(merged), len(flex), len(merged), 100 * len(merged) / max(len(flex), 1),
    )

    if len(merged) < MIN_JOINED_ROWS and not args.test:
        log.error(
            "joined rows %d below threshold %d. STOP — write join diagnostic.",
            len(merged), MIN_JOINED_ROWS,
        )
        diag_lines = [
            f"Joined rows: {len(merged)} (threshold {MIN_JOINED_ROWS})",
            f"flex_l2s_labels rows: {len(flex)}",
            f"joint_l1p5 rows after platform=='flex' & !is_artifact: {len(joint_flex)}",
            "",
            "Sample flex cell_ids: " + ", ".join(flex["cell_id"].head(3).astype(str).tolist()),
            "Sample joint cell_ids: " + ", ".join(
                joint_flex["cell_id"].head(3).astype(str).tolist()
            ),
        ]
        diag_path = args.project_root / "coordination" / "reports" / "CP_supp7_4_join_diagnostic.md"
        diag_path.parent.mkdir(parents=True, exist_ok=True)
        diag_path.write_text(
            "# CP — Supp 7-4 join diagnostic\n\n"
            "Schema-verification gate tripped during pre-render schema check.\n\n"
            "## Diagnostic\n\n```\n" + "\n".join(diag_lines) + "\n```\n"
        )
        raise SystemExit(2)

    # ---- Build crosstab + column-normalize ----
    # Column-normalized: each L1.5 column sums to 1; cell value = P(L2S | L1.5).
    # Reads as "for each L1.5 label, which L2S labels contributed to it?"
    ct = pd.crosstab(merged["l2s_label"], merged["l1p5_short"])
    log.info("crosstab shape: %d L2S labels x %d L1.5 labels", ct.shape[0], ct.shape[1])
    props = ct.div(ct.sum(axis=0), axis=1)
    # Two-axis diagonal sort: sort columns first (by their dominant row), then
    # sort rows (by their dominant column in the new column order). diagonal_sort
    # operates on rows by default; apply via transpose to sort columns.
    props = diagonal_sort(props.T).T
    props = diagonal_sort(props)

    # ---- L1.5 -> compartment map (from joint_l1p5 dominant compartment per L1.5) ----
    l1p5_to_compartment = (
        joint_flex.groupby("l1p5_short")["compartment"]
        .agg(lambda s: s.value_counts().idxmax())
        .to_dict()
    )

    # ---- Theme + dimensions ----
    plt.rcParams.update(get_matplotlib_theme(aes))
    heat_dims = get_dimensions("heatmap_tile", aes)
    base_w = float(heat_dims["width"])  # 7.2
    base_h = float(heat_dims["height"])  # 5.0
    n_rows = props.shape[0]
    extra_per_row = 0.12
    scaled_h = base_h + max(0, (n_rows - 25) * extra_per_row)
    max_h_in = float(aes["dimensions"]["max_height_mm"]) / 25.4
    fig_h = min(scaled_h, max_h_in)
    fig_w = base_w
    log.info("figure: %.2fin x %.2fin (rows=%d)", fig_w, fig_h, n_rows)

    # ---- Render ----
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(
        props,
        cmap="Blues",
        vmin=0,
        vmax=1,
        ax=ax,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"shrink": 0.4, "label": "P(L2S | L1.5)"},
        annot=n_rows <= 40,
        fmt=".2f",
        annot_kws={"size": 4},
    )
    ax.set_xlabel("Joint L1.5 label")
    ax.set_ylabel("FLEX L2S label")
    # tick_params doesn't expose ha; set rotation+alignment via setp on
    # the rotated tick label artists so the right edge anchors to the column.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    plt.setp(ax.get_yticklabels(), rotation=0)

    # Compartment colour strip above x-axis labels (per 1.S3 pattern, but with
    # framework palette instead of inline L0_COLORS)
    comp_palette = get_palette("compartment", aes)
    fallback_color = "#CCCCCC"
    for i, col_name in enumerate(props.columns):
        comp = l1p5_to_compartment.get(col_name, "Unknown")
        face = comp_palette.get(comp, fallback_color)
        ax.add_patch(
            Rectangle((i, -0.5), 1, 0.5, color=face,
                      transform=ax.transData, clip_on=False)
        )

    # ---- Save ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / "supp7_4_l2s_l1p5_confusion.pdf"
    csv_path = args.out_dir / "supp7_4_l2s_l1p5_confusion_data.csv"

    fig.tight_layout()
    fig.savefig(pdf_path, bbox_inches="tight", dpi=600)
    plt.close(fig)
    log.info("wrote %s", pdf_path)

    props.to_csv(csv_path)
    log.info("wrote %s (%d rows x %d cols)", csv_path, props.shape[0], props.shape[1])

    return 0


if __name__ == "__main__":
    sys.exit(main())
