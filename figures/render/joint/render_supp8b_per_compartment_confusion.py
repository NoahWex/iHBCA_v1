"""render_supp8b_per_compartment_confusion.py — Fig 3 Supp 8b: per-compartment
FLEX-label × joint-label confusion matrix.

For each compartment {Epithelial, Immune, Stromal}, computes the column-
normalized crosstab of FLEX-side L2S labels (from
publication/analysis/annotation/flex/outputs/flex_l2s_labels.csv) against the
joint cascade L1.5 labels (from publication/analysis/annotation/xenium/
outputs/compartment/{C}/cell_annotations.csv) — restricted to FLEX cells in
that compartment. Diagonal mass = the two cascades agree; off-diagonal mass
= where they disagree.

This is the per-compartment counterpart to Supp 7-4's atlas-wide L2S × L1.5
confusion matrix, scoped to per-compartment annotation cascades.

Pattern source:
- publication/figures/render/joint/render_supp7_4_l2s_l1p5_confusion.py
  (full-script adaptation; replace single-axis joint_l1p5.csv read with
  per-compartment cell_annotations.csv loop, and filter L2S by compartment).

No on-plot title, panel letter, or method caption (Illustrator handles those).

Outputs (per compartment):
    supp8b_confusion_{epithelial,immune,stromal}.pdf
    supp8b_confusion_{epithelial,immune,stromal}_data.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp8b_confusion")


COMPARTMENT_FULL = {"epi": "Epithelial", "imm": "Immune", "str": "Stromal"}
COMPARTMENT_TOKEN = {
    "epi": "compartment_epi_cell_annotations",
    "imm": "compartment_imm_cell_annotations",
    "str": "compartment_str_cell_annotations",
}

# Per-compartment artifact reclassifications applied at render time. BMYO-NC
# is reclassified because the cluster is sample-specific (single patient,
# single Xenium sample) and does not represent a generalizable biological
# state. Joint-side artifact labels (canonical + override) form the right-edge
# column block in the confusion matrix; their FLEX-L2S contribution profile
# is the artifact-state evidence (smear vs canonical column diagonals).
ARTIFACT_OVERRIDES = {
    "epi": {"BMYO-NC"},
    "imm": set(),
    "str": set(),
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--compartment", required=True, choices=["epi", "imm", "str"])
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample to 5000 FLEX cells before crosstab.")
    return p.parse_args()


def diagonal_sort(props: pd.DataFrame) -> pd.DataFrame:
    """Sort rows by argmax-col then -row_sum.

    Pattern: render_supp7_4_l2s_l1p5_confusion.py:82-103.
    """
    col_index = {c: i for i, c in enumerate(props.columns)}
    argmax_col = props.idxmax(axis=1)
    sort_key = argmax_col.map(col_index).astype(int)
    row_sums = props.sum(axis=1)
    order = (
        pd.DataFrame({"col_idx": sort_key, "row_sum": -row_sums})
        .sort_values(by=["col_idx", "row_sum"]).index
    )
    return props.loc[order]


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    config_dir = str(project_root / "publication" / "config")

    sys.path.insert(0, config_dir)
    from load_aesthetics import (  # type: ignore
        load_aesthetics, get_matplotlib_theme,
    )
    from load_paths import load_paths, resolve_path  # type: ignore

    paths = load_paths(config_dir=config_dir)
    aes = load_aesthetics(config_dir=config_dir)
    plt.rcParams.update(get_matplotlib_theme(aes))

    compartment_short = args.compartment
    compartment_full = COMPARTMENT_FULL[compartment_short]

    flex_path = Path(resolve_path("flex_l2s_labels", paths))
    ann_path = Path(resolve_path(COMPARTMENT_TOKEN[compartment_short], paths))
    log.info("FLEX L2S labels: %s", flex_path)
    log.info("Joint annotations: %s", ann_path)

    flex = pd.read_csv(flex_path, usecols=["cell_id", "compartment", "l2s_label"])
    flex_comp = flex[flex["compartment"] == compartment_short].copy()
    log.info("FLEX cells in %s: %d (L2S unique: %d)",
             compartment_full, len(flex_comp), flex_comp["l2s_label"].nunique())

    ann = pd.read_csv(ann_path, usecols=["cell_id", "label", "is_artifact"])
    ann["is_artifact"] = ann["is_artifact"].astype(bool)

    # Apply ARTIFACT_OVERRIDES so reclassified labels group with cascade
    # artifacts in the right-edge confusion column block.
    overrides = ARTIFACT_OVERRIDES.get(compartment_short, set())
    if overrides:
        n_pre = int(ann["is_artifact"].sum())
        ann.loc[ann["label"].isin(overrides), "is_artifact"] = True
        n_post = int(ann["is_artifact"].sum())
        log.info("Override applied: %s → +%d cells reclassified as artifact",
                 sorted(overrides), n_post - n_pre)

    log.info("Joint annotations: %d cells (artifacts=%d, %d unique labels)",
             len(ann), int(ann["is_artifact"].sum()), ann["label"].nunique())

    if args.test:
        n_sub = min(5000, len(flex_comp))
        flex_comp = flex_comp.sample(n=n_sub, random_state=42).reset_index(drop=True)
        log.info("--test: subsampled FLEX to %d cells", len(flex_comp))

    merged = flex_comp.merge(ann, on="cell_id", how="inner")
    log.info("Inner-joined: %d FLEX cells with joint labels (retention %.1f%%)",
             len(merged), 100 * len(merged) / max(len(flex_comp), 1))

    if len(merged) == 0:
        sys.exit(f"FATAL: no FLEX cells in {compartment_full} matched joint annotations")

    # Joint-side artifact cells are retained as columns. Their lack of clean
    # alignment to any FLEX L2S identity is the artifact-state evidence for
    # this panel; diagonal sort plus a separator after the canonical block
    # surfaces the smear pattern visually.
    artifact_joint_labels = sorted(set(merged.loc[merged["is_artifact"], "label"]))
    canonical_joint_labels = sorted(set(merged.loc[~merged["is_artifact"], "label"]))
    log.info("Joint labels: %d canonical + %d artifact",
             len(canonical_joint_labels), len(artifact_joint_labels))

    # Crosstab: rows = FLEX L2S labels, cols = joint L1.5 labels (all).
    ct = pd.crosstab(merged["l2s_label"], merged["label"])
    log.info("Crosstab: %d L2S × %d joint labels", *ct.shape)

    # Column-normalize: each joint label column sums to 1.
    # Reads as "for each joint label, which FLEX L2S labels contributed?"
    # For artifact columns: a clean signal would mean they ARE biology
    # (concerning); diffuse / smeared columns are the expected artifact pattern.
    props = ct.div(ct.sum(axis=0), axis=1).fillna(0.0)

    # Diagonal sort within canonical columns first (preserves the bio
    # cascade); artifact columns appended at the right edge in cascade order.
    canonical_cols_present = [c for c in canonical_joint_labels if c in props.columns]
    artifact_cols_present = [c for c in artifact_joint_labels if c in props.columns]
    canonical_props = props[canonical_cols_present] if canonical_cols_present \
                      else pd.DataFrame(index=props.index)
    if not canonical_props.empty:
        canonical_props = diagonal_sort(canonical_props.T).T
        canonical_props = diagonal_sort(canonical_props)
        # Apply the canonical row order to the full props (carries to artifacts).
        props = props.loc[canonical_props.index, canonical_cols_present + artifact_cols_present]
    else:
        props = props[canonical_cols_present + artifact_cols_present]
    log.info("After sort: %d × %d (canonical=%d, artifact=%d)",
             props.shape[0], props.shape[1],
             len(canonical_cols_present), len(artifact_cols_present))

    # ---- Render ----
    n_rows, n_cols = props.shape
    # Per-compartment confusion matrices vary substantially in column count
    # (Epi ~8 cols, Imm ~17 cols, Str ~8 cols). Width scales with n_cols so
    # the column tick labels (joint L1.5 names) don't crowd. Compose stage
    # scales the panel down to fit the 55mm grid cell, but adequate
    # source-panel width preserves text legibility post-scale.
    width_in = max(70.0, 30.0 + n_cols * 5.5) / 25.4
    height_in = max(50.0, 35.0 + n_rows * 4.0) / 25.4

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    sns.heatmap(
        props,
        cmap="Blues",
        vmin=0, vmax=1,
        ax=ax,
        linewidths=0.3,
        linecolor="white",
        cbar_kws={"shrink": 0.4, "label": "P(L2S | joint label)"},
        annot=(n_rows <= 20 and n_cols <= 20),
        fmt=".2f",
        annot_kws={"size": 4},
    )
    # Heavy vertical separator between canonical joint labels and artifact
    # joint labels — surfaces the "artifact columns smear across L2S labels"
    # signal vs the canonical column diagonals.
    if artifact_cols_present and canonical_cols_present:
        ax.axvline(len(canonical_cols_present), color="black", linewidth=0.8)

    ax.set_xlabel("Joint L1.5 label")
    ax.set_ylabel("FLEX L2S label")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor", fontsize=5)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=5)

    # Tag artifact column tick labels with gray to differentiate from canonical
    # at a glance.
    if artifact_cols_present:
        artifact_set = set(artifact_cols_present)
        for tick in ax.get_xticklabels():
            if tick.get_text() in artifact_set:
                tick.set_color("#666666")
                tick.set_fontstyle("italic")

    # ---- Save ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / f"supp8b_confusion_{compartment_full.lower()}.pdf"
    csv_path = args.out_dir / f"supp8b_confusion_{compartment_full.lower()}_data.csv"

    fig.tight_layout()
    fig.savefig(pdf_path, bbox_inches="tight",
                dpi=int(aes.get("rendering", {}).get("dpi", 600)))
    plt.close(fig)
    log.info("Wrote %s", pdf_path)

    props.to_csv(csv_path)
    log.info("Wrote %s (%d × %d)", csv_path, *props.shape)
    return 0


if __name__ == "__main__":
    sys.exit(main())
