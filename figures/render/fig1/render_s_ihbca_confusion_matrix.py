#!/usr/bin/env python3
"""
Fig 1 supp S3 annotation validation: per-study native author cell-type labels
vs the harmonized iHBCA cell-type vocabulary at L1 (11 canonical types) and L2
(42 fine types), rendered as a stacked unified composite with shared author
x-axis. Both row axes (L1, L2) are compartment-blocked (Epithelial → Stromal →
Immune); author x-axis is ordered by L2 cascade within each compartment block.
Bracket lines mark compartment transitions on all three axes.

Substrate pipeline (3 inputs):

  1. obs_names_to_numeric_id.csv (from V1_Annotation/dev/l20_annotation/outputs/)
     numeric_id <-> obs_name (cell_id with sample prefix, e.g. "PM-D_AAACCTGAGAATCTCC").
     Canonical bridge between numeric_id (used in per_study_native_labels) and
     cell_id (used in panel_substrate).

  2. per_study_native_labels.csv (sibling of above)
     Long-form: (numeric_id, native_label, study). One row per (cell, contributing-
     study-author-label). ~2.08M rows (only cells with native labels). Eliminates
     the need to pivot wide native_<study> columns.

  3. publication/figures/data/fig1/panel_substrate.csv
     Per-cell L1 (l1_annotation) + L2 (label) under the canonical iHBCA vocabulary.
     Indexed by cell_id. l1_annotation = Myeloid/Basal-myoepithelial/Luminal hormone
     sensing/etc. — the 11-type canonical L1. label = Macro_C1Q/LASP-major/etc. —
     the L2 fine-type vocabulary.

The L2 -> L1 mapping is derived directly from panel_substrate (groupby label,
modal l1_annotation per L2) — no join required for that mapping itself. Validated
at 97-99%+ purity across all L2 types in the 2026-05-16 audit.

Cross-tab + cascade ordering + GridSpec layout unchanged from prior version.

CLI (same shape as superseded publication script + new substrate args):
  --obs-id-map:  obs_names_to_numeric_id.csv
  --native-long: per_study_native_labels.csv  (numeric_id, native_label, study)
  --panel-substrate: panel_substrate.csv  (cell_id, l1_annotation, label, ...)
  --outdir:      destination directory for PDF + PNG + data CSV
  --level:       L1 | L2 | unified

Substrate locations (HPC, fetched via hpc file cat to /tmp/):
  /share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Annotation/dev/l20_annotation/outputs/obs_names_to_numeric_id.csv
  /share/crsp/lab/dalawson/nwechter/iHBCA_V1/Analysis/stages/V1_Annotation/dev/l20_annotation/outputs/per_study_native_labels.csv
  /share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/figures/data/fig1/panel_substrate.csv  (already in publication/)
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import pandas as pd
import seaborn as sns


# L0 compartment for canonical L1 labels (color-bar mapping).
L1_TO_L0 = {
    "Basal-myoepithelial": "Epithelial",
    "Luminal hormone sensing": "Epithelial",
    "Luminal adaptive secretory precursor": "Epithelial",
    "Lactocyte": "Epithelial",
    "Fibroblast": "Stromal",
    "Perivascular": "Stromal",
    "Vascular endothelial": "Stromal",
    "Lymphatic endothelial": "Stromal",
    "Myeloid": "Immune",
    "T-lymphocyte": "Immune",
    "B-lymphocyte": "Immune",
}

L0_COLORS = {
    "Epithelial": "#7570B3",
    "Stromal": "#D95F02",
    "Immune": "#1B9E77",
    "Other": "#B09C85",
}

STUDY_COLORS = {
    "gray": "#D62728",
    "kumar": "#1F77B4",
    "murrow": "#2CA02C",
    "nee": "#9467BD",
    "pal": "#E377C2",
    "reed": "#8C564B",
    "twigger": "#FF7F0E",
}


def load_substrate(obs_id_map_path, native_long_path, panel_substrate_path):
    """Load + join the three substrates into one long-form table.

    Returns DataFrame with columns: cell_id, numeric_id, study, native_label,
    l1_annotation, label (L2), is_artifact.
    """
    print(f"Loading obs_names_to_numeric_id from {obs_id_map_path} ...")
    obs_map = pd.read_csv(obs_id_map_path)
    obs_map = obs_map.rename(columns={"obs_name": "cell_id"})
    print(f"  {len(obs_map):,} cell_id <-> numeric_id mappings")

    print(f"Loading per_study_native_labels from {native_long_path} ...")
    native_long = pd.read_csv(native_long_path)
    print(f"  {len(native_long):,} (cell, study, native_label) rows; "
          f"{native_long['study'].nunique()} studies")

    print(f"Loading panel_substrate from {panel_substrate_path} ...")
    substrate = pd.read_csv(panel_substrate_path,
                             usecols=["cell_id", "l1_annotation", "label", "is_artifact"])
    print(f"  {len(substrate):,} cells in panel_substrate")

    print("Joining native_long + obs_map -> cell_id ...")
    native_with_cell = native_long.merge(obs_map, on="numeric_id", how="inner")
    print(f"  {len(native_with_cell):,} rows after numeric_id -> cell_id bridge")

    print("Joining + panel_substrate -> l1_annotation + label ...")
    joined = native_with_cell.merge(substrate, on="cell_id", how="inner")
    print(f"  {len(joined):,} rows after panel_substrate join")

    # Drop artifact cells (Doublet, Single nuclei) at substrate level.
    n_before = len(joined)
    joined = joined[~joined["is_artifact"].fillna(False)]
    print(f"  {len(joined):,} rows after dropping is_artifact ({n_before - len(joined):,} removed)")

    return joined


def build_l2_to_l1_from_substrate(joined: pd.DataFrame) -> dict:
    """Derive L2 -> L1 from the joined long-form table.

    Modal l1_annotation per L2 label. Purity audit (2026-05-16): all L2 types
    >= 97% modal purity in their L1 family.
    """
    modal = joined.groupby("label")["l1_annotation"].agg(
        lambda s: s.value_counts().index[0]
    )
    return modal.to_dict()


def build_l2_to_l0_from_substrate(joined: pd.DataFrame) -> dict:
    """Derive L2 -> L0 from the joined long-form table.

    Modal compartment per L2 label via L1_TO_L0 chain. Equivalent: modal l1
    -> L0 lookup.
    """
    l2_to_l1 = build_l2_to_l1_from_substrate(joined)
    return {l2: L1_TO_L0.get(l1, "Other") for l2, l1 in l2_to_l1.items()}


def _cascade_orderings(joined: pd.DataFrame, ihbca_col: str):
    """Compute author-study × iHBCA crosstab + cascade orderings.

    ihbca_col: 'l1_annotation' for L1 panel, 'label' for L2 panel.

    Returns: (ct_prop_author_x_ihbca, ihbca_order, author_order).
    """
    df = joined.copy()
    df[ihbca_col] = df[ihbca_col].astype(str)
    df["native_label"] = df["native_label"].astype(str)
    df["study"] = df["study"].astype(str)
    df = df[~df[ihbca_col].isin(["Doublet", "Single nuclei", "nan"])]
    df["author_study_label"] = df["study"] + ": " + df["native_label"]

    ct = pd.crosstab(df["author_study_label"], df[ihbca_col])
    # Row-normalize in (author × ihbca) orientation: each author row sums to 1.
    # After transpose for display (rows=ihbca, cols=author), each displayed author
    # column sums to 1 = P(iHBCA | author): where each author's cells end up in the
    # iHBCA vocabulary.
    ct_prop = ct.div(ct.sum(axis=1), axis=0)

    dominant_per_author = ct_prop.idxmax(axis=1)
    n_correspondences = dominant_per_author.value_counts()
    ihbca_abundance = df[ihbca_col].value_counts()
    col_score = pd.DataFrame({
        "n_labels": n_correspondences,
        "n_cells": ihbca_abundance,
    }).fillna(0)
    col_score = col_score.sort_values(["n_labels", "n_cells"], ascending=[False, False])
    ihbca_order = [c for c in col_score.index if c in ct_prop.columns]
    ct_prop = ct_prop[ihbca_order]

    dominant_col_idx = ct_prop.values.argmax(axis=1)
    dominant_prop = ct_prop.values.max(axis=1)
    row_meta = pd.DataFrame({
        "author_study_label": ct_prop.index,
        "dominant_col": dominant_col_idx,
        "dominant_prop": dominant_prop,
    })
    row_meta = row_meta.sort_values(["dominant_col", "dominant_prop"], ascending=[True, False])
    ct_prop = ct_prop.loc[row_meta["author_study_label"]]
    return ct_prop, ihbca_order, list(ct_prop.index)


def make_unified_confusion_matrix(joined: pd.DataFrame,
                                  l2_to_l1: dict,
                                  l2_to_l0: dict,
                                  out_path: Path):
    """Stacked L1 + L2 panel; shared author x-axis locked by L2 cascade.

    Both L1 and L2 row axes are compartment-blocked (Epithelial → Stromal → Immune),
    then cascade-ordered within each compartment. Author x-axis ordered by L2 cascade
    so the L2 panel reads as a clean diagonal at the finest available resolution;
    the L1 panel inherits the same x-axis and reads as a coarser-grain version
    of the same compartment-blocked structure.
    """
    COMPARTMENT_ORDER = {"Epithelial": 0, "Stromal": 1, "Immune": 2, "Other": 3}

    # Compute base cascades; we then re-derive author order to match compartment blocking.
    ct_l2_full, l2_order, _ = _cascade_orderings(joined, "label")
    ct_l1_full, l1_order, _ = _cascade_orderings(joined, "l1_annotation")

    # Author x-axis order: compartment-first (via dominant-L2 compartment),
    # then within-compartment by L2 cascade position, then by dominant proportion.
    l2_cascade_pos = {l2: i for i, l2 in enumerate(l2_order)}
    dom_l2_per_author = ct_l2_full.idxmax(axis=1)
    dom_prop_per_author = ct_l2_full.values.max(axis=1)
    author_meta = pd.DataFrame({
        "author": ct_l2_full.index,
        "dom_l2": dom_l2_per_author.values,
        "dom_prop": dom_prop_per_author,
    })
    author_meta["compartment"] = author_meta["dom_l2"].map(lambda lbl: l2_to_l0.get(lbl, "Other"))
    author_meta["comp_pos"] = author_meta["compartment"].map(
        lambda c: COMPARTMENT_ORDER.get(c, len(COMPARTMENT_ORDER))
    )
    author_meta["dom_l2_pos"] = author_meta["dom_l2"].map(
        lambda l: l2_cascade_pos.get(l, len(l2_cascade_pos))
    )
    author_meta = author_meta.sort_values(
        ["comp_pos", "dom_l2_pos", "dom_prop"],
        ascending=[True, True, False],
    )
    author_order = list(author_meta["author"])

    ct_l1 = ct_l1_full.reindex(author_order)
    ct_l2 = ct_l2_full.reindex(author_order)

    plot_l1 = ct_l1.T  # 11 × 208
    plot_l2 = ct_l2.T  # 42 × 208

    # ----- L1 row ordering: compartment-blocked + within-compartment cascade -----
    l1_dom_col_idx = plot_l1.values.argmax(axis=1)
    l1_dom_prop = plot_l1.values.max(axis=1)
    l1_row_meta = pd.DataFrame({
        "ihbca_l1": plot_l1.index,
        "compartment": [L1_TO_L0.get(lbl, "Other") for lbl in plot_l1.index],
        "dom_col": l1_dom_col_idx,
        "dom_prop": -l1_dom_prop,
    })
    l1_row_meta["compartment_pos"] = l1_row_meta["compartment"].map(
        lambda x: COMPARTMENT_ORDER.get(x, len(COMPARTMENT_ORDER))
    )
    l1_row_meta = l1_row_meta.sort_values(["compartment_pos", "dom_col", "dom_prop"])
    plot_l1 = plot_l1.loc[l1_row_meta["ihbca_l1"]]

    # ----- L2 row ordering: compartment-blocked + within-compartment cascade -----
    l2_dom_col_idx = plot_l2.values.argmax(axis=1)
    l2_dom_prop = plot_l2.values.max(axis=1)
    l2_row_meta = pd.DataFrame({
        "ihbca_l2": plot_l2.index,
        "compartment": [l2_to_l0.get(lbl, "Other") for lbl in plot_l2.index],
        "dom_col": l2_dom_col_idx,
        "dom_prop": -l2_dom_prop,
    })
    l2_row_meta["compartment_pos"] = l2_row_meta["compartment"].map(
        lambda x: COMPARTMENT_ORDER.get(x, len(COMPARTMENT_ORDER))
    )
    l2_row_meta = l2_row_meta.sort_values(["compartment_pos", "dom_col", "dom_prop"])
    plot_l2 = plot_l2.loc[l2_row_meta["ihbca_l2"]]

    n_cols = plot_l1.shape[1]
    n_rows_l1 = plot_l1.shape[0]
    n_rows_l2 = plot_l2.shape[0]

    # Layout
    heatmap_width = max(12, n_cols * 0.075)
    height_l1 = max(2.5, n_rows_l1 * 0.3)
    height_l2 = max(5, n_rows_l2 * 0.3)
    compartment_bar_width = 0.25
    prop_cbar_width = 0.3
    legend_width = 2.2
    fig_width = compartment_bar_width + heatmap_width + prop_cbar_width + legend_width + 1.0
    fig_height = height_l1 + height_l2 + 3.5

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(
        2, 4, figure=fig,
        width_ratios=[compartment_bar_width, heatmap_width, prop_cbar_width, legend_width],
        height_ratios=[height_l1, height_l2],
        wspace=0.02, hspace=0.04,
    )

    ax_compartment_l1 = fig.add_subplot(gs[0, 0])
    ax_l1 = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[0, 3])
    ax_compartment_l2 = fig.add_subplot(gs[1, 0])
    ax_l2 = fig.add_subplot(gs[1, 1])
    ax_prop_full = fig.add_subplot(gs[1, 2])
    for empty in [gs[0, 2], gs[1, 3]]:
        fig.add_subplot(empty).axis("off")
    ax_legend.axis("off")
    ax_prop_full.axis("off")

    # L1 heatmap (top)
    sns.heatmap(plot_l1, ax=ax_l1, cmap="Blues", vmin=0, vmax=1,
                linewidths=0.2, linecolor="white", cbar=False,
                xticklabels=False, yticklabels=False)
    ax_l1.set_xlabel("")
    ax_l1.set_ylabel("")
    ax_l1.tick_params(axis="x", bottom=False)
    ax_l1.tick_params(axis="y", left=False)

    # L1 compartment strip
    row_l0_colors_l1 = [L0_COLORS.get(L1_TO_L0.get(r, "Other"), L0_COLORS["Other"]) for r in plot_l1.index]
    ax_compartment_l1.imshow(
        [[i] for i in range(n_rows_l1)], aspect="auto",
        cmap=mpl.colors.ListedColormap(row_l0_colors_l1),
    )
    ax_compartment_l1.set_ylim(n_rows_l1 - 0.5, -0.5)
    ax_compartment_l1.set_xticks([])
    ax_compartment_l1.set_yticks(range(n_rows_l1))
    ax_compartment_l1.set_yticklabels(list(plot_l1.index), fontsize=6)
    ax_compartment_l1.tick_params(axis="y", left=True, right=False, length=3, pad=2)
    ax_compartment_l1.set_ylabel("iHBCA L1", fontsize=9, fontweight="bold", labelpad=4)

    # L2 heatmap (bottom, with shared author x-tick labels)
    sns.heatmap(plot_l2, ax=ax_l2, cmap="Blues", vmin=0, vmax=1,
                linewidths=0.2, linecolor="white", cbar=False,
                xticklabels=True, yticklabels=False)
    ax_l2.set_xlabel("")
    ax_l2.set_ylabel("")
    ax_l2.tick_params(axis="y", left=False)
    xlabels = ax_l2.get_xticklabels()
    for label in xlabels:
        text = label.get_text()
        study = text.split(": ", 1)[0] if ": " in text else ""
        label.set_color(STUDY_COLORS.get(study, "#333333"))
        label.set_fontsize(4.5)
        label.set_rotation(90)
    ax_l2.set_xticklabels(xlabels)

    # L2 compartment strip
    row_l0_colors_l2 = [L0_COLORS.get(l2_to_l0.get(r, "Other"), L0_COLORS["Other"]) for r in plot_l2.index]
    ax_compartment_l2.imshow(
        [[i] for i in range(n_rows_l2)], aspect="auto",
        cmap=mpl.colors.ListedColormap(row_l0_colors_l2),
    )
    ax_compartment_l2.set_ylim(n_rows_l2 - 0.5, -0.5)
    ax_compartment_l2.set_xticks([])
    ax_compartment_l2.set_yticks(range(n_rows_l2))
    ax_compartment_l2.set_yticklabels(list(plot_l2.index), fontsize=5)
    ax_compartment_l2.tick_params(axis="y", left=True, right=False, length=3, pad=2)
    ax_compartment_l2.set_ylabel("iHBCA L2", fontsize=9, fontweight="bold", labelpad=4)

    # Compartment-block bracket lines on both panels
    # Author x-axis boundaries: when dominant-L2's compartment changes
    dominant_l2_per_author = ct_l2.idxmax(axis=1)
    dominant_compartment_per_author = dominant_l2_per_author.map(
        lambda lbl: l2_to_l0.get(lbl, "Other")
    )
    col_boundaries = []
    prev_c = None
    for i, c in enumerate(dominant_compartment_per_author.values):
        if prev_c is not None and c != prev_c:
            col_boundaries.append(i)
        prev_c = c

    # L1 row boundaries: between compartment blocks
    l1_compartment_seq = [L1_TO_L0.get(lbl, "Other") for lbl in plot_l1.index]
    row_boundaries_l1 = []
    prev = None
    for i, p in enumerate(l1_compartment_seq):
        if prev is not None and p != prev:
            row_boundaries_l1.append(i)
        prev = p

    # L2 row boundaries: between compartment blocks
    l2_compartment_seq = [l2_to_l0.get(lbl, "Other") for lbl in plot_l2.index]
    row_boundaries_l2 = []
    prev = None
    for i, p in enumerate(l2_compartment_seq):
        if prev is not None and p != prev:
            row_boundaries_l2.append(i)
        prev = p

    for col_idx in col_boundaries:
        ax_l1.axvline(x=col_idx, color="black", linewidth=0.9, alpha=0.7)
        ax_l2.axvline(x=col_idx, color="black", linewidth=0.9, alpha=0.7)
    for row_idx in row_boundaries_l1:
        ax_l1.axhline(y=row_idx, color="black", linewidth=0.9, alpha=0.7)
    for row_idx in row_boundaries_l2:
        ax_l2.axhline(y=row_idx, color="black", linewidth=0.9, alpha=0.7)

    # Proportion colorbar
    fig.canvas.draw()
    hm_pos = ax_l2.get_position()
    cbar_height_in = min(2.5, height_l2 * 0.5)
    cbar_bottom = hm_pos.y0 + (hm_pos.height - cbar_height_in / fig_height) / 2
    cbar_left = hm_pos.x1 + 0.005
    cbar_w = 0.15 / fig_width
    cbar_h = cbar_height_in / fig_height
    ax_cbar = fig.add_axes([cbar_left, cbar_bottom, cbar_w, cbar_h])
    sm = plt.cm.ScalarMappable(cmap="Blues", norm=mpl.colors.Normalize(vmin=0, vmax=1))
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax_cbar)
    cb.set_label("Proportion", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    # Legends
    l0_patches = [Patch(facecolor=c, label=l) for l, c in L0_COLORS.items() if l != "Other"]
    leg_l0 = ax_legend.legend(
        handles=l0_patches, title="Compartment",
        loc="upper left", bbox_to_anchor=(0.05, 1.0),
        fontsize=7, title_fontsize=8, frameon=True,
    )
    ax_legend.add_artist(leg_l0)
    study_patches = [Patch(facecolor=c, label=s) for s, c in STUDY_COLORS.items()]
    ax_legend.legend(
        handles=study_patches, title="Study",
        loc="upper left", bbox_to_anchor=(0.05, 0.45),
        fontsize=7, title_fontsize=8, frameon=True,
    )

    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved (unified): {out_path} (L1: {n_rows_l1} x {n_cols}; L2: {n_rows_l2} x {n_cols})")
    return ct_l1, ct_l2


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--obs-id-map", required=True,
                        help="obs_names_to_numeric_id.csv (V1_Annotation/dev/l20_annotation/outputs)")
    parser.add_argument("--native-long", required=True,
                        help="per_study_native_labels.csv (numeric_id, native_label, study)")
    parser.add_argument("--panel-substrate", required=True,
                        help="panel_substrate.csv (cell_id, l1_annotation, label, is_artifact)")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--level", required=True, choices=["L1", "L2", "unified"],
                        help="L1 / L2 single-level, or 'unified' (stacked L1+L2 with shared author x-axis)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    joined = load_substrate(args.obs_id_map, args.native_long, args.panel_substrate)

    if args.level == "unified":
        print("Building L2 -> L1 mapping from substrate (modal l1_annotation per L2) ...")
        l2_to_l1 = build_l2_to_l1_from_substrate(joined)
        l2_to_l0 = {l2: L1_TO_L0.get(l1, "Other") for l2, l1 in l2_to_l1.items()}
        out_pdf = outdir / "s1_4_confusion_matrix.pdf"
        out_png = outdir / "s1_4_confusion_matrix.png"
        ct_l1, ct_l2 = make_unified_confusion_matrix(joined, l2_to_l1, l2_to_l0, out_pdf)
        make_unified_confusion_matrix(joined, l2_to_l1, l2_to_l0, out_png)
        ct_l1.to_csv(outdir / "s1_4_confusion_matrix_L1_data.csv")
        ct_l2.to_csv(outdir / "s1_4_confusion_matrix_L2_data.csv")
        print("\nDone.")
        return

    # Single-level path (L1 or L2) — reuses unified logic for the level requested.
    raise NotImplementedError("Single-level rendering not wired up in this DRAFT; "
                              "only --level unified is implemented. Add L1/L2 standalone "
                              "paths post-validation if needed.")


if __name__ == "__main__":
    main()
