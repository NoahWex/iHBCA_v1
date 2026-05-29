"""
render_supp2_flex_contamination_umap.py — Fig 3 Supp 2, Panel 2a.

Cohort step-08 UMAP colored by `cluster_annotation` (consolidated cell-type
labels at leiden_scvi_1.0 from step 15 MiloContamination decision
application). Contamination-flagged annotations (Keratinocyte, Eccrine
Gland, Melanocyte) get saturated callout colors; mammary-canon annotations
use the compartment palette.

Reviewer question (Supp 2 narrative): "Where in the manifold do the
populations we removed live?"

Inputs (consumes-from-dev — bundle into Supp 2 promotion package):
    Spatial_HBCA_preprocessing/outputs/preprocessing/08_ScviIntegration/
        embeddings/umap_2d.csv
    Spatial_HBCA_preprocessing/outputs/preprocessing/15_MiloContamination/
        metadata/cell_compartments_step15.csv
        (cell_id, leiden_scvi_1.0, cluster_annotation, compartment)

Outputs:
    supp2_flex_contamination_umap.pdf
    supp2_flex_contamination_umap_data.csv

Pattern source:
    publication/figures/render/flex/render_supp1p5_flex_celltype_umap.py
    (Option B framework boilerplate, qualitative palette, rasterized scatter)

Reference for the cluster_annotation vocabulary:
    Spatial_HBCA/project/01_Preprocessing/scripts/step_15_milo_contamination/
    4_decision_application.Rmd (lines 545–897 — defines cluster_annotation
    label set and DoHeatmap usage; this script reuses the label vocabulary
    only, not the per-cell heatmap pattern).

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
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
    get_dimensions,
    get_matplotlib_theme,
    load_aesthetics,
)
from load_paths import load_paths  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp2_2a")


SHUFFLE_SEED = 42

# Contamination-flagged cluster_annotation labels per step 15 decisions
# (cluster_annotations.yaml + 4_decision_application.Rmd consolidation)
CONTAMINATION_LABELS = {"Keratinocyte", "Eccrine_Gland", "Melanocyte"}

# Annotation order: contamination first (saturated), then by compartment
ANNOTATION_ORDER = [
    # Contamination
    "Keratinocyte", "Eccrine_Gland", "Melanocyte",
    # Epithelial
    "Basal", "Luminal_HR", "Luminal_Secretory", "Contractile_Myoepithelial",
    "PS_Epithelial",
    # Stromal
    "Fibroblast", "Pericyte", "VSMC", "Adipocyte",
    "Vascular_Endothelium", "Lymphatic_Endothelium", "Schwann_Cell",
    # Immune
    "T_Cell", "B_Cell", "Plasma_cell", "Myeloid", "Mast", "Neutrophil",
    "pDC", "Langerhans_Cell",
    "RBC",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K cells for sanity render.")
    return p.parse_args()


def get_annotation_palette() -> dict:
    """Per-annotation color, contamination saturated + compartment ramps elsewhere.

    No `cluster_annotation` palette token in aesthetics.yaml — flagged for
    coordinator. Contamination set uses saturated callout colors; remaining
    annotations use sequential colormap blocks per compartment family.
    """
    contam = {
        "Keratinocyte":  "#D6604D",
        "Eccrine_Gland": "#E08214",
        "Melanocyte":    "#9B6FAB",
    }
    epi = plt.get_cmap("Blues")(np.linspace(0.45, 0.85, 5))
    str_ = plt.get_cmap("Greens")(np.linspace(0.35, 0.90, 7))
    imm = plt.get_cmap("Oranges")(np.linspace(0.40, 0.92, 8))
    pal = dict(contam)
    epi_labels = ["Basal", "Luminal_HR", "Luminal_Secretory",
                  "Contractile_Myoepithelial", "PS_Epithelial"]
    for lbl, c in zip(epi_labels, epi):
        pal[lbl] = tuple(c)
    str_labels = ["Fibroblast", "Pericyte", "VSMC", "Adipocyte",
                  "Vascular_Endothelium", "Lymphatic_Endothelium", "Schwann_Cell"]
    for lbl, c in zip(str_labels, str_):
        pal[lbl] = tuple(c)
    imm_labels = ["T_Cell", "B_Cell", "Plasma_cell", "Myeloid", "Mast",
                  "Neutrophil", "pDC", "Langerhans_Cell"]
    for lbl, c in zip(imm_labels, imm):
        pal[lbl] = tuple(c)
    pal["RBC"] = "#888888"
    return pal


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    root = Path(paths["_root"])
    umap_path = (root / "Spatial_HBCA_preprocessing" / "outputs" / "preprocessing"
                 / "08_ScviIntegration" / "embeddings" / "umap_2d.csv")
    cc_path = (root / "Spatial_HBCA_preprocessing" / "outputs" / "preprocessing"
               / "15_MiloContamination" / "metadata" / "cell_compartments_step15.csv")

    log.info("Reading UMAP: %s", umap_path)
    umap = pd.read_csv(umap_path, index_col=0)
    umap.index.name = "cell_id"
    log.info("UMAP: %d cells", len(umap))

    log.info("Reading cell_compartments_step15: %s", cc_path)
    cc = pd.read_csv(cc_path, usecols=["cell_id", "cluster_annotation",
                                       "compartment"])
    log.info("cell_compartments: %d cells", len(cc))

    df = umap.join(cc.set_index("cell_id"), how="inner").reset_index()
    df = df.dropna(subset=["cluster_annotation"])
    log.info("Joined: %d cells; %d cluster_annotation values",
             len(df), df["cluster_annotation"].nunique())

    if args.test:
        df = df.sample(n=50_000, random_state=SHUFFLE_SEED).reset_index(drop=True)
        log.info("Subsampled to %d cells (test)", len(df))

    rng = np.random.default_rng(SHUFFLE_SEED)
    perm = rng.permutation(len(df))
    df = df.iloc[perm].reset_index(drop=True)

    pal = get_annotation_palette()
    from matplotlib.colors import to_rgba
    color_arr = np.array([to_rgba(pal.get(lbl, "#B3B3B3"))
                          for lbl in df["cluster_annotation"]])

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    dims = get_dimensions("umap_wide", aes)

    fig, ax = plt.subplots(figsize=(dims["width"], dims["height"]))
    ax.scatter(df["UMAP_1"], df["UMAP_2"], c=color_arr,
               s=0.5, linewidths=0, alpha=0.5, rasterized=True)

    # Centroid labels for contamination clusters only (keep panel clean)
    for lbl in CONTAMINATION_LABELS:
        sub = df[df["cluster_annotation"] == lbl]
        if len(sub) == 0:
            continue
        cx = sub["UMAP_1"].median()
        cy = sub["UMAP_2"].median()
        ax.text(cx, cy, lbl.replace("_", " "), fontsize=6, weight="bold",
                ha="center", va="center", color="black",
                bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                          edgecolor="none", alpha=0.85))

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    handles = []
    for lbl in ANNOTATION_ORDER:
        if lbl in pal:
            handles.append(Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor=pal[lbl], markersize=4,
                                  label=lbl.replace("_", " ")))
    ax.legend(handles=handles, loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=5,
              handletextpad=0.5, labelspacing=0.4, ncol=1)

    out_pdf = args.out_dir / "supp2_flex_contamination_umap.pdf"
    out_csv = args.out_dir / "supp2_flex_contamination_umap_data.csv"

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)

    df[["cell_id", "UMAP_1", "UMAP_2", "cluster_annotation",
        "compartment"]].to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(df))


if __name__ == "__main__":
    main()
