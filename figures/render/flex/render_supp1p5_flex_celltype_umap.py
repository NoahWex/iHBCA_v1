"""
render_supp1p5_flex_celltype_umap.py — Fig 3 Supp 1.5, Panel 1.5a.

Cohort step-08 UMAP feature plot colored by Kumar SingleR cell-type labels
(`HBCATransferredLabels.Kumar_2023`, transferred at step 04). Shows the
mature post-MAD/doublet manifold structure — the heterogeneous problem
region (where multiple cell types co-cluster despite intervening filters)
is visually identifiable here, motivating the BigSur manifold filter.

Reviewer question (Act I of Supp 1.5): "Where do canonical cell types live
on this manifold?"

Inputs (consumes-from-dev):
    Spatial_HBCA_preprocessing/outputs/preprocessing/08_ScviIntegration/
        embeddings/umap_2d.csv (rownames=cell_id, UMAP_1, UMAP_2)
    Spatial_HBCA_preprocessing/outputs/preprocessing/04_LabelTransfer/
        cell_metadata/<sample>_labels.csv (cell_id, HBCATransferredLabels.Kumar_2023)

Outputs:
    supp1p5_flex_celltype_umap.pdf
    supp1p5_flex_celltype_umap_data.csv

Substrate gap:
    No `kumar_l1` palette token in publication/config/aesthetics.yaml. This
    script uses a placeholder mapping anchored to compartment colors:
        Epithelial: Basal, Luminal_HR, Luminal_Secretory
        Stromal:    Fibroblast, Pericyte, Vascular_Endothelium,
                    Lymphatic_Endothelium
        Immune:     Myeloid, T_Cell, B_Cell, Mast
    Each compartment's types use distinct hues from a qualitative ramp,
    flagged inline. Coordinator should add a formal `kumar_l1` palette
    token with vetted colors.

Pattern source:
    publication/figures/render/flex/render_supp1p5_flex_delta_umap.py
    (rasterized scatter scaffolding)

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
    get_palette,
    load_aesthetics,
)
from load_paths import load_paths, resolve_path  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp1p5_1a")


SHUFFLE_SEED = 42
LABEL_COL = "HBCATransferredLabels.Kumar_2023"

LABEL_COMPARTMENT = {
    "Basal":                 "Epithelial",
    "Luminal_HR":            "Epithelial",
    "Luminal_Secretory":     "Epithelial",
    "Fibroblast":            "Stromal",
    "Pericyte":              "Stromal",
    "Vascular_Endothelium":  "Stromal",
    "Lymphatic_Endothelium": "Stromal",
    "Adipocyte":             "Stromal",
    "Myeloid":               "Immune",
    "T_Cell":                "Immune",
    "B_Cell":                "Immune",
    "Mast":                  "Immune",
}

LABEL_ORDER = [
    "Basal", "Luminal_HR", "Luminal_Secretory",
    "Fibroblast", "Pericyte", "Vascular_Endothelium", "Lymphatic_Endothelium",
    "Adipocyte",
    "Myeloid", "T_Cell", "B_Cell", "Mast",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K cells for sanity render.")
    return p.parse_args()


def load_kumar_labels(meta_dir: Path) -> pd.Series:
    csvs = sorted(meta_dir.glob("*_labels.csv"))
    log.info("Loading %d label CSVs from %s", len(csvs), meta_dir)
    frames = []
    for f in csvs:
        df = pd.read_csv(f, index_col=0, usecols=["Unnamed: 0", LABEL_COL])
        df.index.name = "cell_id"
        frames.append(df)
    out = pd.concat(frames)
    log.info("Loaded %d cell labels", len(out))
    return out[LABEL_COL]


def get_kumar_palette(config: dict) -> dict:
    """Kumar L1 palette via aesthetics.yaml `kumar_l1` token."""
    return get_palette("kumar_l1", config)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    umap_path = Path(resolve_path("flex_scvi_n100_umap", paths))
    label_dir = Path(resolve_path("flex_kumar_labels_dir", paths))

    log.info("Reading UMAP: %s", umap_path)
    umap = pd.read_csv(umap_path, index_col=0)
    umap.index.name = "cell_id"
    log.info("UMAP: %d cells", len(umap))

    labels = load_kumar_labels(label_dir)
    df = umap.join(labels.rename("kumar_label"), how="inner").reset_index()
    df = df.dropna(subset=["kumar_label"])
    log.info("Joined: %d cells; label counts: %s",
             len(df), df["kumar_label"].value_counts().to_dict())

    if args.test:
        df = df.sample(n=50_000, random_state=SHUFFLE_SEED).reset_index(drop=True)
        log.info("Subsampled to %d cells (test)", len(df))

    rng = np.random.default_rng(SHUFFLE_SEED)
    perm = rng.permutation(len(df))
    df = df.iloc[perm].reset_index(drop=True)

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    pal = get_kumar_palette(aes)
    dims = get_dimensions("umap_wide", aes)

    color_arr = np.array([pal[lbl] for lbl in df["kumar_label"]])

    fig, ax = plt.subplots(figsize=(dims["width"], dims["height"]))
    ax.scatter(df["UMAP_1"], df["UMAP_2"], c=color_arr,
               s=0.5, linewidths=0, alpha=0.5, rasterized=True)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    handles = []
    for label in LABEL_ORDER:
        if label in pal:
            handles.append(Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor=pal[label], markersize=4,
                                  label=label.replace("_", " ")))
    ax.legend(handles=handles, loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=False, fontsize=5,
              handletextpad=0.5, labelspacing=0.4)

    out_pdf = args.out_dir / "supp1p5_flex_celltype_umap.pdf"
    out_csv = args.out_dir / "supp1p5_flex_celltype_umap_data.csv"

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)

    df[["cell_id", "UMAP_1", "UMAP_2", "kumar_label"]].to_csv(
        out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(df))


if __name__ == "__main__":
    main()
