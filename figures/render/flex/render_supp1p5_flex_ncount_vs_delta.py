"""
render_supp1p5_flex_ncount_vs_delta.py — Fig 3 Supp 1.5, Panel 1.5d.

Per-cell scatter of log10(nCount_RNA) vs Δ (delta_mad_scaled). Demonstrates
that the manifold-incoherence metric is independent of read depth: cells
flagged by step_10b are not preferentially low- or high-count cells, they
span the full nCount range. Cell-level UMI/MAD filtering cannot find them.

Layout:
- Background: hexbin of all 275,542 cells (cohort density in nCount × Δ
  space) at low alpha
- Foreground: highlighted scatter of removed cells colored by removal_reason
  (Global_Incoherence_GMM vs TrashCluster_FinePocket)

Reviewer question (Act II of Supp 1.5): "Does Δ just track depth?"
The visible answer: no — failed cells span 3 orders of magnitude in nCount;
Δ is a depth-orthogonal metric.

Inputs (consumes-from-dev — bundle into Supp 1.5 promotion package):
    Spatial_HBCA_preprocessing/outputs/preprocessing/10b_InteractiveExploration/
        exports/full_decision_matrix.csv
    cols: cell_id, nCount_RNA, delta_mad_scaled, removal_reason

Outputs:
    supp1p5_flex_ncount_vs_delta.pdf
    supp1p5_flex_ncount_vs_delta_data.csv

Pattern source:
    publication/figures/render/flex/render_supp1_flex_filter_umap.py
    (Option B framework boilerplate)

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
log = logging.getLogger("supp1p5_1d")


SHUFFLE_SEED = 42
GRIDSIZE = 80


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K cells for sanity render.")
    return p.parse_args()


def get_status_palette(config: dict) -> dict:
    """Binary kept/removed via qc_status."""
    qc = get_palette("qc_status", config)
    return {
        "kept": qc["kept"],
        "removed": qc["removed_doublet"],
    }


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    dm_path = Path(resolve_path("flex_manifold_decision_matrix", paths))
    log.info("Reading 10b decision matrix: %s", dm_path)
    df = pd.read_csv(dm_path, low_memory=False,
                     usecols=["cell_id", "nCount_RNA", "delta_mad_scaled",
                              "removal_reason"])
    log.info("Loaded %d cells", len(df))

    df = df.dropna(subset=["nCount_RNA", "delta_mad_scaled"])
    df["log10_nCount"] = np.log10(df["nCount_RNA"].clip(lower=1))
    log.info("Removal counts: %s", df["removal_reason"].value_counts().to_dict())

    if args.test:
        df = df.sample(n=50_000, random_state=SHUFFLE_SEED).reset_index(drop=True)
        log.info("Subsampled to %d cells (test mode)", len(df))

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    status_pal = get_status_palette(aes)
    dims = get_dimensions("biplot", aes)

    fig, ax = plt.subplots(figsize=(dims["width"], dims["height"]))

    hb = ax.hexbin(df["log10_nCount"], df["delta_mad_scaled"],
                   gridsize=GRIDSIZE, cmap="Greys", mincnt=1, bins="log",
                   linewidths=0, alpha=0.85, rasterized=True)
    cb = fig.colorbar(hb, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("Cells (log)", fontsize=6)
    cb.ax.tick_params(labelsize=5)

    removed = df[df["removal_reason"] != "Pass"].copy()
    rng = np.random.default_rng(SHUFFLE_SEED)
    perm = rng.permutation(len(removed))
    removed = removed.iloc[perm].reset_index(drop=True)
    ax.scatter(removed["log10_nCount"], removed["delta_mad_scaled"],
               c=status_pal["removed"], s=1.2, linewidths=0, alpha=0.6,
               rasterized=True)

    ax.axhline(0, color="black", linestyle=":", linewidth=0.6, alpha=0.5)

    ax.set_xlabel(r"log$_{10}$(nCount RNA)")
    ax.set_ylabel(r"$\Delta$ (MAD-scaled)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    n_removed = int((df["removal_reason"] != "Pass").sum())
    handles = [
        Line2D([0], [0], marker="h", color="w",
               markerfacecolor="#666666", markersize=5,
               label=f"All cells (n={len(df):,})"),
        Line2D([0], [0], marker="o", color="w",
               markerfacecolor=status_pal["removed"], markersize=4,
               label=f"Removed (n={n_removed:,})"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=5,
              handletextpad=0.5)

    out_pdf = args.out_dir / "supp1p5_flex_ncount_vs_delta.pdf"
    out_csv = args.out_dir / "supp1p5_flex_ncount_vs_delta_data.csv"

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)

    df[["cell_id", "nCount_RNA", "log10_nCount", "delta_mad_scaled",
        "removal_reason"]].to_csv(out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(df))


if __name__ == "__main__":
    main()
