"""
render_supp1p5_flex_delta_umap.py — Fig 3 Supp 1.5, Panel 1.5e.

Cohort step-08 UMAP feature plot colored by Δ (delta_mad_scaled). Shows
where on the integrated manifold the BigSur incoherence signal lights up.
Visually paired with 1.5a (cell-type UMAP) — Δ should localize to the
heterogeneous problem region.

Color scale: divergent blue-orange centered at 0. Negative Δ = artifact-
candidate (smoothed signal collapsed). Positive Δ ≈ neighborhood
reinforcement (no signal of incoherence). Clipped to symmetric range
[-q99, +q99] so single-cell extreme outliers don't dominate.

Reviewer question (Act II of Supp 1.5): "Where on the manifold does Δ
fire?"

Inputs (consumes-from-dev):
    Spatial_HBCA_preprocessing/outputs/preprocessing/08_ScviIntegration/
        embeddings/umap_2d.csv (rownames=cell_id, UMAP_1, UMAP_2)
    Spatial_HBCA_preprocessing/outputs/preprocessing/10b_InteractiveExploration/
        exports/full_decision_matrix.csv (cell_id, delta_mad_scaled)

Outputs:
    supp1p5_flex_delta_umap.pdf
    supp1p5_flex_delta_umap_data.csv

Pattern source:
    publication/figures/render/flex/render_supp1_flex_filter_umap.py
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

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
    get_dimensions,
    get_matplotlib_theme,
    load_aesthetics,
)
from load_paths import load_paths, resolve_path  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp1p5_1e")


SHUFFLE_SEED = 42
CMAP = "RdBu_r"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K cells for sanity render.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    umap_path = Path(resolve_path("flex_scvi_n100_umap", paths))
    dm_path = Path(resolve_path("flex_manifold_decision_matrix", paths))

    log.info("Reading UMAP: %s", umap_path)
    umap = pd.read_csv(umap_path, index_col=0)
    umap.index.name = "cell_id"
    log.info("UMAP: %d cells", len(umap))

    log.info("Reading decision matrix: %s", dm_path)
    dm = pd.read_csv(dm_path, low_memory=False,
                     usecols=["cell_id", "delta_mad_scaled"]).set_index("cell_id")
    log.info("Decision matrix: %d cells", len(dm))

    df = umap.join(dm, how="inner").reset_index()
    df = df.dropna(subset=["delta_mad_scaled"])
    log.info("Joined: %d cells", len(df))

    if args.test:
        df = df.sample(n=50_000, random_state=SHUFFLE_SEED).reset_index(drop=True)
        log.info("Subsampled to %d cells (test)", len(df))

    rng = np.random.default_rng(SHUFFLE_SEED)
    perm = rng.permutation(len(df))
    df = df.iloc[perm].reset_index(drop=True)

    delta = df["delta_mad_scaled"].to_numpy()
    q99 = float(np.nanpercentile(np.abs(delta), 99))
    log.info("Δ symmetric clip: ±%.2f (99th pct of |Δ|)", q99)

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    dims = get_dimensions("umap", aes)

    fig, ax = plt.subplots(figsize=(dims["width"], dims["height"]))
    sc = ax.scatter(df["UMAP_1"], df["UMAP_2"], c=delta, cmap=CMAP,
                    vmin=-q99, vmax=q99, s=0.5, linewidths=0, alpha=0.5,
                    rasterized=True)

    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(r"$\Delta$ (MAD-scaled)", fontsize=6)
    cb.ax.tick_params(labelsize=5)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    out_pdf = args.out_dir / "supp1p5_flex_delta_umap.pdf"
    out_csv = args.out_dir / "supp1p5_flex_delta_umap_data.csv"

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)

    df[["cell_id", "UMAP_1", "UMAP_2", "delta_mad_scaled"]].to_csv(
        out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(df))


if __name__ == "__main__":
    main()
