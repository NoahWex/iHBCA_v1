"""
render_supp1p5_flex_ncount_umap.py — Fig 3 Supp 1.5, Panel 1.5b.

Cohort step-08 UMAP feature plot colored by log10(n_umi). Pairs with 1.5a
(cell-type UMAP) to show that high-nCount cells are not localized to one
type but cluster in the heterogeneous problem region. Cell-level UMI
filtering can't pick them out without removing real biology.

Color scale: viridis on log10 nCount, clipped to 5–99.5th pct so a small
number of extreme outliers don't dominate.

Reviewer question (Act I of Supp 1.5): "Are these cells just high-count
outliers?"
The visible answer: no — the high-nCount density spans the manifold and
concentrates in the heterogeneous region from 1.5a, but cell-level UMI
filtering would either keep them all or remove real biology.

Inputs (consumes-from-dev):
    Spatial_HBCA_preprocessing/outputs/preprocessing/08_ScviIntegration/
        embeddings/umap_2d.csv (rownames=cell_id, UMAP_1, UMAP_2)
    Spatial_HBCA_preprocessing/outputs/preprocessing/02_CellFiltering/
        cell_metadata/<sample>_metadata.csv (cell_id, n_umi)

Outputs:
    supp1p5_flex_ncount_umap.pdf
    supp1p5_flex_ncount_umap_data.csv

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
log = logging.getLogger("supp1p5_1b")


SHUFFLE_SEED = 42
CMAP = "viridis"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K cells for sanity render.")
    return p.parse_args()


def load_n_umi(meta_dir: Path) -> pd.Series:
    csvs = sorted(meta_dir.glob("*_metadata.csv"))
    log.info("Loading %d sample CSVs from %s", len(csvs), meta_dir)
    frames = []
    for f in csvs:
        df = pd.read_csv(f, index_col=0, usecols=["Unnamed: 0", "n_umi"])
        df.index.name = "cell_id"
        frames.append(df)
    out = pd.concat(frames)
    log.info("Loaded n_umi for %d cells", len(out))
    return out["n_umi"]


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    paths = load_paths(config_dir=str(args.project_root / "publication" / "config"))
    umap_path = Path(resolve_path("flex_scvi_n100_umap", paths))
    meta_dir = Path(resolve_path("flex_cell_metadata_dir", paths))

    log.info("Reading UMAP: %s", umap_path)
    umap = pd.read_csv(umap_path, index_col=0)
    umap.index.name = "cell_id"
    log.info("UMAP: %d cells", len(umap))

    n_umi = load_n_umi(meta_dir)

    df = umap.join(n_umi.rename("n_umi"), how="inner").reset_index()
    df["log10_n_umi"] = np.log10(df["n_umi"].clip(lower=1))
    log.info("Joined: %d cells", len(df))

    if args.test:
        df = df.sample(n=50_000, random_state=SHUFFLE_SEED).reset_index(drop=True)
        log.info("Subsampled to %d cells (test)", len(df))

    rng = np.random.default_rng(SHUFFLE_SEED)
    perm = rng.permutation(len(df))
    df = df.iloc[perm].reset_index(drop=True)

    vals = df["log10_n_umi"].to_numpy()
    vmin = float(np.nanpercentile(vals, 5))
    vmax = float(np.nanpercentile(vals, 99.5))
    log.info("log10(n_umi) clip: [%.2f, %.2f]", vmin, vmax)

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    dims = get_dimensions("umap", aes)

    fig, ax = plt.subplots(figsize=(dims["width"], dims["height"]))
    sc = ax.scatter(df["UMAP_1"], df["UMAP_2"], c=vals, cmap=CMAP,
                    vmin=vmin, vmax=vmax, s=0.5, linewidths=0, alpha=0.5,
                    rasterized=True)

    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(r"log$_{10}$(n UMI)", fontsize=6)
    cb.ax.tick_params(labelsize=5)

    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    out_pdf = args.out_dir / "supp1p5_flex_ncount_umap.pdf"
    out_csv = args.out_dir / "supp1p5_flex_ncount_umap_data.csv"

    fig.tight_layout()
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)

    df[["cell_id", "UMAP_1", "UMAP_2", "n_umi", "log10_n_umi"]].to_csv(
        out_csv, index=False)
    log.info("Wrote %s (%d rows)", out_csv, len(df))


if __name__ == "__main__":
    main()
