"""render_supp8a_joint_umap_per_compartment.py — Fig 3 Supp 8a: per-compartment
joint embedding label UMAP.

One panel per compartment (Epithelial / Immune / Stromal). Reads the joint
(FLEX + Xenium) UMAP coordinates and the L1.5 cascade labels produced by the
per-compartment annotation pipeline (publication/analysis/annotation/xenium/
outputs/compartment/{C}/{joint_umap.csv,cell_annotations.csv}). Renders a
rasterized scatter inside a vector PDF, with centroid text labels and
force-directed repulsion (no external legend).

Pattern reference (full-script adaptation):
- publication/figures/render/flex/render_supp5b_label_umap.py — entire script
  is the structural template. Adapted by:
    * Swapping FLEX integration intermediate (umap.csv + obs.csv joined with
      flex_l2s_labels.csv) for the per-compartment joint substrate
      (compartment_{epi,imm,str}_{joint_umap,cell_annotations} tokens).
    * Reading the canonical label column directly from cell_annotations.csv
      (`label`, with `is_artifact` flag) rather than parsing label strings
      for an ARTIFACT_ prefix.
    * Using UMAP column names UMAP1/UMAP2 (no underscore — joint substrate
      schema) instead of UMAP_1/UMAP_2.

No on-plot title, panel letter, or method caption — Illustrator handles those.

Outputs:
    supp8a_joint_umap_{epithelial,immune,stromal}.pdf
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List, Tuple

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp8a")


COMPARTMENT_FULL = {"epi": "Epithelial", "imm": "Immune", "str": "Stromal"}
COMPARTMENT_TOKEN = {
    "epi": ("compartment_epi_joint_umap", "compartment_epi_cell_annotations"),
    "imm": ("compartment_imm_joint_umap", "compartment_imm_cell_annotations"),
    "str": ("compartment_str_joint_umap", "compartment_str_cell_annotations"),
}
ARTIFACT_COLOR = "#CCCCCC"

# Per-compartment artifact reclassifications applied at render time. Labels
# listed here render as artifact-class even when the cascade vocabulary lists
# them under canonical labels. BMYO-NC is reclassified because the cluster is
# sample-specific (single patient, single Xenium sample) and does not
# represent a generalizable biological state.
ARTIFACT_OVERRIDES = {
    "epi": {"BMYO-NC"},
    "imm": set(),
    "str": set(),
}


def _setup_aesthetics(project_root: Path) -> Tuple[dict, dict]:
    """Load aesthetics framework + return (config, dimensions).

    Pattern: render_supp5b_label_umap.py:55-67.
    """
    sys.path.insert(0, str(project_root / "publication" / "config"))
    from load_aesthetics import (  # type: ignore
        load_aesthetics, get_matplotlib_theme, get_dimensions,
    )
    cfg = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(cfg))
    dims = get_dimensions("umap", cfg)
    return cfg, dims


def categorical_palette(n: int) -> List[Tuple[float, float, float]]:
    """Tab20-family palette extended with HSV. Pattern: supp5b:70-80."""
    base = (
        list(colormaps["tab20"].colors)
        + list(colormaps["tab20b"].colors)
        + list(colormaps["tab20c"].colors)
    )
    if n <= len(base):
        return base[:n]
    extra = colormaps["hsv"](np.linspace(0, 1, n - len(base) + 1))[:-1]
    return base + [tuple(c[:3]) for c in extra]


def repel_labels(
    anchors: np.ndarray, x_range: float, y_range: float,
    iters: int = 160, anchor_pull: float = 0.04, sep_frac: float = 0.075,
    repel_frac: float = 0.0055,
) -> np.ndarray:
    """Force-directed repulsion fallback. Pattern: supp5b:83-111."""
    if len(anchors) == 0:
        return anchors.copy()
    diag = (x_range ** 2 + y_range ** 2) ** 0.5
    min_sep = diag * sep_frac
    repel_strength = diag * repel_frac
    positions = anchors.copy()
    for _ in range(iters):
        disp = np.zeros_like(positions)
        for i in range(len(positions)):
            for j in range(len(positions)):
                if i == j:
                    continue
                d = positions[i] - positions[j]
                dist = float(np.linalg.norm(d))
                if dist < min_sep:
                    if dist < 1e-6:
                        d = np.array([1.0, 1.0]) * 1e-3
                        dist = float(np.linalg.norm(d))
                    push = (min_sep - dist) / min_sep
                    disp[i] += (d / dist) * push * repel_strength
            disp[i] += (anchors[i] - positions[i]) * anchor_pull
        positions += disp
    return positions


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--project-root", type=Path, required=True,
                   help="iHBCA_publication root (resolves paths.yaml + aesthetics.yaml)")
    p.add_argument("--compartment", required=True, choices=["epi", "imm", "str"])
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--test", action="store_true",
                   help="Subsample to 5000 cells for fast iteration")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    cfg, _dims = _setup_aesthetics(project_root)

    sys.path.insert(0, str(project_root / "publication" / "config"))
    from load_paths import load_paths, resolve_path  # type: ignore
    paths = load_paths(str(project_root / "publication" / "config"))

    compartment_full = COMPARTMENT_FULL[args.compartment]
    umap_token, ann_token = COMPARTMENT_TOKEN[args.compartment]
    umap_csv = Path(resolve_path(umap_token, paths))
    ann_csv = Path(resolve_path(ann_token, paths))

    log.info("UMAP CSV  : %s", umap_csv)
    log.info("Annot CSV : %s", ann_csv)
    for p in (umap_csv, ann_csv):
        if not p.exists():
            sys.exit(f"FATAL: missing input: {p}")

    umap = pd.read_csv(umap_csv)
    if "UMAP1" not in umap.columns or "UMAP2" not in umap.columns:
        sys.exit(f"FATAL: expected UMAP1/UMAP2 columns in {umap_csv}; got {list(umap.columns)}")
    log.info("UMAP rows: %d", len(umap))

    ann = pd.read_csv(ann_csv, usecols=["cell_id", "label", "is_artifact"])
    ann["is_artifact"] = ann["is_artifact"].astype(bool)
    log.info("Annotation rows: %d (artifacts=%d)",
             len(ann), int(ann["is_artifact"].sum()))

    # Apply ARTIFACT_OVERRIDES so reclassified labels render in the artifact
    # group (gray points + numbered marker outlined rather than filled).
    overrides = ARTIFACT_OVERRIDES.get(args.compartment, set())
    if overrides:
        n_pre = int(ann["is_artifact"].sum())
        ann.loc[ann["label"].isin(overrides), "is_artifact"] = True
        n_post = int(ann["is_artifact"].sum())
        log.info("Override applied: %s → +%d cells reclassified as artifact",
                 sorted(overrides), n_post - n_pre)

    df = umap.merge(ann, on="cell_id", how="inner")
    log.info("Joined: %d cells", len(df))

    if args.test:
        df = df.sample(n=min(5000, len(df)), random_state=0).reset_index(drop=True)
        log.info("--test: subsampled to %d cells", len(df))

    df["is_artifact"] = df["is_artifact"].astype(bool)

    nonart = df[~df["is_artifact"]]
    label_counts = nonart.groupby("label").size().sort_values(ascending=False)
    label_order = label_counts.index.tolist()
    palette = categorical_palette(len(label_order))
    color_map = {lbl: palette[i] for i, lbl in enumerate(label_order)}

    # Scatter aesthetic ported from publication/analysis/annotation/xenium/
    # scripts/02_compartment_analysis.py:41-58 (scatter_cat): s=1.0,
    # alpha=0.3 manages density via low opacity rather than dense paint.
    # Centroid numbered markers + right-side legend match the same upstream
    # legend style (small swatch + label text).
    pt_size = 1.0
    alpha = 0.3
    # Panel: scatter on the left + legend strip on the right. Total ~75mm wide.
    fig = plt.figure(figsize=(75.0 / 25.4, 60.0 / 25.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[3.5, 1], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    legend_ax = fig.add_subplot(gs[0, 1])
    ax.set_aspect("equal")
    ax.axis("off")
    legend_ax.axis("off")

    art = df[df["is_artifact"]]
    if len(art):
        ax.scatter(
            art["UMAP1"], art["UMAP2"],
            c=ARTIFACT_COLOR, s=pt_size, alpha=alpha * 0.5,
            rasterized=True, linewidths=0,
        )
    for lbl in label_order:
        sub = nonart[nonart["label"] == lbl]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP1"], sub["UMAP2"],
            c=[color_map[lbl]], s=pt_size, alpha=alpha,
            rasterized=True, linewidths=0,
        )

    canonical_centroids = (
        nonart.groupby("label", as_index=False)
              .agg(UMAP1=("UMAP1", "median"), UMAP2=("UMAP2", "median"),
                   n_cells=("cell_id", "size"))
    )
    canonical_centroids = canonical_centroids[canonical_centroids["n_cells"] >= 50] \
                              .reset_index(drop=True)
    canonical_centroids = canonical_centroids.sort_values(
        "n_cells", ascending=False).reset_index(drop=True)

    artifact_centroids = (
        art.groupby("label", as_index=False)
           .agg(UMAP1=("UMAP1", "median"), UMAP2=("UMAP2", "median"),
                n_cells=("cell_id", "size"))
    )
    artifact_centroids = artifact_centroids[artifact_centroids["n_cells"] >= 50] \
                              .reset_index(drop=True)
    artifact_centroids = artifact_centroids.sort_values(
        "n_cells", ascending=False).reset_index(drop=True)

    # Number assignment: canonical 1..N, artifacts (N+1)..M.
    legend_entries: List[Tuple[int, str, bool, Tuple[float, float, float]]] = []
    n = 1
    for _, row in canonical_centroids.iterrows():
        col = color_map.get(row["label"], (0.4, 0.4, 0.4))
        legend_entries.append((n, str(row["label"]), False, col))
        n += 1
    for _, row in artifact_centroids.iterrows():
        legend_entries.append((n, str(row["label"]), True,
                               (0.55, 0.55, 0.55)))
        n += 1

    # Centroid numbered markers — small disks with the legend index.
    # Sizing matches upstream legend (markersize=4 on Line2D ≈ scatter s≈40)
    # rather than the prior bulky s=110.
    def _draw_centroid(x: float, y: float, num: int, is_artifact: bool) -> None:
        face = "white"
        edge = "#666666" if is_artifact else "black"
        text_color = "#666666" if is_artifact else "black"
        ax.scatter([x], [y], s=42, marker="o",
                   facecolors=face, edgecolors=edge, linewidths=0.4,
                   zorder=10)
        ax.text(x, y, str(num), fontsize=3.5, ha="center", va="center",
                color=text_color, weight="bold", zorder=11)

    for i, row in canonical_centroids.iterrows():
        _draw_centroid(float(row["UMAP1"]), float(row["UMAP2"]),
                       i + 1, is_artifact=False)
    offset = len(canonical_centroids)
    for i, row in artifact_centroids.iterrows():
        _draw_centroid(float(row["UMAP1"]), float(row["UMAP2"]),
                       offset + i + 1, is_artifact=True)

    # Right-side legend — upstream scatter_cat style (small swatch + label).
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    n_entries = len(legend_entries)
    if n_entries > 0:
        line_h = 1.0 / max(n_entries + 1, 1)
        for i, (num, lbl, is_art, color) in enumerate(legend_entries):
            y = 1.0 - (i + 0.5) * line_h
            legend_ax.scatter([0.05], [y], s=22, marker="o",
                              facecolors=("white" if is_art else color),
                              edgecolors=("#666666" if is_art else "black"),
                              linewidths=0.3, transform=legend_ax.transAxes,
                              clip_on=False)
            legend_ax.text(0.05, y, str(num), fontsize=3.0,
                           ha="center", va="center",
                           color=("#666666" if is_art else "black"),
                           weight="bold", transform=legend_ax.transAxes)
            legend_ax.text(0.13, y, lbl, fontsize=4.5,
                           ha="left", va="center",
                           color=("#666666" if is_art else "black"),
                           style=("italic" if is_art else "normal"),
                           transform=legend_ax.transAxes)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = args.out_dir / f"supp8a_joint_umap_{compartment_full.lower()}.pdf"
    fig.savefig(out_pdf, bbox_inches="tight",
                dpi=int(cfg.get("rendering", {}).get("dpi_umap_raster", 600)))
    plt.close(fig)
    log.info("Wrote %s (%d canonical + %d artifact labels)",
             out_pdf, len(canonical_centroids), len(artifact_centroids))
    return 0


if __name__ == "__main__":
    sys.exit(main())
