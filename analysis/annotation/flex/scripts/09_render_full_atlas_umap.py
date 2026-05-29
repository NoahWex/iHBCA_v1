"""09_render_full_atlas_umap.py — project all FLEX L2S labels onto the integrated atlas UMAP.

Pre-promotion sanity panel: confirms per-compartment L2S labels (concatenated to
flex_l2s_labels.csv) project sensibly onto the joint atlas UMAP — i.e. compartment
boundaries are tight and lineage neighborhoods (LASP near LHS, Vas-arterial near
Vas-vein, etc.) emerge spatially as expected.

Inputs:
  - {full_bundle}/umap.csv          full-atlas UMAP from scvi_n100 integration
  - {full_bundle}/obs.csv           full-atlas obs (used as fallback for cell_id)
  - {labels_csv}                    concatenated FLEX labels (cell_id × l2s_label)

Outputs:
  - {out_dir}/flex_l2s_full_atlas_umap.pdf
  - {out_dir}/flex_l2s_full_atlas_umap.png            (preview raster)
  - {out_dir}/flex_l2s_full_atlas_umap_by_compartment.pdf  (3-panel facet)

Styling matches the per-compartment label_umap.pdf produced by 06_render_annotation.py
(centroid text, force-directed label repulsion via adjustText, categorical palette).
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path
from typing import Dict, List, Tuple

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
log = logging.getLogger("full_atlas_umap")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--full-bundle", type=Path, required=True,
                   help="dir with umap.csv + obs.csv (full atlas scvi_n100)")
    p.add_argument("--labels-csv", type=Path, required=True,
                   help="concatenated flex_l2s_labels.csv")
    p.add_argument("--label-column", default="l2s_label",
                   help="column in labels-csv to color by")
    p.add_argument("--compartment-column", default="compartment")
    p.add_argument("--out-dir", type=Path, required=True)
    return p.parse_args()


def categorical_palette(n: int) -> List[tuple]:
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
    label_pos: List[Tuple[str, float, float, float, float]],
    df: pd.DataFrame,
    n_iter: int = 120,
) -> List[Tuple[str, float, float, float, float]]:
    if not label_pos:
        return label_pos
    anchors = np.array([(x, y) for _, x, y, _, _ in label_pos])
    positions = anchors.copy()
    x_range = float(df["umap_1"].max() - df["umap_1"].min())
    y_range = float(df["umap_2"].max() - df["umap_2"].min())
    diag = (x_range ** 2 + y_range ** 2) ** 0.5
    min_sep = diag * 0.045
    anchor_pull = 0.05
    repel_strength = diag * 0.0035
    for _ in range(n_iter):
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
    out = []
    for i, (lbl, ax_x, ax_y, _, _) in enumerate(label_pos):
        out.append((lbl, ax_x, ax_y, float(positions[i, 0]), float(positions[i, 1])))
    return out


def assign_colors(labels: List[str]) -> Dict[str, tuple]:
    bio = [l for l in labels if not str(l).startswith("ARTIFACT")]
    bio.sort()
    palette = categorical_palette(len(bio))
    color_map: Dict[str, tuple] = dict(zip(bio, palette))
    artifact_grey = (0.78, 0.78, 0.78)
    for l in labels:
        if str(l).startswith("ARTIFACT"):
            color_map[l] = artifact_grey
    return color_map


def render_panel(
    ax,
    df: pd.DataFrame,
    color_map: Dict[str, tuple],
    point_size: float = 1.2,
    alpha: float = 0.85,
    label_min_cells: int = 200,
    title: str | None = None,
) -> None:
    colors = df["label"].map(color_map).tolist()
    ax.scatter(df["umap_1"], df["umap_2"], c=colors, s=point_size, alpha=alpha,
               rasterized=True, linewidths=0)

    label_pos: List[Tuple[str, float, float, float, float]] = []
    for lbl, sub in df.groupby("label"):
        if len(sub) < label_min_cells:
            continue
        x = float(sub["umap_1"].median())
        y = float(sub["umap_2"].median())
        label_pos.append((lbl, x, y, x, y))
    label_pos = repel_labels(label_pos, df)

    for lbl, ax_x, ax_y, lx, ly in label_pos:
        if (lx != ax_x) or (ly != ax_y):
            ax.plot([ax_x, lx], [ax_y, ly], color="gray", lw=0.25, alpha=0.5, zorder=9)
        ax.text(
            lx, ly, lbl,
            fontsize=6, ha="center", va="center", color="black", weight="bold",
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="black", lw=0.5, alpha=0.92),
            zorder=10,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.3)
    if title is not None:
        ax.set_title(title, fontsize=8)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    umap_csv = args.full_bundle / "umap.csv"
    obs_csv = args.full_bundle / "obs.csv"
    log.info("Loading UMAP: %s", umap_csv)
    umap = pd.read_csv(umap_csv)
    umap.columns = [c.lower() for c in umap.columns]
    if "cell_id" not in umap.columns and obs_csv.exists():
        log.info("Falling back to obs.csv for cell_id ordering")
        obs = pd.read_csv(obs_csv)
        umap = obs[["cell_id", "UMAP_1", "UMAP_2"]].copy()
        umap.columns = ["cell_id", "umap_1", "umap_2"]
    if not {"umap_1", "umap_2"}.issubset(umap.columns):
        cs = [c for c in umap.columns if c != "cell_id"]
        umap = umap.rename(columns={cs[0]: "umap_1", cs[1]: "umap_2"})

    log.info("Loading labels: %s", args.labels_csv)
    labels_df = pd.read_csv(args.labels_csv)
    if "cell_id" not in labels_df.columns:
        raise SystemExit("ERROR: labels CSV missing cell_id column")
    if args.label_column not in labels_df.columns:
        raise SystemExit(f"ERROR: labels CSV missing {args.label_column} column "
                         f"(have: {list(labels_df.columns)})")

    log.info("UMAP rows: %d, label rows: %d", len(umap), len(labels_df))
    df = umap.merge(
        labels_df[["cell_id", args.label_column, args.compartment_column]],
        on="cell_id",
        how="inner",
    )
    df = df.rename(columns={args.label_column: "label",
                            args.compartment_column: "compartment"})
    log.info("Merged: %d cells (%d dropped from UMAP, %d dropped from labels)",
             len(df), len(umap) - len(df), len(labels_df) - len(df))

    n_pre = len(df)
    df = df[~df["label"].astype(str).str.startswith("ARTIFACT")].copy()
    log.info("Dropped %d ARTIFACT_* cells from projection (kept %d biological)",
             n_pre - len(df), len(df))

    counts = df["label"].value_counts()
    log.info("L2S label counts (top 15):\n%s", counts.head(15).to_string())
    n_labels = counts.shape[0]
    log.info("Total L2S labels in plot: %d", n_labels)

    labels_in_data = sorted(df["label"].dropna().unique())
    color_map = assign_colors(labels_in_data)

    # ---- main panel: full atlas, all biological labels ----
    fig, ax = plt.subplots(figsize=(12, 12))
    render_panel(ax, df, color_map, point_size=0.8, alpha=0.7,
                 label_min_cells=300, title=None)
    out_pdf = args.out_dir / "flex_l2s_full_atlas_umap.pdf"
    out_png = args.out_dir / "flex_l2s_full_atlas_umap.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)
    log.info("  wrote %s", out_png)

    # ---- 3-panel facet: by compartment (active compartment colored, others grey) ----
    compartments = sorted(df["compartment"].dropna().unique())
    fig, axes = plt.subplots(1, len(compartments), figsize=(6 * len(compartments), 6.5),
                             sharex=True, sharey=True)
    if len(compartments) == 1:
        axes = [axes]
    for ax_i, comp in zip(axes, compartments):
        sub = df[df["compartment"].astype(str) == comp].copy()
        comp_color_map = {l: color_map.get(l, (0.4, 0.4, 0.4)) for l in sub["label"].unique()}
        render_panel(ax_i, sub, comp_color_map, point_size=0.6, alpha=0.7,
                     label_min_cells=250, title=str(comp))
    out_pdf2 = args.out_dir / "flex_l2s_full_atlas_umap_by_compartment.pdf"
    fig.savefig(out_pdf2, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf2)

    # ---- summary table for traceability ----
    summary = (df.groupby(["compartment", "label"]).size()
                 .reset_index(name="n_cells")
                 .sort_values("n_cells", ascending=False))
    summary_path = args.out_dir / "flex_l2s_full_atlas_umap_summary.csv"
    summary.to_csv(summary_path, index=False)
    log.info("  wrote %s", summary_path)


if __name__ == "__main__":
    main()
