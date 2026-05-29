"""render_supp5b_label_umap.py — Fig 3 Supp 5b: per-compartment L2S label UMAP.

One panel per compartment (Epithelial / Immune / Stromal). Reads canonical
FLEX integration intermediate (umap.csv + obs.csv) joined with the Track C
L2S labels (publication/analysis/annotation/flex/outputs/flex_l2s_labels.csv).
Renders a rasterized scatter inside a vector PDF wrapper (Nature spec for
high-density UMAPs), with centroid text labels and force-directed repulsion.

Pattern reference:
- Cell coloring + centroid label placement: publication/analysis/annotation/
  flex/scripts/06_render_annotation.py:415-498 (render_label_umap)
- Framework conformance (load_aesthetics + load_paths): publication/figures/
  render/flex/render_flex_umap_l2s.py
- adjustText fallback (internal force-directed repel): same 06_render_annotation
  block 444-468.

No on-plot title, panel letter, or method caption — this is a publication
panel intended for Illustrator composition.

Outputs:
    supp5b_flex_label_umap_{epi,imm,str}.pdf
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
log = logging.getLogger("supp5b")


COMPARTMENT_FULL = {"epi": "Epithelial", "imm": "Immune", "str": "Stromal"}
ARTIFACT_COLOR = "#CCCCCC"


def _setup_aesthetics(project_root: Path) -> Tuple[dict, dict]:
    """Load load_aesthetics and apply matplotlib rcParams.

    Returns (aesthetics_config, dimensions_dict).
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
    """Tab20-family palette extended with HSV; matches 06_render_annotation:403-412."""
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
    iters: int = 4000, anchor_pull: float = 0.001, sep_frac: float = 0.25,
    repel_frac: float = 0.06,
) -> np.ndarray:
    """Force-directed repulsion of anchors. Same algorithm as 06_render_annotation:443-468."""
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
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
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
    cfg, dims = _setup_aesthetics(project_root)

    sys.path.insert(0, str(project_root / "publication" / "config"))
    from load_paths import load_paths, resolve_path  # type: ignore
    paths = load_paths(str(project_root / "publication" / "config"))

    compartment_full = COMPARTMENT_FULL[args.compartment]
    flex_root = Path(resolve_path("flex_integration_intermediate", paths))
    bundle = flex_root / compartment_full / "scvi_n100"
    umap_csv = bundle / "umap.csv"
    obs_csv = bundle / "obs.csv"
    l2s_csv = Path(resolve_path("flex_l2s_labels", paths))

    log.info("UMAP CSV  : %s", umap_csv)
    log.info("Obs CSV   : %s", obs_csv)
    log.info("L2S CSV   : %s", l2s_csv)
    for p in (umap_csv, obs_csv, l2s_csv):
        if not p.exists():
            sys.exit(f"FATAL: missing input: {p}")

    umap = pd.read_csv(umap_csv)
    if "UMAP_1" not in umap.columns:
        # Fallback for lower-cased schema
        cols = [c for c in umap.columns if c != "cell_id"]
        umap = umap.rename(columns={cols[0]: "UMAP_1", cols[1]: "UMAP_2"})
    log.info("UMAP rows: %d", len(umap))

    # flex_l2s_labels.csv uses lowercase shortnames {epi, imm, str} in the
    # compartment column (not the full Epithelial/Immune/Stromal). Match on the
    # CLI shortname.
    l2s = pd.read_csv(l2s_csv, usecols=["cell_id", "compartment", "l2s_label"])
    l2s = l2s[l2s["compartment"] == args.compartment]
    log.info("L2S rows for %s (compartment=%s): %d",
             compartment_full, args.compartment, len(l2s))

    df = umap.merge(l2s, on="cell_id", how="inner")
    log.info("Joined: %d cells", len(df))

    if args.test:
        df = df.sample(n=min(5000, len(df)), random_state=0).reset_index(drop=True)
        log.info("--test: subsampled to %d cells", len(df))

    label_str = df["l2s_label"].astype(str)
    is_artifact = (
        df["l2s_label"].isna()
        | (label_str.str.lower() == "nan")
        | label_str.str.startswith("ARTIFACT_")
    )
    df["is_artifact"] = is_artifact.values
    log.info("Artifact cells: %d", int(df["is_artifact"].sum()))

    # Order labels by descending cell count for stable color assignment.
    nonart = df[~df["is_artifact"]]
    label_counts = nonart.groupby("l2s_label").size().sort_values(ascending=False)
    label_order = label_counts.index.tolist()
    palette = categorical_palette(len(label_order))
    color_map = {lbl: palette[i] for i, lbl in enumerate(label_order)}

    # Per Noah review (2026-05-07): 5b panels rendered too light at print scale.
    # Override the aesthetics.yaml umap-panel pt_size (0.05) and alpha (0.3)
    # with print-readable values. Width bumped to ~60mm to match the v2 layout
    # cell (3 cells across the 183mm Nature double column with 5mm gutters).
    pt_size = 2.5
    alpha = 0.9
    width_in = 60.0 / 25.4   # 60mm in inches
    height_in = width_in     # square UMAP cell

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    ax.set_aspect("equal")
    ax.axis("off")

    art = df[df["is_artifact"]]
    if len(art):
        ax.scatter(
            art["UMAP_1"], art["UMAP_2"],
            c=ARTIFACT_COLOR, s=pt_size, alpha=alpha * 0.5,
            rasterized=True, linewidths=0,
        )

    for lbl in label_order:
        sub = nonart[nonart["l2s_label"] == lbl]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP_1"], sub["UMAP_2"],
            c=[color_map[lbl]], s=pt_size, alpha=alpha,
            rasterized=True, linewidths=0,
        )

    centroids = (
        nonart.groupby("l2s_label", as_index=False)
              .agg(UMAP_1=("UMAP_1", "median"), UMAP_2=("UMAP_2", "median"),
                   n_cells=("cell_id", "size"))
    )
    # Drop tiny labels (<50 cells) from text placement to avoid clutter; they
    # remain colored on the plot but unlabeled. Same threshold as 06_render_annotation:431.
    centroids = centroids[centroids["n_cells"] >= 50].reset_index(drop=True)

    x_range = float(df["UMAP_1"].max() - df["UMAP_1"].min())
    y_range = float(df["UMAP_2"].max() - df["UMAP_2"].min())
    anchors = centroids[["UMAP_1", "UMAP_2"]].to_numpy(dtype=float)

    use_adjust = False
    final_pos = anchors.copy()
    try:
        from adjustText import adjust_text  # type: ignore
        use_adjust = True
    except ImportError:
        log.info("adjustText unavailable; using internal repel")
        final_pos = repel_labels(anchors, x_range, y_range)

    text_artists = []
    for i, row in centroids.iterrows():
        ax_x, ax_y = float(row["UMAP_1"]), float(row["UMAP_2"])
        lx, ly = float(final_pos[i, 0]), float(final_pos[i, 1])
        if not use_adjust and ((lx != ax_x) or (ly != ax_y)):
            ax.plot([ax_x, lx], [ax_y, ly], color="grey", lw=0.5, alpha=0.7, zorder=9)
        t = ax.text(
            lx, ly, str(row["l2s_label"]),
            fontsize=5, ha="center", va="center", color="black", weight="bold",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="black",
                      lw=0.4, alpha=0.92),
            zorder=10,
        )
        text_artists.append(t)
    if use_adjust:
        adjust_text(
            text_artists, ax=ax,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.4, alpha=0.6),
            expand_text=(2.2, 2.4), expand_points=(1.6, 1.8),
            force_text=3.0, force_points=1.0,
            only_move={"text": "xy"}, iter_lim=1000,
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = args.out_dir / f"supp5b_flex_label_umap_{args.compartment}.pdf"
    fig.savefig(out_pdf, bbox_inches="tight",
                dpi=int(cfg.get("rendering", {}).get("dpi_umap_raster", 600)))
    plt.close(fig)
    log.info("Wrote %s", out_pdf)
    return 0


if __name__ == "__main__":
    sys.exit(main())
