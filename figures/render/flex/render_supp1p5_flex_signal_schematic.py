"""
render_supp1p5_flex_signal_schematic.py — Fig 3 Supp 1.5, Panel 1.5c.

Schematic illustrating the BigSur manifold-incoherence metric. Two side-by-
side cartoons of the same operation (averaging a cell's expression vector
with its k=5 neighbors), one for biology and one for noise:

    LEFT — Valid signal: a real cell surrounded by 5 transcriptionally
    similar neighbors. Their expression vectors all point roughly the same
    way; the smoothed average reinforces the cell's signal. Smooth-R² ≈
    raw-R² → Δ ≈ 0.

    RIGHT — Noise (debris/artifact): a cell whose nearest neighbors in the
    integrated embedding are unrelated (it sits where it does because of
    incoherent counts, not biology). Their expression vectors point in
    different directions; the smoothed average collapses. Smooth-R² → 0
    while raw-R² stays positive → Δ << 0.

Reviewer question (Act II of Supp 1.5): "How does Δ work mechanically?"

This is a stylized matplotlib drawing — no data. Pure compositional
illustration intended to be dropped into the supplemental layout.

Outputs:
    supp1p5_flex_signal_schematic.pdf  (no _data.csv — schematic)

Framework conformance: Option B per CP_supp7_python_framework_gap.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
    get_dimensions,
    get_matplotlib_theme,
    get_palette,
    load_aesthetics,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("supp1p5_1c")


N_NEIGHBORS = 5
CENTER_R = 0.18
NEIGHBOR_R = 0.13
ORBIT = 0.55
ARROW_LEN = 0.22
SEED = 7


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def neighbor_positions(k: int, orbit: float) -> np.ndarray:
    angles = np.linspace(0, 2 * np.pi, k, endpoint=False) + np.pi / 2
    return np.column_stack([orbit * np.cos(angles), orbit * np.sin(angles)])


def draw_cell(ax, x, y, r, color, edgecolor="white", lw=1.0):
    ax.add_patch(mpatches.Circle((x, y), r, facecolor=color, edgecolor=edgecolor,
                                 linewidth=lw, zorder=3))


def draw_arrow(ax, x, y, dx, dy, color, alpha=0.95):
    ax.annotate("", xy=(x + dx, y + dy), xytext=(x, y),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2,
                                shrinkA=0, shrinkB=0, alpha=alpha),
                zorder=4)


def render_panel(ax, signal: bool, palette: dict) -> None:
    rng = np.random.default_rng(SEED if signal else SEED + 1)

    if signal:
        center_color = palette.get("Epithelial", "#4A7BB7")
        neighbor_color = palette.get("Epithelial", "#4A7BB7")
        base_angle = np.pi / 3
        jitter = 0.18
    else:
        center_color = "#A0A0A0"
        neighbor_color = "#C8C8C8"
        base_angle = 0.0
        jitter = 0.0

    pos = neighbor_positions(N_NEIGHBORS, ORBIT)
    for i, (nx, ny) in enumerate(pos):
        if signal:
            angle = base_angle + rng.normal(0, jitter)
        else:
            angle = rng.uniform(0, 2 * np.pi)
        dx = ARROW_LEN * np.cos(angle)
        dy = ARROW_LEN * np.sin(angle)
        draw_cell(ax, nx, ny, NEIGHBOR_R, neighbor_color, lw=0.8)
        draw_arrow(ax, nx, ny, dx, dy, neighbor_color, alpha=0.95)

    if signal:
        center_angle = base_angle + rng.normal(0, jitter)
    else:
        center_angle = rng.uniform(0, 2 * np.pi)
    dx_c = ARROW_LEN * np.cos(center_angle)
    dy_c = ARROW_LEN * np.sin(center_angle)
    draw_cell(ax, 0, 0, CENTER_R, center_color, edgecolor="black", lw=1.4)
    draw_arrow(ax, 0, 0, dx_c, dy_c, "black", alpha=1.0)

    if signal:
        avg_dx = ARROW_LEN * np.cos(base_angle)
        avg_dy = ARROW_LEN * np.sin(base_angle)
        delta_text = r"$\Delta \approx 0$"
        verdict = "Coherent neighborhood"
    else:
        avg_dx = ARROW_LEN * 0.05
        avg_dy = ARROW_LEN * 0.05
        delta_text = r"$\Delta \ll 0$"
        verdict = "Incoherent neighborhood"

    ax.annotate("", xy=(0.0 + 1.5 * avg_dx, -0.95 + 1.5 * avg_dy),
                xytext=(0.0, -0.95),
                arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6,
                                shrinkA=0, shrinkB=0),
                zorder=4)
    ax.text(0, -1.18, "smoothed average", ha="center", va="top", fontsize=6,
            color="black", style="italic")
    ax.text(0, 0.95, verdict, ha="center", va="bottom", fontsize=7,
            color="black", weight="bold")
    ax.text(0.85, -0.95, delta_text, ha="left", va="center", fontsize=8,
            color="black", weight="bold")

    ax.set_xlim(-1.0, 1.0)
    ax.set_ylim(-1.4, 1.1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    aes = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(aes))
    palette = get_palette("compartment", aes)
    dims = get_dimensions("schematic", aes)

    fig, axes = plt.subplots(1, 2, figsize=(dims["width"], dims["height"]))
    render_panel(axes[0], signal=True, palette=palette)
    render_panel(axes[1], signal=False, palette=palette)

    fig.text(0.27, 0.93, "Valid signal (biology)", ha="center", va="bottom",
             fontsize=8, weight="bold")
    fig.text(0.78, 0.93, "Noise (artifact)", ha="center", va="bottom",
             fontsize=8, weight="bold")

    out_pdf = args.out_dir / "supp1p5_flex_signal_schematic.pdf"
    fig.savefig(out_pdf, format="pdf", dpi=1200, bbox_inches="tight")
    log.info("Wrote %s", out_pdf)


if __name__ == "__main__":
    main()
