"""
render_fig1b_l2_umap_v2.py — full-atlas L2 UMAP, compartment-grouped legend.

Replaces the earlier R render. Reads panel_substrate.csv (already joined),
renders with matplotlib using the prior 01_project_l2_onto_full_umap.py
visual conventions (KDE-based L1 dashed outlines, density-mode labels)
plus a compartment-grouped multi-legend layout.

Output: fig1_ihbca_umap_l2_numbered.pdf + .png + fig1b_label_key.csv
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("fig1b_v2")


# Distinct 23-color shared palette — reused per compartment.
# Max within-compartment label count is immune at 23, so 23 colors are required.
# Hand-tuned for maximum hue separation: alternating warm/cool, distinct lightness.
SHARED_PALETTE = [
    "#E6194B",  # 1  red
    "#3CB44B",  # 2  green
    "#4363D8",  # 3  blue
    "#F58231",  # 4  orange
    "#911EB4",  # 5  purple
    "#42D4F4",  # 6  cyan
    "#F032E6",  # 7  magenta
    "#BFEF45",  # 8  lime
    "#FABED4",  # 9  pink
    "#469990",  # 10 teal
    "#9A6324",  # 11 brown
    "#800000",  # 12 maroon
    "#AAFFC3",  # 13 mint
    "#808000",  # 14 olive
    "#000075",  # 15 navy
    "#FFD700",  # 16 gold
    "#A9A9A9",  # 17 dark grey
    "#FF6347",  # 18 tomato
    "#4169E1",  # 19 royal blue
    "#228B22",  # 20 forest green
    "#DA70D6",  # 21 orchid
    "#CD853F",  # 22 peru
    "#000000",  # 23 black
]

ARTIFACT_COLOR = "#CCCCCC"

COMPARTMENT_DISPLAY = {"imm": "Immune", "str": "Stromal", "epi": "Epithelial"}

# Per-L1.0-label nudges in DATA-UNIT increments (UMAP1/UMAP2 range ~30 each).
# Positive y = up, positive x = right. Match keys are substrings (case-insens).
L1_NUDGE = {
    # cumulative nudges from group median (data units; UMAP range ~30 each axis)
    "basal":     {"dy": -3.02},
    "lhs":       {"dy":  3.02},                    # up 10% (was 0)
    "lumhr":     {"dy":  3.02},
    "lasp":      {"dy":  6.03},
    "lumsec":    {"dy":  6.03},
    "b-":        {"dy":  0.60, "dx": -1.54},
    "b cell":    {"dy":  0.60, "dx": -1.54},
    "b lymph":   {"dy":  0.60, "dx": -1.54},
    "myeloid":   {"dy": -1.07, "dx":  2.05},       # down 1% (was -0.77)
    "t-":        {"dy":  0.91, "dx": -2.05},
    "t cell":    {"dy":  0.91, "dx": -2.05},
    "t lymph":   {"dy":  0.91, "dx": -2.05},
    "lymphatic": {"dy":  1.52},
    "vascular":  {"dy": -0.30},                    # down 1%
    "vasc":      {"dy": -0.30},
}

# Wrap L1 display strings if longer than this many chars (split at last space
# inside the limit). Used so multi-word L1 names like "Luminal Secretory" or
# "Vascular endothelial" wrap onto two lines.
L1_DISPLAY_MAX_CHARS = 12


def _wrap_l1_label(s: str) -> str:
    if len(s) <= L1_DISPLAY_MAX_CHARS:
        return s
    cut = s.rfind(" ", 0, L1_DISPLAY_MAX_CHARS + 1)
    if cut <= 0:
        return s
    return s[:cut] + "\n" + s[cut + 1:]


# Per-(compartment, label) L2 nudges in data units. Applied AFTER repel.
L2_NUDGE = {
    ("imm", "IFNg_T"):     {"dx": 0.90, "dy": -0.15},
    ("epi", "LASP-major"): {"dy":  1.51},              # up 5% (#1)
    ("epi", "LHS-major"):  {"dy": -1.51},              # down 5% (#2)
}


def _apply_l1_nudge(label: str, x: float, y: float) -> tuple[float, float]:
    key = label.lower()
    for needle, off in L1_NUDGE.items():
        if needle in key:
            return x + float(off.get("dx", 0.0)), y + float(off.get("dy", 0.0))
    return x, y
# Numbering order: epithelial (1..) → stromal → immune; within compartment by
# descending cell count.
COMPARTMENT_ORDER = ("epi", "str", "imm")


def _repel_labels(
    ax, movable_texts: list, static_texts: list,
    iters: int = 100, step: float = 0.05, min_sep_data_units: float = 0.4,
) -> None:
    """Iterative force-directed repel.

    movable_texts move to avoid each other AND static_texts. After each
    iteration, draw an anchor line from the moved label back to its
    original (centroid) position.
    """
    if not movable_texts:
        return
    origins = [(t.get_position(), t) for t in movable_texts]
    positions = np.array([list(t.get_position()) for t in movable_texts], dtype=float)
    static_positions = (
        np.array([list(t.get_position()) for t in static_texts], dtype=float)
        if static_texts else np.zeros((0, 2))
    )

    for _ in range(iters):
        moved = False
        # Movable vs movable
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                d = positions[j] - positions[i]
                dist = float(np.linalg.norm(d))
                if dist < min_sep_data_units and dist > 1e-9:
                    push = (min_sep_data_units - dist) * 0.5 * step / max(dist, 1e-9)
                    positions[i] -= d * push
                    positions[j] += d * push
                    moved = True
        # Movable vs static (movable yields)
        for i in range(len(positions)):
            for s in static_positions:
                d = s - positions[i]
                dist = float(np.linalg.norm(d))
                if dist < min_sep_data_units * 1.5 and dist > 1e-9:
                    push = (min_sep_data_units * 1.5 - dist) * step / max(dist, 1e-9)
                    positions[i] -= d * push
                    moved = True
        if not moved:
            break

    # Apply new positions; draw thin anchor lines back to original centroid
    for (orig, t), new in zip(origins, positions):
        t.set_position((float(new[0]), float(new[1])))
        d = float(np.linalg.norm(np.array(new) - np.array(list(orig))))
        if d > 0.15:
            ax.plot(
                [orig[0], float(new[0])], [orig[1], float(new[1])],
                color="grey", linewidth=0.3, alpha=0.5, zorder=11,
            )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--substrate", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def draw_l1_outlines(
    ax, umap_df: pd.DataFrame, l1_col: str = "compartment",
    subsample: int = 8000, density_frac: float = 0.98,
    grid_n: int = 200,
    line_color: str = "black", line_width: float = 0.8,
    line_alpha: float = 0.85, label_fontsize: int = 8,
    text_collector: list | None = None,
) -> None:
    """KDE-based dashed isocontour + density-mode label per L1 (compartment) group.

    Mirrors 01_project_l2_onto_full_umap.draw_l1_outlines verbatim except
    label_fontsize default raised from 7 → 8 for the larger composite figure.
    """
    labeled = umap_df[umap_df[l1_col].notna() & (umap_df[l1_col] != "nan")]
    rng = np.random.default_rng(42)

    for l1_type, group in labeled.groupby(l1_col):
        pts = pd.DataFrame(group)[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
        if len(pts) < 10:
            continue
        if len(pts) > subsample:
            idx = rng.choice(len(pts), subsample, replace=False)
            pts = pts[idx]
        try:
            kde = gaussian_kde(pts.T, bw_method="scott")
        except Exception as e:
            log.warning("KDE failed for %s: %s", l1_type, e)
            continue
        pad = 0.05 * (pts.max(axis=0) - pts.min(axis=0))
        xmin, ymin = pts.min(axis=0) - pad
        xmax, ymax = pts.max(axis=0) + pad
        xx, yy = np.meshgrid(np.linspace(xmin, xmax, grid_n), np.linspace(ymin, ymax, grid_n))
        zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(grid_n, grid_n)
        flat = np.sort(zz.ravel())[::-1]
        cumulative = np.cumsum(flat) / flat.sum()
        cut_idx = int(np.searchsorted(cumulative, density_frac))
        threshold = float(flat[min(cut_idx, len(flat) - 1)])
        ax.contour(
            xx, yy, zz, levels=[threshold],
            colors=line_color, linewidths=line_width,
            alpha=line_alpha, linestyles="dashed", zorder=10,
        )
        # L1.0 label: median of group + per-label nudge (deterministic).
        cx = float(np.median(pts[:, 0]))
        cy = float(np.median(pts[:, 1]))
        display = COMPARTMENT_DISPLAY.get(str(l1_type), str(l1_type))
        tx, ty = _apply_l1_nudge(display, cx, cy)
        wrapped = _wrap_l1_label(display)
        t = ax.text(
            tx, ty, wrapped,
            fontsize=label_fontsize, fontweight="bold",
            ha="center", va="center", zorder=11,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="black", alpha=0.88, linewidth=0.5),
        )
        if text_collector is not None:
            text_collector.append(t)


def build_color_map(label_order_per_compartment: dict[str, list[str]]) -> dict[tuple[str, str], str]:
    """Per-(compartment, label) hex from SHARED_PALETTE, restarted per compartment.

    label_order_per_compartment: {compartment: [labels in display order]}.
    Color N within each compartment uses SHARED_PALETTE[N]. Same color appears
    in multiple compartments by design — cluster position + L1.0 outline +
    numbered labels disambiguate.
    """
    color_map: dict[tuple[str, str], str] = {}
    for compartment, labels in label_order_per_compartment.items():
        for i, lbl in enumerate(labels):
            color_map[(compartment, lbl)] = SHARED_PALETTE[i % len(SHARED_PALETTE)]
    return color_map


def main() -> int:
    args = parse_args()

    log.info("loading substrate %s", args.substrate)
    dt = pd.read_csv(args.substrate)
    log.info("  rows: %d", len(dt))

    # Numbering order: compartment block (epi → str → imm) × within-compartment by
    # DESCENDING cell count. So #1 = largest epithelial label, ..., #N = smallest
    # immune label.
    label_order_per_comp: dict[str, list[str]] = {}
    for compartment in COMPARTMENT_ORDER:
        sub = dt[(dt["compartment"] == compartment) & (~dt["is_artifact"])]
        if sub.empty:
            label_order_per_comp[compartment] = []
            continue
        comp_counts = sub.groupby("label").size().sort_values(ascending=False)
        label_order_per_comp[compartment] = comp_counts.index.tolist()

    label_order: list[tuple[str, str]] = []
    for compartment in COMPARTMENT_ORDER:
        label_order.extend((compartment, lbl) for lbl in label_order_per_comp[compartment])
    label_to_num = {lbl: i + 1 for i, lbl in enumerate(label_order)}

    color_map = build_color_map(label_order_per_comp)
    log.info("  unique non-artifact (compartment, label) pairs: %d", len(color_map))

    # Per-row plot color
    dt = dt.copy()
    dt["plot_color"] = ARTIFACT_COLOR
    for (compartment, lbl), hex_color in color_map.items():
        mask = (~dt["is_artifact"]) & (dt["compartment"] == compartment) & (dt["label"] == lbl)
        dt.loc[mask, "plot_color"] = hex_color

    # Centroids per (compartment, label) — non-artifact only
    centroids = (
        dt[~dt["is_artifact"]]
        .groupby(["compartment", "label"], as_index=False)
        .agg(UMAP1=("UMAP1", "median"), UMAP2=("UMAP2", "median"), n_cells=("cell_id", "size"))
    )
    centroids["label_num"] = centroids.apply(lambda r: label_to_num[(r["compartment"], r["label"])], axis=1)
    centroids = centroids.sort_values("label_num").reset_index(drop=True)

    # Figure size — generous for legends on the right
    fig, ax = plt.subplots(figsize=(11.5, 8.0))
    ax.set_aspect("equal")
    ax.axis("off")

    # Filter artifacts entirely from final UMAP (per Pascal feedback 2026-04-30)
    nonart = dt[~dt["is_artifact"]]
    for compartment, lbl in label_order:
        sub = nonart[(nonart["compartment"] == compartment) & (nonart["label"] == lbl)]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP1"], sub["UMAP2"],
            c=color_map[(compartment, lbl)], s=0.08, alpha=0.6,
            rasterized=True, linewidths=0,
        )

    # Numbered centroids — collision-aware via adjustText (ggrepel equivalent).
    # Falls back to fixed placement if adjustText unavailable.
    number_texts = []
    for _, r in centroids.iterrows():
        t = ax.text(
            r["UMAP1"], r["UMAP2"], str(int(r["label_num"])),
            fontsize=6, fontweight="bold", ha="center", va="center", zorder=12,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="black", alpha=0.92, linewidth=0.4),
        )
        number_texts.append(t)

    # L1.0 outlines (level1_annotation from cell_metadata, joined into substrate)
    l1_texts: list = []
    if "l1_annotation" in dt.columns:
        log.info("drawing L1.0 outlines on l1_annotation (%d unique)", int(dt["l1_annotation"].nunique()))
        draw_l1_outlines(ax, dt, l1_col="l1_annotation", label_fontsize=8, text_collector=l1_texts)
    else:
        log.warning("l1_annotation column missing from substrate; outlines skipped")

    # L1 labels stay at their nudged median positions (no repel — explicit).
    # L2 numbers get TARGETED stronger repel within dense regions (myeloid +
    # T-cell). Identify "dense region" L2 numbers as the immune labels that
    # sit within 4 data-units of any T-cell or Myeloid centroid.
    #
    # Build immune-group anchors from the substrate's compartment + l1_annotation
    immune_anchors: list[tuple[float, float]] = []
    if "l1_annotation" in dt.columns:
        for _, group in dt[dt["compartment"] == "imm"].groupby("l1_annotation"):
            key = str(group["l1_annotation"].iloc[0]).lower()
            if "t" in key or "myeloid" in key or "macro" in key:
                immune_anchors.append((float(np.median(group["UMAP1"])), float(np.median(group["UMAP2"]))))

    dense_radius = 4.0  # data units
    dense_idx, sparse_idx = [], []
    for i, t in enumerate(number_texts):
        x, y = t.get_position()
        in_dense = any(((x - ax_)**2 + (y - ay_)**2) ** 0.5 < dense_radius for ax_, ay_ in immune_anchors)
        (dense_idx if in_dense else sparse_idx).append(i)

    log.info("L2 repel: %d dense (myeloid/T-cell) + %d sparse", len(dense_idx), len(sparse_idx))
    if dense_idx:
        _repel_labels(
            ax, [number_texts[i] for i in dense_idx], static_texts=[],
            iters=300, step=0.06, min_sep_data_units=0.85,
        )
    if sparse_idx:
        _repel_labels(
            ax, [number_texts[i] for i in sparse_idx], static_texts=[],
            iters=120, step=0.04, min_sep_data_units=0.45,
        )

    # Apply per-(compartment, label) L2 nudges AFTER repel.
    # Maps label_num back to (compartment, label) via centroids table.
    centroids_indexed = centroids.set_index("label_num")
    for i, t in enumerate(number_texts):
        try:
            num = int(t.get_text())
        except ValueError:
            continue
        if num not in centroids_indexed.index:
            continue
        row = centroids_indexed.loc[num]
        key = (str(row["compartment"]), str(row["label"]))
        nudge = L2_NUDGE.get(key)
        if nudge is None:
            continue
        x, y = t.get_position()
        t.set_position((float(x) + float(nudge.get("dx", 0.0)),
                        float(y) + float(nudge.get("dy", 0.0))))
        log.info("nudged L2 #%d %s: dx=%.2f dy=%.2f", num, key, nudge.get("dx", 0.0), nudge.get("dy", 0.0))

    # Emit final positions for hand-tuning
    import json as _json
    pos_dump = {
        "l1": [{"text": t.get_text(), "x": float(t.get_position()[0]), "y": float(t.get_position()[1])} for t in l1_texts],
        "l2_numbered": [{"text": t.get_text(), "x": float(t.get_position()[0]), "y": float(t.get_position()[1])} for t in number_texts],
    }
    pos_path = args.out_dir / "fig1b_label_positions.json"
    pos_path.write_text(_json.dumps(pos_dump, indent=2))
    log.info("wrote label positions to %s", pos_path)

    # Two-row compartment-grouped legend layout:
    #   Top row:    Stromal (left)   |  Epithelial (right)
    #   Bottom row: Immune (split into 2 columns to handle 23 labels)
    legend_artists = []
    # Final pre-promotion: epi -1% (0.675 → 0.665), str at 0.825,
    # imm +2% up (0.79 → 0.81).
    top_layout = [("epi", 0.665, 1.00, 1), ("str", 0.825, 1.00, 1)]
    bottom_layout = [("imm", 0.695, 0.81, 2)]

    def _make_legend(compartment: str, x: float, y: float, ncol: int):
        comp_centroids = centroids[centroids["compartment"] == compartment].sort_values("label_num")
        if comp_centroids.empty:
            return None
        handles = [
            mpatches.Patch(
                color=color_map[(compartment, r["label"])],
                label=f"{int(r['label_num'])}. {r['label']}  ({int(r['n_cells']):,})",
            )
            for _, r in comp_centroids.iterrows()
        ]
        leg = ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(x, y),
            fontsize=4, frameon=False, markerscale=1.6, handlelength=0.6, handletextpad=0.3,
            title=f"{COMPARTMENT_DISPLAY[compartment]} (n={int(comp_centroids['n_cells'].sum()):,})",
            title_fontsize=5.5,
            borderaxespad=0.0,
            ncol=ncol, columnspacing=0.6,
        )
        leg.get_title().set_fontweight("bold")
        return leg

    for compartment, x, y, ncol in top_layout + bottom_layout:
        leg = _make_legend(compartment, x, y, ncol)
        if leg is not None:
            legend_artists.append(leg)
    # Re-add prior legends so they all render (matplotlib only keeps the last by default)
    for leg in legend_artists[:-1]:
        ax.add_artist(leg)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / "fig1_ihbca_umap_l2_numbered.pdf"
    png_path = args.out_dir / "fig1_ihbca_umap_l2_numbered.png"
    # DPI for rasterized scatter from aesthetics.yaml (rendering.dpi_umap_raster).
    import yaml as _yaml  # local import to avoid touching module-level imports
    _aes_path = Path(__file__).resolve().parents[3] / "config" / "aesthetics.yaml"
    with open(_aes_path) as _fh:
        _raster_dpi = int(_yaml.safe_load(_fh).get("rendering", {}).get("dpi_umap_raster", 600))
    fig.savefig(pdf_path, bbox_inches="tight", dpi=_raster_dpi)
    fig.savefig(png_path, bbox_inches="tight", dpi=_raster_dpi)
    plt.close(fig)
    log.info("wrote %s + PNG", pdf_path)

    key_path = args.out_dir / "fig1b_label_key.csv"
    centroids[["label_num", "compartment", "label", "n_cells", "UMAP1", "UMAP2"]].rename(
        columns={"UMAP1": "UMAP1_centroid", "UMAP2": "UMAP2_centroid"}
    ).to_csv(key_path, index=False)
    log.info("wrote label-num key %s", key_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
