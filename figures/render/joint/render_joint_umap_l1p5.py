"""
render_joint_umap_l1p5.py — Fig 3c: joint FLEX+Xenium UMAP labeled by L1.5.

Reads joint_umap.csv and joint_l1p5.csv. Inner-joins on cell_id; renders
with the fig1b L2-numbered convention adapted for L1.5 vocabulary:
SHARED_PALETTE per L1/L0.5 lineage, KDE-based lineage outlines, force-directed
repel on L1.5 number labels, lineage-grouped legend.

Lineage = `l1p5_lineage` (~13 non-artifact: Basal / Luminal / Fibroblast /
Pericyte / Adipocyte / Endothelial / Myeloid + DC/Granulocyte/Mast /
Lymphoid_T / Lymphoid_B / Lymphoid_NK). The compartment (L0) layer is
intentionally NOT used — joint integration figure operates at L1/L0.5
granularity to mirror the panel-resolvable Xenium tier.

Substrate respects joint_l1p5.csv `is_artifact` flag (cells dropped from
biology-facing analysis are still drawn under, but in artifact gray).

Outputs:
    fig3c_joint_umap_l1p5.pdf + .png + fig3c_joint_label_key.csv
    fig3c_joint_label_positions.json (for hand-tuning)
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
import pandas as pd
from scipy.stats import gaussian_kde

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("fig3c")


SHARED_PALETTE = [
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4",
    "#F032E6", "#BFEF45", "#FABED4", "#469990", "#9A6324", "#800000",
    "#AAFFC3", "#808000", "#000075", "#FFD700", "#A9A9A9", "#FF6347",
    "#4169E1", "#228B22", "#DA70D6", "#CD853F", "#000000",
]


def darken(hex_color: str, amount: float = 0.10) -> str:
    """Mix `hex_color` toward black by `amount` (0.0=no change, 1.0=black)."""
    r = int(hex_color[1:3], 16) * (1.0 - amount)
    g = int(hex_color[3:5], 16) * (1.0 - amount)
    b = int(hex_color[5:7], 16) * (1.0 - amount)
    return f"#{int(r):02X}{int(g):02X}{int(b):02X}"


# Globally darken SHARED_PALETTE by 1 step (~10%) for visual consistency
# with fig3b and to lift light entries (lime, pink, mint) off white background.
SHARED_PALETTE = [darken(c, 0.12) for c in SHARED_PALETTE]

ARTIFACT_COLOR = "#CCCCCC"

# L1 / L0.5 lineage order (general lineage types — the layer at which
# iHBCA L1 and Xenium L0.5 align). Joint figure uses this for dashed outlines
# instead of L0 (compartment, too coarse) or L1.5 (cell types, too fine).
# Numbering of L1.5 entries is CONTINUOUS across L0 (no per-compartment restart).
L1_LINEAGE_ORDER = (
    "Basal", "Luminal",                                    # epithelial
    "Fibroblast", "Pericyte", "Adipocyte", "Endothelial",  # stromal
    "Myeloid", "Lymphoid_T", "Lymphoid_B", "Lymphoid_NK",  # immune
)

# Collapse L1.5-level Myeloid sub-divisions in `l1p5_lineage` back into a
# single Myeloid L1 bucket. Other lineage values pass through unchanged.
LINEAGE_COLLAPSE = {
    "Myeloid": "Myeloid",
    "Myeloid_DC": "Myeloid",
    "Myeloid_Granulocyte": "Myeloid",
    "Myeloid_Mast": "Myeloid",
}

# Display names for L1 lineage outlines (drawn on the UMAP plot).
LINEAGE_DISPLAY = {
    "Basal": "Basal", "Luminal": "Luminal",
    "Fibroblast": "Fibroblast", "Pericyte": "Pericyte",
    "Adipocyte": "Adipocyte", "Endothelial": "Endothelial",
    "Myeloid": "Myeloid",
    "Lymphoid_T": "T cell", "Lymphoid_B": "B / Plasma", "Lymphoid_NK": "NK",
}

# L0 compartment order (legend grouping for 3c). Joint l1p5.csv uses
# capitalized strings (Epithelial/Stromal/Immune) — distinct from FLEX
# 3-letter lowercase (epi/str/imm).
COMPARTMENT_ORDER = ("Epithelial", "Stromal", "Immune")
COMPARTMENT_DISPLAY = {"Epithelial": "Epithelial", "Stromal": "Stromal", "Immune": "Immune"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--joint-umap", required=True, type=Path)
    p.add_argument("--joint-l1p5", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument(
        "--split-by-platform", action="store_true",
        help=(
            "Render a 1x2 side-by-side variant faceted by platform (FLEX | Xenium). "
            "Color map, numbered labels, and lineage convention stay identical across "
            "panels so the same legend applies to both. Output filename suffix: "
            "_split_by_platform."
        ),
    )
    return p.parse_args()


def build_substrate(umap_path: Path, l1p5_path: Path) -> pd.DataFrame:
    log.info("loading joint umap: %s", umap_path)
    umap = pd.read_csv(umap_path)
    log.info("  rows: %d", len(umap))

    log.info("loading joint l1p5: %s", l1p5_path)
    cols = ["cell_id", "platform", "compartment", "l1p5_label", "l1p5_lineage", "is_artifact"]
    l1p5 = pd.read_csv(l1p5_path, usecols=cols)
    log.info("  rows: %d", len(l1p5))

    sub = umap.merge(l1p5, on="cell_id", how="inner")
    log.info("substrate rows: %d", len(sub))
    sub["is_artifact"] = sub["is_artifact"].fillna(False).astype(bool)
    log.info("  artifact rows: %d", int(sub["is_artifact"].sum()))

    # Collapse L1.5-level Myeloid sub-lineages back to L1 (Myeloid_DC →
    # Myeloid, etc.). Other lineage values pass through unchanged.
    sub["l1_lineage"] = sub["l1p5_lineage"].map(lambda v: LINEAGE_COLLAPSE.get(v, v))
    return sub


def repel_labels(ax, texts: list, iters: int = 200, step: float = 0.05, min_sep: float = 0.6) -> None:
    if not texts:
        return
    origins = [(t.get_position(), t) for t in texts]
    pos = np.array([list(t.get_position()) for t in texts], dtype=float)
    for _ in range(iters):
        moved = False
        for i in range(len(pos)):
            for j in range(i + 1, len(pos)):
                d = pos[j] - pos[i]
                dist = float(np.linalg.norm(d))
                if dist < min_sep and dist > 1e-9:
                    push = (min_sep - dist) * 0.5 * step / max(dist, 1e-9)
                    pos[i] -= d * push
                    pos[j] += d * push
                    moved = True
        if not moved:
            break
    for (orig, t), new in zip(origins, pos):
        t.set_position((float(new[0]), float(new[1])))
        d = float(np.linalg.norm(np.array(new) - np.array(list(orig))))
        if d > 0.15:
            ax.plot(
                [orig[0], float(new[0])], [orig[1], float(new[1])],
                color="grey", linewidth=0.3, alpha=0.5, zorder=11,
            )


def draw_lineage_outlines(
    ax, dt: pd.DataFrame, density_frac: float = 0.95, grid_n: int = 200,
    label_fontsize: int = 7,
) -> list:
    """KDE outline + L1 lineage label. Operates on the collapsed `l1_lineage`
    column (post-Myeloid-collapse), so 10 outlines instead of 13."""
    rng = np.random.default_rng(42)
    text_artists: list = []
    labeled = dt[(~dt["is_artifact"]) & dt["l1_lineage"].notna() & (dt["l1_lineage"] != "nan")]
    for lineage, group in labeled.groupby("l1_lineage"):
        if str(lineage) not in L1_LINEAGE_ORDER:
            continue
        pts = group[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
        if len(pts) < 10:
            continue
        if len(pts) > 8000:
            idx = rng.choice(len(pts), 8000, replace=False)
            pts = pts[idx]
        try:
            kde = gaussian_kde(pts.T, bw_method="scott")
        except Exception as e:
            log.warning("KDE failed for %s: %s", lineage, e)
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
            colors="black", linewidths=0.7, alpha=0.7, linestyles="dashed", zorder=10,
        )
        cx = float(np.median(pts[:, 0])); cy = float(np.median(pts[:, 1]))
        display = LINEAGE_DISPLAY.get(str(lineage), str(lineage))
        t = ax.text(
            cx, cy, display,
            fontsize=label_fontsize, fontweight="bold",
            ha="center", va="center", zorder=11,
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="black", alpha=0.85, linewidth=0.4),
        )
        text_artists.append(t)
    return text_artists


def _render_axis(
    ax,
    dt: pd.DataFrame,
    label_order: list,
    label_to_num: dict,
    color_map: dict,
    centroids: pd.DataFrame,
    label_dominant_compartment: dict,
    min_label_n: int = 0,
) -> list:
    """Scatter + numbered L1.5 labels onto a single axis.

    Used both for the combined render (full dt) and per-platform-faceted renders
    (platform-filtered dt). Artifact cells drawn in artifact gray; non-artifact
    cells colored by (dominant_compartment, l1p5_label) lookup. Numbered labels
    placed at the centroid of cells present in this view; labels with fewer than
    ``min_label_n`` cells in this view are skipped (set to 0 to draw all).

    Returns the list of label text artists for downstream repel.
    """
    ax.set_aspect("equal")
    ax.axis("off")

    art = dt[dt["is_artifact"]]
    if len(art):
        ax.scatter(
            art["UMAP1"], art["UMAP2"],
            c=ARTIFACT_COLOR, s=0.05, alpha=0.20,
            rasterized=True, linewidths=0,
        )

    nonart_scatter = dt[~dt["is_artifact"]].copy()
    nonart_scatter["dominant_compartment"] = nonart_scatter["l1p5_label"].map(
        label_dominant_compartment
    )
    label_counts_here = (
        nonart_scatter.groupby(["dominant_compartment", "l1p5_label"]).size().to_dict()
    )
    for comp, lbl in label_order:
        sub = nonart_scatter[
            (nonart_scatter["dominant_compartment"] == comp)
            & (nonart_scatter["l1p5_label"] == lbl)
        ]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP1"], sub["UMAP2"],
            c=color_map[(comp, lbl)], s=0.10, alpha=0.75,
            rasterized=True, linewidths=0,
        )

    number_texts: list = []
    for _, r in centroids.iterrows():
        key = (r["compartment"], r["l1p5_label"])
        n_here = label_counts_here.get(key, 0)
        if n_here < min_label_n:
            continue
        t = ax.text(
            r["UMAP1"], r["UMAP2"], str(int(r["label_num"])),
            fontsize=6, fontweight="bold", ha="center", va="center", zorder=12,
            bbox=dict(
                boxstyle="round,pad=0.18", facecolor="white", edgecolor="black",
                alpha=0.92, linewidth=0.4,
            ),
        )
        number_texts.append(t)
    return number_texts


def main() -> int:
    args = parse_args()

    dt = build_substrate(args.joint_umap, args.joint_l1p5)

    # Assign each L1.5 cell-type label its DOMINANT L0 compartment (mode of
    # `compartment` for non-artifact cells with that label). Drops cross-
    # compartment noise (e.g. cells assigned Fibroblast L1.5 but landing in
    # an immune compartment cluster) and yields one (compartment, l1p5_label)
    # row per label.
    nonart = dt[~dt["is_artifact"]]
    label_dominant_compartment = (
        nonart.groupby("l1p5_label")["compartment"]
        .agg(lambda s: s.value_counts().idxmax())
        .to_dict()
    )

    # Group + number L1.5 labels by their dominant L0 compartment. Numbering is
    # CONTINUOUS across compartments (no per-section restart).
    label_order_per_comp: dict[str, list[str]] = {c: [] for c in COMPARTMENT_ORDER}
    counts_per_label = nonart.groupby("l1p5_label").size().to_dict()
    for lbl, comp in label_dominant_compartment.items():
        if comp in COMPARTMENT_ORDER:
            label_order_per_comp[comp].append(lbl)
    for comp in COMPARTMENT_ORDER:
        label_order_per_comp[comp].sort(key=lambda l: -counts_per_label.get(l, 0))

    label_order: list[tuple[str, str]] = []
    for comp in COMPARTMENT_ORDER:
        label_order.extend((comp, lbl) for lbl in label_order_per_comp[comp])
    label_to_num = {pair: i + 1 for i, pair in enumerate(label_order)}

    # Per-compartment palette offset: shift SHARED_PALETTE start index by
    # compartment ordinal so colors rotate slightly across compartments and
    # don't collide for visually adjacent labels of different compartments.
    PALETTE_OFFSETS = {"Epithelial": 0, "Stromal": 1, "Immune": 2}
    color_map: dict[tuple[str, str], str] = {}
    for comp, labels in label_order_per_comp.items():
        offset = PALETTE_OFFSETS.get(comp, 0)
        for i, lbl in enumerate(labels):
            color_map[(comp, lbl)] = SHARED_PALETTE[(i + offset) % len(SHARED_PALETTE)]
    log.info("unique non-artifact L1.5 labels: %d (across %d L0 compartments)",
             sum(len(v) for v in label_order_per_comp.values()),
             sum(1 for v in label_order_per_comp.values() if v))

    # Centroids computed per L1.5 label (cells assigned their dominant L0
    # compartment); cross-compartment noise rows collapsed. Platform-stratified
    # cell counts (n_flex, n_xenium) added for legend display.
    nonart = nonart.copy()
    nonart["dominant_compartment"] = nonart["l1p5_label"].map(label_dominant_compartment)
    nonart = nonart[nonart["dominant_compartment"].isin(COMPARTMENT_ORDER)]
    centroids = (
        nonart.groupby(["dominant_compartment", "l1p5_label"], as_index=False)
        .agg(
            UMAP1=("UMAP1", "median"),
            UMAP2=("UMAP2", "median"),
            n_cells=("cell_id", "size"),
            n_flex=("platform", lambda s: int((s == "flex").sum())),
            n_xenium=("platform", lambda s: int((s == "xenium").sum())),
        )
        .rename(columns={"dominant_compartment": "compartment"})
    )
    centroids["label_num"] = centroids.apply(
        lambda r: label_to_num[(r["compartment"], r["l1p5_label"])], axis=1
    )
    centroids = centroids.sort_values("label_num").reset_index(drop=True)

    # Two render modes:
    #   - combined (default): single axis over full joint substrate
    #   - split_by_platform: 1x2 axes, FLEX | Xenium, same color/numbering
    if args.split_by_platform:
        # 3-column gridspec: FLEX | Xenium | legend.
        #   - Tiny wspace between the two UMAPs (move them together).
        #   - Spacer + dedicated legend column on the right (margin between
        #     UMAPs and legend comes from gridspec wspace + the legend axis's
        #     own left padding).
        fig = plt.figure(figsize=(22.0, 8.0))
        gs = fig.add_gridspec(
            1, 3, width_ratios=[1.0, 1.0, 0.32],
            wspace=0.04, left=0.03, right=0.99, top=0.95, bottom=0.05,
        )
        ax_flex = fig.add_subplot(gs[0, 0])
        ax_xen = fig.add_subplot(gs[0, 1])
        ax_leg = fig.add_subplot(gs[0, 2])
        ax_leg.axis("off")
        axes = [ax_flex, ax_xen]
        ax = ax_leg  # legends are drawn on the dedicated legend axis
        platform_axes = [("flex", ax_flex, "FLEX"), ("xenium", ax_xen, "Xenium")]
        # In split mode, only number L1.5 labels that have at least a minimal
        # cell complement in that platform's view (avoids a forest of numbered
        # bubbles in regions with near-zero same-platform cells).
        MIN_LABEL_N_SPLIT = 50
        all_texts: list = []
        for pf, pf_ax, pf_title in platform_axes:
            dt_pf = dt[dt["platform"] == pf]
            texts = _render_axis(
                pf_ax, dt_pf, label_order, label_to_num, color_map, centroids,
                label_dominant_compartment, min_label_n=MIN_LABEL_N_SPLIT,
            )
            repel_labels(pf_ax, texts, iters=200, step=0.05, min_sep=0.6)
            pf_ax.set_title(
                f"{pf_title} (n={len(dt_pf):,})",
                fontsize=10, fontweight="bold", pad=6,
            )
            all_texts.extend(texts)
        number_texts = all_texts
    else:
        fig, ax = plt.subplots(figsize=(11.5, 8.0))
        number_texts = _render_axis(
            ax, dt, label_order, label_to_num, color_map, centroids,
            label_dominant_compartment, min_label_n=0,
        )
        # Dashed L1 outlines intentionally omitted — labels are already at L1
        # granularity, so outlines would be redundant with the colored labels.
        repel_labels(ax, number_texts, iters=200, step=0.05, min_sep=0.6)

    import json as _json
    pos_dump = {
        "l1p5_numbered": [
            {"text": t.get_text(), "x": float(t.get_position()[0]), "y": float(t.get_position()[1])}
            for t in number_texts
        ],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "fig3c_joint_label_positions.json").write_text(_json.dumps(pos_dump, indent=2))

    # Compartment-stacked legend on the dedicated right-side legend axis.
    # In split mode, the three compartment legends stack vertically (single
    # column each); in combined mode, the prior 3-block top/bottom layout is
    # preserved on the main UMAP axis.
    legend_artists: list = []

    def _make_legend(target_ax, comp: str, x: float, y: float, ncol: int):
        comp_centroids = centroids[centroids["compartment"] == comp].sort_values("label_num")
        if comp_centroids.empty:
            return None
        handles = [
            mpatches.Patch(
                color=color_map[(comp, r["l1p5_label"])],
                label=f"{int(r['label_num'])}. {r['l1p5_label']}  ({int(r['n_flex']):,} | {int(r['n_xenium']):,})",
            )
            for _, r in comp_centroids.iterrows()
        ]
        n_flex_total = int(comp_centroids["n_flex"].sum())
        n_xen_total = int(comp_centroids["n_xenium"].sum())
        leg = target_ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(x, y),
            fontsize=4.8, frameon=False, markerscale=1.6, handlelength=0.6, handletextpad=0.3,
            title=(r"$\bf{" + COMPARTMENT_DISPLAY[comp] + r"}$" +
                   f"\n(FLEX {n_flex_total:,} | Xenium {n_xen_total:,})"),
            title_fontsize=5.8,
            borderaxespad=0.0,
            ncol=ncol, columnspacing=0.6,
            labelspacing=0.30,
        )
        # Center each line of the multiline title relative to the title block,
        # AND center the title block within the legend box.
        leg.get_title().set_multialignment("center")
        leg._legend_box.align = "center"
        return leg

    if args.split_by_platform:
        # Vertical stack on the dedicated legend axis. Y positions distribute
        # Epithelial → Stromal → Immune top→bottom; small left padding (x=0.10)
        # creates the margin between UMAP and legend that Noah requested.
        stack_layout = [
            ("Epithelial", 0.10, 0.98, 1),
            ("Stromal",    0.10, 0.66, 1),
            ("Immune",     0.10, 0.42, 1),
        ]
        for comp, x, y, ncol in stack_layout:
            leg = _make_legend(ax_leg, comp, x, y, ncol)
            if leg is not None:
                legend_artists.append(leg)
        for leg in legend_artists[:-1]:
            ax_leg.add_artist(leg)
    else:
        # Combined-mode layout unchanged: 3-block top/bottom on main axis.
        top_layout = [("Epithelial", 0.490, 0.975, 1), ("Stromal", 0.725, 0.975, 1)]
        bottom_layout = [("Immune", 0.595, 0.81, 2)]
        for comp, x, y, ncol in top_layout + bottom_layout:
            leg = _make_legend(ax, comp, x, y, ncol)
            if leg is not None:
                legend_artists.append(leg)
        for leg in legend_artists[:-1]:
            ax.add_artist(leg)

    suffix = "_split_by_platform" if args.split_by_platform else ""
    pdf_path = args.out_dir / f"fig3c_joint_umap_l1p5{suffix}.pdf"
    png_path = args.out_dir / f"fig3c_joint_umap_l1p5{suffix}.png"
    # DPI for rasterized scatter from aesthetics.yaml (rendering.dpi_umap_raster).
    import yaml as _yaml
    _aes_path = Path(__file__).resolve().parents[3] / "config" / "aesthetics.yaml"
    with open(_aes_path) as _fh:
        _raster_dpi = int(_yaml.safe_load(_fh).get("rendering", {}).get("dpi_umap_raster", 600))
    fig.savefig(pdf_path, bbox_inches="tight", dpi=_raster_dpi)
    fig.savefig(png_path, bbox_inches="tight", dpi=_raster_dpi)
    plt.close(fig)
    log.info("wrote %s + PNG", pdf_path)

    key_path = args.out_dir / "fig3c_joint_label_key.csv"
    centroids[["label_num", "compartment", "l1p5_label", "n_cells", "UMAP1", "UMAP2"]].rename(
        columns={"UMAP1": "UMAP1_centroid", "UMAP2": "UMAP2_centroid"}
    ).to_csv(key_path, index=False)
    log.info("wrote label-num key %s", key_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
