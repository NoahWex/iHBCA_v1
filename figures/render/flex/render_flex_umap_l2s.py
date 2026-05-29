"""
render_flex_umap_l2s.py — Fig 3b: FLEX UMAP labeled by L2S.

Reads canonical FLEX bundle (umap.csv + obs.csv) and Track C L2S labels.
Builds substrate at runtime (CSVs are small; no h5ad needed). Renders with
the fig1b L2-numbered convention: SHARED_PALETTE per compartment, KDE-based
L1/compartment dashed outlines, force-directed repel on L2S number labels,
compartment-grouped multi-legend.

Outputs:
    fig3b_flex_umap_l2s.pdf + .png + fig3b_flex_label_key.csv
    fig3b_flex_label_positions.json (for hand-tuning)
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
log = logging.getLogger("fig3b")


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


# Globally darken SHARED_PALETTE by 1 step (~10%) so light entries (lime,
# pink, mint, gold) read more clearly against the white panel background.
SHARED_PALETTE = [darken(c, 0.12) for c in SHARED_PALETTE]

ARTIFACT_COLOR = "#CCCCCC"

COMPARTMENT_DISPLAY = {"epi": "Epithelial", "str": "Stromal", "imm": "Immune"}
COMPARTMENT_ORDER = ("epi", "str", "imm")

# V1-harmonized L1 display labels (cross-referenced from iHBCA V1 panel
# substrate l1_annotation: Basal-myoepithelial, Luminal adaptive secretory
# precursor, Luminal hormone sensing, Fibroblast, Perivascular, Vascular
# endothelial, Lymphatic endothelial, Myeloid, T-lymphocyte, B-lymphocyte).
# FLEX-specific additions: Adipocyte (no V1 L1 equivalent at full scale) and
# Plasma (visually distinct subcluster of Lymphoid_B in spatial atlas).
# NK folded into T-lymphocyte per V1 L1 nomenclature.
V1_L1_ORDER = (
    "Basal-myoepithelial",
    "Luminal adaptive secretory precursor",
    "Luminal hormone sensing",
    "Fibroblast",
    "Perivascular",
    "Adipocyte",
    "Vascular endothelial",
    "Lymphatic endothelial",
    "Myeloid",
    "T-lymphocyte",
    "B-lymphocyte",
)

# On-plot wrapped display strings for long L1 labels (fonts ≤7pt; without
# wrapping, "Luminal adaptive secretory precursor" overruns the cluster region).
# Lymphatic endothelial wrapped to match LASP visual weight.
L1_DISPLAY_OVERRIDE = {
    "Luminal adaptive secretory precursor": "Luminal adaptive\nsecretory precursor",
    "Lymphatic endothelial": "Lymphatic\nendothelial",
}


def assign_v1_l1(lineage: str | float, l2s_label: str | float) -> str | None:
    """Map (joint l1p5_lineage, FLEX l2s_label) → V1-harmonized L1 display group.

    Lymphoid_NK collapses into T-lymphocyte per V1 nomenclature. Luminal
    splits into LASP / LHS by L2S label prefix. Endothelial splits into
    Vascular / Lymphatic by L2S label content. Lymphoid_B splits into
    B-lymphocyte / Plasma by L2S label content.
    """
    if pd.isna(lineage):
        return None
    lin = str(lineage)
    l2s = str(l2s_label) if pd.notna(l2s_label) else ""
    if lin == "Basal":
        return "Basal-myoepithelial"
    if lin == "Luminal":
        if l2s.startswith("LASP"):
            return "Luminal adaptive secretory precursor"
        if l2s.startswith("LHS"):
            return "Luminal hormone sensing"
        return None
    if lin == "Fibroblast":
        return "Fibroblast"
    if lin == "Pericyte":
        return "Perivascular"
    if lin == "Adipocyte":
        return "Adipocyte"
    if lin == "Endothelial":
        if "Lym" in l2s:
            return "Lymphatic endothelial"
        return "Vascular endothelial"
    if lin.startswith("Myeloid"):
        return "Myeloid"
    if lin == "Lymphoid_T" or lin == "Lymphoid_NK":
        return "T-lymphocyte"
    if lin == "Lymphoid_B":
        if "Plasma" in l2s:
            return "Plasma"
        return "B-lymphocyte"
    return None


# Per-V1-L1 nudges in DATA-UNIT increments (UMAP1 range ~19, UMAP2 range ~26;
# 2.5% ≈ 0.65 units, 5% ≈ 1.3, 10% ≈ 2.5). Positive y up, positive x right.
# Net positions reflect cumulative user-directed adjustments.
L1_NUDGE = {
    "Basal-myoepithelial":                  {"dy":  1.3},                            # +5% up
    "Luminal adaptive secretory precursor": {"dx": -1.0,  "dy":  0.26},              # -5% left, +1% up
    "Luminal hormone sensing":              {},
    "Fibroblast":                           {"dy":  0.65},                           # +2.5% up
    "Perivascular":                         {"dx": -0.19},                           # -1% left
    "Lymphatic endothelial":                {"dx": -0.19, "dy":  0.91},              # +3.5% up, -1% left
    "Myeloid":                              {"dy":  0.0},                            # net 0
    "T-lymphocyte":                         {"dy":  1.04, "dx": -2.0},               # +4% up, -10% left
    # B-lymphocyte: nudge toward Plasma cluster (centroid weighted by B>>Plasma).
    "B-lymphocyte":                         {"dx": -1.5,  "dy":  2.0},
}

# Per-label-num L2S nudges, applied AFTER repel.
# #2 = LHS-major, #3 = LASP-major, #11 = VSMC (Perivascular),
# #12 = Adipo-lipolysis, #14 = Lym-major (Lymphatic endothelial),
# #21 = Macro_FOLR2 (Myeloid), #23 = B_cell (B-lymphocyte area).
L2_NUDGE_BY_NUM = {
    2:  {"dy": -0.65},  # LHS-2 down 2.5%
    3:  {"dy": -0.26},  # LASP-3 down 1%
    11: {"dy": -0.65},  # PV-11 down 2.5%
    12: {"dy":  0.26},  # Adipo-12 up 1%
    14: {"dy": -0.65},  # LE-14 down 2.5%
    21: {"dy": -1.3},   # myeloid label-21 down 5%
    23: {"dy": -1.04},  # b-lymph label-23 down 4% (was -5%, +1%)
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--flex-umap", required=True, type=Path, help="canonical FLEX umap.csv")
    p.add_argument("--flex-obs", required=True, type=Path, help="canonical FLEX obs.csv")
    p.add_argument("--l2s-labels", required=True, type=Path, help="Track C flex_l2s_labels.csv")
    p.add_argument("--joint-l1p5", required=True, type=Path, help="joint_l1p5.csv (for L1 lineage outlines)")
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def build_substrate(umap_path: Path, obs_path: Path, l2s_path: Path, l1p5_path: Path) -> pd.DataFrame:
    """Inner-join FLEX UMAP × L2S labels; left-join joint l1p5_lineage."""
    log.info("loading FLEX umap.csv: %s", umap_path)
    umap = pd.read_csv(umap_path)
    if "UMAP1" not in umap.columns and "UMAP_1" in umap.columns:
        umap = umap.rename(columns={"UMAP_1": "UMAP1", "UMAP_2": "UMAP2"})
    log.info("  rows: %d", len(umap))

    log.info("loading FLEX obs.csv: %s", obs_path)
    obs = pd.read_csv(obs_path, usecols=lambda c: c in {"cell_id", "patient_id", "sample_id", "position"})
    log.info("  rows: %d", len(obs))

    log.info("loading L2S labels: %s", l2s_path)
    l2s = pd.read_csv(l2s_path, usecols=["cell_id", "compartment", "l2s_label"])
    log.info("  rows: %d", len(l2s))

    log.info("loading joint l1p5 (for L1 lineage outlines): %s", l1p5_path)
    l1p5 = pd.read_csv(l1p5_path, usecols=["cell_id", "l1p5_lineage"])
    log.info("  rows: %d", len(l1p5))

    sub = umap.merge(l2s, on="cell_id", how="inner").merge(obs, on="cell_id", how="left")
    sub = sub.merge(l1p5, on="cell_id", how="left")
    log.info("substrate rows after inner-join (umap × l2s) + left-join l1p5: %d", len(sub))

    sub["v1_l1"] = sub.apply(
        lambda r: assign_v1_l1(r["l1p5_lineage"], r["l2s_label"]), axis=1
    )

    label_str = sub["l2s_label"].astype(str)
    sub["is_artifact"] = (
        sub["l2s_label"].isna()
        | (label_str.str.lower() == "nan")
        | label_str.str.startswith("ARTIFACT_")
    )
    log.info("  artifact rows: %d", int(sub["is_artifact"].sum()))
    log.info("  cells without v1_l1: %d", int(sub["v1_l1"].isna().sum()))
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


def draw_v1_l1_labels(
    ax, dt: pd.DataFrame, label_fontsize: int = 7,
) -> list:
    """Place V1-harmonized L1 group labels at group medians. No KDE outlines.

    Each label position is the median UMAP coordinate of cells assigned to
    that V1 L1 group, then nudged per L1_NUDGE.
    """
    text_artists: list = []
    labeled = dt[(~dt["is_artifact"]) & dt["v1_l1"].notna()]
    for v1_l1, group in labeled.groupby("v1_l1"):
        if str(v1_l1) not in V1_L1_ORDER:
            continue
        pts = pd.DataFrame(group)[["UMAP1", "UMAP2"]].to_numpy(dtype=float)
        if len(pts) < 10:
            continue
        cx = float(np.median(pts[:, 0])); cy = float(np.median(pts[:, 1]))
        nudge = L1_NUDGE.get(str(v1_l1), {})
        cx += float(nudge.get("dx", 0.0))
        cy += float(nudge.get("dy", 0.0))
        display = L1_DISPLAY_OVERRIDE.get(str(v1_l1), str(v1_l1))
        t = ax.text(
            cx, cy, display,
            fontsize=label_fontsize, fontweight="bold",
            ha="center", va="center", zorder=11,
            bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="black", alpha=0.88, linewidth=0.5),
        )
        text_artists.append(t)
    return text_artists


def main() -> int:
    args = parse_args()

    dt = build_substrate(args.flex_umap, args.flex_obs, args.l2s_labels, args.joint_l1p5)

    label_order_per_comp: dict[str, list[str]] = {}
    for compartment in COMPARTMENT_ORDER:
        sub = dt[(dt["compartment"] == compartment) & (~dt["is_artifact"])]
        if sub.empty:
            label_order_per_comp[compartment] = []
            continue
        comp_counts = sub.groupby("l2s_label").size().sort_values(ascending=False)
        label_order_per_comp[compartment] = comp_counts.index.tolist()

    label_order: list[tuple[str, str]] = []
    for compartment in COMPARTMENT_ORDER:
        label_order.extend((compartment, lbl) for lbl in label_order_per_comp[compartment])
    label_to_num = {lbl: i + 1 for i, lbl in enumerate(label_order)}

    color_map: dict[tuple[str, str], str] = {}
    for compartment, labels in label_order_per_comp.items():
        for i, lbl in enumerate(labels):
            color_map[(compartment, lbl)] = SHARED_PALETTE[i % len(SHARED_PALETTE)]
    log.info("unique non-artifact (compartment, label) pairs: %d", len(color_map))

    centroids = (
        dt[~dt["is_artifact"]]
        .groupby(["compartment", "l2s_label"], as_index=False)
        .agg(UMAP1=("UMAP1", "median"), UMAP2=("UMAP2", "median"), n_cells=("cell_id", "size"))
    )
    centroids["label_num"] = centroids.apply(
        lambda r: label_to_num[(r["compartment"], r["l2s_label"])], axis=1
    )
    centroids = centroids.sort_values("label_num").reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(11.5, 8.0))
    ax.set_aspect("equal")
    ax.axis("off")

    art = dt[dt["is_artifact"]]
    if len(art):
        ax.scatter(art["UMAP1"], art["UMAP2"], c=ARTIFACT_COLOR, s=0.05, alpha=0.20, rasterized=True, linewidths=0)

    nonart = dt[~dt["is_artifact"]]
    for compartment, lbl in label_order:
        sub = nonart[(nonart["compartment"] == compartment) & (nonart["l2s_label"] == lbl)]
        if sub.empty:
            continue
        ax.scatter(
            sub["UMAP1"], sub["UMAP2"],
            c=color_map[(compartment, lbl)], s=0.11, alpha=0.95,
            rasterized=True, linewidths=0,
        )

    number_texts = []
    for _, r in centroids.iterrows():
        t = ax.text(
            r["UMAP1"], r["UMAP2"], str(int(r["label_num"])),
            fontsize=6, fontweight="bold", ha="center", va="center", zorder=12,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="black", alpha=0.92, linewidth=0.4),
        )
        number_texts.append(t)

    # Drop KDE outlines per design — labels remain at V1-harmonized L1 medians.
    draw_v1_l1_labels(ax, dt, label_fontsize=7)
    repel_labels(ax, number_texts, iters=200, step=0.05, min_sep=0.6)

    # Apply per-label-num L2 nudges AFTER repel (e.g. label-21 myeloid down 5%).
    for t in number_texts:
        try:
            num = int(t.get_text())
        except ValueError:
            continue
        nudge = L2_NUDGE_BY_NUM.get(num)
        if nudge is None:
            continue
        x, y = t.get_position()
        t.set_position((float(x) + float(nudge.get("dx", 0.0)),
                        float(y) + float(nudge.get("dy", 0.0))))

    import json as _json
    pos_dump = {
        "l2s_numbered": [
            {"text": t.get_text(), "x": float(t.get_position()[0]), "y": float(t.get_position()[1])}
            for t in number_texts
        ],
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "fig3b_flex_label_positions.json").write_text(_json.dumps(pos_dump, indent=2))

    # Single-column legend, shifted slightly left of the previous position.
    # Compartment-section headers as colorless Patches with bold label text;
    # blank Patch inserted between compartments for visual breathing room.
    handles: list = []
    for ci, compartment in enumerate(COMPARTMENT_ORDER):
        comp_centroids = centroids[centroids["compartment"] == compartment].sort_values("label_num")
        if comp_centroids.empty:
            continue
        if ci > 0:
            handles.append(mpatches.Patch(
                facecolor="white", edgecolor="white", label=" ",
            ))
        n_total = int(comp_centroids["n_cells"].sum())
        handles.append(mpatches.Patch(
            facecolor="white", edgecolor="white",
            label=r"$\bf{" + COMPARTMENT_DISPLAY[compartment].replace(" ", r"\ ") +
                  r"}$" + f"  (n={n_total:,})",
        ))
        for _, r in comp_centroids.iterrows():
            handles.append(mpatches.Patch(
                color=color_map[(compartment, r["l2s_label"])],
                label=f"  {int(r['label_num'])}. {r['l2s_label']}  ({int(r['n_cells']):,})",
            ))

    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.945, 1.0),
        fontsize=4.5, frameon=False, markerscale=1.6, handlelength=0.6, handletextpad=0.3,
        borderaxespad=0.0, ncol=1, labelspacing=0.20,
    )

    pdf_path = args.out_dir / "fig3b_flex_umap_l2s.pdf"
    png_path = args.out_dir / "fig3b_flex_umap_l2s.png"
    # DPI for rasterized scatter from aesthetics.yaml (rendering.dpi_umap_raster).
    import yaml as _yaml
    _aes_path = Path(__file__).resolve().parents[3] / "config" / "aesthetics.yaml"
    with open(_aes_path) as _fh:
        _raster_dpi = int(_yaml.safe_load(_fh).get("rendering", {}).get("dpi_umap_raster", 600))
    fig.savefig(pdf_path, bbox_inches="tight", dpi=_raster_dpi)
    fig.savefig(png_path, bbox_inches="tight", dpi=_raster_dpi)
    plt.close(fig)
    log.info("wrote %s + PNG", pdf_path)

    key_path = args.out_dir / "fig3b_flex_label_key.csv"
    centroids[["label_num", "compartment", "l2s_label", "n_cells", "UMAP1", "UMAP2"]].rename(
        columns={"UMAP1": "UMAP1_centroid", "UMAP2": "UMAP2_centroid"}
    ).to_csv(key_path, index=False)
    log.info("wrote label-num key %s", key_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
