"""
render_supp7_1_concord_receipt.py — Fig 3 Supp 7-1.

Concord integration receipt: two facet grids over the joint FLEX + Xenium
UMAP embedding.

    Top row 1x2 (FLEX | Xenium)   coloured by patient
                                  -> patients overlay across platforms,
                                     showing Concord aligned them

    Bottom 5x2 position grid      Xenium cells only
                                  coloured by patient (shared legend)
                                  -> per-position views show patient
                                     overlap within each anatomic site

Inputs:
    paths.joint_v6_umap            cell_id, UMAP1, UMAP2 (joint embedding)
    paths.xenium_joint_l1p5        cell_id, platform, patient_id,
                                   is_artifact
    paths.xenium_cell_annotations  cell_id, position (full anatomic
                                   position string with P3 anterior/
                                   middle/posterior depth split)
    paths.joint_concord_config     Concord training config (referenced
                                   in index.yaml only, not loaded)

Outputs:
    supp7_1_concord_receipt.pdf

Substrate note: paths.joint_v6_umap and paths.joint_v6_obs resolve to
publication/preprocessing/integration/joint_concord/embeddings/ as of
2026-05-07 (substrate promotion logged in publication/manifest.yaml,
audit at coordination/audit/supp7_joint_v6_substrate_audit_20260507.md).
joint_v6_latent stays in the Spatial_HBCA_Xenium dev repo (1.2 GB; only
consumed for re-training, not at render time).

Palette note: this version colours both rows by patient (using the
existing 'patient' palette in aesthetics.yaml). Position is now an
axis-facet rather than a colour, so no 'position' palette is required.

Drawing order: rows are shuffled with a fixed seed before scatter so no
single platform / position dominates the visible top of the rasterised
layer. Artifacts go down first under all categories.

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
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT_DEFAULT = SCRIPT_DIR.parents[3]
sys.path.insert(0, str(PROJECT_ROOT_DEFAULT / "publication" / "config"))
from load_aesthetics import (  # noqa: E402
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
log = logging.getLogger("supp7_1")


QUADRANT_FULL = {
    "LIQ": "Lower Inner",
    "LOQ": "Lower Outer",
    "UIQ": "Upper Inner",
    "UOQ": "Upper Outer",
    "Lower": "Lower",      # P2 sub-quadrant unspecified
    "Upper": "Upper",
}

DEPTH_FULL = {
    "A": "Anterior",
    "M": "Middle",
    "P": "Posterior",
}

ARTIFACT_COLOR = "lightgrey"
UNMAPPED_COLOR = "darkgrey"

PLATFORM_ORDER = ("flex", "xenium")
PATIENT_ORDER = ("Pat1", "Pat2", "UCI604", "UCI220228")

# 19 Xenium anatomic positions, ordered by p_level then depth (A->M->P) then
# quadrant (UO->UI->LO->LI). 5x4 grid in row-major order with one trailing
# empty slot (slot 20).
POSITION_ORDER = (
    "P1",
    "P2 Upper Outer", "P2 Upper Inner", "P2 Lower Outer", "P2 Lower Inner",
    "P2 Upper", "P2 Lower",
    "P3 Anterior Upper Outer", "P3 Anterior Upper Inner",
    "P3 Anterior Lower Outer", "P3 Anterior Lower Inner",
    "P3 Middle Upper Outer", "P3 Middle Upper Inner",
    "P3 Middle Lower Outer", "P3 Middle Lower Inner",
    "P3 Posterior Upper Outer", "P3 Posterior Upper Inner",
    "P3 Posterior Lower Outer", "P3 Posterior Lower Inner",
)

SHUFFLE_SEED = 42


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--test", action="store_true",
                   help="Subsample 50K cells for quick sanity render.")
    return p.parse_args()


def expand_position(pos) -> str:
    """Expand `xenium_cell_annotations.position` codes to full names.

    Examples:
        FLEX             -> 'FLEX'
        P1               -> 'P1'
        P2_LIQ           -> 'P2 Lower Inner'
        P2_Lower         -> 'P2 Lower' (P2 cells with no UQ/LQ split)
        P2_Upper         -> 'P2 Upper'
        P3_A_LIQ         -> 'P3 Anterior Lower Inner'
        P3_M_UOQ         -> 'P3 Middle Upper Outer'
        P3_P_LOQ         -> 'P3 Posterior Lower Outer'
        NaN / ''         -> 'Unknown'
    """
    if pos is None or (isinstance(pos, float) and pd.isna(pos)):
        return "Unknown"
    s = str(pos).strip()
    if s == "" or s.lower() == "nan":
        return "Unknown"
    if s == "FLEX" or s == "P1":
        return s
    parts = s.split("_")
    level = parts[0]
    if len(parts) == 2:
        # P2_LIQ / P2_Upper / P2_Lower
        suffix = QUADRANT_FULL.get(parts[1], parts[1])
        return f"{level} {suffix}"
    if len(parts) == 3:
        # P3_A_LIQ
        depth = DEPTH_FULL.get(parts[1], parts[1])
        quadrant = QUADRANT_FULL.get(parts[2], parts[2])
        return f"{level} {depth} {quadrant}"
    return s


def scatter_subpanel(
    ax,
    dt: pd.DataFrame,
    color_col: str,
    color_map: dict,
    point_size: float,
    title: str,
    full_extent: tuple,
) -> None:
    """One UMAP scatter subpanel.

    Rows are shuffled with SHUFFLE_SEED before drawing so categories
    interleave; artifacts go down first underneath.
    """
    if len(dt) == 0:
        ax.set_xlim(full_extent[0], full_extent[1])
        ax.set_ylim(full_extent[2], full_extent[3])
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlabel(title, fontsize=7)
        return

    dt = dt.sample(frac=1, random_state=SHUFFLE_SEED)
    art_mask = dt["is_artifact"].astype(bool)
    art = dt[art_mask]
    nonart = dt[~art_mask]

    if len(art):
        ax.scatter(
            art["UMAP1"], art["UMAP2"],
            c=ARTIFACT_COLOR, s=point_size, alpha=0.20,
            rasterized=True, linewidths=0,
        )
    if len(nonart):
        colors = [color_map.get(v, UNMAPPED_COLOR) for v in nonart[color_col]]
        ax.scatter(
            nonart["UMAP1"], nonart["UMAP2"],
            c=colors, s=point_size, alpha=0.7,
            rasterized=True, linewidths=0,
        )

    ax.set_xlim(full_extent[0], full_extent[1])
    ax.set_ylim(full_extent[2], full_extent[3])
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel(title, fontsize=7)


def main() -> int:
    args = parse_args()
    config_dir = str(args.project_root / "publication" / "config")
    paths = load_paths(config_dir=config_dir)
    aes = load_aesthetics(config_dir=config_dir)

    umap_path = resolve_path("joint_v6_umap", paths)
    l1p5_path = resolve_path("xenium_joint_l1p5", paths)
    log.info("loading joint UMAP coords: %s", umap_path)
    umap = pd.read_csv(umap_path)
    log.info("  rows: %d", len(umap))

    log.info("loading joint substrate: %s", l1p5_path)
    l1p5 = pd.read_csv(
        l1p5_path,
        usecols=["cell_id", "platform", "patient_id", "is_artifact"],
    )
    log.info("  rows: %d | platforms: %s | patients: %s",
             len(l1p5), sorted(l1p5["platform"].dropna().unique()),
             sorted(l1p5["patient_id"].dropna().unique()))

    pos_path = resolve_path("xenium_cell_annotations", paths)
    log.info("loading position substrate: %s", pos_path)
    pos = pd.read_csv(pos_path, usecols=["cell_id", "position"])
    log.info("  rows: %d | distinct positions: %d", len(pos), pos["position"].nunique())

    dt = umap.merge(l1p5, on="cell_id", how="inner").merge(
        pos, on="cell_id", how="left"
    )
    log.info("3-way join: %d cells", len(dt))

    if args.test:
        n_sub = min(50_000, len(dt))
        dt = dt.sample(n=n_sub, random_state=42)
        log.info("--test subsample: %d", len(dt))

    dt["position_full"] = dt["position"].apply(expand_position)

    # ---- Theme + patient palette ----
    plt.rcParams.update(get_matplotlib_theme(aes))
    patient_palette_full = get_palette("patient", aes)
    actual_patients = sorted(dt["patient_id"].dropna().unique())
    patient_color_map = {
        p: patient_palette_full[p]
        for p in actual_patients if p in patient_palette_full
    }
    missing_patients = [p for p in actual_patients if p not in patient_palette_full]
    if missing_patients:
        log.warning("patients without palette entry: %s", missing_patients)

    xen_only = dt[dt["platform"] == "xenium"]
    xen_nonart = xen_only[~xen_only["is_artifact"].astype(bool)]
    log.info("Xenium non-artifact cells: %d", len(xen_nonart))
    pos_present = set(xen_nonart["position_full"].dropna().unique())
    pos_missing_from_canonical = [p for p in pos_present if p not in POSITION_ORDER]
    if pos_missing_from_canonical:
        log.warning("positions in data but not in POSITION_ORDER: %s",
                    pos_missing_from_canonical)

    # ---- Layout ----
    # Row 0: 1x2 platform sub-panels (left) + patient side-legend (right)
    #        with (FLEX | Xenium) counts -- legend keyed to top plot.
    # Row 1: 5x4 position sub-panels (Xenium only), full width.
    # Row 2: patient colour-key legend strip (full-width horizontal) keyed
    #        to the bottom plot (no counts; per-patient counts already in
    #        each position panel header).
    fig_w = 7.2
    fig_h = 9.6
    fig = plt.figure(figsize=(fig_w, fig_h))
    outer = fig.add_gridspec(
        nrows=3, ncols=2,
        height_ratios=[3.0, 6.0, 0.6],
        width_ratios=[4.8, 2.4],
        hspace=0.18, wspace=0.05,
        left=0.02, right=0.98, top=0.99, bottom=0.01,
    )
    top_grid = outer[0, 0].subgridspec(1, 2, hspace=0.0, wspace=0.10)
    bot_grid = outer[1, :].subgridspec(4, 5, hspace=0.45, wspace=0.10)
    ax_top_legend = fig.add_subplot(outer[0, 1])
    ax_bot_legend = fig.add_subplot(outer[2, :])
    ax_top_legend.axis("off")
    ax_bot_legend.axis("off")
    log.info("figure: %.2fin x %.2fin", fig_w, fig_h)

    # Shared UMAP extent so all sub-panels overlay onto the same coordinate frame
    pad = 0.05
    x_min = dt["UMAP1"].min() - pad * (dt["UMAP1"].max() - dt["UMAP1"].min())
    x_max = dt["UMAP1"].max() + pad * (dt["UMAP1"].max() - dt["UMAP1"].min())
    y_min = dt["UMAP2"].min() - pad * (dt["UMAP2"].max() - dt["UMAP2"].min())
    y_max = dt["UMAP2"].max() + pad * (dt["UMAP2"].max() - dt["UMAP2"].min())
    full_extent = (x_min, x_max, y_min, y_max)

    pt_size = 0.05

    # ---- Top 1x2: platform sub-panels, colour by patient ----
    for c, plat in enumerate(PLATFORM_ORDER):
        ax = fig.add_subplot(top_grid[0, c])
        plat_dt = dt[dt["platform"] == plat]
        if len(plat_dt) == 0:
            ax.axis("off")
            continue
        n_total = int((~plat_dt["is_artifact"].astype(bool)).sum())
        title = "FLEX" if plat == "flex" else "Xenium"
        scatter_subpanel(ax, plat_dt, "patient_id", patient_color_map, pt_size,
                         f"{title} (n = {n_total:,})", full_extent)

    # ---- Top legend (side, keyed to top plot): patient × platform counts ----
    nonart_all = dt[~dt["is_artifact"].astype(bool)]
    top_handles = []
    for p in PATIENT_ORDER:
        if p not in patient_color_map:
            continue
        n_flex = int(((nonart_all["patient_id"] == p) & (nonart_all["platform"] == "flex")).sum())
        n_xen = int(((nonart_all["patient_id"] == p) & (nonart_all["platform"] == "xenium")).sum())
        top_handles.append(
            mpatches.Patch(color=patient_color_map[p],
                           label=f"{p} ({n_flex:,} | {n_xen:,})")
        )
    ax_top_legend.legend(
        handles=top_handles, loc="center", fontsize=7,
        frameon=False, handlelength=1.0, handletextpad=0.4,
        labelspacing=0.5,
        title="Patient (FLEX | Xenium)", title_fontsize=7,
    )

    # ---- Bottom 5x4: position sub-panels, Xenium only, colour by patient ----
    # Sub-panel header is 2-line: "{position}\n({n_pat1:,} | {n_pat2:,} | ...)"
    # listing per-patient counts in PATIENT_ORDER (matching the top legend).
    def _fmt_count(n: int) -> str:
        if n >= 1000:
            return f"{n / 1000:.1f}k".replace(".0k", "k")
        return str(n)

    patient_order_present = [p for p in PATIENT_ORDER if p in patient_color_map]
    for idx, position in enumerate(POSITION_ORDER):
        r, c = divmod(idx, 5)
        ax = fig.add_subplot(bot_grid[r, c])
        pos_dt = xen_only[xen_only["position_full"] == position]
        if len(pos_dt) == 0:
            ax.axis("off")
            continue
        pos_nonart = pos_dt[~pos_dt["is_artifact"].astype(bool)]
        per_pat_counts = [
            int((pos_nonart["patient_id"] == p).sum())
            for p in patient_order_present
        ]
        counts_str = " | ".join(_fmt_count(n) for n in per_pat_counts)
        title = f"{position}\n({counts_str})"
        scatter_subpanel(ax, pos_dt, "patient_id", patient_color_map, pt_size,
                         title, full_extent)

    # blank slot 20 (5x4 = 20, only 19 positions)
    last_idx = len(POSITION_ORDER)
    if last_idx < 20:
        r, c = divmod(last_idx, 5)
        ax_blank = fig.add_subplot(bot_grid[r, c])
        ax_blank.axis("off")

    # ---- Bottom strip: patient colour-key legend keyed to position grid ----
    # Per-patient cell counts within each position are shown in each position
    # sub-panel header; this legend is purely a colour key.
    bot_handles = [
        mpatches.Patch(color=patient_color_map[p], label=p)
        for p in PATIENT_ORDER if p in patient_color_map
    ]
    ax_bot_legend.legend(
        handles=bot_handles, loc="center",
        ncol=len(bot_handles), fontsize=7,
        frameon=False, handlelength=1.2, handletextpad=0.5,
        columnspacing=2.0,
        title=f"Position — Xenium cells only "
              f"(counts shown as: "
              f"{' | '.join(p for p in PATIENT_ORDER if p in patient_color_map)})",
        title_fontsize=7,
    )

    # ---- Save ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / "supp7_1_concord_receipt.pdf"
    _raster_dpi = int(aes.get("rendering", {}).get("dpi_umap_raster", 600))
    fig.savefig(pdf_path, bbox_inches="tight", dpi=_raster_dpi)
    plt.close(fig)
    log.info("wrote %s", pdf_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
