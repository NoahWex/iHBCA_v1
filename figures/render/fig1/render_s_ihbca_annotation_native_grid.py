"""
Render s1.7 annotation-per-compartment supplementary figure.

Per-compartment 2×4 grid (3 compartments stacked top-to-bottom: epi, str, imm):
  Top-left:  L2 reference UMAP (canonical iHBCA L2 labels, this compartment)
  Remaining 7: per-study native author labels with gray backdrop

Same UMAP coords across all 8 sub-panels (the per-compartment scANVI UMAP).
Cells in the current study are colored by native label; out-of-study cells
shown as gray backdrop. Native label centroids labeled via adjustText.

Manual off-target filter (--exclude-yaml) drops native labels per (compartment,
study) that aren't appropriate for this compartment (Noah-curated list).

Substrate (resolved via CLI args; canonical promoted paths):
  --umap-epi    publication/analysis/annotation/scanvi/epi/n_latent_50/epi_scanvi_umap_barcoded.csv.gz
  --umap-str    publication/analysis/annotation/scanvi/str/n_latent_50/str_scanvi_umap_barcoded.csv.gz
  --umap-imm    publication/analysis/annotation/scanvi/imm/n_latent_50/imm_scanvi_umap_barcoded.csv.gz
  --native-csv  publication/analysis/annotation/native_labels.csv.gz
  --labels-csv  publication/analysis/annotation/labels_full.csv
  --exclude-yaml publication/figures/data/fig1/annotation_per_compartment/exclude_curated.yaml

Output (under --out-dir):
  s1_7_annotation_per_compartment.pdf
  s1_7_annotation_per_compartment.png

exclude_yaml schema (optional):
  epi:
    pal: [T_cell, B_cell]      # native labels to drop for Pal in epi
    reed: [Macrophage, ...]
  str: {...}
  imm: {...}
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
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("anno_supp")


# ----------------------------------------------------------------------------
# Locate project root + load aesthetics framework at module import time.
# Native-label palette is loaded from scales.qualitative_23 in aesthetics.yaml;
# label sets are study-specific so a key->color map would be wrong here —
# cycle-by-index is the correct semantic.
# ----------------------------------------------------------------------------

def _find_project_root() -> Path:
    here = Path(__file__).resolve().parent
    for p in [here] + list(here.parents):
        if (p / "publication" / "config" / "aesthetics.yaml").exists():
            return p
    raise FileNotFoundError("Could not locate publication/config/aesthetics.yaml")


_PROJECT_ROOT = _find_project_root()
sys.path.insert(0, str(_PROJECT_ROOT / "publication" / "config"))

from load_aesthetics import load_aesthetics  # noqa: E402

_CONFIG = load_aesthetics(config_dir=str(_PROJECT_ROOT / "publication" / "config"))
PALETTE = list(_CONFIG["scales"]["qualitative_23"])


PANEL_ID = "s_ihbca_annotation_native_grid"
OUTPUT_SLUG = "s1_7_annotation_per_compartment"

BACKDROP_COLOR = "#D8D8D8"
BACKDROP_SIZE = 0.30
BACKDROP_ALPHA = 0.35

FG_SIZE = 0.60
FG_ALPHA = 0.85

CENTROID_MIN_CELLS = 50

COMPARTMENT_DISPLAY = {"epi": "Epithelial", "str": "Stromal", "imm": "Immune"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--umap-epi", required=True, type=Path,
                   help="Epi scANVI UMAP CSV (cell_id + UMAP_1/2)")
    p.add_argument("--umap-str", required=True, type=Path,
                   help="Str scANVI UMAP CSV")
    p.add_argument("--umap-imm", required=True, type=Path,
                   help="Imm scANVI UMAP CSV")
    p.add_argument("--native-csv", required=True, type=Path,
                   help="wide-format native labels (cell_id + native_{study} columns)")
    p.add_argument("--labels-csv", required=True, type=Path,
                   help="labels_full.csv (cell_id, compartment, label, is_artifact)")
    p.add_argument("--exclude-yaml", type=Path, default=None,
                   help="optional manual off-target filter list")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--studies", default="gray,kumar,murrow,nee,pal,reed,twigger",
                   help="Comma-separated study list (fills the 7 per-study panels per row)")
    p.add_argument("--min-cells", type=int, default=CENTROID_MIN_CELLS,
                   help="Minimum cells per native label to draw a centroid label")
    p.add_argument("--adjust-iters", type=int, default=1000,
                   help="adjustText iter_lim (locked-PDF default)")
    p.add_argument("--adjust-force-text", type=float, default=3.0,
                   help="adjustText force_text (locked-PDF default)")
    return p.parse_args()


def load_substrate_for_compartment(compartment: str, umap_path: Path,
                                   labels: pd.DataFrame, native: pd.DataFrame) -> pd.DataFrame:
    log.info("[%s] Loading UMAP: %s", compartment, umap_path)
    umap = pd.read_csv(umap_path)
    rename = {c: c.replace("UMAP1", "UMAP_1").replace("UMAP2", "UMAP_2") for c in umap.columns}
    umap = umap.rename(columns=rename)
    log.info("[%s]   %d UMAP rows", compartment, len(umap))

    comp_labels = labels[labels["compartment"] == compartment]
    log.info("[%s]   %d cells in labels_full (post artifact drop)", compartment, len(comp_labels))

    df = comp_labels.merge(umap, on="cell_id", how="inner")
    df = df.merge(native, on="cell_id", how="left")
    log.info("[%s]   joined: %d cells", compartment, len(df))
    return df


def load_shared_substrates(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load labels_full + native (wide) once; per-compartment UMAPs joined per row."""
    log.info("Loading labels_full: %s", args.labels_csv)
    labels = pd.read_csv(args.labels_csv, usecols=["cell_id", "compartment", "label", "is_artifact"])
    labels = labels[~labels["is_artifact"].fillna(False)]
    labels = labels.rename(columns={"label": "L2"})
    log.info("  %d non-artifact cells", len(labels))

    log.info("Loading native (wide): %s", args.native_csv)
    native = pd.read_csv(args.native_csv, low_memory=False)
    native_cols = [c for c in native.columns if c.startswith("native_")]
    log.info("  %d cells; %d native_* columns", len(native), len(native_cols))
    return labels, native


def apply_exclude_filter(df: pd.DataFrame, compartment: str, exclude_yaml: Path | None,
                          studies: list[str]) -> pd.DataFrame:
    if exclude_yaml is None or not exclude_yaml.exists():
        return df

    cfg = yaml.safe_load(exclude_yaml.read_text()) or {}
    rules = cfg.get(compartment, {}) or {}
    if not rules:
        return df
    log.info("[%s] Exclude rules: %s", compartment, {k: len(v) for k, v in rules.items()})

    for study in studies:
        col = f"native_{study}"
        if col not in df.columns:
            continue
        drops = rules.get(study, [])
        if not drops:
            continue
        mask = df[col].isin(drops)
        n_drop = int(mask.sum())
        if n_drop > 0:
            df.loc[mask, col] = np.nan
            log.info("[%s]   %s: nulled %d cells matching %d excluded labels",
                     compartment, study, n_drop, len(drops))
    return df


def project_centroids(ax, df: pd.DataFrame, lbl_col: str, min_cells: int,
                      iters: int, force_text: float) -> None:
    """Draw label-name text at each label's median position, with adjustText repel."""
    grp = df.groupby(lbl_col).agg(
        UMAP_1=("UMAP_1", "median"),
        UMAP_2=("UMAP_2", "median"),
        n=("cell_id", "size"),
    ).reset_index()
    grp = grp[grp["n"] >= min_cells]
    if grp.empty:
        return

    texts = []
    for _, r in grp.iterrows():
        t = ax.text(
            float(r["UMAP_1"]), float(r["UMAP_2"]), str(r[lbl_col]),
            fontsize=5, fontweight="bold", ha="center", va="center", zorder=10,
            bbox=dict(boxstyle="round,pad=0.18", facecolor="white",
                      edgecolor="black", alpha=0.85, linewidth=0.3),
        )
        texts.append(t)

    try:
        from adjustText import adjust_text
        adjust_text(
            texts, ax=ax,
            expand_text=(2.2, 2.4), expand_points=(1.6, 1.8),
            force_text=force_text, force_points=1.0,
            arrowprops=dict(arrowstyle="-", color="grey", lw=0.3, alpha=0.5),
            iter_lim=iters,
            only_move={"text": "xy"},
        )
    except ImportError:
        log.warning("adjustText not installed; labels not repelled.")


def render_l2_panel(ax, df: pd.DataFrame, compartment: str, min_cells: int,
                    iters: int, force_text: float, show_title: bool = True) -> None:
    """Leftmost L2 reference panel for this row (one row per compartment)."""
    labels = sorted(df["L2"].dropna().unique().tolist())
    color_map = {lbl: PALETTE[i % len(PALETTE)] for i, lbl in enumerate(labels)}
    colors = df["L2"].map(color_map).fillna(BACKDROP_COLOR).values

    ax.scatter(df["UMAP_1"], df["UMAP_2"], c=colors, s=FG_SIZE, alpha=FG_ALPHA,
               rasterized=True, linewidths=0)
    project_centroids(ax, df.assign(_lbl=df["L2"]), "_lbl", min_cells, iters, force_text)

    if show_title:
        ax.set_title("iHBCA L2", fontsize=7, fontweight="bold")
    # Compartment label on the y-axis (row label, publication style)
    ax.set_ylabel(COMPARTMENT_DISPLAY[compartment], fontsize=8, fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.3)


def render_study_panel(ax, df: pd.DataFrame, study: str, min_cells: int,
                       iters: int, force_text: float, show_title: bool = True) -> None:
    """Per-study panel: gray backdrop + colored study cells + centroid labels."""
    col = f"native_{study}"
    if col not in df.columns:
        ax.text(0.5, 0.5, f"{study}\n(column missing)", ha="center", va="center",
                fontsize=6, transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        return

    # Backdrop: all compartment cells
    ax.scatter(df["UMAP_1"], df["UMAP_2"], c=BACKDROP_COLOR, s=BACKDROP_SIZE,
               alpha=BACKDROP_ALPHA, rasterized=True, linewidths=0)

    fg = df[df[col].notna()].copy()
    if len(fg) == 0:
        if show_title:
            ax.set_title(study, fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
        return

    labels = sorted(fg[col].astype(str).unique())
    color_map = {lbl: PALETTE[i % len(PALETTE)] for i, lbl in enumerate(labels)}
    fg_colors = fg[col].astype(str).map(color_map).values

    ax.scatter(fg["UMAP_1"], fg["UMAP_2"], c=fg_colors, s=FG_SIZE, alpha=FG_ALPHA,
               rasterized=True, linewidths=0)
    project_centroids(ax, fg.assign(_lbl=fg[col].astype(str)), "_lbl",
                      min_cells, iters, force_text)

    if show_title:
        ax.set_title(study, fontsize=7, fontweight="bold")
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.3)


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    studies = [s.strip() for s in args.studies.split(",") if s.strip()]
    if len(studies) != 7:
        log.warning("Expected 7 studies, got %d: %s", len(studies), studies)

    labels, native = load_shared_substrates(args)

    compartments = [
        ("epi", args.umap_epi),
        ("str", args.umap_str),
        ("imm", args.umap_imm),
    ]

    # 3 columns (one per compartment), 4 rows × 2 sub-cols each.
    # Per compartment: L2 reference + 7 study panels = 8 panels, 2 per row.
    fig = plt.figure(figsize=(15, 11))
    subfigs = fig.subfigures(1, 3, wspace=0.04)

    for comp_idx, (comp, umap_path) in enumerate(compartments):
        df = load_substrate_for_compartment(comp, umap_path, labels, native)
        df = apply_exclude_filter(df, comp, args.exclude_yaml, studies)

        sf = subfigs[comp_idx]
        sf.suptitle(COMPARTMENT_DISPLAY[comp], fontsize=9, fontweight="bold", y=0.99)
        axs = sf.subplots(4, 2, gridspec_kw={"hspace": 0.25, "wspace": 0.08})
        sub_axes = [axs[r, c] for r in range(4) for c in range(2)]

        render_l2_panel(sub_axes[0], df, comp,
                        args.min_cells, args.adjust_iters, args.adjust_force_text,
                        show_title=True)
        sub_axes[0].set_ylabel("")
        for i, study in enumerate(studies):
            render_study_panel(sub_axes[i + 1], df, study,
                               args.min_cells, args.adjust_iters, args.adjust_force_text,
                               show_title=True)

    pdf_path = args.out_dir / f"{OUTPUT_SLUG}.pdf"
    png_path = args.out_dir / f"{OUTPUT_SLUG}.png"
    fig.savefig(pdf_path, bbox_inches="tight", dpi=600)
    fig.savefig(png_path, bbox_inches="tight", dpi=300)
    plt.close(fig)
    log.info("Wrote %s (%d bytes)", pdf_path, pdf_path.stat().st_size)
    log.info("Wrote %s (%d bytes)", png_path, png_path.stat().st_size)
    return 0


if __name__ == "__main__":
    sys.exit(main())
