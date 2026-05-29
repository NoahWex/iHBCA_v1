"""
00_build_structural_package.py
Stage 2 of annotation_pipeline_v2 — pre-annotation structural evidence.

Produces for a given compartment:
  structural/study_composition.csv         cells per (study, donor)
  structural/study_composition.pdf         per-study barplot
  structural/leiden_umap_grid.pdf          UMAPs across all Leiden resolutions
  structural/native_label_umap_grid.pdf    per-study UMAPs with native labels projected on
  structural/native_label_f1_res{X}.pdf    native-label × cluster F1 heatmap per resolution
  structural/native_label_f1_res{X}.csv    underlying F1 matrix

Inputs (paths via CLI — no hardcoded absolute paths in the script body):
  --compartment {imm,epi,str}
  --leiden-csv       scANVI Leiden CSV (cell_id + leiden_{res}...)
  --umap-csv         scANVI UMAP CSV (cell_id + UMAP_1 + UMAP_2)
  --metadata-csv     per-compartment metadata_enriched.csv (cell_id + dataset + patientID ...)
  --native-labels-csv  all-cells native labels CSV (cell_id + native_{study}...)
  --config           compartment_config.yaml (for native_labels.studies list)
  --out-dir          output directory (structural/ subdir will be created)

Design doc: iHBCA_publication/coordination/plans/annotation_pipeline_v2.md §4.1
Stewardship: pattern-ref to extract_native_labels.sh for container env; fresh
logic for composition/UMAP/F1 plotting (no working pattern in V1_Annotation
that matches all four outputs).
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import List

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("build_structural_package")


STUDY_COLORS = {
    "kumar":   "#1f77b4",
    "reed":    "#ff7f0e",
    "nee":     "#2ca02c",
    "gray":    "#d62728",
    "murrow":  "#9467bd",
    "twigger": "#8c564b",
    "pal":     "#e377c2",
}


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compartment", required=True, choices=["imm", "epi", "str"])
    p.add_argument("--leiden-csv", required=True, type=Path)
    p.add_argument("--umap-csv", required=True, type=Path)
    p.add_argument("--metadata-csv", required=True, type=Path)
    p.add_argument("--native-labels-csv", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path, help="compartment_config.yaml")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument(
        "--subsample",
        type=int,
        default=None,
        help="If set, subsample to N cells stratified by leiden_1.0 for testing. Omit for full run.",
    )
    return p.parse_args()


# =============================================================================
# Loading
# =============================================================================

def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, dict, List[str]]:
    """Return (joined_df, config, resolution_names)."""
    log.info("Loading config: %s", args.config)
    cfg = yaml.safe_load(args.config.read_text())

    log.info("Loading Leiden: %s", args.leiden_csv)
    leiden = pd.read_csv(args.leiden_csv)
    res_cols = [c for c in leiden.columns if c.startswith("leiden_")]
    log.info("  %d cells, %d resolutions: %s", len(leiden), len(res_cols), res_cols)

    log.info("Loading UMAP: %s", args.umap_csv)
    umap = pd.read_csv(args.umap_csv)

    log.info("Loading metadata: %s", args.metadata_csv)
    meta_cols = ["cell_id", "dataset", "patientID"]
    meta = pd.read_csv(args.metadata_csv, usecols=meta_cols, low_memory=False)

    log.info("Loading native labels: %s", args.native_labels_csv)
    native = pd.read_csv(args.native_labels_csv, low_memory=False)
    native_study_cols = [c for c in native.columns if c.startswith("native_")]
    log.info("  native label columns: %s", native_study_cols)

    log.info("Joining by cell_id")
    df = leiden.merge(umap, on="cell_id", how="inner")
    df = df.merge(meta, on="cell_id", how="inner")
    df = df.merge(native, on="cell_id", how="left")  # left: some cells may not have native label
    log.info("Joined rows: %d", len(df))

    if args.subsample:
        log.info("Subsampling to %d cells (stratified by leiden_1.0)", args.subsample)
        # Use leiden_1.0 if present, else first resolution
        strat = "leiden_1.0" if "leiden_1.0" in df.columns else res_cols[0]
        df = (
            df.groupby(strat, group_keys=False)
            .apply(lambda g: g.sample(min(len(g), max(1, args.subsample // df[strat].nunique())), random_state=42))
            .reset_index(drop=True)
        )
        log.info("  post-subsample: %d cells", len(df))

    return df, cfg, res_cols


# =============================================================================
# Output 1: Study composition
# =============================================================================

def emit_study_composition(df: pd.DataFrame, out_dir: Path) -> None:
    log.info("Emitting study composition")
    tab = df.groupby(["dataset", "patientID"], observed=True).size().reset_index(name="n_cells")
    csv_path = out_dir / "study_composition.csv"
    tab.to_csv(csv_path, index=False)
    log.info("  wrote %s", csv_path)

    # Per-study totals for the barplot
    per_study = df["dataset"].value_counts().sort_values(ascending=False)
    per_study_donors = df.groupby("dataset", observed=True)["patientID"].nunique().reindex(per_study.index)

    fig, ax = plt.subplots(1, 1, figsize=(7, 3.5))
    bars = ax.bar(per_study.index, per_study.values, color="#4a7ba6", edgecolor="black", linewidth=0.5)
    for bar, n_donor in zip(bars, per_study_donors.values):
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + per_study.values.max() * 0.02,
            f"{int(h):,}\n({n_donor}d)",
            ha="center",
            va="bottom",
            fontsize=6,
        )
    ax.set_ylabel("cells", fontsize=7)
    ax.set_xlabel("study (dataset)", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_ylim(0, per_study.values.max() * 1.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    pdf_path = out_dir / "study_composition.pdf"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", pdf_path)


# =============================================================================
# Output 2: Leiden UMAP grid (all resolutions, colored by cluster)
# =============================================================================

def _categorical_palette(n: int) -> List[tuple]:
    """Return n distinct RGB tuples cycling through tab20 + tab20b + tab20c."""
    base = (
        list(colormaps["tab20"].colors)
        + list(colormaps["tab20b"].colors)
        + list(colormaps["tab20c"].colors)
    )
    if n <= len(base):
        return base[:n]
    hsv = colormaps["hsv"]
    return [hsv(i / n) for i in range(n)]


def _scatter_clusters(ax, df, cluster_col, title):
    clusters = sorted(df[cluster_col].astype(str).unique(), key=lambda s: (len(s), s))
    palette = _categorical_palette(len(clusters))
    color_map = dict(zip(clusters, palette))
    colors = df[cluster_col].astype(str).map(color_map).values
    ax.scatter(df["UMAP_1"], df["UMAP_2"], c=colors, s=0.5, alpha=0.6, rasterized=True, linewidths=0)
    ax.set_title(title, fontsize=7)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.3)
    # Minimal cluster count annotation
    ax.text(
        0.02,
        0.98,
        f"k={len(clusters)}",
        transform=ax.transAxes,
        fontsize=5,
        va="top",
        ha="left",
        color="black",
        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7),
    )


def emit_leiden_umap_grid(df: pd.DataFrame, res_cols: List[str], out_dir: Path) -> None:
    log.info("Emitting Leiden UMAP grid (%d resolutions)", len(res_cols))
    n = len(res_cols)
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.5, nrows * 2.5))
    axes = np.atleast_2d(axes).flatten()

    for i, col in enumerate(res_cols):
        title = col.replace("leiden_", "res=")
        _scatter_clusters(axes[i], df, col, title)
    for j in range(len(res_cols), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    pdf_path = out_dir / "leiden_umap_grid.pdf"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", pdf_path)


# =============================================================================
# Output 3: Native-label UMAP grid (per-study, labels projected on plot)
# =============================================================================

def _project_label_centroids(ax, df, label_col, min_cells: int = 50) -> None:
    """Annotate ax with median-coord label text, with mild repulsion so labels
    don't stack. Uses adjustText when available; degrades to static placement."""
    sub = df[df[label_col].notna()]
    if sub.empty:
        return
    counts = sub[label_col].value_counts()
    shown = counts[counts >= min_cells].index

    texts = []
    for label in shown:
        cells = sub[sub[label_col] == label]
        x = float(cells["UMAP_1"].median())
        y = float(cells["UMAP_2"].median())
        t = ax.text(
            x,
            y,
            str(label),
            fontsize=5,
            ha="center",
            va="center",
            color="black",
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="black", lw=0.25, alpha=0.85),
            zorder=10,
        )
        texts.append(t)

    if not texts:
        return

    try:
        from adjustText import adjust_text

        adjust_text(
            texts,
            ax=ax,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.25, alpha=0.6),
            expand_text=(1.05, 1.1),
            expand_points=(1.05, 1.1),
            force_text=(0.2, 0.3),
            force_points=(0.05, 0.1),
            only_move={"text": "xy"},
            lim=50,
        )
    except ImportError:
        log.warning("adjustText not installed; label centroids rendered without repulsion")


def emit_native_label_umap_grid(df: pd.DataFrame, cfg: dict, out_dir: Path) -> None:
    studies = cfg["pipeline"]["native_labels"]["studies"]
    prefix = cfg["pipeline"]["native_labels"]["column_prefix"]
    log.info("Emitting native-label UMAP grid (%d studies)", len(studies))

    ncols = min(4, len(studies))
    nrows = int(np.ceil(len(studies) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.2, nrows * 3.2))
    axes = np.atleast_2d(axes).flatten()

    # Background cells (all, gray)
    for i, study in enumerate(studies):
        ax = axes[i]
        ax.scatter(df["UMAP_1"], df["UMAP_2"], c="#d0d0d0", s=0.3, alpha=0.4, rasterized=True, linewidths=0)

        col = f"{prefix}{study}"
        if col not in df.columns:
            ax.set_title(f"{study}\n(column missing)", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        study_cells = df[df[col].notna()]
        n = len(study_cells)
        if n == 0:
            ax.set_title(f"{study}\n(no native labels)", fontsize=7)
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        labels = sorted(study_cells[col].astype(str).unique())
        palette = _categorical_palette(len(labels))
        color_map = dict(zip(labels, palette))
        colors = study_cells[col].astype(str).map(color_map).values
        ax.scatter(
            study_cells["UMAP_1"],
            study_cells["UMAP_2"],
            c=colors,
            s=0.6,
            alpha=0.8,
            rasterized=True,
            linewidths=0,
        )

        _project_label_centroids(ax, df.assign(_lbl=df[col]), "_lbl", min_cells=50)

        ax.set_title(f"{study}  n={n:,}  labels={len(labels)}", fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_linewidth(0.3)

    for j in range(len(studies), len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    pdf_path = out_dir / "native_label_umap_grid.pdf"
    fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", pdf_path)


# =============================================================================
# Output 4: Native-label × cluster F1 heatmap (per resolution)
# =============================================================================

def _compute_f1_matrix(df: pd.DataFrame, native_col: str, cluster_col: str) -> pd.DataFrame:
    """F1(label, cluster) treating both as binary memberships over all cells."""
    sub = df[df[native_col].notna()].copy()
    if sub.empty:
        return pd.DataFrame()

    sub[native_col] = sub[native_col].astype(str)
    sub[cluster_col] = sub[cluster_col].astype(str)

    # Contingency: rows = native label, cols = cluster, values = count of cells in both
    ct = pd.crosstab(sub[native_col], sub[cluster_col])
    label_totals = ct.sum(axis=1)  # cells per native label
    cluster_totals = sub[cluster_col].value_counts().reindex(ct.columns).fillna(0)

    # Vectorized F1: precision = tp/cluster_total, recall = tp/label_total, f1 = 2pr/(p+r)
    tp = ct.values.astype(float)
    prec = np.divide(tp, cluster_totals.values[np.newaxis, :], out=np.zeros_like(tp), where=cluster_totals.values[np.newaxis, :] > 0)
    rec = np.divide(tp, label_totals.values[:, np.newaxis], out=np.zeros_like(tp), where=label_totals.values[:, np.newaxis] > 0)
    denom = prec + rec
    f1 = np.divide(2 * prec * rec, denom, out=np.zeros_like(tp), where=denom > 0)
    return pd.DataFrame(f1, index=ct.index, columns=ct.columns)


def _diagonal_ordering(mat: pd.DataFrame) -> tuple[list, list]:
    """Row/column permutations that push high-F1 cells onto the diagonal.

    Row order: sort rows by argmax-column index, break ties by larger max value first.
    Col order: after rows reordered, each column placed by the first row (top-down)
    where it is that row's argmax; columns never argmax go to the end, sorted by max.
    """
    col_names = list(mat.columns)
    col_to_pos = {c: i for i, c in enumerate(col_names)}

    row_argmax = mat.idxmax(axis=1).map(col_to_pos)
    row_max = mat.max(axis=1)
    row_order = (
        pd.DataFrame({"argmax": row_argmax, "neg_max": -row_max}, index=mat.index)
        .sort_values(["argmax", "neg_max"])
        .index.tolist()
    )

    mat_r = mat.loc[row_order]
    col_first: dict = {}
    for ridx, cname in enumerate(mat_r.idxmax(axis=1).values):
        col_first.setdefault(cname, ridx)

    col_max = mat.max(axis=0)
    cols_used = sorted([c for c in col_names if c in col_first], key=lambda c: col_first[c])
    cols_unused = sorted(
        [c for c in col_names if c not in col_first], key=lambda c: -float(col_max[c])
    )
    return row_order, cols_used + cols_unused


def emit_native_label_f1(df: pd.DataFrame, cfg: dict, res_cols: List[str], out_dir: Path) -> None:
    studies = cfg["pipeline"]["native_labels"]["studies"]
    prefix = cfg["pipeline"]["native_labels"]["column_prefix"]
    log.info("Emitting native-label F1 heatmaps (%d resolutions × %d studies)", len(res_cols), len(studies))

    for res_col in res_cols:
        # Per-study F1 matrices; row labels are "{study}|{native_label}"
        frames = []
        for study in studies:
            native_col = f"{prefix}{study}"
            if native_col not in df.columns:
                continue
            m = _compute_f1_matrix(df, native_col, res_col)
            if m.empty:
                continue
            m.index = [f"{study}|{lbl}" for lbl in m.index]
            frames.append(m)

        if not frames:
            log.warning("  %s: no native-label frames, skipping", res_col)
            continue

        all_clusters = sorted(set().union(*[f.columns for f in frames]), key=lambda s: (len(s), s))
        stacked = pd.concat([f.reindex(columns=all_clusters, fill_value=0.0) for f in frames])

        row_order, col_order = _diagonal_ordering(stacked)
        stacked_ord = stacked.loc[row_order, col_order]

        csv_path = out_dir / f"native_label_f1_{res_col}.csv"
        stacked_ord.to_csv(csv_path)
        log.info("  wrote %s", csv_path)

        # Study color per row
        row_colors = pd.Series(
            [STUDY_COLORS.get(r.split("|", 1)[0], "#888888") for r in stacked_ord.index],
            index=stacked_ord.index,
            name="study",
        )

        # Sizing informed by actual immune F1 dimensions (128 rows fixed; cols 6-100)
        # Tall aspect is unavoidable with 128 labels; balance the axes per resolution.
        height = max(4.0, 0.11 * len(stacked_ord) + 1.0)   # ~14 in for 128 rows
        width = max(8.0, 0.24 * len(stacked_ord.columns) + 4.0)  # 4 in padding for y-labels, cbar, legend

        g = sns.clustermap(
            stacked_ord,
            row_cluster=False,
            col_cluster=False,
            row_colors=row_colors,
            cmap="viridis",
            vmin=0,
            vmax=1,
            figsize=(width, height),
            cbar_kws={"label": "F1"},
            cbar_pos=(0.95, 0.30, 0.010, 0.40),  # right side, very narrow vertical
            linewidths=0,
            xticklabels=True,
            yticklabels=True,
            dendrogram_ratio=0.01,
        )
        # Y-tick labels on the LEFT (clustermap default puts them on the right)
        g.ax_heatmap.yaxis.tick_left()
        g.ax_heatmap.yaxis.set_label_position("left")

        g.ax_heatmap.set_xlabel(f"cluster ({res_col})", fontsize=7)
        g.ax_heatmap.set_ylabel("")
        g.ax_heatmap.tick_params(axis="x", labelsize=5, rotation=90)
        g.ax_heatmap.tick_params(axis="y", labelsize=4)
        g.ax_heatmap.set_title(
            f"Native label × {res_col} (diagonal-ordered; rows color-coded by study)",
            fontsize=7,
        )

        from matplotlib.patches import Patch  # local import to avoid top-level

        used_studies = {r.split("|", 1)[0] for r in stacked_ord.index}
        handles = [
            Patch(facecolor=STUDY_COLORS.get(s, "#888888"), label=s)
            for s in studies
            if s in used_studies
        ]
        # With y-labels now on the left, the right side is free — place legend there
        # in figure coordinates (clustermap layout doesn't cooperate with ax-relative
        # bbox reliably; figure-relative keeps it outside the heatmap).
        g.fig.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(0.99, 0.5),
            fontsize=6,
            title="study",
            title_fontsize=6,
            frameon=False,
        )

        pdf_path = out_dir / f"native_label_f1_{res_col}.pdf"
        g.savefig(pdf_path, dpi=300, bbox_inches="tight")
        plt.close(g.fig)
        log.info("  wrote %s", pdf_path)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()
    structural = args.out_dir / "structural"
    structural.mkdir(parents=True, exist_ok=True)

    df, cfg, res_cols = load_inputs(args)

    emit_study_composition(df, structural)
    emit_leiden_umap_grid(df, res_cols, structural)
    emit_native_label_umap_grid(df, cfg, structural)
    emit_native_label_f1(df, cfg, res_cols, structural)

    log.info("Structural package complete: %s", structural)


if __name__ == "__main__":
    main()
