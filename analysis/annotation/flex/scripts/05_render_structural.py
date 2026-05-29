"""05_render_structural.py — pre-review structural diagnostic package per compartment.

FLEX is single-study, so `study_composition` is replaced with `patient_composition`
(4 patients). FLEX has no per-study native author labels, so the F1 metric uses
`v1_label_f1` — F1 of UCell-positive cells (combined score >= threshold) for each
V1 L2 label vs cluster membership at each resolution. This answers "at which
resolution is each V1 label best resolved by some cluster?"

Outputs to {out_dir}/structural/:
    patient_composition.{csv,pdf}
    leiden_umap_grid.pdf
    v1_label_umap_grid.pdf                   per-V1-label UMAP highlight (UCell-positive)
    v1_label_f1_leiden_{res}.{csv,pdf}       F1 heatmap (V1 label × cluster)
    ucell_score_heatmap_leiden_{res}.pdf     cluster × V1 label combined-score heatmap
    top_marker_heatmap_leiden_{res}.{csv,pdf}  limma top markers × cluster heatmap
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
log = logging.getLogger("render_structural")


# ---------- helpers (pattern: V1 v2 00_build_structural_package.py) -----------

def categorical_palette(n: int) -> List[tuple]:
    base = (
        list(colormaps["tab20"].colors)
        + list(colormaps["tab20b"].colors)
        + list(colormaps["tab20c"].colors)
    )
    if n <= len(base):
        return base[:n]
    hsv = colormaps["hsv"]
    return [hsv(i / n) for i in range(n)]


def diagonal_ordering(mat: pd.DataFrame) -> tuple[list, list]:
    col_names = list(mat.columns)
    col_to_pos = {c: i for i, c in enumerate(col_names)}
    row_argmax = mat.idxmax(axis=1).map(col_to_pos)
    row_max = mat.max(axis=1)
    row_order = (
        pd.DataFrame({"argmax": row_argmax, "neg_max": -row_max}, index=mat.index)
        .sort_values(by=["argmax", "neg_max"])
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


def f1_matrix_from_positives(
    df: pd.DataFrame, positive_col: str, cluster_col: str
) -> pd.DataFrame:
    """F1 matrix: rows = positive_col values (binary {0,1}), cols = cluster.
    Returns single-row DataFrame (positive_col=1 vs cluster).

    For single-binary-label F1 we collapse positives into one row per V1 label;
    callers stack rows for many labels. See `emit_v1_label_f1`.
    """
    raise NotImplementedError("Use compute_v1_label_f1_stack for stacked F1.")


def compute_v1_label_f1_stack(
    df: pd.DataFrame,
    label_cols: list[str],
    cluster_col: str,
    threshold: float,
) -> pd.DataFrame:
    """For each V1 label column (combined UCell score), threshold to positives,
    then compute F1 of each cluster as predictor. Returns DataFrame indexed by
    V1 label, columns = cluster, values = F1.
    """
    cluster_vals = df[cluster_col].astype(str).values
    cluster_unique = np.unique(cluster_vals)
    rows = {}
    for lb in label_cols:
        scores = df[lb].values
        pos = scores >= threshold
        n_pos = int(pos.sum())
        if n_pos == 0:
            continue
        f1_per_cluster = {}
        for cl in cluster_unique:
            in_cluster = cluster_vals == cl
            tp = int((in_cluster & pos).sum())
            pred_pos = int(in_cluster.sum())
            if pred_pos == 0:
                f1_per_cluster[cl] = 0.0
                continue
            precision = tp / pred_pos
            recall = tp / n_pos
            denom = precision + recall
            f1 = (2 * precision * recall / denom) if denom > 0 else 0.0
            f1_per_cluster[cl] = f1
        rows[lb] = f1_per_cluster
    return pd.DataFrame(rows).T.fillna(0.0)


# ---------- CLI / loaders -----------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compartment", required=True, choices=["Epithelial", "Immune", "Stromal"])
    p.add_argument("--clusters-csv", required=True, type=Path)
    p.add_argument("--umap-csv", required=True, type=Path,
                   help="scvi_n100/umap.csv (cell_id + UMAP_1 + UMAP_2 or X_umap dims)")
    p.add_argument("--ucell-per-cell-csv", required=True, type=Path)
    p.add_argument("--limma-dir", required=True, type=Path,
                   help="dir with limma_top_markers_{Compartment}_leiden_{res}.csv")
    p.add_argument("--gene-data", type=Path, default=None,
                   help="optional gene_data.csv with ENSG → symbol mapping")
    p.add_argument("--resolutions", nargs="+", required=True)
    p.add_argument("--weight-identity", type=float, default=2.0)
    p.add_argument("--weight-canonical", type=float, default=1.0)
    p.add_argument("--positive-threshold", type=float, default=0.5)
    p.add_argument("--top-markers-per-cluster", type=int, default=15)
    p.add_argument("--out-dir", required=True, type=Path)
    return p.parse_args()


def load_umap(umap_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(umap_csv)
    if "cell_id" not in df.columns:
        df = df.rename(columns={df.columns[0]: "cell_id"})
    # Allow either UMAP_1/UMAP_2 or first two non-cell_id columns
    other_cols = [c for c in df.columns if c != "cell_id"]
    if "UMAP_1" not in df.columns:
        df = df.rename(columns={other_cols[0]: "UMAP_1", other_cols[1]: "UMAP_2"})
    return df[["cell_id", "UMAP_1", "UMAP_2"]]


def load_ucell_combined(per_cell_csv: Path, w_identity: float, w_canonical: float) -> pd.DataFrame:
    per_cell = pd.read_csv(per_cell_csv)
    sig_cols = [c for c in per_cell.columns if "__" in c]
    labels = sorted(set(c.rsplit("__", 1)[0] for c in sig_cols))
    w_total = w_identity + w_canonical
    out = pd.DataFrame({"cell_id": per_cell["cell_id"]})
    for lb in labels:
        canon = per_cell.get(f"{lb}__canonical")
        ident = per_cell.get(f"{lb}__identity")
        if canon is None: canon = pd.Series(0.0, index=per_cell.index)
        if ident is None: ident = pd.Series(0.0, index=per_cell.index)
        out[lb] = (ident * w_identity + canon * w_canonical) / w_total
    return out


# ---------- output emitters ---------------------------------------------------

def emit_patient_composition(df: pd.DataFrame, out_dir: Path) -> None:
    if "patient_id" not in df.columns:
        log.warning("patient_id missing; skipping patient_composition")
        return
    log.info("Emitting patient_composition")
    tab = df.groupby("patient_id", observed=True).size().reset_index(name="n_cells")
    csv_path = out_dir / "patient_composition.csv"
    tab.to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.bar(tab["patient_id"], tab["n_cells"], color="#4a7ba6", edgecolor="black", linewidth=0.5)
    for i, n in enumerate(tab["n_cells"]):
        ax.text(i, n + tab["n_cells"].max() * 0.02, f"{int(n):,}",
                ha="center", va="bottom", fontsize=6)
    ax.set_xlabel("patient_id", fontsize=7)
    ax.set_ylabel("cells", fontsize=7)
    ax.tick_params(labelsize=6)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_dir / "patient_composition.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_dir / "patient_composition.pdf")


def emit_leiden_umap_grid(df: pd.DataFrame, res_cols: list[str], out_dir: Path) -> None:
    log.info("Emitting leiden_umap_grid")
    n = len(res_cols)
    ncols = min(5, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 2.5, nrows * 2.5))
    axes = np.atleast_2d(axes).flatten()
    for i, col in enumerate(res_cols):
        ax = axes[i]
        clusters = sorted(df[col].astype(str).unique(), key=lambda s: (len(s), s))
        palette = categorical_palette(len(clusters))
        cm = dict(zip(clusters, palette))
        colors = df[col].astype(str).map(cm).values
        ax.scatter(df["UMAP_1"], df["UMAP_2"], c=colors,
                   s=0.5, alpha=0.6, rasterized=True, linewidths=0)
        ax.set_title(col.replace("leiden_", "res="), fontsize=7)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_linewidth(0.3)
        ax.text(0.02, 0.98, f"k={len(clusters)}", transform=ax.transAxes,
                fontsize=5, va="top", ha="left",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
    for j in range(len(res_cols), len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_dir / "leiden_umap_grid.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_dir / "leiden_umap_grid.pdf")


def emit_v1_label_umap_grid(
    df: pd.DataFrame, label_cols: list[str], threshold: float, out_dir: Path
) -> None:
    log.info("Emitting v1_label_umap_grid (%d labels, threshold=%.2f)",
             len(label_cols), threshold)
    n = len(label_cols)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.0, nrows * 3.0))
    axes = np.atleast_2d(axes).flatten()
    for i, lb in enumerate(label_cols):
        ax = axes[i]
        ax.scatter(df["UMAP_1"], df["UMAP_2"], c="#d8d8d8",
                   s=0.3, alpha=0.4, rasterized=True, linewidths=0)
        pos = df[lb].values >= threshold
        ax.scatter(df.loc[pos, "UMAP_1"], df.loc[pos, "UMAP_2"],
                   c="#d62728", s=0.6, alpha=0.8, rasterized=True, linewidths=0)
        ax.set_title(f"{lb}  n+={int(pos.sum()):,}", fontsize=6)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_linewidth(0.3)
    for j in range(len(label_cols), len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_dir / "v1_label_umap_grid.pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_dir / "v1_label_umap_grid.pdf")


def emit_f1_heatmaps(
    df: pd.DataFrame,
    label_cols: list[str],
    res_cols: list[str],
    threshold: float,
    out_dir: Path,
) -> None:
    log.info("Emitting v1_label_f1_leiden_<res> heatmaps (%d resolutions × %d labels)",
             len(res_cols), len(label_cols))
    for res_col in res_cols:
        f1 = compute_v1_label_f1_stack(df, label_cols, res_col, threshold)
        if f1.empty:
            log.warning("  %s: empty F1 matrix; skipping", res_col)
            continue
        all_clusters = sorted(f1.columns, key=lambda s: (len(s), s))
        f1 = f1.reindex(columns=all_clusters, fill_value=0.0)
        row_order, col_order = diagonal_ordering(f1)
        f1_ord = f1.loc[row_order, col_order]

        csv_path = out_dir / f"v1_label_f1_{res_col}.csv"
        f1_ord.to_csv(csv_path)

        h = max(3.0, 0.18 * len(f1_ord) + 1.0)
        w = max(6.0, 0.20 * len(f1_ord.columns) + 2.0)
        fig, ax = plt.subplots(figsize=(w, h))
        sns.heatmap(f1_ord, cmap="Blues", vmin=0, vmax=1,
                    cbar_kws={"label": "F1", "shrink": 0.4},
                    linewidths=0, ax=ax,
                    xticklabels=True, yticklabels=True)
        ax.set_xlabel(f"cluster ({res_col})", fontsize=7)
        ax.set_ylabel("V1 label (UCell-positive)", fontsize=7)
        ax.tick_params(axis="x", labelsize=5, rotation=90)
        ax.tick_params(axis="y", labelsize=5)
        plt.tight_layout()
        pdf_path = out_dir / f"v1_label_f1_{res_col}.pdf"
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        log.info("  wrote %s + %s", csv_path, pdf_path)


def emit_ucell_score_heatmaps(
    df: pd.DataFrame, label_cols: list[str], res_cols: list[str], out_dir: Path
) -> None:
    log.info("Emitting ucell_score_heatmap_leiden_<res> (cluster × V1 label combined-score)")
    for res_col in res_cols:
        cluster_vals = df[res_col].astype(str).values
        clusters = sorted(np.unique(cluster_vals), key=lambda s: (len(s), s))
        mat = pd.DataFrame(0.0, index=clusters, columns=label_cols)
        for cl in clusters:
            mask = cluster_vals == cl
            if not mask.any():
                continue
            mat.loc[cl] = df.loc[mask, label_cols].median(axis=0).values
        row_order, col_order = diagonal_ordering(mat.T)  # transpose → labels × clusters for ordering
        mat_ord = mat.loc[col_order, row_order]
        h = max(3.0, 0.18 * len(mat_ord) + 1.0)
        w = max(6.0, 0.20 * len(mat_ord.columns) + 2.0)
        fig, ax = plt.subplots(figsize=(w, h))
        sns.heatmap(mat_ord, cmap="viridis",
                    cbar_kws={"label": "median UCell combined", "shrink": 0.4},
                    linewidths=0, ax=ax,
                    xticklabels=True, yticklabels=True)
        ax.set_xlabel("V1 label", fontsize=7)
        ax.set_ylabel(f"cluster ({res_col})", fontsize=7)
        ax.tick_params(axis="x", labelsize=5, rotation=90)
        ax.tick_params(axis="y", labelsize=5)
        plt.tight_layout()
        pdf_path = out_dir / f"ucell_score_heatmap_{res_col}.pdf"
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        log.info("  wrote %s", pdf_path)


def emit_top_marker_heatmaps(
    limma_dir: Path,
    res_cols: list[str],
    compartment: str,
    n_top: int,
    gene_data: Path | None,
    out_dir: Path,
) -> None:
    log.info("Emitting top_marker_heatmap_leiden_<res> (%d resolutions)", len(res_cols))
    ensg_to_sym = {}
    if gene_data is not None and gene_data.exists():
        gd = pd.read_csv(gene_data, index_col=0)
        if "gene_symbol" in gd.columns:
            ensg_to_sym = gd["gene_symbol"].to_dict()
    for res_col in res_cols:
        res = res_col.replace("leiden_", "")
        top_csv = limma_dir / f"limma_top_markers_{compartment}_leiden_{res}.csv"
        if not top_csv.exists():
            log.warning("  missing %s — skipping", top_csv)
            continue
        markers = pd.read_csv(top_csv)
        if "cluster" not in markers.columns:
            log.warning("  %s missing cluster column — skipping", top_csv)
            continue
        # Top-N per cluster, sorted by padj asc
        markers["cluster"] = markers["cluster"].astype(str)
        top = (markers.sort_values(by=["cluster", "padj"], ascending=[True, True])
                      .groupby("cluster", group_keys=False)
                      .head(n_top))
        # Wide: rows = cluster, cols = gene, values = log2FoldChange (use NA elsewhere)
        wide = top.pivot_table(index="cluster", columns="gene",
                               values="log2FoldChange", aggfunc="first")
        # Order: cluster numeric, gene by first-appearing cluster
        gene_first_cluster = (top.sort_values(by=["cluster"])
                                  .drop_duplicates("gene")
                                  .set_index("gene")["cluster"])
        gene_order = gene_first_cluster.index.tolist()
        cluster_order = sorted(wide.index, key=lambda s: (len(s), s))
        wide = wide.reindex(index=cluster_order, columns=gene_order)
        csv_path = out_dir / f"top_marker_heatmap_{res_col}.csv"
        wide.to_csv(csv_path)
        # Render
        gene_labels = [ensg_to_sym.get(g, g) for g in gene_order]
        h = max(3.0, 0.18 * len(cluster_order) + 1.0)
        w = max(8.0, 0.14 * len(gene_order) + 2.5)
        fig, ax = plt.subplots(figsize=(w, h))
        sns.heatmap(wide.fillna(0.0), cmap="RdBu_r", center=0, vmin=-3, vmax=3,
                    xticklabels=gene_labels, yticklabels=True,
                    cbar_kws={"label": "log2FC vs rest", "shrink": 0.4},
                    linewidths=0, ax=ax)
        ax.set_xlabel(f"top {n_top} markers per cluster ({len(gene_order)} unique)", fontsize=7)
        ax.set_ylabel(f"cluster ({res_col})", fontsize=7)
        ax.tick_params(axis="x", labelsize=5, rotation=90)
        ax.tick_params(axis="y", labelsize=5)
        plt.tight_layout()
        pdf_path = out_dir / f"top_marker_heatmap_{res_col}.pdf"
        fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        log.info("  wrote %s + %s", csv_path, pdf_path)


# ---------- main --------------------------------------------------------------

def main() -> None:
    args = parse_args()
    out_dir = args.out_dir / "structural"
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading clusters: %s", args.clusters_csv)
    clusters = pd.read_csv(args.clusters_csv)
    res_cols = [f"leiden_{r}" for r in args.resolutions]
    missing = [c for c in res_cols if c not in clusters.columns]
    if missing:
        sys.exit(f"ERROR: missing leiden columns {missing}")

    log.info("Loading UMAP: %s", args.umap_csv)
    umap = load_umap(args.umap_csv)

    log.info("Loading per-cell UCell: %s", args.ucell_per_cell_csv)
    combined = load_ucell_combined(
        args.ucell_per_cell_csv, args.weight_identity, args.weight_canonical
    )
    label_cols = [c for c in combined.columns if c != "cell_id"]
    log.info("  %d V1 labels", len(label_cols))

    df = clusters.merge(umap, on="cell_id", how="inner")
    df = df.merge(combined, on="cell_id", how="inner")
    log.info("Joined: %d cells, %d cols", len(df), df.shape[1])

    emit_patient_composition(df, out_dir)
    emit_leiden_umap_grid(df, res_cols, out_dir)
    emit_v1_label_umap_grid(df, label_cols, args.positive_threshold, out_dir)
    emit_f1_heatmaps(df, label_cols, res_cols, args.positive_threshold, out_dir)
    emit_ucell_score_heatmaps(df, label_cols, res_cols, out_dir)
    emit_top_marker_heatmaps(
        args.limma_dir, res_cols, args.compartment,
        args.top_markers_per_cluster, args.gene_data, out_dir
    )

    log.info("Structural package complete: %s", out_dir)


if __name__ == "__main__":
    main()
