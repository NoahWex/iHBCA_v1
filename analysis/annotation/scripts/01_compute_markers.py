"""
01_compute_markers.py — top-marker heatmaps per Leiden resolution.

Orchestrator. Reads existing limma outputs from `markers/limma/{comp}/leiden_{res}/`;
does NOT regenerate limma (regeneration requires R; deferred to a separate script
if needed).

For each resolution with existing limma markers:
  1. Pick top N markers per cluster targeting 75-90 unique genes total
  2. Decode ENSEMBL IDs to gene symbols via gene_data.csv
  3. Compute cluster-mean log1p expression for selected genes from counts.npz
  4. Render z-scored heatmap (genes × clusters) as PDF
  5. Emit selected-markers CSV with cluster, rank, ensembl, symbol, log2FC, padj

Outputs under out_dir/structural/:
  top_marker_heatmap_res{X}.pdf
  top_markers_selected_res{X}.csv

Resolutions without existing limma outputs are logged as WARN and skipped.

Design doc: iHBCA_publication/coordination/plans/annotation_pipeline_v2.md §4.1
Pattern-ref: 00_build_structural_package.py (same container, same arg style).
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
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from scipy import sparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("compute_markers")


# =============================================================================
# CLI
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--compartment", required=True, choices=["imm", "epi", "str"])
    p.add_argument("--leiden-csv", required=True, type=Path)
    p.add_argument("--counts-npz", required=True, type=Path)
    p.add_argument("--gene-data", required=True, type=Path)
    p.add_argument("--limma-dir", required=True, type=Path, help="markers/limma/{comp}/")
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--target-n-markers", type=int, default=85, help="target unique markers per heatmap")
    p.add_argument("--max-per-cluster", type=int, default=5)
    p.add_argument("--padj-threshold", type=float, default=0.01)
    return p.parse_args()


# =============================================================================
# Loading
# =============================================================================

def load_ensembl_to_symbol(gene_data_path: Path) -> dict:
    log.info("Loading gene_data: %s", gene_data_path)
    gd = pd.read_csv(gene_data_path, index_col=0)
    # Column layout verified 2026-04-24: index = ENSG; gene_symbol column has symbol
    if "gene_symbol" not in gd.columns:
        raise RuntimeError(f"gene_data missing gene_symbol column; cols: {list(gd.columns)}")
    mapping = gd["gene_symbol"].to_dict()
    log.info("  %d ENSG -> symbol mappings", len(mapping))
    return mapping


def load_counts_and_cells(counts_npz: Path, leiden_csv: Path) -> tuple[sparse.csr_matrix, pd.DataFrame]:
    log.info("Loading counts: %s", counts_npz)
    loaded = np.load(counts_npz, allow_pickle=True)
    # Conventional scipy sparse .npz: data, indices, indptr, shape
    if "data" in loaded.files:
        X = sparse.csr_matrix(
            (loaded["data"], loaded["indices"], loaded["indptr"]),
            shape=tuple(loaded["shape"]),
        )
    else:
        raise RuntimeError(f"counts npz has unexpected keys: {loaded.files}")
    log.info("  counts shape: %s (cells × genes)", X.shape)

    log.info("Loading leiden: %s", leiden_csv)
    leiden = pd.read_csv(leiden_csv)
    log.info("  %d cells, resolutions: %s", len(leiden), [c for c in leiden.columns if c.startswith("leiden_")])
    return X, leiden


# =============================================================================
# Marker selection
# =============================================================================

def select_top_markers(
    markers: pd.DataFrame,
    target_n: int,
    max_per_cluster: int,
    padj_threshold: float,
) -> pd.DataFrame:
    """Select top markers per cluster aiming for target_n unique genes, evenly distributed.

    Iteratively increase per-cluster depth until unique-marker count reaches target.
    Sorted within cluster by log2FoldChange descending.
    """
    if "cell_type" in markers.columns:
        markers = markers.rename(columns={"cell_type": "cluster"})
    required = {"cluster", "gene", "log2FoldChange", "padj"}
    missing = required - set(markers.columns)
    if missing:
        raise RuntimeError(f"markers missing columns: {missing}; present: {list(markers.columns)}")

    sig = markers[markers["padj"] < padj_threshold].copy()
    sig = sig.sort_values(["cluster", "log2FoldChange"], ascending=[True, False])
    clusters = sig["cluster"].unique()

    if max_per_cluster < 1:
        raise ValueError("max_per_cluster must be >= 1")

    picks = sig.groupby("cluster", group_keys=False).head(1)
    depth = 1
    for d in range(1, max_per_cluster + 1):
        picks = sig.groupby("cluster", group_keys=False).head(d)
        depth = d
        if picks["gene"].nunique() >= target_n:
            break
    picks = picks.assign(rank=picks.groupby("cluster").cumcount() + 1)
    picks = picks[["cluster", "rank", "gene", "log2FoldChange", "padj"]].reset_index(drop=True)
    log.info("  selected %d rows (%d unique genes) across %d clusters at depth=%d", len(picks), picks["gene"].nunique(), len(clusters), depth)
    return picks


# =============================================================================
# Expression aggregation
# =============================================================================

def cluster_mean_expression(
    X: sparse.csr_matrix, cluster_labels: np.ndarray, gene_indices: np.ndarray
) -> pd.DataFrame:
    """Return DataFrame (cluster × gene) of mean log1p-transformed expression.

    log1p applied element-wise BEFORE averaging (standard Seurat/scanpy norm).
    """
    clusters = np.sort(np.unique(cluster_labels))
    X_sub = X[:, gene_indices]
    # log1p on raw counts — normalized per cell would require size factors; we use
    # log1p of raw since that matches scanpy's log1p-normalize pipeline when
    # counts are already normalized. If counts are raw, this is approximate but
    # adequate for cluster-relative z-score visualization.
    X_log = X_sub.copy()
    X_log.data = np.log1p(X_log.data)
    out = np.zeros((len(clusters), len(gene_indices)))
    for i, c in enumerate(clusters):
        mask = cluster_labels == c
        if mask.sum() == 0:
            continue
        out[i] = np.asarray(X_log[mask].mean(axis=0)).flatten()
    return pd.DataFrame(out, index=clusters.astype(str), columns=gene_indices)


# =============================================================================
# Heatmap rendering
# =============================================================================

def render_heatmap(
    expr: pd.DataFrame,
    ensg_to_symbol: dict,
    picks: pd.DataFrame,
    res_name: str,
    out_pdf: Path,
) -> None:
    # Rename columns (gene idx → ENSG) is done elsewhere; expr columns are ENSG
    # Order genes by their first-appearing cluster (for diagonal pattern)
    gene_first_cluster = picks.sort_values(["cluster", "rank"]).drop_duplicates("gene")[["gene", "cluster"]]
    gene_order = gene_first_cluster["gene"].tolist()
    expr_ordered = expr[gene_order]
    cluster_order = sorted(expr.index.astype(str).unique(), key=lambda s: (len(s), s))
    expr_ordered = expr_ordered.loc[cluster_order]

    # Z-score per gene across clusters
    z = (expr_ordered - expr_ordered.mean(axis=0)) / (expr_ordered.std(axis=0).replace(0, 1))

    # Symbol labels
    symbol_labels = [ensg_to_symbol.get(g, g) for g in gene_order]

    n_genes = len(gene_order)
    n_clusters = len(cluster_order)
    height = max(4.0, 0.15 * n_clusters + 1.5)
    width = max(8.0, 0.14 * n_genes + 2.5)

    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        z,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        xticklabels=symbol_labels,
        yticklabels=True,
        cbar_kws={"label": "z-score", "shrink": 0.4},
        linewidths=0,
        ax=ax,
    )
    ax.set_xlabel(f"gene ({n_genes} unique)", fontsize=7)
    ax.set_ylabel(f"cluster ({res_name})", fontsize=7)
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=5)
    ax.set_title(f"Top markers × clusters ({res_name}, z-scored)", fontsize=7)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()
    structural = args.out_dir / "structural"
    structural.mkdir(parents=True, exist_ok=True)

    ensg_to_symbol = load_ensembl_to_symbol(args.gene_data)
    X, leiden = load_counts_and_cells(args.counts_npz, args.leiden_csv)
    gene_data = pd.read_csv(args.gene_data, index_col=0)
    ensg_to_idx = {ensg: i for i, ensg in enumerate(gene_data.index)}

    res_cols = [c for c in leiden.columns if c.startswith("leiden_")]
    log.info("Processing %d resolutions", len(res_cols))

    for res_col in res_cols:
        res_suffix = res_col.replace("leiden_", "")
        limma_subdir = args.limma_dir / f"leiden_{res_suffix}"
        top_csv = limma_subdir / f"limma_top_markers_leiden_{res_suffix}.csv"
        if not top_csv.exists():
            log.warning("  %s: no limma output at %s, skipping", res_col, top_csv)
            continue

        log.info("Resolution %s", res_col)
        markers = pd.read_csv(top_csv)

        picks = select_top_markers(
            markers, args.target_n_markers, args.max_per_cluster, args.padj_threshold
        )
        picks_with_symbols = picks.assign(symbol=picks["gene"].map(ensg_to_symbol))
        picks_csv = structural / f"top_markers_selected_{res_col}.csv"
        picks_with_symbols.to_csv(picks_csv, index=False)
        log.info("  wrote %s", picks_csv)

        unique_genes = picks["gene"].unique().tolist()
        gene_indices = np.array([ensg_to_idx[g] for g in unique_genes if g in ensg_to_idx], dtype=int)
        missing = [g for g in unique_genes if g not in ensg_to_idx]
        if missing:
            log.warning("  %d selected genes not in gene_data index, dropped", len(missing))

        cluster_labels = leiden[res_col].astype(str).values
        expr = cluster_mean_expression(X, cluster_labels, gene_indices)
        expr.columns = [unique_genes[i] for i, idx in enumerate(gene_indices)]

        out_pdf = structural / f"top_marker_heatmap_{res_col}.pdf"
        render_heatmap(expr, ensg_to_symbol, picks, res_col, out_pdf)

    log.info("Marker heatmaps complete: %s", structural)


if __name__ == "__main__":
    main()
