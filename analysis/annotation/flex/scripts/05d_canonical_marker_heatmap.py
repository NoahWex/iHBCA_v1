"""05d_canonical_marker_heatmap.py — V1 yaml canonical+identity markers × cluster heatmap.

For each resolution, render a heatmap of cluster mean log1p expression for V1 yaml
canonical_markers + identity_markers (rows grouped by V1 label, columns = clusters).

This is the dot-product to the existing top_marker_heatmap (which uses limma de novo
markers). The canonical heatmap shows whether the V1 label markers are actually
expressed where UCell says they are.

Outputs to {out_dir}/structural/:
    canonical_marker_heatmap_leiden_{res}.{pdf,csv}      one per resolution

Inputs:
    --bundle-dir            scvi_n100 bundle (cells.tsv, genes.tsv, counts.mtx.gz)
    --clusters-csv          per-cell cluster IDs (multi-resolution)
    --yaml                  V1 annotation yaml
    --resolutions           space-separated resolutions to render
    --out-dir               compartment output dir (e.g. outputs/Immune/)
    --max-markers-per-label max canonical+identity markers per label (default 8)

Style: no on-plot titles, viridis cmap, 7pt body text, 600 DPI. Follows the
cluster-mean-expression and V1-label marker-grouping logic of the canonical
marker heatmap convention used elsewhere in the publication tree.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("canonical_heatmap")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--bundle-dir", type=Path, required=True)
    p.add_argument("--clusters-csv", type=Path, required=True)
    p.add_argument("--yaml", type=Path, required=True)
    p.add_argument("--resolutions", nargs="+", required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--assignments-dir", type=Path, required=True,
                   help="Dir containing Immune_leiden_{res}_assignments.csv (UCell argmax per cluster)")
    p.add_argument("--compartment", required=True,
                   help="Used for the assignments file prefix (e.g. Immune)")
    p.add_argument("--max-markers-per-label", type=int, default=8)
    p.add_argument("--target-sum", type=float, default=1e4)
    return p.parse_args()


def load_mtx_gz(path: Path) -> sp.csr_matrix:
    with gzip.open(path, "rb") as fh:
        m = sio.mmread(fh)
    if not sp.issparse(m):
        m = sp.csr_matrix(m)
    return m.tocsr()


def normalize_log1p(X: sp.csr_matrix, target_sum: float) -> sp.csr_matrix:
    sums = np.asarray(X.sum(axis=1)).ravel()
    sums[sums == 0] = 1.0
    scale = target_sum / sums
    Xn = X.multiply(scale[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def cluster_means(Xn: sp.csr_matrix, cluster_ids: pd.Series, gene_idx: List[int]) -> pd.DataFrame:
    sub = Xn[:, gene_idx]
    df = pd.DataFrame(sub.toarray(), index=cluster_ids.index)
    df["__cluster"] = cluster_ids.values
    means = df.groupby("__cluster").mean()
    means.columns = [str(g) for g in means.columns]
    return means


def collect_label_markers(yaml_path: Path, max_per_label: int) -> Dict[str, List[str]]:
    with open(yaml_path) as fh:
        cfg = yaml.safe_load(fh)
    markers: Dict[str, List[str]] = {}
    for label, body in (cfg.get("labels") or {}).items():
        if not isinstance(body, dict):
            continue
        canon = body.get("canonical_markers") or []
        ident = body.get("identity_markers") or []
        seen, ordered = set(), []
        for g in list(canon) + list(ident):
            if g and g not in seen:
                seen.add(g)
                ordered.append(g)
            if len(ordered) >= max_per_label:
                break
        if ordered:
            markers[label] = ordered
    return markers


def _column_block_order(label_to_markers: Dict[str, List[str]], label_order: List[str]) -> List[tuple]:
    """Place each canonical marker under its first declaring V1 label.

    Returns list of (gene, label) tuples in column order — true per-label blocks
    even when markers are shared across labels (matches V1 publication pattern).
    """
    placed: Dict[str, str] = {}
    ordered: List[tuple] = []
    for label in label_order:
        for g in label_to_markers.get(label, []):
            if g not in placed:
                placed[g] = label
                ordered.append((g, label))
    return ordered


def render_heatmap(
    means: pd.DataFrame,
    label_to_markers: Dict[str, List[str]],
    label_order: List[str],
    cluster_order: List[int],
    cluster_argmax_labels: Dict[int, str],
    title_panel: str,
    out_pdf: Path,
    out_csv: Path,
) -> None:
    block_order = _column_block_order(label_to_markers, label_order)
    block_order = [(g, lbl) for (g, lbl) in block_order if g in means.columns]
    if not block_order:
        log.warning("[%s] no markers in expression data", title_panel)
        return
    genes = [g for g, _ in block_order]
    col_groups = [lbl for _, lbl in block_order]

    M = means.loc[cluster_order, genes].copy()
    M.index = [f"c{int(c)}" for c in cluster_order]
    Z = (M - M.mean(axis=0)) / M.std(axis=0).replace(0, 1)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    Z.to_csv(out_csv, index_label="cluster")

    n_rows, n_cols = Z.shape
    height = max(3.5, 0.16 * n_rows + 1.2)
    width = max(8.0, 0.13 * n_cols + 3.5)

    fig, ax = plt.subplots(figsize=(width, height), dpi=150)
    import seaborn as sns

    sns.heatmap(
        Z,
        cmap="RdBu_r",
        center=0,
        vmin=-2,
        vmax=2,
        xticklabels=genes,
        yticklabels=[f"{c}  [{cluster_argmax_labels.get(int(c[1:]), '?')}]" for c in Z.index],
        cbar_kws={"label": "z-score", "shrink": 0.3},
        linewidths=0,
        ax=ax,
    )
    ax.set_xlabel("V1 canonical + identity markers (block-ordered by first declaring label)", fontsize=7)
    ax.set_ylabel("cluster  [UCell argmax]", fontsize=7)
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=5)

    boundaries: List[int] = []
    last = None
    for i, lbl in enumerate(col_groups):
        if lbl != last:
            boundaries.append(i)
            last = lbl
    boundaries.append(n_cols)
    for b in boundaries[1:-1]:
        ax.axvline(b, color="black", linewidth=0.4)

    secax = ax.secondary_xaxis("top")
    centers = [(boundaries[k] + boundaries[k + 1]) / 2 for k in range(len(boundaries) - 1)]
    block_labels = [col_groups[boundaries[k]] for k in range(len(boundaries) - 1)]
    secax.set_xticks(centers)
    secax.set_xticklabels(block_labels, fontsize=5, rotation=45, ha="left")
    secax.tick_params(axis="x", length=0)

    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("[%s] %d clusters × %d markers -> %s", title_panel, n_rows, n_cols, out_pdf.name)


def main() -> None:
    args = parse_args()

    cells_tsv = args.bundle_dir / "cells.tsv"
    genes_tsv = args.bundle_dir / "genes.tsv"
    mtx = args.bundle_dir / "counts.mtx.gz"
    for f in (cells_tsv, genes_tsv, mtx, args.yaml, args.clusters_csv):
        if not f.exists():
            sys.exit(f"ERROR: missing input: {f}")

    cells = pd.read_csv(cells_tsv, sep="\t", header=None, names=["cell_id"])
    genes = pd.read_csv(genes_tsv, sep="\t", header=None, names=["gene"])
    log.info("Cells=%d, Genes=%d", len(cells), len(genes))

    X = load_mtx_gz(mtx)
    if X.shape == (len(cells), len(genes)):
        log.info("mtx orientation: cells × genes")
    elif X.shape == (len(genes), len(cells)):
        log.info("mtx orientation: genes × cells (transposing)")
        X = X.T.tocsr()
    else:
        sys.exit(f"ERROR: mtx shape {X.shape}")

    log.info("Normalizing (target_sum=%s, log1p)", args.target_sum)
    Xn = normalize_log1p(X, args.target_sum)

    gene_index = pd.Series(np.arange(len(genes)), index=genes["gene"].values)
    if gene_index.index.has_duplicates:
        gene_index = gene_index[~gene_index.index.duplicated(keep="first")]

    label_to_markers = collect_label_markers(args.yaml, args.max_markers_per_label)
    label_order = list(label_to_markers.keys())
    all_markers = sorted(
        {g for gs in label_to_markers.values() for g in gs if g in gene_index.index}
    )
    gene_idx = [int(gene_index[g]) for g in all_markers]
    log.info("Markers: %d unique across %d labels", len(all_markers), len(label_to_markers))

    clusters = pd.read_csv(args.clusters_csv)

    out_root = args.out_dir / "structural"
    out_root.mkdir(parents=True, exist_ok=True)

    for res in args.resolutions:
        col = f"leiden_{res}"
        if col not in clusters.columns:
            log.warning("Skip res=%s (column %s missing)", res, col)
            continue
        cluster_series = clusters[col].astype(int)
        cluster_series.index = clusters["cell_id"] if "cell_id" in clusters else clusters.index
        if len(cluster_series) != len(cells):
            cluster_series = cluster_series.reindex(cells["cell_id"]).reset_index(drop=True)
        else:
            cluster_series = pd.Series(cluster_series.values, index=range(len(cells)))

        means_full = cluster_means(Xn, cluster_series, gene_idx)
        means_full.columns = all_markers

        cluster_order = sorted(means_full.index.astype(int).tolist())

        cluster_argmax_labels: Dict[int, str] = {}
        assign_csv = args.assignments_dir / f"{args.compartment}_leiden_{res}_assignments.csv"
        if assign_csv.exists():
            adf = pd.read_csv(assign_csv)
            for _, row in adf.iterrows():
                cluster_argmax_labels[int(row["cluster"])] = str(row["L2S_assignment"])
        else:
            log.warning("Assignments CSV missing: %s (row labels will show '?')", assign_csv)

        out_pdf = out_root / f"canonical_marker_heatmap_leiden_{res}.pdf"
        out_csv = out_root / f"canonical_marker_heatmap_leiden_{res}.csv"
        render_heatmap(
            means_full,
            label_to_markers,
            label_order,
            cluster_order,
            cluster_argmax_labels,
            f"leiden_{res}",
            out_pdf,
            out_csv,
        )

    log.info("Done.")


if __name__ == "__main__":
    main()
