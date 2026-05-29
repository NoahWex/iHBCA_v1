"""
03_render_annotation_package.py — post-YAML annotation evidence figures.

Produces §4.2 outputs under outputs/{comp}/annotation/:
  label_umap.pdf                 proposed-label UMAP with labels projected on
  heatmap_canonical.pdf          declared canonical markers per label
  heatmap_top.pdf                top DE markers per label (fresh one-vs-rest)
  heatmap_lineage_tiled.pdf      L1 lineage markers tiled across curated resolutions
  heatmap_artifact_tiled.pdf     artifact markers tiled across curated resolutions
  label_study_composition.pdf    per-label study barplot (validity evidence for multi-study)
  label_patient_composition.pdf  per-label patient barplot (validity evidence for multi-patient)

Inputs: v2 YAML, labels.csv (from 06), leiden CSV, UMAP CSV, counts.npz, gene_data,
metadata_enriched.csv, compartment_config.yaml

Design doc: iHBCA_publication/coordination/plans/annotation_pipeline_v2.md §4.2
Pattern-ref: 00_build_structural_package.py for plot helpers; 01_compute_markers.py
for marker expression aggregation.
"""

from __future__ import annotations

import argparse
import logging
import os
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
from scipy import sparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("render_annotation_package")


# =============================================================================
# CLI + Loading
# =============================================================================

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--v2-yaml", required=True, type=Path)
    p.add_argument("--labels-csv", required=True, type=Path)
    p.add_argument("--leiden-csv", required=True, type=Path)
    p.add_argument("--umap-csv", required=True, type=Path)
    p.add_argument("--counts-npz", required=True, type=Path)
    p.add_argument("--gene-data", required=True, type=Path)
    p.add_argument("--metadata-csv", required=True, type=Path)
    p.add_argument("--config", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--top-n-per-label", type=int, default=5)
    return p.parse_args()


def _categorical_palette(n: int) -> List[tuple]:
    base = list(colormaps["tab20"].colors) + list(colormaps["tab20b"].colors) + list(colormaps["tab20c"].colors)
    if n <= len(base):
        return base[:n]
    hsv = colormaps["hsv"]
    return [hsv(i / n) for i in range(n)]


# =============================================================================
# Cluster-mean expression helper (labels as groups)
# =============================================================================

def cluster_mean_log1p(X: sparse.csr_matrix, group_labels: np.ndarray, gene_indices: np.ndarray) -> pd.DataFrame:
    groups = np.sort(np.unique(group_labels[group_labels != None]))  # noqa: E711
    X_sub = X[:, gene_indices]
    X_log = X_sub.copy()
    X_log.data = np.log1p(X_log.data)
    out = np.zeros((len(groups), len(gene_indices)))
    for i, g in enumerate(groups):
        mask = group_labels == g
        if mask.sum() == 0:
            continue
        out[i] = np.asarray(X_log[mask].mean(axis=0)).flatten()
    return pd.DataFrame(out, index=groups.astype(str), columns=gene_indices)


def _diagonal_ordering(mat: pd.DataFrame) -> tuple[list, list]:
    """Pattern-ref: 00_build_structural_package.py::_diagonal_ordering."""
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


def _lineage_rank(lin: str, lineage_order: list[str]) -> int:
    try:
        return lineage_order.index(lin)
    except ValueError:
        return len(lineage_order) + 1


def _build_lineage_groups(
    yml: dict, labels_in_data: set, include_artifact: bool,
    lineage_order: list[str],
) -> list[tuple[str, list[str]]]:
    """Return [(lineage_name, [labels_in_that_lineage]), ...] in canonical order.
    Order pulled from per-compartment `lineage_order` in the config.
    """
    groups: dict[tuple[int, bool, str], list[str]] = {}
    for lbl, rec in yml["labels"].items():
        if lbl not in labels_in_data:
            continue
        is_art = bool(rec.get("is_artifact", False))
        if is_art and not include_artifact:
            continue
        lin = rec.get("lineage", "") or ""
        key = (_lineage_rank(lin, lineage_order), is_art, lin)
        groups.setdefault(key, []).append(lbl)
    ordered = sorted(groups.items(), key=lambda x: x[0])
    return [(k[2], sorted(v)) for k, v in ordered]


def _row_order_within_groups(
    z: pd.DataFrame, groups: list[tuple[str, list[str]]]
) -> list[str]:
    """Row order: preserve group order; within each group sort by argmax-column position."""
    row_order: list[str] = []
    for _, members in groups:
        members_in = [l for l in members if l in z.index]
        if not members_in:
            continue
        if len(members_in) == 1:
            row_order.extend(members_in)
            continue
        sub = z.loc[members_in]
        col_names = list(sub.columns)
        col_to_pos = {c: i for i, c in enumerate(col_names)}
        argmax_pos = sub.idxmax(axis=1).map(col_to_pos)
        max_val = sub.max(axis=1)
        sub_order = (
            pd.DataFrame({"argmax": argmax_pos, "neg_max": -max_val}, index=sub.index)
            .sort_values(["argmax", "neg_max"])
            .index.tolist()
        )
        row_order.extend(sub_order)
    return row_order


def _canonical_block_col_order(
    row_order: list[str], canonical_by_label: dict[str, list[str]], all_cols: list[str]
) -> list[str]:
    """Column order for canonical-marker heatmap: each gene sits under its first
    declaring label in row order. Declaration order is preserved per label. Genes
    declared only by labels not in row_order (e.g., artifacts in non-artifact view)
    are appended at the end.
    """
    col_order: list[str] = []
    seen: set = set()
    cols_set = set(all_cols)
    for lbl in row_order:
        for sym in canonical_by_label.get(lbl, []):
            if sym in cols_set and sym not in seen:
                col_order.append(sym)
                seen.add(sym)
    for sym in all_cols:
        if sym not in seen:
            col_order.append(sym)
    return col_order


def group_frac_expressed(X: sparse.csr_matrix, group_labels: np.ndarray, gene_indices: np.ndarray) -> pd.DataFrame:
    """Fraction of cells per group with each gene detected (> 0)."""
    groups = np.sort(np.unique(group_labels[group_labels != None]))  # noqa: E711
    X_sub = X[:, gene_indices]
    detected = (X_sub > 0).astype(np.int8)
    out = np.zeros((len(groups), len(gene_indices)))
    for i, g in enumerate(groups):
        mask = group_labels == g
        if mask.sum() == 0:
            continue
        out[i] = np.asarray(detected[mask].mean(axis=0)).flatten()
    return pd.DataFrame(out, index=groups.astype(str), columns=gene_indices)


# =============================================================================
# Output: label UMAP with projection
# =============================================================================

def render_label_umap(df: pd.DataFrame, out_pdf: Path) -> None:
    log.info("Rendering label UMAP")
    labels_present = sorted(df["label"].dropna().unique())
    palette = _categorical_palette(len(labels_present))
    color_map = dict(zip(labels_present, palette))

    fig, ax = plt.subplots(figsize=(9, 9))
    sub = df.dropna(subset=["label"])
    colors = sub["label"].map(color_map).values
    ax.scatter(sub["UMAP_1"], sub["UMAP_2"], c=colors, s=0.4, alpha=0.7, rasterized=True, linewidths=0)

    # Label centroids (project on plot)
    texts = []
    for lbl in labels_present:
        cells = sub[sub["label"] == lbl]
        if len(cells) < 50:
            continue
        x = float(cells["UMAP_1"].median())
        y = float(cells["UMAP_2"].median())
        t = ax.text(
            x, y, lbl,
            fontsize=5, ha="center", va="center", color="black",
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="black", lw=0.3, alpha=0.85),
            zorder=10,
        )
        texts.append(t)
    try:
        from adjustText import adjust_text
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.25, alpha=0.6),
                    expand_text=(1.05, 1.1), expand_points=(1.05, 1.1),
                    force_text=(0.2, 0.3), force_points=(0.05, 0.1),
                    only_move={"text": "xy"}, lim=50)
    except ImportError:
        log.warning("adjustText not installed; labels rendered without repulsion")

    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.3)
    ax.set_title(f"Proposed labels ({len(labels_present)})", fontsize=7)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# =============================================================================
# Output: heatmap_canonical (labels × declared canonical markers)
# =============================================================================

def render_heatmap_canonical(
    yml: dict,
    df: pd.DataFrame,
    X: sparse.csr_matrix,
    ensg_to_idx: dict,
    symbol_to_ensg: dict,
    out_pdf: Path,
    lineage_order: list[str],
    include_artifact: bool = True,
) -> None:
    log.info("Rendering canonical-marker heatmap (include_artifact=%s)", include_artifact)
    # Filter labels by include_artifact; collect declared canonical markers for that set
    candidate_labels = {
        lbl for lbl, rec in yml["labels"].items()
        if include_artifact or not rec.get("is_artifact", False)
    }
    canonical_by_label: dict[str, list[str]] = {}
    all_symbols: list[str] = []
    for label, rec in yml["labels"].items():
        if label not in candidate_labels:
            continue
        syms = list(rec.get("canonical_markers") or [])
        canonical_by_label[label] = syms
        for s in syms:
            if s not in all_symbols:
                all_symbols.append(s)
    if not all_symbols:
        log.warning("  no canonical markers declared; skipping")
        return

    # Resolve to gene indices
    sym_to_gidx = {}
    for s in all_symbols:
        ensg = symbol_to_ensg.get(s)
        if ensg and ensg in ensg_to_idx:
            sym_to_gidx[s] = ensg_to_idx[ensg]
    resolved = [s for s in all_symbols if s in sym_to_gidx]
    dropped = [s for s in all_symbols if s not in sym_to_gidx]
    if dropped:
        log.warning("  %d canonical markers not in gene_data: %s", len(dropped), dropped[:10])
    gene_indices = np.array([sym_to_gidx[s] for s in resolved], dtype=int)

    labels_in_data = set(df["label"].unique()) & candidate_labels
    groups = _build_lineage_groups(yml, labels_in_data, include_artifact=include_artifact,
                                    lineage_order=lineage_order)

    # Compute mean log1p per label × gene; z across labels
    group_labels = df["label"].values
    mean_expr = cluster_mean_log1p(X, group_labels, gene_indices)
    mean_expr.columns = resolved
    mean_expr = mean_expr.loc[mean_expr.index.isin(labels_in_data)]
    z = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0).replace(0, 1))

    # Row order: lineage groups; within-group by argmax column. Column order: each
    # canonical marker placed under its first declaring label → true per-label blocks.
    row_order = _row_order_within_groups(z, groups)
    col_order = _canonical_block_col_order(row_order, canonical_by_label, list(z.columns))
    z = z.loc[row_order, col_order]
    labels_order = row_order
    resolved = list(col_order)

    n_labels = len(labels_order)
    n_genes = len(resolved)
    height = max(4.0, 0.18 * n_labels + 1.5)
    width = max(8.0, 0.16 * n_genes + 3.0)

    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(z, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                xticklabels=resolved, yticklabels=labels_order,
                cbar_kws={"label": "z-score", "shrink": 0.3}, linewidths=0, ax=ax)
    _draw_group_separators(ax, labels_order, groups)
    title_suffix = " + artifacts" if include_artifact else ""
    ax.set_xlabel("canonical markers (declared)", fontsize=7)
    ax.set_ylabel("")
    ax.set_title(f"Declared canonical markers × labels (z-scored){title_suffix}", fontsize=7)
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=5)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


def _draw_group_separators(ax, row_order: list[str], groups: list[tuple[str, list[str]]]) -> None:
    """Thin horizontal lines between lineage groups to reinforce the block structure."""
    seen = 0
    for i, (_, members) in enumerate(groups):
        seen += len([m for m in members if m in row_order])
        if i < len(groups) - 1 and 0 < seen < len(row_order):
            ax.axhline(seen, color="black", lw=0.4)


# =============================================================================
# Output: heatmap_top (labels × top DE markers, one-vs-rest per label)
# =============================================================================

def compute_top_markers_per_label(
    df: pd.DataFrame, X: sparse.csr_matrix, labels_order: list[str],
    gene_data: pd.DataFrame, top_n: int,
    candidate_labels: set | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """One-vs-rest log-fold-change per label. Returns (label × gene z-matrix, gene_order).

    df + X must be row-aligned (full compartment). `candidate_labels` restricts the
    "rest" pool so exclusions (e.g., artifacts) don't inflate rest-mean.
    """
    if candidate_labels is None:
        candidate_labels = set(labels_order)
    try:
        import scanpy as sc
        import anndata as ad
    except Exception as e:
        log.warning("scanpy/anndata import failed (%s: %s); falling back to raw logFC computation",
                    type(e).__name__, str(e))
        return _top_markers_raw(df, X, labels_order, gene_data, top_n, candidate_labels)

    log.info("  running rank_genes_groups (wilcoxon) on labels")
    # Build AnnData from sparse counts + label obs
    adata = ad.AnnData(X=X, obs=df.set_index("cell_id")[["label"]], var=gene_data.reset_index().rename(columns={"index": "ensg"}).set_index("ensg"))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    # Rank using labels
    adata.obs["label"] = adata.obs["label"].astype("category")
    sc.tl.rank_genes_groups(adata, groupby="label", method="wilcoxon", use_raw=False)

    # Select top_n per label
    picks: list[tuple[str, str]] = []
    seen_genes: set = set()
    for lbl in labels_order:
        if lbl not in adata.uns["rank_genes_groups"]["names"].dtype.names:
            continue
        names = adata.uns["rank_genes_groups"]["names"][lbl][:top_n * 3]  # over-collect, dedup
        for g in names:
            if g in seen_genes:
                continue
            picks.append((lbl, g))
            seen_genes.add(g)
            if len([p for p in picks if p[0] == lbl]) >= top_n:
                break

    gene_order = [g for _, g in picks]
    if not gene_order:
        log.warning("  no top markers found")
        return pd.DataFrame(), []

    # Mean log1p per label × gene
    gene_indices = np.array([gene_data.index.get_loc(g) for g in gene_order if g in gene_data.index], dtype=int)
    resolved_genes = [g for g in gene_order if g in gene_data.index]
    mean_expr = cluster_mean_log1p(X, df["label"].values, gene_indices)
    mean_expr.columns = resolved_genes
    mean_expr = mean_expr.reindex(labels_order)
    z = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0).replace(0, 1))
    return z, resolved_genes


def _top_markers_raw(df, X, labels_order, gene_data, top_n, candidate_labels: set | None = None):
    log.info("  computing raw one-vs-rest logFC per label (memory-efficient, 1-pass sums)")
    symbols_series = gene_data["gene_symbol"]
    has_symbol = symbols_series.notna() & (symbols_series.astype(str).str.len() > 0)
    valid_gene_mask = has_symbol.values
    picks: list[tuple[str, str]] = []
    seen_genes: set = set()
    group_labels = df["label"].values
    if candidate_labels is None:
        candidate_labels = set(labels_order)
    in_candidates = np.isin(group_labels, list(candidate_labels))

    # Reconstruct X_log via shared buffers (no second full copy of indices/indptr):
    # only the data array doubles in memory. Saves ~8GB on 1B-nnz matrices.
    X_log_data = np.log1p(X.data)
    X_log = sparse.csr_matrix((X_log_data, X.indices, X.indptr), shape=X.shape)

    # Pre-compute per-label sum (one sparse-mask reduction per label) and total.
    # X_log[mask].sum(axis=0) returns a 1×G dense row; cast to flat array.
    total_sum = np.asarray(X_log[in_candidates].sum(axis=0)).flatten()
    total_n = int(in_candidates.sum())
    label_sums: dict[str, tuple[np.ndarray, int]] = {}
    for lbl in labels_order:
        mask = (group_labels == lbl) & in_candidates
        n = int(mask.sum())
        if n == 0:
            continue
        label_sums[lbl] = (np.asarray(X_log[mask].sum(axis=0)).flatten(), n)

    # logFC = label_mean - rest_mean; rest_sum = total_sum - label_sum
    for lbl in labels_order:
        if lbl not in label_sums:
            continue
        s, n = label_sums[lbl]
        rest_s = total_sum - s
        rest_n = total_n - n
        if rest_n == 0:
            continue
        logfc = (s / n) - (rest_s / rest_n)
        # Exclude genes without a symbol (would render as ENSG)
        logfc_masked = np.where(valid_gene_mask, logfc, -np.inf)
        order = np.argsort(-logfc_masked)
        added = 0
        for idx in order:
            ensg = gene_data.index[idx]
            if ensg in seen_genes:
                continue
            picks.append((lbl, ensg))
            seen_genes.add(ensg)
            added += 1
            if added >= top_n:
                break

    gene_order = [g for _, g in picks]
    if not gene_order:
        return pd.DataFrame(), []
    gene_indices = np.array([gene_data.index.get_loc(g) for g in gene_order], dtype=int)
    mean_expr = cluster_mean_log1p(X, group_labels, gene_indices)
    mean_expr.columns = gene_order
    mean_expr = mean_expr.reindex(labels_order)
    z = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0).replace(0, 1))
    return z, gene_order


def render_heatmap_top(
    yml: dict, df: pd.DataFrame, X: sparse.csr_matrix, gene_data: pd.DataFrame,
    ensg_to_symbol: dict, top_n: int, out_pdf: Path,
    lineage_order: list[str],
    include_artifact: bool = True,
) -> None:
    log.info("Rendering top-markers heatmap (include_artifact=%s)", include_artifact)
    # Pre-order labels by lineage so gene picks follow the canonical biology order
    labels_in_data = {
        l for l in yml["labels"].keys()
        if l in df["label"].unique()
        and (include_artifact or not yml["labels"][l].get("is_artifact", False))
    }
    if not labels_in_data:
        log.warning("  no labels; skipping")
        return
    groups = _build_lineage_groups(yml, labels_in_data, include_artifact=include_artifact,
                                    lineage_order=lineage_order)
    labels_order_pre = [l for _, members in groups for l in members]

    # Pass full df + X (row-aligned); candidate_labels restricts one-vs-rest pool
    z, gene_order = compute_top_markers_per_label(
        df, X, labels_order_pre, gene_data, top_n, candidate_labels=labels_in_data,
    )
    if z.empty:
        log.warning("  skipping (no markers)")
        return

    # Genes were picked in labels_order_pre order (top_n per label, deduped) so
    # gene_order is already label-blocked. Keep that; only drop rows not in data.
    z.columns = [ensg_to_symbol.get(g, g) for g in gene_order]
    row_order = [l for l in labels_order_pre if l in z.index]
    z = z.loc[row_order]

    symbol_labels = list(z.columns)
    n_labels = len(z.index)
    n_genes = len(symbol_labels)
    height = max(4.0, 0.18 * n_labels + 1.5)
    width = max(8.0, 0.14 * n_genes + 3.0)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(z, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                xticklabels=symbol_labels, yticklabels=True,
                cbar_kws={"label": "z-score", "shrink": 0.3}, linewidths=0, ax=ax)
    _draw_group_separators(ax, row_order, groups)
    title_suffix = " + artifacts" if include_artifact else ""
    ax.set_xlabel(f"top {top_n} markers per label (one-vs-rest)", fontsize=7)
    ax.set_ylabel("")
    ax.set_title(f"Top DE markers × labels (z-scored){title_suffix}", fontsize=7)
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=5)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# =============================================================================
# Output: heatmap_lineage_tiled / heatmap_artifact_tiled (clusters × markers tiled across resolutions)
# =============================================================================

def render_tiled_marker_heatmap(
    df: pd.DataFrame,
    X: sparse.csr_matrix,
    symbol_groups: dict[str, list[str]],
    symbol_to_ensg: dict,
    ensg_to_idx: dict,
    res_cols: list[str],
    out_pdf: Path,
    title_prefix: str,
) -> None:
    """Grid of (cluster × gene) z-score heatmaps, one panel per resolution."""
    log.info("Rendering tiled marker heatmap: %s", title_prefix)
    # Flatten markers with group label → pass-through list; keep group boundaries for separators
    flat_symbols: list[str] = []
    boundaries: list[int] = []
    group_labels_x: list[str] = []
    for group, symbols in symbol_groups.items():
        for s in symbols:
            flat_symbols.append(s)
            group_labels_x.append(group)
        boundaries.append(len(flat_symbols))

    gene_indices = []
    resolved_symbols = []
    resolved_groups = []
    for s, g in zip(flat_symbols, group_labels_x):
        ensg = symbol_to_ensg.get(s)
        if ensg and ensg in ensg_to_idx:
            gene_indices.append(ensg_to_idx[ensg])
            resolved_symbols.append(s)
            resolved_groups.append(g)
    if not gene_indices:
        log.warning("  no markers resolved; skipping")
        return
    gene_indices = np.array(gene_indices, dtype=int)

    # 2-row layout: top row holds the first ceil(N/2) (coarse resolutions),
    # bottom row holds the remaining (fine resolutions).
    n = len(res_cols)
    ncols_plot = int(np.ceil(n / 2)) if n > 1 else 1
    nrows_plot = 2 if n > 1 else 1
    fig, axes = plt.subplots(nrows_plot, ncols_plot, figsize=(ncols_plot * 5, nrows_plot * 5))
    axes = np.atleast_2d(axes).flatten()

    for i, res_col in enumerate(res_cols):
        ax = axes[i]
        if res_col not in df.columns:
            ax.set_visible(False)
            continue
        cluster_labels = df[res_col].astype(str).values
        mean_expr = cluster_mean_log1p(X, cluster_labels, gene_indices)
        mean_expr.columns = resolved_symbols
        z = (mean_expr - mean_expr.mean(axis=0)) / (mean_expr.std(axis=0).replace(0, 1))
        # Rotate: genes as rows (few, <50), clusters as cols (up to ~hundreds at high res)
        z_t = z.T
        # Order cluster cols numerically where possible
        try:
            col_order = sorted(z_t.columns, key=lambda c: int(str(c)))
            z_t = z_t[col_order]
        except ValueError:
            pass
        sns.heatmap(z_t, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
                    xticklabels=True, yticklabels=resolved_symbols,
                    cbar=False, linewidths=0, ax=ax)
        ax.set_title(f"{title_prefix} @ {res_col}", fontsize=6)
        ax.set_xlabel("cluster", fontsize=5)
        ax.tick_params(axis="x", labelsize=4, rotation=90)
        ax.tick_params(axis="y", labelsize=4)
        # Horizontal separators between marker groups (rows now)
        for b in boundaries[:-1]:
            if 0 < b < len(resolved_symbols):
                ax.axhline(b, color="black", lw=0.3)

    for j in range(len(res_cols), len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# =============================================================================
# Output: artifact_compromise heatmap (label × artifact category, mean-z)
# =============================================================================

def render_artifact_compromise_heatmap(
    yml: dict,
    df: pd.DataFrame,
    X: sparse.csr_matrix,
    symbol_to_ensg: dict,
    ensg_to_idx: dict,
    artifact_cfg: dict,
    z_threshold: float,
    out_pdf: Path,
    lineage_order: list[str],
) -> None:
    """Label × artifact-category mean-z heatmap.

    For each (label, category): mean log1p of category genes in label cells, compared
    against mean of the same genes across all other cells; reported as z-score using
    the rest pool's std. Mirrors the metric in 04_compute_evidence.py §7.2.

    Ordered by lineage (non-artifact first, then artifact-like), matching the other
    marker heatmaps. Cell values capped to [-z_threshold, +z_threshold] for readability.
    """
    log.info("Rendering artifact-compromise heatmap")
    labels_in_data = set(df["label"].dropna().unique())
    groups = _build_lineage_groups(yml, labels_in_data, include_artifact=True,
                                    lineage_order=lineage_order)
    labels_ordered = [l for _, members in groups for l in members]

    categories = list(artifact_cfg.keys())
    # Resolve gene indices per category once
    cat_gene_indices: dict[str, np.ndarray] = {}
    for cat, symbols in artifact_cfg.items():
        idxs = []
        for s in symbols:
            ensg = symbol_to_ensg.get(s)
            if ensg and ensg in ensg_to_idx:
                idxs.append(ensg_to_idx[ensg])
        cat_gene_indices[cat] = np.array(idxs, dtype=int) if idxs else np.array([], dtype=int)

    group_labels = df["label"].values
    # Build the z-matrix
    z_mat = np.zeros((len(labels_ordered), len(categories)))
    for ci, cat in enumerate(categories):
        gidx = cat_gene_indices[cat]
        if gidx.size == 0:
            z_mat[:, ci] = np.nan
            continue
        # mean log1p across category genes, one value per cell
        X_sub = X[:, gidx]
        X_log = X_sub.copy()
        X_log.data = np.log1p(X_log.data)
        per_cell = np.asarray(X_log.mean(axis=1)).flatten()
        for ri, lbl in enumerate(labels_ordered):
            mask = group_labels == lbl
            rest = ~mask
            if mask.sum() == 0 or rest.sum() == 0:
                z_mat[ri, ci] = np.nan
                continue
            mu_rest = float(per_cell[rest].mean())
            sd_rest = float(per_cell[rest].std(ddof=0))
            if sd_rest <= 0:
                z_mat[ri, ci] = 0.0
                continue
            z_mat[ri, ci] = (per_cell[mask].mean() - mu_rest) / sd_rest

    z_df = pd.DataFrame(z_mat, index=labels_ordered, columns=categories)

    vmax = max(float(z_threshold), 3.0)
    n_labels = len(labels_ordered)
    n_cats = len(categories)
    height = max(4.0, 0.22 * n_labels + 1.5)
    width = max(6.0, 0.6 * n_cats + 4.0)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        z_df, cmap="RdBu_r", center=0, vmin=-vmax, vmax=vmax,
        xticklabels=categories, yticklabels=labels_ordered,
        cbar_kws={"label": f"mean z (clipped ±{vmax:g})", "shrink": 0.3},
        linewidths=0.1, linecolor="white", ax=ax,
        annot=True, fmt=".2f", annot_kws={"fontsize": 4},
    )
    _draw_group_separators(ax, labels_ordered, groups)
    ax.set_xlabel("artifact category", fontsize=7)
    ax.set_ylabel("")
    ax.set_title(f"Label compromise score per artifact category (flag at z ≥ {z_threshold:g})", fontsize=7)
    ax.tick_params(axis="x", labelsize=5, rotation=45)
    ax.tick_params(axis="y", labelsize=5)
    for tick in ax.get_xticklabels():
        tick.set_ha("right")
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# =============================================================================
# Output: label_study_composition / label_patient_composition barplots
# =============================================================================

def render_patient_composition_by_study(
    df: pd.DataFrame, out_pdf: Path, labels_order: list[str],
) -> None:
    """Patient composition tiled per study: one panel per study, bars = patient fractions
    within that study's contribution to each label. Exposes donor-concentration within
    each study's contribution.
    """
    log.info("Rendering patient composition tiled by study")
    studies = sorted(df["dataset"].dropna().unique())
    n = len(studies)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    row_height = max(3.0, 0.18 * len(labels_order) + 1.0)
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4.5, nrows * row_height), sharey=True)
    axes = np.atleast_2d(axes).flatten()

    for si, study in enumerate(studies):
        ax = axes[si]
        sub = df[df["dataset"] == study]
        tab = (
            sub.groupby(["label", "patientID"], observed=True).size().unstack(fill_value=0)
            .reindex(labels_order, fill_value=0)
        )
        row_sums = tab.sum(axis=1)
        frac = tab.div(row_sums.where(row_sums > 0, 1), axis=0).fillna(0)
        categories = list(frac.columns)
        if not categories:
            ax.set_title(f"{study} (no cells)", fontsize=6)
            ax.set_visible(False)
            continue
        palette = _categorical_palette(len(categories))
        bottom = np.zeros(len(labels_order))
        for cat, color in zip(categories, palette):
            vals = frac[cat].values
            ax.barh(range(len(labels_order)), vals, left=bottom, color=color,
                    edgecolor="white", lw=0.15)
            bottom = bottom + vals
        ax.set_yticks(range(len(labels_order)))
        ax.set_yticklabels(labels_order, fontsize=4)
        ax.set_xlabel("fraction of label (within study)", fontsize=5)
        ax.set_xlim(0, 1)
        ax.invert_yaxis()
        ax.set_title(f"{study}  (n_donors={len(categories)}, n_cells={int(row_sums.sum()):,})", fontsize=6)
        ax.tick_params(axis="x", labelsize=4)
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    for j in range(n, len(axes)):
        axes[j].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


def render_composition_barplot(
    df: pd.DataFrame, by: str, title: str, out_pdf: Path, labels_order: list[str],
) -> None:
    log.info("Rendering composition barplot: %s", by)
    tab = (
        df.groupby(["label", by], observed=True).size().unstack(fill_value=0)
        .reindex(labels_order, fill_value=0)
    )
    # Fraction per label
    frac = tab.div(tab.sum(axis=1), axis=0).fillna(0)

    categories = list(frac.columns)
    palette = _categorical_palette(len(categories))

    height = max(4.0, 0.25 * len(labels_order) + 1.5)
    fig, ax = plt.subplots(figsize=(8, height))
    bottom = np.zeros(len(labels_order))
    for cat, color in zip(categories, palette):
        vals = frac[cat].values
        ax.barh(range(len(labels_order)), vals, left=bottom, color=color, label=cat, edgecolor="white", lw=0.2)
        bottom = bottom + vals
    ax.set_yticks(range(len(labels_order)))
    ax.set_yticklabels(labels_order, fontsize=5)
    ax.set_xlabel("fraction of label", fontsize=7)
    ax.set_xlim(0, 1)
    ax.invert_yaxis()
    ax.set_title(title, fontsize=7)
    ax.tick_params(axis="x", labelsize=6)
    ax.legend(fontsize=5, bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False, ncol=1)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()
    annotation_dir = args.out_dir / "annotation"
    annotation_dir.mkdir(parents=True, exist_ok=True)

    yml = yaml.safe_load(args.v2_yaml.read_text())
    cfg = yaml.safe_load(args.config.read_text())
    compartment = yml.get("compartment", "imm")
    comp_cfg = cfg["compartments"][compartment]
    lineage_order: list[str] = comp_cfg.get("lineage_order", [])

    log.info("Loading labels: %s", args.labels_csv)
    labels_df = pd.read_csv(args.labels_csv)
    log.info("Loading leiden: %s", args.leiden_csv)
    leiden = pd.read_csv(args.leiden_csv)
    log.info("Loading umap: %s", args.umap_csv)
    umap = pd.read_csv(args.umap_csv)
    log.info("Loading metadata: %s", args.metadata_csv)
    meta = pd.read_csv(args.metadata_csv, usecols=["cell_id", "dataset", "patientID"], low_memory=False)
    log.info("Loading gene_data: %s", args.gene_data)
    gene_data = pd.read_csv(args.gene_data, index_col=0)
    ensg_to_idx = {ensg: i for i, ensg in enumerate(gene_data.index)}
    ensg_to_symbol = gene_data["gene_symbol"].to_dict()
    symbol_to_ensg = {sym: ensg for ensg, sym in ensg_to_symbol.items() if isinstance(sym, str)}

    log.info("Loading counts: %s", args.counts_npz)
    loaded = np.load(args.counts_npz, allow_pickle=True)
    X = sparse.csr_matrix(
        (loaded["data"], loaded["indices"], loaded["indptr"]),
        shape=tuple(loaded["shape"]),
    )

    df = leiden.merge(umap, on="cell_id", how="inner").merge(meta, on="cell_id", how="inner").merge(
        labels_df[["cell_id", "label", "lineage", "is_artifact"]], on="cell_id", how="inner"
    )
    log.info("Joined rows: %d", len(df))

    # Tile resolutions (D3: curated set; include any used for fine_cluster labels)
    curated = [f"leiden_{r}" for r in cfg["pipeline"]["heatmap_tile_resolutions"]]
    fine_res = {f"leiden_{r.get('source_resolution')}" for r in yml["labels"].values() if r.get("mode") == "fine_cluster" and r.get("source_resolution")}
    tile_resolutions = [c for c in (curated + list(fine_res - set(curated))) if c in df.columns]

    # Labels order (non-artifact first, by lineage)
    labels_order = sorted(
        [l for l in yml["labels"].keys() if l in df["label"].unique()],
        key=lambda lbl: (yml["labels"][lbl].get("is_artifact", False),
                         yml["labels"][lbl].get("lineage", ""),
                         lbl),
    )

    # --- Renders ---
    render_label_umap(df, annotation_dir / "label_umap.pdf")
    # Canonical heatmaps: non-artifact-only (publication) + with-artifact (QA)
    render_heatmap_canonical(yml, df, X, ensg_to_idx, symbol_to_ensg,
                             annotation_dir / "heatmap_canonical.pdf",
                             lineage_order=lineage_order,
                             include_artifact=False)
    render_heatmap_canonical(yml, df, X, ensg_to_idx, symbol_to_ensg,
                             annotation_dir / "heatmap_canonical_with_artifact.pdf",
                             lineage_order=lineage_order,
                             include_artifact=True)
    # Top DE heatmaps: same split
    render_heatmap_top(yml, df, X, gene_data, ensg_to_symbol, args.top_n_per_label,
                       annotation_dir / "heatmap_top.pdf",
                       lineage_order=lineage_order,
                       include_artifact=False)
    render_heatmap_top(yml, df, X, gene_data, ensg_to_symbol, args.top_n_per_label,
                       annotation_dir / "heatmap_top_with_artifact.pdf",
                       lineage_order=lineage_order,
                       include_artifact=True)
    render_tiled_marker_heatmap(
        df, X, comp_cfg["lineage_markers"], symbol_to_ensg, ensg_to_idx,
        tile_resolutions, annotation_dir / "heatmap_lineage_tiled.pdf", "Lineage markers"
    )
    render_tiled_marker_heatmap(
        df, X, comp_cfg["artifact_markers"], symbol_to_ensg, ensg_to_idx,
        tile_resolutions, annotation_dir / "heatmap_artifact_tiled.pdf", "Artifact markers"
    )
    # Artifact-compromise heatmap: composite label × artifact-category z
    z_threshold = float(cfg["pipeline"]["flagging"]["artifact_enrichment"]["z_threshold"])
    render_artifact_compromise_heatmap(
        yml, df, X, symbol_to_ensg, ensg_to_idx,
        comp_cfg["artifact_markers"], z_threshold,
        annotation_dir / "heatmap_artifact_compromise.pdf",
        lineage_order=lineage_order,
    )
    render_composition_barplot(
        df, by="dataset", title="Study composition per label",
        out_pdf=annotation_dir / "label_study_composition.pdf",
        labels_order=labels_order,
    )
    render_composition_barplot(
        df, by="patientID", title="Patient composition per label",
        out_pdf=annotation_dir / "label_patient_composition.pdf",
        labels_order=labels_order,
    )
    render_patient_composition_by_study(
        df, annotation_dir / "label_patient_composition_by_study.pdf", labels_order,
    )

    log.info("Annotation package complete: %s", annotation_dir)


if __name__ == "__main__":
    main()
