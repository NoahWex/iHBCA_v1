"""05c_featureplots.py — UMAP featureplots for V1 yaml markers per compartment.

For each V1 label in the compartment yaml, render a grid of UMAP scatters colored by
per-cell normalized expression of its canonical_markers + identity_markers.

One PDF per V1 label written to {out_dir}/structural/featureplots/{label}.pdf.

Inputs:
    --bundle-dir            scvi_n100 bundle (cells.tsv, genes.tsv, counts.mtx.gz, umap.csv)
    --yaml                  V1 annotation yaml (publication/analysis/annotation/yamls/annotation_v2_*.yaml)
    --out-dir               compartment output dir (e.g. outputs/Immune/)
    --max-markers           max markers per label to plot (default 12)
    --point-size            scatter point size (default 0.5)
    --pct-cap               percentile cap for color scaling (default 99)

Pattern: standalone (not adapted from a single source). Uses the same bundle-loading
approach as 03a_pseudobulk.py (cells.tsv / genes.tsv / counts.mtx.gz orientation check).
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
from matplotlib.backends.backend_pdf import PdfPages
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
log = logging.getLogger("featureplots")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--bundle-dir", type=Path, required=True)
    p.add_argument("--yaml", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--max-markers", type=int, default=12)
    p.add_argument("--point-size", type=float, default=0.5)
    p.add_argument("--pct-cap", type=float, default=99.0)
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


def collect_markers(yaml_path: Path) -> Dict[str, List[str]]:
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
        if ordered:
            markers[label] = ordered
    return markers


def render_label(
    label: str,
    markers: List[str],
    gene_index: pd.Series,
    Xn: sp.csr_matrix,
    umap: np.ndarray,
    out_pdf: Path,
    max_markers: int,
    point_size: float,
    pct_cap: float,
) -> None:
    use = [g for g in markers if g in gene_index.index][:max_markers]
    missing = [g for g in markers if g not in gene_index.index]
    if not use:
        log.warning("[%s] no markers found in genes.tsv (%d requested missing)", label, len(missing))
        return
    n = len(use)
    ncol = min(4, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(2.6 * ncol, 2.6 * nrow), squeeze=False, dpi=150
    )
    fig.suptitle(f"{label} — V1 canonical + identity markers", fontsize=10)
    for ax in axes.flat:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_aspect("equal")
    for i, gene in enumerate(use):
        gi = gene_index[gene]
        col = Xn[:, gi].toarray().ravel()
        cap = np.percentile(col[col > 0], pct_cap) if (col > 0).any() else 1.0
        cap = max(cap, 1e-6)
        ax = axes[i // ncol, i % ncol]
        order = np.argsort(col)
        ax.scatter(
            umap[order, 0],
            umap[order, 1],
            c=np.clip(col[order], 0, cap),
            cmap="viridis",
            s=point_size,
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(gene, fontsize=8)
    for j in range(len(use), nrow * ncol):
        axes[j // ncol, j % ncol].axis("off")
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info(
        "[%s] %d markers plotted (%d missing) -> %s",
        label,
        len(use),
        len(missing),
        out_pdf.name,
    )


def main() -> None:
    args = parse_args()

    cells_tsv = args.bundle_dir / "cells.tsv"
    genes_tsv = args.bundle_dir / "genes.tsv"
    mtx = args.bundle_dir / "counts.mtx.gz"
    umap_csv = args.bundle_dir / "umap.csv"
    for f in (cells_tsv, genes_tsv, mtx, umap_csv, args.yaml):
        if not f.exists():
            sys.exit(f"ERROR: missing input: {f}")

    cells = pd.read_csv(cells_tsv, sep="\t", header=None, names=["cell_id"])
    genes = pd.read_csv(genes_tsv, sep="\t", header=None, names=["gene"])
    umap_df = pd.read_csv(umap_csv)
    umap_df.columns = [c.lower() for c in umap_df.columns]
    if "cell_id" in umap_df.columns:
        umap_df = umap_df.set_index("cell_id").reindex(cells["cell_id"]).reset_index()
    if not {"umap_1", "umap_2"}.issubset(umap_df.columns):
        cols = [c for c in umap_df.columns if c != "cell_id"]
        umap_df = umap_df.rename(columns={cols[0]: "umap_1", cols[1]: "umap_2"})
    umap = umap_df[["umap_1", "umap_2"]].to_numpy()
    log.info("Cells=%d, Genes=%d, UMAP=%s", len(cells), len(genes), umap.shape)

    X = load_mtx_gz(mtx)
    if X.shape == (len(cells), len(genes)):
        log.info("mtx orientation: cells × genes")
    elif X.shape == (len(genes), len(cells)):
        log.info("mtx orientation: genes × cells (transposing)")
        X = X.T.tocsr()
    else:
        sys.exit(f"ERROR: mtx shape {X.shape} does not match cells/genes")

    log.info("Normalizing (target_sum=%s, log1p)", args.target_sum)
    Xn = normalize_log1p(X, args.target_sum)

    gene_index = pd.Series(np.arange(len(genes)), index=genes["gene"].values)
    if gene_index.index.has_duplicates:
        log.warning("Duplicate gene symbols in genes.tsv; keeping first occurrence")
        gene_index = gene_index[~gene_index.index.duplicated(keep="first")]

    markers_by_label = collect_markers(args.yaml)
    log.info("Markers per label loaded for %d labels", len(markers_by_label))

    out_root = args.out_dir / "structural" / "featureplots"
    out_root.mkdir(parents=True, exist_ok=True)
    for label, markers in markers_by_label.items():
        safe = label.replace("/", "_").replace(" ", "_")
        render_label(
            label,
            markers,
            gene_index,
            Xn,
            umap,
            out_root / f"{safe}.pdf",
            max_markers=args.max_markers,
            point_size=args.point_size,
            pct_cap=args.pct_cap,
        )
    log.info("Done. Outputs at %s", out_root)


if __name__ == "__main__":
    main()
