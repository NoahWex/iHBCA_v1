"""06_render_annotation.py — post-yaml evidence package per compartment.

Reads the per-compartment annotation yaml (anchor labels + fine_cluster_overrides),
applies the two-pass label assignment (Pass 1: cluster mode; Pass 2: fine_cluster
mode overrides) to per-cell cluster IDs, and renders the post-review evidence
package:

    label_umap.pdf                  cells colored by assigned L2S label, centroid text
    heatmap_canonical.pdf           labels × declared canonical+identity markers (z-scored)
    heatmap_canonical_with_artifact.pdf
    heatmap_top.pdf                 labels × top one-vs-rest markers (aggregated from
                                    per-cluster limma at the anchor resolution)
    heatmap_top_with_artifact.pdf
    label_summary.csv               n_cells per label

Heatmap styling: column_block_order grouping, RdBu_r z-score, group separators.
Label UMAP uses adjustText for force-directed centroid label repulsion.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Set, Tuple

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import seaborn as sns
import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("render_annotation")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--bundle-dir", type=Path, required=True)
    p.add_argument("--clusters-csv", type=Path, required=True)
    p.add_argument("--yaml-draft", type=Path, required=True, help="annotation_v2s_DRAFT.yaml")
    p.add_argument("--v1-yaml", type=Path, required=True, help="V1 annotation_v2_*.yaml for canonical markers")
    p.add_argument("--limma-dir", type=Path, required=True, help="Dir with limma_markers_{Comp}_leiden_*.csv")
    p.add_argument("--compartment", required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--top-n-per-label", type=int, default=10)
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


# ---------- yaml application -------------------------------------------------

def apply_yaml(yml: dict, clusters: pd.DataFrame) -> pd.Series:
    """Per-cell L2S label from anchor + fine_cluster_overrides."""
    anchor = str(yml["anchor_resolution"])
    anchor_col = f"leiden_{anchor}"
    if anchor_col not in clusters.columns:
        sys.exit(f"ERROR: anchor column {anchor_col} missing from clusters")
    anchor_cells = clusters[anchor_col].astype(int)
    base = yml["clusters_by_resolution"][anchor]
    cluster_to_label: Dict[int, str] = {}
    for ckey, body in base.items():
        cnum = int(str(ckey).lstrip("c"))
        cluster_to_label[cnum] = body.get("assigned_label", "PENDING")
    labels = anchor_cells.map(cluster_to_label).fillna("PENDING")

    for group_name, group_body in (yml.get("fine_cluster_overrides") or {}).items():
        if not isinstance(group_body, dict):
            continue
        for ov in (group_body.get("overrides") or []):
            res = str(ov["resolution"])
            cnum = int(ov["cluster"])
            new_label = str(ov["label"])
            col = f"leiden_{res}"
            if col not in clusters.columns:
                log.warning("[%s] override resolution %s missing in clusters CSV", group_name, res)
                continue
            mask = (clusters[col].astype(int) == cnum).values
            n_changed = int(mask.sum())
            if n_changed == 0:
                continue
            labels = pd.Series(np.where(mask, new_label, labels.values), index=labels.index)
            log.info("Override [%s] res=%s c=%d → %s (%d cells)", group_name, res, cnum, new_label, n_changed)
    return labels


# ---------- shared helpers ---------------------------------------------------

def cluster_mean_log1p(Xn: sp.csr_matrix, group_labels: np.ndarray, gene_indices: np.ndarray) -> pd.DataFrame:
    sub = Xn[:, gene_indices]
    df = pd.DataFrame(sub.toarray())
    df["__g"] = group_labels
    return df.groupby("__g").mean()


def collect_canonical_markers(v1_yaml_path: Path, max_per_label: int) -> Tuple[Dict[str, List[str]], List[str]]:
    with open(v1_yaml_path) as fh:
        cfg = yaml.safe_load(fh)
    markers: Dict[str, List[str]] = {}
    label_order: List[str] = []
    for label, body in (cfg.get("labels") or {}).items():
        if not isinstance(body, dict):
            continue
        canon = body.get("canonical_markers") or []
        ident = body.get("identity_markers") or []
        seen: Set[str] = set()
        ordered: List[str] = []
        for g in list(canon) + list(ident):
            if g and g not in seen:
                seen.add(g)
                ordered.append(g)
            if len(ordered) >= max_per_label:
                break
        if ordered:
            markers[label] = ordered
            label_order.append(label)
    return markers, label_order


def column_block_order(label_to_markers: Dict[str, List[str]], label_order: List[str]) -> List[Tuple[str, str]]:
    placed: Dict[str, str] = {}
    out: List[Tuple[str, str]] = []
    for label in label_order:
        for g in label_to_markers.get(label, []):
            if g not in placed:
                placed[g] = label
                out.append((g, label))
    return out


def render_block_heatmap(
    Z: pd.DataFrame,
    col_groups: List[str],
    out_pdf: Path,
    xlabel: str,
) -> None:
    n_rows, n_cols = Z.shape
    height = max(4.0, 0.18 * n_rows + 1.5)
    width = max(8.0, 0.16 * n_cols + 3.0)
    fig, ax = plt.subplots(figsize=(width, height))
    sns.heatmap(
        Z, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        xticklabels=Z.columns.tolist(), yticklabels=Z.index.tolist(),
        cbar_kws={"label": "z-score", "shrink": 0.3},
        linewidths=0, ax=ax,
    )
    boundaries: List[int] = []
    last = None
    for i, lbl in enumerate(col_groups):
        if lbl != last:
            boundaries.append(i)
            last = lbl
    boundaries.append(n_cols)
    for b in boundaries[1:-1]:
        ax.axvline(b, color="black", linewidth=0.4)
    ax.set_xlabel(xlabel, fontsize=7)
    ax.set_ylabel("")
    ax.tick_params(axis="x", labelsize=5, rotation=90)
    ax.tick_params(axis="y", labelsize=5)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------- canonical heatmap ------------------------------------------------

def render_heatmap_canonical(
    label_per_cell: pd.Series,
    Xn: sp.csr_matrix,
    gene_index: pd.Series,
    label_to_markers: Dict[str, List[str]],
    label_order: List[str],
    out_pdf: Path,
    out_csv: Path,
    include_artifact: bool,
    artifact_labels: Set[str],
    artifact_markers: Dict[str, List[str]] | None = None,
) -> None:
    present = set(label_per_cell.unique())
    eligible = [l for l in label_order if l in present]
    # FLEX-novel labels (not in V1 yaml but present in this annotation): append after V1 order
    extra_markers = artifact_markers or {}
    for l in sorted(present):
        if l in eligible or l in artifact_labels:
            continue
        if l in extra_markers:
            eligible.append(l)
    if include_artifact and artifact_markers:
        for l in artifact_labels:
            if l in present and l in artifact_markers and l not in eligible:
                eligible.append(l)
    eligible = [l for l in eligible if include_artifact or l not in artifact_labels]
    if not eligible:
        log.warning("[canonical] no labels (include_artifact=%s)", include_artifact)
        return
    # Hybrid markers: V1 canonical/identity FIRST (preserves V1 semantic), then
    # limma top-marker supplements deduped against earlier-block claims. Ensures
    # labels like LASP-basal don't end up with only 2 columns when KRT5/14/17
    # get block-claimed by BMYO-basal first; ensures FLEX-novel labels (LHS-apocrine)
    # get a real column block since V1 yaml has no entry for them.
    combined_markers: Dict[str, List[str]] = {}
    for l in eligible:
        merged: List[str] = []
        seen: Set[str] = set()
        for g in (label_to_markers.get(l) or []):
            if g and g not in seen:
                seen.add(g); merged.append(g)
        if extra_markers and l in extra_markers:
            for g in extra_markers[l]:
                if g and g not in seen:
                    seen.add(g); merged.append(g)
        if merged:
            combined_markers[l] = merged
    block = column_block_order({l: combined_markers[l] for l in eligible if l in combined_markers}, eligible)
    block = [(g, lbl) for g, lbl in block if g in gene_index.index]
    if not block:
        return
    genes = [g for g, _ in block]
    col_groups = [lbl for _, lbl in block]
    gene_idx = np.array([int(gene_index[g]) for g in genes])

    mask = label_per_cell.isin(eligible).values
    means = cluster_mean_log1p(Xn[mask], label_per_cell[mask].values, gene_idx)
    means.columns = genes
    means = means.loc[[l for l in eligible if l in means.index]]
    Z = (means - means.mean(axis=0)) / means.std(axis=0).replace(0, 1)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    Z.to_csv(out_csv, index_label="label")
    render_block_heatmap(
        Z, col_groups, out_pdf,
        xlabel="V1 canonical + identity markers (block-ordered by first declaring label)",
    )
    log.info("[canonical%s] %d labels × %d markers -> %s",
             " +artifact" if include_artifact else "", *Z.shape, out_pdf.name)


# ---------- top markers heatmap (FLEX-specific aggregation) ------------------

def label_origin_clusters(yml: dict) -> Dict[str, List[Tuple[str, int]]]:
    """Build label → [(resolution_str, cluster_int)] map honoring fine_cluster_overrides.

    Anchor clusters contribute (anchor_res, cnum) for their assigned_label; each override
    contributes (override_res, override_cluster) for its label. A label can have origins
    at multiple resolutions (e.g. Mono_NC at anchor c21 + 3.0_c22)."""
    anchor_res = str(yml["anchor_resolution"])
    out: Dict[str, List[Tuple[str, int]]] = {}
    for ckey, body in (yml["clusters_by_resolution"][anchor_res] or {}).items():
        if not isinstance(body, dict):
            continue
        cnum = int(str(ckey).lstrip("c"))
        lbl = body.get("assigned_label", "PENDING")
        out.setdefault(lbl, []).append((anchor_res, cnum))
    for group_body in (yml.get("fine_cluster_overrides") or {}).values():
        if not isinstance(group_body, dict):
            continue
        for ov in (group_body.get("overrides") or []):
            res = str(ov["resolution"])
            cnum = int(ov["cluster"])
            lbl = str(ov["label"])
            out.setdefault(lbl, []).append((res, cnum))
    return out


def derive_top_markers_per_label(
    label_origin: Dict[str, List[Tuple[str, int]]],
    limma_dir: Path,
    compartment: str,
    top_n: int,
    candidate_labels: Set[str],
) -> Dict[str, List[str]]:
    """For each label, pull limma markers from every (resolution, cluster) origin tuple,
    aggregate by max -log10(padj), take top N. Honors fine_cluster_overrides — Treg from
    3.0_c20 gets its 3.0 limma even though anchor c7 maps to CD4_Th_like."""
    limma_cache: Dict[str, pd.DataFrame] = {}
    def get_limma(res: str) -> pd.DataFrame | None:
        if res in limma_cache:
            return limma_cache[res]
        path = limma_dir / f"limma_markers_{compartment}_leiden_{res}.csv"
        if not path.exists():
            log.warning("Limma CSV missing: %s", path)
            limma_cache[res] = None  # type: ignore
            return None
        df = pd.read_csv(path)
        df["cluster"] = df["cluster"].astype(int)
        df["gene"] = df["gene"].astype(str).str.strip('"')
        df["score"] = -np.log10(df["padj"].clip(lower=1e-300))
        limma_cache[res] = df
        return df

    out: Dict[str, List[str]] = {}
    for label in candidate_labels:
        origins = label_origin.get(label, [])
        if not origins:
            continue
        frames: List[pd.DataFrame] = []
        for res, cnum in origins:
            df = get_limma(res)
            if df is None:
                continue
            sub = df[(df["cluster"] == cnum) & (df["log2FoldChange"] > 0) & (df["padj"] < 0.05)]
            if not sub.empty:
                frames.append(sub[["gene", "score"]])
        if not frames:
            continue
        pooled = pd.concat(frames, ignore_index=True)
        gene_score = pooled.groupby("gene")["score"].max().sort_values(ascending=False)
        out[label] = gene_score.head(top_n).index.tolist()
    return out


def render_heatmap_top(
    label_per_cell: pd.Series,
    Xn: sp.csr_matrix,
    gene_index: pd.Series,
    top_markers: Dict[str, List[str]],
    label_order: List[str],
    out_pdf: Path,
    out_csv: Path,
    include_artifact: bool,
    artifact_labels: Set[str],
) -> None:
    present = set(label_per_cell.unique())
    extended_order = list(label_order)
    # FLEX-novel labels with top_markers: append after V1 order
    for l in sorted(present):
        if l in extended_order or l in artifact_labels:
            continue
        if l in top_markers:
            extended_order.append(l)
    if include_artifact:
        for l in sorted(artifact_labels):
            if l not in extended_order:
                extended_order.append(l)
    eligible = [
        l for l in extended_order
        if l in present and l in top_markers
        and (include_artifact or l not in artifact_labels)
    ]
    if not eligible:
        log.warning("[top] no labels (include_artifact=%s)", include_artifact)
        return
    block = column_block_order({l: top_markers[l] for l in eligible}, eligible)
    block = [(g, lbl) for g, lbl in block if g in gene_index.index]
    if not block:
        return
    genes = [g for g, _ in block]
    col_groups = [lbl for _, lbl in block]
    gene_idx = np.array([int(gene_index[g]) for g in genes])

    mask = label_per_cell.isin(eligible).values
    means = cluster_mean_log1p(Xn[mask], label_per_cell[mask].values, gene_idx)
    means.columns = genes
    means = means.loc[[l for l in eligible if l in means.index]]
    Z = (means - means.mean(axis=0)) / means.std(axis=0).replace(0, 1)

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    Z.to_csv(out_csv, index_label="label")
    render_block_heatmap(
        Z, col_groups, out_pdf,
        xlabel="top one-vs-rest markers per label (aggregated from anchor-resolution cluster limma)",
    )
    log.info("[top%s] %d labels × %d markers -> %s",
             " +artifact" if include_artifact else "", *Z.shape, out_pdf.name)


# ---------- label UMAP -------------------------------------------------------

def categorical_palette(n: int) -> List[tuple]:
    base = (
        list(colormaps["tab20"].colors)
        + list(colormaps["tab20b"].colors)
        + list(colormaps["tab20c"].colors)
    )
    if n <= len(base):
        return base[:n]
    extra = colormaps["hsv"](np.linspace(0, 1, n - len(base) + 1))[:-1]
    return base + [tuple(c[:3]) for c in extra]


def render_label_umap(label_per_cell: pd.Series, umap: pd.DataFrame, out_pdf: Path) -> None:
    log.info("Rendering label UMAP")
    df = umap.copy()
    df["label"] = label_per_cell.values
    labels = sorted(df["label"].dropna().unique())
    palette = categorical_palette(len(labels))
    color_map = dict(zip(labels, palette))

    fig, ax = plt.subplots(figsize=(10, 10))
    colors = df["label"].map(color_map).tolist()
    ax.scatter(df["umap_1"], df["umap_2"], c=colors, s=1.5, alpha=0.85,
               rasterized=True, linewidths=0)

    label_pos: List[Tuple[str, float, float, float, float]] = []
    for lbl in labels:
        cells = df[df["label"] == lbl]
        if len(cells) < 50:
            continue
        x = float(cells["umap_1"].median())
        y = float(cells["umap_2"].median())
        label_pos.append((lbl, x, y, x, y))

    use_adjust = False
    try:
        from adjustText import adjust_text  # type: ignore
        use_adjust = True
    except ImportError:
        log.info("adjustText not available; using internal force-directed repulsion")
        anchors = np.array([(x, y) for _, x, y, _, _ in label_pos])
        positions = anchors.copy()
        x_range = float(df["umap_1"].max() - df["umap_1"].min())
        y_range = float(df["umap_2"].max() - df["umap_2"].min())
        diag = (x_range ** 2 + y_range ** 2) ** 0.5
        min_sep = diag * 0.045
        anchor_pull = 0.05
        repel_strength = diag * 0.0035
        for _ in range(120):
            disp = np.zeros_like(positions)
            for i in range(len(positions)):
                for j in range(len(positions)):
                    if i == j:
                        continue
                    d = positions[i] - positions[j]
                    dist = float(np.linalg.norm(d))
                    if dist < min_sep:
                        if dist < 1e-6:
                            d = np.array([1.0, 1.0]) * 1e-3
                            dist = float(np.linalg.norm(d))
                        push = (min_sep - dist) / min_sep
                        disp[i] += (d / dist) * push * repel_strength
                disp[i] += (anchors[i] - positions[i]) * anchor_pull
            positions += disp
        for i, (lbl, ax_x, ax_y, _, _) in enumerate(label_pos):
            label_pos[i] = (lbl, ax_x, ax_y, float(positions[i, 0]), float(positions[i, 1]))

    texts = []
    for lbl, ax_x, ax_y, lx, ly in label_pos:
        if (lx != ax_x) or (ly != ax_y):
            ax.plot([ax_x, lx], [ax_y, ly], color="gray", lw=0.25, alpha=0.5, zorder=9)
        t = ax.text(
            lx, ly, lbl,
            fontsize=7, ha="center", va="center", color="black", weight="bold",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="black", lw=0.5, alpha=0.92),
            zorder=10,
        )
        texts.append(t)
    if use_adjust:
        from adjustText import adjust_text  # type: ignore
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="gray", lw=0.25, alpha=0.6),
            expand_text=(1.05, 1.1), expand_points=(1.05, 1.1),
            force_text=(0.2, 0.3), force_points=(0.05, 0.1),
            only_move={"text": "xy"}, lim=80,
        )

    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_linewidth(0.3)
    plt.tight_layout()
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    log.info("  wrote %s", out_pdf)


# ---------- main -------------------------------------------------------------

def main() -> None:
    args = parse_args()

    out_dir = args.out_dir / "annotation"
    out_dir.mkdir(parents=True, exist_ok=True)

    log.info("Loading DRAFT yaml: %s", args.yaml_draft)
    yml = yaml.safe_load(args.yaml_draft.read_text())

    log.info("Loading clusters: %s", args.clusters_csv)
    clusters = pd.read_csv(args.clusters_csv)

    log.info("Applying yaml -> per-cell labels")
    label_per_cell = apply_yaml(yml, clusters)
    label_summary = label_per_cell.value_counts()
    log.info("Labels: %d unique, %d cells\n%s",
             label_summary.shape[0], label_per_cell.shape[0], label_summary.to_string())
    label_summary.rename("n_cells").to_csv(out_dir / "label_summary.csv", header=True)

    cells_tsv = args.bundle_dir / "cells.tsv"
    genes_tsv = args.bundle_dir / "genes.tsv"
    mtx = args.bundle_dir / "counts.mtx.gz"
    umap_csv = args.bundle_dir / "umap.csv"
    cells = pd.read_csv(cells_tsv, sep="\t", header=None, names=["cell_id"])
    genes = pd.read_csv(genes_tsv, sep="\t", header=None, names=["gene"])
    log.info("Cells=%d, Genes=%d", len(cells), len(genes))

    umap = pd.read_csv(umap_csv)
    umap.columns = [c.lower() for c in umap.columns]
    if "cell_id" in umap.columns:
        umap = umap.set_index("cell_id").reindex(cells["cell_id"]).reset_index()
    if not {"umap_1", "umap_2"}.issubset(umap.columns):
        cs = [c for c in umap.columns if c != "cell_id"]
        umap = umap.rename(columns={cs[0]: "umap_1", cs[1]: "umap_2"})

    X = load_mtx_gz(mtx)
    if X.shape != (len(cells), len(genes)):
        if X.shape == (len(genes), len(cells)):
            X = X.T.tocsr()
        else:
            sys.exit(f"ERROR: mtx shape {X.shape}")
    log.info("Normalizing (target_sum=%s)", args.target_sum)
    Xn = normalize_log1p(X, args.target_sum)
    gene_index = pd.Series(np.arange(len(genes)), index=genes["gene"].values)
    if gene_index.index.has_duplicates:
        gene_index = gene_index[~gene_index.index.duplicated(keep="first")]

    artifact_labels = {l for l in label_per_cell.unique() if str(l).startswith("ARTIFACT")}
    log.info("Artifact labels: %s", sorted(artifact_labels))

    label_to_canonical, label_order = collect_canonical_markers(args.v1_yaml, args.max_markers_per_label)

    render_label_umap(label_per_cell, umap, out_dir / "label_umap.pdf")

    candidate_labels = set(label_per_cell.unique()) - {"PENDING"}
    label_origin = label_origin_clusters(yml)
    log.info("Label origins (label → [(res, cluster)]): %s",
             {k: v for k, v in sorted(label_origin.items())})
    top_markers = derive_top_markers_per_label(
        label_origin, args.limma_dir, args.compartment, args.top_n_per_label, candidate_labels
    )
    # extra_markers = limma top markers for EVERY label (V1, FLEX-novel, artifacts).
    # render_heatmap_canonical merges V1 yaml markers + limma top per label so V1 labels
    # whose canonical markers get block-claimed by earlier labels (e.g. LASP-basal losing
    # KRT5/14/17 to BMYO-basal) still get a meaningful column block from limma.
    extra_markers = dict(top_markers)

    render_heatmap_canonical(
        label_per_cell, Xn, gene_index, label_to_canonical, label_order,
        out_dir / "heatmap_canonical.pdf", out_dir / "heatmap_canonical.csv",
        include_artifact=False, artifact_labels=artifact_labels,
        artifact_markers=extra_markers,
    )
    render_heatmap_canonical(
        label_per_cell, Xn, gene_index, label_to_canonical, label_order,
        out_dir / "heatmap_canonical_with_artifact.pdf",
        out_dir / "heatmap_canonical_with_artifact.csv",
        include_artifact=True, artifact_labels=artifact_labels,
        artifact_markers=extra_markers,
    )

    render_heatmap_top(
        label_per_cell, Xn, gene_index, top_markers, label_order,
        out_dir / "heatmap_top.pdf", out_dir / "heatmap_top.csv",
        include_artifact=False, artifact_labels=artifact_labels,
    )
    render_heatmap_top(
        label_per_cell, Xn, gene_index, top_markers, label_order,
        out_dir / "heatmap_top_with_artifact.pdf",
        out_dir / "heatmap_top_with_artifact.csv",
        include_artifact=True, artifact_labels=artifact_labels,
    )

    log.info("Done. Outputs at %s", out_dir)


if __name__ == "__main__":
    main()
