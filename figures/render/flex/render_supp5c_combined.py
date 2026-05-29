"""render_supp5c_combined.py — Fig 3 Supp 5c (revision v6, 2026-05-07): single
merged canonical+artifact heatmap covering all 43 L2S labels across the 3
compartments.

Substrate: per-compartment counts.mtx.gz + cells.tsv + genes.tsv at
publication/preprocessing/chromium_spatial/outputs/integration_intermediate/
{Comp}/scvi_n100/. Joined with publication/analysis/annotation/flex/outputs/
flex_l2s_labels.csv. Marker selection: top-N markers per compartment by
within-compartment F-stat proxy across non-artifact labels, drawn from the
YAML × matrix intersection of V1 canonical+identity markers, no per-label
argmax constraint. Plus 8 artifact-discriminating markers.

v6 changes vs v5:
- Within-block column sort revised: STEP-diagonal (cols grouped by argmax-row
  then sorted by col-max desc within each row's group) instead of v5's
  "first-argmaxer per col + descending col_max tail". Each row now receives
  a column-block spanning the full marker width of its compartment block,
  vs v5's N-col diagonal followed by an unsorted tail.

v5 was: N=40 per compartment, within-block diagonal sort (v5 column algorithm),
        ARTIFACT_lowqc compartment-tagged, 370×125mm.
v4 was: N=25 per compartment, artifact-at-end ordering, standalone 260×100mm.
v3 was: N=17 per compartment, artifact-within-block ordering, standalone 173×88mm.

Compute:
  1. For each compartment, stream counts.mtx.gz + log-normalize (target_sum=1e4).
  2. Per L2S label, mean log1p expression across cells for the selected markers.
  3. Concatenate canonical (Epi → Imm → Str) then ARTIFACT (Epi → Imm → Str).
  4. Z-score PER MARKER ACROSS ALL 43 LABELS (cross-compartment normalization).
  5. Render row-block-grouped heatmap with separator between canonical block
     and ARTIFACT block (v4) plus thin separators between Epi/Imm/Str within
     the canonical block. No compartment color strip (v3 change retained).

Pattern reference:
- Per-cell log1p + group-mean compute: publication/analysis/annotation/flex/
  scripts/06_render_annotation.py:68-127 (load_mtx_gz, normalize_log1p,
  cluster_mean_log1p).
- Heatmap styling (RdBu_r centered at 0, vmin=-2, vmax=+2, block separator
  lines): same script lines 165-196 (render_block_heatmap).

No on-plot title, panel letter, or method caption (Illustrator handles those).

Outputs:
    supp5c_flex_canonical_artifact_combined.pdf
    supp5c_flex_canonical_artifact_combined_data.csv     (43 × N_markers, z-scored)
    supp5c_flex_canonical_artifact_combined_means.csv    (43 × N_markers, raw log1p means)
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
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
log = logging.getLogger("supp5c_combined")


COMPARTMENTS = ("Epithelial", "Immune", "Stromal")
COMP_SHORT = {"Epithelial": "epi", "Immune": "imm", "Stromal": "str"}
COMP_DISPLAY = {"Epithelial": "Epi", "Immune": "Imm", "Stromal": "Str"}

# Compartment color strip — pulled from aesthetics.yaml palette at runtime.

# Provisional marker set per CP_supp5c_marker_set.md (2026-05-07).
# Coordinator may revise; updating MARKER_SETS_DEFAULT updates this script.
MARKER_SETS_DEFAULT = {
    "Epithelial": ["PGR", "PTN", "ESR1", "ALDH1A3", "ELF5", "AGR3", "OXTR"],
    "Immune": ["FGFBP2", "FOXP3", "LILRA4", "KLRF1", "CTSG", "CPA3", "MS4A2"],
    "Stromal": ["CCL21", "TFF3", "PROX1", "SELE", "ACKR1", "SELP", "CA4"],
}
ARTIFACT_MARKERS_DEFAULT = ["HBA2", "HBB", "KRT8", "EPCAM", "ACTA2", "COL1A1", "DCN", "MT-ND4"]


def _setup_aesthetics(project_root: Path) -> dict:
    sys.path.insert(0, str(project_root / "publication" / "config"))
    from load_aesthetics import (  # type: ignore
        load_aesthetics, get_matplotlib_theme, get_palette,
    )
    cfg = load_aesthetics()
    plt.rcParams.update(get_matplotlib_theme(cfg))
    return cfg


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-root", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument(
        "--selection-mode", choices=["fstat", "diversified", "default"], default="fstat",
        help="fstat (v3, default): top-N markers by within-compartment F-stat proxy "
             "(z-matrix range across bio labels), no per-label argmax cap. "
             "diversified (v2): greedy per-label argmax in YAML×matrix intersection. "
             "default (v1): hard-coded top 7 per compartment.",
    )
    p.add_argument("--n-per-compartment", type=int, default=40,
                   help="For --selection-mode fstat: top-N markers per compartment "
                        "(v5 default 40 → 120 bio + ~6 artifact ≈ 126 total; spec 40).")
    p.add_argument("--artifact-at-end", action="store_true", default=True,
                   help="v4 default: place ALL ARTIFACT rows at the very bottom of the "
                        "heatmap, after the 37 canonical rows (Epi → Imm → Str). v3 "
                        "behavior (artifact rows last within each compartment block) "
                        "is available with --no-artifact-at-end.")
    p.add_argument("--no-artifact-at-end", dest="artifact_at_end", action="store_false")
    p.add_argument("--gap-after-canonical", action="store_true", default=True,
                   help="Draw a black separator line between the 37 canonical rows and "
                        "the ARTIFACT block. Only meaningful with --artifact-at-end.")
    p.add_argument("--diagonal-sort", action="store_true", default=True,
                   help="v5 default: apply within-compartment-block diagonal sort to "
                        "rows AND columns (rows by argmax-col then -row_max; cols by "
                        "first-row argmax then -col_max). ARTIFACT rows NOT sorted "
                        "(kept in compartment order). Use --no-diagonal-sort for v4 mode.")
    p.add_argument("--no-diagonal-sort", dest="diagonal_sort", action="store_false")
    p.add_argument("--target-sum", type=float, default=1e4)
    return p.parse_args()


def load_mtx_gz(path: Path, compartment: str | None = None) -> sp.csr_matrix:
    """Read a gzipped MatrixMarket file. Pattern: 06_render_annotation.py:68-73.

    On macOS CRSP mounts, multi-hundred-MB files often present as 0-byte
    locally and bulk reads fail with Errno 5 mid-stream. The script staged
    the 3 compartment MTX files to /tmp/supp5c_mtx_cache/{Comp}_counts.mtx.gz
    via scp before invocation; check there first and fall back to the CRSP
    path only if the cache is absent.
    """
    cache = Path(f"/tmp/supp5c_mtx_cache/{compartment}_counts.mtx.gz")
    if compartment is not None and cache.exists() and cache.stat().st_size > 1024:
        log.info("  using local /tmp cache: %s (%d bytes)", cache, cache.stat().st_size)
        path = cache
    with gzip.open(path, "rb") as fh:
        m = sio.mmread(fh)
    if not sp.issparse(m):
        m = sp.csr_matrix(m)
    return m.tocsr()


def normalize_log1p(X: sp.csr_matrix, target_sum: float) -> sp.csr_matrix:
    """Pattern: 06_render_annotation.py:76-82."""
    sums = np.asarray(X.sum(axis=1)).ravel()
    sums[sums == 0] = 1.0
    scale = target_sum / sums
    Xn = X.multiply(scale[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def collect_canonical_yaml_markers(yaml_path: Path) -> Dict[str, List[str]]:
    """Per-label canonical+identity markers from V1 marker yaml."""
    cfg = yaml.safe_load(yaml_path.read_text())
    out: Dict[str, List[str]] = {}
    for label, body in (cfg.get("labels") or {}).items():
        if not isinstance(body, dict):
            continue
        ms = (body.get("canonical_markers") or []) + (body.get("identity_markers") or [])
        out[label] = list(dict.fromkeys(ms))
    return out


def select_markers_diversified(
    matrix_csv: Path, yaml_markers: Dict[str, List[str]],
) -> List[str]:
    """Greedy per-label argmax in z-scored canonical+artifact matrix.

    For each non-ARTIFACT label, pick its strongest canonical+identity marker
    (max z) that hasn't already been claimed. Ensures every label has at least
    one column where it is the argmax. Used in v2 — left in place for archival
    reproducibility but not invoked by v3.
    """
    Z = pd.read_csv(matrix_csv, index_col=0)
    bio = Z[~Z.index.str.startswith("ARTIFACT")]
    yaml_pool = set()
    for ms in yaml_markers.values():
        yaml_pool.update(ms)
    candidates = [m for m in Z.columns if m in yaml_pool]
    chosen: List[str] = []
    used = set()
    for label in bio.index:
        # Within yaml-canonical markers, sort by this label's z value.
        # Pick the strongest still-available column whose argmax is this label.
        scores = bio.loc[label, candidates].sort_values(ascending=False)
        for m in scores.index:
            if m in used:
                continue
            # Confirm this label is argmax for m (otherwise skip).
            if bio[m].idxmax() != label:
                continue
            chosen.append(m)
            used.add(m)
            break
    log.info("  diversified pick (%d labels → %d cols): %s",
             len(bio), len(chosen), chosen)


def select_markers_top_fstat(
    matrix_csv: Path, yaml_markers: Dict[str, List[str]], n_top: int,
) -> List[str]:
    """Top-N markers by within-compartment F-stat proxy across non-artifact labels.

    Uses the per-compartment z-scored canonical+artifact matrix as the score
    substrate. For z-scored input each marker has variance ≈ 1 across the full
    label set, so variance is uninformative; range across BIO labels (max - min)
    is the discriminative signal — high range markers separate at least one
    label strongly from the rest, low range markers are uniform.

    Restricted to the YAML × matrix intersection so the chosen markers are
    V1-validated identity markers.

    No per-label argmax constraint (v3 change vs v2 diversified). All top-N
    markers regardless of which label they peak for.
    """
    Z = pd.read_csv(matrix_csv, index_col=0)
    bio = Z[~Z.index.str.startswith("ARTIFACT")]
    yaml_pool = set()
    for ms in yaml_markers.values():
        yaml_pool.update(ms)
    candidates = [m for m in Z.columns if m in yaml_pool]
    score = (bio[candidates].max(axis=0) - bio[candidates].min(axis=0))
    score = score.sort_values(ascending=False)
    chosen = score.head(n_top).index.tolist()
    log.info("  F-stat-proxy pick (top %d of %d candidates, %d bio labels): %s",
             n_top, len(candidates), len(bio), chosen)
    return chosen


def diagonal_sort_block(
    Z: pd.DataFrame, row_labels: List[str], col_labels: List[str],
) -> tuple:
    """Within-block step-diagonal sort (v6).

    Step 1 — Row order: rows ordered by argmax-column position then by row-max
    descending (1.S3 pattern, iHBCA_V1/.../1S3_confusion_matrices/scripts/
    render_1S3_confusion_matrices.py:82-97).

    Step 2 — Column order: for each column, find its argmax row in the
    diagonal-sorted row sequence. Group columns by their argmax-row position;
    within each row's column-group, sort columns by column-max descending.
    Each row receives a column-block of all markers that argmax to it,
    yielding a step-diagonal cascade across the FULL marker width (vs v5's
    "first-argmaxer per col" which produced an N-col diagonal then a tail).

    Operates on Z[row_labels, col_labels]; returns (ordered_row_labels,
    ordered_col_labels). Rows/cols outside this block are unaffected.
    """
    sub = Z.loc[row_labels, col_labels]
    if sub.empty:
        return row_labels, col_labels

    # Step 1: row order
    col_to_pos = {c: i for i, c in enumerate(col_labels)}
    row_argmax = sub.idxmax(axis=1).map(col_to_pos)
    row_max = sub.max(axis=1)
    row_order = (
        pd.DataFrame({"argmax": row_argmax, "neg_max": -row_max}, index=sub.index)
          .sort_values(by=["argmax", "neg_max"]).index.tolist()
    )

    # Step 2: column order — group cols by argmax-row in row_order
    row_to_pos = {r: i for i, r in enumerate(row_order)}
    col_argmax_row = sub.idxmax(axis=0)   # for each col, which row is argmax
    col_max = sub.max(axis=0)
    col_order = sorted(
        col_labels,
        key=lambda c: (row_to_pos[col_argmax_row[c]], -float(col_max[c])),
    )
    return row_order, col_order


def per_label_means(
    bundle_dir: Path, l2s: pd.DataFrame, marker_panel: List[str],
    target_sum: float, compartment: str | None = None,
) -> pd.DataFrame:
    """Compute per-label mean log1p for the marker panel within one compartment."""
    cells_tsv = bundle_dir / "cells.tsv"
    genes_tsv = bundle_dir / "genes.tsv"
    mtx = bundle_dir / "counts.mtx.gz"
    cells = pd.read_csv(cells_tsv, sep="\t", header=None, names=["cell_id"])
    genes = pd.read_csv(genes_tsv, sep="\t", header=None, names=["gene"])
    log.info("  cells=%d genes=%d", len(cells), len(genes))

    # cell_id → l2s_label
    j = cells.merge(l2s, on="cell_id", how="left")
    if j["l2s_label"].isna().any():
        n_miss = int(j["l2s_label"].isna().sum())
        log.warning("  %d cells missing L2S label (will be dropped)", n_miss)

    # Gene index for the markers we want
    gene_to_pos = pd.Series(np.arange(len(genes)), index=genes["gene"].astype(str))
    if gene_to_pos.index.has_duplicates:
        gene_to_pos = gene_to_pos[~gene_to_pos.index.duplicated(keep="first")]
    available = [m for m in marker_panel if m in gene_to_pos.index]
    missing = [m for m in marker_panel if m not in gene_to_pos.index]
    if missing:
        log.warning("  markers missing from this compartment's panel: %s", missing)
    gene_idx = np.array([int(gene_to_pos[m]) for m in available])

    log.info("  loading counts.mtx.gz: %s", mtx)
    X = load_mtx_gz(mtx, compartment=compartment)
    if X.shape == (len(genes), len(cells)):
        X = X.T.tocsr()
    elif X.shape != (len(cells), len(genes)):
        sys.exit(f"FATAL: mtx shape {X.shape} mismatch with cells/genes")
    log.info("  normalizing (target_sum=%.0f)", target_sum)
    Xn = normalize_log1p(X, target_sum)
    Xsub = Xn[:, gene_idx]   # cells × len(available)

    # Per-label mean. Drop cells with missing labels.
    label_arr = j["l2s_label"].values
    mask = ~pd.isna(label_arr)
    Xsub = Xsub[mask]
    label_arr = label_arr[mask]

    df = pd.DataFrame(Xsub.toarray(), columns=available)
    df["__label"] = label_arr
    means = df.groupby("__label").mean()
    # Reindex to the full marker panel; missing markers fill with NaN.
    means = means.reindex(columns=marker_panel)
    return means


def build_compartment_color_strip(
    row_labels: List[str], compartment_of: Dict[str, str],
    palette: Dict[str, str],
) -> np.ndarray:
    """RGB color strip array for the row-axis compartment band."""
    from matplotlib.colors import to_rgb
    out = np.zeros((len(row_labels), 1, 3))
    for i, lbl in enumerate(row_labels):
        comp = compartment_of[lbl]
        out[i, 0] = to_rgb(palette[comp])
    return out


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    cfg = _setup_aesthetics(project_root)

    sys.path.insert(0, str(project_root / "publication" / "config"))
    from load_paths import load_paths, resolve_path  # type: ignore
    paths = load_paths(str(project_root / "publication" / "config"))

    # ---- Step 1: substrate paths ----
    flex_root = Path(resolve_path("flex_integration_intermediate", paths))
    l2s_csv = Path(resolve_path("flex_l2s_labels", paths))
    yaml_root = project_root / "publication" / "analysis" / "annotation" / "yamls"
    matrix_root = project_root / "publication" / "analysis" / "annotation" / "flex" / "outputs" / "per_compartment"

    yaml_files = {
        "Epithelial": yaml_root / "annotation_v2_epi.yaml",
        "Immune": yaml_root / "annotation_v2_imm.yaml",
        "Stromal": yaml_root / "annotation_v2_str.yaml",
    }
    matrix_files = {
        c: matrix_root / c / "annotation/heatmap_canonical_with_artifact.csv"
        for c in COMPARTMENTS
    }

    # ---- Step 2: marker selection ----
    per_comp_markers: Dict[str, List[str]] = {}
    if args.selection_mode == "fstat":
        log.info("F-stat marker selection (top-%d per compartment, no per-label cap)",
                 args.n_per_compartment)
        for c in COMPARTMENTS:
            yaml_markers = collect_canonical_yaml_markers(yaml_files[c])
            picks = select_markers_top_fstat(matrix_files[c], yaml_markers, args.n_per_compartment)
            per_comp_markers[c] = picks
    elif args.selection_mode == "diversified":
        log.info("Diversified marker selection (v2 — per-label argmax)")
        for c in COMPARTMENTS:
            yaml_markers = collect_canonical_yaml_markers(yaml_files[c])
            per_comp_markers[c] = select_markers_diversified(matrix_files[c], yaml_markers)
    else:  # default
        log.info("Default v1 marker sets (top 7 per compartment)")
        per_comp_markers = {c: list(MARKER_SETS_DEFAULT[c]) for c in COMPARTMENTS}
    artifact_markers = list(ARTIFACT_MARKERS_DEFAULT)

    # Union, preserving compartment-block order.
    marker_panel: List[str] = []
    for c in COMPARTMENTS:
        for m in per_comp_markers[c]:
            if m not in marker_panel:
                marker_panel.append(m)
    for m in artifact_markers:
        if m not in marker_panel:
            marker_panel.append(m)
    log.info("Final marker panel (%d cols): %s", len(marker_panel), marker_panel)

    # ---- Step 3: per-compartment per-label means ----
    l2s = pd.read_csv(l2s_csv, usecols=["cell_id", "compartment", "l2s_label"])
    log.info("Loaded L2S labels: %d cells", len(l2s))

    per_comp_bio: Dict[str, pd.DataFrame] = {}
    per_comp_art: Dict[str, pd.DataFrame] = {}
    compartment_of_label: Dict[str, str] = {}
    for c in COMPARTMENTS:
        log.info("Computing per-label means for %s ...", c)
        sub = l2s[l2s["compartment"] == COMP_SHORT[c]]
        bundle = flex_root / c / "scvi_n100"
        means = per_label_means(bundle, sub, marker_panel, args.target_sum, compartment=c)
        # ARTIFACT_lowqc appears in all 3 compartments. Tag with compartment
        # short name to disambiguate rows in the merged matrix and avoid the
        # pandas .loc duplicate-label expansion bug downstream. v5.
        means = means.rename(index=lambda lbl: (
            f"{lbl} ({COMP_SHORT[c]})" if str(lbl).startswith("ARTIFACT_lowqc") else lbl
        ))
        bio_rows = sorted([r for r in means.index if not str(r).startswith("ARTIFACT")])
        art_rows = sorted([r for r in means.index if str(r).startswith("ARTIFACT")])
        per_comp_bio[c] = means.loc[bio_rows]
        per_comp_art[c] = means.loc[art_rows]
        for r in means.index:
            compartment_of_label[r] = c
        log.info("  %d labels (%d bio + %d artifact)", len(means), len(bio_rows), len(art_rows))

    # ---- Step 4: assemble row order + z-score across all labels ----
    # v4 default: ALL canonical rows first (Epi → Imm → Str), then ALL
    # ARTIFACT rows at the very bottom (also Epi → Imm → Str). With
    # --no-artifact-at-end, fall back to v3 ordering (artifact within each
    # compartment block).
    if args.artifact_at_end:
        canonical_block = pd.concat([per_comp_bio[c] for c in COMPARTMENTS], axis=0)
        artifact_block = pd.concat([per_comp_art[c] for c in COMPARTMENTS], axis=0)
        full_means = pd.concat([canonical_block, artifact_block], axis=0)
        n_canonical_total = len(canonical_block)
        log.info("Row order: %d canonical (Epi→Imm→Str) → %d ARTIFACT at end",
                 n_canonical_total, len(artifact_block))
    else:
        ordered = []
        for c in COMPARTMENTS:
            ordered.append(per_comp_bio[c])
            ordered.append(per_comp_art[c])
        full_means = pd.concat(ordered, axis=0)
        n_canonical_total = sum(len(per_comp_bio[c]) for c in COMPARTMENTS)
        log.info("Row order: artifact-within-block (v3 mode)")

    log.info("Combined matrix: %d labels × %d markers", *full_means.shape)
    # Z-score per column (marker) across all labels. Mark NaN where any
    # underlying mean was NaN (e.g. marker not in panel for some compartment —
    # not expected since all panel-membership was checked at load time).
    Z = (full_means - full_means.mean(axis=0)) / full_means.std(axis=0).replace(0, 1)

    # ---- Step 4b (v5): within-block diagonal sort ----
    # For each compartment's (row block × col block) submatrix, reorder rows
    # by argmax-col index (then -row_max) and cols by first-row argmax (then
    # -col_max). Cross-block off-diagonals reveal signal leakage. ARTIFACT
    # block (rows) NOT diagonal-sorted — kept in compartment order.
    if args.diagonal_sort:
        ordered_row_labels: List[str] = []
        ordered_col_labels: List[str] = []
        seen_cols: set = set()
        for c in COMPARTMENTS:
            row_block = list(per_comp_bio[c].index)
            # Cols belonging to this compartment block AND not already placed by
            # an earlier compartment (e.g. TFF3/MYH11/CNN1 chosen by both Epi
            # and Stromal go to Epi block; Stromal's slot loses them).
            col_block = [m for m in per_comp_markers[c]
                         if m in Z.columns and m not in seen_cols]
            new_rows, new_cols = diagonal_sort_block(Z, row_block, col_block)
            log.info("  diagonal-sorted %s block: %d rows × %d cols",
                     c, len(new_rows), len(new_cols))
            ordered_row_labels.extend(new_rows)
            ordered_col_labels.extend(new_cols)
            seen_cols.update(new_cols)
        # Append ARTIFACT rows in their existing order (no diagonal sort).
        artifact_rows_in_Z = [r for r in Z.index if r not in set(ordered_row_labels)]
        # Append remaining marker columns (artifact-only set).
        artifact_cols_in_Z = [c for c in Z.columns if c not in seen_cols]
        Z = Z.loc[ordered_row_labels + artifact_rows_in_Z,
                  ordered_col_labels + artifact_cols_in_Z]
        log.info("Within-block diagonal sort applied to row + col axes "
                 "(ARTIFACT rows + cols held at end). Final shape: %d × %d",
                 *Z.shape)

    # ---- Step 5: render ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_csv_z = args.out_dir / "supp5c_flex_canonical_artifact_combined_data.csv"
    out_csv_means = args.out_dir / "supp5c_flex_canonical_artifact_combined_means.csv"
    out_pdf = args.out_dir / "supp5c_flex_canonical_artifact_combined.pdf"

    Z.to_csv(out_csv_z, index_label="l2s_label")
    full_means.to_csv(out_csv_means, index_label="l2s_label")
    log.info("Wrote z-score matrix: %s", out_csv_z)
    log.info("Wrote raw means matrix: %s", out_csv_means)

    n_rows, n_cols = Z.shape
    # v5 layout: panoramic 5c standalone, larger to accommodate ~120-126 cols.
    # Target ~370mm wide × ~125mm tall = ~3:1 aspect. The standalone PDF is the
    # canonical deliverable; composer fits-to-cell when placing into the
    # supp5_layout_v5.pdf row (173 × 92mm), which scales the panel down.
    width_in = 370.0 / 25.4
    height_in = 125.0 / 25.4

    fig, main_ax = plt.subplots(figsize=(width_in, height_in), constrained_layout=False)
    sns.heatmap(
        Z, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        cbar_kws={"label": f"z-score (across {n_rows} labels)", "shrink": 0.4},
        linewidths=0, ax=main_ax,
        xticklabels=Z.columns.tolist(), yticklabels=Z.index.tolist(),
    )

    # Row separators — placement depends on row-ordering mode.
    if args.artifact_at_end:
        # v4: separator between Epi/Imm and Imm/Str within the canonical block
        # (thin lines), then a heavier separator between canonical block and
        # ARTIFACT block at the bottom.
        cum = 0
        for c in COMPARTMENTS[:-1]:
            cum += len(per_comp_bio[c])
            main_ax.axhline(cum, color="black", linewidth=0.4, alpha=0.7)
        if args.gap_after_canonical and n_canonical_total < n_rows:
            main_ax.axhline(n_canonical_total, color="black", linewidth=0.8)
    else:
        # v3 mode: separator after each (bio + artifact) compartment block.
        cum = 0
        for c in COMPARTMENTS[:-1]:
            cum += len(per_comp_bio[c]) + len(per_comp_art[c])
            main_ax.axhline(cum, color="black", linewidth=0.6)

    # Column block separators between Epi/Imm/Str marker groups + artifact.
    col_block_sizes = [len(per_comp_markers[c]) for c in COMPARTMENTS]
    cum_c = 0
    for size in col_block_sizes[:-1]:
        cum_c += size
        main_ax.axvline(cum_c, color="black", linewidth=0.4, alpha=0.7)
    n_bio_cols = sum(col_block_sizes)
    if n_bio_cols < n_cols:
        main_ax.axvline(n_bio_cols, color="black", linewidth=0.5, linestyle="--", alpha=0.7)

    main_ax.set_xlabel("Markers (Epi | Imm | Str | artifact)", fontsize=7)
    main_ax.set_ylabel("")
    main_ax.tick_params(axis="x", labelsize=5, rotation=90)
    main_ax.tick_params(axis="y", labelsize=5)
    fig.savefig(out_pdf, bbox_inches="tight",
                dpi=int(cfg.get("rendering", {}).get("dpi", 600)))
    plt.close(fig)
    log.info("Wrote %s (%.0fmm × %.0fmm)", out_pdf, width_in * 25.4, height_in * 25.4)
    return 0


if __name__ == "__main__":
    sys.exit(main())
