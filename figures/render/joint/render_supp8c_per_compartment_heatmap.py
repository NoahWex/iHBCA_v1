"""render_supp8c_per_compartment_heatmap.py — Fig 3 Supp 8c: per-compartment
canonical+artifact marker expression heatmap.

For each compartment {Epithelial, Immune, Stromal}, computes per-(label) mean
log1p of cascade markers from Xenium nuclear pooled counts (matches the
cascade marker provenance from publication/analysis/annotation/xenium/scripts/
02_compartment_analysis.py:85-221). Z-scored per marker across labels.
Canonical rows + artifact rows; thin separator before artifact block.
Diagonal sort within canonical rows for cascade visual.

Pattern source:
- publication/figures/render/flex/render_supp5c_combined.py (per-compartment
  per-label means + diagonal sort + heatmap render). Adapted by:
    * Substrate: Xenium nuclear pooled counts (single bundle, filter to
      compartment cells via cell_annotations) instead of FLEX per-compartment
      bundles.
    * Markers: parsed from annotation_cascade.yaml vocabulary.{label}.notes
      (regex tokenization + Xenium 280-gene panel intersection) instead of
      annotation_v2_{compartment}.yaml canonical_markers list.
    * Scope: single compartment per invocation (looped 3x via wrapper)
      instead of merged-across-compartments.

No on-plot title, panel letter, or method caption (Illustrator handles those).

Outputs (per compartment):
    supp8c_heatmap_{epithelial,immune,stromal}.pdf
    supp8c_heatmap_{epithelial,immune,stromal}_data.csv  (z-scored)
    supp8c_heatmap_{epithelial,immune,stromal}_means.csv (raw log1p means)
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_config")

import matplotlib
matplotlib.use("Agg")
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
log = logging.getLogger("supp8c")


COMPARTMENT_FULL = {"epi": "Epithelial", "imm": "Immune", "str": "Stromal"}
COMPARTMENT_DIR = {"epi": "epithelial", "imm": "immune", "str": "stromal"}
COMPARTMENT_TOKEN = {
    "epi": ("compartment_epi_cell_annotations", "compartment_epi_joint_obs"),
    "imm": ("compartment_imm_cell_annotations", "compartment_imm_joint_obs"),
    "str": ("compartment_str_cell_annotations", "compartment_str_joint_obs"),
}
NOTES_STOPWORDS = {"FLEX", "Cyto", "Nuclear", "RNA", "DNA", "DC", "Treg",
                   "PCR", "FOXP", "STK", "TF", "IL"}

# Per-compartment artifact reclassifications applied at render time. BMYO-NC
# is reclassified because the cluster is sample-specific (single patient,
# single Xenium sample) and does not represent a generalizable biological
# state. The override moves BMYO-NC into the artifact row block of the
# heatmap (separated from canonical labels by a heavy hline).
ARTIFACT_OVERRIDES = {
    "epi": {"BMYO-NC"},
    "imm": set(),
    "str": set(),
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--project-root", required=True, type=Path)
    p.add_argument("--compartment", required=True, choices=["epi", "imm", "str"])
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--target-sum", type=float, default=1e4)
    p.add_argument("--diagonal-sort", action="store_true", default=True)
    p.add_argument("--no-diagonal-sort", dest="diagonal_sort", action="store_false")
    p.add_argument("--test", action="store_true",
                   help="Subsample to 5000 Xenium cells.")
    return p.parse_args()


def parse_cascade_markers(
    cascade_yaml: Path, gene_set: set,
) -> Tuple[Dict[str, List[str]], Dict[str, List[str]], List[str], List[str]]:
    """Pattern: render_supp8b_joint_heatmap_merged.py (cascade parser).

    Same regex + intersection logic. Returns canonical_markers,
    artifact_markers, final_labels, artifact_labels.
    """
    cfg = yaml.safe_load(cascade_yaml.read_text())
    vocab = cfg.get("vocabulary", {}) or {}
    final_labels = cfg.get("final_labels", []) or []
    artifact_labels = cfg.get("artifact_labels", []) or []

    canonical_markers: Dict[str, List[str]] = {}
    artifact_markers: Dict[str, List[str]] = {}
    token_re = re.compile(r"[A-Z][A-Z0-9]{1,}(?:-[A-Z0-9]+)?")

    def _resolve(notes: str | None) -> List[str]:
        if not notes:
            return []
        toks = token_re.findall(notes)
        out: List[str] = []
        seen: set = set()
        for t in toks:
            if t in NOTES_STOPWORDS or t not in gene_set or t in seen:
                continue
            out.append(t)
            seen.add(t)
        return out

    for label, body in vocab.items():
        if isinstance(body, dict):
            canonical_markers[str(label)] = _resolve(body.get("notes"))

    # Walk the entire cascade YAML for any dict node with both `label` and
    # `notes` fields. Epithelial cascades nest artifact entries under
    # `luminal_detail:` (override blocks) rather than directly under `base:`;
    # recursive walk catches both layouts so contamination markers from
    # all artifact base-cluster notes reach the marker panel.
    base_notes_by_label: Dict[str, List[str]] = {}

    def _walk(node):
        if isinstance(node, dict):
            if "label" in node and "notes" in node:
                lbl = str(node.get("label", ""))
                toks = _resolve(node.get("notes"))
                for t in toks:
                    base_notes_by_label.setdefault(lbl, []).append(t)
            for v in node.values():
                _walk(v)
        elif isinstance(node, list):
            for item in node:
                _walk(item)

    _walk(cfg)
    # Dedupe per label, preserve order.
    for lbl in list(base_notes_by_label):
        seen: set = set()
        out: List[str] = []
        for t in base_notes_by_label[lbl]:
            if t not in seen:
                out.append(t)
                seen.add(t)
        base_notes_by_label[lbl] = out

    # Canonical-marker backfill: when vocabulary notes are sparse (<2 tokens),
    # backfill from base notes. Then truncate per-label at MAX_CANONICAL_PER_LABEL
    # to keep the panel compact (otherwise large compartments with rich
    # vocabularies blow the marker count and dilute z-score variance).
    MAX_CANONICAL_PER_LABEL = 3
    for lbl, toks in canonical_markers.items():
        if len(toks) < 2 and lbl in base_notes_by_label:
            extra = [g for g in base_notes_by_label[lbl] if g not in set(toks)]
            toks = toks + extra[: 5 - len(toks)]
        canonical_markers[lbl] = toks[:MAX_CANONICAL_PER_LABEL]
    # Artifacts: pull from ALL base notes for that artifact label (no cap).
    # The contamination-marker signal is the artifact-state evidence and must
    # be visible in the heatmap. Deduplication against canonical happens at
    # the marker_panel union step; here we collect the full contaminating-
    # transcript set per artifact.
    for lbl in artifact_labels:
        if lbl in base_notes_by_label:
            artifact_markers[lbl] = list(base_notes_by_label[lbl])
        else:
            artifact_markers[lbl] = []

    return canonical_markers, artifact_markers, list(final_labels), list(artifact_labels)


# {X}_artifact in compartment Y => cells assigned compartment Y but expressing
# compartment X canonicals. To make that call evident in the heatmap we pull
# the OTHER compartment's canonical markers as additional artifact-evidence
# columns. Keys are bare artifact-label prefixes; values are compartment_short.
CROSS_COMPARTMENT_BY_ARTIFACT_NAME = {
    "Epithelial_artifact": "epi",
    "Stromal_artifact": "str",
    "Immune_artifact": "imm",
}


def load_cross_compartment_artifact_evidence(
    artifact_labels: List[str],
    current_compartment: str,
    project_root: Path,
    gene_set: set,
    max_per_artifact: int = 6,
) -> Dict[str, List[str]]:
    """For each cross-compartment-named artifact label (Immune_artifact /
    Stromal_artifact / Epithelial_artifact), pull a few top discriminating
    markers from the OTHER compartment's cascade. Capped at `max_per_artifact`
    per artifact so the artifact-evidence block stays compact and the per-
    marker z-score normalisation isn't diluted by hundreds of weak markers."""
    extra: Dict[str, List[str]] = {}
    for al in artifact_labels:
        target_short = CROSS_COMPARTMENT_BY_ARTIFACT_NAME.get(al)
        if target_short is None or target_short == current_compartment:
            continue
        target_dir = COMPARTMENT_DIR[target_short]
        cascade_yaml = (project_root / "publication" / "analysis" / "annotation"
                        / "xenium" / "outputs" / "compartment" / target_dir
                        / "annotation_cascade.yaml")
        if not cascade_yaml.exists():
            continue
        other_canonical, _, _, _ = parse_cascade_markers(cascade_yaml, gene_set)
        # Pick top markers across the other compartment's canonical labels:
        # one or two from each, prioritising the most-distinctive (vocab)
        # tokens, until we hit max_per_artifact.
        picked: List[str] = []
        seen: set = set()
        # Round-robin one token per label to spread coverage across cell types.
        for round_idx in range(5):
            for _lbl, toks in other_canonical.items():
                if round_idx < len(toks):
                    t = toks[round_idx]
                    if t not in seen:
                        picked.append(t)
                        seen.add(t)
                        if len(picked) >= max_per_artifact:
                            break
            if len(picked) >= max_per_artifact:
                break
        extra[al] = picked
    return extra


def normalize_log1p_subset(
    X: sp.csr_matrix, marker_idx: np.ndarray, target_sum: float,
) -> sp.csr_matrix:
    """Compute scale on full matrix, then normalize subset of marker columns.

    Memory-efficient: avoids creating a full normalized copy. Pattern: extracted
    from render_supp8b_joint_heatmap_merged.py per_label_means_flex.
    """
    sums = np.asarray(X.sum(axis=1)).ravel()
    sums[sums == 0] = 1.0
    scale = target_sum / sums
    X_sub = X[:, marker_idx].tocsr()
    Xn = X_sub.multiply(scale[:, None]).tocsr()
    Xn.data = np.log1p(Xn.data)
    return Xn


def diagonal_sort_block(
    Z: pd.DataFrame, row_labels: List[str], col_labels: List[str],
) -> Tuple[List[str], List[str]]:
    """Step-diagonal sort. Pattern: render_supp5c_combined.py:247-287."""
    sub = Z.loc[row_labels, col_labels]
    if sub.empty:
        return row_labels, col_labels
    col_to_pos = {c: i for i, c in enumerate(col_labels)}
    row_argmax = sub.idxmax(axis=1).map(col_to_pos)
    row_max = sub.max(axis=1)
    row_order = (
        pd.DataFrame({"argmax": row_argmax, "neg_max": -row_max}, index=sub.index)
          .sort_values(by=["argmax", "neg_max"]).index.tolist()
    )
    row_to_pos = {r: i for i, r in enumerate(row_order)}
    col_argmax_row = sub.idxmax(axis=0)
    col_max = sub.max(axis=0)
    col_order = sorted(
        col_labels,
        key=lambda c: (row_to_pos[col_argmax_row[c]], -float(col_max[c])),
    )
    return row_order, col_order


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    config_dir = str(project_root / "publication" / "config")

    sys.path.insert(0, config_dir)
    from load_aesthetics import (  # type: ignore
        load_aesthetics, get_matplotlib_theme,
    )
    from load_paths import load_paths, resolve_path  # type: ignore

    paths = load_paths(config_dir=config_dir)
    aes = load_aesthetics(config_dir=config_dir)
    plt.rcParams.update(get_matplotlib_theme(aes))

    compartment_short = args.compartment
    compartment_full = COMPARTMENT_FULL[compartment_short]
    compartment_dir = COMPARTMENT_DIR[compartment_short]

    # ---- Substrate ----
    bundle_dir = Path(resolve_path("xenium_nuclear_bundle", paths))
    cascade_yaml = (project_root / "publication" / "analysis" / "annotation"
                    / "xenium" / "outputs" / "compartment" / compartment_dir
                    / "annotation_cascade.yaml")
    ann_token, obs_token = COMPARTMENT_TOKEN[compartment_short]
    ann_csv = Path(resolve_path(ann_token, paths))
    obs_csv = Path(resolve_path(obs_token, paths))

    log.info("Compartment: %s", compartment_full)
    log.info("Cascade YAML: %s", cascade_yaml)

    # ---- Genes + cascade markers ----
    genes = pd.read_csv(bundle_dir / "xenium_genes.tsv", sep="\t",
                        header=None, names=["gene"])["gene"].tolist()
    gene_set = set(genes)
    canonical_markers, artifact_markers, final_labels, artifact_labels = \
        parse_cascade_markers(cascade_yaml, gene_set)

    # Render-time artifact override: move overridden labels from canonical to
    # artifact lists; their cascade markers still seed the marker_panel for
    # cross-row visibility.
    overrides = ARTIFACT_OVERRIDES.get(compartment_short, set())
    if overrides:
        n_pre = len(artifact_labels)
        for ovr in overrides:
            if ovr in final_labels and ovr not in artifact_labels:
                artifact_labels.append(ovr)
            if ovr in canonical_markers and ovr not in artifact_markers:
                artifact_markers[ovr] = canonical_markers[ovr]
        log.info("Override applied: %s → artifact_labels grew %d → %d",
                 sorted(overrides), n_pre, len(artifact_labels))

    log.info("Cascade: %d canonical labels (%d resolve markers), %d artifact labels",
             len(canonical_markers),
             len([v for v in canonical_markers.values() if v]),
             len(artifact_labels))

    # Marker panel: canonical (final_labels order) + artifact-discriminating.
    marker_panel: List[str] = []
    seen: set = set()
    for lbl in final_labels:
        if lbl in set(artifact_labels):
            continue
        for m in canonical_markers.get(lbl, []):
            if m not in seen:
                marker_panel.append(m)
                seen.add(m)
    for lbl in artifact_labels:
        for m in artifact_markers.get(lbl, []):
            if m not in seen:
                marker_panel.append(m)
                seen.add(m)
    log.info("Marker panel: %d markers", len(marker_panel))

    # ---- Cell metadata + filter Xenium ----
    ann = pd.read_csv(ann_csv, usecols=["cell_id", "label", "is_artifact"])
    obs = pd.read_csv(obs_csv, usecols=["cell_id", "platform"])
    meta = ann.merge(obs, on="cell_id", how="inner")
    meta["is_artifact"] = meta["is_artifact"].astype(bool)
    meta = meta[meta["platform"] == "xenium"].reset_index(drop=True)
    if overrides:
        n_pre = int(meta["is_artifact"].sum())
        meta.loc[meta["label"].isin(overrides), "is_artifact"] = True
        n_post = int(meta["is_artifact"].sum())
        log.info("Override applied to cells: +%d reclassified", n_post - n_pre)
    log.info("Xenium cells: %d (artifact=%d)",
             len(meta), int(meta["is_artifact"].sum()))

    if args.test:
        meta = meta.sample(n=min(5000, len(meta)), random_state=0).reset_index(drop=True)
        log.info("--test: subsampled to %d cells", len(meta))

    # ---- Load + slice counts ----
    log.info("Loading MTX from %s", bundle_dir / "xenium_nuclear_counts.mtx.gz")
    cells_full = pd.read_csv(bundle_dir / "xenium_cells.tsv", sep="\t",
                             header=None, names=["cell_id"])["cell_id"].tolist()
    with gzip.open(bundle_dir / "xenium_nuclear_counts.mtx.gz", "rb") as fh:
        m = sio.mmread(fh)
    if not sp.issparse(m):
        m = sp.csr_matrix(m)
    m = m.tocsr()
    if m.shape == (len(genes), len(cells_full)):
        m = m.T.tocsr()
    elif m.shape != (len(cells_full), len(genes)):
        sys.exit(f"FATAL: MTX shape {m.shape} mismatch")

    cell_to_row = pd.Series(np.arange(len(cells_full)), index=cells_full)
    if cell_to_row.index.has_duplicates:
        cell_to_row = cell_to_row[~cell_to_row.index.duplicated(keep="first")]
    matched = meta["cell_id"].map(cell_to_row).dropna().astype(int)
    n_dropped = len(meta) - len(matched)
    if n_dropped:
        log.warning("%d cells dropped (no MTX row)", n_dropped)
    meta = meta.loc[matched.index].reset_index(drop=True)

    # Normalize and aggregate across the FULL Xenium panel (all 280 genes).
    # We then z-score per gene across labels and pick artifact-discriminating
    # markers from the data rather than relying on cascade-YAML notes alone.
    # This makes the artifact-evidence columns reflect actual differential
    # expression rather than literature-curated marker mentions.
    log.info("Normalizing %d cells × %d (full panel) genes", len(matched), len(genes))
    Xn_full = normalize_log1p_subset(
        m[matched.values, :], np.arange(len(genes)), args.target_sum
    )
    df_full = pd.DataFrame(Xn_full.toarray(), columns=genes)
    df_full["__label"] = meta["label"].astype(str).values
    means_full = df_full.groupby("__label").mean(numeric_only=True)
    means_full.index = means_full.index.set_names(["label"])

    # ---- Row order: canonical (cascade order) + artifact ----
    canonical_rows = [lbl for lbl in final_labels
                      if lbl not in set(artifact_labels)
                      and lbl in set(means_full.index)]
    artifact_rows = [lbl for lbl in artifact_labels if lbl in set(means_full.index)]

    means_full = means_full.loc[canonical_rows + artifact_rows]
    # Full-panel z-score per gene across labels — basis for data-driven
    # marker selection.
    Z_full = (means_full - means_full.mean(axis=0)) / means_full.std(axis=0).replace(0, 1)

    # ---- Data-driven artifact marker selection ----
    # For each artifact label, pick the top-N genes by z-score (markers where
    # this artifact cluster expresses more strongly than the canonical labels
    # — direct evidence for the artifact call).
    TOP_N_PER_ARTIFACT = 4
    data_artifact_markers: List[str] = []
    seen: set = set()
    for art_lbl in artifact_rows:
        top = Z_full.loc[art_lbl].sort_values(ascending=False).head(TOP_N_PER_ARTIFACT)
        picked = [g for g in top.index.tolist() if g not in seen]
        for g in picked:
            seen.add(g)
            data_artifact_markers.append(g)
        log.info("Data-driven artifact markers for %s: %s",
                 art_lbl, list(top.index[:TOP_N_PER_ARTIFACT]))

    # Compose final marker_panel: canonical cascade markers + data-driven
    # artifact markers. Drop the cascade artifact markers from the panel
    # since data-driven selection supersedes them.
    final_marker_panel: List[str] = []
    seen = set()
    for lbl in canonical_rows:
        for mk in canonical_markers.get(lbl, []):
            if mk not in seen and mk in Z_full.columns:
                final_marker_panel.append(mk)
                seen.add(mk)
    for mk in data_artifact_markers:
        if mk not in seen:
            final_marker_panel.append(mk)
            seen.add(mk)
    marker_panel = final_marker_panel
    log.info("Final marker panel: %d canonical + %d artifact-evidence = %d total",
             len(seen) - len(data_artifact_markers),
             len([g for g in data_artifact_markers if g in seen]),
             len(marker_panel))

    full_means = means_full[marker_panel]
    Z = Z_full[marker_panel]
    log.info("Matrix: %d canonical + %d artifact rows × %d markers",
             len(canonical_rows), len(artifact_rows), len(marker_panel))

    # ---- Diagonal sort within canonical block (artifacts held at end) ----
    new_cols: List[str] = []
    if args.diagonal_sort and canonical_rows:
        canonical_cols = [m for m in marker_panel
                          if m in [g for lbl in canonical_rows
                                   for g in canonical_markers.get(lbl, [])]]
        new_rows, new_cols = diagonal_sort_block(Z, canonical_rows, canonical_cols)
        log.info("Diagonal-sorted canonical block: %d × %d",
                 len(new_rows), len(new_cols))
        leftover_cols = [m for m in marker_panel if m not in set(new_cols)]
        Z = Z.loc[new_rows + artifact_rows, new_cols + leftover_cols]
        full_means = full_means.loc[Z.index, Z.columns]

    # Transpose so markers index rows and labels index columns. Canonical
    # labels occupy the left columns; artifact labels are at the right.
    # Canonical markers occupy the top rows (diagonal-sorted); artifact-
    # discriminating markers (`leftover_cols` from above, pulled from
    # cascade artifact-label notes via _walk()) occupy the bottom rows.
    Zt = Z.T
    n_rows, n_cols = Zt.shape
    n_canonical = len(canonical_rows)
    n_canonical_markers = len(new_cols)

    # ---- Render ----
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = args.out_dir / f"supp8c_heatmap_{compartment_dir}.pdf"
    out_z_csv = args.out_dir / f"supp8c_heatmap_{compartment_dir}_data.csv"
    out_means_csv = args.out_dir / f"supp8c_heatmap_{compartment_dir}_means.csv"

    Z.to_csv(out_z_csv, index_label="label")
    full_means.to_csv(out_means_csv, index_label="label")

    # Per-compartment heatmap dimensions. Anisotropic cells: label columns
    # are wider than marker rows are tall, so the panel aspect (W/H) lands
    # close to 1.0–1.6 across compartments. This matches the layout
    # heatmap cell (78 × 55 mm, aspect ~1.4) and minimizes whitespace
    # after aspect-preserve scaling in compose.
    label_cell_mm = 6.0
    marker_cell_mm = 2.0
    width_in = max(40.0, 20.0 + n_cols * label_cell_mm) / 25.4
    height_in = max(50.0, 25.0 + n_rows * marker_cell_mm) / 25.4

    fig, ax = plt.subplots(figsize=(width_in, height_in))
    sns.heatmap(
        Zt, cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        cbar_kws={"label": "z-score", "shrink": 0.4},
        linewidths=0, ax=ax,
        xticklabels=Zt.columns.tolist(),
        yticklabels=Zt.index.tolist(),
    )
    # Heavy vline before artifact-label block (artifact labels on the right).
    if n_canonical < n_cols:
        ax.axvline(n_canonical, color="black", linewidth=0.8)
    # Heavy hline before artifact-marker block (artifact-discriminating
    # markers — pulled from cascade artifact-label notes — on the bottom).
    # Makes the canonical/artifact marker split visually parallel to the
    # canonical/artifact label split, so the artifact-evidence relationship
    # (artifact label column + artifact marker row → off-diagonal block) reads
    # clearly.
    if 0 < n_canonical_markers < n_rows:
        ax.axhline(n_canonical_markers, color="black", linewidth=0.8)

    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor", fontsize=5)
    ax.tick_params(axis="y", labelsize=5)
    # Bold y-tick labels for artifact-marker rows so the bottom block stands
    # out even at thumbnail size in the layout composite.
    for tick_lbl in ax.get_yticklabels()[n_canonical_markers:]:
        tick_lbl.set_fontweight("bold")

    fig.tight_layout()
    fig.savefig(out_pdf, bbox_inches="tight",
                dpi=int(aes.get("rendering", {}).get("dpi", 600)))
    plt.close(fig)
    log.info("Wrote %s (%dmm × %dmm, %d markers × %d labels)",
             out_pdf, int(width_in * 25.4), int(height_in * 25.4), n_rows, n_cols)
    return 0


if __name__ == "__main__":
    sys.exit(main())
