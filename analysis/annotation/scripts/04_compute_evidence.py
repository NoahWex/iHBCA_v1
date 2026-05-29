"""
04_compute_evidence.py — per-label evidence report with flagging.

Computes the four-axis evidence from §7.2 and emits `evidence_report.csv`. Does
NOT auto-fail on scientific predicates (D1 resolved: flagging only; user signs off).

Hard failures (§7.1): handled upstream by 06_emit_labels_csv.py. This script
assumes labels.csv is valid.

Evidence metrics per label:
  - Multi-study: n_studies, top_study, top_study_frac, cells-per-study JSON
  - Multi-patient: n_donors, top_donor_frac
  - Artifact marker enrichment: z-score per category (out_of_compartment_*, proliferation, dissociation_stress)
  - Canonical marker presence: declared, found (frac_expressed >= threshold), missing, required_gene_present

Flags (NOTE column, comma-separated):
  - single-study: label entirely one study
  - single-donor: label entirely one donor
  - possible-artifact: any artifact-category z >= threshold AND not is_artifact
  - sparse-markers: fewer than min_count declared canonical markers "present"

Signature scoring (deferred): per-label UCell scoring NOT implemented in this
script; `scoring_signature` YAML field is read but not scored. See design doc
§7 note.

Design doc: iHBCA_publication/coordination/plans/annotation_pipeline_v2.md §7
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache")

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("compute_evidence")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--v2-yaml", required=True, type=Path)
    p.add_argument("--labels-csv", required=True, type=Path, help="Output of 06_emit_labels_csv.py")
    p.add_argument("--counts-npz", required=True, type=Path)
    p.add_argument("--gene-data", required=True, type=Path)
    p.add_argument("--metadata-csv", required=True, type=Path, help="compartment metadata_enriched.csv")
    p.add_argument("--config", required=True, type=Path, help="compartment_config.yaml")
    p.add_argument("--out-csv", required=True, type=Path)
    p.add_argument("--detection-threshold", type=float, default=0.2,
                   help="Fraction of label cells that must express a canonical marker for it to count as 'present'")
    return p.parse_args()


# =============================================================================
# Helpers
# =============================================================================

def safe_z(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    if sigma <= 0:
        return np.zeros_like(x)
    return (x - mu) / sigma


def frac_expressed(X: sparse.csr_matrix, cell_mask: np.ndarray, gene_idx: int) -> float:
    """Fraction of cells in mask with gene expressed (count > 0)."""
    if cell_mask.sum() == 0:
        return 0.0
    col = X[:, gene_idx]
    expressed = (col.toarray().flatten() > 0)
    return float(expressed[cell_mask].mean())


def mean_log1p(X: sparse.csr_matrix, cell_mask: np.ndarray, gene_indices: np.ndarray) -> np.ndarray:
    """Mean log1p expression across given gene indices, one value per cell in mask."""
    if cell_mask.sum() == 0:
        return np.array([])
    sub = X[np.where(cell_mask)[0]][:, gene_indices]
    sub_log = sub.copy()
    sub_log.data = np.log1p(sub_log.data)
    return np.asarray(sub_log.mean(axis=1)).flatten()


# =============================================================================
# Per-label metrics
# =============================================================================

def evaluate_label(
    label_name: str,
    rec: dict,
    label_cells: pd.DataFrame,
    all_cells_meta: pd.DataFrame,
    X: sparse.csr_matrix,
    gene_data: pd.DataFrame,
    ensg_to_idx: dict,
    symbol_to_ensg: dict,
    artifact_cfg: dict,
    flagging: dict,
    detection_threshold: float,
) -> dict:
    n = len(label_cells)
    row: dict = {
        "label": label_name,
        "n_cells": n,
        "is_artifact": bool(rec.get("is_artifact", False)),
        "lineage": rec.get("lineage", ""),
        "mode": rec.get("mode", "cluster"),
        "source_resolution": rec.get("source_resolution"),
    }
    flags: list[str] = []
    if n == 0:
        row["NOTE"] = "empty-label"
        return row

    # --- Multi-study ---
    studies = label_cells["dataset"].value_counts()
    row["n_studies"] = int(studies.size)
    top_study = studies.index[0]
    top_study_n = int(studies.iloc[0])
    row["top_study"] = top_study
    row["top_study_frac"] = round(top_study_n / n, 4)
    row["cells_per_study"] = json.dumps(studies.to_dict())
    if flagging["multi_study"].get("flag_if_single_study") if False else top_study_n == n:
        flags.append(flagging["multi_study"]["flag_label"])

    # --- Multi-patient ---
    donors = label_cells["patientID"].value_counts()
    row["n_donors"] = int(donors.size)
    row["top_donor_frac"] = round(int(donors.iloc[0]) / n, 4) if donors.size else 0.0
    if int(donors.iloc[0]) == n if donors.size else True:
        flags.append(flagging["multi_patient"]["flag_label"])

    # --- Canonical markers ---
    canonical_declared = list(rec.get("canonical_markers") or [])
    found = []
    missing = []
    label_cell_idx = np.asarray(label_cells.index, dtype=int)
    for sym in canonical_declared:
        ensg = symbol_to_ensg.get(sym)
        if ensg is None or ensg not in ensg_to_idx:
            missing.append(sym + ":notInGeneData")
            continue
        f = frac_expressed(X, _mask_from_indices(label_cell_idx, X.shape[0]), ensg_to_idx[ensg])
        if f >= detection_threshold:
            found.append(f"{sym}:{f:.2f}")
        else:
            missing.append(f"{sym}:{f:.2f}")
    row["canonical_declared_n"] = len(canonical_declared)
    row["canonical_found_n"] = len(found)
    row["canonical_found"] = "|".join(found)
    row["canonical_missing"] = "|".join(missing)
    req = rec.get("required_gene")
    row["required_gene"] = req or ""
    if req:
        ensg = symbol_to_ensg.get(req)
        if ensg and ensg in ensg_to_idx:
            f = frac_expressed(X, _mask_from_indices(label_cell_idx, X.shape[0]), ensg_to_idx[ensg])
            row["required_gene_frac_expressed"] = round(f, 4)
        else:
            row["required_gene_frac_expressed"] = -1.0  # not in gene_data (should have been caught by 06)
    else:
        row["required_gene_frac_expressed"] = None

    if row["canonical_found_n"] < flagging["canonical_marker_count"]["min_count"]:
        if not row["is_artifact"]:
            flags.append(flagging["canonical_marker_count"]["flag_label"])

    # --- Artifact marker enrichment ---
    z_threshold = float(flagging["artifact_enrichment"]["z_threshold"])
    artifact_z_by_category: dict = {}
    any_artifact_positive = False
    label_mask = _mask_from_indices(label_cell_idx, X.shape[0])
    rest_mask = ~label_mask
    for category, symbols in artifact_cfg.items():
        gene_indices = []
        for sym in symbols:
            ensg = symbol_to_ensg.get(sym)
            if ensg and ensg in ensg_to_idx:
                gene_indices.append(ensg_to_idx[ensg])
        if not gene_indices:
            artifact_z_by_category[category] = None
            continue
        gene_indices = np.array(gene_indices, dtype=int)
        label_vals = mean_log1p(X, label_mask, gene_indices)
        rest_vals = mean_log1p(X, rest_mask, gene_indices)
        if label_vals.size == 0 or rest_vals.size == 0:
            artifact_z_by_category[category] = None
            continue
        mu = float(rest_vals.mean())
        sigma = float(rest_vals.std(ddof=0))
        z = safe_z(label_vals, mu, sigma)
        mean_z = float(z.mean())
        artifact_z_by_category[category] = round(mean_z, 3)
        if mean_z >= z_threshold:
            any_artifact_positive = True
    row["artifact_z_by_category"] = json.dumps(artifact_z_by_category)
    # Explicit per-category columns (easier to scan than JSON)
    for category, zv in artifact_z_by_category.items():
        row[f"artifact_z__{category}"] = zv
    # Composite compromise: max positive z and count of categories at/above threshold
    numeric_z = [v for v in artifact_z_by_category.values() if isinstance(v, (int, float))]
    row["artifact_max_z"] = round(max(numeric_z), 3) if numeric_z else None
    row["artifact_n_flagged_categories"] = int(sum(1 for v in numeric_z if v >= z_threshold))
    if any_artifact_positive and not row["is_artifact"]:
        flags.append(flagging["artifact_enrichment"]["flag_label"])

    # --- Signature scoring (deferred) ---
    row["scoring_signature_declared"] = bool(rec.get("scoring_signature"))
    row["signature_score_median"] = None
    row["signature_score_q25"] = None
    row["signature_score_q75"] = None

    row["NOTE"] = ",".join(flags) if flags else ""
    row["reviewer_signoff"] = ""  # user fills this in during review
    return row


def _mask_from_indices(indices: np.ndarray, n_total: int) -> np.ndarray:
    mask = np.zeros(n_total, dtype=bool)
    mask[indices] = True
    return mask


# =============================================================================
# Main
# =============================================================================

def main() -> None:
    args = parse_args()

    yml = yaml.safe_load(args.v2_yaml.read_text())
    cfg = yaml.safe_load(args.config.read_text())
    compartment = yml.get("compartment", "imm")
    artifact_cfg = cfg["compartments"][compartment]["artifact_markers"]
    flagging = cfg["pipeline"]["flagging"]

    log.info("Loading labels: %s", args.labels_csv)
    labels = pd.read_csv(args.labels_csv)

    log.info("Loading gene_data: %s", args.gene_data)
    gene_data = pd.read_csv(args.gene_data, index_col=0)
    ensg_to_idx = {ensg: i for i, ensg in enumerate(gene_data.index)}
    symbol_to_ensg = {sym: ensg for ensg, sym in gene_data["gene_symbol"].items() if isinstance(sym, str)}

    log.info("Loading metadata: %s", args.metadata_csv)
    meta = pd.read_csv(args.metadata_csv, usecols=["cell_id", "dataset", "patientID"], low_memory=False)

    log.info("Loading counts: %s", args.counts_npz)
    loaded = np.load(args.counts_npz, allow_pickle=True)
    X = sparse.csr_matrix(
        (loaded["data"], loaded["indices"], loaded["indptr"]),
        shape=tuple(loaded["shape"]),
    )
    log.info("  counts shape: %s", X.shape)

    # Join labels + metadata on cell_id; preserve row order aligned with counts matrix.
    # Counts rows align with cell order in leiden_csv / metadata_csv — enforce.
    combined = meta.merge(labels, on="cell_id", how="inner").reset_index(drop=True)
    log.info("Joined rows: %d", len(combined))
    if len(combined) != X.shape[0]:
        log.warning("joined rows %d != counts rows %d — row alignment relies on meta's native order", len(combined), X.shape[0])

    rows = []
    for label_name, rec in yml["labels"].items():
        label_cells = combined[combined["label"] == label_name]
        log.info("Evaluating %s (n=%d)", label_name, len(label_cells))
        rows.append(
            evaluate_label(
                label_name, rec, label_cells, combined, X, gene_data,
                ensg_to_idx, symbol_to_ensg, artifact_cfg, flagging,
                args.detection_threshold,
            )
        )

    out = pd.DataFrame(rows)
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    # Header with config fingerprint (§7.3)
    header_comment = (
        f"# evidence_report generated {pd.Timestamp.now().isoformat()}\n"
        f"# detection_threshold={args.detection_threshold}\n"
        f"# flagging_thresholds={json.dumps(flagging)}\n"
    )
    with args.out_csv.open("w") as f:
        f.write(header_comment)
        out.to_csv(f, index=False)
    log.info("Wrote evidence report: %s (%d labels)", args.out_csv, len(out))
    n_flagged = int((out["NOTE"] != "").sum())
    log.info("Flagged labels (NOTE non-empty): %d / %d", n_flagged, len(out))


if __name__ == "__main__":
    main()
