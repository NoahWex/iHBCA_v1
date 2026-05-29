"""
06_emit_labels_csv.py — derive labels.csv (cell_id → label) from v2 YAML.

Deterministic assignment algorithm (§6 precedence rule in design doc):
  1. Every cell assigned its cluster-mode label at the declared source resolution.
  2. fine_cluster labels override their parent_label on the cluster cells.
  3. expression-mode labels override parent_label for cells passing the signature.
  4. Parent-residue cells retain the parent label (Q1 resolved: no residue suffix).
  5. Cells fully unassigned after steps 1-4 → hard failure (§7.1 unassigned_cells).

Hard-fail conditions (§7.1):
  - required_gene_missing: label declares `required_gene` but gene not in gene_data
  - schema_violation: bad mode, missing required fields, precedence conflict
  - unassigned_cells: any cell not covered by a label at source_resolution

Output: labels.csv with cell_id, label, lineage, is_artifact
Design doc: iHBCA_publication/coordination/plans/annotation_pipeline_v2.md §6, §7.1
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("emit_labels_csv")

VALID_MODES = {"cluster", "fine_cluster", "expression"}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--v2-yaml", required=True, type=Path)
    p.add_argument("--leiden-csv", required=True, type=Path)
    p.add_argument("--gene-data", required=True, type=Path)
    p.add_argument("--counts-npz", type=Path, default=None, help="Required if any expression-mode labels")
    p.add_argument("--out-csv", required=True, type=Path)
    return p.parse_args()


# =============================================================================
# Validation
# =============================================================================

def validate_schema(yml: dict) -> None:
    if yml.get("schema_version") != 2:
        raise RuntimeError(f"schema_violation: expected schema_version=2, got {yml.get('schema_version')}")
    if "labels" not in yml or "clusters_by_resolution" not in yml:
        raise RuntimeError("schema_violation: missing top-level labels/clusters_by_resolution")

    for label_name, rec in yml["labels"].items():
        mode = rec.get("mode")
        if mode not in VALID_MODES:
            raise RuntimeError(f"schema_violation: label {label_name} has invalid mode {mode!r}")
        if mode == "cluster" or mode == "fine_cluster":
            if not rec.get("source_resolution"):
                raise RuntimeError(f"schema_violation: label {label_name} missing source_resolution")
            if not rec.get("source_cluster_ids"):
                raise RuntimeError(f"schema_violation: label {label_name} missing source_cluster_ids")
        if mode == "fine_cluster" and not rec.get("parent_label"):
            raise RuntimeError(f"schema_violation: fine_cluster label {label_name} missing parent_label")
        if mode == "expression":
            if not rec.get("parent_label"):
                raise RuntimeError(f"schema_violation: expression label {label_name} missing parent_label")
            sig = rec.get("expression_signature") or {}
            if not sig.get("required_positive"):
                raise RuntimeError(f"schema_violation: expression label {label_name} missing required_positive")


def validate_required_genes(yml: dict, gene_symbols: set[str]) -> None:
    for label_name, rec in yml["labels"].items():
        req = rec.get("required_gene")
        if req and req not in gene_symbols:
            raise RuntimeError(f"required_gene_missing: label {label_name} requires {req!r} which is not in gene_data")


# =============================================================================
# Assignment
# =============================================================================

def normalize_cluster_id(cid) -> str:
    """Strip optional 'c' prefix used in v1 YAMLs (c0, c3, ...). v2 YAMLs may use
    raw cluster integers as strings. Always return a string matching the Leiden
    column value."""
    s = str(cid)
    if s.startswith("c") and s[1:].isdigit():
        return s[1:]
    return s


def assign_cluster_mode(
    cell_df: pd.DataFrame, yml: dict, resolutions_used: list[str]
) -> pd.DataFrame:
    """Step 1: assign every cell to its cluster-mode label at source_resolution."""
    cell_df = cell_df.copy()
    cell_df["label"] = None
    cell_df["lineage"] = None
    cell_df["is_artifact"] = False
    cell_df["source_mode"] = None

    for label_name, rec in yml["labels"].items():
        if rec.get("mode") != "cluster":
            continue
        res = str(rec["source_resolution"])
        col = f"leiden_{res}"
        if col not in cell_df.columns:
            raise RuntimeError(f"schema_violation: leiden column {col} not in leiden_csv for label {label_name}")
        cluster_ids = {normalize_cluster_id(c) for c in rec["source_cluster_ids"]}
        mask = cell_df[col].astype(str).isin(cluster_ids)
        existing = cell_df.loc[mask, "label"].notna()
        if existing.any():
            conflict = cell_df.loc[mask & cell_df["label"].notna(), "label"].unique()
            raise RuntimeError(
                f"schema_violation: cluster-mode label {label_name} overlaps with previously assigned label(s) {list(conflict)}"
            )
        cell_df.loc[mask, "label"] = label_name
        cell_df.loc[mask, "lineage"] = rec.get("lineage", "unknown")
        cell_df.loc[mask, "is_artifact"] = bool(rec.get("is_artifact", False))
        cell_df.loc[mask, "source_mode"] = "cluster"
    return cell_df


def apply_fine_cluster_overrides(cell_df: pd.DataFrame, yml: dict) -> pd.DataFrame:
    """Step 2: fine_cluster labels override parent_label on their cluster cells."""
    for label_name, rec in yml["labels"].items():
        if rec.get("mode") != "fine_cluster":
            continue
        res = str(rec["source_resolution"])
        col = f"leiden_{res}"
        cluster_ids = {normalize_cluster_id(c) for c in rec["source_cluster_ids"]}
        parent = rec["parent_label"]
        mask = (cell_df[col].astype(str).isin(cluster_ids)) & (cell_df["label"] == parent)
        n_override = int(mask.sum())
        if n_override == 0:
            log.warning(
                "fine_cluster %s: 0 cells overridden — parent %r may not cover %s at res %s",
                label_name, parent, cluster_ids, res,
            )
        cell_df.loc[mask, "label"] = label_name
        cell_df.loc[mask, "lineage"] = rec.get("lineage", "unknown")
        cell_df.loc[mask, "is_artifact"] = bool(rec.get("is_artifact", False))
        cell_df.loc[mask, "source_mode"] = "fine_cluster"
        log.info("fine_cluster %s: overrode %d parent cells", label_name, n_override)
    return cell_df


def apply_expression_overrides(
    cell_df: pd.DataFrame,
    yml: dict,
    X: sparse.csr_matrix | None,
    ensg_to_idx: dict,
    symbol_to_ensg: dict,
) -> pd.DataFrame:
    """Step 3: expression-mode labels override parent_label where signature matches."""
    for label_name, rec in yml["labels"].items():
        if rec.get("mode") != "expression":
            continue
        if X is None:
            raise RuntimeError(f"expression-mode label {label_name} requires --counts-npz")
        sig = rec["expression_signature"]
        threshold = float(sig.get("threshold", 0.5))
        pos_genes = sig["required_positive"]
        parent = rec["parent_label"]

        gene_indices = []
        for g in pos_genes:
            ensg = symbol_to_ensg.get(g, g if g in ensg_to_idx else None)
            if ensg is None:
                raise RuntimeError(f"required_gene_missing: expression signature for {label_name} references {g!r} not in gene_data")
            gene_indices.append(ensg_to_idx[ensg])
        gene_indices = np.array(gene_indices, dtype=int)

        parent_mask = cell_df["label"] == parent
        if parent_mask.sum() == 0:
            log.warning("expression %s: parent label %s has 0 cells — skipping", label_name, parent)
            continue
        # log1p-transformed mean expression over required_positive genes per cell
        X_sub = X[:, gene_indices]
        X_log = X_sub.copy()
        X_log.data = np.log1p(X_log.data)
        per_cell_mean = np.asarray(X_log.mean(axis=1)).flatten()
        express_mask = parent_mask & (per_cell_mean >= threshold)
        n_override = int(express_mask.sum())
        cell_df.loc[express_mask, "label"] = label_name
        cell_df.loc[express_mask, "lineage"] = rec.get("lineage", "unknown")
        cell_df.loc[express_mask, "is_artifact"] = bool(rec.get("is_artifact", False))
        cell_df.loc[express_mask, "source_mode"] = "expression"
        log.info("expression %s: overrode %d parent cells (threshold=%.2f over %d genes)", label_name, n_override, threshold, len(gene_indices))
    return cell_df


def main() -> None:
    args = parse_args()
    yml = yaml.safe_load(args.v2_yaml.read_text())
    validate_schema(yml)

    log.info("Loading gene_data: %s", args.gene_data)
    gene_data = pd.read_csv(args.gene_data, index_col=0)
    gene_symbols = set(gene_data["gene_symbol"].dropna().astype(str))
    ensg_to_idx = {ensg: i for i, ensg in enumerate(gene_data.index)}
    symbol_to_ensg = {sym: ensg for ensg, sym in gene_data["gene_symbol"].items() if isinstance(sym, str)}
    validate_required_genes(yml, gene_symbols)

    log.info("Loading leiden: %s", args.leiden_csv)
    leiden = pd.read_csv(args.leiden_csv)

    resolutions_used = [str(r) for r in yml.get("resolutions_used", [])]
    log.info("resolutions_used: %s; labels: %d", resolutions_used, len(yml["labels"]))

    X = None
    needs_counts = any(rec.get("mode") == "expression" for rec in yml["labels"].values())
    if needs_counts:
        log.info("Loading counts (expression-mode labels present): %s", args.counts_npz)
        loaded = np.load(args.counts_npz, allow_pickle=True)
        X = sparse.csr_matrix(
            (loaded["data"], loaded["indices"], loaded["indptr"]),
            shape=tuple(loaded["shape"]),
        )

    cell_df = assign_cluster_mode(leiden, yml, resolutions_used)
    cell_df = apply_fine_cluster_overrides(cell_df, yml)
    cell_df = apply_expression_overrides(cell_df, yml, X, ensg_to_idx, symbol_to_ensg)

    unassigned = cell_df["label"].isna().sum()
    if unassigned > 0:
        # Detail which clusters have unassigned cells
        for res in resolutions_used:
            col = f"leiden_{res}"
            if col not in cell_df.columns:
                continue
            unassigned_clusters = cell_df.loc[cell_df["label"].isna(), col].value_counts().head(20)
            log.error("unassigned_cells at %s: %d total. Top clusters: %s", col, unassigned, unassigned_clusters.to_dict())
        raise RuntimeError(f"unassigned_cells: {unassigned} cells have no label (hard failure per §7.1)")

    out = cell_df[["cell_id", "label", "lineage", "is_artifact", "source_mode"]].copy()
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    log.info("Wrote %d rows to %s", len(out), args.out_csv)
    log.info("Label summary:")
    for (label, is_art), n in out.groupby(["label", "is_artifact"], observed=True).size().items():
        flag = " [ARTIFACT]" if is_art else ""
        log.info("  %-35s %7d%s", label, n, flag)


if __name__ == "__main__":
    main()
