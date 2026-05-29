#!/usr/bin/env python3
"""Compare two h5ad files for structural equivalence.

Used to validate regenerated h5ads against Spatial_HBCA originals.
Reports PASS/FAIL per check with details on mismatches.

Usage:
    python3 diff_h5ads.py \
        --new path/to/regenerated.h5ad \
        --reference path/to/original.h5ad \
        [--output report.yaml]

Part of plan: Activation/regen_pilot
"""

import argparse
import sys
from collections import OrderedDict
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import yaml


def compare_h5ads(new_path: str, ref_path: str) -> dict:
    """Compare two h5ad files. Returns structured report."""
    report = OrderedDict()
    report["timestamp"] = datetime.now().isoformat()
    report["new_file"] = new_path
    report["reference_file"] = ref_path
    report["checks"] = OrderedDict()
    all_pass = True

    print(f"Loading new: {new_path}")
    new = ad.read_h5ad(new_path)
    print(f"Loading ref: {ref_path}")
    ref = ad.read_h5ad(ref_path)

    # --- 1. Shape ---
    check = OrderedDict()
    check["new"] = {"cells": new.n_obs, "genes": new.n_vars}
    check["ref"] = {"cells": ref.n_obs, "genes": ref.n_vars}
    check["pass"] = new.shape == ref.shape
    if not check["pass"]:
        check["detail"] = f"Shape mismatch: {new.shape} vs {ref.shape}"
        all_pass = False
    report["checks"]["shape"] = check
    print(f"  Shape: {'PASS' if check['pass'] else 'FAIL'} ({new.shape} vs {ref.shape})")

    # --- 2. obs columns ---
    check = OrderedDict()
    new_cols = set(new.obs.columns)
    ref_cols = set(ref.obs.columns)
    check["new_only"] = sorted(new_cols - ref_cols)
    check["ref_only"] = sorted(ref_cols - new_cols)
    check["shared"] = len(new_cols & ref_cols)
    check["pass"] = len(check["new_only"]) == 0 and len(check["ref_only"]) == 0
    if not check["pass"]:
        all_pass = False
    report["checks"]["obs_columns"] = check
    print(f"  obs columns: {'PASS' if check['pass'] else 'FAIL'} "
          f"(shared={check['shared']}, new_only={len(check['new_only'])}, ref_only={len(check['ref_only'])})")

    # --- 3. obs dtypes (shared columns) ---
    shared_cols = sorted(new_cols & ref_cols)
    dtype_mismatches = []
    for col in shared_cols:
        nd = str(new.obs[col].dtype)
        rd = str(ref.obs[col].dtype)
        if nd != rd:
            dtype_mismatches.append({"column": col, "new": nd, "ref": rd})
    check = OrderedDict()
    check["mismatches"] = dtype_mismatches[:20]
    check["total_mismatches"] = len(dtype_mismatches)
    check["pass"] = len(dtype_mismatches) == 0
    if not check["pass"]:
        all_pass = False
    report["checks"]["obs_dtypes"] = check
    print(f"  obs dtypes: {'PASS' if check['pass'] else 'FAIL'} ({len(dtype_mismatches)} mismatches)")

    # --- 4. var index (gene IDs) ---
    check = OrderedDict()
    new_genes = set(new.var_names)
    ref_genes = set(ref.var_names)
    check["new_only"] = len(new_genes - ref_genes)
    check["ref_only"] = len(ref_genes - new_genes)
    check["shared"] = len(new_genes & ref_genes)
    check["order_match"] = list(new.var_names) == list(ref.var_names)
    check["pass"] = check["new_only"] == 0 and check["ref_only"] == 0
    if not check["pass"]:
        all_pass = False
        check["new_only_examples"] = sorted(new_genes - ref_genes)[:10]
        check["ref_only_examples"] = sorted(ref_genes - new_genes)[:10]
    report["checks"]["var_index"] = check
    print(f"  var index: {'PASS' if check['pass'] else 'FAIL'} "
          f"(shared={check['shared']}, order_match={check['order_match']})")

    # --- 5. var columns ---
    check = OrderedDict()
    new_vcols = set(new.var.columns)
    ref_vcols = set(ref.var.columns)
    check["new_only"] = sorted(new_vcols - ref_vcols)
    check["ref_only"] = sorted(ref_vcols - new_vcols)
    check["pass"] = len(check["new_only"]) == 0 and len(check["ref_only"]) == 0
    if not check["pass"]:
        all_pass = False
    report["checks"]["var_columns"] = check
    print(f"  var columns: {'PASS' if check['pass'] else 'FAIL'}")

    # --- 6. obsm keys ---
    check = OrderedDict()
    new_obsm = set(new.obsm.keys())
    ref_obsm = set(ref.obsm.keys())
    check["new_only"] = sorted(new_obsm - ref_obsm)
    check["ref_only"] = sorted(ref_obsm - new_obsm)
    check["shared"] = sorted(new_obsm & ref_obsm)
    # new_only = additions (INFO), ref_only = data loss (FAIL)
    check["pass"] = len(check["ref_only"]) == 0
    if not check["pass"]:
        all_pass = False
    # Check shapes of shared embeddings
    shape_mismatches = []
    for key in check["shared"]:
        ns = new.obsm[key].shape
        rs = ref.obsm[key].shape
        if ns != rs:
            shape_mismatches.append({"key": key, "new": list(ns), "ref": list(rs)})
    check["shape_mismatches"] = shape_mismatches
    report["checks"]["obsm"] = check
    status = "PASS" if check["pass"] else "FAIL"
    detail = f"shared={check['shared']}"
    if check["new_only"]:
        detail += f", new_only={check['new_only']}"
    if check["ref_only"]:
        detail += f", ref_only={check['ref_only']}"
    print(f"  obsm: {status} ({detail})")

    # --- 7. uns keys ---
    check = OrderedDict()
    new_uns = set(new.uns.keys())
    ref_uns = set(ref.uns.keys())
    check["new_only"] = sorted(new_uns - ref_uns)
    check["ref_only"] = sorted(ref_uns - new_uns)
    check["shared"] = sorted(new_uns & ref_uns)
    # uns differences are expected (schema_version, title may change)
    check["pass"] = True
    report["checks"]["uns"] = check
    print(f"  uns: INFO (new_only={check['new_only']}, ref_only={check['ref_only']})")

    # --- 8. X format and stats ---
    check = OrderedDict()
    from scipy import sparse
    check["new_format"] = type(new.X).__name__
    check["ref_format"] = type(ref.X).__name__
    check["new_dtype"] = str(new.X.dtype)
    check["ref_dtype"] = str(ref.X.dtype)
    if sparse.issparse(new.X) and sparse.issparse(ref.X):
        check["new_nnz"] = int(new.X.nnz)
        check["ref_nnz"] = int(ref.X.nnz)
        check["nnz_match"] = check["new_nnz"] == check["ref_nnz"]
    check["format_match"] = check["new_format"] == check["ref_format"]
    check["dtype_match"] = check["new_dtype"] == check["ref_dtype"]
    check["pass"] = check.get("nnz_match", True) and check["format_match"]
    if not check["pass"]:
        all_pass = False
    report["checks"]["X_matrix"] = check
    print(f"  X matrix: {'PASS' if check['pass'] else 'FAIL'} "
          f"(format={check['new_format']}, nnz={'match' if check.get('nnz_match', 'N/A') else 'MISMATCH'})")

    # --- 9. raw layer ---
    check = OrderedDict()
    check["new_has_raw"] = new.raw is not None
    check["ref_has_raw"] = ref.raw is not None
    check["pass"] = check["new_has_raw"] == check["ref_has_raw"]
    if not check["pass"]:
        all_pass = False
    report["checks"]["raw"] = check
    print(f"  raw: {'PASS' if check['pass'] else 'FAIL'} "
          f"(new={check['new_has_raw']}, ref={check['ref_has_raw']})")

    # --- 10. Spot-check obs values (first 5 cells, shared columns) ---
    check = OrderedDict()
    sample_cols = [c for c in ["donor_id", "cell_type", "tissue_ontology_term_id",
                                "assay_ontology_term_id", "sex_ontology_term_id"]
                   if c in shared_cols]
    spot_mismatches = []
    for col in sample_cols:
        n_vals = new.obs[col].head(5).tolist()
        r_vals = ref.obs[col].head(5).tolist()
        if n_vals != r_vals:
            spot_mismatches.append({
                "column": col,
                "new_sample": [str(v) for v in n_vals],
                "ref_sample": [str(v) for v in r_vals],
            })
    check["columns_checked"] = sample_cols
    check["mismatches"] = spot_mismatches
    check["pass"] = len(spot_mismatches) == 0
    if not check["pass"]:
        all_pass = False
    report["checks"]["obs_spot_check"] = check
    print(f"  obs spot check: {'PASS' if check['pass'] else 'FAIL'} "
          f"({len(sample_cols)} columns, {len(spot_mismatches)} mismatches)")

    report["overall_pass"] = all_pass
    print(f"\nOVERALL: {'PASS' if all_pass else 'FAIL'}")
    return report


def main():
    parser = argparse.ArgumentParser(description="Compare two h5ad files")
    parser.add_argument("--new", required=True, help="Path to regenerated h5ad")
    parser.add_argument("--reference", required=True, help="Path to original h5ad")
    parser.add_argument("--output", help="Output YAML report path (default: stdout)")
    args = parser.parse_args()

    for path in [args.new, args.reference]:
        if not Path(path).exists():
            print(f"ERROR: File not found: {path}", file=sys.stderr)
            sys.exit(1)

    report = compare_h5ads(args.new, args.reference)

    # Custom YAML representer for OrderedDict
    def represent_ordereddict(dumper, data):
        return dumper.represent_mapping("tag:yaml.org,2002:map", data.items())
    yaml.add_representer(OrderedDict, represent_ordereddict)

    yaml_str = yaml.dump(report, default_flow_style=False, allow_unicode=True, width=120)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            f.write(yaml_str)
        print(f"\nReport written to: {args.output}")
    else:
        print("\n--- YAML Report ---")
        print(yaml_str)

    sys.exit(0 if report["overall_pass"] else 1)


if __name__ == "__main__":
    main()
