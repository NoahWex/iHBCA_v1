#!/usr/bin/env python3
"""Compare current h5ads against baseline manifest.

Reports structural changes (cell counts, gene counts, obs columns, obsm keys,
var columns, uns keys) between the baseline manifest and current h5ad files.
SHA256 differences are reported but expected when code changes affect outputs.

Usage:
    python3 diff_against_baseline.py \
        --manifest /path/to/baseline_manifest.yaml \
        --datasets all-breast-cells breast-epithelial-lineage  # optional filter
"""

import argparse
import os
import time

import yaml


def diff_dataset(label, baseline, current_path):
    """Compare one dataset against its baseline entry. Returns list of diffs."""
    import anndata

    diffs = []

    if not os.path.exists(current_path):
        return [f"MISSING: {current_path}"]

    t0 = time.time()
    adata = anndata.read_h5ad(current_path, backed="r")

    # Cell count
    if adata.n_obs != baseline["n_cells"]:
        diffs.append(f"n_cells: {baseline['n_cells']} -> {adata.n_obs}")

    # Gene count
    if adata.n_vars != baseline["n_genes"]:
        diffs.append(f"n_genes: {baseline['n_genes']} -> {adata.n_vars}")

    # Obs columns
    current_obs = sorted(adata.obs.columns.tolist())
    baseline_obs = baseline["obs_columns"]
    added_obs = set(current_obs) - set(baseline_obs)
    removed_obs = set(baseline_obs) - set(current_obs)
    if added_obs:
        diffs.append(f"obs columns added: {sorted(added_obs)}")
    if removed_obs:
        diffs.append(f"obs columns removed: {sorted(removed_obs)}")

    # Obs dtypes (only for columns present in both)
    if "obs_dtypes" in baseline:
        common_cols = set(current_obs) & set(baseline_obs)
        dtype_changes = []
        for col in sorted(common_cols):
            cur_dtype = str(adata.obs[col].dtype)
            base_dtype = baseline["obs_dtypes"].get(col, "")
            if cur_dtype != base_dtype:
                dtype_changes.append(f"  {col}: {base_dtype} -> {cur_dtype}")
        if dtype_changes:
            diffs.append("obs dtype changes:\n" + "\n".join(dtype_changes))

    # Var columns
    current_var = sorted(adata.var.columns.tolist())
    baseline_var = baseline.get("var_columns", [])
    added_var = set(current_var) - set(baseline_var)
    removed_var = set(baseline_var) - set(current_var)
    if added_var:
        diffs.append(f"var columns added: {sorted(added_var)}")
    if removed_var:
        diffs.append(f"var columns removed: {sorted(removed_var)}")

    # Obsm keys
    current_obsm = sorted(list(adata.obsm.keys()))
    baseline_obsm = baseline.get("obsm_keys", [])
    added_obsm = set(current_obsm) - set(baseline_obsm)
    removed_obsm = set(baseline_obsm) - set(current_obsm)
    if added_obsm:
        diffs.append(f"obsm keys added: {sorted(added_obsm)}")
    if removed_obsm:
        diffs.append(f"obsm keys removed: {sorted(removed_obsm)}")

    # Uns keys
    current_uns = sorted(list(adata.uns.keys()))
    baseline_uns = baseline.get("uns_keys", [])
    added_uns = set(current_uns) - set(baseline_uns)
    removed_uns = set(baseline_uns) - set(current_uns)
    if added_uns:
        diffs.append(f"uns keys added: {sorted(added_uns)}")
    if removed_uns:
        diffs.append(f"uns keys removed: {sorted(removed_uns)}")

    elapsed = time.time() - t0
    adata.file.close()

    return diffs


def main():
    parser = argparse.ArgumentParser(description="Diff h5ads against baseline manifest")
    parser.add_argument("--manifest", required=True, help="Path to baseline_manifest.yaml")
    parser.add_argument("--repo-root", required=True, help="Path to repo root")
    parser.add_argument("--datasets", nargs="*", help="Specific datasets to check (default: all)")
    args = parser.parse_args()

    with open(args.manifest) as f:
        manifest = yaml.safe_load(f)

    repo = args.repo_root
    source_dir = os.path.join(repo, "publication", "outputs", "source_datasets")
    integrated_dir = os.path.join(repo, "publication", "outputs", "integrated_objects")

    datasets = manifest["datasets"]
    if args.datasets:
        datasets = {k: v for k, v in datasets.items() if k in args.datasets}

    print(f"Comparing {len(datasets)} datasets against baseline from {manifest['generated_at']}")
    print("=" * 70)

    n_pass = 0
    n_diff = 0
    n_missing = 0

    for label, baseline in datasets.items():
        # Determine path
        if label.startswith("breast-") or label == "all-breast-cells":
            current_path = os.path.join(integrated_dir, f"{label}.h5ad")
        else:
            current_path = os.path.join(source_dir, f"{label}.h5ad")

        print(f"\n{label}:")
        diffs = diff_dataset(label, baseline, current_path)

        if not diffs:
            print(f"  PASS (structure identical)")
            n_pass += 1
        elif diffs == [f"MISSING: {current_path}"]:
            print(f"  MISSING: {current_path}")
            n_missing += 1
        else:
            print(f"  DIFF ({len(diffs)} change(s)):")
            for d in diffs:
                print(f"    - {d}")
            n_diff += 1

    print("\n" + "=" * 70)
    print(f"SUMMARY: {n_pass} pass, {n_diff} diff, {n_missing} missing "
          f"(of {len(datasets)} total)")

    if n_diff == 0 and n_missing == 0:
        print("GATE: PASS — all structures match baseline")
    else:
        print("GATE: REVIEW — structural changes detected")


if __name__ == "__main__":
    main()
