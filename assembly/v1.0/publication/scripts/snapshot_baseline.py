#!/usr/bin/env python3
"""Snapshot baseline manifest for regression verification.

Records SHA256, cell counts, gene counts, obs columns/dtypes, and obsm keys
for all 11 h5ads (7 source + 4 integrated). Output is a YAML file that
diff_h5ads.py (or manual inspection) can compare against after refactoring.

Usage:
    python3 snapshot_baseline.py --repo-root /path/to/iHBCAv1_upload \
                                 --output /path/to/baseline_manifest.yaml
"""

import argparse
import hashlib
import os
import time

import yaml


def sha256_file(path, chunk_size=1 << 20):
    """Compute SHA256 of a file without loading it entirely into memory."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def snapshot_h5ad(path):
    """Return a dict of properties for one h5ad file."""
    import anndata

    info = {
        "path": path,
        "size_bytes": os.path.getsize(path),
        "sha256": sha256_file(path),
    }

    print(f"  Reading {os.path.basename(path)} ...", flush=True)
    t0 = time.time()
    adata = anndata.read_h5ad(path, backed="r")

    info["n_cells"] = adata.n_obs
    info["n_genes"] = adata.n_vars

    # obs columns and dtypes
    info["obs_columns"] = sorted(adata.obs.columns.tolist())
    info["obs_dtypes"] = {
        col: str(adata.obs[col].dtype) for col in sorted(adata.obs.columns)
    }

    # var columns
    info["var_columns"] = sorted(adata.var.columns.tolist())

    # obsm keys
    info["obsm_keys"] = sorted(list(adata.obsm.keys()))

    # uns keys (top-level only)
    info["uns_keys"] = sorted(list(adata.uns.keys()))

    elapsed = time.time() - t0
    print(
        f"    {adata.n_obs:,} cells x {adata.n_vars:,} genes, "
        f"{len(info['obs_columns'])} obs cols, "
        f"{len(info['obsm_keys'])} obsm keys "
        f"({elapsed:.1f}s)",
        flush=True,
    )

    adata.file.close()
    return info


def main():
    parser = argparse.ArgumentParser(description="Snapshot h5ad baseline manifest")
    parser.add_argument("--repo-root", required=True, help="Path to iHBCAv1_upload root")
    parser.add_argument("--output", required=True, help="Output YAML path")
    args = parser.parse_args()

    repo = args.repo_root
    source_dir = os.path.join(repo, "publication", "outputs", "source_datasets")
    integrated_dir = os.path.join(repo, "publication", "outputs", "integrated_objects")

    source_files = [
        "gray2022.h5ad",
        "kumar2023.h5ad",
        "murrow2022.h5ad",
        "nee2023.h5ad",
        "twigger2022.h5ad",
        "reed2024.h5ad",
        "pal2021.h5ad",
    ]

    integrated_files = [
        "all-breast-cells.h5ad",
        "breast-epithelial-lineage.h5ad",
        "breast-stromal-lineage.h5ad",
        "breast-immune-lineage.h5ad",
    ]

    manifest = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "purpose": "Pre-refactor baseline for pipeline_rebuild gated migration",
        "datasets": {},
    }

    # Source datasets
    print("=== Source datasets ===", flush=True)
    for fname in source_files:
        path = os.path.join(source_dir, fname)
        if not os.path.exists(path):
            print(f"  SKIP (not found): {path}")
            continue
        study = fname.replace(".h5ad", "")
        manifest["datasets"][study] = snapshot_h5ad(path)

    # Integrated objects
    print("\n=== Integrated objects ===", flush=True)
    for fname in integrated_files:
        path = os.path.join(integrated_dir, fname)
        if not os.path.exists(path):
            print(f"  SKIP (not found): {path}")
            continue
        label = fname.replace(".h5ad", "")
        manifest["datasets"][label] = snapshot_h5ad(path)

    # Summary
    found = len(manifest["datasets"])
    expected = len(source_files) + len(integrated_files)
    print(f"\n=== Summary: {found}/{expected} datasets snapshotted ===", flush=True)

    # Write manifest
    with open(args.output, "w") as f:
        yaml.dump(manifest, f, default_flow_style=False, sort_keys=False, width=120)

    print(f"Manifest written to {args.output}")


if __name__ == "__main__":
    main()
