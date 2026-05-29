#!/usr/bin/env python3
"""
Regenerate reed_cells.csv using pure Python (no reticulate).

Reed's source is h5ad, which map_study_generic.R loads via reticulate.
Reticulate fails on HPC singularity containers. This script replaces
the R mapping for reed only.

Mapping type: identity (cell IDs match 1:1 between reed.h5ad and iHBCA inventory).

Output format matches map_study_generic.R exactly:
  ihbca_cell_id, reed_component_cell_id, reed_mapping_method, reed_{obs_col}...
"""

import sys
import argparse
from pathlib import Path

import anndata
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Reed cell mapping (Python-only)")
    parser.add_argument(
        "--h5ad",
        default="${SOURCE_COMPONENT_STUDIES}/reed.h5ad",
        help="Path to reed.h5ad",
    )
    parser.add_argument(
        "--inventory",
        default=None,
        help="Path to ihbca_cell_inventory.csv (auto-detected from repo root)",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Output directory (default: harmonization/studies/reed/outputs/)",
    )
    parser.add_argument(
        "--repo-root",
        default="${PROJECT_IHBCAV1_UPLOAD}",
        help="Repository root",
    )
    args = parser.parse_args()

    repo_root = Path(args.repo_root)

    if args.inventory is None:
        inventory_path = (
            repo_root
            / "external_studies/harmonization/outputs/ihbca_reference/ihbca_cell_inventory.csv"
        )
    else:
        inventory_path = Path(args.inventory)

    if args.output_dir is None:
        output_dir = (
            repo_root / "external_studies/harmonization/studies/reed/outputs"
        )
    else:
        output_dir = Path(args.output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Load iHBCA inventory — filter to reed cells
    print("Loading iHBCA inventory...")
    inv = pd.read_csv(inventory_path, usecols=["cell_id", "dataset"])
    reed_inv = inv[inv["dataset"] == "reed"]
    ihbca_cell_ids = set(reed_inv["cell_id"])
    print(f"  Total iHBCA cells: {len(inv):,}")
    print(f"  Reed cells in inventory: {len(ihbca_cell_ids):,}")

    # Load reed.h5ad
    print(f"\nLoading reed.h5ad: {args.h5ad}")
    adata = anndata.read_h5ad(args.h5ad, backed="r")
    obs = adata.obs.copy()
    component_ids = obs.index.tolist()
    print(f"  Cells: {len(component_ids):,}")
    print(f"  Obs columns: {len(obs.columns)}")

    # Identity mapping
    print("\nRunning identity mapping...")
    mapped_ids = [cid for cid in component_ids if cid in ihbca_cell_ids]
    unmapped_ids = [cid for cid in component_ids if cid not in ihbca_cell_ids]
    print(f"  Mapped: {len(mapped_ids):,}")
    print(f"  Unmapped: {len(unmapped_ids):,}")
    print(f"  Match rate: {100 * len(mapped_ids) / len(component_ids):.1f}%")

    if len(unmapped_ids) > 0:
        print(f"  WARNING: {len(unmapped_ids)} unmapped cells")

    # Build output DataFrame
    # Filter obs to mapped cells
    mapped_obs = obs.loc[mapped_ids].copy()
    mapped_obs.index.name = None

    # Add mapping columns
    result = pd.DataFrame(
        {
            "ihbca_cell_id": mapped_ids,
            "component_cell_id": mapped_ids,
            "mapping_method": "identity",
        }
    )

    # Merge with metadata (obs columns)
    # The R script merges on component_cell_id == _cell_id (obs index)
    mapped_obs["_cell_id"] = mapped_obs.index
    result = result.merge(mapped_obs, left_on="component_cell_id", right_on="_cell_id")
    result.drop(columns=["_cell_id"], inplace=True)

    # Prefix all columns except ihbca_cell_id with reed_
    rename_map = {}
    for col in result.columns:
        if col == "ihbca_cell_id":
            continue
        if col.startswith("reed_"):
            continue
        rename_map[col] = f"reed_{col}"
    result.rename(columns=rename_map, inplace=True)

    # Reorder: ihbca_cell_id first
    cols = ["ihbca_cell_id"] + [c for c in result.columns if c != "ihbca_cell_id"]
    result = result[cols]

    # Write output
    cells_path = output_dir / "reed_cells.csv"
    print(f"\nWriting: {cells_path}")
    result.to_csv(cells_path, index=False)
    print(f"  Rows: {len(result):,}")
    print(f"  Columns: {len(result.columns)}")

    # Verify
    print(f"\n{'=' * 50}")
    print(f"VERIFICATION")
    print(f"  Expected: 803,283 cells")
    print(f"  Got:      {len(result):,} cells")
    print(f"  Columns:  {len(result.columns)}")
    print(f"  Match:    {'PASS' if len(result) == 803283 else 'CHECK'}")
    print(f"{'=' * 50}")


if __name__ == "__main__":
    main()
