#!/usr/bin/env python3
"""
Build Unified Cell Metadata (Raw UNION) for IntegrationAssessment.

This script merges:
1. iHBCA reference (67 columns from author_share) - prefixed with ihbca_
2. Phase A per-study mappings (9 studies) - already prefixed with {study}_

The result is a raw union with NO harmonization - just all columns combined
with clear source attribution via prefixes.

Outputs:
- unified_cell_metadata.csv: Full cell inventory (~2.1M cells)
- unified_cell_metadata.parquet: Same in columnar format

Memory estimate: ~8-10GB peak. Designed for HPC execution.

Why do programmers prefer dark mode? Because light attracts bugs.
"""

import argparse
import gc
import os
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


# =============================================================================
# CONFIGURATION
# =============================================================================

# iHBCA source file (67 columns, ~2.1M cells)
IHBCA_SOURCE = "${SOURCE_AUTHOR_SHARE}/ihbca_level1.5_annotations.csv"

# Phase A per-study mapping files (relative to BASE_DIR)
STUDIES = [
    'gray',
    'kumar',
    'murrow',
    'nee',
    'twigger',
    'pal_norm_epi',
    'pal_norm_total',
    'pal_norm_b1',
    'reed'
]

# Base directory (will be set from args or default)
DEFAULT_BASE_DIR = "${SOURCE_IHBCAV1_HARMONIZATION%/}"


# =============================================================================
# MAIN LOGIC
# =============================================================================

def load_ihbca_reference(ihbca_path: Path) -> pd.DataFrame:
    """Load iHBCA reference and prefix all columns with ihbca_."""
    print(f"\n=== Loading iHBCA Reference ===")
    print(f"Source: {ihbca_path}")

    if not ihbca_path.exists():
        raise FileNotFoundError(f"iHBCA source not found: {ihbca_path}")

    # Read the full CSV
    print("Reading CSV (this may take a moment)...")
    ihbca = pd.read_csv(ihbca_path, low_memory=False)

    print(f"Loaded: {len(ihbca):,} cells x {len(ihbca.columns)} columns")

    # Rename cellID to ihbca_cell_id (this is the join key)
    if 'cellID' not in ihbca.columns:
        raise ValueError("Expected 'cellID' column in iHBCA source")

    # Prefix ALL columns with ihbca_ (except cellID which becomes ihbca_cell_id)
    new_columns = {}
    for col in ihbca.columns:
        if col == 'cellID':
            new_columns[col] = 'ihbca_cell_id'
        else:
            new_columns[col] = f'ihbca_{col}'

    ihbca = ihbca.rename(columns=new_columns)

    print(f"Columns prefixed with 'ihbca_'")
    print(f"Sample columns: {list(ihbca.columns[:5])}")

    return ihbca


def load_study_cells(base_dir: Path, study: str) -> pd.DataFrame:
    """Load Phase A per-study cell mapping."""
    cells_path = base_dir / 'studies' / study / 'outputs' / f'{study}_cells.csv'

    if not cells_path.exists():
        print(f"  WARNING: {study}_cells.csv not found at {cells_path}")
        return None

    df = pd.read_csv(cells_path, low_memory=False)
    print(f"  {study}: {len(df):,} cells x {len(df.columns)} columns")

    # Verify ihbca_cell_id is present
    if 'ihbca_cell_id' not in df.columns:
        print(f"  WARNING: {study} missing ihbca_cell_id column!")
        return None

    return df


def merge_study(unified: pd.DataFrame, study_df: pd.DataFrame, study: str) -> pd.DataFrame:
    """Left-join a study's cells to the unified dataframe."""
    # Get columns from study_df that aren't already in unified (except join key)
    existing_cols = set(unified.columns)
    new_cols = [c for c in study_df.columns if c not in existing_cols or c == 'ihbca_cell_id']

    if len(new_cols) <= 1:  # Only ihbca_cell_id
        print(f"  {study}: No new columns to add")
        return unified

    study_subset = study_df[new_cols]

    # Left join on ihbca_cell_id
    unified = unified.merge(study_subset, on='ihbca_cell_id', how='left')

    added_cols = len(new_cols) - 1  # Exclude join key
    print(f"  {study}: Added {added_cols} columns -> total {len(unified.columns)}")

    return unified


def verify_output(df: pd.DataFrame, studies: list):
    """Verify the unified output has expected structure."""
    print("\n=== Verification ===")

    # Check cell count
    print(f"Total cells: {len(df):,}")

    # Check column prefixes
    ihbca_cols = [c for c in df.columns if c.startswith('ihbca_')]
    print(f"iHBCA columns: {len(ihbca_cols)}")

    for study in studies:
        study_cols = [c for c in df.columns if c.startswith(f'{study}_')]
        if study_cols:
            print(f"{study} columns: {len(study_cols)}")
        else:
            print(f"{study} columns: 0 (WARNING - no columns found!)")

    # Check for unprefixed columns
    prefixes = ['ihbca_'] + [f'{s}_' for s in studies]
    unprefixed = [c for c in df.columns if not any(c.startswith(p) for p in prefixes)]
    if unprefixed:
        print(f"Unprefixed columns: {unprefixed}")

    # Check coverage per study (non-null in first study column)
    print("\nPer-study coverage (cells with non-null data):")
    for study in studies:
        study_cols = [c for c in df.columns if c.startswith(f'{study}_')]
        if study_cols:
            # Check first column for coverage
            sample_col = study_cols[0]
            coverage = df[sample_col].notna().sum()
            pct = 100 * coverage / len(df)
            print(f"  {study}: {coverage:,} cells ({pct:.1f}%)")


def generate_column_stats(df: pd.DataFrame, output_path: Path):
    """Generate column inventory with statistics."""
    print("\n=== Generating Column Inventory ===")

    column_info = {}
    for i, col in enumerate(df.columns):
        info = {
            'position': i + 1,
            'dtype': str(df[col].dtype),
            'n_null': int(df[col].isna().sum()),
            'n_non_null': int(df[col].notna().sum()),
            'coverage_pct': round(100 * df[col].notna().sum() / len(df), 1)
        }

        # Identify source
        if col == 'ihbca_cell_id':
            info['source'] = 'primary_key'
        elif col.startswith('ihbca_'):
            info['source'] = 'ihbca_reference'
            info['original_name'] = col[6:]  # Strip ihbca_ prefix
        else:
            for study in STUDIES:
                if col.startswith(f'{study}_'):
                    info['source'] = f'study:{study}'
                    info['original_name'] = col[len(study)+1:]
                    break

        # For categorical columns, get unique count and sample values
        if df[col].dtype == 'object' or df[col].nunique() < 100:
            n_unique = df[col].nunique()
            info['n_unique'] = int(n_unique)

            # Get top values if reasonable number
            if n_unique <= 20:
                value_counts = df[col].value_counts(dropna=False).head(20)
                info['value_counts'] = {
                    str(k) if pd.notna(k) else 'NULL': int(v)
                    for k, v in value_counts.items()
                }

        column_info[col] = info

    # Save as YAML
    inventory = {
        'generated_at': datetime.now().isoformat(),
        'total_cells': len(df),
        'total_columns': len(df.columns),
        'columns': column_info
    }

    with open(output_path, 'w') as f:
        yaml.dump(inventory, f, default_flow_style=False, sort_keys=False)

    print(f"Column inventory saved to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Build unified cell metadata (raw union)"
    )
    parser.add_argument(
        '--base-dir',
        type=str,
        default=DEFAULT_BASE_DIR,
        help="Base directory for harmonization"
    )
    parser.add_argument(
        '--ihbca-source',
        type=str,
        default=IHBCA_SOURCE,
        help="Path to iHBCA level1.5_annotations.csv"
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help="Output directory (default: {base_dir}/outputs/unified_metadata)"
    )
    parser.add_argument(
        '--test',
        action='store_true',
        help="Test mode: sample 10000 rows from iHBCA"
    )
    parser.add_argument(
        '--skip-parquet',
        action='store_true',
        help="Skip parquet output (useful for testing)"
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    ihbca_path = Path(args.ihbca_source)
    output_dir = Path(args.output_dir) if args.output_dir else base_dir / 'outputs' / 'unified_metadata'

    print("=" * 70)
    print("BUILD UNIFIED CELL METADATA (RAW UNION)")
    print("=" * 70)
    print(f"Base directory: {base_dir}")
    print(f"iHBCA source: {ihbca_path}")
    print(f"Output directory: {output_dir}")
    print(f"Test mode: {args.test}")
    print()

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Load iHBCA reference
    ihbca = load_ihbca_reference(ihbca_path)

    if args.test:
        print(f"\nTEST MODE: Sampling 10000 cells from iHBCA")
        ihbca = ihbca.sample(n=min(10000, len(ihbca)), random_state=42)
        print(f"Sample size: {len(ihbca):,} cells")

    # Step 2: Load and merge each study
    print(f"\n=== Loading and Merging Studies ===")
    unified = ihbca.copy()
    del ihbca  # Free memory
    gc.collect()

    for study in STUDIES:
        study_df = load_study_cells(base_dir, study)
        if study_df is not None:
            unified = merge_study(unified, study_df, study)
            del study_df
            gc.collect()

    # Step 3: Verify output
    verify_output(unified, STUDIES)

    # Step 4: Generate column inventory
    inventory_path = output_dir / 'column_inventory.yaml'
    generate_column_stats(unified, inventory_path)

    # Step 5: Save outputs
    print(f"\n=== Saving Outputs ===")

    csv_path = output_dir / 'unified_cell_metadata.csv'
    print(f"Saving CSV to: {csv_path}")
    unified.to_csv(csv_path, index=False)
    print(f"CSV saved: {csv_path.stat().st_size / (1024*1024):.1f} MB")

    if not args.skip_parquet:
        parquet_path = output_dir / 'unified_cell_metadata.parquet'
        print(f"Saving Parquet to: {parquet_path}")
        unified.to_parquet(parquet_path, index=False)
        print(f"Parquet saved: {parquet_path.stat().st_size / (1024*1024):.1f} MB")

    # Final summary
    print("\n" + "=" * 70)
    print("BUILD COMPLETE")
    print("=" * 70)
    print(f"Total cells: {len(unified):,}")
    print(f"Total columns: {len(unified.columns)}")
    print(f"\nOutputs:")
    print(f"  CSV: {csv_path}")
    if not args.skip_parquet:
        print(f"  Parquet: {parquet_path}")
    print(f"  Column inventory: {inventory_path}")


if __name__ == "__main__":
    main()
