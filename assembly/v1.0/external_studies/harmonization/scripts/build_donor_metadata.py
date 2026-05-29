#!/usr/bin/env python3
"""
Build Unified Donor Metadata (Raw UNION) for IntegrationAssessment.

Phase 3: Uses verified donor_id_mapping.yaml files to join iHBCA donors
with study reference CSVs.

Mapping files verified via cell count comparison (2026-01-22):
- pal: 22/22 donors (underscore→hyphen transform)
- reed: 55/55 donors (identity - source study)
- twigger: 18/18 donors (prefix substitution, cell count verified)
- nee: 22/22 donors (complex mapping, cell count verified)
- gray: direct match (no mapping file needed)
- murrow: direct match (no mapping file needed)

Outputs:
- unified_donor_metadata.csv
- unified_donor_metadata.parquet

Fun fact: This script has been rewritten more times than a PhD thesis abstract.
"""

import argparse
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_BASE_DIR = "${SOURCE_IHBCAV1_HARMONIZATION%/}"
IHBCA_SOURCE = "${SOURCE_AUTHOR_SHARE}/ihbca_level1.5_annotations.csv"

# Reference CSV configuration per study
# Only for studies that need external reference metadata (not reed - source study)
REFERENCE_CONFIG = {
    'gray': {
        'file': 'studies/gray/reference/GRAY_scRNA-seq.csv',
        'id_column': 'scRNA-seq',
        'skip_rows': 0,
        'mapping_file': None,  # Direct match
    },
    'murrow': {
        'file': 'studies/murrow/reference/MURROW_Table_S1.csv',
        'id_column': 'Sample ID',
        'skip_rows': 1,
        'mapping_file': None,  # Direct match
    },
    'pal': {
        'file': 'studies/pal_shared/reference/pal_table_ev1_extracted.csv',
        'id_column': 'specimen_id',
        'skip_rows': 0,
        'mapping_file': 'studies/pal_shared/outputs/donor_id_mapping.yaml',
    },
    'twigger': {
        'file': 'studies/twigger/reference/twigger_suppl_fig4.csv',
        'id_column': 'Sample',
        'skip_rows': 11,
        'mapping_file': 'studies/twigger/outputs/donor_id_mapping.yaml',
    },
    'nee': {
        'file': 'studies/nee/reference/NEE_Table_1.csv',
        'id_column': 'Patient ID',
        'skip_rows': 2,
        'mapping_file': 'studies/nee/outputs/donor_id_mapping.yaml',
        # Multi-section CSV: filter to scRNAseq rows only
        'filter_column': 'Application',
        'filter_value': 'scRNAseq',
    },
    'reed': {
        'file': 'studies/reed/reference/REED_Supplementary_Table_1.csv',
        'id_column': 'BCN_donor_id',
        'skip_rows': 1,
        'mapping_file': 'studies/reed/outputs/donor_id_mapping.yaml',
        'aggregate_to_donor': True,  # Multiple samples per donor (161 rows / 55 donors)
    },
    'kumar': {
        'file': 'studies/kumar/reference/KUMAR_table1_complete.csv',
        'id_column': 'Patient_ID',
        'skip_rows': 7,  # Skip comment header lines
        'mapping_file': None,  # Direct match (P01 = P01)
    },
}


# =============================================================================
# MAPPING LOADERS
# =============================================================================

def load_donor_mapping(mapping_path: Path) -> dict:
    """Load donor_id_mapping.yaml and return as lookup dict."""
    if not mapping_path.exists():
        return None

    with open(mapping_path) as f:
        data = yaml.safe_load(f)

    # Build lookup: ihbca_donor_id -> reference_donor_id
    lookup = {}
    for m in data.get('mappings', []):
        ihbca_id = m['ihbca_donor_id']
        ref_id = m['reference_donor_id']
        lookup[ihbca_id] = ref_id

    return lookup


def get_unique_donors_from_ihbca(ihbca_path: Path) -> pd.DataFrame:
    """Extract unique donors from iHBCA reference with key metadata columns."""
    print(f"\n=== Extracting Donors from iHBCA ===")
    print(f"Source: {ihbca_path}")

    # Read donor-relevant columns
    cols_to_read = [
        'patientID', 'dataset',
        'age', 'parous', 'parity', 'risk_status',
        'tissue_origin', 'FACS_status', 'sample_type',
        'brca_status', 'ethnicity', 'BMI'
    ]

    df = pd.read_csv(ihbca_path, usecols=cols_to_read, low_memory=False)

    # Get unique donors (one row per donor)
    donors = df.drop_duplicates(subset=['patientID', 'dataset'])

    # Rename to ihbca_ prefix
    donors = donors.rename(columns={
        'patientID': 'ihbca_donor_id',
        'dataset': 'ihbca_dataset'
    })

    # Prefix remaining columns
    rename_map = {c: f'ihbca_{c}' for c in donors.columns
                  if c not in ['ihbca_donor_id', 'ihbca_dataset']}
    donors = donors.rename(columns=rename_map)

    print(f"Unique donors: {len(donors)}")
    print("Per dataset:")
    for ds, count in donors['ihbca_dataset'].value_counts().items():
        print(f"  {ds}: {count}")

    return donors


# =============================================================================
# REFERENCE CSV LOADERS
# =============================================================================

def load_reference_csv(base_dir: Path, study: str, config: dict) -> pd.DataFrame:
    """Load and clean a study's reference CSV."""
    if config.get('file') is None:
        return None

    csv_path = base_dir / config['file']
    if not csv_path.exists():
        print(f"  {study}: Reference CSV not found: {csv_path}")
        return None

    skip_rows = config.get('skip_rows', 0)
    df = pd.read_csv(csv_path, skiprows=skip_rows, low_memory=False)
    df.columns = df.columns.str.strip()

    # Apply row filter if specified (for multi-section CSVs like Nee)
    filter_col = config.get('filter_column')
    filter_val = config.get('filter_value')
    if filter_col and filter_val:
        if filter_col in df.columns:
            before_count = len(df)
            df = df[df[filter_col] == filter_val].copy()
            print(f"  {study}: Filtered {filter_col}=={filter_val}: {before_count} -> {len(df)} rows")
        else:
            print(f"  {study}: Filter column '{filter_col}' not found")

    # Identify ID column and clean values
    id_col = config.get('id_column')
    if id_col and id_col in df.columns:
        # Strip whitespace from ID values (e.g., Nee has trailing spaces)
        df[id_col] = df[id_col].astype(str).str.strip()
        df = df.rename(columns={id_col: f'{study}_ref_donor_id'})
    else:
        print(f"  {study}: ID column '{id_col}' not found in {list(df.columns)[:5]}")
        return None

    # Prefix all other columns
    rename_map = {c: f'{study}_{c}' for c in df.columns
                  if c != f'{study}_ref_donor_id'}
    df = df.rename(columns=rename_map)

    # Aggregate to donor level if needed (e.g., Reed has multiple samples per donor)
    if config.get('aggregate_to_donor'):
        ref_id = f'{study}_ref_donor_id'
        before = len(df)
        # Take first non-null value per donor for each column
        df = df.groupby(ref_id, as_index=False).first()
        print(f"  {study}: Aggregated {before} rows -> {len(df)} donors (first non-null per column)")

    print(f"  {study}: Loaded {len(df)} rows x {len(df.columns)} columns")
    return df


# =============================================================================
# MAIN MERGE LOGIC
# =============================================================================

def merge_study_metadata(
    donors: pd.DataFrame,
    base_dir: Path,
    study: str,
    config: dict
) -> pd.DataFrame:
    """
    Merge reference metadata for a single study using verified mapping.
    """
    print(f"\n--- {study} ---")

    # Get donors for this study
    # Handle pal_* variants
    if study == 'pal':
        mask = donors['ihbca_dataset'].str.startswith('pal')
    else:
        mask = donors['ihbca_dataset'] == study

    study_donors = donors[mask]['ihbca_donor_id'].unique()
    print(f"iHBCA donors: {len(study_donors)}")

    # Skip if no reference file
    if config.get('file') is None:
        print(f"No reference CSV - skipping")
        return donors

    # Load mapping if exists
    mapping_path = base_dir / config['mapping_file'] if config.get('mapping_file') else None
    mapping = load_donor_mapping(mapping_path) if mapping_path else None

    if mapping:
        print(f"Using verified mapping: {len(mapping)} entries")
    else:
        print(f"Using direct ID match")

    # Load reference CSV
    ref_df = load_reference_csv(base_dir, study, config)
    if ref_df is None:
        return donors

    # Build join key
    ref_id_col = f'{study}_ref_donor_id'

    # For each iHBCA donor in this study, find the reference donor ID
    # Note: iHBCA patientIDs may have leading/trailing whitespace (e.g., Nee Controls)
    join_keys = []
    for ihbca_id in study_donors:
        lookup_id = ihbca_id.strip() if isinstance(ihbca_id, str) else ihbca_id
        if mapping:
            ref_id = mapping.get(lookup_id)
        else:
            ref_id = ihbca_id  # Direct match

        join_keys.append({
            'ihbca_donor_id': ihbca_id,
            ref_id_col: ref_id
        })

    join_df = pd.DataFrame(join_keys)

    # Count matches
    ref_ids = set(ref_df[ref_id_col].dropna().astype(str))
    matched = sum(1 for k in join_keys if k[ref_id_col] in ref_ids)
    print(f"Matched: {matched}/{len(study_donors)}")

    # Merge join keys with reference
    join_df = join_df.merge(ref_df, on=ref_id_col, how='left')

    # Merge back to donors
    donors = donors.merge(
        join_df,
        on='ihbca_donor_id',
        how='left',
        suffixes=('', '_dup')
    )

    # Drop duplicates
    donors = donors[[c for c in donors.columns if not c.endswith('_dup')]]

    return donors


def main():
    parser = argparse.ArgumentParser(
        description="Build unified donor metadata using verified mappings"
    )
    parser.add_argument(
        '--base-dir', type=str, default=DEFAULT_BASE_DIR,
        help="Base directory for harmonization"
    )
    parser.add_argument(
        '--ihbca-source', type=str, default=IHBCA_SOURCE,
        help="Path to iHBCA level1.5_annotations.csv"
    )
    parser.add_argument(
        '--output-dir', type=str, default=None,
        help="Output directory"
    )
    parser.add_argument(
        '--skip-parquet', action='store_true',
        help="Skip parquet output"
    )
    args = parser.parse_args()

    base_dir = Path(args.base_dir)
    ihbca_path = Path(args.ihbca_source)
    output_dir = Path(args.output_dir) if args.output_dir else base_dir / 'outputs' / 'unified_metadata'

    print("=" * 70)
    print("BUILD UNIFIED DONOR METADATA (PHASE 3)")
    print("=" * 70)
    print(f"Base: {base_dir}")
    print(f"iHBCA: {ihbca_path}")
    print(f"Output: {output_dir}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Get iHBCA donors
    donors = get_unique_donors_from_ihbca(ihbca_path)
    initial_count = len(donors)

    # Step 2: Merge each study's reference metadata
    print("\n=== Merging Study Reference Metadata ===")

    for study, config in REFERENCE_CONFIG.items():
        donors = merge_study_metadata(donors, base_dir, study, config)

    # Step 3: Verify
    print("\n=== Verification ===")
    final_count = len(donors)
    print(f"Initial donors: {initial_count}")
    print(f"Final donors: {final_count}")

    if final_count != initial_count:
        print(f"WARNING: Donor count changed during merge!")

    # Check for duplicates
    dups = donors[donors.duplicated(subset=['ihbca_donor_id'], keep=False)]
    if len(dups) > 0:
        print(f"WARNING: {len(dups)} duplicate ihbca_donor_id entries!")
        print(dups[['ihbca_donor_id', 'ihbca_dataset']].drop_duplicates())

    print(f"\nTotal columns: {len(donors.columns)}")

    # Step 4: Save
    print("\n=== Saving Outputs ===")

    csv_path = output_dir / 'unified_donor_metadata.csv'
    donors.to_csv(csv_path, index=False)
    print(f"CSV: {csv_path} ({csv_path.stat().st_size / 1024:.1f} KB)")

    if not args.skip_parquet:
        parquet_path = output_dir / 'unified_donor_metadata.parquet'
        donors.to_parquet(parquet_path, index=False)
        print(f"Parquet: {parquet_path} ({parquet_path.stat().st_size / 1024:.1f} KB)")

    # Summary
    summary = {
        'generated_at': datetime.now().isoformat(),
        'total_donors': len(donors),
        'total_columns': len(donors.columns),
        'donors_per_dataset': donors['ihbca_dataset'].value_counts().to_dict(),
        'studies_merged': list(REFERENCE_CONFIG.keys()),
    }

    summary_path = output_dir / 'donor_metadata_summary.yaml'
    with open(summary_path, 'w') as f:
        yaml.dump(summary, f, default_flow_style=False)

    print(f"\nSummary: {summary_path}")
    print("\n" + "=" * 70)
    print("PHASE 3 COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
