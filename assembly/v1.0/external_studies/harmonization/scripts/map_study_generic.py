#!/usr/bin/env python3
"""
Generic Cell ID Mapping Script for Phase A.

Maps any component study to iHBCA reference based on study_manifest.yaml config.
Handles multiple mapping types: identity, batch_aware, identity_with_normalization, custom_pal.

Usage:
    python map_study_generic.py --study kumar
    python map_study_generic.py --study gray
    python map_study_generic.py --study pal_norm_epi

Outputs per study:
    - {study}_cells.csv: All mapped cells with {study}_* prefixed columns
    - mapping_findings.yaml: Comprehensive mapping report
"""

import argparse
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


def load_manifest(manifest_path: str) -> dict:
    """Load study manifest configuration."""
    with open(manifest_path, 'r') as f:
        return yaml.safe_load(f)


def load_ihbca_inventory(path: str, dataset_filter: str) -> pd.DataFrame:
    """Load iHBCA inventory and filter by dataset."""
    print(f"Loading iHBCA inventory...")
    df = pd.read_csv(path)
    filtered = df[df['dataset'] == dataset_filter].copy()
    print(f"  Total iHBCA cells: {len(df):,}")
    print(f"  {dataset_filter} cells: {len(filtered):,}")
    return filtered


def load_component_study(path: str, file_format: str):
    """Load component study and return cell IDs + metadata."""
    print(f"Loading component study: {path}")

    if file_format == 'h5ad':
        import scanpy as sc
        adata = sc.read_h5ad(path)
        cell_ids = list(adata.obs_names)
        metadata = adata.obs.copy()
        metadata['_cell_id'] = metadata.index
        print(f"  Cells: {len(cell_ids):,}")
        print(f"  Metadata columns: {len(metadata.columns)}")
        return cell_ids, metadata

    elif file_format == 'rds':
        # For RDS, we need to use R via rpy2 or read pre-extracted data
        # For now, implement a fallback that reads from Track 2 extraction
        # or use subprocess to call R
        raise NotImplementedError(
            f"RDS loading not yet implemented in generic mapper. "
            f"Use study-specific script or add rpy2 support."
        )

    else:
        raise ValueError(f"Unknown format: {file_format}")


def map_identity(component_ids: list, ihbca_df: pd.DataFrame) -> tuple:
    """
    Identity mapping: cell IDs should match exactly.

    Returns: (mapped_list, unmapped_list)
    """
    ihbca_ids = set(ihbca_df['cell_id'].values)

    mapped = []
    unmapped = []

    for cell_id in component_ids:
        if cell_id in ihbca_ids:
            mapped.append({
                'ihbca_cell_id': cell_id,
                'component_cell_id': cell_id,
                'mapping_method': 'identity'
            })
        else:
            unmapped.append({
                'component_cell_id': cell_id,
                'reason': 'not_in_ihbca'
            })

    return mapped, unmapped


def map_identity_with_normalization(component_ids: list, ihbca_df: pd.DataFrame,
                                     reversed_batches: list) -> tuple:
    """
    Identity mapping with format normalization for some cells.

    Some cells have reversed format (BARCODE-BATCH instead of BATCH_BARCODE).
    Try identity first, then try normalized format.
    """
    ihbca_ids = set(ihbca_df['cell_id'].values)

    mapped = []
    unmapped = []

    for cell_id in component_ids:
        # Try identity first
        if cell_id in ihbca_ids:
            mapped.append({
                'ihbca_cell_id': cell_id,
                'component_cell_id': cell_id,
                'mapping_method': 'identity'
            })
            continue

        # Try normalized format (for reversed cells)
        # Twigger reversed: BARCODE-BATCH should become BATCH_BARCODE-1
        # Check if this looks like a reversed format
        normalized = None
        if '-' in cell_id and '_' not in cell_id:
            # Might be reversed: AAACCCAAGAGGTCAC-RB5 -> RB5_AAACCCAAGAGGTCAC-1
            parts = cell_id.rsplit('-', 1)
            if len(parts) == 2:
                barcode_part, batch_part = parts
                # Check if batch_part is a known reversed batch
                if any(batch_part.startswith(rb) for rb in reversed_batches):
                    normalized = f"{batch_part}_{barcode_part}-1"

        if normalized and normalized in ihbca_ids:
            mapped.append({
                'ihbca_cell_id': normalized,
                'component_cell_id': cell_id,
                'mapping_method': 'normalized_reversed'
            })
        else:
            unmapped.append({
                'component_cell_id': cell_id,
                'reason': 'not_in_ihbca_after_normalization'
            })

    return mapped, unmapped


def map_batch_aware(component_ids: list, ihbca_df: pd.DataFrame,
                    cell_id_regex: str) -> tuple:
    """
    Two-pass batch-aware mapping (as implemented for Gray).

    Pass 1: Discover batch→patient mapping via majority vote
    Pass 2: Map all cells using discovered mapping
    """
    # Build (patient, barcode) → ihbca_cell_id lookup
    ihbca_patient_barcode_lookup = {}
    ihbca_patients = set()

    for _, row in ihbca_df.iterrows():
        cell_id = row['cell_id']
        patient = row['patient_id']
        # Extract barcode (last 16 chars before any suffix)
        parts = cell_id.split('_')
        if len(parts) >= 2:
            barcode_part = parts[-1]
            # Remove -1 suffix if present
            barcode = barcode_part.rstrip('-1').replace('-1', '')
            if len(barcode) >= 16:
                barcode = barcode[:16] if len(barcode) > 16 else barcode
            key = (patient, barcode)
            ihbca_patient_barcode_lookup[key] = cell_id
            ihbca_patients.add(patient)

    print(f"  iHBCA patients: {sorted(ihbca_patients)}")

    # Parse component cell IDs
    pattern = re.compile(cell_id_regex)
    parsed_cells = []

    for cell_id in component_ids:
        match = pattern.match(cell_id)
        if match:
            groups = match.groups()
            # Assume format: (status, batch, barcode) for Gray-like
            if len(groups) == 3:
                status, batch, barcode = groups
                batch_key = f"{status}-{batch}"
                parsed_cells.append({
                    'cell_id': cell_id,
                    'batch_key': batch_key,
                    'barcode': barcode
                })

    # Pass 1: Discover batch→patient mapping
    batch_to_patient_votes = defaultdict(Counter)

    for cell in parsed_cells:
        barcode = cell['barcode']
        batch_key = cell['batch_key']

        for patient in ihbca_patients:
            if (patient, barcode) in ihbca_patient_barcode_lookup:
                batch_to_patient_votes[batch_key][patient] += 1

    batch_to_patient = {}
    print("\n  Batch → Patient mapping:")
    for batch_key in sorted(batch_to_patient_votes.keys()):
        votes = batch_to_patient_votes[batch_key]
        if votes:
            majority_patient, count = votes.most_common(1)[0]
            batch_to_patient[batch_key] = majority_patient
            print(f"    {batch_key} → {majority_patient} ({count} cells)")

    # Pass 2: Map all cells
    mapped = []
    unmapped = []

    for cell in parsed_cells:
        batch_key = cell['batch_key']
        barcode = cell['barcode']
        cell_id = cell['cell_id']

        if batch_key not in batch_to_patient:
            unmapped.append({
                'component_cell_id': cell_id,
                'reason': 'batch_not_mapped'
            })
            continue

        expected_patient = batch_to_patient[batch_key]

        if (expected_patient, barcode) in ihbca_patient_barcode_lookup:
            ihbca_cell_id = ihbca_patient_barcode_lookup[(expected_patient, barcode)]
            mapped.append({
                'ihbca_cell_id': ihbca_cell_id,
                'component_cell_id': cell_id,
                'mapping_method': 'batch_aware',
                'batch_key': batch_key,
                'patient': expected_patient
            })
        else:
            unmapped.append({
                'component_cell_id': cell_id,
                'batch_key': batch_key,
                'expected_patient': expected_patient,
                'reason': 'not_in_ihbca_for_patient'
            })

    return mapped, unmapped, batch_to_patient


def generate_report(study: str, config: dict, mapped: list, unmapped: list,
                    metadata_columns: list, extra_info: dict = None) -> dict:
    """Generate comprehensive mapping report."""
    total = len(mapped) + len(unmapped)
    match_rate = 100 * len(mapped) / total if total > 0 else 0
    expected_rate = config.get('expected_match_rate', 100.0)

    # Analyze unmapped reasons
    unmapped_reasons = Counter(u.get('reason', 'unknown') for u in unmapped)

    report = {
        'study': study,
        'extraction_date': datetime.now().isoformat(),
        'source_file': config['source_path'],
        'mapping_type': config['mapping_type'],
        'summary': {
            'total_component_cells': total,
            'mapped': len(mapped),
            'unmapped': len(unmapped),
            'match_rate_pct': round(match_rate, 2),
            'expected_match_rate_pct': expected_rate,
            'rate_within_tolerance': abs(match_rate - expected_rate) < 0.5
        },
        'validation': {
            'expected_issues': config.get('expected_issues', []),
            'verified': [],
            'unexpected_gaps': []
        },
        'unmapped_by_reason': dict(unmapped_reasons),
        'metadata_columns': metadata_columns
    }

    if extra_info:
        report.update(extra_info)

    if unmapped:
        report['unmapped_sample'] = unmapped[:20]

    # Check for unexpected gaps
    if match_rate < expected_rate - 0.5:
        report['validation']['unexpected_gaps'].append({
            'issue': 'match_rate_below_expected',
            'expected': expected_rate,
            'actual': round(match_rate, 2),
            'gap': round(expected_rate - match_rate, 2)
        })

    return report


def print_report(report: dict):
    """Print formatted report to stdout."""
    print("\n" + "=" * 70)
    print(f"{report['study'].upper()} CELL ID MAPPING REPORT")
    print("=" * 70)

    s = report['summary']
    print(f"\nSUMMARY")
    print(f"  Total component cells: {s['total_component_cells']:,}")
    print(f"  Mapped:                {s['mapped']:,}")
    print(f"  Unmapped:              {s['unmapped']:,}")
    print(f"  Match rate:            {s['match_rate_pct']}%")
    print(f"  Expected rate:         {s['expected_match_rate_pct']}%")

    if s['rate_within_tolerance']:
        print(f"  Status:                ✓ WITHIN TOLERANCE")
    else:
        print(f"  Status:                ✗ BELOW EXPECTED")

    if report.get('unmapped_by_reason'):
        print(f"\nUNMAPPED BY REASON")
        for reason, count in report['unmapped_by_reason'].items():
            print(f"  {reason}: {count:,}")

    if report.get('batch_to_patient_mapping'):
        print(f"\nBATCH → PATIENT MAPPING")
        for batch, patient in sorted(report['batch_to_patient_mapping'].items()):
            print(f"  {batch} → {patient}")

    v = report['validation']
    if v['expected_issues']:
        print(f"\nEXPECTED ISSUES")
        for issue in v['expected_issues']:
            status = "✓ verified" if issue in v['verified'] else "? pending"
            print(f"  {issue}: {status}")

    if v['unexpected_gaps']:
        print(f"\n⚠️  UNEXPECTED GAPS")
        for gap in v['unexpected_gaps']:
            print(f"  {gap['issue']}: expected {gap['expected']}%, got {gap['actual']}%")

    print("\n" + "=" * 70)


def main():
    parser = argparse.ArgumentParser(description="Generic cell ID mapping")
    parser.add_argument("--study", required=True, help="Study name from manifest")
    parser.add_argument(
        "--manifest",
        default="${SOURCE_IHBCAV1_HARMONIZATION%/}/config/study_manifest.yaml",
        help="Path to study manifest"
    )
    parser.add_argument("--dry-run", action="store_true", help="Print config only, don't map")
    args = parser.parse_args()

    # Load manifest
    manifest = load_manifest(args.manifest)
    shared = manifest['shared']

    if args.study not in manifest['studies']:
        print(f"ERROR: Study '{args.study}' not in manifest")
        print(f"Available: {list(manifest['studies'].keys())}")
        sys.exit(1)

    config = manifest['studies'][args.study]

    print(f"=== Mapping {args.study} ===")
    print(f"Source: {config['source_path']}")
    print(f"Format: {config['format']}")
    print(f"Mapping type: {config['mapping_type']}")
    print(f"Expected cells: {config['expected_cells']:,}")
    print(f"Expected match rate: {config['expected_match_rate']}%")
    print()

    if args.dry_run:
        print("DRY RUN - exiting without mapping")
        return

    # Load iHBCA inventory
    ihbca_df = load_ihbca_inventory(shared['ihbca_inventory'], config['ihbca_dataset'])

    # Load component study
    try:
        component_ids, metadata = load_component_study(config['source_path'], config['format'])
    except NotImplementedError as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    # Map based on type
    mapping_type = config['mapping_type']
    extra_info = {}

    if mapping_type == 'identity':
        mapped, unmapped = map_identity(component_ids, ihbca_df)

    elif mapping_type == 'identity_with_normalization':
        reversed_batches = config.get('normalization', {}).get('reversed_batches', [])
        mapped, unmapped = map_identity_with_normalization(
            component_ids, ihbca_df, reversed_batches
        )

    elif mapping_type == 'batch_aware':
        regex = config.get('cell_id_regex')
        if not regex:
            print(f"ERROR: batch_aware mapping requires cell_id_regex")
            sys.exit(1)
        mapped, unmapped, batch_mapping = map_batch_aware(component_ids, ihbca_df, regex)
        extra_info['batch_to_patient_mapping'] = batch_mapping

    elif mapping_type == 'custom_pal':
        print(f"ERROR: custom_pal mapping not yet implemented in generic mapper")
        sys.exit(1)

    else:
        print(f"ERROR: Unknown mapping type: {mapping_type}")
        sys.exit(1)

    print(f"\nMapping results:")
    print(f"  Mapped: {len(mapped):,}")
    print(f"  Unmapped: {len(unmapped):,}")

    # Build output dataframe
    print("\nBuilding output dataframe...")
    mapped_df = pd.DataFrame(mapped)

    # Add metadata
    metadata['_cell_id'] = metadata.index
    merged = mapped_df.merge(
        metadata,
        left_on='component_cell_id',
        right_on='_cell_id',
        how='left'
    )

    # Prefix columns with study name
    rename_cols = {}
    for col in merged.columns:
        if col == 'ihbca_cell_id':
            continue
        elif col.startswith(f'{args.study}_'):
            continue
        else:
            rename_cols[col] = f'{args.study}_{col}'

    merged = merged.rename(columns=rename_cols)

    # Reorder: ihbca_cell_id first
    cols = ['ihbca_cell_id'] + [c for c in merged.columns if c != 'ihbca_cell_id']
    merged = merged[cols]

    # Write outputs
    output_dir = Path(shared['output_base']) / args.study / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)

    cells_path = output_dir / f'{args.study}_cells.csv'
    print(f"\nWriting cells to: {cells_path}")
    merged.to_csv(cells_path, index=False)

    # Generate and write report
    report = generate_report(
        args.study, config, mapped, unmapped,
        list(metadata.columns), extra_info
    )

    findings_path = output_dir / 'mapping_findings.yaml'
    print(f"Writing findings to: {findings_path}")
    with open(findings_path, 'w') as f:
        yaml.dump(report, f, default_flow_style=False, sort_keys=False)

    # Print report
    print_report(report)

    # Exit with error if unexpected gaps
    if report['validation']['unexpected_gaps']:
        print("\n⚠️  FAILED: Unexpected gaps found. Investigation required.")
        sys.exit(1)

    print(f"\n✓ {args.study} mapping complete")


if __name__ == "__main__":
    main()
