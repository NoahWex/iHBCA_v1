#!/usr/bin/env python3
"""
Map Gray component study cells to iHBCA reference.

Two-pass approach:
1. First pass: Discover batch→patient mapping by matching unique barcodes
2. Second pass: Use mapping to resolve all cells including duplicates

Key insight: Barcodes are only unique within a batch/patient, not globally.
Gray uses format `Human-{BRCA_status}-{batch}_{barcode}` where batch maps to patient.

Expected issues to verify (from knowledge.md):
- gray-rm-donor-mismatch: RM-* donors should appear, ~19K WT cells
"""

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import scanpy as sc
import pandas as pd
import yaml


def extract_barcode_from_gray(cell_id: str) -> tuple:
    """
    Extract barcode and batch info from Gray component cell ID.

    Format: Human-{BRCA_status}-{batch}_{barcode}
    Example: Human-BRCA1-A_AAACCTGAGACAAGCC

    Returns: (barcode, batch_letter, brca_status)
    """
    match = re.match(r'^Human-(BRCA1|BRCA2|WT)-([A-Z])_([ACGT]{16})$', cell_id)
    if match:
        brca_status = match.group(1)
        batch_letter = match.group(2)
        barcode = match.group(3)
        return (barcode, batch_letter, brca_status)
    return (None, None, None)


def extract_barcode_from_ihbca(cell_id: str) -> tuple:
    """
    Extract barcode and patient from iHBCA cell ID.

    Format: {patient}_{barcode}
    Example: PM-A_AAACCTGAGCAATATG

    Returns: (barcode, patient)
    """
    parts = cell_id.split('_')
    if len(parts) == 2:
        patient = parts[0]
        barcode = parts[1]
        if len(barcode) == 16 and all(c in 'ACGT' for c in barcode):
            return (barcode, patient)
    return (None, None)


def main():
    parser = argparse.ArgumentParser(description="Map Gray cells to iHBCA")
    parser.add_argument(
        "--gray-h5ad",
        default="${SOURCE_COMPONENT_STUDIES}/gray.h5ad",
        help="Gray component study h5ad"
    )
    parser.add_argument(
        "--ihbca-inventory",
        default="${SOURCE_IHBCAV1_HARMONIZATION%/}/outputs/ihbca_reference/ihbca_cell_inventory.csv",
        help="iHBCA cell inventory CSV"
    )
    parser.add_argument(
        "--output-dir",
        default="${SOURCE_IHBCAV1_HARMONIZATION%/}/studies/gray/outputs",
        help="Output directory"
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Gray Cell ID Mapping (Two-Pass) ===")
    print(f"Gray h5ad: {args.gray_h5ad}")
    print(f"iHBCA inventory: {args.ihbca_inventory}")
    print(f"Output: {output_dir}")
    print()

    # Step 1: Load iHBCA inventory for Gray cells
    print("Loading iHBCA cell inventory...")
    ihbca_inventory = pd.read_csv(args.ihbca_inventory)
    ihbca_gray = ihbca_inventory[ihbca_inventory['dataset'] == 'gray'].copy()
    print(f"  iHBCA Gray cells: {len(ihbca_gray):,}")

    # Build (patient, barcode) → ihbca_cell_id lookup
    ihbca_patient_barcode_lookup = {}  # (patient, barcode) → cell_id
    ihbca_patients = set()
    for _, row in ihbca_gray.iterrows():
        barcode, patient = extract_barcode_from_ihbca(row['cell_id'])
        if barcode and patient:
            key = (patient, barcode)
            ihbca_patient_barcode_lookup[key] = row['cell_id']
            ihbca_patients.add(patient)

    print(f"  Unique (patient, barcode) pairs: {len(ihbca_patient_barcode_lookup):,}")
    print(f"  Patients: {sorted(ihbca_patients)}")

    # Step 2: Load Gray component study
    print("\nLoading Gray component study...")
    adata = sc.read_h5ad(args.gray_h5ad)
    print(f"  Cells: {adata.n_obs:,}")
    print(f"  Genes: {adata.n_vars:,}")
    print(f"  Metadata columns: {list(adata.obs.columns)}")

    # Step 3: First pass - discover batch→patient mapping
    print("\nPass 1: Discovering batch→patient mapping...")

    # For each Gray batch, find which iHBCA patient(s) its barcodes match
    batch_to_patient_votes = defaultdict(Counter)  # batch_key → {patient: count}

    for cell_id in adata.obs_names:
        barcode, batch_letter, brca_status = extract_barcode_from_gray(cell_id)
        if barcode is None:
            continue

        batch_key = f"{brca_status}-{batch_letter}"

        # Check which patient(s) have this barcode
        for patient in ihbca_patients:
            if (patient, barcode) in ihbca_patient_barcode_lookup:
                batch_to_patient_votes[batch_key][patient] += 1

    # Determine majority patient for each batch
    batch_to_patient = {}
    print("\nBatch → Patient mapping:")
    for batch_key in sorted(batch_to_patient_votes.keys()):
        votes = batch_to_patient_votes[batch_key]
        if votes:
            majority_patient, count = votes.most_common(1)[0]
            total_votes = sum(votes.values())
            batch_to_patient[batch_key] = majority_patient

            # Check for ambiguity
            if len(votes) > 1:
                second_patient, second_count = votes.most_common(2)[1]
                if second_count > 10:  # Significant secondary
                    print(f"  {batch_key} → {majority_patient} ({count}/{total_votes}) "
                          f"[NOTE: {second_count} votes for {second_patient}]")
                else:
                    print(f"  {batch_key} → {majority_patient} ({count} cells)")
            else:
                print(f"  {batch_key} → {majority_patient} ({count} cells)")

    # Step 4: Second pass - map all cells using discovered mapping
    print("\nPass 2: Mapping all cells using batch→patient mapping...")

    mapped = []
    unmapped = []

    for cell_id in adata.obs_names:
        barcode, batch_letter, brca_status = extract_barcode_from_gray(cell_id)

        if barcode is None:
            unmapped.append({'gray_cell_id': cell_id, 'reason': 'parse_failed'})
            continue

        batch_key = f"{brca_status}-{batch_letter}"

        if batch_key not in batch_to_patient:
            unmapped.append({
                'gray_cell_id': cell_id,
                'gray_barcode': barcode,
                'gray_batch_letter': batch_letter,
                'gray_brca_status': brca_status,
                'reason': 'batch_not_mapped'
            })
            continue

        # Get expected patient for this batch
        expected_patient = batch_to_patient[batch_key]

        # Construct expected iHBCA cell ID
        expected_ihbca_id = f"{expected_patient}_{barcode}"

        # Verify it exists
        if (expected_patient, barcode) in ihbca_patient_barcode_lookup:
            actual_ihbca_id = ihbca_patient_barcode_lookup[(expected_patient, barcode)]
            mapped.append({
                'ihbca_cell_id': actual_ihbca_id,
                'gray_cell_id': cell_id,
                'gray_barcode': barcode,
                'gray_batch_letter': batch_letter,
                'gray_brca_status': brca_status,
                'gray_batch_key': batch_key,
                'ihbca_patient': expected_patient
            })
        else:
            unmapped.append({
                'gray_cell_id': cell_id,
                'gray_barcode': barcode,
                'gray_batch_letter': batch_letter,
                'gray_brca_status': brca_status,
                'gray_batch_key': batch_key,
                'expected_patient': expected_patient,
                'reason': 'not_in_ihbca_for_patient'
            })

    print(f"  Mapped: {len(mapped):,}")
    print(f"  Unmapped: {len(unmapped):,}")

    match_rate = 100 * len(mapped) / adata.n_obs
    print(f"  Match rate: {match_rate:.1f}%")

    # Step 5: Build output dataframe with all metadata
    print("\nBuilding output dataframe...")

    mapped_df = pd.DataFrame(mapped)

    # Add Gray metadata columns (prefixed)
    gray_obs = adata.obs.copy()
    gray_obs['gray_cell_id'] = gray_obs.index

    # Merge mapping with metadata
    output_df = mapped_df.merge(gray_obs, on='gray_cell_id', how='left')

    # Prefix all columns except ihbca_cell_id
    rename_cols = {}
    for col in output_df.columns:
        if col == 'ihbca_cell_id':
            continue
        elif col.startswith('gray_'):
            continue
        else:
            rename_cols[col] = f'gray_{col}'

    output_df = output_df.rename(columns=rename_cols)

    # Reorder columns: ihbca_cell_id first
    cols = ['ihbca_cell_id'] + [c for c in output_df.columns if c != 'ihbca_cell_id']
    output_df = output_df[cols]

    print(f"  Output columns: {len(output_df.columns)}")
    print(f"  Output rows: {len(output_df):,}")

    # Step 6: Write outputs
    cells_path = output_dir / 'gray_cells.csv'
    print(f"\nWriting cells to: {cells_path}")
    output_df.to_csv(cells_path, index=False)

    # Analyze unmapped reasons
    unmapped_reasons = Counter(u.get('reason', 'unknown') for u in unmapped)

    # Count cells by status
    wt_mapped = sum(1 for m in mapped if m['gray_brca_status'] == 'WT')
    brca1_mapped = sum(1 for m in mapped if m['gray_brca_status'] == 'BRCA1')
    brca2_mapped = sum(1 for m in mapped if m['gray_brca_status'] == 'BRCA2')

    # Write mapping findings
    findings = {
        'study': 'gray',
        'extraction_date': datetime.now().isoformat(),
        'source_file': args.gray_h5ad,
        'method': 'two_pass_batch_aware',
        'summary': {
            'total_component_cells': int(adata.n_obs),
            'total_ihbca_cells': len(ihbca_gray),
            'mapped': len(mapped),
            'unmapped': len(unmapped),
            'match_rate_pct': round(match_rate, 2)
        },
        'batch_to_patient_mapping': batch_to_patient,
        'unmapped_by_reason': dict(unmapped_reasons),
        'metadata_columns': list(adata.obs.columns),
        'output_columns': list(output_df.columns),
        'verification': {
            'expected_issue': 'gray-rm-donor-mismatch',
            'rm_patients_found': sorted([p for p in batch_to_patient.values() if p.startswith('RM-')]),
            'pm_patients_found': sorted([p for p in batch_to_patient.values() if p.startswith('PM-')]),
            'wt_cells_mapped': wt_mapped,
            'brca1_cells_mapped': brca1_mapped,
            'brca2_cells_mapped': brca2_mapped
        }
    }

    if unmapped:
        findings['unmapped_sample'] = unmapped[:20]

    findings_path = output_dir / 'mapping_findings.yaml'
    print(f"Writing findings to: {findings_path}")
    with open(findings_path, 'w') as f:
        yaml.dump(findings, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    # Print summary
    print("\n" + "=" * 60)
    print("MAPPING COMPLETE")
    print("=" * 60)
    print(f"\nTotal cells: {adata.n_obs:,}")
    print(f"Mapped: {len(mapped):,} ({match_rate:.1f}%)")
    print(f"Unmapped: {len(unmapped):,}")

    if unmapped_reasons:
        print(f"\nUnmapped by reason:")
        for reason, count in unmapped_reasons.most_common():
            print(f"  {reason}: {count:,}")

    print(f"\nPatients by type:")
    rm_patients = sorted(set(p for p in batch_to_patient.values() if p.startswith('RM-')))
    pm_patients = sorted(set(p for p in batch_to_patient.values() if p.startswith('PM-')))
    print(f"  PM (prophylactic mastectomy): {len(pm_patients)} ({', '.join(pm_patients)})")
    print(f"  RM (reduction mammoplasty): {len(rm_patients)} ({', '.join(rm_patients)})")

    print(f"\nCells by BRCA status (mapped):")
    print(f"  WT: {wt_mapped:,}")
    print(f"  BRCA1: {brca1_mapped:,}")
    print(f"  BRCA2: {brca2_mapped:,}")

    # Verify expected issue
    print(f"\n--- Verification: gray-rm-donor-mismatch ---")
    print(f"Expected: RM-* donors with ~19K WT cells")
    if rm_patients and wt_mapped > 15000:
        print(f"VERIFIED: Found {len(rm_patients)} RM donors with {wt_mapped:,} WT cells")
    else:
        print(f"CHECK: RM donors: {rm_patients}, WT cells: {wt_mapped:,}")

    print(f"\nOutputs:")
    print(f"  Cells: {cells_path}")
    print(f"  Findings: {findings_path}")

    # Why did the bioinformatician refuse to debug their code?
    # Because they were already in too many batches of trouble.


if __name__ == "__main__":
    main()
