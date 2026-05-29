#!/usr/bin/env python3
"""
Extract iHBCA Reference Table from author_share annotations.

This is the iHBCA - the joint integration of all 7 component studies (~2.1M cells).
The output becomes the reference for mapping component study cell IDs.

Outputs:
- ihbca_reference_summary.yaml: Structure documentation
- ihbca_cell_inventory.csv: Full cell inventory with key columns

Why this script exists: The iHBCA annotations CSV from author_share is the
authoritative source for cell ID → study mapping. Component studies will be
mapped against this reference.
"""

import argparse
import csv
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import yaml


def main():
    parser = argparse.ArgumentParser(description="Extract iHBCA reference from author_share")
    parser.add_argument(
        "--input",
        default="${SOURCE_AUTHOR_SHARE}/ihbca_level1.5_annotations.csv",
        help="Path to ihbca_level1.5_annotations.csv"
    )
    parser.add_argument(
        "--output-dir",
        default="${SOURCE_IHBCAV1_HARMONIZATION%/}/outputs/ihbca_reference",
        help="Output directory"
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== iHBCA Reference Extraction ===")
    print(f"Input: {input_path}")
    print(f"Output: {output_dir}")
    print()

    # Verify input exists
    if not input_path.exists():
        print(f"ERROR: Input file not found: {input_path}")
        sys.exit(1)

    file_size_mb = input_path.stat().st_size / (1024 * 1024)
    print(f"Input file size: {file_size_mb:.1f} MB")

    # Read and process
    print("\nReading annotations...")

    # Containers for analysis
    column_info = {}  # col -> {dtype, n_unique, sample_values}
    dataset_counts = Counter()
    dataset_patient_counts = defaultdict(set)
    cell_ids = []
    sample_rows = []  # First 100 rows for inspection

    # Key columns to preserve in inventory
    key_columns = [
        'cellID', 'dataset', 'patientID', 'sampleID', 'batch',
        'level0', 'level1', 'level2',
        'level0_annotation', 'level1_annotation', 'level1.5_annotation',
        'age', 'parous', 'parity', 'risk_status', 'risk_status2'
    ]

    n_rows = 0
    with open(input_path, 'r') as f:
        reader = csv.DictReader(f)
        columns = reader.fieldnames

        # Initialize column tracking
        for col in columns:
            column_info[col] = {
                'values': Counter(),
                'n_null': 0
            }

        for row in reader:
            n_rows += 1

            # Track cell IDs and datasets
            cell_id = row.get('cellID', '')
            dataset = row.get('dataset', '')
            patient_id = row.get('patientID', '')

            cell_ids.append({
                'cell_id': cell_id,
                'dataset': dataset,
                'patient_id': patient_id,
                'batch': row.get('batch', '')
            })

            dataset_counts[dataset] += 1
            if patient_id:
                dataset_patient_counts[dataset].add(patient_id)

            # Track column values (capped to prevent memory explosion)
            for col in columns:
                val = row.get(col, '')
                if val == '' or val is None:
                    column_info[col]['n_null'] += 1
                else:
                    if len(column_info[col]['values']) < 1000:
                        column_info[col]['values'][val] += 1

            # Sample rows for inspection
            if n_rows <= 100:
                sample_rows.append({k: row.get(k, '') for k in key_columns if k in columns})

            # Progress
            if n_rows % 500000 == 0:
                print(f"  Processed {n_rows:,} rows...")

    print(f"  Total: {n_rows:,} cells")

    # Build summary
    print("\nBuilding summary...")

    # Column documentation
    column_docs = {}
    for col in columns:
        info = column_info[col]
        n_unique = len(info['values'])
        coverage_pct = 100 * (1 - info['n_null'] / n_rows) if n_rows > 0 else 0

        col_doc = {
            'position': columns.index(col) + 1,
            'n_unique': n_unique,
            'coverage_pct': round(coverage_pct, 1),
            'n_null': info['n_null']
        }

        # Include all values if < 50, else top 20
        if n_unique <= 50:
            col_doc['values'] = dict(info['values'].most_common())
        else:
            col_doc['sample_values'] = dict(info['values'].most_common(20))
            col_doc['note'] = f"Showing top 20 of {n_unique} unique values"

        column_docs[col] = col_doc

    # Dataset summary
    dataset_summary = {}
    for ds in sorted(dataset_counts.keys()):
        dataset_summary[ds] = {
            'n_cells': dataset_counts[ds],
            'n_patients': len(dataset_patient_counts[ds]),
            'patients': sorted(list(dataset_patient_counts[ds]))
        }

    # Generate summary YAML
    summary = {
        'extraction_metadata': {
            'source_file': str(input_path),
            'extracted_at': datetime.now().isoformat(),
            'file_size_mb': round(file_size_mb, 1),
            'total_cells': n_rows,
            'n_columns': len(columns),
            'n_datasets': len(dataset_counts)
        },
        'dataset_summary': dataset_summary,
        'columns': column_docs,
        'sample_rows': sample_rows[:10]  # Just 10 for inspection
    }

    # Write summary YAML
    summary_path = output_dir / 'ihbca_reference_summary.yaml'
    print(f"\nWriting summary to: {summary_path}")
    with open(summary_path, 'w') as f:
        yaml.dump(summary, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    # Write cell inventory CSV (key columns only for efficiency)
    # This gives us the ihbca_cell_id → dataset mapping needed for component study reconciliation
    inventory_path = output_dir / 'ihbca_cell_inventory.csv'
    print(f"Writing cell inventory to: {inventory_path}")

    with open(inventory_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['cell_id', 'dataset', 'patient_id', 'batch'])
        writer.writeheader()
        writer.writerows(cell_ids)

    # Print summary to stdout
    print("\n" + "=" * 60)
    print("EXTRACTION COMPLETE")
    print("=" * 60)
    print(f"\nTotal cells: {n_rows:,}")
    print(f"Total columns: {len(columns)}")
    print(f"\nCells per dataset:")
    for ds in sorted(dataset_counts.keys()):
        pct = 100 * dataset_counts[ds] / n_rows
        n_patients = len(dataset_patient_counts[ds])
        print(f"  {ds:20s}: {dataset_counts[ds]:>10,} ({pct:5.1f}%) - {n_patients} patients")

    print(f"\nOutputs:")
    print(f"  Summary: {summary_path}")
    print(f"  Cell inventory: {inventory_path}")

    # A small Easter egg hidden in the logs, as per the project rule
    # "Every new script gets one stupid joke hidden in the comments somewhere."
    # It's in this comment: Why do bioinformaticians hate karaoke?
    # Because they can't handle all the single-cell tracks.


if __name__ == "__main__":
    main()
