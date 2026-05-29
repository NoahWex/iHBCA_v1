#!/usr/bin/env python3
"""
Step 6: Map Harmonized Donor Metadata to Cells (Phase A.2)

Joins harmonized_donor_metadata.csv to ihbca_cell_inventory.csv
to create cell-level metadata for DA analysis.

Input:
  - harmonized_donor_metadata.csv (287 donors)
  - ihbca_cell_inventory.csv (2.12M cells)

Output:
  - harmonized_cell_metadata.csv (2.12M cells with conditions)

# Every cell deserves to know where it came from. Existentially.
"""

import pandas as pd
from pathlib import Path


def main():
    # Paths
    base_dir = Path(__file__).parent.parent
    donor_path = base_dir / 'outputs' / 'harmonized_metadata' / 'harmonized_donor_metadata.csv'
    cell_path = base_dir / 'outputs' / 'ihbca_reference' / 'ihbca_cell_inventory.csv'
    output_dir = base_dir / 'outputs' / 'harmonized_metadata'

    # Load donor metadata
    print(f"Loading {donor_path}")
    donors = pd.read_csv(donor_path)
    print(f"  Loaded {len(donors)} donors")

    # Load cell inventory
    print(f"Loading {cell_path}")
    cells = pd.read_csv(cell_path)
    print(f"  Loaded {len(cells):,} cells")

    # Join on patient_id = ihbca_donor_id
    print("Joining donor metadata to cells...")
    merged = cells.merge(
        donors,
        left_on='patient_id',
        right_on='ihbca_donor_id',
        how='left'
    )

    # Validate join
    n_matched = merged['ihbca_donor_id'].notna().sum()
    n_unmatched = merged['ihbca_donor_id'].isna().sum()
    print(f"  Matched: {n_matched:,} cells ({100*n_matched/len(cells):.2f}%)")
    print(f"  Unmatched: {n_unmatched:,} cells ({100*n_unmatched/len(cells):.2f}%)")

    if n_unmatched > 0:
        # Show unmatched donors
        unmatched_patients = cells.loc[~cells['patient_id'].isin(donors['ihbca_donor_id']), 'patient_id'].unique()
        print(f"  Unmatched patient_ids: {list(unmatched_patients)[:10]}")

    # Select output columns
    output_cols = [
        'cell_id',
        'patient_id',  # keep original for reference
        'dataset',
        'batch',
        'study',
        'age_continuous',
        'age_binary',
        'parity_binary',
        'risk_status_binary',
        'menopausal_status_binary'
    ]
    output = merged[output_cols].copy()

    # Write output
    output_path = output_dir / 'harmonized_cell_metadata.csv'
    print(f"Writing {output_path}")
    output.to_csv(output_path, index=False)

    # Summary
    print("\n=== SUMMARY ===")
    print(f"Total cells: {len(output):,}")
    print(f"Unique donors: {output['patient_id'].nunique()}")

    print("\nCells by study:")
    for study, count in output['study'].value_counts().items():
        print(f"  {study}: {count:,}")

    print("\nCoverage per condition:")
    for col in ['age_continuous', 'age_binary', 'parity_binary', 'risk_status_binary', 'menopausal_status_binary']:
        n_with = output[col].notna().sum()
        print(f"  {col}: {n_with:,} cells ({100*n_with/len(output):.1f}%)")

    print(f"\nOutput saved to {output_path}")


if __name__ == '__main__':
    main()
