#!/usr/bin/env python3
"""
Step 3: Coverage Analysis for Donor Metadata Harmonization (Phase A.2)

Compute cross-study coverage matrix for each semantic cluster.
Determine which conditions have sufficient data for DA testing.

Output: coverage_matrix.yaml

# Coverage: the forgotten middle child of completeness and accuracy.
"""

import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from collections import defaultdict
import re


def is_missing(value) -> bool:
    """Check if value is effectively missing."""
    if pd.isna(value):
        return True
    str_val = str(value).strip().lower()
    return str_val in ['', 'nan', 'none', 'na', 'n/a', '-', 'unknown']


def get_study_from_donor(donor_id: str, df: pd.DataFrame) -> str:
    """Determine which study a donor belongs to based on ihbca_dataset."""
    row = df[df['ihbca_donor_id'] == donor_id]
    if len(row) > 0:
        dataset = row['ihbca_dataset'].values[0]
        if pd.notna(dataset):
            return dataset.lower()
    return 'unknown'


def parse_nee_parity(value: str) -> int | None:
    """Parse GXPX notation to extract parity (P value)."""
    if is_missing(value):
        return None
    # Match patterns like G0P0M0A0, G3P3M0A0, G2P1
    match = re.search(r'P(\d+)', str(value))
    if match:
        return int(match.group(1))
    return None


def analyze_age_coverage(df: pd.DataFrame) -> dict:
    """Analyze age coverage by study."""
    studies = ['gray', 'murrow', 'nee', 'twigger', 'pal', 'kumar', 'reed']
    coverage = {}

    for study in studies:
        study_donors = df[df['ihbca_dataset'] == study]
        n_total = len(study_donors)

        if n_total == 0:
            continue

        # Continuous age
        age_col = f'{study}_Age' if study != 'ihbca' else 'ihbca_age'
        if study == 'pal':
            age_col = 'pal_patient_age'
        elif study == 'twigger':
            age_col = 'twigger_Maternal_age'
        elif study == 'reed':
            age_col = 'ihbca_age'

        if age_col in df.columns:
            non_missing = study_donors[age_col].apply(lambda x: not is_missing(x))

            # Check if continuous or binary
            if study == 'kumar':
                # Kumar has binary only (Y/O)
                coverage[study] = {
                    'n_donors': int(n_total),
                    'continuous_available': False,
                    'binary_available': True,
                    'n_with_age': int(non_missing.sum()),
                    'coverage': round(non_missing.sum() / n_total, 3),
                    'encoding': 'binary (Y=<50, O>=50)'
                }
            else:
                # Try to parse as numeric
                values = study_donors.loc[non_missing, age_col]
                try:
                    numeric_values = pd.to_numeric(values, errors='coerce')
                    valid_numeric = numeric_values.dropna()
                    coverage[study] = {
                        'n_donors': int(n_total),
                        'continuous_available': True,
                        'binary_available': True,
                        'n_with_age': int(len(valid_numeric)),
                        'coverage': round(len(valid_numeric) / n_total, 3),
                        'age_range': [float(valid_numeric.min()), float(valid_numeric.max())] if len(valid_numeric) > 0 else None,
                        'encoding': 'numeric'
                    }
                except:
                    coverage[study] = {
                        'n_donors': int(n_total),
                        'continuous_available': False,
                        'binary_available': False,
                        'n_with_age': 0,
                        'coverage': 0.0,
                        'encoding': 'unknown'
                    }
        else:
            coverage[study] = {
                'n_donors': int(n_total),
                'continuous_available': False,
                'binary_available': False,
                'n_with_age': 0,
                'coverage': 0.0,
                'encoding': 'not collected'
            }

    # Summary
    total_donors = sum(c['n_donors'] for c in coverage.values())
    continuous_donors = sum(c['n_with_age'] for c in coverage.values() if c.get('continuous_available', False))
    binary_donors = sum(c['n_with_age'] for c in coverage.values() if c.get('binary_available', False))

    return {
        'per_study': coverage,
        'summary': {
            'total_donors': total_donors,
            'continuous_donors': continuous_donors,
            'binary_donors': binary_donors,
            'continuous_studies': [s for s, c in coverage.items() if c.get('continuous_available', False)],
            'binary_studies': [s for s, c in coverage.items() if c.get('binary_available', False)]
        }
    }


def analyze_parity_coverage(df: pd.DataFrame) -> dict:
    """Analyze parity coverage by study."""
    studies = ['gray', 'murrow', 'nee', 'twigger', 'pal', 'kumar', 'reed']
    coverage = {}

    for study in studies:
        study_donors = df[df['ihbca_dataset'] == study]
        n_total = len(study_donors)

        if n_total == 0:
            continue

        # Different parity columns by study
        parity_col = None
        if study == 'gray':
            parity_col = 'gray_Births'  # or gray_Pregnancies
        elif study == 'murrow':
            parity_col = 'murrow_Parity'
        elif study == 'nee':
            parity_col = 'nee_Parity'  # GXPX format
        elif study == 'twigger':
            parity_col = 'twigger_Parity'
        elif study == 'pal':
            parity_col = 'pal_parity'  # Nulliparous/Parous
        elif study == 'kumar':
            parity_col = 'kumar_Parity'
        elif study == 'reed':
            # Reed uses ihbca_parity from the master column
            parity_col = 'ihbca_parity'

        if parity_col and parity_col in df.columns:
            non_missing = study_donors[parity_col].apply(lambda x: not is_missing(x))
            values = study_donors.loc[non_missing, parity_col]

            # Check encoding
            if study == 'nee':
                # Parse GXPX
                parsed = values.apply(parse_nee_parity)
                n_valid = parsed.notna().sum()
                encoding = 'GXPX notation'
            elif study == 'pal':
                encoding = 'binary (Nulliparous/Parous)'
                n_valid = non_missing.sum()
            else:
                encoding = 'numeric/categorical'
                n_valid = non_missing.sum()

            coverage[study] = {
                'n_donors': int(n_total),
                'n_with_parity': int(n_valid),
                'coverage': round(n_valid / n_total, 3),
                'encoding': encoding,
                'confounded': study.startswith('pal')  # Pal is confounded with menopausal
            }
        else:
            coverage[study] = {
                'n_donors': int(n_total),
                'n_with_parity': 0,
                'coverage': 0.0,
                'encoding': 'not collected',
                'confounded': False
            }

    # Summary
    total_donors = sum(c['n_donors'] for c in coverage.values())
    with_parity = sum(c['n_with_parity'] for c in coverage.values() if not c.get('confounded', False))

    return {
        'per_study': coverage,
        'summary': {
            'total_donors': total_donors,
            'donors_with_parity': with_parity,
            'testable_studies': [s for s, c in coverage.items()
                                if c.get('n_with_parity', 0) > 0 and not c.get('confounded', False)],
            'confounded_studies': [s for s, c in coverage.items() if c.get('confounded', False)]
        }
    }


def analyze_risk_coverage(df: pd.DataFrame) -> dict:
    """Analyze risk/BRCA status coverage by study."""
    studies = ['gray', 'murrow', 'nee', 'twigger', 'pal', 'kumar', 'reed']
    coverage = {}

    for study in studies:
        study_donors = df[df['ihbca_dataset'] == study]
        n_total = len(study_donors)

        if n_total == 0:
            continue

        # Use ihbca_risk_status as master
        risk_values = study_donors['ihbca_risk_status']
        non_missing = risk_values.apply(lambda x: not is_missing(x))
        values = risk_values[non_missing]

        # Count AR vs HR
        n_ar = (values == 'AR').sum()
        n_hr = values.str.startswith('HR').sum() if len(values) > 0 else 0

        # Check if study has contrast
        has_contrast = n_ar > 0 and n_hr > 0

        coverage[study] = {
            'n_donors': int(n_total),
            'n_with_risk': int(non_missing.sum()),
            'n_ar': int(n_ar),
            'n_hr': int(n_hr),
            'coverage': round(non_missing.sum() / n_total, 3),
            'has_contrast': has_contrast,
            'testable': has_contrast
        }

    return {
        'per_study': coverage,
        'summary': {
            'total_donors': sum(c['n_donors'] for c in coverage.values()),
            'testable_studies': [s for s, c in coverage.items() if c.get('testable', False)],
            'ar_only_studies': [s for s, c in coverage.items()
                               if c.get('n_ar', 0) > 0 and c.get('n_hr', 0) == 0],
            'hr_only_studies': [s for s, c in coverage.items()
                               if c.get('n_ar', 0) == 0 and c.get('n_hr', 0) > 0]
        }
    }


def analyze_menopausal_coverage(df: pd.DataFrame) -> dict:
    """Analyze menopausal status coverage by study."""
    studies = ['gray', 'murrow', 'nee', 'twigger', 'pal', 'kumar', 'reed']
    coverage = {}

    for study in studies:
        study_donors = df[df['ihbca_dataset'] == study]
        n_total = len(study_donors)

        if n_total == 0:
            continue

        # Murrow: all donors premenopausal (age 19-47), no raw column
        if study == 'murrow':
            coverage[study] = {
                'n_donors': int(n_total),
                'n_with_status': int(n_total),
                'n_pre': int(n_total),
                'n_post': 0,
                'coverage': 1.0,
                'has_contrast': False,
                'testable': False,
                'notes': 'Inferred from age (all 19-47). No contrast.'
            }
            continue

        # Different columns by study
        meno_col = None
        if study == 'gray':
            meno_col = 'gray_Menopause'
        elif study == 'kumar':
            meno_col = 'kumar_Menopause'
        elif study == 'nee':
            meno_col = 'nee_Menstrual Status'
        elif study == 'pal':
            meno_col = 'pal_menopausal_status'
        elif study == 'reed':
            meno_col = 'reed_Menopause_status'

        if meno_col and meno_col in df.columns:
            non_missing = study_donors[meno_col].apply(lambda x: not is_missing(x))
            values = study_donors.loc[non_missing, meno_col]

            # Count pre vs post
            values_lower = values.astype(str).str.lower()
            n_pre = values_lower.str.contains('pre').sum()
            n_post = (values_lower.str.contains('post') | values_lower.str.contains('surgical') |
                     values_lower.str.contains('menopausal')).sum()

            has_contrast = n_pre > 0 and n_post > 0

            coverage[study] = {
                'n_donors': int(n_total),
                'n_with_status': int(non_missing.sum()),
                'n_pre': int(n_pre),
                'n_post': int(n_post),
                'coverage': round(non_missing.sum() / n_total, 3),
                'has_contrast': has_contrast,
                'testable': has_contrast
            }
        else:
            coverage[study] = {
                'n_donors': int(n_total),
                'n_with_status': 0,
                'n_pre': 0,
                'n_post': 0,
                'coverage': 0.0,
                'has_contrast': False,
                'testable': False
            }

    return {
        'per_study': coverage,
        'summary': {
            'total_donors': sum(c['n_donors'] for c in coverage.values()),
            'testable_studies': [s for s, c in coverage.items() if c.get('testable', False)]
        }
    }


def main():
    # Paths
    base_dir = Path(__file__).parent.parent
    input_path = base_dir / 'outputs' / 'unified_metadata' / 'unified_donor_metadata.csv'
    output_dir = base_dir / 'outputs' / 'harmonization_audit'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load data
    print(f"Loading {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded: {len(df)} donors × {len(df.columns)} columns")

    # Analyze each concept
    print("Analyzing coverage...")

    coverage_matrix = {
        'age': analyze_age_coverage(df),
        'parity': analyze_parity_coverage(df),
        'risk_status': analyze_risk_coverage(df),
        'menopausal_status': analyze_menopausal_coverage(df)
    }

    # Testability summary
    testability = {}
    for concept, data in coverage_matrix.items():
        summary = data.get('summary', {})
        testability[concept] = {
            'native_da': summary.get('testable_studies', summary.get('continuous_studies', [])),
            'pooled_da': len(summary.get('testable_studies', summary.get('binary_studies', []))) >= 2,
            'notes': []
        }

        # Add condition-specific notes
        if concept == 'age':
            testability[concept]['notes'].append("Kumar has binary only (Y/O) - excluded from continuous tests")
        elif concept == 'parity':
            if summary.get('confounded_studies'):
                testability[concept]['notes'].append(f"CONFOUNDED in {summary['confounded_studies']} - skip DA")
        elif concept == 'risk_status':
            if summary.get('ar_only_studies'):
                testability[concept]['notes'].append(f"AR-only (no contrast): {summary['ar_only_studies']}")

    # Build output
    output = {
        'summary': {
            'n_donors': len(df),
            'studies': df['ihbca_dataset'].value_counts().to_dict()
        },
        'coverage_matrix': coverage_matrix,
        'testability': testability
    }

    # Write YAML
    output_path = output_dir / 'coverage_matrix.yaml'
    print(f"Writing {output_path}")

    with open(output_path, 'w') as f:
        yaml.dump(output, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"Done. Coverage matrix saved to {output_path}")

    # Print summary
    print("\n=== TESTABILITY SUMMARY ===")
    for concept, info in testability.items():
        print(f"\n{concept.upper()}:")
        print(f"  Native DA: {info['native_da']}")
        print(f"  Pooled DA possible: {info['pooled_da']}")
        for note in info['notes']:
            print(f"  Note: {note}")


if __name__ == '__main__':
    main()
