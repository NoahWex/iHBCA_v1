#!/usr/bin/env python3
"""
Step 1: Value Inventory for Donor Metadata Harmonization (Phase A.2)

Enumerate all columns in unified_donor_metadata.csv with:
- Unique values, frequencies, dtypes, missing counts
- Per-study breakdown where column has study prefix

Output: column_inventory.yaml

# Don't audit the auditor, the auditor audits you.
"""

import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from collections import Counter
import re


def get_study_prefix(col_name: str) -> str | None:
    """Extract study prefix from column name."""
    prefixes = ['ihbca_', 'gray_', 'murrow_', 'pal_', 'twigger_', 'nee_', 'kumar_']
    for prefix in prefixes:
        if col_name.startswith(prefix):
            return prefix.rstrip('_')
    return None


def analyze_column(df: pd.DataFrame, col_name: str) -> dict:
    """Analyze a single column for the inventory."""
    series = df[col_name]

    # Basic stats
    n_total = len(series)
    n_missing = series.isna().sum()
    # Also count empty strings and 'nan' strings as missing
    string_missing = series.apply(lambda x: str(x).strip().lower() in ['', 'nan', 'none', 'na', 'n/a'] if pd.notna(x) else False).sum()
    n_effective_missing = n_missing + string_missing

    # Get non-missing values
    non_missing = series.dropna()
    non_missing_clean = non_missing[~non_missing.apply(lambda x: str(x).strip().lower() in ['', 'nan', 'none', 'na', 'n/a'])]

    # Unique values
    unique_values = non_missing_clean.unique()
    n_unique = len(unique_values)

    # Coverage
    coverage = (n_total - n_effective_missing) / n_total if n_total > 0 else 0.0

    # Value frequencies (top 20)
    value_counts = Counter(non_missing_clean.astype(str))
    top_values = dict(value_counts.most_common(20))

    # Sample values (up to 10 distinct)
    sample_values = [str(v) for v in list(unique_values)[:10]]

    # Infer dtype
    dtype_raw = str(series.dtype)
    dtype_inferred = infer_semantic_dtype(non_missing_clean, dtype_raw)

    # Study prefix
    study_prefix = get_study_prefix(col_name)

    result = {
        'dtype_raw': dtype_raw,
        'dtype_inferred': dtype_inferred,
        'n_unique': int(n_unique),
        'n_missing': int(n_missing),
        'n_effective_missing': int(n_effective_missing),
        'coverage': round(coverage, 3),
        'sample_values': sample_values,
        'value_frequencies': top_values,
    }

    if study_prefix:
        result['study_prefix'] = study_prefix
        # Get concept name (column name without prefix)
        result['concept'] = col_name[len(study_prefix) + 1:]

    return result


def infer_semantic_dtype(series: pd.Series, raw_dtype: str) -> str:
    """Infer semantic data type from values."""
    if len(series) == 0:
        return 'empty'

    # Try numeric conversion
    try:
        numeric = pd.to_numeric(series, errors='raise')
        if all(numeric.dropna() == numeric.dropna().astype(int)):
            return 'integer'
        return 'numeric'
    except (ValueError, TypeError):
        pass

    # Check for boolean-like
    str_values = set(series.astype(str).str.lower().str.strip())
    bool_like = {'true', 'false', 'yes', 'no', '0', '1', 't', 'f', 'y', 'n'}
    if str_values.issubset(bool_like | {'nan', '', 'none', 'na'}):
        return 'boolean'

    # Check for binary categories (2 main values)
    value_counts = series.value_counts()
    if len(value_counts) == 2:
        return 'binary_categorical'
    elif len(value_counts) <= 10:
        return 'categorical'

    # Otherwise string/mixed
    return 'string'


def group_by_concept(inventory: dict) -> dict:
    """Group columns by semantic concept across studies."""
    concept_map = {}

    for col_name, col_info in inventory.items():
        if 'concept' in col_info:
            concept = col_info['concept'].lower()
            study = col_info['study_prefix']

            if concept not in concept_map:
                concept_map[concept] = {'columns': {}, 'studies': []}

            concept_map[concept]['columns'][col_name] = {
                'study': study,
                'coverage': col_info['coverage'],
                'dtype_inferred': col_info['dtype_inferred'],
                'n_unique': col_info['n_unique'],
                'sample_values': col_info['sample_values'][:5]
            }
            if study not in concept_map[concept]['studies']:
                concept_map[concept]['studies'].append(study)

    # Sort by number of studies (most shared concepts first)
    return dict(sorted(concept_map.items(), key=lambda x: -len(x[1]['studies'])))


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

    # Analyze each column
    print("Analyzing columns...")
    inventory = {}
    for col in df.columns:
        inventory[col] = analyze_column(df, col)

    # Group by concept
    print("Grouping by concept...")
    concept_groups = group_by_concept(inventory)

    # Summary statistics
    summary = {
        'n_donors': int(len(df)),
        'n_columns': int(len(df.columns)),
        'studies': ['ihbca', 'gray', 'murrow', 'pal', 'twigger', 'nee', 'kumar'],
        'columns_by_study': {},
        'shared_concepts': list(concept_groups.keys())[:20]  # top 20
    }

    # Count columns by study
    for study in summary['studies']:
        summary['columns_by_study'][study] = len([c for c in df.columns if c.startswith(f'{study}_')])

    # Build output
    output = {
        'summary': summary,
        'columns': inventory,
        'concept_groups': concept_groups
    }

    # Write YAML
    output_path = output_dir / 'column_inventory.yaml'
    print(f"Writing {output_path}")

    with open(output_path, 'w') as f:
        yaml.dump(output, f, default_flow_style=False, sort_keys=False, allow_unicode=True)

    print(f"Done. Inventory saved to {output_path}")

    # Print summary
    print("\n=== SUMMARY ===")
    print(f"Total columns: {len(df.columns)}")
    for study, count in summary['columns_by_study'].items():
        print(f"  {study}: {count} columns")
    print(f"\nShared concepts across studies: {len(concept_groups)}")
    for concept, info in list(concept_groups.items())[:10]:
        print(f"  {concept}: {len(info['studies'])} studies ({', '.join(info['studies'])})")


if __name__ == '__main__':
    main()
