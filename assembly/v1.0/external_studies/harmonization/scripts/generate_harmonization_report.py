#!/usr/bin/env python3
"""
Generate HTML Report for Donor Metadata Harmonization (Phase A.2)

Creates a comprehensive visual summary showing:
- Per-study original column values
- Mapping logic applied
- Harmonized output distributions
- Interactive distribution plots (Chart.js - client-side rendering)
- Cross-study comparison

Output: harmonization_report.html

# Documentation: because future you will have no idea what past you was thinking.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime
from collections import Counter
import html as html_lib
import json


COLORS = {
    'young': '#4ecca3',
    'old': '#e94560',
    'nulliparous': '#4ecca3',
    'parous': '#f9a825',
    'AR': '#4ecca3',
    'HR': '#e94560',
    'pre': '#4ecca3',
    'post': '#f9a825',
}

STUDY_COLORS = {
    'gray': '#e94560',
    'murrow': '#f9a825',
    'nee': '#4ecca3',
    'twigger': '#00d9ff',
    'pal': '#a855f7',
    'kumar': '#ff6b6b',
    'reed': '#74b9ff',
}


def generate_chart_data(df):
    """Generate all chart data as JSON for Chart.js."""
    data = {}

    # Donor counts by study
    counts = df['study'].value_counts().sort_values(ascending=True)
    data['donorCounts'] = {
        'labels': list(counts.index),
        'values': [int(v) for v in counts.values],
        'colors': [STUDY_COLORS.get(s, '#888') for s in counts.index]
    }

    # Age continuous - histogram data per study
    age_data = {}
    for study in sorted(df['study'].unique()):
        study_ages = df[df['study'] == study]['age_continuous'].dropna().tolist()
        if study_ages:
            age_data[study] = [float(a) for a in study_ages]
    data['ageContinuous'] = {
        'studies': age_data,
        'colors': {s: STUDY_COLORS.get(s, '#888') for s in age_data.keys()}
    }

    # Age binary by study
    binary_data = df[df['age_binary'].notna()].groupby(['study', 'age_binary']).size().unstack(fill_value=0)
    if not binary_data.empty:
        data['ageBinary'] = {
            'labels': list(binary_data.index),
            'young': [int(binary_data.loc[s, 'young']) if 'young' in binary_data.columns else 0 for s in binary_data.index],
            'old': [int(binary_data.loc[s, 'old']) if 'old' in binary_data.columns else 0 for s in binary_data.index]
        }

    # Parity by study
    parity_data = df[df['parity_binary'].notna()].groupby(['study', 'parity_binary']).size().unstack(fill_value=0)
    if not parity_data.empty:
        data['parity'] = {
            'labels': list(parity_data.index),
            'nulliparous': [int(parity_data.loc[s, 'nulliparous']) if 'nulliparous' in parity_data.columns else 0 for s in parity_data.index],
            'parous': [int(parity_data.loc[s, 'parous']) if 'parous' in parity_data.columns else 0 for s in parity_data.index]
        }

    # Risk status by study
    risk_data = df[df['risk_status_binary'].notna()].groupby(['study', 'risk_status_binary']).size().unstack(fill_value=0)
    if not risk_data.empty:
        data['risk'] = {
            'labels': list(risk_data.index),
            'AR': [int(risk_data.loc[s, 'AR']) if 'AR' in risk_data.columns else 0 for s in risk_data.index],
            'HR': [int(risk_data.loc[s, 'HR']) if 'HR' in risk_data.columns else 0 for s in risk_data.index]
        }

    # Menopausal status by study
    meno_data = df[df['menopausal_status_binary'].notna()].groupby(['study', 'menopausal_status_binary']).size().unstack(fill_value=0)
    if not meno_data.empty:
        data['menopausal'] = {
            'labels': list(meno_data.index),
            'pre': [int(meno_data.loc[s, 'pre']) if 'pre' in meno_data.columns else 0 for s in meno_data.index],
            'post': [int(meno_data.loc[s, 'post']) if 'post' in meno_data.columns else 0 for s in meno_data.index]
        }

    # Coverage matrix - expanded for all key columns
    conditions = [
        'age_binary', 'parity_binary', 'risk_status_binary', 'menopausal_status_binary',  # Legacy
        'brca_genotype', 'cancer_history', 'tissue_indication', 'risk_genotype_only',  # Orthogonal risk
        'ethnicity_grouped', 'bmi_category',  # Demographics
    ]
    studies = sorted(df['study'].unique())
    coverage_matrix = []
    for study in studies:
        study_df = df[df['study'] == study]
        row = []
        for cond in conditions:
            if cond in df.columns:
                # For orthogonal dimensions, count non-'unknown' values as coverage
                if cond in ['brca_genotype', 'cancer_history', 'tissue_indication']:
                    coverage = ((study_df[cond].notna()) & (study_df[cond] != 'unknown')).mean() * 100
                else:
                    coverage = study_df[cond].notna().mean() * 100
            else:
                coverage = 0.0
            row.append(round(coverage, 1))
        coverage_matrix.append(row)
    data['coverage'] = {
        'studies': studies,
        'conditions': [c.replace('_', ' ').title() for c in conditions],
        'matrix': coverage_matrix
    }

    return json.dumps(data)


def escape(text):
    """HTML escape text."""
    if pd.isna(text):
        return '<span class="na">NA</span>'
    return html_lib.escape(str(text))


def value_distribution_table(series, max_rows=15):
    """Create HTML table of value distribution."""
    counts = series.value_counts(dropna=False)
    total = len(series)

    rows = []
    for i, (val, count) in enumerate(counts.items()):
        if i >= max_rows:
            remaining = len(counts) - max_rows
            rows.append(f'<tr class="truncated"><td colspan="3">... and {remaining} more values</td></tr>')
            break
        pct = 100 * count / total
        val_display = '<span class="na">NA/Missing</span>' if pd.isna(val) else escape(val)
        rows.append(f'<tr><td>{val_display}</td><td class="num">{count:,}</td><td class="num">{pct:.1f}%</td></tr>')

    return f'''
    <table class="dist-table">
        <thead><tr><th>Value</th><th>Count</th><th>%</th></tr></thead>
        <tbody>{''.join(rows)}</tbody>
    </table>
    '''


def mapping_box(source_col, transform, warning=None):
    """Create mapping logic display box."""
    warning_html = f'<div class="warning">⚠️ {escape(warning)}</div>' if warning else ''
    return f'''
    <div class="mapping-box">
        <div class="source-col">Source: <code>{escape(source_col)}</code></div>
        <div class="transform">→ {escape(transform)}</div>
        {warning_html}
    </div>
    '''


def generate_study_section(study, unified_df, harmonized_df, mappings):
    """Generate HTML section for a single study."""
    study_unified = unified_df[unified_df['ihbca_dataset'] == study]
    study_harmonized = harmonized_df[harmonized_df['study'] == study]
    n_donors = len(study_unified)

    # Conditions config - organized by category
    condition_categories = {
        'Age': [
            {
                'name': 'age_continuous',
                'title': 'Age (Continuous)',
                'source_cols': {
                    'gray': 'gray_Age', 'murrow': 'murrow_Age', 'nee': 'nee_Age',
                    'twigger': 'twigger_Maternal_age', 'pal': 'pal_patient_age',
                    'kumar': 'kumar_Age', 'reed': 'ihbca_age'
                },
                'transforms': {
                    'gray': 'cast to float', 'murrow': 'cast to float', 'nee': 'cast to float',
                    'twigger': 'cast to float', 'pal': 'cast to float',
                    'kumar': 'EXCLUDED (binary only: Y/O)', 'reed': 'cast to float (from ihbca_age)'
                }
            },
            {
                'name': 'age_binary',
                'title': 'Age (Binary: young/old)',
                'source_cols': {
                    'gray': 'gray_Age', 'murrow': 'murrow_Age', 'nee': 'nee_Age',
                    'twigger': 'twigger_Maternal_age', 'pal': 'pal_patient_age',
                    'kumar': 'kumar_Age', 'reed': 'ihbca_age'
                },
                'transforms': {
                    'gray': 'age < 50 → young, else → old', 'murrow': 'age < 50 → young, else → old',
                    'nee': 'age < 50 → young, else → old', 'twigger': 'age < 50 → young, else → old',
                    'pal': 'age < 50 → young, else → old', 'kumar': 'Y → young, O → old',
                    'reed': 'age < 50 → young, else → old (from ihbca_age)'
                },
                'warnings': {'murrow': 'All donors <40, no contrast'}
            },
        ],
        'Parity': [
            {
                'name': 'parity_count',
                'title': 'Parity (Numeric Count)',
                'source_cols': {
                    'gray': 'gray_Births', 'murrow': 'murrow_Parity', 'nee': 'nee_Parity',
                    'twigger': 'twigger_Parity', 'pal': 'pal_parity',
                    'kumar': 'kumar_Parity', 'reed': 'ihbca_parity'
                },
                'transforms': {
                    'gray': 'cast to int', 'murrow': 'cast to int (unknown → NA)',
                    'nee': 'Parse GXPX → P value', 'twigger': 'cast to int',
                    'pal': 'binary only → NA', 'kumar': 'binary only → NA',
                    'reed': 'cast to int'
                }
            },
            {
                'name': 'parity_binary',
                'title': 'Parity (Binary: nulliparous/parous)',
                'source_cols': {
                    'gray': 'gray_Births', 'murrow': 'murrow_Parity', 'nee': 'nee_Parity',
                    'twigger': 'twigger_Parity', 'pal': 'pal_parity',
                    'kumar': 'kumar_Parity', 'reed': 'ihbca_parous'
                },
                'transforms': {
                    'gray': '0 → nulliparous, ≥1 → parous',
                    'murrow': '0 → nulliparous, 1/2/3 → parous, unknown → NA',
                    'nee': 'Parse GXPX: P0 → nulliparous, P≥1 → parous',
                    'twigger': '0 → nulliparous, ≥1 → parous',
                    'pal': 'Nulliparous/Parous directly',
                    'kumar': '0 → nulliparous, 1 → parous, unknown → NA',
                    'reed': 'True/1 → parous, False/0 → nulliparous'
                },
                'warnings': {'pal': 'CONFOUNDED with menopausal status - skip DA'}
            },
            {
                'name': 'age_at_first_birth',
                'title': 'Age at First Birth',
                'source_cols': {
                    'gray': None, 'murrow': None, 'nee': 'nee_Age at First Pregnancy',
                    'twigger': None, 'pal': None,
                    'kumar': None, 'reed': 'reed_Age_at_first_live_birth'
                },
                'transforms': {
                    'gray': 'NOT COLLECTED', 'murrow': 'NOT COLLECTED',
                    'nee': 'cast to int', 'twigger': 'NOT COLLECTED',
                    'pal': 'NOT COLLECTED', 'kumar': 'NOT COLLECTED',
                    'reed': 'cast to int'
                }
            },
        ],
        'Risk Dimensions (Orthogonal)': [
            {
                'name': 'brca_genotype',
                'title': 'BRCA Genotype',
                'source_cols': {
                    'gray': 'ihbca_risk_status', 'murrow': None,
                    'nee': 'nee_BRCA1 Mutation', 'twigger': None,
                    'pal': 'pal_brca1_mutation', 'kumar': None, 'reed': 'reed_brca_status'
                },
                'transforms': {
                    'gray': 'HR-BR1 → BRCA1, HR-BR2 → BRCA2, HR-RAD51C → RAD51C, AR → negative',
                    'murrow': 'All unknown (not tested)',
                    'nee': 'Confirmed carrier → BRCA1, No (confirmed) → negative',
                    'twigger': 'All unknown (not tested)',
                    'pal': 'Yes → BRCA1, No/Wild type → negative',
                    'kumar': 'All unknown (not tested)',
                    'reed': 'BRCA1/BRCA2/None/Wild-type mapping'
                }
            },
            {
                'name': 'cancer_history',
                'title': 'Cancer History',
                'source_cols': {
                    'gray': 'ihbca_risk_status', 'murrow': None,
                    'nee': 'nee_Cancer History', 'twigger': None,
                    'pal': None, 'kumar': 'kumar_Tissue_Source', 'reed': 'ihbca_risk_status'
                },
                'transforms': {
                    'gray': 'HR-*-c suffix → yes, else → no',
                    'murrow': 'All no (reduction mammoplasty)',
                    'nee': 'Yes/No/Unknown directly',
                    'twigger': 'All no (lactating donors)',
                    'pal': 'All no (prophylactic/reduction)',
                    'kumar': 'Cancer Mastectomy → yes (contralateral), else → no',
                    'reed': 'HR-*-c suffix → yes, else → no'
                }
            },
            {
                'name': 'tissue_indication',
                'title': 'Tissue Indication',
                'source_cols': {
                    'gray': 'gray_Surgery', 'murrow': None,
                    'nee': 'nee_Tissue Source', 'twigger': None,
                    'pal': 'ihbca_risk_status', 'kumar': 'kumar_Tissue_Source', 'reed': 'reed_surgery_type'
                },
                'transforms': {
                    'gray': 'Reduction/Prophylactic/Mastectomy mapping',
                    'murrow': 'All reduction (mammoplasty)',
                    'nee': 'Reduction/Prophylactic/Mastectomy mapping',
                    'twigger': 'All unknown',
                    'pal': 'AR → reduction, HR-* → prophylactic',
                    'kumar': 'Reduction/Prophylactic/Cancer Mastectomy → contralateral',
                    'reed': 'Reduction/Prophylactic/Mastectomy mapping'
                }
            },
        ],
        'Risk Composite': [
            {
                'name': 'risk_status_binary',
                'title': 'Risk Status (Binary: AR/HR) - DERIVED',
                'source_cols': {
                    'gray': 'DERIVED', 'murrow': 'DERIVED',
                    'nee': 'DERIVED', 'twigger': 'DERIVED',
                    'pal': 'DERIVED', 'kumar': 'DERIVED', 'reed': 'DERIVED'
                },
                'transforms': {
                    'gray': 'HR if: BRCA+ OR cancer_history=yes OR tissue≠reduction',
                    'murrow': 'AR (all reduction, no BRCA, no cancer)',
                    'nee': 'HR if: BRCA+ OR cancer_history=yes OR tissue≠reduction',
                    'twigger': 'AR (lactating, general population)',
                    'pal': 'HR if: BRCA+ OR tissue=prophylactic',
                    'kumar': 'HR if: cancer_history=yes OR tissue≠reduction',
                    'reed': 'HR if: BRCA+ OR cancer_history=yes OR tissue≠reduction'
                },
                'warnings': {
                    'murrow': 'AR-only study - no contrast',
                    'twigger': 'AR-only study - no contrast'
                }
            },
            {
                'name': 'risk_genotype_only',
                'title': 'Risk (Genotype Only: strict BRCA)',
                'source_cols': {
                    'gray': 'DERIVED', 'murrow': 'DERIVED',
                    'nee': 'DERIVED', 'twigger': 'DERIVED',
                    'pal': 'DERIVED', 'kumar': 'DERIVED', 'reed': 'DERIVED'
                },
                'transforms': {
                    'gray': 'BRCA1/BRCA2/RAD51C → HR, negative → AR, else → None',
                    'murrow': 'All None (unknown genotype)',
                    'nee': 'BRCA1 → HR, negative → AR',
                    'twigger': 'All None (unknown genotype)',
                    'pal': 'BRCA1 → HR, negative → AR',
                    'kumar': 'All None (unknown genotype)',
                    'reed': 'BRCA1/BRCA2 → HR, negative → AR'
                },
                'warnings': {
                    'murrow': 'All None - not testable',
                    'twigger': 'All None - not testable',
                    'kumar': 'All None - not testable'
                }
            },
        ],
        'Menopausal': [
            {
                'name': 'menopausal_status_detailed',
                'title': 'Menopausal Status (Detailed)',
                'source_cols': {
                    'gray': 'gray_Menopause', 'murrow': None,
                    'nee': 'nee_Menstrual Status', 'twigger': None,
                    'pal': 'pal_menopausal_status', 'kumar': 'kumar_Menopause',
                    'reed': 'reed_Menopause_status'
                },
                'transforms': {
                    'gray': 'Premenopausal → pre, Surgical menopause → post_surgical',
                    'murrow': 'Inferred from age → pre',
                    'nee': 'Premenopause → pre, Menopausal → post_natural',
                    'twigger': 'NOT COLLECTED',
                    'pal': 'pre/peri/post mapping',
                    'kumar': 'pre/post/unknown mapping',
                    'reed': 'Pre/Peri/Post/Post_(surgically_induced) → pre/peri/post_natural/post_surgical'
                }
            },
            {
                'name': 'menopausal_status_binary',
                'title': 'Menopausal Status (Binary: pre/post)',
                'source_cols': {
                    'gray': 'DERIVED', 'murrow': 'DERIVED',
                    'nee': 'DERIVED', 'twigger': 'DERIVED',
                    'pal': 'DERIVED', 'kumar': 'DERIVED', 'reed': 'DERIVED'
                },
                'transforms': {
                    'gray': 'pre → pre, post_surgical → post',
                    'murrow': 'pre → pre (no contrast)',
                    'nee': 'pre → pre, post_natural → post',
                    'twigger': 'NOT COLLECTED',
                    'pal': 'pre/peri → pre, post → post',
                    'kumar': 'pre → pre, post → post',
                    'reed': 'pre/peri → pre, post_natural/post_surgical → post'
                },
                'warnings': {'gray': 'Includes surgical menopause'}
            },
        ],
        'Demographics': [
            {
                'name': 'ethnicity_verbatim',
                'title': 'Ethnicity (Verbatim)',
                'source_cols': {
                    'gray': 'gray_Race', 'murrow': None,
                    'nee': 'nee_Ethnicity', 'twigger': None,
                    'pal': None, 'kumar': 'kumar_Race', 'reed': 'reed_Race'
                },
                'transforms': {
                    'gray': 'Direct copy', 'murrow': 'NOT COLLECTED',
                    'nee': 'Direct copy', 'twigger': 'NOT COLLECTED',
                    'pal': 'NOT COLLECTED', 'kumar': 'Direct copy', 'reed': 'Direct copy'
                }
            },
            {
                'name': 'ethnicity_grouped',
                'title': 'Ethnicity (Grouped)',
                'source_cols': {
                    'gray': 'DERIVED', 'murrow': 'DERIVED',
                    'nee': 'DERIVED', 'twigger': 'DERIVED',
                    'pal': 'DERIVED', 'kumar': 'DERIVED', 'reed': 'DERIVED'
                },
                'transforms': {
                    'gray': 'White/Black/Jewish → white/black/jewish',
                    'murrow': 'NOT COLLECTED',
                    'nee': 'White/Black/Hispanic → white/black/hispanic',
                    'twigger': 'NOT COLLECTED',
                    'pal': 'NOT COLLECTED',
                    'kumar': 'Caucasian/African/Asian/Filipino → white/black/asian',
                    'reed': 'White/Black/Hispanic/Asian/Mixed → grouped'
                }
            },
            {
                'name': 'bmi_continuous',
                'title': 'BMI (Continuous)',
                'source_cols': {
                    'gray': None, 'murrow': 'murrow_BMI',
                    'nee': 'nee_BMI', 'twigger': None,
                    'pal': None, 'kumar': 'kumar_BMI', 'reed': 'reed_BMI'
                },
                'transforms': {
                    'gray': 'NOT COLLECTED', 'murrow': 'cast to float',
                    'nee': 'cast to float', 'twigger': 'NOT COLLECTED',
                    'pal': 'NOT COLLECTED', 'kumar': 'cast to float', 'reed': 'cast to float'
                }
            },
            {
                'name': 'bmi_category',
                'title': 'BMI (Category)',
                'source_cols': {
                    'gray': 'DERIVED', 'murrow': 'DERIVED',
                    'nee': 'DERIVED', 'twigger': 'DERIVED',
                    'pal': 'DERIVED', 'kumar': 'DERIVED', 'reed': 'DERIVED'
                },
                'transforms': {
                    'gray': 'NOT COLLECTED',
                    'murrow': '<18.5 underweight, 18.5-24.9 normal, 25-29.9 overweight, ≥30 obese',
                    'nee': '<18.5 underweight, 18.5-24.9 normal, 25-29.9 overweight, ≥30 obese',
                    'twigger': 'NOT COLLECTED',
                    'pal': 'NOT COLLECTED',
                    'kumar': '<18.5 underweight, 18.5-24.9 normal, 25-29.9 overweight, ≥30 obese',
                    'reed': '<18.5 underweight, 18.5-24.9 normal, 25-29.9 overweight, ≥30 obese'
                }
            },
        ],
        'Sample Metadata': [
            {
                'name': 'sample_preservation',
                'title': 'Sample Preservation',
                'source_cols': {
                    'gray': None, 'murrow': None,
                    'nee': None, 'twigger': 'twigger_Tissue_state',
                    'pal': None, 'kumar': None, 'reed': None
                },
                'transforms': {
                    'gray': 'NOT COLLECTED', 'murrow': 'NOT COLLECTED',
                    'nee': 'NOT COLLECTED', 'twigger': 'Fresh/Frozen → fresh/frozen',
                    'pal': 'NOT COLLECTED', 'kumar': 'NOT COLLECTED', 'reed': 'NOT COLLECTED'
                }
            },
            {
                'name': 'sample_type',
                'title': 'Sample Type',
                'source_cols': {
                    'gray': 'ihbca_sample_type', 'murrow': 'ihbca_sample_type',
                    'nee': 'ihbca_sample_type', 'twigger': 'ihbca_sample_type',
                    'pal': 'ihbca_sample_type', 'kumar': 'ihbca_sample_type', 'reed': 'ihbca_sample_type'
                },
                'transforms': {
                    'gray': 'mixed', 'murrow': 'mixed',
                    'nee': 'mixed', 'twigger': 'mixed',
                    'pal': 'mixed', 'kumar': 'mixed',
                    'reed': 'mixed/supernatant/organoid/organoid_lp'
                }
            },
            {
                'name': 'facs_status',
                'title': 'FACS Status',
                'source_cols': {
                    'gray': 'ihbca_FACS_status', 'murrow': 'ihbca_FACS_status',
                    'nee': 'ihbca_FACS_status', 'twigger': 'ihbca_FACS_status',
                    'pal': 'ihbca_FACS_status', 'kumar': 'ihbca_FACS_status', 'reed': 'ihbca_FACS_status'
                },
                'transforms': {
                    'gray': 'not_sorted', 'murrow': 'live_sorted/cell_type_sorted',
                    'nee': 'cell_type_sorted', 'twigger': 'no_sort',
                    'pal': 'not_sorted/cell_type_sorted', 'kumar': 'live_sorted',
                    'reed': 'not_sorted/cell_type_sorted'
                }
            },
            {
                'name': 'dissociation_minutes',
                'title': 'Dissociation Time (minutes)',
                'source_cols': {
                    'gray': None, 'murrow': None,
                    'nee': None, 'twigger': None,
                    'pal': None, 'kumar': None, 'reed': 'reed_dissociation_minutes'
                },
                'transforms': {
                    'gray': 'NOT COLLECTED', 'murrow': 'NOT COLLECTED',
                    'nee': 'NOT COLLECTED', 'twigger': 'NOT COLLECTED',
                    'pal': 'NOT COLLECTED', 'kumar': 'NOT COLLECTED',
                    'reed': '5/7/10 min protocols'
                }
            },
        ],
    }

    condition_sections = []
    for category, conditions in condition_categories.items():
        # Add category header
        condition_sections.append(f'<h4 class="category-divider">{category}</h4>')

        for cond in conditions:
            source_col = cond['source_cols'].get(study)
            transform = cond['transforms'].get(study, 'N/A')
            warning = cond.get('warnings', {}).get(study)

            # Handle DERIVED columns specially
            if source_col == 'DERIVED':
                orig_dist = '<div class="derived-box">DERIVED from other columns</div>'
                mapping_html = mapping_box('Multiple columns', transform, warning)
            elif source_col and source_col in study_unified.columns:
                orig_dist = value_distribution_table(study_unified[source_col])
                mapping_html = mapping_box(source_col, transform, warning)
            elif source_col:
                orig_dist = f'<div class="na-box">Column <code>{source_col}</code> not found</div>'
                mapping_html = mapping_box(source_col, transform, warning)
            else:
                orig_dist = '<div class="na-box">Not collected for this study</div>'
                mapping_html = '<div class="na-box">Not collected</div>'

            # Check if harmonized column exists
            if cond['name'] in study_harmonized.columns:
                harm_dist = value_distribution_table(study_harmonized[cond['name']])
            else:
                harm_dist = '<div class="na-box">Column not in output</div>'

            condition_sections.append(f'''
            <div class="condition-block">
                <h5>{cond['title']}</h5>
                <div class="condition-grid">
                    <div class="original">
                        <h6>Original Values</h6>
                        {orig_dist}
                    </div>
                    <div class="mapping">
                        <h6>Mapping Logic</h6>
                        {mapping_html}
                    </div>
                    <div class="harmonized">
                        <h6>Harmonized Output</h6>
                        {harm_dist}
                    </div>
                </div>
            </div>
            ''')

    return f'''
    <div class="study-section" id="study-{study}">
        <h3>{study.upper()} <span class="donor-count">({n_donors} donors)</span></h3>
        {''.join(condition_sections)}
    </div>
    '''


def generate_cross_study_summary(harmonized_df):
    """Generate cross-study comparison tables."""
    conditions = [
        # Legacy columns
        'age_continuous', 'age_binary', 'parity_binary', 'risk_status_binary', 'menopausal_status_binary',
        # Orthogonal risk dimensions
        'brca_genotype', 'cancer_history', 'tissue_indication', 'risk_genotype_only',
        # New detailed columns
        'menopausal_status_detailed', 'parity_count',
        # Demographics
        'ethnicity_grouped', 'bmi_category',
        # Sample metadata
        'sample_type', 'facs_status',
    ]
    # Filter to columns that exist
    conditions = [c for c in conditions if c in harmonized_df.columns]
    studies = harmonized_df['study'].unique()

    tables = []
    for cond in conditions:
        rows = []
        for study in sorted(studies):
            study_data = harmonized_df[harmonized_df['study'] == study]
            n_total = len(study_data)
            n_valid = study_data[cond].notna().sum()
            pct = 100 * n_valid / n_total if n_total > 0 else 0

            if cond == 'age_continuous':
                vals = study_data[cond].dropna()
                dist = f'range: {vals.min():.0f}-{vals.max():.0f}' if len(vals) > 0 else '-'
            else:
                dist_counts = study_data[cond].value_counts()
                dist = ', '.join(f'{k}={v}' for k, v in dist_counts.items())

            rows.append(f'<tr><td>{study}</td><td class="num">{n_valid}/{n_total}</td><td class="num">{pct:.1f}%</td><td>{dist}</td></tr>')

        tables.append(f'''
        <div class="summary-table-block">
            <h4>{cond}</h4>
            <table class="summary-table">
                <thead><tr><th>Study</th><th>Coverage</th><th>%</th><th>Distribution</th></tr></thead>
                <tbody>{''.join(rows)}</tbody>
            </table>
        </div>
        ''')

    return f'''
    <div class="cross-study-summary">
        <h2>Cross-Study Summary</h2>
        {''.join(tables)}
    </div>
    '''


def generate_testability_matrix(harmonized_df):
    """Generate testability matrix showing which conditions can be tested where.

    Updated to include all 22 harmonized columns grouped by category.
    """
    # Define matrix by category for better organization
    matrix_categories = {
        'Age': {
            'age_continuous': {
                'gray': ('✓', None), 'murrow': ('✓', 'no contrast (<50)'), 'nee': ('✓', None),
                'twigger': ('✓', None), 'pal': ('✓', None),
                'kumar': ('✗', 'binary only'), 'reed': ('✓', None)
            },
            'age_binary': {
                'gray': ('✓', None), 'murrow': ('✗', 'no contrast'), 'nee': ('✓', None),
                'twigger': ('✓', None), 'pal': ('✓', None),
                'kumar': ('✓', None), 'reed': ('✓', None)
            },
        },
        'Parity': {
            'parity_count': {
                'gray': ('✓', None), 'murrow': ('✓', '32% coverage'), 'nee': ('✓', '64% coverage'),
                'twigger': ('✓', None), 'pal': ('✗', 'binary only'),
                'kumar': ('✗', 'binary only'), 'reed': ('✓', None)
            },
            'parity_binary': {
                'gray': ('✓', None), 'murrow': ('✓', None), 'nee': ('✓', None),
                'twigger': ('✓', None), 'pal': ('✗', 'confounded'),
                'kumar': ('✓', None), 'reed': ('✓', None)
            },
            'age_at_first_birth': {
                'gray': ('✗', 'not collected'), 'murrow': ('✗', 'not collected'), 'nee': ('✓', '27% coverage'),
                'twigger': ('✗', 'not collected'), 'pal': ('✗', 'not collected'),
                'kumar': ('✗', 'not collected'), 'reed': ('✓', '35% coverage')
            },
        },
        'Risk Dimensions (Orthogonal)': {
            'brca_genotype': {
                'gray': ('✓', 'BR1/BR2/RAD51C/neg'), 'murrow': ('✗', 'all unknown'), 'nee': ('✓', 'BR1/neg'),
                'twigger': ('✗', 'all unknown'), 'pal': ('✓', 'BR1/neg'),
                'kumar': ('✗', 'all unknown'), 'reed': ('✓', 'BR1/BR2/neg')
            },
            'cancer_history': {
                'gray': ('✓', 'yes/no'), 'murrow': ('✗', 'all no'), 'nee': ('✓', 'yes/no/unk'),
                'twigger': ('✗', 'all no'), 'pal': ('✗', 'all no'),
                'kumar': ('✓', 'yes/no'), 'reed': ('✓', 'yes/no')
            },
            'tissue_indication': {
                'gray': ('✓', 'red/proph/contra'), 'murrow': ('✗', 'all reduction'), 'nee': ('✓', 'red/proph/contra'),
                'twigger': ('✗', 'all unknown'), 'pal': ('✓', 'red/proph'),
                'kumar': ('✓', 'red/proph/contra'), 'reed': ('✓', 'red/proph/contra')
            },
        },
        'Risk Composite': {
            'risk_status_binary': {
                'gray': ('✓', None), 'murrow': ('✗', 'no contrast'), 'nee': ('✓', None),
                'twigger': ('✗', 'no contrast'), 'pal': ('✓', 'small HR n=4'),
                'kumar': ('✓', None), 'reed': ('✓', None)
            },
            'risk_genotype_only': {
                'gray': ('✓', 'strict BRCA'), 'murrow': ('✗', 'all None'), 'nee': ('✓', 'strict BRCA'),
                'twigger': ('✗', 'all None'), 'pal': ('✓', 'strict BRCA'),
                'kumar': ('✗', 'all None'), 'reed': ('✓', 'strict BRCA')
            },
        },
        'Menopausal': {
            'menopausal_status_detailed': {
                'gray': ('✓', 'pre/post_surg'), 'murrow': ('✗', 'all pre'), 'nee': ('✓', 'pre/post'),
                'twigger': ('✗', 'not collected'), 'pal': ('✓', 'pre/peri/post'),
                'kumar': ('✓', 'pre/post/unk'), 'reed': ('✓', 'pre/peri/post_nat/surg')
            },
            'menopausal_status_binary': {
                'gray': ('✓', 'includes surgical'), 'murrow': ('✗', 'no contrast'), 'nee': ('✓', None),
                'twigger': ('✗', 'not collected'), 'pal': ('✓', None),
                'kumar': ('✓', None), 'reed': ('✓', '67% coverage')
            },
        },
        'Demographics': {
            'ethnicity_verbatim': {
                'gray': ('✓', '4 groups'), 'murrow': ('✗', 'not collected'), 'nee': ('✓', '5 groups'),
                'twigger': ('✗', 'not collected'), 'pal': ('✗', 'not collected'),
                'kumar': ('✓', '4 groups'), 'reed': ('✓', '6 groups')
            },
            'ethnicity_grouped': {
                'gray': ('✓', 'white/black/jewish'), 'murrow': ('✗', 'not collected'), 'nee': ('✓', 'white/black/hisp'),
                'twigger': ('✗', 'not collected'), 'pal': ('✗', 'not collected'),
                'kumar': ('✓', 'white/black/asian'), 'reed': ('✓', '5 groups')
            },
            'bmi_continuous': {
                'gray': ('✗', 'not collected'), 'murrow': ('✓', '54% coverage'), 'nee': ('✓', '68% coverage'),
                'twigger': ('✗', 'not collected'), 'pal': ('✗', 'not collected'),
                'kumar': ('✓', '94% coverage'), 'reed': ('✓', '71% coverage')
            },
            'bmi_category': {
                'gray': ('✗', 'not collected'), 'murrow': ('✓', 'normal/over/obese'), 'nee': ('✓', 'normal/over/obese'),
                'twigger': ('✗', 'not collected'), 'pal': ('✗', 'not collected'),
                'kumar': ('✓', 'all 4 categories'), 'reed': ('✓', 'all 4 categories')
            },
        },
        'Sample Metadata': {
            'sample_preservation': {
                'gray': ('✗', 'not collected'), 'murrow': ('✗', 'not collected'), 'nee': ('✗', 'not collected'),
                'twigger': ('✓', 'fresh/frozen'), 'pal': ('✗', 'not collected'),
                'kumar': ('✗', 'not collected'), 'reed': ('✗', 'not collected')
            },
            'sample_type': {
                'gray': ('✗', 'all mixed'), 'murrow': ('✗', 'all mixed'), 'nee': ('✗', 'all mixed'),
                'twigger': ('✗', 'all mixed'), 'pal': ('✗', 'all mixed'),
                'kumar': ('✗', 'all mixed'), 'reed': ('✓', 'mixed/sup/org/org_lp')
            },
            'facs_status': {
                'gray': ('✗', 'all not_sorted'), 'murrow': ('✓', 'live/cell_type'), 'nee': ('✗', 'all cell_type'),
                'twigger': ('✗', 'all no_sort'), 'pal': ('✓', 'not/cell_type'),
                'kumar': ('✗', 'all live_sorted'), 'reed': ('✓', 'not/cell_type')
            },
            'dissociation_minutes': {
                'gray': ('✗', 'not collected'), 'murrow': ('✗', 'not collected'), 'nee': ('✗', 'not collected'),
                'twigger': ('✗', 'not collected'), 'pal': ('✗', 'not collected'),
                'kumar': ('✗', 'not collected'), 'reed': ('✓', '5/7/10 min')
            },
        },
    }

    studies = ['gray', 'murrow', 'nee', 'twigger', 'pal', 'kumar', 'reed']

    # Build HTML with category headers
    header = '<tr><th>Condition</th>' + ''.join(f'<th class="study-header">{s}</th>' for s in studies) + '</tr>'

    rows = []
    for category, conditions in matrix_categories.items():
        # Category header row
        rows.append(f'<tr class="category-row"><td colspan="{len(studies)+1}" class="category-header">{category}</td></tr>')

        for cond, study_data in conditions.items():
            cells = []
            for study in studies:
                status, note = study_data.get(study, ('?', None))
                cell_class = 'testable' if status == '✓' else 'not-testable' if status == '✗' else 'unknown'
                note_html = f'<span class="note" title="{note}">*</span>' if note else ''
                cells.append(f'<td class="{cell_class}">{status}{note_html}</td>')
            rows.append(f'<tr><td class="cond-name">{cond}</td>{"".join(cells)}</tr>')

    return f'''
    <div class="testability-matrix">
        <h2>DA Testability Matrix</h2>
        <p class="matrix-note">✓ = testable (has contrast), ✗ = not testable, * = see note (hover)</p>
        <table class="matrix-table">
            <thead>{header}</thead>
            <tbody>{''.join(rows)}</tbody>
        </table>
    </div>
    '''


def generate_html_report(unified_df, harmonized_df):
    """Generate complete HTML report with interactive Chart.js plots."""
    studies = ['gray', 'murrow', 'nee', 'twigger', 'pal', 'kumar', 'reed']

    # Generate chart data as JSON
    print("  Generating chart data...")
    chart_data = generate_chart_data(harmonized_df)

    nav_items = ''.join(f'<a href="#study-{s}">{s.upper()}</a>' for s in studies)
    study_sections = ''.join(generate_study_section(s, unified_df, harmonized_df, {}) for s in studies)
    cross_study = generate_cross_study_summary(harmonized_df)
    testability = generate_testability_matrix(harmonized_df)

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    html_content = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Phase A.2: Donor Metadata Harmonization Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
    <style>
        :root {{
            --bg: #1a1a2e;
            --card-bg: #16213e;
            --text: #e8e8e8;
            --text-muted: #a0a0a0;
            --accent: #0f3460;
            --success: #4ecca3;
            --warning: #f9a825;
            --error: #e94560;
            --border: #2a2a4a;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: var(--bg);
            color: var(--text);
            line-height: 1.6;
            padding: 2rem;
        }}
        h1 {{ color: var(--success); margin-bottom: 0.5rem; }}
        h2 {{ color: var(--text); margin: 2rem 0 1rem; border-bottom: 2px solid var(--accent); padding-bottom: 0.5rem; }}
        h3 {{ color: var(--success); margin-bottom: 1rem; }}
        h4 {{ color: var(--text); margin-bottom: 0.5rem; font-size: 1rem; }}
        h5 {{ color: var(--text-muted); font-size: 0.85rem; margin-bottom: 0.5rem; }}
        .header {{ margin-bottom: 2rem; }}
        .timestamp {{ color: var(--text-muted); font-size: 0.9rem; }}
        .donor-count {{ color: var(--text-muted); font-weight: normal; font-size: 0.9rem; }}
        nav {{
            background: var(--card-bg);
            padding: 1rem;
            border-radius: 8px;
            margin-bottom: 2rem;
            position: sticky;
            top: 0;
            z-index: 100;
        }}
        nav a {{
            color: var(--text);
            text-decoration: none;
            padding: 0.5rem 1rem;
            margin-right: 0.5rem;
            border-radius: 4px;
            background: var(--accent);
            display: inline-block;
        }}
        nav a:hover {{ background: var(--success); color: var(--bg); }}
        .overview-section {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }}
        .plot-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 1rem;
            margin-top: 1rem;
        }}
        .plot-container {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 0.5rem;
        }}
        .plot-full {{
            margin-top: 1rem;
        }}
        .study-section {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }}
        .condition-block {{
            background: var(--bg);
            border-radius: 6px;
            padding: 1rem;
            margin-bottom: 1rem;
        }}
        .condition-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 1rem;
        }}
        .original, .mapping, .harmonized {{
            background: var(--card-bg);
            padding: 1rem;
            border-radius: 6px;
        }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.85rem; }}
        th, td {{ padding: 0.4rem 0.6rem; text-align: left; border-bottom: 1px solid var(--border); }}
        th {{ color: var(--text-muted); font-weight: 600; }}
        .num {{ text-align: right; font-family: monospace; }}
        .dist-table {{ margin-top: 0.5rem; }}
        .truncated {{ color: var(--text-muted); font-style: italic; }}
        .na {{ color: var(--text-muted); font-style: italic; }}
        .na-box {{ color: var(--text-muted); font-style: italic; padding: 1rem; text-align: center; }}
        .mapping-box {{
            background: var(--bg);
            padding: 0.75rem;
            border-radius: 4px;
            font-size: 0.9rem;
        }}
        .source-col {{ margin-bottom: 0.5rem; }}
        .source-col code {{ background: var(--accent); padding: 0.2rem 0.4rem; border-radius: 3px; }}
        .transform {{ color: var(--success); }}
        .warning {{ color: var(--warning); margin-top: 0.5rem; font-size: 0.85rem; }}
        .cross-study-summary, .testability-matrix {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 1.5rem;
            margin-bottom: 2rem;
        }}
        .summary-table-block {{ margin-bottom: 1.5rem; }}
        .summary-table {{ margin-top: 0.5rem; }}
        .matrix-table {{ font-size: 0.8rem; }}
        .matrix-table th, .matrix-table td {{ text-align: center; padding: 0.5rem; }}
        .study-header {{ writing-mode: vertical-rl; text-orientation: mixed; }}
        .cond-name {{ text-align: left !important; font-weight: 600; }}
        .testable {{ background: rgba(78, 204, 163, 0.2); color: var(--success); }}
        .not-testable {{ background: rgba(233, 69, 96, 0.2); color: var(--error); }}
        .note {{ color: var(--warning); cursor: help; }}
        .matrix-note {{ color: var(--text-muted); font-size: 0.85rem; margin-bottom: 1rem; }}
        .category-row {{ background: var(--accent); }}
        .category-header {{ font-weight: 700; color: var(--success); text-align: left !important; padding: 0.6rem !important; }}
        .category-divider {{ color: var(--success); border-bottom: 1px solid var(--accent); padding-bottom: 0.3rem; margin: 1.5rem 0 0.5rem 0; font-size: 1.1rem; }}
        .derived-box {{ color: var(--warning); font-style: italic; padding: 0.5rem; text-align: center; background: rgba(249, 168, 37, 0.1); border-radius: 4px; }}
        h6 {{ color: var(--text-muted); font-size: 0.8rem; margin-bottom: 0.3rem; font-weight: 500; }}
        @media (max-width: 1200px) {{
            .condition-grid {{ grid-template-columns: 1fr; }}
            .plot-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Phase A.2: Donor Metadata Harmonization</h1>
        <p class="timestamp">Generated: {timestamp}</p>
        <p>Input: {len(unified_df)} donors × {len(unified_df.columns)} columns → Output: {len(harmonized_df)} donors × 5 conditions</p>
    </div>

    <nav>
        <strong>Navigation:</strong>
        <a href="#overview">Overview</a>
        {nav_items}
        <a href="#cross-study">Summary</a>
        <a href="#testability">Testability</a>
    </nav>

    <div class="overview-section" id="overview">
        <h2>Distribution Overview</h2>
        <p style="color: var(--text-muted); margin-bottom: 1rem; font-size: 0.9rem;">
            Interactive charts: hover for details, click legend to toggle series
        </p>
        <div class="plot-grid">
            <div class="plot-container"><canvas id="donorChart"></canvas></div>
            <div class="plot-container"><canvas id="coverageChart"></canvas></div>
        </div>
        <div class="plot-grid">
            <div class="plot-container"><canvas id="ageContinuousChart"></canvas></div>
            <div class="plot-container"><canvas id="ageBinaryChart"></canvas></div>
        </div>
        <div class="plot-grid">
            <div class="plot-container"><canvas id="parityChart"></canvas></div>
            <div class="plot-container"><canvas id="riskChart"></canvas></div>
        </div>
        <div class="plot-full plot-container"><canvas id="menopausalChart"></canvas></div>
    </div>

    <script>
    const chartData = {chart_data};

    // Common chart options
    const darkTheme = {{
        color: '#e8e8e8',
        borderColor: '#2a2a4a',
        backgroundColor: '#16213e'
    }};

    Chart.defaults.color = '#e8e8e8';
    Chart.defaults.borderColor = '#2a2a4a';

    // Donor counts chart
    new Chart(document.getElementById('donorChart'), {{
        type: 'bar',
        data: {{
            labels: chartData.donorCounts.labels,
            datasets: [{{
                label: 'Donors',
                data: chartData.donorCounts.values,
                backgroundColor: chartData.donorCounts.colors,
                borderWidth: 0
            }}]
        }},
        options: {{
            indexAxis: 'y',
            plugins: {{ legend: {{ display: false }}, title: {{ display: true, text: 'Donors per Study' }} }},
            scales: {{ x: {{ grid: {{ color: '#2a2a4a' }} }}, y: {{ grid: {{ color: '#2a2a4a' }} }} }}
        }}
    }});

    // Coverage heatmap (using bar chart approximation)
    const coverageDatasets = chartData.coverage.conditions.map((cond, i) => ({{
        label: cond,
        data: chartData.coverage.matrix.map(row => row[i]),
        backgroundColor: `hsl(${{i * 60}}, 70%, 50%)`
    }}));
    new Chart(document.getElementById('coverageChart'), {{
        type: 'bar',
        data: {{
            labels: chartData.coverage.studies,
            datasets: coverageDatasets
        }},
        options: {{
            plugins: {{ title: {{ display: true, text: 'Coverage by Study (%)' }}, legend: {{ position: 'bottom', labels: {{ boxWidth: 12 }} }} }},
            scales: {{ x: {{ grid: {{ color: '#2a2a4a' }} }}, y: {{ max: 100, grid: {{ color: '#2a2a4a' }} }} }}
        }}
    }});

    // Age continuous (box plot approximation - showing scatter)
    const ageScatterData = [];
    Object.entries(chartData.ageContinuous.studies).forEach(([study, ages]) => {{
        ages.forEach(age => ageScatterData.push({{ x: study, y: age }}));
    }});
    new Chart(document.getElementById('ageContinuousChart'), {{
        type: 'scatter',
        data: {{
            datasets: Object.entries(chartData.ageContinuous.studies).map(([study, ages]) => ({{
                label: study,
                data: ages.map((age, i) => ({{ x: i * 0.1 - 0.2, y: age }})),
                backgroundColor: chartData.ageContinuous.colors[study],
                pointRadius: 4
            }}))
        }},
        options: {{
            plugins: {{ title: {{ display: true, text: 'Age Distribution (Years)' }}, legend: {{ position: 'bottom', labels: {{ boxWidth: 12 }} }} }},
            scales: {{ x: {{ display: false }}, y: {{ title: {{ display: true, text: 'Age' }}, grid: {{ color: '#2a2a4a' }} }} }}
        }}
    }});

    // Age binary chart
    if (chartData.ageBinary) {{
        new Chart(document.getElementById('ageBinaryChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.ageBinary.labels,
                datasets: [
                    {{ label: 'Young (<50)', data: chartData.ageBinary.young, backgroundColor: '#4ecca3' }},
                    {{ label: 'Old (≥50)', data: chartData.ageBinary.old, backgroundColor: '#e94560' }}
                ]
            }},
            options: {{
                plugins: {{ title: {{ display: true, text: 'Age Binary by Study' }}, legend: {{ position: 'bottom' }} }},
                scales: {{ x: {{ grid: {{ color: '#2a2a4a' }} }}, y: {{ grid: {{ color: '#2a2a4a' }} }} }}
            }}
        }});
    }}

    // Parity chart
    if (chartData.parity) {{
        new Chart(document.getElementById('parityChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.parity.labels,
                datasets: [
                    {{ label: 'Nulliparous', data: chartData.parity.nulliparous, backgroundColor: '#4ecca3' }},
                    {{ label: 'Parous', data: chartData.parity.parous, backgroundColor: '#f9a825' }}
                ]
            }},
            options: {{
                plugins: {{ title: {{ display: true, text: 'Parity by Study' }}, legend: {{ position: 'bottom' }} }},
                scales: {{ x: {{ grid: {{ color: '#2a2a4a' }} }}, y: {{ grid: {{ color: '#2a2a4a' }} }} }}
            }}
        }});
    }}

    // Risk chart
    if (chartData.risk) {{
        new Chart(document.getElementById('riskChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.risk.labels,
                datasets: [
                    {{ label: 'AR (Average Risk)', data: chartData.risk.AR, backgroundColor: '#4ecca3' }},
                    {{ label: 'HR (High Risk)', data: chartData.risk.HR, backgroundColor: '#e94560' }}
                ]
            }},
            options: {{
                plugins: {{ title: {{ display: true, text: 'Risk Status by Study' }}, legend: {{ position: 'bottom' }} }},
                scales: {{ x: {{ grid: {{ color: '#2a2a4a' }} }}, y: {{ grid: {{ color: '#2a2a4a' }} }} }}
            }}
        }});
    }}

    // Menopausal chart
    if (chartData.menopausal) {{
        new Chart(document.getElementById('menopausalChart'), {{
            type: 'bar',
            data: {{
                labels: chartData.menopausal.labels,
                datasets: [
                    {{ label: 'Pre-menopausal', data: chartData.menopausal.pre, backgroundColor: '#4ecca3' }},
                    {{ label: 'Post-menopausal', data: chartData.menopausal.post, backgroundColor: '#f9a825' }}
                ]
            }},
            options: {{
                plugins: {{ title: {{ display: true, text: 'Menopausal Status by Study' }}, legend: {{ position: 'bottom' }} }},
                scales: {{ x: {{ grid: {{ color: '#2a2a4a' }} }}, y: {{ grid: {{ color: '#2a2a4a' }} }} }}
            }}
        }});
    }}
    </script>

    {study_sections}

    <div id="cross-study">
        {cross_study}
    </div>

    <div id="testability">
        {testability}
    </div>

    <footer style="color: var(--text-muted); text-align: center; margin-top: 3rem; padding: 1rem;">
        Phase A.2 Harmonization Report • Spatial HBCA Project
    </footer>
</body>
</html>
'''
    return html_content


def main():
    base_dir = Path(__file__).parent.parent
    unified_path = base_dir / 'outputs' / 'unified_metadata' / 'unified_donor_metadata.csv'
    harmonized_path = base_dir / 'outputs' / 'harmonized_metadata' / 'harmonized_donor_metadata.csv'
    output_path = base_dir / 'outputs' / 'harmonization_audit' / 'harmonization_report.html'

    print(f"Loading {unified_path}")
    unified_df = pd.read_csv(unified_path)
    print(f"  Loaded {len(unified_df)} donors × {len(unified_df.columns)} columns")

    print(f"Loading {harmonized_path}")
    harmonized_df = pd.read_csv(harmonized_path)
    print(f"  Loaded {len(harmonized_df)} donors × {len(harmonized_df.columns)} columns")

    print("Generating HTML report...")
    html_content = generate_html_report(unified_df, harmonized_df)

    print(f"Writing {output_path}")
    with open(output_path, 'w') as f:
        f.write(html_content)

    print(f"Done. Report saved to {output_path}")


if __name__ == '__main__':
    main()
