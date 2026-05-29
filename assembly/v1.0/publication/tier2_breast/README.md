# Breast Tier 2 Metadata

HCA Tier 2 breast-specific metadata field proposals for the Atlas Tracker.

Tier 2 fields extend the standard HCA Tier 1 schema with tissue-specific
covariates relevant to breast biology (parity, menopausal status, BRCA
genotype, etc.).

## Files

| File | Purpose |
|------|---------|
| `breast_tier2_draft.csv` | Initial draft of proposed breast Tier 2 fields |
| `breast_tier2_conservative.csv` | Conservative subset (fields with >70% coverage) |
| `breast_tier2_idealized.csv` | Full idealized field list |
| `breast_tier2_prop.csv` | Final proposal submitted for review |
| `breast_tier2_original.csv` | Original field list from EBI |
| `breast_tier2_field_audit.csv` | Per-field coverage audit across 7 studies |
| `breast_tier2_combined_audit.csv` | Combined coverage summary |

## Context

Requested by Ida Zucchi (EBI) as part of the HCA Atlas Tracker Tier 2
metadata standardization effort. The breast proposal follows the format
established by the Lung Tier 2 template (68 fields).

Fields were cross-referenced against the harmonized donor metadata
(287 donors x 23 covariates) to assess population feasibility.

## Contributors

Kessenbrock Lab (UCI), Antony Rose (Newcastle), Walid Khaled (Cambridge)
