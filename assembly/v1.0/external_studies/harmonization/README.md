# Harmonization

Cross-study donor metadata harmonization pipeline. Standardizes clinical
covariates, demographics, and sample metadata across all 7 source studies
into a unified schema (287 donors x 23 covariates).

## Harmonized covariates

| Covariate | Description | Coverage |
|-----------|-------------|----------|
| age_continuous | Age in years | 56% (Kumar binary only) |
| age_binary | <50 vs >=50 | 92% |
| parity_binary | Nulliparous vs Parous | 78% |
| risk_status_binary | Average risk vs High risk (BRCA) | 65% |
| menopausal_status_binary | Pre vs Post | 44% |

Additional covariates: parity count, age at first birth, BRCA genotype,
cancer history, tissue indication, menopausal status (detailed), ethnicity
(verbatim + grouped), BMI, sample preservation, sample type, FACS status,
dissociation time.

## Directory structure

| Directory | Purpose |
|-----------|---------|
| `scripts/` | R and Python harmonization code |
| `run/` | SLURM job wrappers |
| `config/` | Schemas, manifests, extraction format definitions |
| `common/` | Cell ID canonicalization config, process guides |
| `studies/` | Per-study config, extracted data, reference tables, edge case docs |

## Per-study structure (`studies/{study}/`)

| Subdirectory | Purpose |
|--------------|---------|
| `config.yaml` | Study-specific harmonization configuration |
| `STUDY.md` | Study overview and harmonization notes |
| `extracted/` | Structured YAML extractions from supplemental tables |
| `reference/` | Published supplemental tables (CSV, XLSX) |
| `knowledge/` | Edge case documentation and data quality notes |
| `scripts/` | Study-specific extraction scripts (if needed) |
| `run/` | Study-specific SLURM wrappers (if needed) |

## Studies

| Study | Donors | Cells | Notes |
|-------|--------|-------|-------|
| Gray | 10 | 52,681 | WT + BRCA1/BRCA2 carriers |
| Kumar | 126 | 714,331 | Largest study, binary age only |
| Murrow | 16 | 86,136 | All average risk, all premenopausal |
| Nee | 22 | 230,100 | BRCA1 vs controls |
| Twigger | 13 | 110,744 | Lactating + non-lactating tissue |
| Pal (3 sub-cohorts) | 28 | 131,288 | Epithelial, total, BRCA1 |
| Reed | 72 | 803,283 | iHBCA core reference study |

## Reproduction

Pipeline steps (run in order):

```bash
BASE=publication/run

# 1. Extract iHBCA reference cell inventory
sbatch $BASE/../external_studies/harmonization/run/extract_ihbca_reference.sh

# 2. Build unified donor metadata (merge all study metadata)
sbatch $BASE/../external_studies/harmonization/run/build_donor_metadata.sh

# 3. Audit metadata columns
sbatch $BASE/../external_studies/harmonization/run/run_audit_metadata.sh

# 4. Compute coverage matrix
sbatch $BASE/../external_studies/harmonization/run/run_coverage_matrix.sh

# 5. Build harmonized donor table
sbatch $BASE/../external_studies/harmonization/run/run_build_harmonized.sh

# 6. Map harmonized metadata to cells
sbatch $BASE/../external_studies/harmonization/run/run_map_cells.sh

# 7. Generate HTML report
sbatch $BASE/../external_studies/harmonization/run/run_harmonization_report.sh
```

Steps 3-5 depend on Step 2. Steps 6-7 depend on Step 5.

## Output

`outputs/harmonized_metadata/harmonized_donor_metadata.csv` (gitignored) --
consumed by `publication/scripts/assemble_h5ad.py` and entry sheet generators.
