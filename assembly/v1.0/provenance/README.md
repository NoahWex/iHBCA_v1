# Provenance

Original construction code for the iHBCA v1.0 integrated atlas.

## original_construction/

Austin Marckx's scripts that built the published iHBCA CELLxGENE deposit.
Archived here for reproducibility and provenance documentation.

| Directory | Purpose |
|-----------|---------|
| `iHBCA_build/` | Per-study data preparation (`10{a-f}_*_data_preparation.{R,py}`), atlas merging (`merge_iHBCA.py`), CAP object construction |
| `integration/` | scVI, Harmony, and scPoli integration runs |
| `integration_metrics1/` | Integration quality metrics (scIB-style) |
| `integration_metrics2/` | CellTypist leiden analysis and comparison |
| `integration_metrics3/` | Confusion matrices, pseudobulk PCR, embedding analysis |
| `misc/` | Additional analysis scripts |

These scripts are not part of the active pipeline. The upload assembly code
in `publication/scripts/` takes the outputs of this construction pipeline
as its starting inputs.

## Frozen state

Files under `original_construction/` are preserved exactly as published. Hardcoded paths will not resolve in the publication container. Do not re-run — the canonical pipeline is `publication/scripts/`.

## Unpromoted validation utilities

Five post-hoc utilities live in `iHBCAv1_upload/scripts/` (`fix_gray_barcodes.py`, `check_gene_overlap.py`, `check_mh0023_drop.py`, `check_unmatched.py`, `check_ethnicity_by_study.py`). Redundant with checks integrated into `external_studies/harmonization/`. Not promoted.
