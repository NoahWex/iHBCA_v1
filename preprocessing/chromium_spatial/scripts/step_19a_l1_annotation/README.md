# step_19a_l1_annotation — L1 annotation (FROZEN)

> **FROZEN — see `FROZEN.md` in this directory.** Scripts are retained for
> provenance only and are not part of the canonical chromium_spatial
> preprocessing flow. Track C (FLEX annotation pipeline, planned) supersedes
> this stage.

The reference content below describes the historical pipeline as it existed
when frozen, for reviewers tracing the project's evolution.

## Pipeline order (historical reference)

| Step | Script | Wrapper | Output |
|------|--------|---------|--------|
| 09 | `09_build_l1_reference.py` | `run_build_l1_reference.sh` | iHBCA reference SCE/h5ad |
| 09a | `09a_build_l1_seurat.R` | `run_build_l1_seurat.sh` | per-comp `seurat.rds` |
| 09b | `09b_singler_l1.Rmd` | `run_singler_l1.sh` | `singler_l1_predictions_{C}.csv` + `.html` |
| 09c | `09c_l1_harmonize.Rmd` | `run_l1_harmonize.sh` | `harmonized_l1_labels_{C}.csv` (per-comp, leiden-voted) |
| 09d | `09d_concat_l1.R` | `run_concat_l1.sh` | `all/harmonized_l1_labels.csv` (concat, ~270K rows) |
| 09e | `09e_l1_confusogram.Rmd` | `run_l1_confusogram.sh` | confusogram HTML |
| 09f | `09f_smooth_full_l1.py` | `run_smooth_full_l1.sh` | `all/harmonized_l1_labels_full.csv` |

## Why 09f exists (historical context)

The canonical 5-way QC filter retains 258,317 cells. Per-compartment Seurat
builds (09a) drop ~19,560 of these — cells that pass canonical QC but fail
compartment classification or per-comp Seurat assembly. They have full-object
integration (latent + leiden) but no L1 label.

09f voted within `leiden_2.0` clusters on the full object to assign L1 to
these gap cells. Output: `harmonized_l1_labels_full.csv`, the canonical
full-object L1 label source. Coverage 99.9%, with `L1.0_source`
distinguishing `compartment_harmonize` (238,757 from 09c) vs `leiden_vote`
(19,294 from 09f) vs `unsmoothed` (266 cells in clusters too small to vote).
