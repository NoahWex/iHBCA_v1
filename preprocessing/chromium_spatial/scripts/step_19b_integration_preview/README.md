# Step 19b — Integration Preview

Per-compartment QC report of the integrated atlas. Reassembles a Seurat
object from the step 17 sweep winner sidecars + step 14 compartment
labels + step 15 contamination flags + raw H5 counts, runs per-cluster
marker discovery at each Leiden resolution, and renders an HTML preview
with UMAP composition, retention tracking, and marker tables.

## Inputs

- `$CFG_CANONICAL_SWEEP_WINNER/compartments/{Epithelial,Stromal,Immune}/{metadata_with_umap.csv, umap.csv, leiden_multi.csv}` — locked step 17 sweep winner sidecars
- `$CFG_PREPROCESSING_STEP_14/metadata/cell_compartments.csv` — compartment overlay (validation)
- `$CFG_PREPROCESSING_CONTAMINATION_CELLS` — step 15 contamination overlay (validation)
- raw 10x H5 counts per sample (looked up from raw-data manifest)

The Rmd reads Leiden resolutions from `sweep_configs.yaml` rather than
from the legacy `step_17_scvi_integration.algorithm_params` config block.

## Outputs

Per compartment, written under `$CFG_PREPROCESSING_STEP_19B/{COMPARTMENT}/`:

| Subdirectory | Contents | Promoted? |
|---|---|---|
| `reports/` | `integration_preview_${COMPARTMENT}.html` | yes (per-stage QC) |
| `markers_by_resolution/` | `markers_res_*.csv`, `marker_quality_sweep.csv` | yes |
| `metadata/` | `unified_metadata.csv` | yes |
| `seurat_objects/` | `integrated_seurat.rds` | **no — stage-internal scratch** |
| `logs/` | `wrapper_${TIMESTAMP}.log` | n/a |

`seurat_objects/integrated_seurat.rds` is a stage-internal Seurat object
used by the Rmd's marker / plotting chunks. It is regenerable on re-run
from the inputs listed above and is **not** consumed by any other promoted
stage. Per `output-storage.md`, it is not a canonical promoted output.
Downstream tracks (B / C / D) consume the canonical `scvi_n100` bundle
produced by step 20 finalize (CSV + MTX), not this RDS.

## Wrapper

- `wrapper.R` — canonical-path wrapper with validation preamble + dry-run
  flag. Consumes `CFG_CANONICAL_SWEEP_WINNER`, `CFG_PREPROCESSING_STEP_14`,
  `CFG_PREPROCESSING_CONTAMINATION_CELLS`, `CFG_PREPROCESSING_STEP_19B`.
- `../../run/run_step_19b_integration_preview.sh` — SLURM wrapper, Tier 4
  (8 cpus / 64 GB / 3h, 3-task array over compartments), sources
  `run/_common.sh`.

## Dry-run support

The wrapper's `--dry-run` flag (also reachable via `DRY_RUN=1` env on the
SLURM wrapper) validates input paths, sweep winner sidecar presence, and
output directory writability without rendering the notebook. Useful as a
fast pre-flight before launching the Tier 4 render array.
