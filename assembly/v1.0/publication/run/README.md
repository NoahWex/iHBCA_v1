# Run

SLURM job wrappers for executing pipeline scripts on UCI HPC3.

## Configuration

All scripts source `load_config.sh` which reads `publication/config/pipeline.yaml`
and exports path variables. `REPO_ROOT` must be set before sourcing.

## Key scripts

| Script | SLURM Type | Purpose |
|--------|-----------|---------|
| `reassemble_all_source.sh` | Array (0-6) | Assemble all 7 source h5ads |
| `assemble_integrated.sh` | Single | Assemble all-breast-cells.h5ad |
| `enrich_source.sh` | Array (0-6) | Enrich all 7 source h5ads |
| `enrich_integrated.sh` | Single | Enrich integrated h5ad |
| `validate_all.sh` | Single | Run 3 validators on all 8 h5ads |
| `B1_package_external.sh` | Orchestrator | Full extract + assemble pipeline |

## Supporting scripts

| Script | Purpose |
|--------|---------|
| `load_config.sh` | Export pipeline.yaml paths as shell variables |
| `build_sketch.sh` | Generate sketch h5ad |
| `build_cl_crosswalk.sh` | Build CL term crosswalk tables |
| `build_sra_run_tables.sh` | Extract SRA run tables |
| `build_level15_mapping.sh` | Build level 1.5 to CL mapping |
| `diff_all_source.sh` | Diff all source h5ads against baseline |
| `diff_integrated.sh` | Diff integrated h5ad against baseline |
| `generate_tier1_donor.sh` | Generate donor entry sheets |
| `submit_upload.sh` | Submit to HCA via CAP |

## Resource requirements

| Job | Memory | Time | Notes |
|-----|--------|------|-------|
| Source assembly (reed) | 128 GB | 2h | Largest dataset (803K cells) |
| Integrated assembly | 180 GB | 2h | 2.12M cells |
| Integrated enrichment | 256 GB | 6h | Includes log-normalized layer |
| Validation (all) | 180 GB | 6h | Sequential, 3 validators x 8 files |
