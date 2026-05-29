# Publication

h5ad assembly, enrichment, validation, and upload preparation for iHBCA v1.0
submission to the HCA Data Portal.

## Subdirectories

| Directory | Purpose |
|-----------|---------|
| `scripts/` | Python/R assembly, enrichment, validation, and entry sheet scripts |
| `run/` | SLURM job wrappers for HPC execution |
| `config/` | Pipeline configuration, dataset registry, metadata definitions |
| `mappings/` | Ontology lookup tables (gene IDs, HANCESTRO, EFO, CL, SRA) |
| `tier2_breast/` | HCA Tier 2 breast-specific metadata proposals |

## Pipeline stages

1. **Assembly** (`assemble_h5ad.py`) -- build per-study CxG-compliant h5ads
2. **Integrated assembly** (`assemble_integrated.py`) -- build all-breast-cells.h5ad
3. **Enrichment** (`enrich_h5ads.py`) -- add var annotations, obs fields, uns metadata
4. **Validation** (`run_hca_validator.py`) -- CAP + HCA schema checks
5. **Entry sheets** (`generate_tier1_*.py`) -- donor/sample metadata for tracker registration

## Configuration

All I/O paths are declared in `config/pipeline.yaml`. Scripts load this
config at startup via `scripts/ihbca/loaders.py`. See `config/pipeline.yaml`
for the complete path inventory.
