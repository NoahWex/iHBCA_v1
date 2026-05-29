# Scripts

Assembly, enrichment, and validation code for iHBCA h5ad objects.

## Core pipeline scripts

| Script | Purpose |
|--------|---------|
| `assemble_h5ad.py` | Build one CxG-compliant h5ad per source study from extracted intermediates |
| `assemble_integrated.py` | Build all-breast-cells.h5ad (2.12M cells) from author share primitives |
| `enrich_h5ads.py` | Post-assembly enrichment: var annotations, obs fields, uns metadata, layers |
| `run_hca_validator.py` | Run HCA-specific schema validation checks |
| `validate_config.py` | Validate pipeline.yaml config integrity and cross-references |

## Supporting scripts

| Script | Purpose |
|--------|---------|
| `build_sketch.py` | Generate sketch (subsampled) h5ad for fast computation |
| `build_cl_crosswalk.py` | Build cell type ontology crosswalk tables for non-CxG studies |
| `build_sra_run_tables.py` | Extract SRA run accessions per study from NCBI metadata |
| `build_gtf_gene_mapping.py` | Parse GENCODE GTF to build gene symbol/Ensembl mapping |
| `build_level15_mapping.py` | Build level 1.5 annotation to CL ontology term mapping |
| `diff_h5ads.py` | Compare two h5ads for validation (cell counts, obs, var, X) |
| `diff_against_baseline.py` | Diff regenerated h5ads against baseline snapshots |
| `snapshot_baseline.py` | Save baseline h5ad checksums for regression testing |

## Entry sheet generation

| Script | Purpose |
|--------|---------|
| `generate_entry_sheets.py` | Generate dataset-level entry sheets for HCA tracker |
| `generate_tier1_donor_sheet.py` | Per-study donor metadata CSVs for tracker registration |
| `generate_tier1_sample_sheet.py` | Per-study sample metadata CSVs for tracker registration |

## Shared library (`ihbca/`)

| Module | Purpose |
|--------|---------|
| `loaders.py` | Config-driven data loading (mappings, configs, donor metadata) |
| `constants.py` | Schema field lists, ontology mappings, config-driven study helpers |
| `hca_fields.py` | HCA-required obs field population and dataset metadata utilities |
| `ethnicity.py` | HANCESTRO term normalization and downgrade logic |
| `validation.py` | Tier 1 field validation, Ensembl coverage checks, UMAP validation |

## Patch scripts

One-time fixes applied during assembly iteration:

| Script | Purpose |
|--------|---------|
| `patch_obs_fields.py` | Patch obs columns in assembled h5ads |
| `patch_obs_h5py.py` | Low-level h5py obs patching (avoids full anndata load) |
| `patch_uns_title.py` | Fix uns title fields |
| `patch_var_uns_h5py.py` | Patch var and uns via h5py |
| `patch_remove_feature_is_filtered_raw.py` | Remove feature_is_filtered from raw layer |
| `patch_remove_schema_version.py` | Remove stale schema_version from uns |
| `restructure_raw_x.py` | Restructure raw X layer for schema compliance |
