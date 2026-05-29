# Config

Pipeline configuration and metadata definitions.

## Files

| File | Purpose |
|------|---------|
| `pipeline.yaml` | **Primary config** -- single source of truth for all I/O paths, study definitions, and mapping file locations. Loaded by all scripts at startup. |
| `source_dataset_registry.yaml` | Dataset inventory: output filenames, input paths, cell counts, packaging status per study |
| `dataset_metadata.yaml` | Per-study dataset-level metadata (PI, DOI, GEO accession, alignment software, etc.) populated into h5ad `uns` |
| `ontology_mappings.yaml` | Ontology mapping coverage tracker (HANCESTRO, HsapDv, EFO, CL status) |
| `metadata_gap_tracker.yaml` | Per-study x per-field gap analysis for HCA compliance |

## metadata_stages/

Per-study methods extraction notes and the harmonized donor metadata staging CSV:

| File | Purpose |
|------|---------|
| `L1_harmonized_donor.csv` | Harmonized donor metadata (287 donors x 23 covariates) broadcast into h5ad obs |
| `{study}_tier1_methods.md` | Per-study notes on alignment software, gene annotation version, and other dataset-level fields extracted from publication methods sections |
