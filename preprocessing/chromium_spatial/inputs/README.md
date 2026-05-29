# Raw input contract

This directory holds the input boundary for the FLEX preprocessing pipeline.

## raw_data_manifest.template.tsv

Template manifest file. The canonical published atlas was produced from a manifest with 63 sample rows pointing at lab-specific `cellranger multi` H5 files. The canonical manifest itself is not promoted — it contains lab-internal absolute paths to the CRSP share that no external reviewer can access.

To replicate the pipeline against your own data, fill in this template with rows mapping `sample_id` → absolute path of your `cellranger multi` per-sample H5 (`sample_filtered_feature_bc_matrix.h5`).

## Schema

```
sample_id<TAB>h5_path
```

| Column | Description |
|---|---|
| `sample_id` | Unique identifier per sample. Convention: `{patient_id}_{position_id}` (e.g., `Pat1_P1`, `Pat1_P2_Lower`). The pipeline derives `patient_id` (first underscore-delimited token) and `position_id` (remainder) from this string. |
| `h5_path` | Absolute path to per-sample `cellranger multi` H5 output. Must be readable by the SLURM compute nodes. |

Header row required. Tab-separated. UTF-8.

## CellRanger version

Per-sample H5s were produced by `cellranger multi` (10x Genomics, fixed-RNA / FLEX panel). The H5 schema is the standard `sample_filtered_feature_bc_matrix.h5` produced by `cellranger multi` v8.x or compatible. Reference build, panel selection, and exact CellRanger version are documented in the methods section of the published paper.

## Path resolution

The `flex_raw_data_manifest` token in `publication/config/paths.yaml` resolves to the manifest TSV location. Override at runtime via the `--raw-data-manifest` CLI flag accepted by every consumer wrapper, or by replacing this file in a local fork.
