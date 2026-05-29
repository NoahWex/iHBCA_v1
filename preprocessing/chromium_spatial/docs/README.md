# Chromium FLEX Preprocessing

22-stage preprocessing pipeline for the Spatial HBCA Chromium FLEX dataset
(258K cells, 4 patients × 17 anatomic positions per breast). Produces the
canonical scvi_n100 bundle (latent + UMAP + obs + counts) consumed by the
joint preprocessing, FLEX annotation, and FLEX-only DA tracks.

## Pipeline overview

| # | Stage | Language | Status |
|---|-------|----------|--------|
| 01 | baseline_vfs | R | promoted |
| 02 | cell_filtering (per-cluster MAD) | Python | promoted |
| 03 | doublet_detection (scDblFinder) | R | promoted |
| 04 | label_transfer (SingleR + harmonize) | Python / R | promoted |
| 05 | post_qc_vfs (union VFs) | R | promoted |
| 06 | filtered_preview (AnnData) | Python | promoted |
| 07 | final_report (per-sample QC HTML) | Python | promoted |
| 08 | scvi_integration (canonical, 50D latent) | Python | promoted |
| 09 | integration_preview (UMAP, clusters) | R | promoted |
| 10a | add_residual_normalized_data | R | promoted |
| 10b | interactive_exploration | R | **FROZEN** (interactive Rmd; canonical manifest is the frozen output) |
| 11 | filtered_vfs (compartment basis) | R | promoted |
| 12 | filtered_scvi_integration (per-compartment) | Python | promoted |
| 13 | filtered_integration_preview | R | promoted |
| 14 | compartment_classification (multinomial) | R | promoted |
| 15 | milo_contamination | R | promoted (substages 1–2); **FROZEN** (substages 3–4 interactive) |
| 16 | bigsur_vfs | R | promoted |
| 17 | integration_sweep (9 configs × 4 targets, scIB-scored) | Python | promoted |
| 18 | deseq (per-compartment, per-resolution) | R | promoted |
| 19a | l1_annotation | R / Python | **FROZEN** (parallel pathway; superseded by Track C) |
| 19b | integration_preview (Rmd, sweep-winner) | R | promoted |
| 20 | finalize (canonical scvi_n100 bundle producer) | Python | promoted |

Frozen stages are retained for provenance and audit. `FROZEN.md` in each
stage explains what the stage does, why it's frozen, what consumes its
outputs, and where the canonical alternative lives.

## Raw input contract

Step 01 reads per-sample H5 files referenced by a sample manifest TSV. The
manifest holds paths to the lab's CRSP share and is therefore not promoted as
a runnable artifact. Its schema is documented here so a reviewer or replicator
can supply an equivalent file for their own data.

### Manifest schema

```
sample_id<TAB>h5_path
```

| Column | Type | Description |
|---|---|---|
| `sample_id` | string | Unique identifier per sample. Convention: `{patient_id}_{position_id}` (e.g., `Pat1_P1`, `Pat1_P2_Lower`, `UCI604_P3_A_LIQ`). |
| `h5_path` | absolute path | Path to the per-sample `cellranger multi` output H5 (specifically `sample_filtered_feature_bc_matrix.h5` under `per_sample_outs/{position}/count/`). |

Header row required. Tab-separated. UTF-8.

### Example row (anonymized)

```
Pat1_P1<TAB>/share/lab/HBCA_Phase2/fixed_scRNAseq/Pat1/CellRangerMulti/per_sample_outs/P1/count/sample_filtered_feature_bc_matrix.h5
```

The canonical manifest used for the published atlas resolves via the
`raw_data_sample_manifest` token in `config/paths.yaml`.

### CellRanger version

Per-sample H5s were produced by `cellranger multi` (10x Genomics, fixed-RNA / FLEX panel). Reference build, panel selection, and CellRanger version are documented in the methods section of the published paper. The H5 schema is the standard `sample_filtered_feature_bc_matrix.h5` output of `cellranger multi` v8.x or compatible.

### Path resolution for replication

To reproduce the pipeline against a different dataset:
1. Produce per-sample H5s with `cellranger multi` (matching panel + reference build).
2. Build a manifest TSV following the schema above.
3. Update `raw_data_sample_manifest` in `config/paths.yaml` (or pass via the `--raw-data-manifest` CLI flag accepted by every consumer wrapper).

## Directory layout

```
publication/preprocessing/chromium_spatial/
├── README.md                       (this file)
├── config/
│   └── paths.yaml                  (74 internal pipeline tokens; references publication paths.yaml)
├── run/                            (45 SLURM wrappers + _common.sh)
└── scripts/
    ├── step_01_baseline_vfs/
    ├── step_02_cell_filtering/
    ├── ...
    ├── step_19a_l1_annotation/     (FROZEN)
    ├── step_19b_integration_preview/
    └── step_20_finalize/
```

`run/_common.sh` sources `publication/config/load_paths.sh`, bridges the
publication-level path tokens (`USER_ROOT`, `SHARED_DATASETS`,
`HPC_LOGS_DIR`, `BIND_MOUNT_LIST`) into the pipeline's internal namespace
(`CFG_*`), and exposes a `run_r_singularity` / `run_py_singularity` helper
that the wrappers use to invoke containerized commands with the right
bind mounts.

## How to run a stage

All stages run via SLURM wrappers in `run/`. Sourcing
`publication/config/load_paths.sh` first is mandatory — the wrapper
`#SBATCH --output=${HPC_LOGS_DIR}/...` directives are substituted at submit
time, not by SLURM.

```bash
# 1. Source publication path tokens
source publication/config/load_paths.sh

# 2. Submit a stage (example: step 8 scVI integration)
sbatch publication/preprocessing/chromium_spatial/run/run_step_08_scvi_integration.sh

# 3. Monitor
squeue -u $USER
```

Each wrapper:
- Sources `_common.sh` to populate `CFG_*` tokens from `paths.yaml`.
- Resolves the underlying script via `${CFG_SCRIPTS_DIR}/step_NN_*/`.
- Runs the script in the appropriate Singularity container
  (`r_spatial`, `python_scvi`, etc., per `containers.yaml`).
- Writes outputs to the stage's canonical directory under
  `${CFG_PREPROCESSING_ROOT}/`.

## Pipeline order and dependencies

The 22 stages are nominally linear (01 → 20), but a few stages have
non-obvious dependencies:

- Stage 10b (interactive) is the canonical L2 cell gate. Its frozen output
  (`step10b_manifest_update.csv`) is consumed by stages 11 / 12 / 13 / 14
  as a 5-way QC filter. Re-running 10b invalidates downstream outputs.
- Stage 15 substages 1–2 produce inputs for substages 3–4 (interactive),
  which produce the frozen `contamination_specific_cells.csv` consumed by
  stages 16, 17, and 19b.
- Stage 17 outputs are scored by the sweep aggregator; the winner sidecars
  (`${CFG_CANONICAL_SWEEP_WINNER}/...`) are inputs to step 19b and step 20.
- Stage 18 (DESeq) and step 19b operate on the same step 17 sweep winner
  but produce independent outputs and can run in parallel after step 17.
- Stage 20 finalize requires step 17 winner sidecars + step 18 DESeq
  outputs + step 14 compartment labels + step 15 contamination flags to
  build the canonical `scvi_n100` bundle.

## Output policy

Canonical outputs follow the project `output-storage.md` rule: CSV + MTX
intermediates only, no h5ad in canonical paths. Working AnnData/Seurat
objects are reconstructed at runtime from the bundle (`step_20_finalize/README.md`
describes the bundle schema).

## Downstream consumers

| Track | Reads from this pipeline | Where |
|-------|--------------------------|-------|
| Track B (joint preprocessing) | `step_20_finalize` bundle (full + per-compartment) | `publication/preprocessing/joint/` |
| Track C (FLEX annotation) | `step_20_finalize` bundle | `publication/preprocessing/chromium_spatial/annotation/` |
| Track D (FLEX-only DA) | `step_20_finalize` bundle + Track C labels | `publication/analysis/abundance/spatial/flex/` |
| Existing analysis tracks | `step_20_finalize` bundle | `publication/analysis/abundance/spatial/`, `publication/analysis/spatial/` |

## Provenance

- Source repo: `Spatial_HBCA_preprocessing/` (22-stage pipeline; promoted as Track A)
- Promotion date: 2026-04-27
- Promotion log: `publication/manifest.yaml`
