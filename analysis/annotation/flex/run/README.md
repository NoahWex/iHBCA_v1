# FLEX L2S annotation — pipeline runners

Per-cell L2S labels for the Spatial HBCA Chromium FLEX cohort, produced by an
anchor-resolution + fine-cluster-override review pipeline. The atlas-wide
deliverable is `outputs/flex_l2s_labels.csv` (270,035 cells x 6 columns).

## Pipeline shape

The pipeline has three phases:

1. **Pre-review structural diagnostics** (Steps 1-5d). Computes Leiden over a
   resolution sweep, scores V1 marker signatures (UCell), runs limma-voom de
   novo markers, computes assignment + concordance + novelty, and renders the
   structural diagnostic package per compartment.
2. **Manual review**. The reviewer reads `outputs/{Compartment}/structural/`
   plus per-V1-label feature plots and canonical marker heatmaps and writes
   `annotation_yamls/{compartment}.yaml` declaring the anchor labels and any
   fine_cluster overrides.
3. **Post-review label emission** (Steps 6-9). Renders the per-compartment
   evidence package, emits per-cell labels CSVs, concatenates them into the
   atlas table, and renders the atlas-wide L2S UMAP.

## Wrapper inventory

| Wrapper | Shape | Tier | Description |
|---|---|---|---|
| `run_01_cluster.sh` | array 0-2 (per compartment) | 2.5 (24G/4CPU/1.5h) | Extract precomputed Leiden + recompute {1.5, 3.0} via leidenalg+igraph |
| `run_02_ucell.sh` | array 0-2 | 3 (32G/4CPU/1.5h) | UCell rank-based scoring per cell + per-cluster aggregation, all 7 resolutions |
| `run_03a_pseudobulk.sh` | array 0-20 (compartment x resolution) | 2.5 (24G/4CPU/1h) | Pseudobulk counts to (cluster x library_id) per (compartment, resolution) |
| `run_03b_limma.sh` | array 0-20 | 3 (32G/4CPU/2h) | Limma-voom one-vs-rest markers per cluster, ~ library_id + condition |
| `run_04_assignment.sh` | array 0-20 | 1.5 (8G/2CPU/30min) | Best-match L2S + novelty heuristics per (compartment, resolution) |
| `run_04b_concordance.sh` | array 0-2 | 2 (16G/4CPU/45min) | Cross-resolution F1 concordance per V1 label |
| `run_05_structural.sh` | array 0-2 | 3 (32G/4CPU/1.5h) | Pre-review structural diagnostic package |
| `run_05c_featureplots.sh` | array 0-2 | 3 (32G/4CPU/1.5h) | Per-V1-label UMAP feature plots |
| `run_05d_canonical_heatmap.sh` | array 0-2 | 3 (32G/4CPU/1.5h) | V1 canonical+identity markers x cluster heatmap per resolution |
| `run_06_annotation.sh` | array 0-2 | 3 (32G/4CPU/1.5h) | Post-review evidence package per compartment (requires yaml) |
| `run_07_resolve_l2s_labels.sh` | single | 2 (8G/2CPU/30min) | Resolve per-cell L2S labels for all 3 compartments from the reviewed yamls + clusters CSVs; writes the canonical `flex_l2s_labels.csv` |
| `run_09_full_atlas_umap.sh` | single | 2 (16G/2CPU/30min) | Project atlas labels onto integrated UMAP |

## Array index encoding

For Steps 3a, 3b, 4 (21 tasks total):

```
COMPARTMENTS=(Epithelial Immune Stromal)             # 3
RESOLUTIONS=(0.3 0.5 0.8 1.0 1.5 3.0 5.0)            # 7
C_IDX=$((SLURM_ARRAY_TASK_ID / 7))
R_IDX=$((SLURM_ARRAY_TASK_ID % 7))
```

## Submission order

```bash
JOB1=$(sbatch --parsable run_01_cluster.sh)
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 run_02_ucell.sh)
JOB3a=$(sbatch --parsable --dependency=afterok:$JOB1 run_03a_pseudobulk.sh)
JOB3b=$(sbatch --parsable --dependency=afterok:$JOB3a run_03b_limma.sh)
JOB4=$(sbatch --parsable --dependency=afterok:$JOB2,$JOB3b run_04_assignment.sh)
JOB4b=$(sbatch --parsable --dependency=afterok:$JOB2 run_04b_concordance.sh)
JOB5=$(sbatch --parsable --dependency=afterok:$JOB4 run_05_structural.sh)
JOB5c=$(sbatch --parsable --dependency=afterok:$JOB1 run_05c_featureplots.sh)
JOB5d=$(sbatch --parsable --dependency=afterok:$JOB4 run_05d_canonical_heatmap.sh)
# --- Manual review: produces annotation_yamls/{compartment}.yaml ---
JOB6=$(sbatch --parsable run_06_annotation.sh)
JOB7=$(sbatch --parsable --dependency=afterok:$JOB6 run_07_resolve_l2s_labels.sh)
JOB9=$(sbatch --parsable --dependency=afterok:$JOB7 run_09_full_atlas_umap.sh)
```

## Path resolution

All wrappers source `_common.sh`, which sources `publication/config/load_paths.sh`
(located by walking up from the wrapper's directory). `_common.sh` resolves:

- `CFG_ANNOTATION_ROOT` — this directory's parent (the annotation tree)
- `CFG_SCRIPTS_DIR`, `CFG_OUTPUTS_ROOT`, `CFG_ANNOTATION_YAMLS` — relative to `CFG_ANNOTATION_ROOT`
- `CFG_INTEGRATION_INTERMEDIATE` — FLEX scvi_n100 bundles per compartment, from
  `SOURCE_FLEX_INTEGRATION_INTERMEDIATE` in publication paths.yaml
- `CFG_V1_ANNOTATION_YAMLS` — V1 marker signature yamls at
  `${PROJECT_PUBLICATION}/publication/analysis/annotation/yamls`
- Container, R libs, and bind mounts inherited from publication `load_paths.sh`

## Inputs

Per-compartment scvi_n100 bundles produced upstream by the FLEX preprocessing
pipeline:

```
${CFG_INTEGRATION_INTERMEDIATE}/{Compartment}/scvi_n100/
    cells.tsv           cell_id index
    genes.tsv           gene index
    counts.mtx.gz       sparse counts (cells x genes)
    latent.csv          50-D scVI latent
    umap.csv            2-D UMAP
    obs.csv             cell metadata (sample_id, library_id, patient_id, ...)
```

Plus V1 marker signature yamls (`annotation_v2_{epi,imm,str}.yaml`) defining
canonical and identity markers per V1 L2 label.

## Outputs

```
outputs/
    flex_l2s_labels.csv                          # canonical atlas table
    clusters/{Compartment}_clusters.csv
    ucell/{Compartment}_ucell_per_{cell,cluster}.csv
    pseudobulk/{Compartment}_leiden_{res}_pseudobulk_{counts,meta}.csv
    limma/limma_top_markers_{Compartment}_leiden_{res}.csv
    assignments/{Compartment}_leiden_{res}_assignments.csv
    concordance/{Compartment}_resolution_concordance.{csv,pdf}
    {Compartment}/
        structural/                              # pre-review diagnostics
        annotation/                              # post-review evidence (heatmaps, label_summary, label_umap)
```

Step 7 (`run_07_resolve_l2s_labels.sh`) reads the reviewed
`annotation_yamls/{compartment}.yaml` files plus
`outputs/clusters/{Compartment}_clusters.csv` and writes a single
concatenated `flex_l2s_labels.csv` directly to `outputs/`. There is no
intermediate per-compartment `labels.csv` step.

The promotion bundle keeps `outputs/flex_l2s_labels.csv` and
`outputs/per_compartment/{Compartment}/annotation/` only. Intermediate tables
(clusters, ucell, pseudobulk, limma, assignments, concordance, structural) are
regenerated by the pipeline if needed.
