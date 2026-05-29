# Per-Compartment Integration (scVI / scANVI)

Per-compartment integration of the iHBCA v1 atlas (2.13M cells), producing the
canonical per-compartment latent embeddings, UMAPs, and predicted L1 labels
that downstream annotation, abundance, and communication analyses consume.

## Overview

Each compartment (Immune / Epithelial / Stromal) is integrated independently:

1. **scVI** (`01_compartment_scvi_full.py`) trains a variational autoencoder
   with negative binomial likelihood on the 4000 most highly variable genes,
   conditioning on `dataset` (study) as the batch covariate. Output: 50D
   latent representation `X_scVI`.
2. **scANVI** (`02_compartment_scanvi.py`) initializes from the scVI weights
   and fine-tunes as a semi-supervised classifier using L1 cell-type labels
   (`level1_annotation`). Output: 50D latent `X_scANVI`, predicted labels,
   UMAP, and Leiden clusterings at 10 resolutions.
3. **Comparative baselines** — Harmony (`03_compartment_harmony.py`) and
   Seurat v5 (`04_compartment_seurat.R`, with sketch / rPCA / Harmony
   methods) are included to ground the scVI / scANVI choice in the scIB
   benchmark.
4. **scIB benchmarking** — per-metric sidecar scripts (`05*`) compute NMI,
   ARI, ASW_label, ASW_batch, graph_connectivity, kBET, and PCR; aggregators
   (`06`, `06b`) join the per-config CSVs into a single ranked wide table; an
   SVD diagnostic (`09`) catches single-axis collapse.

## Canonical configuration

| Parameter | Value | Notes |
|---|---|---|
| `n_latent` | **50** | LOCKED — Track A's L2 annotation pipeline depends on these latents |
| `n_layers` | 3 | scVI default |
| `n_hvg` | 4000 | per-batch HVG selection (`batch_key=dataset`) |
| `max_epochs` (scVI) | 300 | early-stopping patience 20 |
| `max_epochs` (scANVI) | 30 | fine-tune; early-stopping patience 10 |
| `n_samples_per_label` | 100 | scANVI semi-supervised sampling |
| `seed` | 42 | numpy + torch + scvi |
| `batch_size` | 128 | both stages |

`n_latent = 75` and `n_latent = 100` were run as comparative sweeps and are
visible in the scIB benchmark, but **n_latent = 50 is the canonical / promoted
configuration**. The two larger sweeps did not improve composite scIB scores
sufficiently to justify the additional latent dimensions.

## Pipeline stages

| Stage | Script | Container | Inputs | Outputs |
|---|---|---|---|---|
| scVI training | `scripts/01_compartment_scvi_full.py` | `python_scvi_gpu_2025Q3` | `{comp}_counts.npz`, `gene_data.csv`, `{comp}_metadata_enriched.csv` | `scvi_model/` (with embedded `adata.h5ad`), `{comp}_leiden.csv`, `{comp}_scvi_umap.csv` |
| scANVI fine-tune | `scripts/02_compartment_scanvi.py` | `python_scvi_gpu_2025Q3` | scVI model + L1 labels | `scanvi_model/model.pt`, `{comp}_leiden.csv`, `{comp}_scanvi_umap.csv`, `{comp}_scanvi_predictions.csv` |
| Harmony baseline | `scripts/03_compartment_harmony.py` | `python_scvi_gpu_2025Q3` | NPZ + metadata | `{comp}_harmony.h5ad` (benchmark only) |
| Seurat baseline | `scripts/04_compartment_seurat.R` | `r_spatial_4.3.3` | 10x HDF5 + metadata | `{comp}_seurat_{method}.rds` (benchmark only) |
| scIB sidecars (per metric) | `scripts/05b/05d/05e/05f/05g/05h_*.py` | `python_scvi_gpu_2025Q3` | integrated h5ad + metadata | per-config CSV per metric |
| scIB aggregator (long) | `scripts/06_scib_aggregate.py` | `python_scvi_gpu_2025Q3` | scib output dir | `scib_summary.csv`, `scib_winners.csv` |
| scIB aggregator (wide) | `scripts/06b_build_wide_csv.py` | `python_scvi_gpu_2025Q3` | scib output dir | `*_aggregated_wide.csv` |
| SVD diagnostic | `scripts/09_metric_svd_diagnostic.py` | `python_scvi_gpu_2025Q3` | wide CSV | stdout report + optional JSON |
| (Deprecated monolith) | `scripts/05_scib_benchmark.py` | — | — | Retained for provenance — do not run |

## Composite scIB scores (canonical n_latent = 50)

| Compartment | Cells | scANVI score | Gap over scVI |
|---|---|---|---|
| Immune | 162,626 | 0.7226 | +0.3234 |
| Epithelial | 991,379 | 0.7761 | +0.3967 |
| Stromal | 974,500 | 0.7584 | +0.3946 |

Composite weighting: `Overall = 0.4 * batch_correction + 0.6 * bio_conservation`
(scIB paper convention).

## Reproducibility

### Path resolution

All inputs / outputs resolve through `publication/config/paths.yaml` via
`load_paths.sh`. The single relevant token is `ihbcav1_assembly_*`:

| Token | Path |
|---|---|
| `ihbcav1_assembly_root` | `iHBCAv1_upload/.dev/assembly_rebuild_20260405/outputs/` |
| `ihbcav1_assembly_components` | `${ihbcav1_assembly_root}/components/` |
| `ihbcav1_assembly_scanvi` | `${ihbcav1_assembly_root}/scanvi/` |

The scIB benchmarking outputs live at `${ihbcav1_assembly_root}/scib/` (peer
to `scanvi/`).

### Submitting the canonical run

```bash
# Source env vars so #SBATCH directives can substitute ${HPC_LOGS_DIR}
source publication/config/load_paths.sh

# 1. scVI training (3-task array, ~24h on 1 A100 each)
sbatch publication/preprocessing/integration/run/run_scvi_full.sh

# 2. scANVI fine-tune (after scVI completes; ~1h per compartment)
sbatch publication/preprocessing/integration/run/run_scanvi_full.sh
```

`N_LATENT` defaults to 50; override via env var to run sweep variants.

### Resources

| Stage | Tier | Cores / Memory / Walltime | GPU |
|---|---|---|---|
| scVI training | 5+GPU | 8 / 128G / 36h | 1 A100 |
| scANVI fine-tune | 5+GPU | 8 / 128G / 8h | 1 A100 |
| scIB sidecars | 3 (graph_conn / nmi_ari) – 4 (asw_*) | 4–8 / 32–64G / 1–6h | none |

## Outputs

Canonical outputs (consumed downstream, referenced via `ihbcav1_assembly_scanvi`):

```
${ihbcav1_assembly_scanvi}/{imm,epi,str}/n_latent_50/
├── {comp}_leiden.csv               # cells × 10 Leiden resolutions
├── {comp}_scanvi_predictions.csv   # cell × (predicted_label, true_label)
├── {comp}_scanvi_umap.csv          # cells × (UMAP_1, UMAP_2)
└── scanvi_model/model.pt           # PyTorch weights (model exception per output-storage.md)
```

Output retention follows `output-storage.md`: only derived intermediates
(CSV / model weights) are tracked as canonical outputs. The full integrated
h5ads are written by the training scripts but are treated as ephemeral —
downstream consumers reconstruct working AnnData objects from CSV + counts at
runtime.

## Consumer contract

Downstream consumers of these outputs (current as of 2026-04-26):

| Consumer | Path | Reads |
|---|---|---|
| Track A annotation pipeline | `publication/analysis/annotation/scripts/00_build_structural_package.py` | `--leiden-csv`, `--umap-csv` |
| Annotation render scripts | `publication/analysis/annotation/run/run_{imm,epi,str}_annotation_pkg.sh` | latents + UMAP + counts |

## Provenance

- Source repo (scripts): `iHBCA_V1/Analysis/stages/V1_Annotation/dev/compartment_integration_20260405/scripts/`
- Source repo (HPC runners): `iHBCAv1_upload/.dev/assembly_rebuild_20260405/run/`
- Promoted: 2026-04-26 (Track 2)
- Out of scope (deferred): `06_full_atlas_markers.py`, `07_singler_train.R`,
  `08_singler_classify.R` belong to a downstream marker / SingleR
  sub-promotion; documented in `SOURCES.md`.
