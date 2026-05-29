# iHBCA: An integrated Human Breast Cell Atlas with spatial and anatomical resolution

Code repository for the iHBCA Nature paper submission.

## Methods-to-code mapping

| Methods section | Directory |
|----------------|-----------|
| 2.2 FLEX preprocessing | `preprocessing/chromium_spatial/` |
| 2.3 Xenium pipeline | `preprocessing/xenium/` |
| 3 Annotation | `analysis/annotation/` |
| 4 iHBCA assembly | `assembly/` |
| 5.1 Clinical DA | `analysis/abundance/clinical/` |
| 5.2 Joint spatial DA | `analysis/abundance/spatial/joint/` |
| 5.3 FLEX anatomic DA | `analysis/abundance/spatial/flex/` |
| 6 NMF motifs | `analysis/spatial/motifs/` |
| Config | `config/` |

## Directory contents

- `preprocessing/` -- Chromium FLEX and Xenium preprocessing pipelines (scripts, configs, documentation)
- `analysis/` -- Annotation, differential abundance, spatial motif analysis (scripts, configs, documentation)
- `assembly/` -- Frozen reference assembly (v1.0 and v1), including external study harmonization, provenance, and publication upload pipeline
- `config/` -- Shared configuration: paths, aesthetics, annotation hierarchy, containers
- `tables/` -- Producer scripts for supplemental data tables
- `manifest.yaml` -- Promotion manifest tracking source provenance for all scripts

## Data availability

This repository contains analysis code only. Data outputs (count matrices, embeddings, labels, figures) are available through the data repositories described in the manuscript.

## Requirements

Analysis scripts depend on:
- R 4.3+ with Seurat v5, miloR, SingleR, edgeR, BiocParallel
- Python 3.10+ with scanpy, scvi-tools, pandas, anndata
- HPC execution via SLURM with Singularity containers (see `config/containers.yaml`)

Path resolution uses `config/paths.yaml`. Adjust paths to match your data locations before execution.
