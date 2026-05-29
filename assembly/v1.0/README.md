# iHBCA v1.0 — Upload Code Repository

Code and configuration for assembling, validating, and uploading the
**integrated Human Breast Cell Atlas (iHBCA) v1.0** to the
[Human Cell Atlas Data Portal](https://data.humancellatlas.org/).

> **Status: frozen reference.** This tree is a methodologically-inspectable archive,
> not an executable pipeline from the publication tier. The canonical h5ad lives on
> CRSP and is referenced downstream via `paths.yaml`. End-to-end re-execution
> requires the upstream development environment. See `SOURCES.md` §"Reproducibility
> status" for full context.

## Atlas overview

The iHBCA v1.0 integrates **2,119,033 cells** from **287 donors** across
7 published single-cell RNA-seq studies of human breast tissue. All objects
are formatted to [CxG schema 5.3.2](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/5.3.0/schema.md)
and registered in the HCA Atlas Tracker.

## Source studies

| Study | Citation | Cells | Donors | Output file |
|-------|----------|------:|-------:|-------------|
| Gray | Gray et al. 2022, *Dev Cell* | 52,681 | 16 | `gray2022.h5ad` |
| Kumar | Kumar et al. 2023, *Nature* | 714,331 | 126 | `kumar2023.h5ad` |
| Murrow | Murrow et al. 2022, *Cell Systems* | 86,136 | 28 | `murrow2022.h5ad` |
| Nee | Nee et al. 2023, *Nat Genet* | 230,100 | 22 | `nee2023.h5ad` |
| Twigger | Twigger et al. 2022, *Nat Commun* | 110,744 | 16 | `twigger2022.h5ad` |
| Reed | Reed et al. 2024, *Nat Genet* | 803,283 | 55 | `reed2024.h5ad` |
| Pal | Pal et al. 2021, *EMBO J* | 131,288 | 22 | `pal2021.h5ad` |
| **Integrated** | — | **2,119,033** | **287** | `all-breast-cells.h5ad` |

## Repository structure

```
iHBCA_upload/
├── external_studies/              # Study ingestion and donor metadata harmonization
│   ├── config/                    # Study definitions
│   ├── preprocessing/             # Count matrix extraction from Seurat objects
│   └── harmonization/             # Cross-study metadata harmonization pipeline
│       ├── scripts/               # R and Python harmonization code
│       ├── run/                   # SLURM job wrappers
│       ├── config/                # Schemas and manifests
│       └── studies/               # Per-study config, reference tables, extraction notes
│
├── publication/                   # h5ad assembly and upload preparation
│   ├── scripts/                   # Assembly, enrichment, validation, entry sheet generation
│   ├── run/                       # SLURM job wrappers
│   ├── config/                    # Dataset registry and metadata configuration
│   ├── mappings/                  # Ontology mapping tables (HANCESTRO, EFO, CL, gene IDs, SRA)
│   └── tier2_breast/              # HCA Tier 2 breast-specific metadata proposals
│
└── provenance/                    # Original atlas construction code
    └── original_construction/     # Austin Marckx's iHBCA build and integration scripts
```

## Pipeline overview

### 1. External study ingestion (`external_studies/`)

Each source study's published data (Seurat RDS objects, supplemental tables)
is extracted into standardized count matrices, metadata CSVs, and UMAP
coordinates. Donor-level metadata is harmonized across all 7 studies into a
unified schema covering demographics, clinical covariates, and sample
processing details.

### 2. h5ad assembly (`publication/scripts/assemble_h5ad.py`)

Per-study h5ad files are assembled from extracted intermediates. The assembly
pipeline maps gene symbols to Ensembl IDs (GENCODE v24), populates CxG Tier 1
ontology fields (organism, assay, tissue, disease, sex, development stage,
ethnicity, cell type), and structures the AnnData object to schema requirements.

### 3. Integrated object assembly (`publication/scripts/assemble_integrated.py`)

The integrated object (`all-breast-cells.h5ad`, 2.12M cells) augments the
published iHBCA CELLxGENE deposit with refined annotations and enriched
metadata from the harmonization pipeline.

### 4. Enrichment and validation

`publication/scripts/enrich_h5ads.py` adds dataset-level metadata (`uns`),
embedding coordinates (`obsm`), and additional obs fields required by the
HCA tracker. Validation is performed with the
[CAP upload validator](https://pypi.org/project/cellxgene-census-upload-validator/)
and HCA-specific schema checks.

### 5. Entry sheet generation

`publication/scripts/generate_tier1_donor_sheet.py` and
`generate_tier1_sample_sheet.py` produce the donor- and sample-level
metadata sheets required for HCA Atlas Tracker registration.

## Original construction provenance

The `provenance/original_construction/` directory contains Austin Marckx's
original scripts for constructing the iHBCA:

- `10{a-f}_*_data_preparation.{R,py}` — per-study data preparation
- `merge_iHBCA.py` / `merge_iHBCA_simple.py` — atlas merging
- `iHBCA_CAP_object.py` / `iHBCA_CAP_meta_NEW.py` — CxG CAP object construction
- `integration/` — scVI, Harmony, scPoli integration runs
- `integration_metrics{1-3}/` — integration quality assessment

## Data access

- **CxG deposit**: [iHBCA collection on CELLxGENE](https://cellxgene.cziscience.com/collections/48259aa8-f168-4bf5-b797-af8e88da6637)
- **CellTypist models**: [Zenodo 10.5281/zenodo.10044650](https://doi.org/10.5281/zenodo.10044650)
- **HCA Data Portal**: Submission in progress

## Execution environment

Assembly and harmonization scripts were executed on the UCI HPC3 cluster
using Singularity containers. Shell scripts in `*/run/` directories are
SLURM job wrappers with resource specifications.

**Containers**:
- R 4.3.3 (Seurat v5, miloR): `Rocky8_jupyter_base_R4.3.3_Spatial.sif`
- Python (scanpy, anndata): `Jupyter_R_4.4.2_Giotto_Spatial_Python_2025Q2.sif`

## License

TBD

## Citation

TBD
