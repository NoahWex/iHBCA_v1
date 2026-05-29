# External Studies

Ingestion and harmonization of the 7 published scRNA-seq studies that comprise
the iHBCA v1.0.

## Subdirectories

| Directory | Purpose |
|-----------|---------|
| `config/` | Study manifest (`studies.yaml`) defining all 7 source studies |
| `preprocessing/` | Count matrix extraction from Seurat RDS objects |
| `harmonization/` | Cross-study donor metadata harmonization pipeline |

## Data flow

```
Published data (Seurat RDS, supplemental tables)
        |
        v
  preprocessing/        -> count matrices (MTX), feature lists, barcodes
        |
        v
  harmonization/        -> standardized donor metadata, cell ID mappings,
                           UMAP coordinates, per-study published/ directories
        |
        v
  outputs/ (gitignored) -> per-study published/ dirs consumed by
                           publication/scripts/assemble_h5ad.py
```

## Studies

| Study | Cells | Donors | Citation |
|-------|------:|-------:|----------|
| Gray | 52,681 | 16 | Gray et al. 2022, *Dev Cell* |
| Kumar | 714,331 | 126 | Kumar et al. 2023, *Nature* |
| Murrow | 86,136 | 28 | Murrow et al. 2022, *Cell Systems* |
| Nee | 230,100 | 22 | Nee et al. 2023, *Nat Genet* |
| Twigger | 110,744 | 16 | Twigger et al. 2022, *Nat Commun* |
| Reed | 803,283 | 55 | Reed et al. 2024, *Nat Genet* |
| Pal | 131,288 | 22 | Pal et al. 2021, *EMBO J* |
