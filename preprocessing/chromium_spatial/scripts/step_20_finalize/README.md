# step_20_finalize — Canonical scvi_n100 bundle producer

## What this stage does

Materializes the canonical `scvi_n100` bundle for the FLEX preprocessing
atlas — the handoff artifact consumed by all downstream tracks. Three
scripts run in sequence:

1. `20a_concat_full_object.py` — concatenates per-compartment counts +
   obs into a full-object `integration_intermediate/full/scvi_n100/`
   directory (CSV + MTX format, no h5ad).
2. `20b_promote_full_scvi_n100.py` — extracts latent + UMAP + obs from the
   step 17 sweep winner h5ad and writes the minimal schema consumed by
   Milo build. Counts matrix is not re-written here; 20a handles that.
3. `20_build_manifest.py` — walks the assembled directory and writes
   `outputs/manifest.yaml` recording every file. Fails loudly on missing
   expected files.

## Bundle contents (per compartment + full object)

```
integration_intermediate/{Epithelial,Immune,Stromal,full}/scvi_n100/
├── latent.csv       (cell_id × n_latent=100)
├── umap.csv         (cell_id × 2)
├── obs.csv          (cell_id × metadata columns)
├── counts.mtx.gz    (sparse counts, genes × cells)
├── genes.tsv        (gene_id, gene_symbol per row of counts)
├── cells.tsv        (cell_id per column of counts)
└── manifest.json    (per-bundle integrity record)
```

Plus `outputs/manifest.yaml` at the directory root recording all bundle
files with relative paths.

## Cell counts

- Per-compartment sum: 256,858 cells (from the 5-way QC filter, partitioned
  into Epithelial/Immune/Stromal)
- Full object: 258,317 cells (from the 5-way QC filter, including ~1,459
  cells dropped by the per-compartment Seurat assemblies)

The full-object bundle is built by 20b directly from the sweep winner
h5ad to preserve the 258K count; 20a's concat path produces only the
256,858-cell partitioned view.

## Downstream consumers

| Track | What it reads | Why |
|-------|---------------|-----|
| Track B (joint preprocessing) | full + per-compartment latent + counts | Joint Concord 50D integration with Xenium |
| Track C (FLEX annotation) | full latent + obs + counts | Per-cell L1.5 + L2 label production |
| Track D (FLEX-only DA) | full latent + obs (parameterized on label_set from C) | Milo DA per contrast |
| `publication/analysis/abundance/spatial/*` | latent + obs + labels | Existing DA pipelines |

## Output stability contract

The bundle file names and column conventions above are the canonical
handoff. Downstream tracks pin against these names. Changes to the bundle
schema invalidate Track B / C / D consumers.

## Comments on h5ad

Per project policy (output-storage rule), no h5ad files are emitted as
canonical pipeline output. The bundle is CSV + MTX. Downstream tracks
that need a working AnnData object reconstruct it at runtime from the
bundle.
