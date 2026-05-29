# step_17_integration_sweep — Multi-config integration benchmark

## What this stage does

Runs a benchmark of 9 integration configurations on the FLEX preprocessed
atlas (full object + 3 compartments × 9 configs = 36 sweep runs), scores
each against scIB metrics, and selects a single winner whose latent + UMAP
+ leiden sidecars become the canonical input for step_19b and step_20_finalize.

## Sweep architecture

Single orchestrator (`wrapper.py`) parameterized on:

- `--target {full, epi, str, imm, all}` — which atlas to integrate
- `--config-idx {0..8}` — which config in `source/configs/sweep_configs.yaml`

The SLURM array (`run/run_step_17_sweep.sh`) enumerates target × config_idx
pairs to fill the 36-run grid. Postprocess (`run/run_step_17_postprocess_winner.sh`)
selects the winner by composite scIB score and materializes the canonical
sidecars under `winner/`.

## Sweep config grid (9 configs)

Defined in `source/configs/sweep_configs.yaml`:

- scVI: n_latent ∈ {20, 50, 100, 200}
- Harmony (n_top_genes ∈ {2000, 5000})
- scANVI (n_latent = 50, with iHBCA L1 reference)
- BBKNN (default)
- ComBat (default)

The scoring formula in the yaml is pinned (see file header). Editing it
invalidates the canonical winner manifest.

## Winner sidecar contract

After postprocess, the winner directory contains:

```
${CFG_CANONICAL_SWEEP_WINNER}/
├── latent.csv              (cell_id × n_latent)
├── umap.csv                (cell_id × 2)
├── leiden_multi.csv        (cell_id × {leiden_0.5, leiden_1.0, leiden_2.0, ...})
├── metadata_with_umap.csv  (full obs table including the 2D coords)
└── compartments/
    ├── Epithelial/
    │   ├── latent.csv
    │   ├── umap.csv
    │   └── metadata_with_umap.csv
    ├── Stromal/   (same files)
    └── Immune/    (same files)
```

Step 19b consumes the per-compartment `metadata_with_umap.csv` files to
render the integration preview. Step 20_finalize consumes the full-object
`latent.csv` + `metadata_with_umap.csv` to assemble the canonical scvi_n100
bundle.

## Source modules

- `wrapper.py` — orchestrator, SLURM array index → (target, config_idx)
- `source/integration_sweep.py` — runs one (target, config) integration
- `source/scib_metric.py` — scIB metric computation (per config)
- `source/scib_aggregate.py` — composite score across configs
- `source/postprocess_winner.py` — winner selection + sidecar materialization
- `source/svd_diagnostic.py` — pre-integration SVD diagnostic (rank check)
- `source/recluster.py` — Leiden re-clustering at multiple resolutions
- `source/build_merged_counts.py` — pre-sweep merged count matrix builder
- `source/prep_kumar_labels.py` — Kumar reference label preparation
- `source/configs/sweep_configs.yaml` — sweep config grid (9 configs)

## Resource tiering

| Wrapper | Tier | Notes |
|---------|------|-------|
| `run_step_17_build_merged.sh` | 4 | Builds merged count matrix once |
| `run_step_17_sweep.sh` | 5 | 36-run array (target × config); per-task Tier 5 |
| `run_step_17_scib_array.sh` | 3 | Per-config scIB metric (36 tasks) |
| `run_step_17_postprocess_winner.sh` | 3 | Single-task winner selection |
| `run_svd_only.sh` | 3 | Standalone SVD diagnostic |
