# publication/analysis/

Analysis modules promoted from development repos. Each module produces
canonical, manifest-tracked artifacts for the manuscript.

## Modules

| Module | Purpose |
|---|---|
| `abundance/` | Differential abundance testing (iHBCA V1 clinical + spatial FLEX/joint Milo) |
| `annotation/` | Label-system development (L1 harmonization, L2.0c/L2.0s formalization, Xenium L1.5) |
| `communication/` | Cell-cell communication (CellChat, CellCharter, L-R co-localization) |
| `markers/` | Within-L1 and all-vs-all marker tables |
| `spatial/` | Spatial-specific analyses (cell-cell, niche, motif decomposition) |

## Module layout convention

Every analysis module under `publication/analysis/` adopts this structure:

```
{module}/
├── scripts/         executable analysis code (.py / .R / .qmd)
├── run/             SLURM wrappers (sourced from canonical resolver first)
├── config/          module-local config (palette, hyperparameters, internal layout)
├── outputs/         all produced artifacts (single root; subdirs for organization)
└── README.md        module-specific entry-point (optional)
```

The `outputs/` consolidation is intentional: every produced artifact lives
under one root, with named subdirectories where granularity helps (e.g.
`outputs/validation/`, `outputs/motif_units/`). This replaces split
dev-pipeline layouts like `data/nmf/` + `data/validation/` + `reports/`,
which forced consumers to remember three location conventions.

## Hybrid config pattern (Option B)

Modules that need filesystem path resolution use a **hybrid resolver**:

1. **Canonical resolver** (`publication/config/load_paths.{sh,py,R}`)
   handles global discovery: HPC/local roots, containers, bind mounts,
   R/Python library paths, cross-module substrate tokens (e.g.
   `${SOURCE_MOTIF_BASIS_CSV}`, `${SOURCE_XENIUM_JOINT_L1P5}`).

2. **Module-local resolver** (`{module}/config/load_paths.{sh,py}`,
   optional) handles only internal output layout: where this module
   writes its artifacts within its own tree. Sources canonical first to
   inherit `SOURCE_*`, `CONTAINER_*`, `BIND_MOUNTS`, then layers internal
   exports.

A SLURM wrapper sources both, in order:

```bash
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"   # canonical
PIPELINE_ROOT="${SOURCE_<MODULE>_PIPELINE_ROOT}"
source "${PIPELINE_ROOT}/config/load_paths.sh"            # module-local
```

This avoids forking the resolved values (one source of truth for what
`/share/crsp/...` is), while letting each module model its own internal
output layout without forcing the canonical config to track every
subdir.

**No module duplicates the canonical `roots:` block.** Modules that
don't need a local resolver omit the `config/load_paths.*` files entirely
and rely only on canonical.

## Promotion target

All scripts in `publication/analysis/` are promoted from dev repos via
the protocol in `iHBCA_publication/CLAUDE.md` (§ Promotion Protocol).
Every file is logged in `publication/manifest.yaml` with source path,
validation gate, and validation result. Two automated gates enforce
parameterization + reviewer-facing comments (see
`.claude/rules/promotion-enforcement.md`).
