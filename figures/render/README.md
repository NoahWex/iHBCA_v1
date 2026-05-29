# Figure render scripts

## Naming conventions

- `render_*.{py,R}` — produces a figure PDF (single panel or composite)
- `compose_*.py` — assembles multi-panel composites from standalone panel PDFs (typically via pikepdf / pypdf form-XObject preservation)
- `run_*.sh` — SLURM wrapper that sources `publication/config/load_paths.sh` and invokes the corresponding render script with paths injected from `publication/config/paths.yaml`

## Subdirectory layout (by figure)

Render scripts live under `render/figN/` (where N = 1, 2, 3, 4) co-located with their SLURM wrappers. This layout was adopted 2026-05-22 per the Fig 1 promotion session decision and supersedes the source-based layout (`ihbca/`, `flex/`, `joint/`, `xenium/`).

- `fig1/` — iHBCA atlas overview panels + supps (1b, 1c, s1.1–s1.7); see `submission/figures/fig1/INDEX.md` for the panel catalog
- `fig2/` — risk-associated cell-state changes (will populate as the Fig 2 promotion session migrates its scripts)
- `fig3/` — anatomically resolved spatial atlas (will populate as the Fig 3 promotion session migrates its scripts)
- `fig4/` — synthesis schematic (Illustrator only; no render scripts expected)

## Wrapper co-location

SLURM wrappers (`run_*.sh`) live alongside their render scripts under `render/figN/`, not in a separate `figures/run/` tree. One panel = one directory with both producer and wrapper.

`publication/figures/run/` remains as a holding pen for un-migrated wrappers during the transition and will empty as figure-owner sessions migrate. (Note: `publication/run/` — top-level — is the SHARED SLURM utilities dir and is distinct from `figures/run/`.)

## Legacy platform-based subdirectories

The source-based subdirs are deprecating, not deleted. They remain populated with un-migrated scripts until their figure-owner sessions migrate them:

- `flex/` — Chromium FLEX preprocessing + analysis renders (mostly Fig 3 supps)
- `ihbca/` — iHBCA atlas renders not yet attributed to a figure session (mostly Fig 2-attributed remaining)
- `joint/` — Joint FLEX + Xenium integration renders (Figs 3 + 4 + joint-aware supps)
- `xenium/` — Xenium-specific renders

When a session promotes a script under the figure-based convention, that script moves to `render/figN/` and the source-based copy is removed. The source-based dirs will be deleted once empty.

## Substrate-prep scripts live in `data_prep/`, not here

Scripts that build/extract/compute intermediate substrate (CSVs, parquets, NPZs consumed by render scripts) live in `publication/figures/data_prep/figN/` — a sibling directory mirroring this one's figure-based layout. Naming:

- `build_*` — joins/composes substrate from multiple upstream sources
- `compute_*` — derives a single numerical quantity (e.g. mean expression per cluster)
- `extract_*` — pulls a slice from a canonical upstream object

These were separated from render scripts on 2026-05-17 to clarify the data-flow boundary: a render script should consume publication-ready substrate and produce a publication-ready PDF; substrate construction is a distinct stage.

## Superseded scripts: `_archived/`

Scripts no longer in use go to `_archived/` with a note on what replaced them. Currently:
- `_archived/render_s_ihbca_confusion_matrix_L2.py` — superseded 2026-05-16 by the unified L1+L2 variant (now at `fig1/render_s_ihbca_confusion_matrix.py` post Atom-1 migration; broken substrate join in L2 variant produced 20-40% L2→L1 purity)

## Per-platform notes

- `joint/.archived_slurm_logs/` (if present) holds historical SLURM stdout files from 2026-05-11 → 2026-05-14 render iterations; preserved for provenance trail but not actively read.
