#!/bin/bash
# Motifs pipeline internal-layout resolver (hybrid config, Option B).
#
# This script MUST be sourced AFTER publication/config/load_paths.sh, which
# exports SOURCE_*, CONTAINER_*, PYTHON_LIBS_*, BIND_MOUNTS, etc. (global
# resolution). This script only layers pipeline-internal output paths.
#
# Required envvars in scope when sourcing:
#   PIPELINE_ROOT — set by the wrapper to ${SOURCE_MOTIF_PIPELINE_ROOT}
#   SOURCE_MOTIF_PIPELINE_ROOT — exported by canonical load_paths.sh
#
# Usage in a wrapper:
#   source "${PUB_CONFIG}/load_paths.sh"                            # canonical
#   PIPELINE_ROOT="${SOURCE_MOTIF_PIPELINE_ROOT}"
#   source "${PIPELINE_ROOT}/config/load_paths.sh"                  # module-local

set -euo pipefail

: "${PIPELINE_ROOT:?PIPELINE_ROOT must be set before sourcing motifs/config/load_paths.sh}"
: "${SOURCE_MOTIF_PIPELINE_ROOT:?Source publication/config/load_paths.sh first to populate SOURCE_*}"

# Module-internal output layout. Matches publication/analysis/spatial/motifs/
# convention: single outputs/ root with named subdirs where granularity helps.
export MOTIF_OUTPUTS_ROOT="${PIPELINE_ROOT}/outputs"
export MOTIF_VALIDATION_DIR="${PIPELINE_ROOT}/outputs/validation"
export MOTIF_UNITS_DIR="${PIPELINE_ROOT}/outputs/motif_units"

# SLURM logs land in the canonical publication logs dir (HPC_LOGS_DIR is
# exported by the canonical resolver). Wrappers mkdir this directory at
# job start; we don't mkdir here so this file is safe to source from a
# read-only local CRSP mount.
export LOGS_ROOT="${HPC_LOGS_DIR}"
