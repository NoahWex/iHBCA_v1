#!/usr/bin/env bash
# _common.sh — Shared config loader for FLEX L2S annotation pipeline wrappers.
#
# Sources publication load_paths.sh to obtain USER_ROOT, SOURCE_*, CONTAINER_*,
# R_LIBS_*, BIND_MOUNTS, and HPC_LOGS_DIR, then defines a small set of
# CFG_* tokens that the wrappers reference for the annotation tree.
#
# Wrappers must set CONTAINER_TYPE before sourcing this file:
#   CONTAINER_TYPE=python_spatial_2025Q2 | r_spatial_4.3.3
#
# Usage from a wrapper:
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/_common.sh"

set -euo pipefail

# ============================================================
# 1. Resolve publication-tier roots via load_paths.sh
# ============================================================
# Walk up from the current run/ directory looking for
# publication/config/load_paths.sh. Works from either the staging tree
# (coordination/staging/track_c/run/) or the post-promotion tree
# (publication/preprocessing/chromium_spatial/annotation/run/).

COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANNOTATION_ROOT_LOCAL="$(cd "$COMMON_DIR/.." && pwd)"

if [ -z "${USER_ROOT:-}" ] || [ -z "${BIND_MOUNTS:-}" ]; then
    _SEARCH_DIR="${ANNOTATION_ROOT_LOCAL}"
    PUBLICATION_CONFIG_DIR=""
    while [ "$_SEARCH_DIR" != "/" ]; do
        if [ -f "${_SEARCH_DIR}/publication/config/load_paths.sh" ]; then
            PUBLICATION_CONFIG_DIR="${_SEARCH_DIR}/publication/config"
            break
        fi
        _SEARCH_DIR="$(dirname "$_SEARCH_DIR")"
    done
    if [ -z "$PUBLICATION_CONFIG_DIR" ]; then
        echo "FATAL: could not locate publication/config/load_paths.sh by walking up from $ANNOTATION_ROOT_LOCAL" >&2
        exit 2
    fi

    # publication load_paths.sh:23 needs python3 with PyYAML on compute nodes.
    module load anaconda/2024.06 2>/dev/null || module load python/3.10.0 2>/dev/null || true

    # shellcheck source=/dev/null
    source "${PUBLICATION_CONFIG_DIR}/load_paths.sh"
fi

# ============================================================
# 2. Annotation-tree CFG_* tokens
# ============================================================
# These resolve relative to the running tree (staging or promoted), not via
# paths.yaml. The pipeline self-contains its scripts/ and outputs/ subtree.

CFG_ANNOTATION_ROOT="${ANNOTATION_ROOT_LOCAL}"
CFG_SCRIPTS_DIR="${CFG_ANNOTATION_ROOT}/scripts"
CFG_OUTPUTS_ROOT="${CFG_ANNOTATION_ROOT}/outputs"
CFG_ANNOTATION_YAMLS="${CFG_ANNOTATION_ROOT}/annotation_yamls"

# Override for SLURM execution: CFG_ANNOTATION_ROOT may differ from where the
# script literally lives if the wrapper was submitted with --chdir to a
# different staging path. Honor it when set explicitly.
if [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${CFG_PROJECT_ROOT:-}" ] && [ "${CFG_PROJECT_ROOT}" != "${CFG_ANNOTATION_ROOT}" ]; then
    CFG_ANNOTATION_ROOT="${CFG_PROJECT_ROOT}"
    CFG_SCRIPTS_DIR="${CFG_ANNOTATION_ROOT}/scripts"
    CFG_OUTPUTS_ROOT="${CFG_ANNOTATION_ROOT}/outputs"
    CFG_ANNOTATION_YAMLS="${CFG_ANNOTATION_ROOT}/annotation_yamls"
fi

# ============================================================
# 3. Upstream tokens from publication paths.yaml
# ============================================================
# FLEX integration intermediate is the per-compartment scvi_n100 bundle
# (counts.mtx.gz, latent.csv, umap.csv, obs.csv) that this pipeline reads.
# Resolved by load_paths.sh as SOURCE_FLEX_INTEGRATION_INTERMEDIATE (pending
# Noah approval of the Track 3 paths.yaml batch). Falls back to
# ${PROJECT_SPATIAL_HBCA_PREPROCESSING}/outputs/integration_intermediate while
# Track 3 is unpromoted.

if [ -n "${SOURCE_FLEX_INTEGRATION_INTERMEDIATE:-}" ]; then
    CFG_INTEGRATION_INTERMEDIATE="${SOURCE_FLEX_INTEGRATION_INTERMEDIATE}"
else
    CFG_INTEGRATION_INTERMEDIATE="${USER_ROOT}/Spatial_HBCA_preprocessing/outputs/integration_intermediate"
fi

# V1 annotation yamls (Track A promoted Apr 26).
CFG_V1_ANNOTATION_YAMLS="${PROJECT_PUBLICATION}/publication/analysis/annotation/yamls"

# ============================================================
# 4. Container resolution (mirrors publication CONTAINER_* exports)
# ============================================================

if [ -z "${CONTAINER_TYPE:-}" ]; then
    echo "FATAL: CONTAINER_TYPE not set before sourcing _common.sh" >&2
    echo "  Set CONTAINER_TYPE=python_spatial_2025Q2 | r_spatial_4.3.3" >&2
    exit 2
fi

case "${CONTAINER_TYPE}" in
    python_spatial_2025Q2)
        CONTAINER_PATH="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
        ;;
    r_spatial_4.3.3|r_spatial_4_3_3)
        CONTAINER_PATH="${CONTAINER_R_SPATIAL_4_3_3}"
        ;;
    *)
        echo "FATAL: unknown CONTAINER_TYPE='${CONTAINER_TYPE}'" >&2
        exit 2
        ;;
esac

# ============================================================
# 5. Singularity bind mounts
# ============================================================
# BIND_MOUNTS is exported by load_paths.sh as a pre-formatted "--bind ..." string.

# ============================================================
# 6. Per-job temp directory with cleanup
# ============================================================

TEMP_ROOT="${USER_ROOT}/tmp"
mkdir -p "${TEMP_ROOT}"
JOB_TEMP_DIR=$(mktemp -d "${TEMP_ROOT}/trackc_${SLURM_JOB_ID:-local}_XXXXXX")
trap 'rm -rf "$JOB_TEMP_DIR"' EXIT

# ============================================================
# 7. Helper: load singularity module
# ============================================================

_load_singularity() {
    module load singularity 2>/dev/null || true
}

# ============================================================
# 8. Startup banner
# ============================================================

echo "=== FLEX L2S Annotation Pipeline ==="
echo "Annotation root: ${CFG_ANNOTATION_ROOT}"
echo "Container:       ${CONTAINER_TYPE} (${CONTAINER_PATH})"
echo "Job ID:          ${SLURM_JOB_ID:-local}"
echo "Task ID:         ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Temp:            ${JOB_TEMP_DIR}"
echo "Start:           $(date)"
echo "====================================="
