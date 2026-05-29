#!/bin/bash
# _common.sh — Shared config loader for FLEX preprocessing pipeline wrappers
#
# Usage: source "$(dirname "$0")/_common.sh"
#
# Sets up: paths from config, container environment, bind mounts,
#          temp directories, DRY_RUN handling, helper functions.
#
# Wrappers must set CONTAINER_TYPE before sourcing this file:
#   CONTAINER_TYPE=scgpt_gpu | giotto_spatial | r_spatial
#
# This file sources publication-tier load_paths.sh first to obtain
# USER_ROOT and SHARED_DATASETS (resolved per environment), then expands
# the local config/paths.yaml token namespace via _expand_refs.

set -euo pipefail

# ============================================================
# 1. Resolve publication-tier roots via load_paths.sh
# ============================================================
# load_paths.sh exports USER_ROOT, SHARED_DATASETS, HPC_HOME, BIND_MOUNTS,
# BIND_MOUNT_LIST, HPC_LOGS_DIR, plus PROJECT_*, SOURCE_*, CONTAINER_*.
#
# Note: SLURM #SBATCH --output directives in the wrappers reference
# ${HPC_LOGS_DIR}; the user must source load_paths.sh BEFORE invoking
# sbatch for those substitutions to resolve. This _common.sh inherits
# already-exported vars from the calling shell.

COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT_LOCAL="$(cd "$COMMON_DIR/.." && pwd)"
CONFIG_DIR="${PROJECT_ROOT_LOCAL}/config"

# Source the publication path resolver (only if not already loaded — under
# SLURM execution USER_ROOT and BIND_MOUNTS are typically unset because
# `hpc submit` does not pre-source publication/config/load_paths.sh in the
# remote shell).
#
# _common.sh may be sourced from either staging
# (coordination/staging/track_3/run/) or post-promotion
# (publication/preprocessing/chromium_spatial/run/) — different depths from
# the iHBCA_publication root. Walk up looking for
# publication/config/load_paths.sh rather than hardcoding a relative depth.
# Fails loud if outside any iHBCA_publication/ tree.
if [ -z "${USER_ROOT:-}" ] || [ -z "${BIND_MOUNTS:-}" ]; then
    _SEARCH_DIR="${PROJECT_ROOT_LOCAL}"
    PUBLICATION_CONFIG_DIR=""
    while [ "$_SEARCH_DIR" != "/" ]; do
        if [ -f "${_SEARCH_DIR}/publication/config/load_paths.sh" ]; then
            PUBLICATION_CONFIG_DIR="${_SEARCH_DIR}/publication/config"
            break
        fi
        _SEARCH_DIR="$(dirname "$_SEARCH_DIR")"
    done
    if [ -z "$PUBLICATION_CONFIG_DIR" ]; then
        echo "FATAL: could not locate publication/config/load_paths.sh by walking up from $PROJECT_ROOT_LOCAL" >&2
        exit 2
    fi

    # Ensure python has yaml available — load_paths.sh:40 invokes a python -c
    # block that imports yaml. Default compute-node python3 lacks the module.
    # This dependency belongs in load_paths.sh (which owns the python invocation),
    # but is patched here as a Track A staging-level workaround pending a
    # publication-side fix authorized by the project lead.
    module load anaconda/2024.06 2>/dev/null || module load python/3.10.0 2>/dev/null || true

    # shellcheck source=/dev/null
    source "${PUBLICATION_CONFIG_DIR}/load_paths.sh"
fi

# Bridge publication-tier exports into _expand_refs's CFG_* namespace.
# _expand_refs resolves ${user_root} by looking up CFG_USER_ROOT, etc.
CFG_USER_ROOT="${USER_ROOT}"
CFG_SHARED_DATASETS="${SHARED_DATASETS}"

# Verify local config files exist
for cfg in paths.yaml containers.yaml methods.yaml; do
    if [ ! -f "${CONFIG_DIR}/${cfg}" ]; then
        echo "FATAL: Config not found: ${CONFIG_DIR}/${cfg}" >&2
        exit 2
    fi
done

# ============================================================
# 2. Parse paths.yaml (key extraction + ${var} expansion)
# ============================================================

_yaml_get() {
    # Extract a top-level key from a YAML file (no nesting support)
    local file="$1" key="$2"
    grep -E "^${key}:" "$file" | head -1 | sed "s/^${key}:[[:space:]]*//"
}

PATHS_YAML="${CONFIG_DIR}/paths.yaml"

# Expand ${var} references in a value using already-loaded CFG_* vars.
# Fails loud if a ref cannot be resolved.
_expand_refs() {
    local val="$1"
    while [[ "$val" =~ \$\{([a-z_0-9]+)\} ]]; do
        local key="${BASH_REMATCH[1]}"
        # Portable uppercase (bash 3.2 compatible — ${x^^} is bash 4+)
        local upper
        upper="CFG_$(printf '%s' "$key" | tr '[:lower:]' '[:upper:]')"
        local replacement="${!upper:-}"
        if [ -z "$replacement" ]; then
            echo "FATAL: unresolved ref \${${key}} in paths.yaml (tried $upper)" >&2
            exit 2
        fi
        val="${val//\$\{${key}\}/${replacement}}"
    done
    printf '%s\n' "$val"
}

# Project root must be resolved first — every other token references it
CFG_PROJECT_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" project_root)")"
CFG_OUTPUTS_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" outputs_root)")"

# Upstream (immutable inputs)
CFG_UPSTREAM_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" upstream_root)")"
CFG_UPSTREAM_MERGED_COUNTS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" upstream_merged_counts)")"
CFG_UPSTREAM_CELL_METADATA="$(_expand_refs "$(_yaml_get "$PATHS_YAML" upstream_cell_metadata)")"
CFG_UPSTREAM_CONTAMINATION="$(_expand_refs "$(_yaml_get "$PATHS_YAML" upstream_contamination)")"
CFG_UPSTREAM_SINGLER_LABELS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" upstream_singler_labels)")"
CFG_UPSTREAM_PASSING_CELLS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" upstream_passing_cells)")"
CFG_COMPARTMENT_CELL_LISTS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" compartment_cell_lists)")"
CFG_SAMPLE_LIST="$(_expand_refs "$(_yaml_get "$PATHS_YAML" sample_list)")"
CFG_SCRIPTS_DIR="$(_expand_refs "$(_yaml_get "$PATHS_YAML" scripts_dir)")"
CFG_LOGS_DIR="$(_expand_refs "$(_yaml_get "$PATHS_YAML" logs_dir)")"
CFG_TEMP_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" temp_root)")"
CFG_ACCOUNT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" account)")"
CFG_L1_PSEUDOBULK="$(_expand_refs "$(_yaml_get "$PATHS_YAML" l1_pseudobulk_csv)")"

# Override CFG_PROJECT_ROOT (and dependent tokens) when the wrapper is
# running from a location that differs from paths.yaml's canonical
# project_root. Under SLURM, PIPELINE_ROOT (derived from pwd at the
# wrapper's `--chdir` target) is authoritative for where the wrapper
# actually lives — use it to keep CFG_SCRIPTS_DIR pointed at the correct
# on-disk scripts/ subtree. This handles staging-tree dry runs (where
# PIPELINE_ROOT = coordination/staging/track_3 but paths.yaml resolves
# project_root to publication/preprocessing/chromium_spatial). Post-
# promotion the two values match and this override is a no-op.
if [ -n "${SLURM_JOB_ID:-}" ] && [ "${PIPELINE_ROOT}" != "${CFG_PROJECT_ROOT}" ]; then
    CFG_PROJECT_ROOT="${PIPELINE_ROOT}"
    CFG_SCRIPTS_DIR="${PIPELINE_ROOT}/scripts"
fi

# L1 reference build inputs (legacy + current)
CFG_IHBCA_AUTHOR_SHARE_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_author_share_root)")"
CFG_IHBCA_COUNTS_NPZ="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_counts_npz)")"
CFG_IHBCA_LEVEL15_ANNOTATIONS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_level15_annotations)")"
CFG_IHBCA_GENE_DATA="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_gene_data)")"
CFG_GENE_SYMBOL_TO_ENSEMBL="$(_expand_refs "$(_yaml_get "$PATHS_YAML" gene_symbol_to_ensembl)")"
CFG_IHBCA_ASSEMBLY_COMPONENTS_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_assembly_components_root)")"
CFG_IHBCA_H5_EPI="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_h5_epi)")"
CFG_IHBCA_H5_STR="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_h5_str)")"
CFG_IHBCA_H5_IMM="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_h5_imm)")"
CFG_IHBCA_ASSEMBLY_CELL_METADATA="$(_expand_refs "$(_yaml_get "$PATHS_YAML" ihbca_assembly_cell_metadata)")"

# Canonical preprocessing root (referenced by per-step tokens)
CFG_PREPROCESSING_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_root)")"
CFG_MODULE_CONFIGS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" module_configs)")"

# L1 per-step (steps 01-07)
CFG_PREPROCESSING_STEP_01="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_01)")"
CFG_PREPROCESSING_STEP_02="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_02)")"
CFG_PREPROCESSING_STEP_03="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_03)")"
CFG_PREPROCESSING_STEP_04="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_04)")"
CFG_PREPROCESSING_STEP_05="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_05)")"
CFG_PREPROCESSING_STEP_06="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_06)")"
CFG_PREPROCESSING_STEP_07="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_07)")"

# L1 aggregation (gates everything downstream)
CFG_PREPROCESSING_QC_STATUS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_qc_status)")"
CFG_PREPROCESSING_CENTRAL_CELL_STATUS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_central_cell_status)")"

# L2 per-step (steps 08-13)
CFG_PREPROCESSING_STEP_08="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_08)")"
CFG_PREPROCESSING_STEP_09="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_09)")"
CFG_PREPROCESSING_STEP_10A="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_10a)")"
CFG_PREPROCESSING_STEP_10B="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_10b)")"
CFG_PREPROCESSING_STEP_11="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_11)")"
CFG_PREPROCESSING_STEP_12="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_12)")"
CFG_PREPROCESSING_STEP_13="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_13)")"

# L3 per-step (steps 14-19b)
CFG_PREPROCESSING_STEP_14="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_14)")"
CFG_PREPROCESSING_STEP_15="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_15)")"
CFG_PREPROCESSING_STEP_16="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_16)")"
CFG_PREPROCESSING_STEP_17="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_17)")"
CFG_PREPROCESSING_STEP_18_DESEQ="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_18_deseq)")"
CFG_PREPROCESSING_STEP_19B="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_19b)")"
CFG_INTEGRATION_INTERMEDIATE="$(_expand_refs "$(_yaml_get "$PATHS_YAML" integration_intermediate)")"

# Frozen canonical handoff artifacts
CFG_PREPROCESSING_CELL_RETENTION_LIST="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_cell_retention_list)")"
CFG_PREPROCESSING_STEP_08_MODEL="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_08_model)")"
CFG_PREPROCESSING_STEP_12_MODEL="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step_12_model)")"
CFG_PREPROCESSING_CELL_COMPARTMENTS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_cell_compartments)")"
CFG_PREPROCESSING_EPIDERMAL_SAMPLES="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_epidermal_samples)")"
CFG_PREPROCESSING_CONTAMINATION_CELLS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_contamination_cells)")"
CFG_PREPROCESSING_DA_TESTING_RESULTS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_da_testing_results)")"
CFG_PREPROCESSING_STEP16_VF_DIR="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step16_vf_dir)")"
CFG_PREPROCESSING_STEP16_RETENTION="$(_expand_refs "$(_yaml_get "$PATHS_YAML" preprocessing_step16_retention)")"

# Canonical sweep root (references step_17)
CFG_CANONICAL_SWEEP_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_sweep_root)")"
CFG_CANONICAL_SWEEP_FULL_OBJECT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_sweep_full_object)")"
CFG_CANONICAL_SWEEP_COMPARTMENTS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_sweep_compartments)")"
CFG_CANONICAL_SWEEP_SCIB="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_sweep_scib)")"
CFG_CANONICAL_SWEEP_WINNER="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_sweep_winner)")"
CFG_CANONICAL_SWEEP_SVD="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_sweep_svd)")"

# Final atlas sidecars
CFG_CANONICAL_LATENT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_latent)")"
CFG_CANONICAL_UMAP="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_umap)")"
CFG_CANONICAL_LEIDEN="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_leiden)")"
CFG_CANONICAL_METADATA="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_metadata)")"
CFG_CANONICAL_SCVI_MODEL="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_scvi_model)")"
CFG_CANONICAL_L1_LABELS="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_l1_labels)")"
CFG_CANONICAL_MANIFEST="$(_expand_refs "$(_yaml_get "$PATHS_YAML" canonical_manifest)")"

# Raw data
CFG_RAW_DATA_ROOT="$(_expand_refs "$(_yaml_get "$PATHS_YAML" raw_data_root)")"
CFG_RAW_DATA_SAMPLE_MANIFEST="$(_expand_refs "$(_yaml_get "$PATHS_YAML" raw_data_sample_manifest)")"

# ============================================================
# 3. Parse containers.yaml for selected container
# ============================================================

CONTAINERS_YAML="${CONFIG_DIR}/containers.yaml"

if [ -z "${CONTAINER_TYPE:-}" ]; then
    echo "FATAL: CONTAINER_TYPE not set before sourcing _common.sh" >&2
    echo "  Set CONTAINER_TYPE=scgpt_gpu|giotto_spatial|r_spatial" >&2
    exit 2
fi

_yaml_get_nested() {
    # Extract a value under a specific container block
    local file="$1" container="$2" key="$3"
    awk "/^  ${container}:/{found=1; next} /^  [a-z]/{if(found) exit} found && /^    ${key}:/{gsub(/^    ${key}:[[:space:]]*/,\"\"); print; exit}" "$file"
}

CONTAINER_PATH="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" path)"

if [ -z "$CONTAINER_PATH" ]; then
    echo "FATAL: Container type '${CONTAINER_TYPE}' not found in containers.yaml" >&2
    exit 2
fi

# Container-specific lib paths
R_LIBS_USER_HOST="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" r_libs_user_host)"
R_LIBS_USER_CONTAINER="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" r_libs_user_container)"
PYTHONUSERBASE_HOST="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" pythonuserbase_host)"
PYTHONUSERBASE_CONTAINER="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" pythonuserbase_container)"
PANDOC_PATH="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" pandoc_path)"
CONTAINER_GPU="$(_yaml_get_nested "$CONTAINERS_YAML" "$CONTAINER_TYPE" gpu)"

# ============================================================
# 4. Singularity bind mounts (from publication-tier load_paths.sh)
# ============================================================
# BIND_MOUNT_LIST is a space-separated list of host:container:mode triples.
# Build BIND_FLAGS as repeated --bind args for singularity exec.

BIND_FLAGS=""
for mount in ${BIND_MOUNT_LIST}; do
    BIND_FLAGS="${BIND_FLAGS} --bind ${mount}"
done

# ============================================================
# 5. Per-job temp directory with cleanup
# ============================================================

mkdir -p "${CFG_TEMP_ROOT}"
JOB_TEMP_DIR=$(mktemp -d "${CFG_TEMP_ROOT}/job_${SLURM_JOB_ID:-local}_XXXXXX")
trap 'rm -rf "$JOB_TEMP_DIR"' EXIT

# ============================================================
# 6. DRY_RUN handling
# ============================================================
# Wrappers pass DRY_RUN=1 to enable dry-run mode.
# Scripts receive --dry-run flag.

DRY_RUN="${DRY_RUN:-0}"
DRY_RUN_FLAG=""
if [ "$DRY_RUN" = "1" ]; then
    DRY_RUN_FLAG="--dry-run"
    echo "=== DRY RUN MODE ==="
fi

# ============================================================
# 7. Helper functions
# ============================================================

_load_singularity() {
    module load singularity 2>/dev/null || true
}

run_python_singularity() {
    local script="$1"; shift
    _load_singularity

    local nv_flag=""
    if [ "${CONTAINER_GPU}" = "true" ]; then
        nv_flag="--nv"
    fi

    local env_flags=""
    if [ -n "${PYTHONUSERBASE_HOST}" ]; then
        env_flags="--env PYTHONUSERBASE=${PYTHONUSERBASE_CONTAINER}"
    fi

    # Per-task numba/matplotlib caches under JOB_TEMP_DIR avoid concurrent
    # array tasks racing on shared /tmp/numba_cache (numba's FunctionCache
    # locator dies on contended writes).
    mkdir -p "${JOB_TEMP_DIR}/numba_cache" "${JOB_TEMP_DIR}/mpl_config"

    singularity exec ${nv_flag} \
        ${BIND_FLAGS} \
        --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
        --env "NUMBA_CACHE_DIR=${JOB_TEMP_DIR}/numba_cache" \
        --env "MPLCONFIGDIR=${JOB_TEMP_DIR}/mpl_config" \
        ${env_flags} \
        "${CONTAINER_PATH}" \
        python3 "${script}" "$@"
}

run_r_singularity() {
    local script="$1"; shift
    _load_singularity

    local r_bind=""
    local r_env=""
    if [ -n "${R_LIBS_USER_HOST}" ]; then
        r_bind="--bind ${R_LIBS_USER_HOST}:${R_LIBS_USER_CONTAINER}:ro"
        r_env="--env R_LIBS_USER=${R_LIBS_USER_CONTAINER}"
    fi

    local pandoc_env=""
    if [ -n "${PANDOC_PATH}" ]; then
        pandoc_env="--env RSTUDIO_PANDOC=${PANDOC_PATH}"
    fi

    singularity exec \
        --no-mount bind-paths \
        ${BIND_FLAGS} \
        ${r_bind} \
        --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
        ${r_env} \
        ${pandoc_env} \
        --env "TMPDIR=${JOB_TEMP_DIR}" \
        --env "TMP=${JOB_TEMP_DIR}" \
        --env "TEMP=${JOB_TEMP_DIR}" \
        --env "LANG=en_US.UTF-8" \
        --env "LC_ALL=en_US.UTF-8" \
        "${CONTAINER_PATH}" \
        Rscript "${script}" "$@"
}

run_rmd_render() {
    local rmd="$1"
    local output_file="$2"
    shift 2
    # Remaining args are param=value pairs
    local params_r="list("
    local first=true
    for param in "$@"; do
        local key="${param%%=*}"
        local val="${param#*=}"
        if [ "$first" = true ]; then
            first=false
        else
            params_r="${params_r}, "
        fi
        params_r="${params_r}${key} = '${val}'"
    done
    params_r="${params_r})"

    run_r_singularity -e "rmarkdown::render('${rmd}', params = ${params_r}, output_file = '${output_file}', knit_root_dir = '$(dirname "$rmd")')"
}

# ============================================================
# 8. Startup banner
# ============================================================

echo "=== FLEX Pipeline ==="
echo "Project:   ${CFG_PROJECT_ROOT}"
echo "Container: ${CONTAINER_TYPE} (${CONTAINER_PATH})"
echo "Job ID:    ${SLURM_JOB_ID:-local}"
echo "Task ID:   ${SLURM_ARRAY_TASK_ID:-N/A}"
echo "Temp:      ${JOB_TEMP_DIR}"
echo "DRY_RUN:   ${DRY_RUN}"
echo "Start:     $(date)"
echo "====================="
