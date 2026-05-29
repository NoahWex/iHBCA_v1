#!/bin/bash
# _common.sh — Shared config loader for Track B (Xenium pipeline) wrappers.
# Wrappers must set CONTAINER_TYPE before sourcing this file:
#   CONTAINER_TYPE=python_scvi_gpu_2025Q3 | python_spatial_2025Q4E2 | r_spatial_4.3.3
# Pattern source: coordination/staging/track_3/run/_common.sh (post-CP5.7.10).

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
COMMON_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRACK_B_ROOT="$(cd "$COMMON_DIR/.." && pwd)"

# Walk up from TRACK_B_ROOT to find publication/config/load_paths.sh
if [ -z "${USER_ROOT:-}" ] || [ -z "${BIND_MOUNTS:-}" ]; then
    _SEARCH_DIR="${TRACK_B_ROOT}"
    PUBLICATION_CONFIG_DIR=""
    while [ "$_SEARCH_DIR" != "/" ]; do
        if [ -f "${_SEARCH_DIR}/publication/config/load_paths.sh" ]; then
            PUBLICATION_CONFIG_DIR="${_SEARCH_DIR}/publication/config"
            break
        fi
        _SEARCH_DIR="$(dirname "$_SEARCH_DIR")"
    done
    if [ -z "$PUBLICATION_CONFIG_DIR" ]; then
        echo "FATAL: could not locate publication/config/load_paths.sh by walking up from $TRACK_B_ROOT" >&2
        exit 2
    fi
    module load anaconda/2024.06 2>/dev/null || module load python/3.10.0 2>/dev/null || true
    # Source publication's load_paths.sh BEFORE setting our own PATHS_YAML
    # (load_paths.sh reads PATHS_YAML env var if set; need it to use its own default)
    unset PATHS_YAML
    source "${PUBLICATION_CONFIG_DIR}/load_paths.sh"
fi

export CFG_USER_ROOT="${USER_ROOT}"
export CFG_SHARED_LAB="${SHARED_LAB:-/share/crsp/lab/dalawson/share}"
# Also export the upstream paths the Python parser may need to seed substitutions.
export CFG_UPSTREAM_FLEX_PREP="${CFG_USER_ROOT}/Spatial_HBCA_preprocessing"

# Now set OUR paths.yaml
TRACK_B_PATHS_YAML="${TRACK_B_ROOT}/config/paths.yaml"
if [ ! -f "$TRACK_B_PATHS_YAML" ]; then
    echo "FATAL: track_b paths.yaml not found at $TRACK_B_PATHS_YAML" >&2
    exit 2
fi

# Parse track_b paths.yaml and export CFG_* (substituting ${var} refs against CFG_*)
eval "$(TRACK_B_PATHS_YAML="${TRACK_B_PATHS_YAML}" python3 << 'PYEOF'
import yaml, os, re
with open(os.environ['TRACK_B_PATHS_YAML']) as f:
    cfg = yaml.safe_load(f) or {}
seed = {k[len('CFG_'):].lower(): v for k, v in os.environ.items() if k.startswith('CFG_')}
resolved = dict(seed)
def resolve(val):
    if not isinstance(val, str): return val
    out = val
    for m in re.findall(r'\$\{([a-z_0-9]+)\}', val):
        if m in resolved:
            out = out.replace('${' + m + '}', str(resolved[m]))
        else:
            return None
    return out
for _ in range(10):
    progress = False
    for k, v in cfg.items():
        if k in resolved or not isinstance(v, str): continue
        r = resolve(v)
        if r is not None and '${' not in r:
            resolved[k] = r
            progress = True
    if not progress: break
for k, v in resolved.items():
    if k in seed: continue
    if isinstance(v, str):
        safe = v.replace('"', '\\"')
        print('export CFG_' + k.upper() + '="' + safe + '"')
PYEOF
)"

for required in CFG_PROJECT_ROOT CFG_PREPROCESSING_XENIUM_ROOT CFG_XENIUM_MANIFEST CFG_XENIUM_COUNTS_ROOT; do
    if [ -z "${!required:-}" ]; then
        echo "FATAL: required token $required did not resolve from paths.yaml" >&2
        exit 2
    fi
done

if [ -z "${CONTAINER_TYPE:-}" ]; then
    echo "FATAL: CONTAINER_TYPE not set before sourcing _common.sh" >&2
    exit 2
fi
CONTAINER_VAR="CONTAINER_$(echo "$CONTAINER_TYPE" | tr '[:lower:].' '[:upper:]_')"
CONTAINER_PATH="${!CONTAINER_VAR:-}"
if [ -z "$CONTAINER_PATH" ]; then
    echo "FATAL: container '$CONTAINER_TYPE' (var $CONTAINER_VAR) not exported by load_paths.sh" >&2
    exit 2
fi

CONTAINER_GPU="false"
case "$CONTAINER_TYPE" in
    *gpu*|python_scvi_gpu_*) CONTAINER_GPU="true" ;;
esac

R_LIBS_USER_HOST=""
R_LIBS_USER_CONTAINER="/home/jovyan/R/library"
case "$CONTAINER_TYPE" in
    r_spatial_4.3.3)
        R_LIBS_USER_HOST="${R_LIBS_R_SPATIAL_4_3_3:-}"
        ;;
esac

PYTHONUSERBASE_HOST=""
PYTHONUSERBASE_CONTAINER=""
case "$CONTAINER_TYPE" in
    python_scvi_gpu_2025Q3)
        PYTHONUSERBASE_HOST="/pub/nwechter/python_user_packages/scGPT_GPU_2025Q3"
        PYTHONUSERBASE_CONTAINER="$PYTHONUSERBASE_HOST"
        ;;
esac

BIND_FLAGS=""
for mount in ${BIND_MOUNT_LIST:-}; do
    BIND_FLAGS="${BIND_FLAGS} --bind ${mount}"
done
[[ "$BIND_FLAGS" != *"/share/crsp/lab/dalawson"* ]] && BIND_FLAGS="${BIND_FLAGS} --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw"
[[ "$BIND_FLAGS" != *"/pub/nwechter"* ]] && BIND_FLAGS="${BIND_FLAGS} --bind /pub/nwechter:/pub/nwechter:rw"
[[ "$BIND_FLAGS" != *"/dfs7"* ]] && BIND_FLAGS="${BIND_FLAGS} --bind /dfs7:/dfs7:ro"
[[ "$BIND_FLAGS" != *"/dfs8"* ]] && BIND_FLAGS="${BIND_FLAGS} --bind /dfs8:/dfs8:ro"

JOB_TEMP_BASE="${CFG_HPC_LOGS_DIR:-/tmp}/job_temp"
mkdir -p "${JOB_TEMP_BASE}" 2>/dev/null || JOB_TEMP_BASE="/tmp"
JOB_TEMP_DIR=$(mktemp -d "${JOB_TEMP_BASE}/job_${SLURM_JOB_ID:-local}_XXXXXX")
trap 'rm -rf "$JOB_TEMP_DIR"' EXIT

_load_singularity() { module load singularity 2>/dev/null || true; }

run_python_singularity() {
    local script="$1"; shift
    _load_singularity
    local nv_flag=""
    if [ "${CONTAINER_GPU}" = "true" ]; then nv_flag="--nv"; fi
    local env_flags=""
    if [ -n "${PYTHONUSERBASE_HOST}" ]; then
        env_flags="--env PYTHONUSERBASE=${PYTHONUSERBASE_CONTAINER}"
    fi
    mkdir -p "${JOB_TEMP_DIR}/numba_cache" "${JOB_TEMP_DIR}/mpl_config"
    singularity exec ${nv_flag} --no-home \
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
    singularity exec --no-home \
        ${BIND_FLAGS} \
        ${r_bind} \
        --bind "${JOB_TEMP_DIR}:${JOB_TEMP_DIR}:rw" \
        ${r_env} \
        --env "TMPDIR=${JOB_TEMP_DIR}" \
        --env "LANG=en_US.UTF-8" \
        --env "LC_ALL=en_US.UTF-8" \
        "${CONTAINER_PATH}" \
        Rscript "${script}" "$@"
}

echo "=== Track B Xenium Pipeline ==="
echo "Project:    ${CFG_PROJECT_ROOT}"
echo "Container:  ${CONTAINER_TYPE} (${CONTAINER_PATH})"
echo "Job ID:     ${SLURM_JOB_ID:-local}"
echo "Temp:       ${JOB_TEMP_DIR}"
echo "Start:      $(date)"
echo "==============================="
