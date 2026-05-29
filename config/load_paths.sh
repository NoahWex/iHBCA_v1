#!/usr/bin/env bash
# =============================================================================
# load_paths.sh — Shell loader for publication path resolution
# =============================================================================
# Source this from SLURM job scripts to get resolved paths as shell variables.
# Uses inline Python to parse paths.yaml (same pattern as load_config.sh).
#
# Usage:
#   source "path/to/publication/config/load_paths.sh"
#
# Exports:
#   USER_ROOT, SHARED_DATASETS, HPC_HOME
#   PROJECT_PUBLICATION, PROJECT_IHBCAV1_UPLOAD, PROJECT_SPATIAL_HBCA, ...
#   CONTAINER_R_SPATIAL_4_3_3, CONTAINER_PYTHON_SPATIAL_2025Q2, ...
#   R_LIBS_R_SPATIAL_4_3_3, PYTHON_LIBS_SPATIAL_2025Q4E2, ...
#   BIND_MOUNTS (space-separated --bind args)
#   SLURM_CFG_ACCOUNT, SLURM_CFG_PARTITION
# =============================================================================

# Ensure python3 has yaml available. Default compute-node python lacks PyYAML;
# anaconda/2024.06 provides it (verified 2026-04-27 on UCI HPC3, yaml 6.0.1).
# Load only on HPC compute (login nodes / local already have working python).
if [ -n "${SLURM_JOB_ID:-}" ] && command -v module &>/dev/null; then
    module load anaconda/2024.06 2>/dev/null || module load python/3.10.0 2>/dev/null || true
fi

# Locate paths.yaml relative to this script
_LOAD_PATHS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PATHS_YAML="${PATHS_YAML:-${_LOAD_PATHS_DIR}/paths.yaml}"

if [ ! -f "${PATHS_YAML}" ]; then
    echo "ERROR: paths.yaml not found at ${PATHS_YAML}" >&2
    return 1 2>/dev/null || exit 1
fi

# Detect environment: IHBCA_ENV override or hostname-based
if [ -z "${IHBCA_ENV:-}" ]; then
    _hostname="$(hostname)"
    case "${_hostname}" in
        *hpc3*|login-*|compute-*) IHBCA_ENV="hpc" ;;
        *) IHBCA_ENV="local" ;;
    esac
fi
export IHBCA_ENV

# Parse paths.yaml and export all variables via inline Python
eval "$(python3 -c "
import yaml, os

env = os.environ.get('IHBCA_ENV', 'hpc')

with open('${PATHS_YAML}') as f:
    cfg = yaml.safe_load(f)

root = cfg['roots'][env]
print(f'export USER_ROOT=\"{root}\"')

# Additional roots (absolute, environment-independent)
for key in ('shared_datasets', 'hpc_home'):
    val = cfg['roots'].get(key, '')
    if val:
        var_name = key.upper()
        print(f'export {var_name}=\"{val}\"')

# Projects (resolved)
for key, val in cfg.get('projects', {}).items():
    resolved = val if val.startswith('/') else f'{root}/{val}'
    var_name = f'PROJECT_{key.upper()}'
    print(f'export {var_name}=\"{resolved}\"')

# Sources (resolved)
for key, val in cfg.get('sources', {}).items():
    resolved = val if val.startswith('/') else f'{root}/{val}'
    var_name = f'SOURCE_{key.upper()}'
    print(f'export {var_name}=\"{resolved}\"')

# Containers
for key, val in cfg.get('containers', {}).items():
    var_name = f'CONTAINER_{key.upper().replace(\".\", \"_\")}'
    print(f'export {var_name}=\"{val}\"')

# R library paths
for key, val in cfg.get('r_libs', {}).items():
    var_name = f'R_LIBS_{key.upper().replace(\".\", \"_\")}'
    print(f'export {var_name}=\"{val}\"')
# Legacy single r_libs_user
rlu = cfg.get('r_libs_user', '')
if rlu:
    print(f'export R_LIBS_USER_DEFAULT=\"{rlu}\"')

# Python library paths
for key, val in cfg.get('python_libs', {}).items():
    var_name = f'PYTHON_LIBS_{key.upper().replace(\".\", \"_\")}'
    print(f'export {var_name}=\"{val}\"')

# Bind mounts as repeated --bind flags
mounts = cfg.get('bind_mounts', [])
bind_args = ' '.join(f'--bind {m}' for m in mounts)
print(f'export BIND_MOUNTS=\"{bind_args}\"')
# Also as a plain list for iteration
mount_list = ' '.join(mounts)
print(f'export BIND_MOUNT_LIST=\"{mount_list}\"')

# Top-level scalar paths (e.g., hpc_logs_dir). Exported in UPPER_CASE.
# IMPORTANT: SLURM SBATCH --output= directives can substitute these
# variables ONLY if the user sources this file BEFORE invoking sbatch.
# Pattern: source publication/config/load_paths.sh THEN sbatch wrapper.sh
for key in ('hpc_logs_dir',):
    val = cfg.get(key, '')
    if val:
        print(f'export {key.upper()}=\"{val}\"')

# SLURM defaults (prefixed to avoid collision with SLURM env vars)
slurm = cfg.get('slurm', {})
for key, val in slurm.items():
    var_name = f'SLURM_CFG_{key.upper()}'
    print(f'export {var_name}=\"{val}\"')
" 2>&1)"

# Verify the critical variable was set
if [ -z "${USER_ROOT:-}" ]; then
    echo "ERROR: load_paths.sh failed to resolve USER_ROOT from ${PATHS_YAML}" >&2
    return 1 2>/dev/null || exit 1
fi
