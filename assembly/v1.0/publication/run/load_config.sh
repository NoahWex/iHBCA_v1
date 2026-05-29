#!/usr/bin/env bash
# =============================================================================
# load_config.sh — Export pipeline.yaml paths as shell variables
# =============================================================================
# Source this from SLURM job scripts after setting REPO_ROOT:
#
#   REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"  # from load_paths.sh
#   source "${REPO_ROOT}/publication/run/load_config.sh"
#
# Exports:
#   PIPELINE_CONFIG          — absolute path to pipeline.yaml
#   OUTPUT_SOURCE_DATASETS   — source h5ad output dir
#   OUTPUT_INTEGRATED_OBJECTS — integrated h5ad output dir
#   OUTPUT_VALIDATION_REPORTS — validation report dir
#   OUTPUT_ENTRY_SHEETS      — entry sheet base dir
#   MAPPINGS_DIR             — publication/mappings
#   CONFIG_DIR               — publication/config
#   HARMONIZED_METADATA      — harmonized donor metadata CSV
#   STUDIES                  — space-separated study names
#   STUDY_FILENAMES          — parallel array of output filenames
# =============================================================================

if [ -z "${REPO_ROOT}" ]; then
    echo "ERROR: REPO_ROOT must be set before sourcing load_config.sh" >&2
    return 1 2>/dev/null || exit 1
fi

PIPELINE_CONFIG="${REPO_ROOT}/publication/config/pipeline.yaml"

if [ ! -f "${PIPELINE_CONFIG}" ]; then
    echo "ERROR: pipeline.yaml not found at ${PIPELINE_CONFIG}" >&2
    return 1 2>/dev/null || exit 1
fi

# Parse pipeline.yaml and export paths using inline Python
eval "$(python3 -c "
import yaml, sys

with open('${PIPELINE_CONFIG}') as f:
    cfg = yaml.safe_load(f)

repo = '${REPO_ROOT}'

# Output directories
for k, v in cfg['paths']['outputs'].items():
    name = 'OUTPUT_' + k.upper()
    print(f'export {name}=\"{repo}/{v}\"')

# Mapping and config directories
print(f'export MAPPINGS_DIR=\"{repo}/publication/mappings\"')
print(f'export CONFIG_DIR=\"{repo}/publication/config\"')

# Harmonized metadata
hm = cfg['paths']['inputs']['harmonized_metadata']
print(f'export HARMONIZED_METADATA=\"{repo}/{hm}\"')

# Study names and filenames as parallel arrays
studies = cfg['studies']
names = list(studies.keys())
filenames = [s['output_filename'] for s in studies.values()]
print(f'export STUDIES=\"{\" \".join(names)}\"')
print(f'export STUDY_FILENAMES=\"{\" \".join(filenames)}\"')
" 2>&1)"

# Convert STUDIES and STUDY_FILENAMES to bash arrays for convenience
# shellcheck disable=SC2206
STUDIES_ARRAY=(${STUDIES})
# shellcheck disable=SC2206
FILENAMES_ARRAY=(${STUDY_FILENAMES})
export STUDIES_ARRAY FILENAMES_ARRAY
