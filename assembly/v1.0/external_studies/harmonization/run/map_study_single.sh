#!/bin/bash
#SBATCH --job-name=map_study
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# =============================================================================
# Single Study Cell ID Mapping
# =============================================================================
# Maps a single component study to iHBCA reference using study_manifest.yaml.
#
# Usage:
#   sbatch --export=STUDY=kumar run/map_study_single.sh
#   sbatch --export=STUDY=murrow run/map_study_single.sh
#
# The STUDY variable selects which study to process from the manifest.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

# Verify STUDY is set
if [ -z "${STUDY:-}" ]; then
    echo "ERROR: STUDY environment variable not set"
    echo "Usage: sbatch --export=STUDY=kumar run/map_study_single.sh"
    exit 1
fi

echo "=== Map ${STUDY} Cells to iHBCA ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Study: ${STUDY}"
echo "Started: $(date)"
echo

# Paths
PROJECT_DIR="${PROJECT_IHBCAV1_UPLOAD}"
PHASE_A_DIR="${PROJECT_DIR}/external_studies/harmonization"
MANIFEST="${PHASE_A_DIR}/config/study_manifest.yaml"
R_SCRIPT="${PHASE_A_DIR}/scripts/map_study_generic.R"
PY_SCRIPT="${PHASE_A_DIR}/scripts/map_study_generic.py"

# Containers
R_CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
PY_CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

# R library path (for Seurat, argparse, etc.)
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

# Load modules
module load singularity

# Verify manifest exists
if [ ! -f "${MANIFEST}" ]; then
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
fi

# Get file format from manifest to determine which script to use
# Use Python to parse YAML and extract format
FORMAT=$(singularity exec \
    ${BIND_MOUNTS} \
    "${PY_CONTAINER}" \
    python3 -c "
import yaml
with open('${MANIFEST}') as f:
    m = yaml.safe_load(f)
if '${STUDY}' not in m['studies']:
    print('STUDY_NOT_FOUND')
else:
    print(m['studies']['${STUDY}']['format'])
")

echo "File format: ${FORMAT}"

if [ "${FORMAT}" == "STUDY_NOT_FOUND" ]; then
    echo "ERROR: Study '${STUDY}' not found in manifest"
    exit 1
fi

# Run appropriate mapper based on format
if [ "${FORMAT}" == "rds" ]; then
    echo "Using R mapper for RDS file..."
    echo

    singularity exec \
        ${BIND_MOUNTS} \
        --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
        --env "R_LIBS_USER=/home/jovyan/R/library" \
        "${R_CONTAINER}" \
        Rscript "${R_SCRIPT}" \
            --study "${STUDY}" \
            --manifest "${MANIFEST}"

elif [ "${FORMAT}" == "h5ad" ]; then
    echo "Using Python mapper for h5ad file..."
    echo

    singularity exec \
        ${BIND_MOUNTS} \
        --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
        --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
        "${PY_CONTAINER}" \
        python3 "${PY_SCRIPT}" \
            --study "${STUDY}" \
            --manifest "${MANIFEST}"

else
    echo "ERROR: Unknown format '${FORMAT}' for study '${STUDY}'"
    exit 1
fi

echo
echo "=== ${STUDY} Mapping Complete ==="
echo "Finished: $(date)"
