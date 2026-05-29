#!/bin/bash
#SBATCH --job-name=map_pal
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00

# =============================================================================
# Pal Substudy Cell ID Mapping
# =============================================================================
# Maps a single Pal substudy to iHBCA reference using custom_pal mapping.
#
# Usage:
#   sbatch --export=STUDY=pal_norm_epi studies/pal/run/map_pal_single.sh
#   sbatch --export=STUDY=pal_norm_total studies/pal/run/map_pal_single.sh
#   sbatch --export=STUDY=pal_norm_b1 studies/pal/run/map_pal_single.sh
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../../../config/load_paths.sh"


if [ -z "${STUDY:-}" ]; then
    echo "ERROR: STUDY environment variable not set"
    exit 1
fi

case "${STUDY}" in
    pal_norm_epi|pal_norm_total|pal_norm_b1) ;;
    *) echo "ERROR: Invalid STUDY: ${STUDY}"; exit 1 ;;
esac

echo "=== Map ${STUDY} Cells to iHBCA ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Started: $(date)"
echo

PROJECT_DIR="${PROJECT_IHBCAV1_UPLOAD}"
PHASE_A_DIR="${PROJECT_DIR}/external_studies/harmonization"
R_SCRIPT="${PHASE_A_DIR}/scripts/map_study_generic.R"
MANIFEST="${PHASE_A_DIR}/config/study_manifest.yaml"

R_CONTAINER="${CONTAINER_R_SPATIAL_4_3_3}"
R_LIBS_USER="${R_LIBS_R_SPATIAL_4_3_3}"

module load singularity

singularity exec \
    --bind ${SHARED_LAB}:${SHARED_LAB}:ro \
    --bind ${USER_ROOT}:${USER_ROOT}:rw \
    ${BIND_MOUNTS} \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    "${R_CONTAINER}" \
    Rscript "${R_SCRIPT}" \
        --study "${STUDY}" \
        --manifest "${MANIFEST}"

echo
echo "=== ${STUDY} Mapping Complete ==="
echo "Finished: $(date)"
