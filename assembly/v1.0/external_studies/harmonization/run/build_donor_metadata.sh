#!/bin/bash
#SBATCH --job-name=build_donor_meta
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Phase 3: Build unified donor metadata using verified mappings

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

echo "==========================================="
echo "Phase 3: Build Unified Donor Metadata"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Started: $(date)"
echo "==========================================="

module load singularity

CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"
BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}"
SCRIPT="${BASE_DIR}/scripts/build_donor_metadata.py"

singularity exec \
    ${BIND_MOUNTS} \
    "${CONTAINER}" \
    python "${SCRIPT}" \
        --base-dir "${BASE_DIR}"

echo ""
echo "==========================================="
echo "Completed: $(date)"
echo "==========================================="
