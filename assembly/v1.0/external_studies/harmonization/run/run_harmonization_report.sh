#!/bin/bash
#SBATCH --job-name=harm_report
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:15:00

# Generate HTML harmonization report

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start: $(date)"

BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

module load singularity/3.11.3

singularity exec \
    ${BIND_MOUNTS} \
    ${CONTAINER} \
    python3 ${BASE_DIR}/scripts/generate_harmonization_report.py

echo "Done: $(date)"
