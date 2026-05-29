#!/bin/bash
#SBATCH --job-name=audit_meta
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Phase A.2 Step 1: Column Inventory Audit
# Analyzes unified_donor_metadata.csv to create column_inventory.yaml

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start: $(date)"

# Paths
BASE_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

module load singularity/3.11.3

singularity exec \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    ${CONTAINER} \
    python3 ${BASE_DIR}/scripts/audit_donor_metadata.py

echo "Done: $(date)"
