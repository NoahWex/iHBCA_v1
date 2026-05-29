#!/bin/bash
#SBATCH --job-name=xen_qc
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/08_assembly/qc_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/08_assembly/qc_%j.err

# Phase 8 — per-sample QC bar chart + CSV (n_cells, passes_all rate). Reads
# cell_annotations.csv (must exist first).

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

JOB_TMP=/tmp/xen_qc_${SLURM_JOB_ID}
mkdir -p "${JOB_TMP}"
trap 'rm -rf "${JOB_TMP}"' EXIT

module load singularity
singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/pipeline/scripts/08_assembly/05_per_sample_qc.py" \
        --project-root "${PROJECT_ROOT}"

echo "Done: $(date)"
