#!/bin/bash
#SBATCH --job-name=joint_leiden
#SBATCH --account=dalawson_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=01:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/03_cluster_alignment/joint_leiden_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/03_cluster_alignment/joint_leiden_%j.err

# Phase 3 — Leiden clustering on the full FLEX + Xenium joint latent (joint_v6).
# Produces joint_leiden_assignments.csv at resolutions {0.3, 0.5, 1.0, 1.5}.
# Used by 04_joint_markers.R to compute per-cluster markers across three
# expression sources.

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

JOB_TMP=/tmp/joint_leiden_${SLURM_JOB_ID}
mkdir -p "${JOB_TMP}"
trap 'rm -rf "${JOB_TMP}"' EXIT

module load singularity
singularity exec --nv \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "NUMBA_CACHE_DIR=${JOB_TMP}/numba" \
    --env "CUPY_CACHE_DIR=${JOB_TMP}/cupy" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/pipeline/scripts/03_cluster_alignment/03_joint_leiden.py" \
        --latent   "${JOINT_LATENT}" \
        --obs      "${JOINT_OBS}" \
        --out-path "${JOINT_LEIDEN}" \
        --k 30 \
        --resolutions 0.3 0.5 1.0 1.5

echo "Done: $(date)"
