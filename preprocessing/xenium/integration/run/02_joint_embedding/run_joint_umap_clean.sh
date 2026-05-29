#!/bin/bash
#SBATCH --job-name=joint_umap_clean
#SBATCH --account=dalawson_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/02_joint_embedding/umap_clean_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/02_joint_embedding/umap_clean_%j.err

# Phase 2 — RAPIDS GPU UMAP on the joint_nuc_n100_clean latent. Diagnostic
# UMAP used at Gate G2; the canonical UMAP is re-derived from joint_v6 in
# Phase 7.

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

WIN=${PROJECT_ROOT}/.dev/02_integration/joint_sweep/joint_nuc_n100_clean
JOB_TMP=/tmp/umap_clean_${SLURM_JOB_ID}
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
    python3 "${PROJECT_ROOT}/pipeline/scripts/02_joint_embedding/03_joint_umap.py" \
        --latent "${WIN}/outputs/joint_latent.csv" \
        --obs "${WIN}/outputs/joint_obs.csv" \
        --out-dir "${WIN}/outputs"

echo "Done umap_clean: $(date)"
