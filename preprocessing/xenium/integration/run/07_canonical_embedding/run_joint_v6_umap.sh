#!/bin/bash
#SBATCH --job-name=joint_v6_umap
#SBATCH --account=dalawson_lab
#SBATCH --partition=free-gpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=192G
#SBATCH --gres=gpu:1
#SBATCH --time=03:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/07_canonical_embedding/umap_v6_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/07_canonical_embedding/umap_v6_%j.err

# Phase 7a — RAPIDS GPU UMAP on the joint_v6 latent. Canonical UMAP
# coordinates exported in cell_annotations.csv (UMAP1, UMAP2 columns).
# Reuses the Phase 2 UMAP script; no code difference from phase 2, just
# a different latent input.

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

OUT=${PROJECT_ROOT}/pipeline/outputs/embedding
JOB_TMP=/tmp/umap_v6_${SLURM_JOB_ID}
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
        --latent "${OUT}/joint_latent.csv" \
        --obs "${OUT}/joint_obs.csv" \
        --out-dir "${OUT}"

echo "Done joint_v6 UMAP: $(date)"
