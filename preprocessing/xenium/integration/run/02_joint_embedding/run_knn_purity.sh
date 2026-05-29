#!/bin/bash
#SBATCH --job-name=knn_purity
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/02_joint_embedding/knn_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/02_joint_embedding/knn_%j.err

# Phase 2c benchmark — FLEX-in-joint-latent kNN purity across the 6 joint
# Concord configurations plus the FLEX-only Track A baseline. Produces the
# winner selection evidence used at Gate G2 (winner: joint_nuc_n100).

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

SWEEP=${PROJECT_ROOT}/.dev/02_integration/joint_sweep
BENCH=${PROJECT_ROOT}/.dev/02_integration/benchmark
FLEX_ONLY=${PROJECT_ROOT}/.dev/_legacy/01_flex_concord_panel/outputs/full_sweep/concord_n50
L0P5=${PROJECT_ROOT}/.dev/_legacy/02_flex_l0p5_nmp/outputs/l0p5_per_cell.csv
OUT=${BENCH}/outputs

JOB_TMP=/tmp/knn_${SLURM_JOB_ID}
mkdir -p "${JOB_TMP}" "${OUT}"
trap 'rm -rf "${JOB_TMP}"' EXIT

CONFIGS=(
  "nuc_n50=${SWEEP}/joint_nuc_n50/outputs/joint_latent.csv:${SWEEP}/joint_nuc_n50/outputs/joint_obs.csv"
  "nuc_n100=${SWEEP}/joint_nuc_n100/outputs/joint_latent.csv:${SWEEP}/joint_nuc_n100/outputs/joint_obs.csv"
  "whl_n50=${SWEEP}/joint_whl_n50/outputs/joint_latent.csv:${SWEEP}/joint_whl_n50/outputs/joint_obs.csv"
  "whl_n100=${SWEEP}/joint_whl_n100/outputs/joint_latent.csv:${SWEEP}/joint_whl_n100/outputs/joint_obs.csv"
  "cyto_n50=${SWEEP}/joint_cyto_n50/outputs/joint_latent.csv:${SWEEP}/joint_cyto_n50/outputs/joint_obs.csv"
  "cyto_n100=${SWEEP}/joint_cyto_n100/outputs/joint_latent.csv:${SWEEP}/joint_cyto_n100/outputs/joint_obs.csv"
  "flex_only=${FLEX_ONLY}/latent.csv:${FLEX_ONLY}/obs.csv"
)

module load singularity
singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "NUMBA_CACHE_DIR=${JOB_TMP}/numba" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/pipeline/scripts/02_joint_embedding/benchmark/01_knn_purity.py" \
        --configs "${CONFIGS[@]}" \
        --l0p5-csv "${L0P5}" \
        --out-dir "${OUT}" \
        -k 30

echo "Done: $(date)"
