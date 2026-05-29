#!/bin/bash
#SBATCH --job-name=smooth_full_l1
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/smooth_full_l1_%j.out
#SBATCH --error=.slurm_stubs/smooth_full_l1_%j.err
#
# Cluster-vote smoothing for full-object L1: fills the ~19K cells that passed
# canonical 5-way QC but were dropped from per-compartment Seurat builds.
#
# Container:   scgpt_gpu (pandas only, light)
# Tier 1:      2cpu/8G/20min
# Pattern src: 09f_smooth_full_l1.py

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=scgpt_gpu
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/smooth_full_l1_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/smooth_full_l1_${SLURM_JOB_ID}.err"
OBS="${PIPELINE_ROOT}/outputs/integration_intermediate/full/scvi_n100/obs.csv"
LABELS="${PIPELINE_ROOT}/outputs/annotation/l1/all/harmonized_l1_labels.csv"
OUT="${PIPELINE_ROOT}/outputs/annotation/l1/all/harmonized_l1_labels_full.csv"
VOTE_LEIDEN="${VOTE_LEIDEN:-leiden_2.0}"
LABEL_COL="${LABEL_COL:-L1.0}"

echo "obs:         $OBS"
echo "labels:      $LABELS"
echo "vote leiden: $VOTE_LEIDEN"
echo "out:         $OUT"

SCRIPT="${CFG_SCRIPTS_DIR}/step_19a_l1_annotation/09f_smooth_full_l1.py"
run_python_singularity \
    "$SCRIPT" \
    --obs "$OBS" \
    --labels "$LABELS" \
    --label-col "$LABEL_COL" \
    --vote-leiden "$VOTE_LEIDEN" \
    --out "$OUT"

echo "Exit: $?, End: $(date)"
