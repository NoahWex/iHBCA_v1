#!/bin/bash
#SBATCH --job-name=step_04_harmonize
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_04_harmonize_%j.out
#SBATCH --error=.slurm_stubs/step_04_harmonize_%j.err

# Step 04 (harmonize): Aggregate all per-sample singler_labels.csv files and
# produce harmonized_l1_labels.csv via SC/SN reconciliation + Leiden vote.
# Runs AFTER all per-sample classify array tasks complete.
#
# Resources (Tier 2): 2 CPU / 16G / 30min
#   CSV aggregation of 62 samples is CPU-light, memory-light.
#
# Required env vars (or pass via sbatch --export=ALL,LEIDEN_CSV=...):
#   LEIDEN_CSV — path to multi-resolution Leiden CSV from the intra-patient merge
#               (per-patient concatenation + Leiden, no batch correction — patients
#                are the batches so a straight merge is sufficient at step 04 time)
# Optional:
#   SINGLER_ROOT — defaults to CFG_PREPROCESSING_STEP_04 (where per-sample dirs live)
#
# Pattern source: run/run_aggregate_l2.sh (Tier 2, single-task, _common.sh pattern)

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=scgpt_gpu  # python container with pandas
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/step_04_harmonize_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_04_harmonize_${SLURM_JOB_ID}.err"
LEIDEN_CSV="${LEIDEN_CSV:-}"
if [ -z "$LEIDEN_CSV" ]; then
  echo "FATAL: LEIDEN_CSV must be set (path to multi-resolution Leiden CSV from sweep winner)" >&2
  exit 1
fi

SCRIPT="${CFG_SCRIPTS_DIR}/step_04_label_transfer/03c_harmonize_labels.py"
SINGLER_ROOT="${SINGLER_ROOT:-${CFG_PREPROCESSING_STEP_04}}"
OUTPUT_DIR="${CFG_PREPROCESSING_STEP_04}/harmonized"

echo "Step:         step_04_harmonize"
echo "Singler root: ${SINGLER_ROOT}"
echo "Leiden CSV:   ${LEIDEN_CSV}"
echo "Output:       ${OUTPUT_DIR}"

run_python_singularity "${SCRIPT}" \
    --singler-root "${SINGLER_ROOT}" \
    --leiden-csv   "${LEIDEN_CSV}" \
    --output-dir   "${OUTPUT_DIR}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
