#!/bin/bash
#SBATCH --job-name=prep_kumar_labels
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=.slurm_stubs/prep_kumar_labels_%j.out
#SBATCH --error=.slurm_stubs/prep_kumar_labels_%j.err

# Prepare per-compartment Kumar label CSVs for scIB bio scoring.
# See scripts/step_17_integration_sweep/source/prep_kumar_labels.py for details.

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
exec >> "${HPC_LOGS_DIR}/prep_kumar_labels_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/prep_kumar_labels_${SLURM_JOB_ID}.err"
CELL_ANN="${CFG_PREPROCESSING_ROOT}/preprocessing_wrapup/cell_annotations.csv"
LABELS_OUT="${CFG_PREPROCESSING_ROOT}/preprocessing_wrapup/labels_for_scib"

run_python_singularity \
    "${CFG_SCRIPTS_DIR}/step_17_integration_sweep/source/prep_kumar_labels.py" \
    --cell-annotations "$CELL_ANN" \
    --output-dir "$LABELS_OUT"
