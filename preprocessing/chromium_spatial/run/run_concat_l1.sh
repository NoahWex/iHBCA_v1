#!/bin/bash
#SBATCH --job-name=l1_concat
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/l1_concat_%j.out
#SBATCH --error=.slurm_stubs/l1_concat_%j.err

# Concatenate per-compartment harmonized L1 CSVs into a single 270,035-row
# file. Runs validation gates (row count, duplicate check).
#
# Container:   r_spatial
# Tier 1:      2cpu/4G/30min — reads 3 CSVs (~270K rows total), rbinds, writes
# Pattern src: run_compartment_integration_report.sh (config/container load)

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=r_spatial
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/l1_concat_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/l1_concat_${SLURM_JOB_ID}.err"
SCRIPT="${CFG_SCRIPTS_DIR}/step_19a_l1_annotation/09d_concat_l1.R"
L1_ROOT="${CFG_PROJECT_ROOT}/outputs/annotation/l1"

EPI_CSV="${L1_ROOT}/Epithelial/harmonized_l1_labels_Epithelial.csv"
STR_CSV="${L1_ROOT}/Stromal/harmonized_l1_labels_Stromal.csv"
IMM_CSV="${L1_ROOT}/Immune/harmonized_l1_labels_Immune.csv"

OUT_DIR="${L1_ROOT}/all"
OUT_CSV="${OUT_DIR}/harmonized_l1_labels.csv"
mkdir -p "${OUT_DIR}"

echo "Epithelial:  ${EPI_CSV}"
echo "Stromal:     ${STR_CSV}"
echo "Immune:      ${IMM_CSV}"
echo "Output:      ${OUT_CSV}"

run_r_singularity "${SCRIPT}" \
    --epithelial "${EPI_CSV}" \
    --stromal    "${STR_CSV}" \
    --immune     "${IMM_CSV}" \
    --out        "${OUT_CSV}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
