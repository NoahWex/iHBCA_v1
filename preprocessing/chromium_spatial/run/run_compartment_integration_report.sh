#!/bin/bash
#SBATCH --job-name=comp_int_rep
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/comp_int_rep_%A_%a.out
#SBATCH --error=.slurm_stubs/comp_int_rep_%A_%a.err

# Render the per-compartment integration report (scripts/08_compartment_integration_report.Rmd).
# One task per compartment, parallel. Tier 4 (8cpu/48G/2h).
# Sources run/_common.sh for config + container/bind setup; calls promoted
# scripts/08_compartment_integration_report.Rmd against the scvi_n100
# winner config.
#   - per-task intermediates_dir from _common.sh JOB_TEMP_DIR
#
# Task layout:
#   0  Epithelial  102,952 cells
#   1  Stromal     119,141 cells
#   2  Immune       47,942 cells
#
# DRY_RUN=1 sbatch ... runs the Rmd's validation preamble and knit_exit.

# SLURM copies the script to a spool dir, so BASH_SOURCE points away from the
# real run/ directory. Hardcode the project root and source _common.sh from
# there. Update this if the pipeline moves.
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
exec >> "${HPC_LOGS_DIR}/comp_int_rep_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/comp_int_rep_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
CONFIG_LABEL=scvi_n100

RMD="${CFG_PROJECT_ROOT}/scripts/step_18_integration_preview/08_compartment_integration_report.Rmd"
INTERMEDIATE_DIR="${CFG_PROJECT_ROOT}/outputs/integration_intermediate/${COMPARTMENT}/${CONFIG_LABEL}"
OUTDIR="${CFG_PROJECT_ROOT}/outputs/integration_reports/${COMPARTMENT}/${CONFIG_LABEL}"
OUTPUT_FILE="phase1c_${COMPARTMENT}_${CONFIG_LABEL}_integration.html"

mkdir -p "${OUTDIR}"

# DRY_RUN flag is consumed by the Rmd's validation preamble (knit_exit on dry).
# Translate the bash DRY_RUN var into an R logical for the params list.
if [ "${DRY_RUN}" = "1" ]; then
    DRY_R="TRUE"
else
    DRY_R="FALSE"
fi

echo "Compartment:  ${COMPARTMENT}"
echo "Config label: ${CONFIG_LABEL}"
echo "Rmd:          ${RMD}"
echo "Output:       ${OUTDIR}/${OUTPUT_FILE}"
echo "Dry run:      ${DRY_R}"

run_r_singularity -e "rmarkdown::render('${RMD}', output_dir = '${OUTDIR}', output_file = '${OUTPUT_FILE}', intermediates_dir = '${JOB_TEMP_DIR}', knit_root_dir = '${CFG_PROJECT_ROOT}', params = list(project_root = '${CFG_PROJECT_ROOT}', compartment = '${COMPARTMENT}', config_label = '${CONFIG_LABEL}', upstream_root = '${CFG_UPSTREAM_ROOT}', intermediate_dir = '${INTERMEDIATE_DIR}', outdir = '${CFG_PROJECT_ROOT}/outputs/integration_reports', dry_run = ${DRY_R}))"

echo "Exit: $?, End: $(date)"
