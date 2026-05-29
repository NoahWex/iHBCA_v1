#!/bin/bash
#SBATCH --job-name=l1_harmonize
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/l1_harmonize_%A_%a.out
#SBATCH --error=.slurm_stubs/l1_harmonize_%A_%a.err

# Per-compartment harmonization of SingleR L1 predictions via cluster_annotation
# overrides + leiden_n30_r1.0 vote consolidation.
#
# Container:   r_spatial
# Tier 4:      8cpu/48G/2h — matches the pattern tier; actual workload is
#              small (tabular operations only) but Seurat load dominates
# Pattern src: run_compartment_integration_report.sh
#
# Tasks:
#   0  Epithelial
#   1  Stromal
#   2  Immune

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
exec >> "${HPC_LOGS_DIR}/l1_harmonize_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/l1_harmonize_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}

RMD="${CFG_SCRIPTS_DIR}/step_19a_l1_annotation/09c_l1_harmonize.Rmd"
SINGLER_CSV="${CFG_PROJECT_ROOT}/outputs/annotation/l1/${COMPARTMENT}/singler_l1_predictions_${COMPARTMENT}.csv"
SEURAT_RDS="${CFG_PROJECT_ROOT}/outputs/annotation/l1/${COMPARTMENT}/seurat.rds"
OUTDIR="${CFG_PROJECT_ROOT}/outputs/annotation/l1/${COMPARTMENT}"
OUTPUT_FILE="l1_harmonize_${COMPARTMENT}.html"

mkdir -p "${OUTDIR}"

if [ "${DRY_RUN}" = "1" ]; then
    DRY_R="TRUE"
else
    DRY_R="FALSE"
fi

echo "Compartment:  ${COMPARTMENT}"
echo "Rmd:          ${RMD}"
echo "SingleR CSV:  ${SINGLER_CSV}"
echo "Seurat RDS:   ${SEURAT_RDS}"
echo "Output dir:   ${OUTDIR}"
echo "Dry run:      ${DRY_R}"

run_r_singularity -e "rmarkdown::render('${RMD}', output_dir = '${OUTDIR}', output_file = '${OUTPUT_FILE}', intermediates_dir = '${JOB_TEMP_DIR}', knit_root_dir = '${CFG_PROJECT_ROOT}', params = list(compartment = '${COMPARTMENT}', singler_csv = '${SINGLER_CSV}', seurat_rds = '${SEURAT_RDS}', outdir = '${OUTDIR}', dry_run = ${DRY_R}))"

echo "Exit: $?, End: $(date)"
