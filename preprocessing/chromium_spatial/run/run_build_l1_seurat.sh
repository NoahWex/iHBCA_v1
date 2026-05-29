#!/bin/bash
#SBATCH --job-name=l1_build_seurat
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/l1_build_seurat_%A_%a.out
#SBATCH --error=.slurm_stubs/l1_build_seurat_%A_%a.err

# Build per-compartment Seurat objects for L1 SingleR classification.
# Lightweight reassembly from the flat intermediate (08a output) + patches +
# contamination CSV + cell_metadata.csv.
#
# Container:   r_spatial (Seurat, Matrix, argparse)
# Tier 4:      8cpu/48G/2h — matches run_compartment_integration_report.sh
#              (same pattern, same working-set scale per compartment)
# Pattern src: run/run_compartment_integration_report.sh (hardcoded
#              PIPELINE_ROOT, _common.sh sourcing, 3-task compartment array)
#
# Tasks:
#   0  Epithelial  102,952 cells
#   1  Stromal     119,141 cells
#   2  Immune       47,942 cells

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
exec >> "${HPC_LOGS_DIR}/l1_build_seurat_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/l1_build_seurat_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
CONFIG_LABEL=scvi_n100

SCRIPT="${CFG_SCRIPTS_DIR}/step_19a_l1_annotation/09a_build_l1_seurat.R"
INTERMEDIATE_DIR="${CFG_PROJECT_ROOT}/outputs/integration_intermediate/${COMPARTMENT}/${CONFIG_LABEL}"
PATCHES_DIR="${INTERMEDIATE_DIR}/patches"
OUTDIR="${CFG_PROJECT_ROOT}/outputs/annotation/l1/${COMPARTMENT}"
OUT_RDS="${OUTDIR}/seurat.rds"

mkdir -p "${OUTDIR}"

echo "Compartment:     ${COMPARTMENT}"
echo "Intermediate:    ${INTERMEDIATE_DIR}"
echo "Patches:         ${PATCHES_DIR}"
echo "Contamination:   ${CFG_UPSTREAM_CONTAMINATION}"
echo "Cell metadata:   ${CFG_UPSTREAM_CELL_METADATA}"
echo "Output RDS:      ${OUT_RDS}"

run_r_singularity "${SCRIPT}" \
    --intermediate-dir   "${INTERMEDIATE_DIR}" \
    --patches-dir        "${PATCHES_DIR}" \
    --contamination-csv  "${CFG_UPSTREAM_CONTAMINATION}" \
    --cell-metadata-csv  "${CFG_UPSTREAM_CELL_METADATA}" \
    --compartment        "${COMPARTMENT}" \
    --out                "${OUT_RDS}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
