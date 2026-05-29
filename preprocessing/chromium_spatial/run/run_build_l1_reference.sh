#!/bin/bash
#SBATCH --job-name=l1_ref_build
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=12
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/l1_ref_build_%j.out
#SBATCH --error=.slurm_stubs/l1_ref_build_%j.err

# Rebuild the iHBCA L1 pseudobulk reference, restricted to the FLEX 18K gene
# panel. One-shot, ~2.12M iHBCA cells collapsed to 11 L1 types.
#
# Container:   scgpt_gpu (Python + numpy/pandas/scipy sparse)
# Tier 2:      4cpu/16G/30min — sparse NPZ load is <8G, dense 11x18K output
#              fits easily
# Pattern src: run_compartment_integration_report.sh (common.sh sourcing,
#              PIPELINE_ROOT hardcode for SLURM spool safety)
#
# DRY_RUN=1 sbatch ... runs the script's --dry-run preamble (validates
# inputs, exits before reading counts).

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
exec >> "${HPC_LOGS_DIR}/l1_ref_build_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/l1_ref_build_${SLURM_JOB_ID}.err"
SCRIPT="${CFG_SCRIPTS_DIR}/step_19a_l1_annotation/09_build_l1_reference.py"

# The FLEX panel is identical across all compartments — any genes.tsv works
FLEX_PANEL="${CFG_PROJECT_ROOT}/outputs/integration_intermediate/Epithelial/scvi_n100/genes.tsv"

OUTDIR="${CFG_PROJECT_ROOT}/outputs/annotation/l1/reference"
mkdir -p "${OUTDIR}"

echo "Script:        ${SCRIPT}"
echo "h5 Epi:        ${CFG_IHBCA_H5_EPI}"
echo "h5 Stromal:    ${CFG_IHBCA_H5_STR}"
echo "h5 Immune:     ${CFG_IHBCA_H5_IMM}"
echo "cell_metadata: ${CFG_IHBCA_ASSEMBLY_CELL_METADATA}"
echo "FLEX panel:    ${FLEX_PANEL}"
echo "Output dir:    ${OUTDIR}"

run_python_singularity "${SCRIPT}" \
    --h5-epi        "${CFG_IHBCA_H5_EPI}" \
    --h5-str        "${CFG_IHBCA_H5_STR}" \
    --h5-imm        "${CFG_IHBCA_H5_IMM}" \
    --cell-metadata "${CFG_IHBCA_ASSEMBLY_CELL_METADATA}" \
    --flex-panel    "${FLEX_PANEL}" \
    --gene-mapping  "${CFG_GENE_SYMBOL_TO_ENSEMBL}" \
    --outdir        "${OUTDIR}" \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
