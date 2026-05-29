#!/bin/bash
#SBATCH --job-name=step_04_train_singler
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/step_04_train_%j.out
#SBATCH --error=.slurm_stubs/step_04_train_%j.err

# Step 04 (train): Build SingleR models from Kumar 2023 SC + SN references.
# Runs ONCE before the per-sample classify array (run_step_04_array.sh).
#
# Resources (Tier 5): 4 CPU / 128G / 4h
#   SC reference: 714K cells, MTX load + pseudobulk peaks ~80-100G.
#   SN reference: smaller Seurat RDS, ~20G.
#
# Required env vars (set before sbatch or export=ALL):
#   SC_REFERENCE_DIR  — path to Kumar SC converted MTX dir (not in paths.yaml)
#   SN_REFERENCE_RDS  — path to annotated_hbca_nuclei_celltype.rds (not in paths.yaml)
# Optional:
#   GENE_MAP_CSV      — two-column (ensembl,symbol) map for SN namespace conversion
#                       defaults to CFG_GENE_SYMBOL_TO_ENSEMBL if set
#
# Pattern source: run/run_step_04_array.sh (SLURM header + _common.sh pattern)

set -euo pipefail

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
exec >> "${HPC_LOGS_DIR}/step_04_train_${SLURM_JOB_ID}.out" 2>> "${HPC_LOGS_DIR}/step_04_train_${SLURM_JOB_ID}.err"
SC_REF_DIR="${SC_REFERENCE_DIR:-}"
SN_REF="${SN_REFERENCE_RDS:-}"
GENE_MAP="${GENE_MAP_CSV:-${CFG_GENE_SYMBOL_TO_ENSEMBL:-}}"

if [ -z "$SC_REF_DIR" ] || [ -z "$SN_REF" ]; then
  echo "FATAL: SC_REFERENCE_DIR and SN_REFERENCE_RDS must be set before running." >&2
  echo "  These are Kumar 2023 reference paths; set them via sbatch --export or shell vars." >&2
  exit 1
fi

SCRIPT="${CFG_SCRIPTS_DIR}/step_04_label_transfer/03a_train_singler_ref.R"
OUTPUT_DIR="${CFG_PREPROCESSING_STEP_04}/models/singler_models"

echo "Step:       step_04_train_singler"
echo "SC ref dir: ${SC_REF_DIR}"
echo "SN ref:     ${SN_REF}"
echo "Output:     ${OUTPUT_DIR}"

GENE_MAP_FLAG=""
if [ -n "$GENE_MAP" ]; then
  GENE_MAP_FLAG="--gene_map ${GENE_MAP}"
fi

run_r_singularity "${SCRIPT}" \
    --sc_ref_dir "${SC_REF_DIR}" \
    --sn_ref     "${SN_REF}" \
    --output_dir "${OUTPUT_DIR}" \
    ${GENE_MAP_FLAG} \
    ${DRY_RUN_FLAG}

echo "Exit: $?, End: $(date)"
