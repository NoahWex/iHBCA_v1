#!/bin/bash
#
# 4c-annotation step 9.2: per-compartment Leiden + markers + diagnostic UMAPs.
#
# Per-compartment Leiden clustering at multiple resolutions (0.1 / 0.3 / 0.5 /
# 0.7 / 1.0 / 1.5), Wilcoxon marker calling, and diagnostic UMAP panels for
# review during cascade YAML curation.
#
# Tier 4: A100 GPU / 8 CPUs / 64G / 2h per task (RAPIDS on GPU for Leiden).

#SBATCH --job-name=xen_comp_analysis
#SBATCH --partition=free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --gres=gpu:A100:1
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=.slurm_stubs/comp_analysis_%A_%a.out
#SBATCH --error=.slurm_stubs/comp_analysis_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/comp_analysis_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/comp_analysis_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
COMP_LC=$(echo "$COMPARTMENT" | tr '[:upper:]' '[:lower:]')
COMP_DIR="${CFG_XENIUM_COMPARTMENT_DIR}/${COMP_LC}"
mkdir -p "${COMP_DIR}/panels"

# FLEX-side joint Leiden (from phase 03 cluster alignment)
FLEX_LEIDEN_JOINT="$(dirname "${CFG_XENIUM_JOINT_LATENT}")/flex_leiden_joint.csv"

run_python_singularity \
    "${PIPELINE_ROOT}/4c-annotation/scripts/02_compartment_analysis.py" \
        --compartment "${COMPARTMENT}" \
        --latent      "${COMP_DIR}/joint_latent.csv" \
        --obs         "${COMP_DIR}/joint_obs.csv" \
        --umap        "${COMP_DIR}/joint_umap.csv" \
        --xenium-bundle "${CFG_XENIUM_POOLED_NUCLEAR}" \
        --l0p5        "${CFG_XENIUM_L0P5}" \
        --flex-leiden "${FLEX_LEIDEN_JOINT}" \
        --flex-dir    "${CFG_FLEX_CONCORD_DIR}" \
        --flex-qc     "${CFG_FLEX_QC}" \
        --nmp         "${CFG_XENIUM_NMP_SCORED}" \
        --out-dir     "${COMP_DIR}/panels" \
        -k 30

echo "Done comp_analysis ${COMPARTMENT}: $(date)"
