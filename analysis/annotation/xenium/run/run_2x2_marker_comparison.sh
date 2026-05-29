#!/bin/bash
#
# 4c-annotation step 9.3: per-compartment 2x2 marker comparison evidence document.
#
# For each compartment, produces a 2x2 panel of marker plots comparing Xenium
# nuclear markers (top-DE in Xenium clusters) and FLEX markers (top-DE in FLEX
# clusters) against both Xenium and FLEX expression. Used during cascade YAML
# curation to identify clusters where Xenium and FLEX agree (canonical L1.5
# label) vs disagree (artifact / unknown). Also produces FLEX-context cluster
# coverage stats.
#
# Tier 4: 8 CPUs / 64G / 4h per task (parallel Wilcoxon DE across resolutions).

#SBATCH --job-name=xen_2x2_marker
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --array=0-2
#SBATCH --output=.slurm_stubs/2x2_marker_%A_%a.out
#SBATCH --error=.slurm_stubs/2x2_marker_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/2x2_marker_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/2x2_marker_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
COMP_LC=$(echo "$COMPARTMENT" | tr '[:upper:]' '[:lower:]')
COMP_DIR="${CFG_XENIUM_COMPARTMENT_DIR}/${COMP_LC}"
PANELS_DIR="${COMP_DIR}/panels2x2"
mkdir -p "${PANELS_DIR}"

FLEX_CLUSTER_SUMMARY="$(dirname "${CFG_XENIUM_JOINT_LATENT}")/cluster_summary_res0.1.csv"
FLEX_LEIDEN_JOINT="$(dirname "${CFG_XENIUM_JOINT_LATENT}")/flex_leiden_joint.csv"

run_python_singularity \
    "${PIPELINE_ROOT}/4c-annotation/scripts/03_2x2_marker_comparison.py" \
        --compartment "${COMPARTMENT}" \
        --latent      "${COMP_DIR}/joint_latent.csv" \
        --obs         "${COMP_DIR}/joint_obs.csv" \
        --umap        "${COMP_DIR}/joint_umap.csv" \
        --nuc-bundle  "${CFG_XENIUM_POOLED_NUCLEAR}" \
        --cyto-bundle "${CFG_XENIUM_POOLED_CYTOPLASMIC}" \
        --flex-dir    "${CFG_FLEX_CONCORD_DIR}" \
        --flex-qc     "${CFG_FLEX_QC}" \
        --l0p5        "${CFG_XENIUM_L0P5}" \
        --flex-cluster-summary "${FLEX_CLUSTER_SUMMARY}" \
        --flex-leiden-joint    "${FLEX_LEIDEN_JOINT}" \
        --nmp         "${CFG_XENIUM_NMP_SCORED}" \
        --out-dir     "${PANELS_DIR}" \
        --html        "${COMP_DIR}/${COMP_LC}_2x2.html" \
        -k 30 --res 0.3 --n-top 4

echo "Done 2x2_marker ${COMPARTMENT}: $(date)"
