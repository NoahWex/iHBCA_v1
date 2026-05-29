#!/bin/bash
#
# 4a-counts step 1: per-position boundary counts.
# Reads xenium-ranger transcripts.parquet + cells.parquet for each sample in
# xenium_manifest.tsv and emits three count matrices per sample (whole,
# nuclear, cytoplasmic) split on overlaps_nucleus, plus per-cell QC metrics.
# Tier 2: 4 CPUs / 32G / 1h per task. Output: per xenium_id, counts_*.mtx.gz
# + obs.csv + cells.tsv + genes.tsv. Methodological core: load_counts.py.

#SBATCH --job-name=xen_03a_counts
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --array=1-65
#SBATCH --output=.slurm_stubs/03a_counts_%A_%a.out
#SBATCH --error=.slurm_stubs/03a_counts_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_scvi_gpu_2025Q3
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/03a_counts_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/03a_counts_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

MANIFEST="${CFG_XENIUM_MANIFEST}"
LINE=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR==i+1' "${MANIFEST}")
XENIUM_ID=$(echo "${LINE}" | cut -f1)
PATIENT_ID=$(echo "${LINE}" | cut -f2)
TX_PATH=$(echo "${LINE}" | cut -f3)
SAMPLE_PATH=$(dirname "${TX_PATH}")

echo "task=${SLURM_ARRAY_TASK_ID} xenium_id=${XENIUM_ID} sample=${SAMPLE_PATH}"
mkdir -p "${CFG_XENIUM_COUNTS_PER_POSITION}"

run_python_singularity \
    "${PIPELINE_ROOT}/4a-counts/scripts/03a_build_boundary_counts.py" \
        --sample-path "${SAMPLE_PATH}" \
        --xenium-id   "${XENIUM_ID}" \
        --patient-id  "${PATIENT_ID}" \
        --out-dir     "${CFG_XENIUM_COUNTS_PER_POSITION}"

echo "Done task=${SLURM_ARRAY_TASK_ID}: $(date)"
