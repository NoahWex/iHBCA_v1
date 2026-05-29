#!/bin/bash
#
# 4a-counts step 4: apply locked Phase 1 threshold gate.
#
# Writes pass_qc_whole = (nFeature_Whole > 5) AND (tx_per_gene_Whole > 1.25)
# in-place on xenium_qc.csv. Cohort pass rate ~93%. Threshold cut splits the
# defensible shoulder (>1.5) and permissive peak-intact cut (>1.0).
#
# Tier 1: 2 CPUs / 4G / 15min.

#SBATCH --job-name=xen_finalize_thr
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=.slurm_stubs/finalize_thr_%j.out
#SBATCH --error=.slurm_stubs/finalize_thr_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/finalize_thr_${SLURM_JOB_ID}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/finalize_thr_${SLURM_JOB_ID}.err"

run_python_singularity \
    "${PIPELINE_ROOT}/4a-counts/scripts/04_finalize_threshold.py" \
        --xqc   "${CFG_XENIUM_QC}" \
        --fqc   "${CFG_FLEX_QC}" \
        --audit "${CFG_COHORT_QC_AUDIT}"

echo "Done finalize_threshold: $(date)"
