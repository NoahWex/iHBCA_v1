#!/bin/bash
#
# 4b-joint phase 4 (step c): final L0.5 arbitration across kNN / SingleR / FICTURE.
#
# Combines three independent L0.5 prediction sources (kNN vote on joint latent,
# SingleR pseudobulk reference, FICTURE topic-based). Produces single
# l0p5_xenium.csv with per-cell label + arbitration evidence (which source
# won, confidence margin).
#
# Tier 2: 4 CPUs / 16G / 1h.

#SBATCH --job-name=xen_p4_arb
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=.slurm_stubs/p4_arb_%j.out
#SBATCH --error=.slurm_stubs/p4_arb_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/p4_arb_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/p4_arb_${SLURM_JOB_ID}.err"

KNN_VOTE="$(dirname "${CFG_XENIUM_L0P5}")/xenium_knn_vote_clean.csv"
LEGACY_ARB="${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/../legacy_arbitration.csv"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/04_label_transfer/03_final_arbitration.py" \
        --legacy-arb "${LEGACY_ARB}" \
        --knn-vote   "${KNN_VOTE}" \
        --out-csv    "${CFG_XENIUM_L0P5}"

echo "Done arbitration: $(date)"
