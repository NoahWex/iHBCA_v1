#!/bin/bash
#
# 4b-joint phase 5: dual nuclear/cytoplasmic NMP scoring.
#
# Negative Marker Probability per Xenium cell, computed independently for
# nuclear and cytoplasmic count layers. Cells with high NMP (many transcripts
# from genes that should not co-occur in the assigned cell type) are flagged
# as likely segmentation artifacts.
#
# Tier 3: 4 CPUs / 24G / 2h.

#SBATCH --job-name=xen_p5_nmp
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=.slurm_stubs/p5_nmp_%j.out
#SBATCH --error=.slurm_stubs/p5_nmp_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/p5_nmp_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/p5_nmp_${SLURM_JOB_ID}.err"

# Negative-marker tables from legacy FLEX L0.5 NMP work (intake C8)
LEGACY_NMP="${CFG_XENIUM_LEGACY_PANEL_INTERSECTION}/../../../02_flex_l0p5_nmp/outputs"

run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/05_nmp_scoring/01_xenium_nmp_dual.py" \
        --nuclear-bundle "${CFG_XENIUM_POOLED_NUCLEAR}" \
        --cyto-bundle    "${CFG_XENIUM_POOLED_CYTOPLASMIC}" \
        --l0p5-xenium    "${CFG_XENIUM_L0P5}" \
        --neg-markers    "${LEGACY_NMP}/flex_negative_markers_l0p5.yaml" \
        --families       "${LEGACY_NMP}/marker_families.yaml" \
        --out-csv        "${CFG_XENIUM_NMP_SCORED}"

echo "Done phase5_nmp: $(date)"
