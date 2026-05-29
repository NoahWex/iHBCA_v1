#!/bin/bash
#SBATCH --job-name=ihbca_anno_imm_evidence
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=${HPC_LOGS_DIR}/%x_%j.out
#SBATCH --error=${HPC_LOGS_DIR}/%x_%j.err

# =============================================================================
# Per-label evidence report (multi-study, multi-patient, marker-support
# tabulation) for the immune compartment.
#
# Tier 3 (32G / 4 CPU / 1h) per hpc-resource-rules.md.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/load_paths.sh"

COMPARTMENT="imm"
OUT_DIR="${PROJECT_PUBLICATION}/analysis/annotation/${COMPARTMENT}"
STAGE_OUT="${OUT_DIR}/evidence/"

V2_YAML="${SCRIPT_DIR}/../yamls/annotation_v2_${COMPARTMENT}.yaml"
LABELS="${OUT_DIR}/labels.csv"
COUNTS="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}${COMPARTMENT}_counts.npz"
GENE_DATA="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}gene_data.csv"
META="${SOURCE_IHBCAV1_ASSEMBLY_COMPONENTS}${COMPARTMENT}_metadata_enriched.csv"
CONFIG="${SCRIPT_DIR}/../config/compartment_config.yaml"
OUT_CSV="${STAGE_OUT}evidence_report.csv"
SCRIPT="${SCRIPT_DIR}/../scripts/04_compute_evidence.py"

mkdir -p "${STAGE_OUT}" "${PROJECT_PUBLICATION}/logs"

for f in "$V2_YAML" "$LABELS" "$COUNTS" "$GENE_DATA" "$META" "$CONFIG" "$SCRIPT"; do
    [ -f "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done

module load singularity

singularity exec --no-home \
    ${BIND_MOUNTS} \
    --env "PYTHONUSERBASE=${PYTHON_LIBS_SCVI_GPU_2025Q3_USERBASE}" \
    --env "PYTHONUNBUFFERED=1" \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    "${CONTAINER_PYTHON_SCVI_GPU_2025Q3}" \
    python "$SCRIPT" \
        --v2-yaml "$V2_YAML" \
        --labels-csv "$LABELS" \
        --counts-npz "$COUNTS" \
        --gene-data "$GENE_DATA" \
        --metadata-csv "$META" \
        --config "$CONFIG" \
        --out-csv "$OUT_CSV"

echo "=== Complete ==="
