#!/bin/bash
#SBATCH --job-name=enrich_prereqs
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Phase 0 prerequisites: parse GTF + re-stage L1 metadata

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "============================================="
echo "Phase 0: Prerequisites"
echo "  Date: $(date)"
echo "  Node: $(hostname)"
echo "============================================="

module load singularity

# 0b. Parse GTF
echo ""
echo "--- Step 1: Parse GENCODE v24 GTF ---"
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    "$CONTAINER" \
    python3 "${SCRIPTS}/parse_gtf.py" \
        --gtf "${REPO_ROOT}/publication/mappings/gencode.v24.annotation.gtf" \
        --output "${REPO_ROOT}/publication/mappings/gencode_v24_gene_annotations.tsv"

echo ""
echo "Verifying GTF parse output:"
wc -l "${REPO_ROOT}/publication/mappings/gencode_v24_gene_annotations.tsv"
head -3 "${REPO_ROOT}/publication/mappings/gencode_v24_gene_annotations.tsv"

# 0c. Re-stage L1 metadata
echo ""
echo "--- Step 2: Re-stage L1 metadata (18 -> 21 obs columns) ---"
singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    "$CONTAINER" \
    python3 "${SCRIPTS}/stage_l1_metadata.py" \
        --repo-root "$REPO_ROOT"

echo ""
echo "Verifying L1 output:"
head -1 "${REPO_ROOT}/publication/config/metadata_stages/L1_harmonized_donor.csv" | tr ',' '\n' | wc -l
echo "columns in L1 CSV"

echo ""
echo "============================================="
echo "Phase 0 COMPLETE"
echo "Date: $(date)"
echo "============================================="
