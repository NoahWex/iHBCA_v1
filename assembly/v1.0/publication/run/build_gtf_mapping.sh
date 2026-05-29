#!/bin/bash
#SBATCH --job-name=gtf_gene_map
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
REPO_ROOT="${PROJECT_IHBCAV1_UPLOAD}"
SCRIPTS="${REPO_ROOT}/publication/scripts"
MAPPINGS="${REPO_ROOT}/publication/mappings"
INTERMEDIATES="${REPO_ROOT}/publication/outputs/source_datasets/intermediates"
GTF="${SOURCE_SPACERANGER_REF}genes/genes.gtf"
CONTAINER="${CONTAINER_PYTHON_SPATIAL_2025Q2}"

echo "============================================="
echo "Build comprehensive GTF-based gene mapping"
echo "============================================="

module load singularity

singularity exec \
    --no-mount bind-paths \
    ${BIND_MOUNTS} \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    "$CONTAINER" \
    python "${SCRIPTS}/build_gtf_gene_mapping.py" \
        --gtf "${GTF}" \
        --output "${MAPPINGS}/gene_symbol_to_ensembl_full.tsv" \
        --existing-mapping "${MAPPINGS}/gene_symbol_to_ensembl.tsv" \
        --hgnc "${MAPPINGS}/hgnc_complete_set.txt" \
        --test-features \
            "${INTERMEDIATES}/murrow/features.tsv.gz" \
            "${INTERMEDIATES}/nee/features.tsv.gz" \
            "${INTERMEDIATES}/pal_norm_epi/features.tsv.gz" \
            "${INTERMEDIATES}/gray/features.tsv.gz" \
            "${INTERMEDIATES}/kumar/features.tsv.gz" \
            "${INTERMEDIATES}/twigger/features.tsv.gz" \
            "${INTERMEDIATES}/reed/features.tsv.gz"

echo ""
echo "COMPLETE"
