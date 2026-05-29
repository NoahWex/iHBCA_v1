#!/bin/bash
#SBATCH --job-name=joint_markers
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/03_cluster_alignment/joint_markers_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/03_cluster_alignment/joint_markers_%j.err

# Phase 3 — presto::wilcoxauc markers per joint Leiden cluster, for three
# expression sources:
#   (1) Xenium nuclear panel   (Xenium cells, 280 genes)
#   (2) FLEX on Xenium panel   (FLEX cells, 280 genes)
#   (3) FLEX full transcriptome (FLEX cells, ~18K genes)

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

JOB_TMP=/tmp/joint_markers_${SLURM_JOB_ID}
mkdir -p "${JOB_TMP}"
trap 'rm -rf "${JOB_TMP}"' EXIT

module load singularity
singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --bind "${R_LIBS_USER}:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${R_CONTAINER}" \
    Rscript "${PROJECT_ROOT}/pipeline/scripts/03_cluster_alignment/04_joint_markers.R" \
        --joint-leiden "${JOINT_LEIDEN}" \
        --joint-obs    "${JOINT_OBS}" \
        --leiden-res   1.0 \
        --nuc-bundle   "${POOLED_NUCLEAR}" \
        --flex-dir     "${UPSTREAM_FLEX_PREP}/outputs/integration_intermediate" \
        --xenium-panel "${XENIUM_PANEL_GENES}" \
        --out-dir      "${PROJECT_ROOT}/pipeline/outputs/annotations"

echo "Done: $(date)"
