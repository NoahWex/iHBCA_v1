#!/bin/bash
#SBATCH --job-name=joint_heatmaps
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/03_cluster_alignment/joint_heatmaps_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/03_cluster_alignment/joint_heatmaps_%j.err

# Phase 3 — render three joint-cluster marker heatmaps from the markers CSVs
# produced by run_joint_markers.sh, into pipeline/outputs/previews/.

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

JOB_TMP=/tmp/joint_heatmaps_${SLURM_JOB_ID}
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
    Rscript "${PROJECT_ROOT}/pipeline/scripts/03_cluster_alignment/05_render_joint_heatmaps.R" \
        --joint-leiden "${JOINT_LEIDEN}" \
        --joint-obs    "${JOINT_OBS}" \
        --leiden-res   1.0 \
        --markers-dir  "${PROJECT_ROOT}/pipeline/outputs/annotations" \
        --nuc-bundle   "${POOLED_NUCLEAR}" \
        --flex-dir     "${UPSTREAM_FLEX_PREP}/outputs/integration_intermediate" \
        --xenium-panel "${XENIUM_PANEL_GENES}" \
        --out-dir      "${PREVIEWS_DIR}" \
        --top-n 5

echo "Done: $(date)"
