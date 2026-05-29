#!/bin/bash
#SBATCH --job-name=xen_panels
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/08_assembly/panels_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium/pipeline/logs/08_assembly/panels_%j.err

# Phase 8 — render diagnostic UMAP panels (platform, patient, compartment,
# L0.5, position, clinical, NMP, UOQ, compartment bg/fg). Reads pipeline/
# outputs/embedding/ + annotations/. Run after run_assemble.sh.

set -euo pipefail
PROJECT_ROOT=/share/crsp/lab/dalawson/nwechter/Spatial_HBCA_Xenium
source "${PROJECT_ROOT}/pipeline/config/load_paths.sh"

PIPE=${PROJECT_ROOT}/pipeline/outputs
PREV=${PIPE}/previews

JOB_TMP=/tmp/xen_panels_${SLURM_JOB_ID}
mkdir -p "${JOB_TMP}" "${PREV}"
trap 'rm -rf "${JOB_TMP}"' EXIT

module load singularity

# Standard UMAP panels
singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/pipeline/scripts/08_assembly/02_render_umap_panels.py" \
        --umap  "${PIPE}/embedding/joint_umap.csv" \
        --obs   "${PIPE}/embedding/joint_obs.csv" \
        --l0p5  "${PIPE}/annotations/l0p5_xenium.csv" \
        --nmp   "${PIPE}/annotations/three_axis_filter.csv" \
        --out-dir "${PREV}"

# UOQ panel
singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/pipeline/scripts/08_assembly/03_render_uoq_umap.py" \
        --umap  "${PIPE}/embedding/joint_umap.csv" \
        --obs   "${PIPE}/embedding/joint_obs.csv" \
        --out-dir "${PREV}"

# Compartment bg/fg on joint UMAP
singularity exec \
    --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
    --bind /pub/nwechter:/pub/nwechter:rw \
    --bind /tmp:/tmp:rw \
    --env "PYTHONUSERBASE=${PYTHONUSERBASE}" \
    --env "MPLCONFIGDIR=${JOB_TMP}/mpl" \
    "${CONTAINER}" \
    python3 "${PROJECT_ROOT}/pipeline/scripts/08_assembly/04_render_compartment_umaps.py" \
        --umap  "${PIPE}/embedding/joint_umap.csv" \
        --obs   "${PIPE}/embedding/joint_obs.csv" \
        --l0p5  "${PIPE}/annotations/l0p5_xenium.csv" \
        --out-dir "${PREV}"

echo "Done: $(date)"
