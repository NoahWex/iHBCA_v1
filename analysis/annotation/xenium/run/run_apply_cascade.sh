#!/bin/bash
#
# 4c-annotation step 9.4: apply cascade YAML to leiden assignments → cell_annotations.
#
# Reads the hand-curated annotation_cascade YAML for each compartment and
# applies it to leiden_assignments.csv to produce per-compartment
# cell_annotations.csv with L1.5 labels (label, label_short, lineage,
# is_artifact). Cascade YAMLs are version-controlled at config/cascades/.
#
# Tier 1: 2 CPUs / 8G / 15min per task.

#SBATCH --job-name=xen_apply_cascade
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --array=0-2
#SBATCH --output=.slurm_stubs/apply_cascade_%A_%a.out
#SBATCH --error=.slurm_stubs/apply_cascade_%A_%a.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/apply_cascade_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" \
    2>> "${CFG_HPC_LOGS_DIR}/apply_cascade_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"

COMPARTMENTS=(Epithelial Stromal Immune)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}
COMP_LC=$(echo "$COMPARTMENT" | tr '[:upper:]' '[:lower:]')
COMP_DIR="${CFG_XENIUM_COMPARTMENT_DIR}/${COMP_LC}"

run_python_singularity \
    "${PIPELINE_ROOT}/4c-annotation/scripts/05_apply_cascade.py" \
        "${CFG_XENIUM_CASCADES_DIR}/${COMP_LC}_cascade.yaml" \
        "${COMP_DIR}/leiden_assignments.csv" \
        "${COMP_DIR}/cell_annotations.csv"

echo "Done apply_cascade ${COMPARTMENT}: $(date)"
