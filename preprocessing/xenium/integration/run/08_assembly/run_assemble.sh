#!/bin/bash
#
# 4b-joint phase 8: assemble master cell_annotations.csv (L0.5 cohort table).
#
# Joins joint embedding (joint_v6) + L0.5 labels (phase 4) + three-axis filter
# (phase 6) + NMP scores (phase 5) + clinical metadata into the master per-cell
# cohort table. ~1.3M rows (FLEX + Xenium cells).
#
# Also writes the provenance manifest (manifest.json) recording version, cohort,
# and row counts.
#
# Tier 3: 4 CPUs / 32G / 1h.

#SBATCH --job-name=xen_assemble
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=.slurm_stubs/assemble_%j.out
#SBATCH --error=.slurm_stubs/assemble_%j.err

if [ -n "${SLURM_JOB_ID:-}" ]; then SCRIPT_DIR="$(pwd)"; else SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"; fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
CONTAINER_TYPE=python_spatial_2025Q4E2
source "${PIPELINE_ROOT}/run/_common.sh"

mkdir -p "${CFG_HPC_LOGS_DIR}"
exec >> "${CFG_HPC_LOGS_DIR}/assemble_${SLURM_JOB_ID}.out" 2>> "${CFG_HPC_LOGS_DIR}/assemble_${SLURM_JOB_ID}.err"

# 4b-joint root (post-promotion: publication/preprocessing/xenium/integration)
# Sub-track scripts read paths via `os.path.join(args.project_root, "config/paths.yaml")`.
SUBTRACK_ROOT="${PIPELINE_ROOT}/4b-joint"

# Step 8.1: assemble cell_annotations
run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/08_assembly/01_assemble_annotations.py" \
        --project-root "${SUBTRACK_ROOT}"

# Step 8.6: write provenance manifest
run_python_singularity \
    "${PIPELINE_ROOT}/4b-joint/scripts/08_assembly/06_write_manifest.py" \
        --project-root "${SUBTRACK_ROOT}"

echo "Done assemble: $(date)"
