#!/bin/bash
# =============================================================================
# submit_metadata_pipeline.sh — Chain metadata harmonization pipeline
# =============================================================================
# Submits: build_donor_metadata → build_harmonized → map_cells → coverage+report
# Each step depends on the previous completing successfully.
#
# Like a relay race, but the baton is a CSV and nobody's running.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../../config/load_paths.sh"

RUN_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}/run"
LOG_DIR="${SOURCE_IHBCAV1_HARMONIZATION%/}/logs"

echo "=== Metadata Harmonization Pipeline ==="
echo "$(date)"
echo ""

# Step 1: Build unified donor metadata (loads reference CSVs, merges)
JOB1=$(sbatch --parsable "${RUN_DIR}/build_donor_metadata.sh")
echo "Step 1 - build_donor_metadata:    job ${JOB1}"

# Step 2: Build harmonized donor metadata (applies harmonization logic)
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} "${RUN_DIR}/run_build_harmonized.sh")
echo "Step 2 - build_harmonized_donors: job ${JOB2} (after ${JOB1})"

# Step 3: Map donor metadata to cells
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} "${RUN_DIR}/run_map_cells.sh")
echo "Step 3 - map_donor_to_cells:      job ${JOB3} (after ${JOB2})"

# Step 4a+4b: Coverage matrix + HTML report (independent, both depend on step 3)
JOB4a=$(sbatch --parsable --dependency=afterok:${JOB3} "${RUN_DIR}/run_coverage_matrix.sh")
JOB4b=$(sbatch --parsable --dependency=afterok:${JOB3} "${RUN_DIR}/run_harmonization_report.sh")
echo "Step 4a - coverage_matrix:        job ${JOB4a} (after ${JOB3})"
echo "Step 4b - harmonization_report:   job ${JOB4b} (after ${JOB3})"

echo ""
echo "=== Pipeline Submitted ==="
echo "Monitor: squeue -j ${JOB1},${JOB2},${JOB3},${JOB4a},${JOB4b}"
echo ""
echo "Tail step 1:"
echo "  while [ ! -f ${LOG_DIR}/*${JOB1}* ]; do sleep 1; done && tail -f ${LOG_DIR}/*${JOB1}*"
