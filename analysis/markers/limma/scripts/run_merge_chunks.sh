#!/bin/bash
#SBATCH --job-name=limma_merge_chunks
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%x_%j.out
#SBATCH --error=slurm-%x_%j.err

# =============================================================================
# Merge per-chunk limma CSVs produced by run_14b_limma_{epi,imm,str}_chunked.sh
# into canonical limma_markers_{level}.csv + limma_top_markers_{level}.csv per
# (compartment, resolution) combo.
#
# Iterates all 3 compartments × 7 resolutions. Skips (compartment, resolution)
# pairs without chunk files (e.g., if the chunked wrapper was only run for a
# subset of resolutions). Merged outputs are written in-place to:
#   ${PROJECT_PUBLICATION}/publication/analysis/markers/limma/outputs/
#     {compartment}/leiden_{res}/
#       limma_markers_leiden_{res}.csv
#       limma_top_markers_leiden_{res}.csv
#
# No separate "promote" cp step: merged outputs land directly in the
# canonical publication path. Only the 3 T2-feeding combos
# (epi/leiden_1.5, str/leiden_1.0, imm/leiden_1.0) are tracked in
# publication/manifest.yaml; the other 18 mergeable combos are derived
# artifacts producible by re-running this wrapper.
# =============================================================================

set -euo pipefail

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$PWD"
else
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"

OUTPUT_ROOT="${PROJECT_PUBLICATION}/publication/analysis/markers/limma/outputs"
RESOLUTIONS=(0.1 0.2 0.3 0.5 0.8 1.0 1.5)
COMPARTMENTS=(epithelial immune stromal)

module load singularity

R_SING_EXEC="singularity exec --no-home \
  ${BIND_MOUNTS} \
  --bind ${R_LIBS_R_SPATIAL_4_3_3}:/home/jovyan/R/userlib:ro \
  --env R_LIBS_USER=/home/jovyan/R/userlib \
  ${CONTAINER_R_SPATIAL_4_3_3}"

echo "=== Merge limma marker chunks ==="
echo "Start: $(date)"

for COMP in "${COMPARTMENTS[@]}"; do
  for RES in "${RESOLUTIONS[@]}"; do
    LEVEL="leiden_${RES}"
    OUTDIR="${OUTPUT_ROOT}/${COMP}/${LEVEL}"

    if [ ! -d "$OUTDIR" ]; then
      echo "[SKIP] $COMP $LEVEL — output dir not found"
      continue
    fi

    N_CHUNKS=$(ls "${OUTDIR}/limma_markers_${LEVEL}_c"*.csv 2>/dev/null | wc -l)
    if [ "$N_CHUNKS" -eq 0 ]; then
      echo "[SKIP] $COMP $LEVEL — no chunk files found"
      continue
    fi

    echo ""
    echo "--- Merging $COMP $LEVEL ($N_CHUNKS chunk files) ---"
    $R_SING_EXEC Rscript "${SCRIPT_DIR}/merge_14b_limma_chunks.R" \
      --output-dir "$OUTDIR" \
      --level-tag  "$LEVEL" \
      --n-top 20
  done
done

echo ""
echo "=== Merge complete ==="
ls -lh "${OUTPUT_ROOT}"/*/*/limma_markers_leiden_*.csv 2>/dev/null \
  | grep -v '_c[0-9]' | awk '{print $5, $9}'
echo "Done: $(date)"
