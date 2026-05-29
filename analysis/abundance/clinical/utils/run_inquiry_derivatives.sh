#!/bin/bash
#SBATCH --job-name=inq_derivatives
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_publication/logs/inq_derivatives_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/iHBCA_publication/logs/inq_derivatives_%j.err
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
set -euo pipefail
#
# Build long-form derivatives for both V1 clinical DA inquiries:
#   parity_findings   (8 parity contrasts)
#   risk_main_effects (3 risk-main contrasts + 1 sensitivity)

PUB="/share/crsp/lab/dalawson/nwechter/iHBCA_publication"
PRODUCER="$PUB/publication/analysis/abundance/clinical/utils/build_inquiry_derivatives.R"

mkdir -p "$PUB/logs"
module load singularity

R_CONTAINER="/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif"
R_LIBS_USER_PATH="/pub/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library"

for INQUIRY in parity_findings risk_main_effects; do
  echo "========================================"
  echo "Inquiry: $INQUIRY"
  echo "========================================"
  singularity exec --no-home \
      --bind "$R_LIBS_USER_PATH:/home/jovyan/R/library:ro" \
      --bind /share/crsp/lab/dalawson:/share/crsp/lab/dalawson:rw \
      --bind /dfs7:/dfs7:ro \
      --bind /pub/nwechter:/pub/nwechter:ro \
      --bind /tmp:/tmp:rw \
      --env "R_LIBS_USER=/home/jovyan/R/library" --env "HOME=/home/jovyan" \
      --env "TMPDIR=/tmp" \
      "$R_CONTAINER" Rscript "$PRODUCER" \
          --project-root "$PUB" \
          --inquiry "$INQUIRY"
done

echo "=== Done ==="
