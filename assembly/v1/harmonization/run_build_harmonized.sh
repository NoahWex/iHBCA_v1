#!/bin/bash
#SBATCH --job-name=build_harmonized
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/iHBCA_publication/coordination/daily/build_harmonized_%j.log

module load singularity

SIF=/dfs7/singularity_containers/rcic/JHUB3/Rocky8_jupyter_base_R4.3.3_Spatial.sif
HARMONIZATION=/share/crsp/lab/dalawson/nwechter/iHBCA_publication/publication/assembly/v1/harmonization
UNIFIED=/share/crsp/lab/dalawson/nwechter/iHBCAv1_upload/external_studies/harmonization/outputs/unified_metadata/unified_donor_metadata.csv

singularity exec \
  --bind /share/crsp/lab/dalawson/nwechter:/share/crsp/lab/dalawson/nwechter:rw \
  --bind /dfs7:/dfs7:ro \
  --env NUMBA_CACHE_DIR=/tmp/numba_cache \
  $SIF \
  python3 $HARMONIZATION/build_harmonized_donors.py \
    --input $UNIFIED \
    --output-dir $HARMONIZATION
