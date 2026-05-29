#!/bin/bash
#SBATCH --job-name=tier1_donor
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00

# =============================================================================
# Generate Tier 1 Donor Metadata CSVs for HCA Tracker entry sheets
#
# Produces 7 CSVs (one per study) in:
#   publication/outputs/entry_sheets/tier1_donor/
#
# Each CSV contains 8 Tier 1 Donor Metadata columns per CxG schema 5.3.2:
#   donor_id, dataset_id, organism_ontology_term_id, manner_of_death,
#   disease_ontology_term_id, sex_ontology_term_id,
#   development_stage_ontology_term_id, self_reported_ethnicity_ontology_term_id
#
# Expected donor counts: gray=16, kumar=126, murrow=28, nee=22,
#                        twigger=18, reed=55, pal=22  (total=287)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/load_paths.sh"
ROOT="${PROJECT_IHBCAV1_UPLOAD}"

echo "=== Generate Tier 1 Donor Entry Sheets ==="
echo "Date: $(date)"
echo "Host: $(hostname)"
echo ""

module load python/3.10.2

python "$ROOT/publication/scripts/generate_tier1_donor_sheet.py" \
    --project-root "$ROOT"

echo ""
echo "Output: $ROOT/publication/outputs/entry_sheets/tier1_donor/"
echo "Done: $(date)"
