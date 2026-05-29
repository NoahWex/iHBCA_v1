#!/bin/bash
#SBATCH --job-name=l1_confusogram
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --array=0-3
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=.slurm_stubs/l1_confusogram_%A_%a.out
#SBATCH --error=.slurm_stubs/l1_confusogram_%A_%a.err

# L1 confusogram rendering — three matrices per task:
#   cluster_annotation -> SingleR L1 raw -> Harmonized L1.0
#
# Container:   r_spatial
# Tier 4:      8cpu/48G/2h — matches the pattern tier; actual workload is
#              pheatmap rendering on small xtabs (~seconds), Seurat load
#              dominates for per-compartment tasks
# Pattern src: run_compartment_integration_report.sh (3-task array) plus an
#              extra task for the 'all' concatenated view
#
# Tasks:
#   0  Epithelial  — per-compartment view (reads Seurat for cluster_anno)
#   1  Stromal     — per-compartment view
#   2  Immune      — per-compartment view
#   3  all         — uses the concatenated harmonized + contamination CSV
#                    direct (no per-compartment Seurat), consumes ALL SingleR
#                    prediction files merged upstream via an ad hoc cat below.

if [ -n "${SLURM_JOB_ID:-}" ]; then
    SCRIPT_DIR="$(pwd)"
else
    SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
fi
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONTAINER_TYPE=r_spatial
source "${PIPELINE_ROOT}/run/_common.sh"


# Redirect stdout/stderr to canonical publication logs location.
# SBATCH directives above land tiny stubs in .slurm_stubs/; real output goes
# to the publication-tier ${HPC_LOGS_DIR}/ which is resolved by _common.sh.
mkdir -p "${HPC_LOGS_DIR}"
exec >> "${HPC_LOGS_DIR}/l1_confusogram_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.out" 2>> "${HPC_LOGS_DIR}/l1_confusogram_${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID}}_${SLURM_ARRAY_TASK_ID:-0}.err"
COMPARTMENTS=(Epithelial Stromal Immune all)
COMPARTMENT=${COMPARTMENTS[$SLURM_ARRAY_TASK_ID]}

RMD="${CFG_SCRIPTS_DIR}/step_19a_l1_annotation/09e_l1_confusogram.Rmd"
L1_ROOT="${CFG_PROJECT_ROOT}/outputs/annotation/l1"

if [ "${DRY_RUN}" = "1" ]; then
    DRY_R="TRUE"
else
    DRY_R="FALSE"
fi

if [ "${COMPARTMENT}" = "all" ]; then
    # Concatenate the 3 per-compartment singler CSVs into a temp merged CSV
    # so the Rmd can consume the full 270K rows in one pass. Use the Rmd's
    # cluster_annotation_csv branch (reads cell_compartments.csv directly).
    MERGED_SINGLER="${JOB_TEMP_DIR}/singler_l1_predictions_all.csv"

    HARMONIZED_CSV="${L1_ROOT}/all/harmonized_l1_labels.csv"
    CLUSTER_ANNO_CSV="${CFG_UPSTREAM_CONTAMINATION}"
    SEURAT_RDS=""
    OUTDIR="${L1_ROOT}/all"
    OUTPUT_FILE="confusogram.html"

    if [ "${DRY_RUN}" = "1" ]; then
        # In dry run, skip the merge — the Rmd's dry-run preamble will check
        # required files; create an empty placeholder so its existence check
        # passes. The preamble exits before the file is read.
        touch "${MERGED_SINGLER}"
    else
        mkdir -p "$(dirname "${MERGED_SINGLER}")"
        HEAD_CSV="${L1_ROOT}/Epithelial/singler_l1_predictions_Epithelial.csv"
        head -n 1 "${HEAD_CSV}" > "${MERGED_SINGLER}"
        for c in Epithelial Stromal Immune; do
            tail -n +2 "${L1_ROOT}/${c}/singler_l1_predictions_${c}.csv" >> "${MERGED_SINGLER}"
        done
        echo "Merged SingleR CSV (rows): $(wc -l < ${MERGED_SINGLER})"
    fi

    SINGLER_CSV="${MERGED_SINGLER}"
else
    SINGLER_CSV="${L1_ROOT}/${COMPARTMENT}/singler_l1_predictions_${COMPARTMENT}.csv"
    HARMONIZED_CSV="${L1_ROOT}/${COMPARTMENT}/harmonized_l1_labels_${COMPARTMENT}.csv"
    SEURAT_RDS="${L1_ROOT}/${COMPARTMENT}/seurat.rds"
    CLUSTER_ANNO_CSV=""
    OUTDIR="${L1_ROOT}/${COMPARTMENT}"
    OUTPUT_FILE="confusogram.html"
fi

mkdir -p "${OUTDIR}"

echo "Task:                 ${COMPARTMENT}"
echo "Rmd:                  ${RMD}"
echo "SingleR CSV:          ${SINGLER_CSV}"
echo "Harmonized CSV:       ${HARMONIZED_CSV}"
echo "Seurat RDS:           ${SEURAT_RDS}"
echo "cluster_anno CSV:     ${CLUSTER_ANNO_CSV}"
echo "Output:               ${OUTDIR}/${OUTPUT_FILE}"
echo "Dry run:              ${DRY_R}"

run_r_singularity -e "rmarkdown::render('${RMD}', output_dir = '${OUTDIR}', output_file = '${OUTPUT_FILE}', intermediates_dir = '${JOB_TEMP_DIR}', knit_root_dir = '${CFG_PROJECT_ROOT}', params = list(compartment = '${COMPARTMENT}', singler_csv = '${SINGLER_CSV}', harmonized_csv = '${HARMONIZED_CSV}', seurat_rds = '${SEURAT_RDS}', cluster_annotation_csv = '${CLUSTER_ANNO_CSV}', outdir = '${OUTDIR}', dry_run = ${DRY_R}))"

echo "Exit: $?, End: $(date)"
