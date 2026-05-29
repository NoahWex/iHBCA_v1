#!/usr/bin/env Rscript

# ##########################################################################
# ##                                                                      ##
# ##   STEP 04: LABEL TRANSFER -- SUPERSEDED FOR CANONICAL PIPELINE FLOW  ##
# ##                                                                      ##
# ##   This script is PORTED FOR PROVENANCE ONLY.                         ##
# ##                                                                      ##
# ##   The canonical L1 label transfer path is:                           ##
# ##       scripts/25b_singler_l1.Rmd                                     ##
# ##   (SingleR against the rebuilt FLEX-panel L1 pseudobulk reference).  ##
# ##                                                                      ##
# ##   Do NOT include this step in canonical dry-run matrices and do NOT  ##
# ##   add Kumar 2023 HBCA reference paths to paths.yaml. If you need to  ##
# ##   reproduce the historical step-04 outputs for debugging step 14     ##
# ##   compartment classification, fill in the SC/SN reference paths      ##
# ##   below and run with --dry-run first to confirm everything resolves. ##
# ##                                                                      ##
# ##########################################################################
#
# Purpose (historical)
#   Dual-reference label transfer from Kumar 2023 HBCA SC and SN atlases
#   onto each sample's cells via Seurat TransferData, followed by a
#   rule-based reconciliation of SC vs SN labels using a mapping file.
#   The output was a per-cell HBCATransferredLabels.Kumar_2023 column
#   used by step 14 compartment classification.
#
# Why superseded
#   The canonical flow now uses SingleR against a pseudobulk L1 reference
#   built directly from the iHBCA v1 atlas (the same atlas this dataset
#   is being integrated into), which is more consistent and requires no
#   external Kumar objects. See coordination/reports/migration_audit/
#   L1_spec.md for the decision trail.
#
# Note: Kumar 2023 reference paths are CLI-required arguments, not stored
# in paths.yaml. Pass them per-invocation via the wrapper's --kumar-* flags.
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
})

parser <- ArgumentParser(description = "Step 04: Label Transfer (PROVENANCE ONLY)")
parser$add_argument("--sample-id", type = "character", default = NULL)
parser$add_argument("--sample-index", type = "integer", default = NULL)
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
parser$add_argument("--central-cell-status", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""))
parser$add_argument("--step-01-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_01", ""))
parser$add_argument("--step-02-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_02", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_04", ""))
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS",
                                         file.path(Sys.getenv("CFG_PROJECT_ROOT", ""),
                                                   "config", "module_configs.yaml")))
parser$add_argument("--sc-reference", type = "character", default = "",
                    help = "Historical Kumar 2023 SC reference RDS (intentionally not in paths.yaml).")
parser$add_argument("--sn-reference", type = "character", default = "",
                    help = "Historical Kumar 2023 SN reference RDS (intentionally not in paths.yaml).")
parser$add_argument("--label-mapping", type = "character", default = "",
                    help = "Historical SC/SN label reconciliation map.")
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

cat("##\n")
cat("## WARNING: step_04_label_transfer is SUPERSEDED by 25b_singler_l1.Rmd\n")
cat("## This wrapper is kept for historical reproduction only.\n")
cat("##\n")

log_event <- function(action, value = NULL, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(value)) cat(sprintf("[%s] %s step_04 %s\n", ts, level, action))
  else cat(sprintf("[%s] %s step_04 %s=%s\n", ts, level, action, as.character(value)))
}
fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, level = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("dry_run", args$dry_run)

# Resolve sample
if (identical(args$raw_data_manifest, "")) fail_loud("raw_data_manifest is empty")
manifest <- read.table(args$raw_data_manifest, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
if (is.null(args$sample_id) && is.null(args$sample_index)) {
  fail_loud("must provide --sample-id or --sample-index")
}
if (!is.null(args$sample_index)) {
  idx <- args$sample_index + 1L
  if (idx < 1L || idx > nrow(manifest)) {
    log_event("graceful_skip", sprintf("index %d out of manifest range", args$sample_index),
              level = "WARN")
    quit(save = "no", status = 0)
  }
  row <- manifest[idx, ]
} else {
  row <- manifest[manifest$sample_id == args$sample_id, ]
  if (nrow(row) != 1L) fail_loud(sprintf("sample-id %s not unique", args$sample_id))
}
sample_id <- row$sample_id
patient_id <- row$patient_id
position_id <- row$position_id
h5_path <- row$h5_path

# Historical references are required for execution but NOT for dry-run.
# Dry-run still validates the CLI surface and output planning so the
# provenance port is at least shape-valid.
output_base <- args$output_root
if (identical(output_base, "")) fail_loud("output-root is empty (CFG_PREPROCESSING_STEP_04).")

cell_metadata_path <- file.path(output_base, "cell_metadata", paste0(sample_id, "_labels.csv"))
sample_manifest_path <- file.path(output_base, "sample_manifests", paste0(sample_id, ".yaml"))
log_path <- file.path(output_base, "logs", paste0(sample_id, ".log"))
viz_dir <- file.path(output_base, "visualizations")
for (d in c(dirname(cell_metadata_path), dirname(sample_manifest_path),
            dirname(log_path), viz_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_04_label_transfer (SUPERSEDED)\n")
  cat("NOTE: this step is excluded from the canonical V1 dry-run matrix.\n")
  cat(sprintf("Sample: %s\n", sample_id))
  cat(sprintf("Output planned: %s\n", cell_metadata_path))
  cat("VALIDATION PASSED (provenance port)\n")
  quit(save = "no", status = 0)
}

# Real execution requires the Kumar references to be passed in explicitly
if (!nzchar(args$sc_reference) || !nzchar(args$sn_reference) ||
    !nzchar(args$label_mapping)) {
  fail_loud(
    "step 04 execution requires --sc-reference, --sn-reference, and --label-mapping",
    "this step is superseded; references are intentionally not in paths.yaml."
  )
}
if (!file.exists(args$sc_reference)) fail_loud(sprintf("SC reference not found: %s", args$sc_reference))
if (!file.exists(args$sn_reference)) fail_loud(sprintf("SN reference not found: %s", args$sn_reference))
if (!file.exists(args$label_mapping)) fail_loud(sprintf("label mapping not found: %s", args$label_mapping))
if (!file.exists(h5_path)) fail_loud(sprintf("raw H5 not found: %s", h5_path))

# Load modules
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
SCRIPT_DIR <- dirname(normalizePath(script_path, mustWork = FALSE))
for (f in c("load_query_object.R", "sc_label_transfer.R", "sn_label_transfer.R",
            "reconcile_labels.R", "visualizations.R", "quality_metrics.R")) {
  source(file.path(SCRIPT_DIR, "source", f))
}

# Load module configs to get transfer params (if present)
mod_cfg <- yaml::read_yaml(args$module_configs)
params <- mod_cfg$modules$step_04_label_transfer$algorithm_params
if (is.null(params)) {
  params <- list(
    transfer_dims = 30, query_assay = "RNA",
    sc_reference_assay = "RNA", sc_label_column = "CellType",
    sn_reference_assay = "RNA", sn_label_column = "CellType",
    use_step02_clusters = TRUE, reconciliation_strategy = "use_mapping_file",
    random_seed = 42
  )
}

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
log_event("execute_start")

step02_meta_path <- file.path(args$step_02_root, "cell_metadata",
                              paste0(sample_id, "_metadata.csv"))
step01_vfs_path <- file.path(args$step_01_root, "features",
                             paste0(sample_id, "_vfs.txt"))
if (!file.exists(step02_meta_path)) fail_loud(sprintf("step-02 metadata missing: %s", step02_meta_path))
if (!file.exists(step01_vfs_path)) fail_loud(sprintf("step-01 VFs missing: %s", step01_vfs_path))

query_obj <- load_query_object(
  h5_path = h5_path,
  sample_id = sample_id,
  central_manifest_path = args$central_cell_status,
  step02_metadata_path = step02_meta_path,
  step04_vfs_path = step01_vfs_path,
  min_cells = 2,
  config = params
)
gc(verbose = FALSE)

query_obj <- perform_sc_label_transfer(query_obj = query_obj,
                                       sc_ref_path = args$sc_reference,
                                       config = params)
gc(verbose = FALSE)

query_obj <- perform_sn_label_transfer(query_obj = query_obj,
                                       sn_ref_path = args$sn_reference,
                                       config = params)
gc(verbose = FALSE)

query_obj <- reconcile_labels(query_obj = query_obj,
                              mapping_path = args$label_mapping,
                              config = params)

metrics <- calculate_transfer_metrics(query_obj)
viz_paths <- generate_label_plots(query_obj = query_obj,
                                  output_dir = viz_dir,
                                  sample_id = sample_id)

metadata_output <- query_obj@meta.data[, c(
  "predicted.celltype.SC", "prediction.score.SC",
  "predicted.celltype.SN", "prediction.score.SN",
  "HBCATransferredLabels.Kumar_2023", "cluster_coarse"
)]
tmp <- paste0(cell_metadata_path, ".tmp.", Sys.getpid())
write.csv(metadata_output, tmp, row.names = TRUE)
file.rename(tmp, cell_metadata_path)
log_event("wrote", cell_metadata_path)

sample_data <- list(
  sample_info = list(
    patient_id = patient_id, position_id = position_id, sample_id = sample_id,
    position_status = "step_04_complete_superseded"
  ),
  step_04_label_transfer = list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    status = "pass",
    outputs = list(cell_metadata = cell_metadata_path, log = log_path),
    summary = list(total_cells = metrics$total_cells),
    note = "SUPERSEDED by 25b_singler_l1.Rmd -- not used in canonical flow."
  )
)
tmp_mfst <- paste0(sample_manifest_path, ".tmp.", Sys.getpid())
yaml::write_yaml(sample_data, tmp_mfst)
file.rename(tmp_mfst, sample_manifest_path)
log_event("execute_end")
quit(save = "no", status = 0)
