#!/usr/bin/env Rscript

# ============================================================================
# Step 05: Post-QC Variable Features - Wrapper
# ============================================================================
# Purpose
#   Re-run BigSur variance-based VF discovery on the post-QC cell set
#   (cells passing all of step 01, step 02, and step 03). The resulting
#   "filtered VFs" are the per-sample feature list consumed by the
#   integration sweep (step 17) and the filtered preview (step 06).
#
# Why this step exists
#   BigSur's variance-based features depend on the per-cell count
#   distribution. Removing debris/doublets changes that distribution
#   enough that baseline-step-01 VFs can contain technically-driven
#   features. Re-running on the filtered cell set produces a more
#   biologically clean feature list AND a retention metric (how much
#   of the baseline VF list survives QC), which is a useful sample-
#   level quality indicator flagged in step 07.
#
# Cell-set-changing? No. Step 05 reads the existing final cell set from
# the central manifest; it only changes the FEATURE set, not the cells.
#
# Inputs
#   - CFG_PREPROCESSING_CENTRAL_CELL_STATUS (requires step 01/02/03 pass columns)
#   - CFG_PREPROCESSING_STEP_01/features/{sample}_vfs.txt (baseline VF list)
#   - Raw H5 from CFG_RAW_DATA_SAMPLE_MANIFEST
#   - Algorithm params from CFG_MODULE_CONFIGS
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_05}/features/filtered_vfs/{sample}_vfs_filtered.txt
#   - ${CFG_PREPROCESSING_STEP_05}/features/feature_ranks/{sample}_feature_ranks.rds
#   - ${CFG_PREPROCESSING_STEP_05}/retention_metrics/{sample}_retention.csv
#   - ${CFG_PREPROCESSING_STEP_05}/normalized_matrices/{sample}_bigsur_normalized.rds
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
})

parser <- ArgumentParser(description = "Step 05: Post-QC Variable Features")
parser$add_argument("--sample-id", type = "character", default = NULL)
parser$add_argument("--sample-index", type = "integer", default = NULL)
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
parser$add_argument("--central-cell-status", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""))
parser$add_argument("--step-01-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_01", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_05", ""))
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS",
                                         file.path(Sys.getenv("CFG_PROJECT_ROOT", ""),
                                                   "config", "module_configs.yaml")))
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

log_event <- function(action, value = NULL, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(value)) cat(sprintf("[%s] %s step_05 %s\n", ts, level, action))
  else cat(sprintf("[%s] %s step_05 %s=%s\n", ts, level, action, as.character(value)))
}
fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, level = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("dry_run", args$dry_run)

for (lbl in c("raw_data_manifest", "central_cell_status", "step_01_root",
              "output_root", "module_configs")) {
  v <- args[[lbl]]
  if (identical(v, "") || is.null(v)) fail_loud(sprintf("%s is empty", lbl))
  log_event(lbl, v)
}

manifest <- read.table(args$raw_data_manifest, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
if (is.null(args$sample_id) && is.null(args$sample_index)) {
  fail_loud("must provide --sample-id or --sample-index")
}
if (!is.null(args$sample_index)) {
  idx <- args$sample_index + 1L
  if (idx < 1L || idx > nrow(manifest)) {
    log_event("graceful_skip", sprintf("index %d out of range", args$sample_index), level = "WARN")
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
log_event("sample_id", sample_id)
if (!file.exists(h5_path)) fail_loud(sprintf("raw H5 not found: %s", h5_path))

vf_path <- file.path(args$step_01_root, "features", paste0(sample_id, "_vfs.txt"))
if (!file.exists(vf_path)) fail_loud(sprintf("step-01 VFs not found: %s", vf_path))
vfs_baseline <- readLines(vf_path)
log_event("n_vfs_baseline", length(vfs_baseline))

mod_cfg <- yaml::read_yaml(args$module_configs)
params <- mod_cfg$modules$step_05_post_qc_vfs$algorithm_params
if (is.null(params)) fail_loud("step_05 params missing in module_configs.yaml")

central_manifest <- read.csv(args$central_cell_status, stringsAsFactors = FALSE)
for (col in c("step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass")) {
  if (!(col %in% colnames(central_manifest))) {
    fail_loud(sprintf("central_cell_status missing %s", col),
              "run step 01/02/03 aggregators before step 05.")
  }
  central_manifest[[col]] <- as.logical(central_manifest[[col]])
}

sample_cells_step01 <- central_manifest[
  central_manifest$sample_id == sample_id & central_manifest$step_01_umi_pass == TRUE,
]
n_cells_baseline <- nrow(sample_cells_step01)
log_event("n_cells_step01", n_cells_baseline)

# Row-count assertion: step-01 sidecar metadata should have the same
# number of rows as the central manifest says it should.
step01_meta_path <- file.path(args$step_01_root, "cell_metadata",
                              paste0(sample_id, "_metadata.csv"))
if (!file.exists(step01_meta_path)) {
  fail_loud(sprintf("step-01 metadata not found: %s", step01_meta_path))
}
step01_meta <- read.csv(step01_meta_path, stringsAsFactors = FALSE)
if (nrow(step01_meta) != n_cells_baseline) {
  # step 01 writes all cells (passing + failing); filter to umi-passing
  n_umi_pass <- sum(step01_meta$flat_threshold_pass)
  if (n_umi_pass != n_cells_baseline) {
    fail_loud(sprintf(
      "row count mismatch: step-01 metadata has %d UMI-pass rows but central manifest has %d",
      n_umi_pass, n_cells_baseline),
      "this usually means the step-01 aggregator was not re-run after a step-01 rerun.")
  }
}

sample_cells_final <- central_manifest[
  central_manifest$sample_id == sample_id &
  central_manifest$step_01_umi_pass == TRUE &
  central_manifest$step_02_mad_pass == TRUE &
  central_manifest$step_03_doublet_pass == TRUE,
]
n_cells_final <- nrow(sample_cells_final)
log_event("n_cells_final", n_cells_final)
if (n_cells_final == 0L) {
  fail_loud(sprintf("no cells passed all three QC steps for %s", sample_id))
}
retention_cells_pct <- 100 * n_cells_final / max(n_cells_baseline, 1)
log_event("cell_retention_pct", sprintf("%.2f", retention_cells_pct))

# Strip sample_id prefix for H5 barcode compatibility
final_passing_cells_prefixed <- sample_cells_final$cell_id
final_passing_cells <- gsub(paste0("^", sample_id, "_"), "",
                             final_passing_cells_prefixed)

# Output paths
output_base <- args$output_root
vf_out_dir     <- file.path(output_base, "features", "filtered_vfs")
rank_out_dir   <- file.path(output_base, "features", "feature_ranks")
metrics_out_dir<- file.path(output_base, "retention_metrics")
matrix_out_dir <- file.path(output_base, "normalized_matrices")
mfst_out_dir   <- file.path(output_base, "sample_manifests")
log_out_dir    <- file.path(output_base, "logs")
for (d in c(vf_out_dir, rank_out_dir, metrics_out_dir,
            matrix_out_dir, mfst_out_dir, log_out_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
vf_output_path       <- file.path(vf_out_dir, paste0(sample_id, "_vfs_filtered.txt"))
ranks_output_path    <- file.path(rank_out_dir, paste0(sample_id, "_feature_ranks.rds"))
metrics_output_path  <- file.path(metrics_out_dir, paste0(sample_id, "_retention.csv"))
matrix_output_path   <- file.path(matrix_out_dir, paste0(sample_id, "_bigsur_normalized.rds"))
sample_manifest_path <- file.path(mfst_out_dir, paste0(sample_id, ".yaml"))
log_path             <- file.path(log_out_dir, paste0(sample_id, ".log"))

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_05_post_qc_vfs\n")
  cat("Container: r_spatial\n")
  cat("Inputs validated:\n")
  cat(sprintf("  %s: OK\n", h5_path))
  cat(sprintf("  %s: OK (%d baseline VFs)\n", vf_path, length(vfs_baseline)))
  cat(sprintf("  %s: OK\n", args$central_cell_status))
  cat("Outputs planned:\n")
  cat(sprintf("  %s\n", vf_output_path))
  cat(sprintf("  %s\n", metrics_output_path))
  cat(sprintf("  %s\n", matrix_output_path))
  cat("Resources: 4 CPUs, 16 GB, 45m (SLURM array per sample)\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(BigSur)
  library(dplyr)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
SCRIPT_DIR <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(SCRIPT_DIR, "source", "post_qc_vfs.R"))

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
log_event("execute_start")

results <- run_post_qc_vfs(
  h5_path = h5_path,
  sample_id = sample_id,
  vfs_baseline = vfs_baseline,
  final_passing_cells = final_passing_cells,
  bigsur_min_cells = params$bigsur_min_cells
)

log_event("n_vfs_filtered", results$retention_metrics$n_vfs_filtered)
log_event("vf_retention_pct", sprintf("%.2f", results$retention_metrics$retention_pct))

# Writes
writeLines(results$vfs_filtered, vf_output_path)
if (ncol(results$retention_metrics) != 12L) {
  fail_loud(sprintf("retention_metrics schema: expected 12 cols, got %d",
                    ncol(results$retention_metrics)))
}
write.csv(results$retention_metrics, metrics_output_path, row.names = FALSE)
saveRDS(results$normalized_matrix, matrix_output_path, compress = TRUE)
saveRDS(results$feature_ranks, ranks_output_path, compress = TRUE)

sample_data <- list(
  sample_info = list(
    patient_id = patient_id, position_id = position_id, sample_id = sample_id,
    position_status = "step_05_complete"
  ),
  step_05_post_qc_vfs = list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    status = "pass",
    outputs = list(
      vfs_filtered = vf_output_path,
      retention_metrics = metrics_output_path,
      normalized_matrix = matrix_output_path,
      feature_ranks = ranks_output_path,
      log = log_path
    ),
    summary = list(
      n_vfs_baseline = as.integer(results$retention_metrics$n_vfs_baseline),
      n_vfs_filtered = as.integer(results$retention_metrics$n_vfs_filtered),
      retention_pct = round(results$retention_metrics$retention_pct, 2),
      n_cells_final = as.integer(n_cells_final),
      cell_retention_pct = round(retention_cells_pct, 2)
    )
  )
)
tmp_mfst <- paste0(sample_manifest_path, ".tmp.", Sys.getpid())
yaml::write_yaml(sample_data, tmp_mfst)
file.rename(tmp_mfst, sample_manifest_path)
log_event("execute_end")
quit(save = "no", status = 0)
