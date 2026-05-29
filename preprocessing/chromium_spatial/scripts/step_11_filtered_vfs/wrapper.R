#!/usr/bin/env Rscript

# ============================================================================
# Step 11: Post-Artifact Filtered Variable Features - Wrapper
# ============================================================================
# NOTE: An earlier sidecar_inventory.md suggested step 11 produced no
# outputs and "did not run". The I3 audit of the current 01_Preprocessing
# execution found $CFG_PREPROCESSING_STEP_11 populated with 285 files /
# ~16.7 GB (normalized_matrices, features, retention, manifests), so step
# 11 DID run in the canonical pipeline. This wrapper is ported accordingly
# and is part of the V1 dry-run matrix.
#
# Purpose
#   Per-sample recomputation of BigSur variable features and normalized
#   matrices on the FOUR-WAY filtered cell set (step_01 ∩ step_02 ∩ step_03
#   ∩ step_10b_manifold_pass). Produces post-artifact VFs consumed by
#   step 12's filtered scVI re-integration.
#
# Methodology
#   1. Resolve sample from --sample-id or --sample-index (SLURM array).
#   2. Load central_cell_status.csv; if it lacks step_10b_manifold_pass,
#      merge in the frozen retention list at
#      $CFG_PREPROCESSING_CELL_RETENTION_LIST.
#   3. Apply FOUR-WAY filter for this sample, strip sample prefix from
#      cell_ids for H5 compatibility, and call run_post_qc_vfs() from
#      source/post_qc_vfs.R (same pure algorithm used in step 05, but
#      applied to the cleaned cell population).
#   4. Write per-sample VF list, feature ranks, retention metrics (12-col
#      schema), normalized matrix, and sample manifest.
#
# Inputs
#   - $CFG_PREPROCESSING_CENTRAL_CELL_STATUS
#   - $CFG_PREPROCESSING_CELL_RETENTION_LIST (frozen step 10b; only used as
#     fallback when step_10b_manifold_pass is not in central_cell_status)
#   - $CFG_RAW_DATA_SAMPLE_MANIFEST (h5_path per sample)
#   - $CFG_MODULE_CONFIGS (step_11_filtered_vfs.algorithm_params)
#   - $CFG_PREPROCESSING_STEP_01/features/{sample}_vfs.txt (baseline VF
#     for retention comparison)
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_11}/features/filtered_vfs/{sample}_vfs_filtered.txt
#   - ${CFG_PREPROCESSING_STEP_11}/features/feature_ranks/{sample}_feature_ranks.rds
#   - ${CFG_PREPROCESSING_STEP_11}/retention_metrics/{sample}_retention.csv
#   - ${CFG_PREPROCESSING_STEP_11}/normalized_matrices/{sample}_bigsur_normalized.rds
#   - ${CFG_PREPROCESSING_STEP_11}/sample_manifests/{sample}.yaml
#   - ${CFG_PREPROCESSING_STEP_11}/logs/{sample}.log
#
# Note: sample discovery uses the raw-data manifest TSV + --sample-index
# (mirrors step 01). The central_cell_status CSV + raw-data TSV are the
# two sources of truth.
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
  library(dplyr)
})

parser <- ArgumentParser(description = "Step 11: Post-artifact filtered VFs")
parser$add_argument("--project-root", type = "character",
                    default = Sys.getenv("CFG_PROJECT_ROOT", ""))
parser$add_argument("--sample-id", type = "character", default = NULL)
parser$add_argument("--sample-index", type = "integer", default = NULL,
                    help = "0-based index into the raw-data manifest TSV (SLURM array).")
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
parser$add_argument("--central-manifest", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""))
parser$add_argument("--cell-retention-list", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_CELL_RETENTION_LIST", ""))
parser$add_argument("--step-01-dir", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_01", ""))
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_11", ""))
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

log_event <- function(a, v = NULL, l = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(v)) cat(sprintf("[%s] %s step_11 %s\n", ts, l, a))
  else cat(sprintf("[%s] %s step_11 %s=%s\n", ts, l, a, as.character(v)))
}

fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, l = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("dry_run", args$dry_run)

if (identical(args$raw_data_manifest, "") || !file.exists(args$raw_data_manifest)) {
  fail_loud(sprintf("raw data manifest missing: %s", args$raw_data_manifest))
}
manifest <- read.table(args$raw_data_manifest, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
for (col in c("patient_id", "position_id", "sample_id", "h5_path")) {
  if (!col %in% colnames(manifest)) {
    fail_loud(sprintf("raw data manifest missing column: %s", col))
  }
}

if (is.null(args$sample_id) && is.null(args$sample_index)) {
  fail_loud("must provide --sample-id or --sample-index",
            "SLURM array wrappers pass --sample-index $SLURM_ARRAY_TASK_ID")
}
if (!is.null(args$sample_index)) {
  idx <- args$sample_index + 1L
  if (idx < 1L || idx > nrow(manifest)) {
    log_event("graceful_skip",
              sprintf("sample-index %d outside manifest", args$sample_index),
              l = "WARN")
    quit(save = "no", status = 0)
  }
  row <- manifest[idx, ]
} else {
  sel <- manifest[manifest$sample_id == args$sample_id, ]
  if (nrow(sel) == 0L) fail_loud(sprintf("sample-id %s not in manifest", args$sample_id))
  row <- sel[1L, ]
}
sample_id <- row$sample_id
patient_id <- row$patient_id
position_id <- row$position_id
h5_path <- row$h5_path
log_event("sample_id", sample_id)
log_event("patient_id", patient_id)
log_event("position_id", position_id)
log_event("h5_path", h5_path)
if (!file.exists(h5_path)) fail_loud(sprintf("raw H5 not found: %s", h5_path))

if (identical(args$central_manifest, "") || !file.exists(args$central_manifest)) {
  fail_loud(sprintf("central manifest missing: %s", args$central_manifest))
}
if (identical(args$cell_retention_list, "") || !file.exists(args$cell_retention_list)) {
  fail_loud(sprintf("step 10b retention list missing: %s", args$cell_retention_list),
            "canonical frozen file; set --cell-retention-list or CFG_PREPROCESSING_CELL_RETENTION_LIST")
}

if (identical(args$step_01_dir, "") || !dir.exists(args$step_01_dir)) {
  fail_loud(sprintf("step 01 dir missing: %s", args$step_01_dir))
}
baseline_vf_path <- file.path(args$step_01_dir, "features", paste0(sample_id, "_vfs.txt"))
if (!file.exists(baseline_vf_path)) {
  fail_loud(sprintf("baseline VF file missing: %s", baseline_vf_path),
            "step 01 must run for this sample before step 11")
}

if (identical(args$module_configs, "") || !file.exists(args$module_configs)) {
  fail_loud(sprintf("module_configs.yaml missing: %s", args$module_configs))
}
mod_cfg <- yaml::read_yaml(args$module_configs)
params <- mod_cfg$modules$step_11_filtered_vfs$algorithm_params
if (is.null(params)) fail_loud("step_11_filtered_vfs.algorithm_params missing")
log_event("bigsur_min_cells", params$bigsur_min_cells)

if (identical(args$output_root, "")) {
  fail_loud("output-root empty", "set --output-root or CFG_PREPROCESSING_STEP_11")
}
for (sub in c("features/filtered_vfs", "features/feature_ranks",
              "retention_metrics", "normalized_matrices",
              "sample_manifests", "logs")) {
  dir.create(file.path(args$output_root, sub), recursive = TRUE, showWarnings = FALSE)
}
if (file.access(args$output_root, mode = 2L) != 0L) {
  fail_loud(sprintf("output root not writable: %s", args$output_root))
}

vf_output_path <- file.path(args$output_root, "features/filtered_vfs",
                            paste0(sample_id, "_vfs_filtered.txt"))
ranks_output_path <- file.path(args$output_root, "features/feature_ranks",
                               paste0(sample_id, "_feature_ranks.rds"))
metrics_output_path <- file.path(args$output_root, "retention_metrics",
                                 paste0(sample_id, "_retention.csv"))
matrix_output_path <- file.path(args$output_root, "normalized_matrices",
                                paste0(sample_id, "_bigsur_normalized.rds"))
sample_manifest_path <- file.path(args$output_root, "sample_manifests",
                                  paste0(sample_id, ".yaml"))
log_path <- file.path(args$output_root, "logs", paste0(sample_id, ".log"))

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_11_filtered_vfs\n")
  cat("Container: r_spatial\n")
  cat("Inputs validated:\n")
  cat(sprintf("  %s: OK\n", h5_path))
  cat(sprintf("  %s: OK\n", args$central_manifest))
  cat(sprintf("  %s: OK (frozen)\n", args$cell_retention_list))
  cat(sprintf("  %s: OK\n", baseline_vf_path))
  cat(sprintf("  %s: OK\n", args$module_configs))
  cat("Outputs planned:\n")
  cat(sprintf("  %s\n", vf_output_path))
  cat(sprintf("  %s\n", ranks_output_path))
  cat(sprintf("  %s\n", metrics_output_path))
  cat(sprintf("  %s\n", matrix_output_path))
  cat(sprintf("  %s\n", sample_manifest_path))
  cat("Resources: Tier 2 (4 CPUs, 16 GB, 45min, SLURM array 0-61)\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(BigSur)
  library(Matrix)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "source", "post_qc_vfs.R"))

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

log_event("execute_start")

vfs_baseline <- readLines(baseline_vf_path)
log_event("baseline_vfs", length(vfs_baseline))

central <- read.csv(args$central_manifest, stringsAsFactors = FALSE)
for (col in c("step_01_umi_pass", "step_02_mad_pass", "step_03_doublet_pass")) {
  if (!col %in% colnames(central)) fail_loud(sprintf("central manifest missing: %s", col))
  central[[col]] <- as.logical(central[[col]])
}
if ("step_10b_manifold_pass" %in% colnames(central)) {
  central$step_10b_manifold_pass <- as.logical(central$step_10b_manifold_pass)
  log_event("step_10b_source", "central_manifest")
} else {
  retention <- read.csv(args$cell_retention_list, stringsAsFactors = FALSE)
  retention$step_10b_manifold_pass <- as.logical(retention$step_10b_manifold_pass)
  merge_cols <- intersect(
    c("cell_id", "step_10b_manifold_pass", "step_10b_removal_reason"),
    colnames(retention)
  )
  central <- central %>% left_join(retention[, merge_cols], by = "cell_id")
  log_event("step_10b_source", "retention_list_merged")
}

sample_cells <- central[
  central$sample_id == sample_id &
  central$step_01_umi_pass & central$step_02_mad_pass &
  central$step_03_doublet_pass &
  !is.na(central$step_10b_manifold_pass) & central$step_10b_manifold_pass,
]
n_final <- nrow(sample_cells)
log_event("cells_four_way_pass", n_final)
if (n_final == 0L) {
  fail_loud(sprintf("no cells pass FOUR-WAY filter for sample %s", sample_id))
}

# Central manifest uses prefixed cell_ids; H5 uses unprefixed barcodes
final_passing_cells <- sub(paste0("^", sample_id, "_"), "", sample_cells$cell_id)

results <- run_post_qc_vfs(
  h5_path = h5_path,
  sample_id = sample_id,
  vfs_baseline = vfs_baseline,
  final_passing_cells = final_passing_cells,
  bigsur_min_cells = params$bigsur_min_cells
)

if (ncol(results$retention_metrics) != 12L) {
  fail_loud(sprintf("retention_metrics schema: expected 12 cols, got %d",
                    ncol(results$retention_metrics)))
}

writeLines(results$vfs_filtered, vf_output_path)
log_event("wrote_vfs", vf_output_path)
log_event("n_vfs_filtered", length(results$vfs_filtered))

write.csv(results$retention_metrics, metrics_output_path,
          row.names = FALSE, quote = FALSE)
log_event("wrote_retention", metrics_output_path)

saveRDS(results$normalized_matrix, matrix_output_path, compress = TRUE)
log_event("wrote_matrix_mb",
          round(file.size(matrix_output_path) / 1024^2, 2))

saveRDS(results$feature_ranks, ranks_output_path, compress = TRUE)
log_event("wrote_feature_ranks", ranks_output_path)

sample_data <- list(
  sample_info = list(
    patient_id = patient_id,
    position_id = position_id,
    sample_id = sample_id,
    position_status = "step_11_complete"
  ),
  step_11_filtered_vfs = list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    status = "pass",
    filter_type = "four_way",
    filter_description = "step_01 AND step_02 AND step_03 AND step_10b_manifold_pass",
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
      retention_pct = round(results$retention_metrics$retention_pct, 1),
      n_cells_four_way_pass = as.integer(n_final)
    )
  )
)
tmp <- paste0(sample_manifest_path, ".tmp.", Sys.getpid())
yaml::write_yaml(sample_data, tmp)
file.rename(tmp, sample_manifest_path)
log_event("wrote", sample_manifest_path)
log_event("done")
quit(save = "no", status = 0)
