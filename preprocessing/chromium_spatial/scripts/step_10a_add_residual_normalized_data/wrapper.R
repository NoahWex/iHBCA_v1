#!/usr/bin/env Rscript

# ============================================================================
# Step 10a: Add BigSur Pearson-Residual Normalization - Wrapper
# ============================================================================
# Purpose
#   Attach BigSur Pearson-residual normalized expression as a new "BigSur"
#   assay on the step 09 integrated Seurat object, and compute discordance
#   metrics between the raw/smoothed counts and the BigSur-normalized
#   residuals. The BigSur assay provides a variance-stabilized matrix used
#   downstream at step 10b for manifold artifact discovery.
#
# Methodology
#   1. Load the step 09 integrated Seurat object (reassembled from scVI
#      sidecars by step 09).
#   2. Per sample, load the step 05 per-sample BigSur normalized matrix
#      (step 05 normalized_matrices/{sample}_bigsur_normalized.rds), use
#      rank-aggregation to pick ~2000 consensus features across samples,
#      and vertically slice each matrix to those features.
#   3. Create per-sample Seurat objects from the sliced matrices, merge +
#      JoinLayers (Seurat v5 pattern), extract the resulting assay and
#      attach it as "BigSur" on the integrated object.
#   4. Run compute_discordance_metrics() across 5 metrics × 4 conditions
#      (20 CSVs) measuring per-cell/per-feature discordance between the
#      raw RNA assay and the BigSur residuals.
#
# Inputs
#   - ${CFG_PREPROCESSING_STEP_09}/seurat_objects/integrated_seurat.rds
#   - ${CFG_PREPROCESSING_STEP_05}/normalized_matrices/{sample}_bigsur_normalized.rds
#   - ${CFG_PREPROCESSING_STEP_05}/features/feature_ranks/{sample}_feature_ranks.rds
#   - $CFG_MODULE_CONFIGS (step_10a_add_residual_normalized_data.algorithm_params)
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_10A}/seurat_objects/integrated_with_bigsur.rds
#   - ${CFG_PREPROCESSING_STEP_10A}/discordance/metrics/{metric}/{condition}.csv  (20 CSVs)
#   - ${CFG_PREPROCESSING_STEP_10A}/logs/step_10a.log
#   - ${CFG_PREPROCESSING_STEP_10A}/completion_manifest.yaml
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
})

parser <- ArgumentParser(description = "Step 10a: Add BigSur normalization")
parser$add_argument("--project-root", type = "character",
                    default = Sys.getenv("CFG_PROJECT_ROOT", ""))
parser$add_argument("--step-09-dir", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_09", ""))
parser$add_argument("--step-05-dir", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_05", ""))
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_10A", ""))
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

log_event <- function(a, v = NULL, l = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(v)) cat(sprintf("[%s] %s step_10a %s\n", ts, l, a))
  else cat(sprintf("[%s] %s step_10a %s=%s\n", ts, l, a, as.character(v)))
}

fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, l = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("dry_run", args$dry_run)

if (identical(args$step_09_dir, "") || !dir.exists(args$step_09_dir)) {
  fail_loud(sprintf("step 09 dir missing: %s", args$step_09_dir),
            "set --step-09-dir or CFG_PREPROCESSING_STEP_09")
}
integrated_seurat_path <- file.path(args$step_09_dir, "seurat_objects", "integrated_seurat.rds")
if (!file.exists(integrated_seurat_path)) {
  fail_loud(sprintf("step 09 integrated Seurat missing: %s", integrated_seurat_path),
            "step 09 must complete before step 10a")
}
log_event("integrated_seurat", integrated_seurat_path)

if (identical(args$step_05_dir, "") || !dir.exists(args$step_05_dir)) {
  fail_loud(sprintf("step 05 dir missing: %s", args$step_05_dir),
            "set --step-05-dir or CFG_PREPROCESSING_STEP_05")
}
bigsur_dir <- file.path(args$step_05_dir, "normalized_matrices")
ranks_dir <- file.path(args$step_05_dir, "features", "feature_ranks")
if (!dir.exists(bigsur_dir)) {
  fail_loud(sprintf("step 05 normalized_matrices/ missing: %s", bigsur_dir))
}
if (!dir.exists(ranks_dir)) {
  fail_loud(sprintf("step 05 features/feature_ranks/ missing: %s", ranks_dir))
}
log_event("bigsur_dir", bigsur_dir)
log_event("ranks_dir", ranks_dir)

if (identical(args$raw_data_manifest, "") || !file.exists(args$raw_data_manifest)) {
  fail_loud(sprintf("raw data manifest missing: %s", args$raw_data_manifest))
}
manifest <- read.table(args$raw_data_manifest, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE)
if (!all(c("sample_id") %in% colnames(manifest))) {
  fail_loud("raw_data_manifest must have sample_id column")
}
sample_ids <- manifest$sample_id
log_event("sample_count", length(sample_ids))

# Discover BigSur matrices + feature_ranks per sample
bigsur_matrix_paths <- list()
feature_ranks_paths <- list()
for (sid in sample_ids) {
  bp <- file.path(bigsur_dir, paste0(sid, "_bigsur_normalized.rds"))
  rp <- file.path(ranks_dir, paste0(sid, "_feature_ranks.rds"))
  if (file.exists(bp)) bigsur_matrix_paths[[sid]] <- bp
  else log_event("bigsur_matrix_missing", sid, l = "WARN")
  if (file.exists(rp)) feature_ranks_paths[[sid]] <- rp
  else log_event("feature_ranks_missing", sid, l = "WARN")
}
if (length(bigsur_matrix_paths) == 0) {
  fail_loud("no step 05 BigSur matrices found",
            "step 05 must write normalized_matrices/{sample}_bigsur_normalized.rds")
}
if (length(feature_ranks_paths) == 0) {
  fail_loud("no step 05 feature_ranks files found",
            "step 05 must write features/feature_ranks/{sample}_feature_ranks.rds")
}
log_event("bigsur_matrices_found", length(bigsur_matrix_paths))
log_event("feature_ranks_found", length(feature_ranks_paths))

if (identical(args$module_configs, "") || !file.exists(args$module_configs)) {
  fail_loud(sprintf("module_configs.yaml missing: %s", args$module_configs))
}
mod_cfg <- yaml::read_yaml(args$module_configs)
algo_params <- mod_cfg$modules$step_10a_add_residual_normalized_data$algorithm_params
if (is.null(algo_params)) {
  fail_loud("step_10a_add_residual_normalized_data.algorithm_params missing from module_configs.yaml")
}
n_consensus_features <- algo_params$n_consensus_features
if (is.null(n_consensus_features)) n_consensus_features <- 2000
log_event("n_consensus_features", n_consensus_features)

if (identical(args$output_root, "")) {
  fail_loud("output-root empty", "set --output-root or CFG_PREPROCESSING_STEP_10A")
}
for (sub in c("seurat_objects", "logs", "discordance/metrics/slope",
              "discordance/metrics/pearson_r", "discordance/metrics/pearson_r2",
              "discordance/metrics/spearman_r", "discordance/metrics/spearman_r2")) {
  dir.create(file.path(args$output_root, sub), recursive = TRUE, showWarnings = FALSE)
}
if (file.access(args$output_root, mode = 2L) != 0L) {
  fail_loud(sprintf("output root not writable: %s", args$output_root))
}

output_path <- file.path(args$output_root, "seurat_objects", "integrated_with_bigsur.rds")
log_path <- file.path(args$output_root, "logs", "step_10a.log")

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_10a_add_residual_normalized_data\n")
  cat("Container: r_spatial\n")
  cat("Inputs validated:\n")
  cat(sprintf("  %s: OK\n", integrated_seurat_path))
  cat(sprintf("  %s: OK (%d matrices)\n", bigsur_dir, length(bigsur_matrix_paths)))
  cat(sprintf("  %s: OK (%d rank files)\n", ranks_dir, length(feature_ranks_paths)))
  cat(sprintf("  %s: OK\n", args$module_configs))
  cat("Outputs planned:\n")
  cat(sprintf("  %s\n", output_path))
  cat(sprintf("  %s/discordance/metrics/*.csv (20 files)\n", args$output_root))
  cat("Resources: Tier 3 (8 CPUs, 64 GB, 2h) — r_spatial container\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))

source(file.path(script_dir, "source", "add_bigsur_normalization.R"))
source(file.path(script_dir, "source", "compute_discordance.R"))

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

log_event("execute_start")
result <- add_bigsur_to_integrated_object(
  integrated_seurat_path = integrated_seurat_path,
  bigsur_matrix_paths = bigsur_matrix_paths,
  feature_ranks_paths = feature_ranks_paths,
  output_path = output_path,
  n_consensus_features = n_consensus_features,
  verbose = TRUE
)
summary_stats <- result$summary_stats
log_event("cells_integrated", summary_stats$n_cells)
log_event("bigsur_genes", summary_stats$n_genes_bigsur)
log_event("object_mb", round(summary_stats$object_size_mb, 1))

log_event("discordance_start")
seu <- result$updated_seurat
discordance_result <- compute_discordance_metrics(
  seurat_object = seu,
  output_base = args$output_root,
  rna_assay = "RNA",
  bigsur_assay = "BigSur",
  scvi_reduction = "scvi",
  n_scvi_dims = 50,
  k_neighbors = 20,
  verbose = TRUE
)
log_event("discordance_cells", discordance_result$summary_stats$n_cells)

completion <- list(
  step_info = list(
    step = "step_10a_add_residual_normalized_data",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    status = "pass"
  ),
  outputs = list(
    integrated_with_bigsur = output_path,
    log = log_path,
    discordance_root = file.path(args$output_root, "discordance", "metrics")
  ),
  summary = summary_stats
)
yaml::write_yaml(completion, file.path(args$output_root, "completion_manifest.yaml"))
log_event("wrote_completion_manifest")
log_event("done")
quit(save = "no", status = 0)
