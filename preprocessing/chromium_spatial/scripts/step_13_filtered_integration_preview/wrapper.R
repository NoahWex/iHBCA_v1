#!/usr/bin/env Rscript

# ============================================================================
# Step 13: Filtered Integration Preview - Wrapper
# ============================================================================
# Purpose
#   R Markdown QC + marker preview on the step 12 FOUR-WAY filtered scVI
#   integration sidecars. Structural twin of step 09, differing only in
#   the input directory (step 12 outputs instead of step 08) and the
#   module-config key used to read the expected clustering resolutions.
#
# Methodology
#   Same as step 09: validate step 12 sidecars, render
#   source/filtered_integration_preview.Rmd, emit HTML report + per-
#   resolution marker CSVs + thin Seurat object.
#
# Inputs
#   - $CFG_PREPROCESSING_STEP_12 (latent/umap/clusters/metadata sidecars)
#   - $CFG_MODULE_CONFIGS (reads clustering_resolutions from
#       step_12_filtered_scvi_integration.algorithm_params.clustering_resolutions)
#   - $CFG_RAW_DATA_SAMPLE_MANIFEST
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_13}/reports/filtered_integration_preview_{ts}.html
#   - ${CFG_PREPROCESSING_STEP_13}/markers_by_resolution/markers_res_{res}.csv
#   - ${CFG_PREPROCESSING_STEP_13}/seurat_objects/integrated_seurat.rds
#   - ${CFG_PREPROCESSING_STEP_13}/metadata/unified_metadata.csv
#
# Note: this wrapper is the step-12 analog of step 09's preview — same
# render mechanics, different inputs (step-12 outputs, step-12 config
# keys, output file prefix "filtered_integration_preview_").
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
})

parser <- ArgumentParser(description = "Step 13: Filtered Integration Preview")
parser$add_argument("--project-root", type = "character",
                    default = Sys.getenv("CFG_PROJECT_ROOT", ""))
parser$add_argument("--step-12-dir", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_12", ""))
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS", ""))
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_13", ""))
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

log_event <- function(a, v = NULL, l = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(v)) cat(sprintf("[%s] %s step_13 %s\n", ts, l, a))
  else cat(sprintf("[%s] %s step_13 %s=%s\n", ts, l, a, as.character(v)))
}

fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, l = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("dry_run", args$dry_run)

if (identical(args$step_12_dir, "") || !dir.exists(args$step_12_dir)) {
  fail_loud(sprintf("step 12 dir missing: %s", args$step_12_dir),
            "set --step-12-dir or CFG_PREPROCESSING_STEP_12")
}
required_inputs <- c(
  file.path(args$step_12_dir, "embeddings", "latent_50d.csv"),
  file.path(args$step_12_dir, "embeddings", "umap_2d.csv"),
  file.path(args$step_12_dir, "metadata", "cell_metadata_integrated.csv")
)
for (p in required_inputs) {
  if (!file.exists(p)) fail_loud(sprintf("step 12 sidecar missing: %s", p))
  log_event("input_ok", p)
}

if (identical(args$module_configs, "") || !file.exists(args$module_configs)) {
  fail_loud(sprintf("module_configs.yaml missing: %s", args$module_configs))
}
mod_cfg <- yaml::read_yaml(args$module_configs)
resolutions <- mod_cfg$modules$step_12_filtered_scvi_integration$algorithm_params$clustering_resolutions
if (is.null(resolutions) || length(resolutions) == 0) {
  fail_loud("clustering_resolutions missing from step_12_filtered_scvi_integration.algorithm_params")
}
log_event("expected_resolutions", paste(resolutions, collapse = ","))

cluster_dir <- file.path(args$step_12_dir, "clusters")
if (!dir.exists(cluster_dir)) fail_loud(sprintf("step 12 clusters dir missing: %s", cluster_dir))
missing_res <- c()
for (res in resolutions) {
  res_str <- format(res, nsmall = 1)
  cf <- file.path(cluster_dir, sprintf("leiden_res_%s.csv", res_str))
  if (!file.exists(cf)) missing_res <- c(missing_res, res)
}
if (length(missing_res) > 0) {
  fail_loud(sprintf("missing cluster files for resolutions: %s",
                    paste(missing_res, collapse = ", ")))
}

if (identical(args$raw_data_manifest, "") || !file.exists(args$raw_data_manifest)) {
  fail_loud(sprintf("raw data manifest missing: %s", args$raw_data_manifest))
}
if (identical(args$output_root, "")) {
  fail_loud("output-root empty", "set --output-root or CFG_PREPROCESSING_STEP_13")
}

reports_dir <- file.path(args$output_root, "reports")
markers_dir <- file.path(args$output_root, "markers_by_resolution")
seurat_dir  <- file.path(args$output_root, "seurat_objects")
metadata_dir <- file.path(args$output_root, "metadata")
log_dir <- file.path(args$output_root, "logs")
for (d in c(reports_dir, markers_dir, seurat_dir, metadata_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
if (file.access(args$output_root, mode = 2L) != 0L) {
  fail_loud(sprintf("output root not writable: %s", args$output_root))
}

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
notebook <- file.path(script_dir, "source", "filtered_integration_preview.Rmd")
renderer <- file.path(script_dir, "source", "render_notebook.R")
if (!file.exists(notebook)) fail_loud(sprintf("notebook source missing: %s", notebook))
if (!file.exists(renderer)) fail_loud(sprintf("render_notebook.R missing: %s", renderer))

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_html <- file.path(reports_dir,
                         sprintf("filtered_integration_preview_%s.html", timestamp))

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_13_filtered_integration_preview\n")
  cat("Container: r_spatial\n")
  cat("Filter: FOUR-WAY (inherited from step 12)\n")
  cat("Inputs validated:\n")
  for (p in required_inputs) cat(sprintf("  %s: OK\n", p))
  cat(sprintf("  %s: OK (%d resolutions)\n", cluster_dir, length(resolutions)))
  cat(sprintf("  %s: OK\n", args$module_configs))
  cat(sprintf("  %s: OK\n", args$raw_data_manifest))
  cat("Outputs planned:\n")
  cat(sprintf("  %s\n", output_html))
  cat(sprintf("  %s/markers_res_*.csv\n", markers_dir))
  cat(sprintf("  %s/integrated_seurat.rds\n", seurat_dir))
  cat(sprintf("  %s/unified_metadata.csv\n", metadata_dir))
  cat("Resources: Tier 4 (8 CPUs, 64 GB, 2h) — r_spatial container\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

log_event("render_start")
log_event("notebook", notebook)
log_event("output_html", output_html)

rmarkdown::render(
  input = notebook,
  output_file = basename(output_html),
  output_dir = reports_dir,
  knit_root_dir = args$project_root,
  params = list(
    step_12_dir = args$step_12_dir,
    raw_data_manifest = args$raw_data_manifest,
    module_configs = args$module_configs,
    output_root = args$output_root,
    resolutions = resolutions
  )
)

if (!file.exists(output_html)) {
  fail_loud(sprintf("render completed but output HTML not found: %s", output_html))
}
log_event("wrote", output_html)
log_event("bytes", file.info(output_html)$size)
log_event("done")
quit(save = "no", status = 0)
