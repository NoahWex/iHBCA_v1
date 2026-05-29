#!/usr/bin/env Rscript

# ============================================================================
# Step 01: Baseline Variable Features - Wrapper
# ============================================================================
# Purpose
#   Per-sample quality control and variable feature discovery on raw
#   CellRanger output. This is the foundation of the L1 preprocessing
#   layer: it defines the starting cell set (flat UMI filter) and the
#   per-sample variable feature list that downstream steps use for
#   clustering context and integration input.
#
# Methodology
#   1. Flat UMI threshold (>= 150 UMI/cell). Cells below this threshold
#      are dropped before Seurat object creation to keep peak memory
#      within the tier-2 budget for a 4-CPU / 24G node. The threshold is
#      deliberate: below ~150 UMI a cell cannot support meaningful
#      clustering or scDblFinder scoring downstream, so there is no
#      benefit to carrying it into step 02/03.
#   2. BigSur variance-based variable feature discovery (NOT
#      correlation-based). Variance-based VFs are more stable with
#      low-coverage cells than Pearson-residual ranks, which matters
#      because the spatial dataset spans a wide coverage range across
#      its 62 samples. Sample passes QC iff n_vfs >= 100.
#   3. Per-cell metrics written for ALL cells (passing + failing). The
#      failing cells remain in the central manifest so later steps can
#      audit removal rates by patient/position.
#
# Inputs
#   - $CFG_RAW_DATA_SAMPLE_MANIFEST : TSV with columns
#       patient_id, position_id, sample_id, h5_path
#     (built once during migration from raw_data_manifest.yaml)
#   - The H5 file referenced by the sample's row.
#   - Algorithm parameters from $CFG_MODULE_CONFIGS
#     (step_01_baseline_vfs.algorithm_params).
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_01}/features/{sample_id}_vfs.txt
#       BigSur variance VF list, one gene per line.
#   - ${CFG_PREPROCESSING_STEP_01}/cell_metadata/{sample_id}_metadata.csv
#       Per-cell QC metrics (cell_id, sample_id, patient_id, position,
#       flat_threshold_pass, n_umi, n_genes, sample_bigsur_vf_count,
#       sample_qc_pass). Contributes to central_cell_status.csv via the
#       step-01 aggregator.
#   - ${CFG_PREPROCESSING_STEP_01}/sample_manifests/{sample_id}.yaml
#       Per-sample summary (n_vfs, status, cell counts). Aggregated into
#       the step-01 manifest.
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

parser <- ArgumentParser(description = "Step 01: Baseline Variable Features")
parser$add_argument("--sample-id", type = "character", default = NULL,
                    help = "Sample ID to process (e.g. Pat1_P1). Required unless --sample-index is given.")
parser$add_argument("--sample-index", type = "integer", default = NULL,
                    help = "0-based sample index into the raw-data manifest TSV. Used by SLURM array.")
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""),
                    help = "TSV with columns patient_id, position_id, sample_id, h5_path.")
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS",
                                         file.path(Sys.getenv("CFG_PROJECT_ROOT", ""),
                                                   "config", "module_configs.yaml")),
                    help = "Path to module_configs.yaml (holds algorithm_params).")
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_01", ""),
                    help = "Step 01 output root (defaults to CFG_PREPROCESSING_STEP_01).")
parser$add_argument("--qc-status-dir", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_QC_STATUS", ""),
                    help = "Directory for qc_status/temp manifest shards (defaults to CFG_PREPROCESSING_QC_STATUS).")
parser$add_argument("--dry-run", action = "store_true", default = FALSE,
                    help = "Validate inputs and print the plan without running BigSur.")

args <- parser$parse_args()

# ---------------------------------------------------------------------------
# Structured logger
# ---------------------------------------------------------------------------

log_event <- function(action, value = NULL, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(value)) {
    cat(sprintf("[%s] %s step_01 %s\n", ts, level, action))
  } else {
    cat(sprintf("[%s] %s step_01 %s=%s\n", ts, level, action, as.character(value)))
  }
}

fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, level = "ERROR")
  if (!is.null(hint)) {
    cat(sprintf("  hint: %s\n", hint))
  }
  quit(save = "no", status = 2)
}

# ---------------------------------------------------------------------------
# Validation preamble
# ---------------------------------------------------------------------------

log_event("start")
log_event("dry_run", args$dry_run)

# Container / environment sanity
for (env_var in c("CFG_PROJECT_ROOT", "CFG_PREPROCESSING_STEP_01")) {
  if (identical(Sys.getenv(env_var), "")) {
    log_event("env_missing", env_var, level = "WARN")
  }
}

# Raw data manifest
if (identical(args$raw_data_manifest, "") || is.null(args$raw_data_manifest)) {
  fail_loud("raw-data manifest path is empty",
            "set --raw-data-manifest or CFG_RAW_DATA_SAMPLE_MANIFEST")
}
if (!file.exists(args$raw_data_manifest)) {
  fail_loud(sprintf("raw-data manifest not found: %s", args$raw_data_manifest),
            "build the TSV at $CFG_RAW_DATA_SAMPLE_MANIFEST before submitting step 01")
}
log_event("raw_data_manifest", args$raw_data_manifest)

manifest <- tryCatch(
  read.table(args$raw_data_manifest, header = TRUE, sep = "\t",
             stringsAsFactors = FALSE, check.names = FALSE),
  error = function(e) fail_loud(sprintf("failed to parse raw_data_manifest TSV: %s", e$message))
)
required_cols <- c("sample_id", "h5_path")
missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols) > 0) {
  fail_loud(sprintf("raw_data_manifest missing required columns: %s",
                    paste(missing_cols, collapse = ", ")),
            "TSV must declare sample_id and h5_path. patient_id and position_id are derived from sample_id ({patient_id}_{position_id}) when not present.")
}
# Derive patient_id and position_id from sample_id when not in the manifest.
# Convention: sample_id = {patient_id}_{position_id} where patient_id is the
# first underscore-delimited token (e.g. Pat1_P2_Lower -> Pat1, P2_Lower).
if (!"patient_id" %in% colnames(manifest)) {
  manifest$patient_id <- sub("_.*$", "", manifest$sample_id)
}
if (!"position_id" %in% colnames(manifest)) {
  manifest$position_id <- sub("^[^_]+_", "", manifest$sample_id)
}
log_event("manifest_rows", nrow(manifest))

# Resolve the row to process
if (is.null(args$sample_id) && is.null(args$sample_index)) {
  fail_loud("must provide --sample-id or --sample-index",
            "SLURM array wrappers should pass --sample-index $SLURM_ARRAY_TASK_ID")
}
if (!is.null(args$sample_index)) {
  idx <- args$sample_index + 1L  # 0-based -> 1-based
  if (idx < 1L || idx > nrow(manifest)) {
    # Graceful skip: SLURM arrays are sized to allow trailing slots that may
    # not correspond to a manifest row (e.g. a sample registered without
    # chromium data). An out-of-range index is treated as an intentionally
    # empty array slot, not a failure.
    log_event("graceful_skip",
              sprintf("sample-index %d outside manifest (no row to process)",
                      args$sample_index),
              level = "WARN")
    quit(save = "no", status = 0)
  }
  row <- manifest[idx, ]
  sample_id <- row$sample_id
  patient_id <- row$patient_id
  position_id <- row$position_id
  h5_path <- row$h5_path
} else {
  sel <- manifest[manifest$sample_id == args$sample_id, ]
  if (nrow(sel) == 0L) {
    fail_loud(sprintf("sample-id %s not found in manifest", args$sample_id))
  }
  if (nrow(sel) > 1L) {
    fail_loud(sprintf("sample-id %s is ambiguous in manifest (%d rows)",
                      args$sample_id, nrow(sel)))
  }
  sample_id <- sel$sample_id
  patient_id <- sel$patient_id
  position_id <- sel$position_id
  h5_path <- sel$h5_path
}

log_event("sample_id", sample_id)
log_event("patient_id", patient_id)
log_event("position_id", position_id)
log_event("h5_path", h5_path)

if (!file.exists(h5_path)) {
  fail_loud(sprintf("raw H5 not found: %s", h5_path),
            "check the h5_path column in the raw-data manifest")
}
h5_size <- file.info(h5_path)$size
log_event("h5_size_bytes", h5_size)
if (h5_size <= 0L) {
  fail_loud(sprintf("raw H5 is empty (0 bytes): %s", h5_path))
}

# Module configs
if (!file.exists(args$module_configs)) {
  fail_loud(sprintf("module_configs.yaml not found: %s", args$module_configs),
            "port config/module_configs.yaml to the canonical config dir")
}
mod_cfg <- yaml::read_yaml(args$module_configs)
params <- mod_cfg$modules$step_01_baseline_vfs$algorithm_params
if (is.null(params)) {
  fail_loud("step_01_baseline_vfs.algorithm_params missing in module_configs.yaml")
}
log_event("umi_threshold", params$umi_threshold)
log_event("bigsur_vf_threshold", params$bigsur_vf_threshold)

# Output root
if (identical(args$output_root, "") || is.null(args$output_root)) {
  fail_loud("output root is empty",
            "set --output-root or CFG_PREPROCESSING_STEP_01")
}
output_base <- args$output_root
vf_dir <- file.path(output_base, "features")
metadata_dir <- file.path(output_base, "cell_metadata")
manifest_dir <- file.path(output_base, "sample_manifests")
log_dir <- file.path(output_base, "logs")

vf_path <- file.path(vf_dir, paste0(sample_id, "_vfs.txt"))
metadata_csv <- file.path(metadata_dir, paste0(sample_id, "_metadata.csv"))
sample_manifest_path <- file.path(manifest_dir, paste0(sample_id, ".yaml"))
log_path <- file.path(log_dir, paste0(sample_id, ".log"))

log_event("output_root", output_base)
log_event("vf_path", vf_path)
log_event("metadata_csv", metadata_csv)

for (d in c(vf_dir, metadata_dir, manifest_dir, log_dir)) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}
if (file.access(output_base, mode = 2L) != 0L) {
  fail_loud(sprintf("output root not writable: %s", output_base))
}

# ---------------------------------------------------------------------------
# Dry-run exit
# ---------------------------------------------------------------------------

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat(sprintf("Step: step_01_baseline_vfs\n"))
  cat(sprintf("Container: r_spatial\n"))
  cat("Inputs validated:\n")
  cat(sprintf("  %s: OK (%d bytes)\n", h5_path, h5_size))
  cat(sprintf("  %s: OK\n", args$module_configs))
  cat(sprintf("  %s: OK\n", args$raw_data_manifest))
  cat("Outputs planned:\n")
  cat(sprintf("  %s\n", vf_path))
  cat(sprintf("  %s\n", metadata_csv))
  cat(sprintf("  %s\n", sample_manifest_path))
  cat("Resources: 4 CPUs, 24 GB, 2h (SLURM array per sample)\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------------------
# Algorithm execution
# ---------------------------------------------------------------------------

# Load heavy libraries only after validation passes (keeps dry-run fast
# and container-agnostic).
suppressPackageStartupMessages({
  library(Seurat)
  library(BigSur)
  library(Matrix)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
SCRIPT_DIR <- dirname(normalizePath(script_path, mustWork = FALSE))
source_path <- file.path(SCRIPT_DIR, "source", "baseline_vf.R")
if (!file.exists(source_path)) {
  fail_loud(sprintf("algorithm source not found: %s", source_path),
            "ensure scripts/step_01_baseline_vfs/source/baseline_vf.R is present")
}
source(source_path)

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

log_event("execute_start")
log_event("cells_pre_umi_filter_note",
          "see per-cell metadata for flat_threshold_pass counts")

results <- run_baseline_vf(
  h5_path = h5_path,
  sample_id = sample_id,
  patient_id = patient_id,
  position_id = position_id,
  umi_threshold = params$umi_threshold,
  bigsur_vf_threshold = params$bigsur_vf_threshold,
  bigsur_min_cells = params$bigsur_min_cells,
  bigsur_correlations = params$bigsur_correlations,
  bigsur_variable_features = params$bigsur_variable_features,
  bigsur_log_file = params$bigsur_log_file
)

log_event("cells_total", nrow(results$metrics))
log_event("cells_passing_umi", sum(results$metrics$flat_threshold_pass))
log_event("n_vfs", results$n_vfs)
log_event("sample_qc_pass", results$sample_qc_pass)

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

write.csv(results$metrics, metadata_csv, row.names = FALSE)
log_event("wrote", metadata_csv)
log_event("bytes", file.info(metadata_csv)$size)

writeLines(results$vf_genes, vf_path)
log_event("wrote", vf_path)
log_event("bytes", file.info(vf_path)$size)

# Write qc_status temp shard (6 columns) that aggregate_step_01 will
# concat into central_cell_status.csv. Atomic write via temp + rename.
qc_status_dir <- args$qc_status_dir
if (identical(qc_status_dir, "") || is.null(qc_status_dir)) {
  qc_status_dir <- file.path(output_base, "..", "qc_status")
}
qc_temp_dir <- file.path(qc_status_dir, "temp")
dir.create(qc_temp_dir, recursive = TRUE, showWarnings = FALSE)
qc_temp_csv <- file.path(qc_temp_dir, paste0(sample_id, "_step01.csv"))
manifest_data <- data.frame(
  cell_id = results$metrics$cell_id,
  sample_id = results$metrics$sample_id,
  patient_id = results$metrics$patient_id,
  position = results$metrics$position,
  step_01_umi_pass = results$metrics$flat_threshold_pass,
  step_01_n_umi = results$metrics$n_umi,
  stringsAsFactors = FALSE,
  row.names = NULL
)
tmp_qc <- paste0(qc_temp_csv, ".tmp.", Sys.getpid())
write.csv(manifest_data, tmp_qc, row.names = FALSE)
file.rename(tmp_qc, qc_temp_csv)
log_event("wrote", qc_temp_csv)

# Sample manifest (atomic write: temp + rename)
status <- if (results$sample_qc_pass) "pass" else "fail"
position_status <- if (status == "pass") "step_01_complete" else "step_01_failed"
sample_data <- list(
  sample_info = list(
    patient_id = patient_id,
    position_id = position_id,
    sample_id = sample_id,
    position_status = position_status
  ),
  step_01_baseline_vfs = list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    status = status,
    outputs = list(
      vfs_list = vf_path,
      cell_metadata = metadata_csv,
      log = log_path
    ),
    summary = list(
      n_vfs = as.integer(results$n_vfs),
      qc_pass = results$sample_qc_pass,
      n_cells_total = nrow(results$metrics),
      n_cells_passing = sum(results$metrics$flat_threshold_pass)
    )
  )
)
tmp_manifest <- paste0(sample_manifest_path, ".tmp.", Sys.getpid())
yaml::write_yaml(sample_data, tmp_manifest)
file.rename(tmp_manifest, sample_manifest_path)
log_event("wrote", sample_manifest_path)

log_event("execute_end")
if (results$sample_qc_pass) {
  quit(save = "no", status = 0)
} else {
  # Sample below VF threshold -- wrapper succeeded, sample excluded from
  # downstream by its QC status. Exit 0 so the SLURM array task is not
  # marked failed; the sample manifest carries the status.
  quit(save = "no", status = 0)
}
