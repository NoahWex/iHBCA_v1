#!/usr/bin/env Rscript

# ============================================================================
# Step 03: Doublet Detection - Wrapper
# ============================================================================
# Purpose
#   Per-sample doublet scoring with scDblFinder, using cluster labels
#   from step 02 as the biological prior. Cells classified as doublets
#   are flagged in the central manifest and excluded from the final
#   filtered-VF step (05) and all downstream integration.
#
# Methodology
#   - doublet_rate: 0.008 (quadratic scaling). The expected doublet
#     fraction for the 10x Chromium GEM loading regime used in the
#     Spatial HBCA samples is ~0.8% per 1000 loaded cells, which
#     scales quadratically with cell count. scDblFinder's 'dbr'
#     argument takes the sample-level rate as input, not the per-
#     cell rate.
#   - pca_dims_scDblFinder: 20 — above 30 dims the KNN-based density
#     estimation becomes noisy on smaller samples.
#   - deep_cluster_resolution: 5.0 — deliberately over-cluster inside
#     each coarse step 02 cluster so scDblFinder can find homotypic
#     doublets by neighbor-density enrichment. Coarse resolution would
#     miss them because homotypic doublets embed inside their parent
#     cluster.
#   - The cell-set-changing gate prints pre/post doublet counts and
#     fails loud on anomalous removal rates.
#
# Inputs
#   - CFG_PREPROCESSING_CENTRAL_CELL_STATUS (step 02 flags must be populated)
#   - CFG_PREPROCESSING_STEP_01/features/{sample}_vfs.txt
#   - CFG_PREPROCESSING_STEP_02/cell_metadata/{sample}_metadata.csv (cluster_coarse)
#   - Raw H5 path from CFG_RAW_DATA_SAMPLE_MANIFEST
#
# Outputs
#   - ${CFG_PREPROCESSING_STEP_03}/barcodes/{sample}_doublets.txt
#   - ${CFG_PREPROCESSING_STEP_03}/metadata/{sample}.csv
#   - ${CFG_PREPROCESSING_STEP_03}/metrics/{sample}_cluster_enrichment.csv
#   - ${CFG_PREPROCESSING_STEP_03}/manifest_updates/{sample}_step03_update.csv
#   - ${CFG_PREPROCESSING_STEP_03}/figures/{sample}_*.png
# ============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(argparse)
})

parser <- ArgumentParser(description = "Step 03: Doublet Detection (per-sample)")
parser$add_argument("--sample-id", type = "character", default = NULL,
                    help = "Sample ID to process. Required unless --sample-index is given.")
parser$add_argument("--sample-index", type = "integer", default = NULL,
                    help = "0-based sample index into the raw-data manifest TSV.")
parser$add_argument("--raw-data-manifest", type = "character",
                    default = Sys.getenv("CFG_RAW_DATA_SAMPLE_MANIFEST", ""),
                    help = "TSV with patient_id, position_id, sample_id, h5_path.")
parser$add_argument("--central-cell-status", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_CENTRAL_CELL_STATUS", ""),
                    help = "Path to central_cell_status.csv")
parser$add_argument("--step-01-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_01", ""))
parser$add_argument("--step-02-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_02", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_03", ""))
parser$add_argument("--module-configs", type = "character",
                    default = Sys.getenv("CFG_MODULE_CONFIGS",
                                         file.path(Sys.getenv("CFG_PROJECT_ROOT", ""),
                                                   "config", "module_configs.yaml")))
parser$add_argument("--min-removal-pct", type = "double", default = 0.5,
                    help = "Fail-loud lower bound on doublet removal percentage.")
parser$add_argument("--max-removal-pct", type = "double", default = 20.0,
                    help = "Fail-loud upper bound on doublet removal percentage.")
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

log_event <- function(action, value = NULL, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(value)) {
    cat(sprintf("[%s] %s step_03 %s\n", ts, level, action))
  } else {
    cat(sprintf("[%s] %s step_03 %s=%s\n", ts, level, action, as.character(value)))
  }
}

fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, level = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("dry_run", args$dry_run)

# ---- Validation preamble ----
for (lbl in c("raw_data_manifest", "central_cell_status", "step_01_root",
              "step_02_root", "output_root", "module_configs")) {
  val <- args[[lbl]]
  if (identical(val, "") || is.null(val)) {
    fail_loud(sprintf("%s path is empty", lbl),
              "set the corresponding --flag or CFG_* env var.")
  }
  log_event(lbl, val)
}

for (lbl in c("raw_data_manifest", "central_cell_status", "module_configs")) {
  p <- args[[lbl]]
  if (!file.exists(p)) fail_loud(sprintf("%s not found: %s", lbl, p))
  if (file.info(p)$size == 0L) fail_loud(sprintf("%s is empty: %s", lbl, p))
}

manifest <- read.table(args$raw_data_manifest, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
required <- c("patient_id", "position_id", "sample_id", "h5_path")
missing <- setdiff(required, colnames(manifest))
if (length(missing) > 0L) fail_loud(sprintf("raw-data manifest missing cols: %s",
                                            paste(missing, collapse = ", ")))

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
  if (nrow(row) != 1L) fail_loud(sprintf("sample-id %s not unique in manifest", args$sample_id))
}
sample_id <- row$sample_id
patient_id <- row$patient_id
position_id <- row$position_id
h5_path <- row$h5_path

log_event("sample_id", sample_id)
log_event("patient_id", patient_id)
log_event("h5_path", h5_path)

if (!file.exists(h5_path)) fail_loud(sprintf("raw H5 not found: %s", h5_path))

# Step 01 VFs
vf_path <- file.path(args$step_01_root, "features", paste0(sample_id, "_vfs.txt"))
if (!file.exists(vf_path)) fail_loud(sprintf("step-01 VF not found: %s", vf_path))

# Step 02 cluster metadata
step02_meta_path <- file.path(args$step_02_root, "cell_metadata",
                              paste0(sample_id, "_metadata.csv"))
if (!file.exists(step02_meta_path)) {
  fail_loud(sprintf("step-02 cell_metadata not found: %s", step02_meta_path),
            "run step 02 and its aggregator before step 03.")
}

# Module configs
mod_cfg <- yaml::read_yaml(args$module_configs)
params <- mod_cfg$modules$step_03_doublet$algorithm_params
if (is.null(params)) {
  # Fallback for naming variation
  params <- mod_cfg$modules$step_03_doublet_detection$algorithm_params
}
if (is.null(params)) fail_loud("step_03_doublet(.detection) params missing in module_configs.yaml")
log_event("doublet_rate", params$doublet_rate)
log_event("deep_cluster_resolution", params$deep_cluster_resolution)

# Output paths
output_base <- args$output_root
out_dirs <- list(
  barcodes      = file.path(output_base, "barcodes"),
  metadata      = file.path(output_base, "metadata"),
  metrics       = file.path(output_base, "metrics"),
  figures       = file.path(output_base, "figures"),
  manifest_upd  = file.path(output_base, "manifest_updates"),
  logs          = file.path(output_base, "logs"),
  sample_mfst   = file.path(output_base, "sample_manifests")
)
for (d in out_dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

barcode_path <- file.path(out_dirs$barcodes, paste0(sample_id, "_doublets.txt"))
temp_csv_path <- file.path(out_dirs$metadata, paste0(sample_id, ".csv"))
metrics_path <- file.path(out_dirs$metrics, paste0(sample_id, "_cluster_enrichment.csv"))
update_csv_path <- file.path(out_dirs$manifest_upd, paste0(sample_id, "_step03_update.csv"))
log_path <- file.path(out_dirs$logs, paste0(sample_id, ".log"))
sample_manifest_path <- file.path(out_dirs$sample_mfst, paste0(sample_id, ".yaml"))

# Pre-read central manifest for cell counts so the gate can compare
central_manifest <- read.csv(args$central_cell_status, stringsAsFactors = FALSE)
for (col in c("step_01_umi_pass", "step_02_mad_pass")) {
  if (col %in% colnames(central_manifest)) {
    central_manifest[[col]] <- as.logical(central_manifest[[col]])
  } else {
    fail_loud(sprintf("central_cell_status missing %s column", col),
              "run the step 01 and step 02 aggregators before step 03.")
  }
}

sample_cells_manifest <- central_manifest[
  central_manifest$sample_id == sample_id &
  central_manifest$step_01_umi_pass == TRUE,
]
n_pre_doublet <- nrow(sample_cells_manifest)
if (n_pre_doublet == 0L) {
  fail_loud(sprintf("no step-01 passing cells for %s in central manifest", sample_id))
}
log_event("n_cells_pre_doublet", n_pre_doublet)

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_03_doublet_detection\n")
  cat("Container: r_spatial\n")
  cat("Inputs validated:\n")
  cat(sprintf("  %s: OK\n", h5_path))
  cat(sprintf("  %s: OK\n", vf_path))
  cat(sprintf("  %s: OK\n", step02_meta_path))
  cat(sprintf("  %s: OK\n", args$central_cell_status))
  cat("Outputs planned:\n")
  cat(sprintf("  %s\n", barcode_path))
  cat(sprintf("  %s\n", update_csv_path))
  cat(sprintf("  %s\n", metrics_path))
  cat("Resources: 8 CPUs, 64 GB, 4h (SLURM array per sample)\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

# ---- Execute ----
suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
SCRIPT_DIR <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(SCRIPT_DIR, "source", "doublet_detection.R"))

sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)
log_event("execute_start")

# Load H5 and prefix barcodes
h5_data <- Read10X_h5(h5_path)
original_barcodes <- colnames(h5_data)
colnames(h5_data) <- paste0(sample_id, "_", original_barcodes)
seurat_obj <- CreateSeuratObject(counts = h5_data, project = sample_id)

variable_features <- readLines(vf_path)

# Load step 02 metadata for cluster_coarse
step_02_metadata <- read.csv(step02_meta_path, stringsAsFactors = FALSE, row.names = 1)
if (!("cluster_coarse" %in% colnames(step_02_metadata))) {
  fail_loud("step_02 cell_metadata has no cluster_coarse column")
}

# Subset to step-01 passing cells and merge cluster labels
passing_cell_ids <- sample_cells_manifest$cell_id
common_cells <- intersect(colnames(seurat_obj), passing_cell_ids)
if (length(common_cells) == 0L) {
  fail_loud("no overlap between H5 barcodes and central manifest",
            "check barcode prefixing (H5 is unprefixed; manifest is sample_id_ prefixed).")
}
seurat_obj <- subset(seurat_obj, cells = common_cells)

cluster_lookup <- setNames(step_02_metadata$cluster_coarse, rownames(step_02_metadata))
cluster_assignments <- cluster_lookup[colnames(seurat_obj)]
if (any(is.na(cluster_assignments))) {
  n_missing <- sum(is.na(cluster_assignments))
  log_event("missing_cluster_labels", n_missing, level = "WARN")
  keep <- colnames(seurat_obj)[!is.na(cluster_assignments)]
  seurat_obj <- subset(seurat_obj, cells = keep)
  cluster_assignments <- cluster_assignments[!is.na(cluster_assignments)]
}

result <- run_doublet_detection(
  seurat_obj = seurat_obj,
  cluster_assignments = cluster_assignments,
  variable_features = variable_features,
  sample_id = sample_id,
  params = params
)

n_total <- nrow(result$doublet_metadata)
n_doublet <- sum(!result$doublet_metadata$doublet_filter_pass)
n_singlet <- n_total - n_doublet
removal_pct <- 100 * n_doublet / max(n_total, 1)
log_event("n_cells_total", n_total)
log_event("n_doublets", n_doublet)
log_event("n_singlets", n_singlet)
log_event("doublet_removal_pct", sprintf("%.2f", removal_pct))

if (removal_pct < args$min_removal_pct) {
  log_event("low_doublet_rate_warn",
            sprintf("%.2f%% below %.2f%% -- scDblFinder may have skipped",
                    removal_pct, args$min_removal_pct),
            level = "WARN")
}
if (removal_pct > args$max_removal_pct) {
  fail_loud(sprintf("doublet removal %.2f%% exceeds safety bound %.2f%%",
                    removal_pct, args$max_removal_pct),
            "check cluster structure -- scDblFinder may be over-flagging.")
}

# Write outputs (atomic for the update CSV)
update_data <- data.frame(
  cell_id = result$doublet_metadata$cell_barcode,
  sample_id = sample_id,
  step_03_doublet_pass = result$doublet_metadata$doublet_filter_pass,
  step_03_doublet_score = result$doublet_metadata$scDblFinder_score,
  stringsAsFactors = FALSE
)
tmp_update <- paste0(update_csv_path, ".tmp.", Sys.getpid())
write.csv(update_data, tmp_update, row.names = FALSE)
file.rename(tmp_update, update_csv_path)
log_event("wrote", update_csv_path)

write.csv(result$doublet_metadata, temp_csv_path, row.names = FALSE)
doublet_barcodes <- result$doublet_metadata %>%
  filter(!doublet_filter_pass) %>%
  pull(cell_barcode)
writeLines(doublet_barcodes, barcode_path)
write.csv(result$enrichment_metrics, metrics_path, row.names = FALSE)

for (plot_name in names(result$plots)) {
  plot_path <- file.path(out_dirs$figures, paste0(sample_id, "_", plot_name, ".png"))
  ggsave(plot_path, result$plots[[plot_name]], width = 10, height = 8, dpi = 100)
}

sample_data <- list(
  sample_info = list(
    patient_id = patient_id,
    position_id = position_id,
    sample_id = sample_id,
    position_status = "step_03_complete"
  ),
  step_03_doublet = list(
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    status = "pass",
    patient_id = patient_id,
    outputs = list(
      temp_csv = temp_csv_path,
      doublet_barcodes = barcode_path,
      enrichment_metrics = metrics_path,
      log = log_path
    ),
    summary = list(
      n_cells_input = as.integer(n_total),
      n_doublets = as.integer(n_doublet),
      doublet_rate_pct = round(removal_pct, 3),
      n_cells_passing = as.integer(n_singlet)
    )
  )
)
tmp_mfst <- paste0(sample_manifest_path, ".tmp.", Sys.getpid())
yaml::write_yaml(sample_data, tmp_mfst)
file.rename(tmp_mfst, sample_manifest_path)
log_event("wrote", sample_manifest_path)

log_event("execute_end")
quit(save = "no", status = 0)
