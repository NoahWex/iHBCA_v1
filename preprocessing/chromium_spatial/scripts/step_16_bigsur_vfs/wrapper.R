#!/usr/bin/env Rscript
#
# Step 16 — BigSur Variable Features on the FIVE-WAY-filtered cell set
#
# Runs BigSur per sample on the final clean cell set (cells that pass all
# of step 01 UMI, step 02 MAD, step 03 doublet, step 10b manifold, AND
# step 15 Milo contamination filters) and produces the canonical VF lists
# that the Phase B scVI sweep (step 17) will consume. This is descriptive:
# there is no QC threshold applied — BigSur just identifies the variable
# features remaining after all upstream filters, so the resulting VF set
# is the one the atlas integration sees.
#
# Two things distinguish step 16 from step 05's baseline VF discovery:
#
#   1. FIVE-WAY filter. Step 05 used step_01..step_03 + step_10b (four-way);
#      step 16 adds step_15_pass derived from contamination_specific_cells.
#      The extra filter removes ~9.9K keratinocyte/melanocyte cells that
#      would otherwise bias the variable-feature ranks toward epidermal
#      signature genes.
#
#   2. Central manifest backfill. The first time step 16 runs, the
#      central_cell_status.csv produced by L1 aggregation does NOT yet
#      have a step_15_pass column. Step 16 merges that column in from
#      contamination_specific_cells.csv and writes the manifest back
#      atomically (temp file + file.rename on the same filesystem). This
#      is critical because central_cell_status.csv is read concurrently
#      by every sample's array task — a partial write would corrupt the
#      gate for every downstream step.
#
# Atomic manifest update: the column-merge logic uses dplyr left_join in
# memory, then writes via temp-file + rename pair so concurrent array tasks
# cannot see a half-written manifest.
#
# Output: $CFG_PREPROCESSING_STEP16_VF_DIR/{sample}_vfs_filtered.txt — this
# path is exactly what the Phase B sweep wrapper reads.

suppressPackageStartupMessages({
    library(yaml)
    library(Seurat)
    library(BigSur)
    library(dplyr)
    library(optparse)
})

# ============================================================
# CLI
# ============================================================

option_list <- list(
    make_option("--sample-index", type = "integer", default = NULL,
                help = "0-based SLURM array index (required unless --dry-run)"),
    make_option("--project-root", type = "character", default = NULL),
    make_option("--step-16-dir", type = "character", default = NULL),
    make_option("--step-01-dir", type = "character", default = NULL),
    make_option("--central-manifest", type = "character", default = NULL),
    make_option("--contamination-cells", type = "character", default = NULL),
    make_option("--raw-data-root", type = "character", default = NULL),
    make_option("--module-config", type = "character", default = NULL),
    make_option("--dry-run", action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(option_list = option_list))

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
get_env <- function(key, default = NULL) {
    val <- Sys.getenv(key, unset = NA_character_)
    if (is.na(val) || val == "") return(default)
    val
}

PROJECT_ROOT        <- opts$`project-root` %||% get_env("CFG_PROJECT_ROOT")
STEP_16_DIR         <- opts$`step-16-dir` %||% get_env("CFG_PREPROCESSING_STEP_16")
STEP_16_VF_DIR      <- get_env("CFG_PREPROCESSING_STEP16_VF_DIR",
                               if (!is.null(STEP_16_DIR))
                                   file.path(STEP_16_DIR, "features", "filtered_vfs"))
STEP_01_DIR         <- opts$`step-01-dir` %||% get_env("CFG_PREPROCESSING_STEP_01")
CENTRAL_MANIFEST    <- opts$`central-manifest` %||% get_env("CFG_PREPROCESSING_CENTRAL_CELL_STATUS")
CONTAMINATION_CELLS <- opts$`contamination-cells` %||% get_env("CFG_PREPROCESSING_CONTAMINATION_CELLS")
RAW_DATA_ROOT       <- opts$`raw-data-root` %||% get_env("CFG_RAW_DATA_ROOT")
SAMPLE_MANIFEST     <- get_env("CFG_RAW_DATA_SAMPLE_MANIFEST")
MODULE_CONFIG       <- opts$`module-config` %||%
    file.path(PROJECT_ROOT %||% ".", "config", "module_configs.yaml")

if (is.null(PROJECT_ROOT) || !nzchar(PROJECT_ROOT)) {
    stop("FATAL: CFG_PROJECT_ROOT not set. Invoke via run/run_step_16_bigsur_vfs.sh.",
         call. = FALSE)
}

cat("============================================================\n")
cat("Step 16 — BigSur VFs (FIVE-WAY filtered)\n")
cat("============================================================\n")
cat("PROJECT_ROOT:        ", PROJECT_ROOT, "\n")
cat("STEP_16_DIR:         ", STEP_16_DIR, "\n")
cat("STEP_16_VF_DIR:      ", STEP_16_VF_DIR, "\n")
cat("STEP_01_DIR:         ", STEP_01_DIR, "\n")
cat("CENTRAL_MANIFEST:    ", CENTRAL_MANIFEST, "\n")
cat("CONTAMINATION_CELLS: ", CONTAMINATION_CELLS, "\n")
cat("RAW_DATA_ROOT:       ", RAW_DATA_ROOT, "\n")
cat("SAMPLE_MANIFEST:     ", SAMPLE_MANIFEST, "\n")
cat("Dry-run:             ", opts$`dry-run`, "\n")
cat("============================================================\n\n")

# ============================================================
# Validation preamble
# ============================================================

assert_path <- function(label, path, must_be_dir = FALSE) {
    if (is.null(path) || is.na(path) || !nzchar(path)) {
        stop(paste0("FATAL: ", label, " not set"), call. = FALSE)
    }
    if (!file.exists(path)) {
        stop(paste0("FATAL: ", label, " does not exist: ", path), call. = FALSE)
    }
    if (must_be_dir && !dir.exists(path)) {
        stop(paste0("FATAL: ", label, " exists but is not a directory: ", path),
             call. = FALSE)
    }
    info <- file.info(path)
    cat(sprintf("  OK  %-22s %s (%s)\n", label, path,
                if (info$isdir) "dir" else sprintf("%.1f MB", info$size / 1024^2)))
}

cat("--- Input validation ---\n")
assert_path("PROJECT_ROOT",        PROJECT_ROOT, must_be_dir = TRUE)
assert_path("STEP_01_DIR",         STEP_01_DIR,  must_be_dir = TRUE)
assert_path("module_configs",      MODULE_CONFIG)
assert_path("CENTRAL_MANIFEST",    CENTRAL_MANIFEST)
assert_path("CONTAMINATION_CELLS", CONTAMINATION_CELLS)
assert_path("RAW_DATA_ROOT",       RAW_DATA_ROOT, must_be_dir = TRUE)

if (!is.null(SAMPLE_MANIFEST) && nzchar(SAMPLE_MANIFEST) && file.exists(SAMPLE_MANIFEST)) {
    cat(sprintf("  OK  SAMPLE_MANIFEST       %s\n", SAMPLE_MANIFEST))
}

# Ensure output dirs
for (d in c(STEP_16_DIR, STEP_16_VF_DIR,
            file.path(STEP_16_DIR, "retention_metrics"),
            file.path(STEP_16_DIR, "sample_manifests"),
            file.path(STEP_16_DIR, "logs"),
            file.path(STEP_16_DIR, "features/feature_ranks"),
            file.path(STEP_16_DIR, "normalized_matrices"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
cat(sprintf("  OK  Outputs writable      %s\n", STEP_16_DIR))

module_config <- read_yaml(MODULE_CONFIG)
params <- module_config$modules$step_16_bigsur_vfs$algorithm_params
if (is.null(params)) {
    stop(paste0("FATAL: step_16_bigsur_vfs.algorithm_params not found in ",
                MODULE_CONFIG), call. = FALSE)
}
bigsur_min_cells <- params$bigsur_min_cells %||% 2L
cat(sprintf("  OK  BigSur min.cells      %d\n", bigsur_min_cells))

# ============================================================
# Dry-run exit
# ============================================================

if (isTRUE(opts$`dry-run`)) {
    cat("\n=== DRY RUN MODE ===\n")
    cat("Step: step_16_bigsur_vfs\n")
    cat("Container: r_spatial\n")
    cat("Array: 0..61 (per sample)\n")
    cat("Outputs planned per sample:\n")
    cat("  ", file.path(STEP_16_VF_DIR, "{sample}_vfs_filtered.txt"), "\n")
    cat("  ", file.path(STEP_16_DIR, "retention_metrics/{sample}_retention.csv"), "\n")
    cat("Central manifest backfill (atomic): ", CENTRAL_MANIFEST, "\n")
    cat("Resources: 4 CPUs, 16 GB, 45 min per task\n")
    cat("VALIDATION PASSED\n")
    quit(status = 0)
}

# ============================================================
# Sample selection (array task)
# ============================================================

sample_index <- opts$`sample-index`
if (is.null(sample_index)) {
    slurm_idx <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
    if (!nzchar(slurm_idx)) {
        stop("FATAL: --sample-index not provided and SLURM_ARRAY_TASK_ID not set",
             call. = FALSE)
    }
    sample_index <- as.integer(slurm_idx)
}

# Build sample list from central manifest (samples with step 01..03 + step 10b passes)
cat("\n--- Building sample list from central manifest ---\n")
central_manifest <- read.csv(CENTRAL_MANIFEST, stringsAsFactors = FALSE)

for (bool_col in c("step_01_umi_pass", "step_02_mad_pass",
                   "step_03_doublet_pass", "step_10b_manifold_pass")) {
    if (bool_col %in% colnames(central_manifest)) {
        central_manifest[[bool_col]] <- as.logical(central_manifest[[bool_col]])
    }
}

if (!"sample_id" %in% colnames(central_manifest)) {
    stop("FATAL: central manifest missing sample_id column", call. = FALSE)
}

samples_list <- sort(unique(central_manifest$sample_id[
    central_manifest$step_01_umi_pass %in% TRUE &
    central_manifest$step_02_mad_pass %in% TRUE &
    central_manifest$step_03_doublet_pass %in% TRUE
]))
cat(sprintf("Samples passing step 01..03: %d\n", length(samples_list)))

if (sample_index < 0 || sample_index >= length(samples_list)) {
    stop(paste0("FATAL: invalid sample index ", sample_index,
                " (valid range: 0..", length(samples_list) - 1, ")"),
         call. = FALSE)
}
sample_id <- samples_list[sample_index + 1L]

cat(sprintf("\n--- Processing sample %d: %s ---\n", sample_index, sample_id))

# ============================================================
# Atomic step_15_pass backfill into central manifest
# ============================================================

if (!"step_15_pass" %in% colnames(central_manifest)) {
    cat("\nCentral manifest missing step_15_pass — backfilling from contamination_specific_cells.csv\n")

    step_15_data <- read.csv(CONTAMINATION_CELLS, stringsAsFactors = FALSE)
    if (!all(c("cell_id", "is_flagged") %in% colnames(step_15_data))) {
        stop("FATAL: contamination CSV missing required columns cell_id/is_flagged",
             call. = FALSE)
    }
    step_15_data$step_15_pass <- !as.logical(step_15_data$is_flagged)

    merged <- central_manifest %>%
        dplyr::left_join(
            step_15_data %>% dplyr::select(cell_id, step_15_pass),
            by = "cell_id"
        )
    # Cells not in step 15 output are upstream of Milo — default to pass
    merged$step_15_pass[is.na(merged$step_15_pass)] <- TRUE

    # ATOMIC WRITE: write to a temp file in the same directory then rename.
    # rename() is atomic on local and most network filesystems as long as
    # source and destination are on the same FS.
    manifest_dir <- dirname(CENTRAL_MANIFEST)
    tmp_path <- tempfile(tmpdir = manifest_dir, pattern = ".central_cell_status.",
                         fileext = ".csv.tmp")
    on.exit(if (file.exists(tmp_path)) file.remove(tmp_path), add = TRUE)

    write.csv(merged, tmp_path, row.names = FALSE)
    ok <- file.rename(tmp_path, CENTRAL_MANIFEST)
    if (!ok) {
        stop(paste0("FATAL: atomic rename failed: ", tmp_path, " -> ", CENTRAL_MANIFEST),
             call. = FALSE)
    }
    cat(sprintf("  step_15_pass merged into central manifest (%d rows)\n", nrow(merged)))
    central_manifest <- merged
} else {
    central_manifest$step_15_pass <- as.logical(central_manifest$step_15_pass)
}

# ============================================================
# Build FIVE-WAY-filtered cell set for this sample
# ============================================================

sample_cells <- central_manifest[
    central_manifest$sample_id == sample_id &
    central_manifest$step_01_umi_pass %in% TRUE &
    central_manifest$step_02_mad_pass %in% TRUE &
    central_manifest$step_03_doublet_pass %in% TRUE &
    !is.na(central_manifest$step_10b_manifold_pass) &
    central_manifest$step_10b_manifold_pass %in% TRUE &
    !is.na(central_manifest$step_15_pass) &
    central_manifest$step_15_pass %in% TRUE,
]

n_final <- nrow(sample_cells)
cat(sprintf("FIVE-WAY filtered cells for %s: %d\n", sample_id, n_final))

if (n_final == 0) {
    stop(paste0("FATAL: no cells passed FIVE-WAY filter for sample ", sample_id),
         call. = FALSE)
}

# Strip sample_id prefix from cell barcodes (H5 matrix uses unprefixed barcodes)
final_passing_cells <- sub(paste0("^", sample_id, "_"), "",
                           sample_cells$cell_id, perl = TRUE)

# Baseline step 01 VFs
baseline_vf_path <- file.path(STEP_01_DIR, "features", paste0(sample_id, "_vfs.txt"))
if (!file.exists(baseline_vf_path)) {
    stop(paste0("FATAL: step 01 baseline VF file not found: ", baseline_vf_path),
         call. = FALSE)
}
vfs_baseline <- readLines(baseline_vf_path)
cat(sprintf("Baseline VFs (step 01): %d\n", length(vfs_baseline)))

# H5 path resolution: SAMPLE_MANIFEST (TSV) first, then RAW_DATA_ROOT/{sample}.h5
h5_path <- NA_character_
if (!is.null(SAMPLE_MANIFEST) && nzchar(SAMPLE_MANIFEST) && file.exists(SAMPLE_MANIFEST)) {
    mf <- read.delim(SAMPLE_MANIFEST, stringsAsFactors = FALSE)
    if (all(c("sample_id", "h5_path") %in% colnames(mf))) {
        row <- mf[mf$sample_id == sample_id, , drop = FALSE]
        if (nrow(row) == 1L) h5_path <- row$h5_path[1]
    }
}
if (is.na(h5_path) || !nzchar(h5_path)) {
    h5_path <- file.path(RAW_DATA_ROOT, paste0(sample_id, ".h5"))
}
if (!file.exists(h5_path)) {
    stop(paste0("FATAL: raw H5 not found for ", sample_id, ": ", h5_path),
         call. = FALSE)
}
cat(sprintf("Raw H5: %s\n", h5_path))

# ============================================================
# Run BigSur
# ============================================================

SCRIPT_DIR <- local({
    cmd_args <- commandArgs(trailingOnly = FALSE)
    script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
    if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path, mustWork = FALSE))
})
source(file.path(SCRIPT_DIR, "source", "bigsur_vfs.R"))

results <- run_post_qc_vfs(
    h5_path             = h5_path,
    sample_id           = sample_id,
    vfs_baseline        = vfs_baseline,
    final_passing_cells = final_passing_cells,
    bigsur_min_cells    = bigsur_min_cells
)

# ============================================================
# Write outputs
# ============================================================

vf_output        <- file.path(STEP_16_VF_DIR, paste0(sample_id, "_vfs_filtered.txt"))
ranks_output     <- file.path(STEP_16_DIR, "features/feature_ranks", paste0(sample_id, "_feature_ranks.rds"))
metrics_output   <- file.path(STEP_16_DIR, "retention_metrics",      paste0(sample_id, "_retention.csv"))
matrix_output    <- file.path(STEP_16_DIR, "normalized_matrices",    paste0(sample_id, "_bigsur_normalized.rds"))
manifest_output  <- file.path(STEP_16_DIR, "sample_manifests",       paste0(sample_id, ".yaml"))

writeLines(results$vfs_filtered, vf_output)
cat(sprintf("VFs written:             %s (%d features)\n",
            vf_output, length(results$vfs_filtered)))

if (ncol(results$retention_metrics) != 12L) {
    stop(paste0("FATAL: retention metrics schema drift — expected 12 columns, got ",
                ncol(results$retention_metrics)), call. = FALSE)
}
write.csv(results$retention_metrics, metrics_output, row.names = FALSE, quote = FALSE)
cat(sprintf("Retention metrics:       %s\n", metrics_output))

saveRDS(results$normalized_matrix, matrix_output, compress = TRUE)
cat(sprintf("Normalized matrix:       %s (%.1f MB)\n",
            matrix_output, file.info(matrix_output)$size / 1024^2))

saveRDS(results$feature_ranks, ranks_output, compress = TRUE)
cat(sprintf("Feature ranks:           %s\n", ranks_output))

sample_yaml <- list(
    sample_info = list(sample_id = sample_id, position_status = "step_16_complete"),
    step_16_bigsur_vfs = list(
        timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        status = "pass",
        filter_type = "five_way",
        outputs = list(
            vfs_filtered       = vf_output,
            retention_metrics  = metrics_output,
            normalized_matrix  = matrix_output,
            feature_ranks      = ranks_output
        ),
        summary = list(
            n_vfs_baseline   = as.integer(results$retention_metrics$n_vfs_baseline),
            n_vfs_filtered   = as.integer(results$retention_metrics$n_vfs_filtered),
            retention_pct    = round(results$retention_metrics$retention_pct, 1),
            n_cells_final    = as.integer(n_final)
        )
    )
)
write_yaml(sample_yaml, manifest_output)
cat(sprintf("Sample manifest:         %s\n", manifest_output))

cat("\n============================================================\n")
cat(sprintf("Step 16 — %s PASSED\n", sample_id))
cat("============================================================\n")
cat(sprintf("  VF retention: %d -> %d (%.1f%%)\n",
            results$retention_metrics$n_vfs_baseline,
            results$retention_metrics$n_vfs_filtered,
            results$retention_metrics$retention_pct))
cat(sprintf("  Final cells:  %d\n", n_final))
cat("============================================================\n")
quit(status = 0)
