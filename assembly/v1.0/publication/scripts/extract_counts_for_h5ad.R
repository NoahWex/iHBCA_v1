#!/usr/bin/env Rscript
# =============================================================================
# EXTRACT COUNTS FOR H5AD — R Phase
# =============================================================================
# Extracts raw count matrix from Seurat RDS and writes 10x-format intermediates
# for downstream Python h5ad assembly.
#
# Outputs:
#   counts.mtx.gz    - sparse count matrix (genes x cells, Market Matrix)
#   features.tsv.gz  - gene names (one per line)
#   barcodes.tsv.gz  - cell IDs (one per line)
#   metadata.csv     - cell metadata (only if --extract-metadata)
#
# Usage:
#   Rscript extract_counts_for_h5ad.R \
#     --study gray \
#     --input-rds /path/to/seurat.rds \
#     --output-dir /path/to/intermediates/
#
# Plan: B1_source_datasets_external
# =============================================================================

# Seurat objects are like onions — layers, layers, and they make you cry.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# -----------------------------------------------------------------------------
# Argument parsing (manual — argparse not guaranteed in container)
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(args, flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx + 1 > length(args)) stop(paste("Missing value for", flag))
  args[idx + 1]
}

has_flag <- function(args, flag) {
  flag %in% args
}

study <- parse_arg(args, "--study")
input_rds <- parse_arg(args, "--input-rds")
output_dir <- parse_arg(args, "--output-dir")
extract_meta <- has_flag(args, "--extract-metadata")

if (is.null(study) || is.null(input_rds) || is.null(output_dir)) {
  stop("Required args: --study <name> --input-rds <path> --output-dir <path>")
}

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------
message("=" , strrep("=", 69))
message("EXTRACT COUNTS FOR H5AD")
message("Study: ", study)
message("Input: ", input_rds)
message("Output: ", output_dir)
message("Extract metadata: ", extract_meta)
message("=", strrep("=", 69))

if (!file.exists(input_rds)) {
  stop("FATAL: Input RDS not found: ", input_rds)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Use SLURM $TMPDIR (per-job, more space) or fall back to /tmp
# CRSP stale file handle workaround: always write to local, then copy
tmp_base <- Sys.getenv("TMPDIR", "/tmp")
tmp_dir <- file.path(tmp_base, paste0("extract_", study, "_", Sys.getpid()))
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

# Log file
log_file <- file.path(output_dir, "extraction_log.txt")
log_con <- file(log_file, open = "wt")
on.exit(close(log_con), add = TRUE)

log_msg <- function(...) {
  msg <- paste0(...)
  message(msg)
  writeLines(msg, log_con)
  flush(log_con)
}

# -----------------------------------------------------------------------------
# Load RDS
# -----------------------------------------------------------------------------
log_msg("[1/5] Loading RDS...")
t0 <- Sys.time()
srt <- readRDS(input_rds)
t1 <- Sys.time()
log_msg("  Loaded in ", round(difftime(t1, t0, units = "mins"), 2), " minutes")

# v3 -> v5 compat
tryCatch({
  srt <- UpdateSeuratObject(srt)
  log_msg("  UpdateSeuratObject applied")
}, error = function(e) {
  log_msg("  UpdateSeuratObject not needed or failed: ", conditionMessage(e))
})

n_cells <- ncol(srt)
default_assay <- DefaultAssay(srt)
available_assays <- Assays(srt)

log_msg("  Cells: ", format(n_cells, big.mark = ","))
log_msg("  Default assay: ", default_assay)
log_msg("  Available assays: ", paste(available_assays, collapse = ", "))

# Find the raw-count assay: RNA > originalexp > first non-derived assay
raw_assay <- NULL
for (candidate in c("RNA", "originalexp")) {
  if (candidate %in% available_assays) {
    raw_assay <- candidate
    break
  }
}
if (is.null(raw_assay)) {
  non_derived <- setdiff(available_assays, c("integrated", "SCT", "prediction.score"))
  if (length(non_derived) > 0) raw_assay <- non_derived[1]
  else raw_assay <- available_assays[1]
}

if (raw_assay != default_assay) {
  log_msg("  Switching default assay from '", default_assay, "' to '", raw_assay, "'")
  DefaultAssay(srt) <- raw_assay
}

# For Assay5 with split layers (e.g. counts.1, counts.2 from merge), join first
assay_obj <- srt[[raw_assay]]
if (inherits(assay_obj, "Assay5")) {
  avail_layers <- Layers(assay_obj)
  log_msg("  Assay5 layers: ", paste(avail_layers, collapse = ", "))
  if (any(grepl("^counts\\.", avail_layers))) {
    log_msg("  Joining split count layers...")
    srt[[raw_assay]] <- JoinLayers(srt[[raw_assay]])
  }
}

n_features <- nrow(srt)
log_msg("  Features (", raw_assay, "): ", format(n_features, big.mark = ","))

# -----------------------------------------------------------------------------
# Extract raw counts
# -----------------------------------------------------------------------------
log_msg("\n[2/5] Extracting raw count matrix...")

# Strategy: try layer="counts" (Assay5), check for empty result, then
# fall back to slot="counts" (v3 Assay). layer= on v3 Assay returns 0x0
# WITHOUT error, so we must check dimensions explicitly.
counts <- NULL

# Approach 1: layer="counts" (Assay5)
tryCatch({
  m <- GetAssayData(srt, assay = raw_assay, layer = "counts")
  if (prod(dim(m)) > 0) {
    counts <- m
    log_msg("  Retrieved via layer='counts'")
  }
}, error = function(e) {
  log_msg("  layer='counts' not available: ", conditionMessage(e))
})

# Approach 2: slot="counts" (v3 Assay)
if (is.null(counts)) {
  tryCatch({
    m <- GetAssayData(srt, assay = raw_assay, slot = "counts")
    if (prod(dim(m)) > 0) {
      counts <- m
      log_msg("  Retrieved via slot='counts'")
    }
  }, error = function(e) {
    log_msg("  slot='counts' not available: ", conditionMessage(e))
  })
}

# Approach 3: data slot (log-normalized — last resort)
if (is.null(counts)) {
  tryCatch({
    m <- GetAssayData(srt, assay = raw_assay, slot = "data")
    if (prod(dim(m)) > 0) {
      counts <- m
      log_msg("  WARNING: Using 'data' slot — values may be log-normalized, not raw counts")
    }
  }, error = function(e) {
    log_msg("  slot='data' not available: ", conditionMessage(e))
  })
}

if (is.null(counts)) {
  stop("FATAL: Could not extract count matrix from assay '", raw_assay,
       "'. Tried layer='counts', slot='counts', slot='data'.")
}

if (!inherits(counts, "dgCMatrix") && !inherits(counts, "dgTMatrix")) {
  log_msg("  Converting to sparse matrix (class was: ", class(counts)[1], ")")
  counts <- as(counts, "dgCMatrix")
}

log_msg("  Count matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")
log_msg("  Non-zero entries: ", format(nnzero(counts), big.mark = ","))
log_msg("  Sparsity: ", round(1 - nnzero(counts) / (as.numeric(nrow(counts)) * ncol(counts)), 4) * 100, "%")

# Verify genes x cells orientation (R/Seurat convention)
stopifnot(nrow(counts) == n_features)
stopifnot(ncol(counts) == n_cells)

# -----------------------------------------------------------------------------
# Write counts to /tmp, then copy (CRSP workaround)
# -----------------------------------------------------------------------------
log_msg("\n[3/5] Writing count matrix...")

tmp_mtx <- file.path(tmp_dir, "counts.mtx")
writeMM(counts, tmp_mtx)
log_msg("  Wrote: counts.mtx (", file.size(tmp_mtx), " bytes)")

# Gzip
system2("gzip", c("-f", tmp_mtx))
tmp_mtx_gz <- paste0(tmp_mtx, ".gz")
log_msg("  Compressed: counts.mtx.gz (", file.size(tmp_mtx_gz), " bytes)")

# Copy to output
final_mtx_gz <- file.path(output_dir, "counts.mtx.gz")
file.copy(tmp_mtx_gz, final_mtx_gz, overwrite = TRUE)
log_msg("  Copied to: ", final_mtx_gz)

# -----------------------------------------------------------------------------
# Write features and barcodes
# -----------------------------------------------------------------------------
log_msg("\n[4/5] Writing features and barcodes...")

# Features (gene names, one per line)
gene_names <- rownames(counts)
tmp_feat <- file.path(tmp_dir, "features.tsv")
writeLines(gene_names, tmp_feat)
system2("gzip", c("-f", tmp_feat))
file.copy(paste0(tmp_feat, ".gz"), file.path(output_dir, "features.tsv.gz"), overwrite = TRUE)
log_msg("  Features: ", length(gene_names), " genes written")

# Barcodes (cell IDs, one per line)
cell_ids <- colnames(counts)
tmp_bc <- file.path(tmp_dir, "barcodes.tsv")
writeLines(cell_ids, tmp_bc)
system2("gzip", c("-f", tmp_bc))
file.copy(paste0(tmp_bc, ".gz"), file.path(output_dir, "barcodes.tsv.gz"), overwrite = TRUE)
log_msg("  Barcodes: ", length(cell_ids), " cells written")

# -----------------------------------------------------------------------------
# Extract metadata (Kumar only)
# -----------------------------------------------------------------------------
if (extract_meta) {
  log_msg("\n[5/5] Extracting metadata from Seurat object...")
  meta <- srt@meta.data
  meta$cell_id <- rownames(meta)

  tmp_meta <- file.path(tmp_dir, "metadata.csv")
  write.csv(meta, tmp_meta, row.names = FALSE)
  file.copy(tmp_meta, file.path(output_dir, "metadata.csv"), overwrite = TRUE)
  log_msg("  Metadata: ", nrow(meta), " cells x ", ncol(meta), " columns")
  log_msg("  Columns: ", paste(colnames(meta), collapse = ", "))
} else {
  log_msg("\n[5/5] Metadata extraction skipped (not requested)")
}

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
log_msg("\n", strrep("=", 70))
log_msg("EXTRACTION COMPLETE")
log_msg("  Study: ", study)
log_msg("  Cells: ", format(n_cells, big.mark = ","))
log_msg("  Genes: ", format(n_features, big.mark = ","))
log_msg("  Assay: ", raw_assay)
log_msg("  Output: ", output_dir)
log_msg("  Files: counts.mtx.gz, features.tsv.gz, barcodes.tsv.gz",
        if (extract_meta) ", metadata.csv" else "")
log_msg(strrep("=", 70))
