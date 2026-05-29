#!/usr/bin/env Rscript
# =============================================================================
# construct_rds.R - Phase 2: Build Seurat RDS from extracted components
# =============================================================================
# Takes intermediate files from convert_h5ad_to_rds.py and constructs
# a proper Seurat v5 object with PCA embedding.
#
# Inputs (from --input-dir):
#   counts.mtx.gz, barcodes.tsv.gz, features.tsv.gz, metadata.csv, pca.csv
#
# Output:
#   reed.rds (Seurat v5 object)
#
# Usage:
#   Rscript construct_rds.R --input-dir /tmp/reed_extract --output reed.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(argparse)
})

# Why do bioinformaticians prefer RDS over H5AD?
# Because they like their objects well-serialized.

main <- function(input_dir, output_file) {
  message("\n", strrep("=", 60))
  message("PHASE 2: SEURAT CONSTRUCTION")
  message(strrep("=", 60))

  # [1] Load count matrix
  message("\n[1] Loading count matrix...")
  mtx_file <- file.path(input_dir, "counts.mtx.gz")
  if (!file.exists(mtx_file)) stop("counts.mtx.gz not found in ", input_dir)

  counts <- readMM(mtx_file)
  message("  Raw matrix: ", nrow(counts), " x ", ncol(counts))

  # [2] Load barcodes and features
  message("\n[2] Loading barcodes and features...")
  barcodes <- readLines(gzfile(file.path(input_dir, "barcodes.tsv.gz")))
  barcodes <- barcodes[barcodes != ""]  # drop empty lines

  features_raw <- read.delim(gzfile(file.path(input_dir, "features.tsv.gz")),
                             header = FALSE, stringsAsFactors = FALSE)
  genes <- features_raw$V1

  message("  Barcodes: ", length(barcodes))
  message("  Features: ", length(genes))

  # Matrix is genes x cells (transposed in Python)
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  message("  Count matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

  # [3] Create Seurat object
  message("\n[3] Creating Seurat object...")
  obj <- CreateSeuratObject(counts = counts, project = "reed")
  message("  Created: ", ncol(obj), " cells x ", nrow(obj), " genes")

  # [4] Add metadata
  message("\n[4] Adding metadata...")
  meta_file <- file.path(input_dir, "metadata.csv")
  if (file.exists(meta_file)) {
    meta <- read.csv(meta_file, row.names = 1, stringsAsFactors = FALSE)
    # Only add columns not already in Seurat
    new_cols <- setdiff(names(meta), names(obj@meta.data))
    if (length(new_cols) > 0) {
      common <- intersect(rownames(meta), colnames(obj))
      obj <- AddMetaData(obj, metadata = meta[common, new_cols, drop = FALSE])
      message("  Added ", length(new_cols), " metadata columns")
    }
  }

  # [5] Add PCA embedding
  message("\n[5] Adding PCA embedding...")
  pca_file <- file.path(input_dir, "pca.csv")
  if (file.exists(pca_file)) {
    pca <- read.csv(pca_file, row.names = 1)
    pca_mat <- as.matrix(pca)
    # Align to Seurat cell order
    pca_mat <- pca_mat[colnames(obj), , drop = FALSE]
    obj[["pca"]] <- CreateDimReducObject(embeddings = pca_mat, key = "PC_",
                                          assay = DefaultAssay(obj))
    message("  PCA: ", nrow(pca_mat), " cells x ", ncol(pca_mat), " dims")
  } else {
    message("  No PCA file found — will need to compute post-hoc")
  }

  # [6] Save
  message("\n[6] Saving RDS...")
  # Write to /tmp first if output is on CRSP (stale file handle workaround)
  tmp_file <- file.path("/tmp", basename(output_file))
  saveRDS(obj, tmp_file)
  file.copy(tmp_file, output_file, overwrite = TRUE)
  file.remove(tmp_file)
  message("  Saved: ", output_file)
  message("  Size: ", round(file.info(output_file)$size / 1e9, 2), " GB")

  # [7] Validation
  message("\n[7] Validation...")
  message("  Cells: ", ncol(obj))
  message("  Genes: ", nrow(obj))
  message("  Reductions: ", paste(Reductions(obj), collapse = ", "))
  if ("pca" %in% Reductions(obj)) {
    message("  PCA dims: ", ncol(Embeddings(obj, "pca")))
  }

  message("\n", strrep("=", 60))
  message("PHASE 2 COMPLETE: ", output_file)
  message(strrep("=", 60), "\n")
}

# CLI
if (!interactive()) {
  parser <- ArgumentParser(description = "Construct Seurat RDS from H5AD intermediates")
  parser$add_argument("--input-dir", required = TRUE,
                      help = "Directory with extracted H5AD components")
  parser$add_argument("--output", required = TRUE,
                      help = "Output RDS file path")
  args <- parser$parse_args()
  main(args$input_dir, args$output)
}
