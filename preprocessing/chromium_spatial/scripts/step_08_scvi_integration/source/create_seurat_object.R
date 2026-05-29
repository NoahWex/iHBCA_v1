#!/usr/bin/env Rscript

# Step 08: scVI Integration - Seurat Object Creation
#
# PURPOSE: Convert scVI integration results to Seurat object
# INPUT: Intermediate H5/CSV files from Python wrapper
# OUTPUT: Seurat .rds object with scVI embeddings
#
# ARCHITECTURE:
# - Reads intermediate files prepared by Python wrapper
# - Creates Seurat object with raw counts
# - Adds scVI embeddings as DimReduc object
# - Adds UMAP coordinates as DimReduc object
# - Preserves all metadata from integration
#
# Author: Agent 4 (I/O Coordinator)
# Date: 2025-11-03

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: create_seurat_object.R <input_dir> <output_dir>")
}

input_dir <- args[1]
output_dir <- args[2]

cat(paste0("=", paste(rep("=", 80), collapse = ""), "\n"))
cat("Step 08: Creating Seurat Object from scVI Integration\n")
cat(paste0("=", paste(rep("=", 80), collapse = ""), "\n"))
cat(sprintf("Input directory:  %s\n", input_dir))
cat(sprintf("Output directory: %s\n", output_dir))
cat("\n")

# ============================================================================
# Load Counts Matrix (H5 Format with COO Components)
# ============================================================================

cat("Loading counts matrix...\n")

# Check for H5 reading capability
if (requireNamespace("rhdf5", quietly = TRUE)) {
  library(rhdf5)
  use_rhdf5 <- TRUE
  cat("  Using rhdf5 package for H5 loading\n")
} else if (requireNamespace("hdf5r", quietly = TRUE)) {
  library(hdf5r)
  use_rhdf5 <- FALSE
  cat("  Using hdf5r package for H5 loading\n")
} else {
  stop(sprintf(
    "[create_seurat_object] Neither 'rhdf5' nor 'hdf5r' package available\n  Cannot load H5 files\n  Hint: Install rhdf5 from Bioconductor: BiocManager::install('rhdf5')"
  ))
}

h5_file <- file.path(input_dir, "counts_matrix.h5")

if (!file.exists(h5_file)) {
  stop(sprintf(
    "[create_seurat_object] H5 file not found\n  Expected path: %s\n  Hint: Check that Python wrapper completed successfully",
    h5_file
  ))
}

cat("  Found H5 format (COO sparse matrix)\n")

# Load H5 file using appropriate library
if (use_rhdf5) {
  # rhdf5 implementation
  tryCatch({
    # Read COO components from H5
    coo_data <- h5read(h5_file, "counts/data")
    coo_row <- h5read(h5_file, "counts/row")
    coo_col <- h5read(h5_file, "counts/col")

    # Read dimensions from attributes
    h5_attrs <- h5readAttributes(h5_file, "/")
    n_genes <- h5_attrs$n_genes
    n_cells <- h5_attrs$n_cells
    nnz <- h5_attrs$nnz

    # Close H5 file
    H5close()

  }, error = function(e) {
    stop(sprintf(
      "[create_seurat_object] Failed to read H5 file with rhdf5\n  Path: %s\n  Error: %s\n  Hint: Check H5 file structure matches specification",
      h5_file, e$message
    ))
  })

} else {
  # hdf5r implementation
  tryCatch({
    # Open H5 file
    h5 <- H5File$new(h5_file, mode = "r")

    # Read COO components
    coo_data <- h5[["counts/data"]][]
    coo_row <- h5[["counts/row"]][]
    coo_col <- h5[["counts/col"]][]

    # Read dimensions from attributes
    n_genes <- h5$attr_open("n_genes")$read()
    n_cells <- h5$attr_open("n_cells")$read()
    nnz <- h5$attr_open("nnz")$read()

    # Close H5 file
    h5$close_all()

  }, error = function(e) {
    stop(sprintf(
      "[create_seurat_object] Failed to read H5 file with hdf5r\n  Path: %s\n  Error: %s\n  Hint: Check H5 file structure matches specification",
      h5_file, e$message
    ))
  })
}

# Validate dimensions
if (n_genes <= 0 || n_cells <= 0 || nnz <= 0) {
  stop(sprintf(
    "[create_seurat_object] Invalid dimensions in H5 file\n  Genes: %d\n  Cells: %d\n  Non-zero: %d\n  Hint: All dimensions must be positive",
    n_genes, n_cells, nnz
  ))
}

# Validate COO component lengths
if (length(coo_data) != nnz || length(coo_row) != nnz || length(coo_col) != nnz) {
  stop(sprintf(
    "[create_seurat_object] COO component length mismatch\n  Expected nnz: %d\n  data length: %d\n  row length: %d\n  col length: %d\n  Hint: All COO components must have same length as nnz",
    nnz, length(coo_data), length(coo_row), length(coo_col)
  ))
}

# CRITICAL: Convert Python 0-based indices to R 1-based indices
coo_row_r <- as.integer(coo_row) + 1L
coo_col_r <- as.integer(coo_col) + 1L

# Validate indices are in bounds
if (min(coo_row_r) < 1 || max(coo_row_r) > n_genes) {
  stop(sprintf(
    "[create_seurat_object] Row indices out of bounds after conversion\n  Min: %d (expected >= 1)\n  Max: %d (expected <= %d)\n  Hint: Index conversion error",
    min(coo_row_r), max(coo_row_r), n_genes
  ))
}
if (min(coo_col_r) < 1 || max(coo_col_r) > n_cells) {
  stop(sprintf(
    "[create_seurat_object] Column indices out of bounds after conversion\n  Min: %d (expected >= 1)\n  Max: %d (expected <= %d)\n  Hint: Index conversion error",
    min(coo_col_r), max(coo_col_r), n_cells
  ))
}

# Reconstruct sparse matrix (dgCMatrix format for Seurat compatibility)
cat("  Reconstructing sparse matrix from COO components...\n")
tryCatch({
  counts <- sparseMatrix(
    i = coo_row_r,
    j = coo_col_r,
    x = as.numeric(coo_data),
    dims = c(n_genes, n_cells)
  )
}, error = function(e) {
  stop(sprintf(
    "[create_seurat_object] Failed to reconstruct sparse matrix\n  Dimensions: %d × %d\n  Non-zero elements: %d\n  Error: %s",
    n_genes, n_cells, nnz, e$message
  ))
})

# Validate reconstructed matrix dimensions
if (nrow(counts) != n_genes || ncol(counts) != n_cells) {
  stop(sprintf(
    "[create_seurat_object] Reconstructed matrix dimension mismatch\n  Expected: %d × %d\n  Actual: %d × %d",
    n_genes, n_cells, nrow(counts), ncol(counts)
  ))
}

# Load gene names and cell IDs
gene_names_file <- file.path(input_dir, "gene_names.csv")
cell_ids_file <- file.path(input_dir, "cell_ids.csv")

if (!file.exists(gene_names_file)) {
  stop(sprintf(
    "[create_seurat_object] Gene names file not found\n  Expected path: %s",
    gene_names_file
  ))
}
if (!file.exists(cell_ids_file)) {
  stop(sprintf(
    "[create_seurat_object] Cell IDs file not found\n  Expected path: %s",
    cell_ids_file
  ))
}

gene_names <- read.csv(gene_names_file)$gene_name
cell_ids <- read.csv(cell_ids_file)$cell_id

# Validate gene and cell name lengths
if (length(gene_names) != n_genes) {
  stop(sprintf(
    "[create_seurat_object] Gene names length mismatch\n  Expected: %d\n  Actual: %d\n  Hint: CSV must contain exactly n_genes entries",
    n_genes, length(gene_names)
  ))
}
if (length(cell_ids) != n_cells) {
  stop(sprintf(
    "[create_seurat_object] Cell IDs length mismatch\n  Expected: %d\n  Actual: %d\n  Hint: CSV must contain exactly n_cells entries",
    n_cells, length(cell_ids)
  ))
}

# Check for NAs in names
if (anyNA(gene_names)) {
  stop(sprintf(
    "[create_seurat_object] NA values found in gene names\n  Number of NAs: %d",
    sum(is.na(gene_names))
  ))
}
if (anyNA(cell_ids)) {
  stop(sprintf(
    "[create_seurat_object] NA values found in cell IDs\n  Number of NAs: %d",
    sum(is.na(cell_ids))
  ))
}

# Assign dimnames
rownames(counts) <- gene_names
colnames(counts) <- cell_ids

# Report matrix statistics
cat(sprintf("  Loaded counts matrix: %d genes × %d cells (%d non-zero)\n",
            n_genes, n_cells, nnz))
cat(sprintf("  Sparsity: %.2f%%\n", 100 * (1 - nnz/(n_genes * n_cells))))
cat(sprintf("  Matrix class: %s\n", class(counts)))

# ============================================================================
# Load Embeddings and Metadata
# ============================================================================

cat("Loading embeddings and metadata...\n")

# Load scVI embeddings
scvi_emb <- read.csv(file.path(input_dir, "scvi_embeddings.csv"), row.names = 1)
cat(sprintf("  scVI embeddings: %d cells × %d dimensions\n", nrow(scvi_emb), ncol(scvi_emb)))

# Load UMAP coordinates
umap_coords <- read.csv(file.path(input_dir, "umap_coords.csv"), row.names = 1)
cat(sprintf("  UMAP coordinates: %d cells × 2 dimensions\n", nrow(umap_coords)))

# Load cell metadata
metadata <- read.csv(file.path(input_dir, "metadata.csv"), row.names = 1)
cat(sprintf("  Cell metadata: %d cells × %d columns\n", nrow(metadata), ncol(metadata)))

# Verify dimensions match
if (ncol(counts) != nrow(scvi_emb)) {
  stop(sprintf("Dimension mismatch: counts has %d cells but embeddings have %d cells",
               ncol(counts), nrow(scvi_emb)))
}

# ============================================================================
# Create Seurat Object
# ============================================================================

cat("\nCreating Seurat object...\n")

# Ensure cell order matches between counts and metadata
cell_order <- colnames(counts)
metadata <- metadata[cell_order, , drop = FALSE]
scvi_emb <- scvi_emb[cell_order, , drop = FALSE]
umap_coords <- umap_coords[cell_order, , drop = FALSE]

# Create Seurat object with raw counts
seurat_obj <- CreateSeuratObject(
  counts = counts,
  meta.data = metadata,
  project = "scVI_Integration",
  min.cells = 0,  # Don't filter - already QC'd
  min.features = 0  # Don't filter - already QC'd
)

cat(sprintf("  Created Seurat object: %d cells × %d genes\n",
            ncol(seurat_obj), nrow(seurat_obj)))

# ============================================================================
# Add scVI Embeddings as Reduction
# ============================================================================

cat("Adding scVI embeddings as reduction...\n")

# Convert to matrix and ensure correct row/column names
scvi_matrix <- as.matrix(scvi_emb)
rownames(scvi_matrix) <- colnames(seurat_obj)

# Ensure column names are properly formatted
if (!all(grepl("^scvi_", colnames(scvi_matrix)))) {
  colnames(scvi_matrix) <- paste0("scvi_", 1:ncol(scvi_matrix))
}

# Create DimReduc object for scVI
seurat_obj[["scvi"]] <- CreateDimReducObject(
  embeddings = scvi_matrix,
  key = "scvi_",
  assay = "RNA"
)

cat(sprintf("  Added 'scvi' reduction with %d dimensions\n", ncol(scvi_matrix)))

# ============================================================================
# Add UMAP as Reduction
# ============================================================================

cat("Adding UMAP coordinates as reduction...\n")

# Convert to matrix and ensure correct row/column names
umap_matrix <- as.matrix(umap_coords)
rownames(umap_matrix) <- colnames(seurat_obj)

# Ensure column names are UMAP_1 and UMAP_2
if (ncol(umap_matrix) == 2) {
  colnames(umap_matrix) <- c("UMAP_1", "UMAP_2")
}

# Create DimReduc object for UMAP
seurat_obj[["umap_scvi"]] <- CreateDimReducObject(
  embeddings = umap_matrix,
  key = "UMAPscvi_",
  assay = "RNA"
)

cat(sprintf("  Added 'umap_scvi' reduction\n"))

# ============================================================================
# Add Cluster Information to Metadata
# ============================================================================

# Check if cluster columns exist in metadata
cluster_cols <- grep("^leiden_scvi_", names(metadata), value = TRUE)
if (length(cluster_cols) > 0) {
  cat(sprintf("Found %d cluster resolution columns:\n", length(cluster_cols)))
  for (col in cluster_cols) {
    # Cluster labels should already be in metadata from Python
    cat(sprintf("  - %s: %d clusters\n", col,
                length(unique(seurat_obj@meta.data[[col]]))))
  }
} else {
  cat("No cluster columns found in metadata\n")
}

# ============================================================================
# Add Integration Summary to Misc
# ============================================================================

cat("\nAdding integration summary to misc slot...\n")

seurat_obj@misc$scvi_integration <- list(
  date = Sys.Date(),
  n_cells = ncol(seurat_obj),
  n_genes = nrow(seurat_obj),
  n_latent_dims = ncol(scvi_matrix),
  patients = unique(seurat_obj@meta.data$patient_id),
  samples = unique(seurat_obj@meta.data$sample_id),
  cluster_resolutions = cluster_cols
)

# ============================================================================
# Save Seurat Object
# ============================================================================

output_file <- file.path(output_dir, "integrated_seurat.rds")

cat(sprintf("\nSaving Seurat object to: %s\n", output_file))

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save RDS
saveRDS(seurat_obj, output_file)

# Verify file was created
if (file.exists(output_file)) {
  file_size <- file.info(output_file)$size / 1024^2  # Convert to MB
  cat(sprintf("  Saved successfully (%.1f MB)\n", file_size))
} else {
  stop("Failed to save Seurat object")
}

# ============================================================================
# Summary Report
# ============================================================================

cat("\n")
cat(paste0("=", paste(rep("=", 80), collapse = ""), "\n"))
cat("SEURAT OBJECT CREATION COMPLETE\n")
cat(paste0("=", paste(rep("=", 80), collapse = ""), "\n"))
cat(sprintf("Cells:      %d\n", ncol(seurat_obj)))
cat(sprintf("Genes:      %d\n", nrow(seurat_obj)))
cat(sprintf("Reductions: %s\n", paste(names(seurat_obj@reductions), collapse = ", ")))
cat(sprintf("Metadata:   %d columns\n", ncol(seurat_obj@meta.data)))

if (length(cluster_cols) > 0) {
  cat(sprintf("Clusters:   %d resolutions\n", length(cluster_cols)))
}

cat(sprintf("Output:     %s\n", output_file))
cat(paste0("=", paste(rep("=", 80), collapse = ""), "\n"))
