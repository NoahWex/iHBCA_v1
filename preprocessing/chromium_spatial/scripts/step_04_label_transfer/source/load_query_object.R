# Step 07 Label Transfer - Load Query Object
#
# Part of: Spatial HBCA Preprocessing Pipeline
# Step: 07_LabelTransfer
# Type: Algorithm (Pure R function, no I/O except reading inputs)
#
# Purpose: Construct Seurat query object from raw H5 file using Step 02-04
#          filtered data. Follows Step 04 wrapper pattern (lines 186-236) exactly.
#
# Dependencies:
#   - Seurat (>= 5.0) **CRITICAL: Seurat v5**
#   - dplyr
#
# Author: Agent 3 (Algorithm Developer)
# Date: 2025-11-03
# Pattern: TRACE principles (wrapper + source separation)
# Seurat Version: v5

#' Load and Process Query Seurat Object
#'
#' Constructs a Seurat query object from raw H5 file using Central Cell Status
#' Manifest to determine which cells to process. Step 04 is READ-ONLY for the
#' manifest and processes ALL cells with step_01_umi_pass == TRUE.
#'
#' @param h5_path Character. Path to raw CellRanger H5 file (from raw_data_manifest)
#' @param sample_id Character. Sample identifier for filtering manifest
#' @param central_manifest_path Character. Path to Central Cell Status Manifest CSV
#' @param step02_metadata_path Character. Path to Step 02 cell metadata CSV (for cluster_coarse)
#' @param step04_vfs_path Character. Path to Step 01 baseline VF list
#' @param min_cells Integer. Minimum cells per gene for CreateSeuratObject (default: 2)
#' @param config List. Algorithm parameters from module_configs.yaml
#'   - transfer_dims: PCA dimensions for UMAP (default: 30)
#'   - random_seed: Random seed for reproducibility (default: 42)
#'
#' @return Seurat object with:
#'   - ALL cells with step_01_umi_pass == TRUE (not filtered on Step 02/03)
#'   - Step 01 baseline VFs set as VariableFeatures
#'   - Step 02 cluster_coarse in metadata (for context)
#'   - Normalized, scaled, PCA, UMAP computed
#'
#' @details
#' Step 04 (Label Transfer) is READ-ONLY for biological annotation.
#' It processes ALL step_01_umi_pass cells regardless of Step 02/03 status.
#' This enables post-hoc QC analysis using cell type labels.
#'
#' Processing steps:
#' 1. Load Central Manifest → filter to step_01_umi_pass == TRUE for this sample
#' 2. Load Step 02 metadata → extract cluster_coarse for algorithmic context
#' 3. Load Step 01 baseline VFs
#' 4. Load raw H5, create Seurat object with Step 01 passing cells
#' 5. Set Step 01 baseline VFs as variable features
#' 6. Add Step 02 cluster_coarse to metadata (no re-clustering)
#' 7. Run Seurat pipeline: NormalizeData, ScaleData, RunPCA, RunUMAP
#'
#' @examples
#' query_obj <- load_query_object(
#'   h5_path = "/path/to/raw.h5",
#'   sample_id = "Pat1_P1",
#'   central_manifest_path = "/path/to/central_cell_status.csv",
#'   step02_metadata_path = "/path/to/step02_metadata.csv",
#'   step04_vfs_path = "/path/to/step01_baseline_vfs.txt",
#'   min_cells = 2,
#'   config = list(transfer_dims = 30, random_seed = 42)
#' )
load_query_object <- function(
  h5_path,
  sample_id,
  central_manifest_path,
  step02_metadata_path,
  step04_vfs_path,
  min_cells = 2,
  config = list()
) {

  # Extract config parameters with defaults
  transfer_dims <- config$transfer_dims %||% 30
  random_seed <- config$random_seed %||% 42

  # Set random seed for reproducibility
  set.seed(random_seed)

  cat("\n=== Loading Query Object ===\n")

  # ============================================================================
  # STEP 1: Load Central Manifest → Filter to step_01_umi_pass for This Sample
  # ============================================================================

  cat("\n[1/5] Loading cells from Central Cell Status Manifest...\n")

  # Defensive check
  if (!file.exists(central_manifest_path)) {
    stop(sprintf("Central Cell Status Manifest not found: %s", central_manifest_path))
  }

  # Load Central Manifest
  central_manifest <- read.csv(central_manifest_path, stringsAsFactors = FALSE)
  cat(sprintf("  Loaded Central Manifest: %d total cells\n", nrow(central_manifest)))

  # Convert boolean columns from string to logical
  if ("step_01_umi_pass" %in% colnames(central_manifest)) {
    central_manifest$step_01_umi_pass <- as.logical(central_manifest$step_01_umi_pass)
  }

  # Filter to THIS sample AND Step 01 passing cells
  # CRITICAL: Step 04 processes ALL step_01_umi_pass cells (not filtering on Step 02/03)
  sample_cells_manifest <- central_manifest[
    central_manifest$sample_id == sample_id & central_manifest$step_01_umi_pass == TRUE,
  ]

  cat(sprintf("  Cells for %s (step_01_umi_pass == TRUE): %d\n",
              sample_id, nrow(sample_cells_manifest)))

  if (nrow(sample_cells_manifest) == 0) {
    stop(sprintf("ERROR: No Step 01 passing cells found for %s in Central Manifest", sample_id))
  }

  # Extract prefixed barcodes for Seurat object
  final_passing_cells <- sample_cells_manifest$cell_id

  # ============================================================================
  # STEP 2: Load Step 02 Metadata → Extract cluster_coarse for Context
  # ============================================================================

  cat("\n[2/5] Loading Step 02 metadata for cluster assignments...\n")

  # Defensive check
  if (!file.exists(step02_metadata_path)) {
    stop(sprintf("Step 02 metadata not found: %s", step02_metadata_path))
  }

  # Read CSV (row names = cell barcodes, prefixed)
  metadata_step_02 <- read.csv(step02_metadata_path, stringsAsFactors = FALSE, row.names = 1)
  cat(sprintf("  Loaded Step 02 metadata: %s cells\n", nrow(metadata_step_02)))

  # Validate cluster_coarse column exists
  if (!"cluster_coarse" %in% colnames(metadata_step_02)) {
    stop("ERROR: cluster_coarse column missing from Step 02 metadata")
  }

  # Create lookup for cluster assignments (using prefixed barcodes)
  coarse_clusters <- setNames(metadata_step_02$cluster_coarse, rownames(metadata_step_02))

  cat(sprintf("  Coarse clusters: %d unique clusters\n",
              length(unique(coarse_clusters))))

  # ============================================================================
  # STEP 3: Load Step 01 Baseline Variable Features
  # ============================================================================

  cat("\n[3/5] Loading Step 01 baseline variable features...\n")

  # Defensive check
  if (!file.exists(step04_vfs_path)) {
    stop(sprintf("Step 04 filtered VFs not found: %s", step04_vfs_path))
  }

  # Read VF list (plain text, one gene per line)
  vfs_filtered <- readLines(step04_vfs_path)

  cat(sprintf("  Loaded %d baseline VFs\n", length(vfs_filtered)))

  # ============================================================================
  # STEP 4: Load Raw H5 and Create Seurat Object
  # ============================================================================

  cat("\n[4/5] Loading raw H5 and creating Seurat object...\n")

  # Defensive check
  if (!file.exists(h5_path)) {
    stop(sprintf("H5 file not found: %s", h5_path))
  }

  # Load raw counts
  raw_counts <- Seurat::Read10X_h5(h5_path)

  # CRITICAL: Cells in manifest are PREFIXED (sample_id + "_" + barcode)
  # H5 files use unprefixed barcodes, need to strip prefix before subsetting
  h5_barcodes <- gsub(paste0("^", sample_id, "_"), "", final_passing_cells)

  # Create Seurat object with unprefixed barcodes (to match H5)
  query_obj <- Seurat::CreateSeuratObject(
    counts = raw_counts,
    cells = h5_barcodes,
    min.cells = min_cells
  )

  # Re-prefix cell names to match manifest
  query_obj <- Seurat::RenameCells(
    query_obj,
    new.names = paste0(sample_id, "_", Seurat::Cells(query_obj))
  )

  cat(sprintf("  Created Seurat object: %d cells, %d genes\n",
              ncol(query_obj), nrow(query_obj)))

  # ============================================================================
  # STEP 5: Set Variable Features and Run Seurat Pipeline
  # ============================================================================

  cat("\n[5/5] Setting variable features and running Seurat pipeline...\n")

  # Set Step 01 baseline VFs
  Seurat::VariableFeatures(query_obj) <- vfs_filtered

  cat(sprintf("  Set %d variable features\n", length(Seurat::VariableFeatures(query_obj))))

  # Add coarse clusters to metadata (for algorithmic context)
  # Subset coarse_clusters to cells in query_obj (maintain order)
  query_obj$cluster_coarse <- coarse_clusters[Seurat::Cells(query_obj)]

  # Verify assignment
  cat(sprintf("  Assigned coarse clusters: %d unique clusters\n",
              length(unique(query_obj$cluster_coarse))))

  # Normalization (log-transform)
  query_obj <- Seurat::NormalizeData(
    query_obj,
    normalization.method = "LogNormalize",
    scale.factor = 1e4
  )
  cat("  Normalization complete\n")

  # Scaling (z-score, only variable features)
  query_obj <- Seurat::ScaleData(
    query_obj,
    features = vfs_filtered,
    vars.to.regress = NULL,
    do.scale = TRUE,
    do.center = TRUE,
    scale.max = 10
  )
  cat("  Scaling complete\n")

  # PCA
  query_obj <- Seurat::RunPCA(
    query_obj,
    features = vfs_filtered,
    npcs = 50,
    verbose = FALSE
  )
  cat("  PCA complete (50 PCs)\n")

  # UMAP
  query_obj <- Seurat::RunUMAP(
    query_obj,
    dims = 1:transfer_dims,
    n.neighbors = 30,
    min.dist = 0.3,
    verbose = FALSE
  )
  cat(sprintf("  UMAP complete (using %d PCs)\n", transfer_dims))

  cat("\n=== Query Object Construction Complete ===\n")
  cat(sprintf("Final object: %d cells, %d genes, %d variable features\n",
              ncol(query_obj), nrow(query_obj), length(Seurat::VariableFeatures(query_obj))))

  return(query_obj)
}
