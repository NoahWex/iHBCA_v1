# ==============================================================================
# STEP 04: Post-QC Variable Features - Core Algorithm
# ==============================================================================
# Purpose:
#   - Re-run BigSur on final passing cells (post-Step 02+03)
#   - Calculate VF retention: baseline → filtered
#   - Purely descriptive (no QC thresholds)
#
# Design:
#   - Standard BigSur function with zero-count gene filtering
#   - Single comparison: Step 01 VFs vs filtered VFs
#   - Returns data structures, wrapper handles I/O
#
# Author: Noah Wechter
# Date: 2025-10-30
# ==============================================================================

#' Standard BigSur execution with zero-count gene filtering
#'
#' @description
#' Executes BigSur variable feature discovery on a count matrix with
#' mandatory zero-count gene filtering prior to Seurat object creation.
#' This function implements the standard BigSur workflow used throughout
#' the preprocessing pipeline.
#'
#' @details
#' Processing steps:
#' 1. Filter genes with zero counts across all cells (Matrix::rowSums)
#' 2. Create Seurat object with min.cells parameter
#' 3. Run BigSur with variance-based feature discovery
#' 4. Extract and return variable feature names
#'
#' Critical design decisions:
#' - Zero-count filtering BEFORE CreateSeuratObject (memory optimization)
#' - BigSur parameters: correlations=FALSE, variable.features=TRUE, log.file=FALSE
#' - No QC thresholds applied (purely descriptive)
#'
#' @param counts_matrix Sparse matrix (dgCMatrix). Count matrix with genes as
#'   rows and cells as columns. Typically output from Seurat::Read10X_h5 or
#'   subset thereof.
#'
#' @param min_cells Integer. Minimum number of cells expressing a gene for
#'   inclusion in Seurat object (passed to CreateSeuratObject min.cells parameter).
#'   Default: 2 (from module_configs.yaml).
#'
#' @return Character vector of variable feature gene names identified by BigSur.
#'   Length varies based on biological heterogeneity in dataset.
#'
#' @examples
#' \dontrun{
#' # Load count matrix
#' counts <- Seurat::Read10X_h5("/path/to/filtered_feature_bc_matrix.h5")
#'
#' # Run BigSur with standard parameters
#' vfs <- run_bigsur_standard(counts, min_cells = 2)
#'
#' # Inspect results
#' length(vfs)
#' head(vfs)
#' }
#'
#' @export
run_bigsur_standard <- function(counts_matrix, min_cells = 2) {

  # ============================================================================
  # STEP 1: Filter zero-count genes
  # CRITICAL: Must occur BEFORE CreateSeuratObject for memory efficiency
  # ============================================================================

  message("\n--- Filtering zero-count genes ---")

  gene_counts <- Matrix::rowSums(counts_matrix)
  genes_keep <- names(gene_counts)[gene_counts > 0]
  counts_filtered <- counts_matrix[genes_keep, , drop = FALSE]

  n_genes_removed <- nrow(counts_matrix) - nrow(counts_filtered)

  message(sprintf("Removed %d zero-count genes", n_genes_removed))
  message(sprintf("Filtered matrix: %d genes × %d cells",
                  nrow(counts_filtered), ncol(counts_filtered)))

  # ============================================================================
  # STEP 2: Create Seurat object
  # ============================================================================

  message(sprintf("\n--- Creating Seurat object (min.cells = %d) ---", min_cells))

  seu <- Seurat::CreateSeuratObject(
    counts = counts_filtered,
    min.cells = min_cells
  )

  message(sprintf("Seurat object: %d cells × %d features",
                  ncol(seu), nrow(seu)))

  # Free memory
  rm(counts_filtered)
  gc(verbose = FALSE)

  # ============================================================================
  # STEP 3: Run BigSur variable feature discovery
  # ============================================================================

  message("\n--- Running BigSur ---")
  message("BigSur parameters:")
  message("  correlations: FALSE")
  message("  variable.features: TRUE")
  message("  log.file: FALSE")

  seu <- BigSur::BigSur(
    seu,
    correlations = FALSE,
    variable.features = TRUE,
    log.file = FALSE
  )

  message("BigSur complete")

  # ============================================================================
  # STEP 4: Extract variable features, normalized matrix, and feature ranks
  # ============================================================================

  vfs <- Seurat::VariableFeatures(seu)
  message(sprintf("BigSur identified %d variable features", length(vfs)))

  # Extract normalized data layer (Pearson-residual normalized expression)
  # This is used in Step 10a for debris detection via normalization comparison
  message("\n--- Extracting normalized data layer ---")
  normalized_matrix <- Seurat::GetAssayData(seu, layer = "data", assay = "RNA")
  message(sprintf("Normalized matrix: %d genes × %d cells",
                  nrow(normalized_matrix), ncol(normalized_matrix)))

  # Extract feature ranks (for Step 10a consensus feature selection)
  message("\n--- Extracting feature ranks ---")
  feature_metadata <- seu@assays$RNA@meta.data
  feature_ranks <- data.frame(
    gene = rownames(feature_metadata),
    rank = feature_metadata$var.features.rank,
    variance = feature_metadata$mcfanos,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Keep only variable features (rows with non-NA ranks)
  feature_ranks <- feature_ranks[!is.na(feature_ranks$rank), ]
  message(sprintf("Extracted ranks for %d variable features", nrow(feature_ranks)))

  # Free memory
  rm(seu)
  gc(verbose = FALSE)

  return(list(
    vfs = vfs,
    normalized_matrix = normalized_matrix,
    feature_ranks = feature_ranks
  ))
}


#' Calculate variable feature retention metrics with expression analysis
#'
#' @description
#' Compares baseline variable features (from Step 01) to filtered variable
#' features (from Step 04 post-QC) and calculates retention statistics and
#' expression metrics.
#'
#' @details
#' Retention metrics calculated:
#' - n_vfs_baseline: VF count from Step 01 (all cells passing UMI threshold)
#' - n_vfs_filtered: VF count from Step 04 (cells passing Steps 02+03)
#' - n_vfs_overlap: VFs present in both baseline and filtered sets
#' - n_vfs_lost: VFs in baseline but not in filtered (lost during QC)
#' - n_vfs_gained: VFs in filtered but not in baseline (gained during QC)
#' - retention_pct: Percentage of baseline VFs retained (overlap/baseline * 100)
#'
#' Expression metrics calculated (when counts_matrix and final_passing_cells provided):
#' - mean_expr_baseline: Mean expression of baseline VFs in baseline cells
#' - mean_expr_filtered: Mean expression of filtered VFs in filtered cells
#' - mean_expr_retained: Mean expression of retained VFs (overlap) in filtered cells
#' - mean_expr_lost: Mean expression of lost VFs in filtered cells
#' - mean_expr_gained: Mean expression of gained VFs in filtered cells
#'
#' Interpretation:
#' - High retention (>90%): QC filtering had minimal impact on VF composition
#' - Moderate retention (70-90%): Expected range for typical QC filtering
#' - Low retention (<70%): Substantial biological signal removed by QC
#' - Gained VFs: May indicate that QC removed confounding variation
#' - Expression metrics: Quantify expression levels of different VF categories
#'
#' @param vfs_baseline Character vector. Variable feature gene names from
#'   Step 01 baseline analysis (pre-QC). Typically loaded from Step 01
#'   sample_output.rds.
#'
#' @param vfs_filtered Character vector. Variable feature gene names from
#'   Step 04 post-QC analysis. Output from run_bigsur_standard().
#'
#' @param sample_id Character. Sample identifier (e.g., "Pat1_P1").
#'   Must match sample_id in manifests.
#'
#' @param counts_matrix Sparse matrix (dgCMatrix), optional. Count matrix with genes as
#'   rows and cells as columns. Required for expression metrics calculation.
#'   Should contain final passing cells only. Default: NULL.
#'
#' @param final_passing_cells Character vector, optional. Cell barcodes that passed
#'   QC (intersection of Steps 02+03). Used to verify matrix columns.
#'   Default: NULL.
#'
#' @return data.frame with 12 columns (or 7 if expression data not provided) and 1 row:
#'   \item{sample_id}{Character. Sample identifier}
#'   \item{n_vfs_baseline}{Integer. VF count from Step 01}
#'   \item{n_vfs_filtered}{Integer. VF count from Step 04}
#'   \item{n_vfs_overlap}{Integer. VFs in both sets}
#'   \item{n_vfs_lost}{Integer. VFs lost during QC}
#'   \item{n_vfs_gained}{Integer. VFs gained during QC}
#'   \item{retention_pct}{Numeric. Percentage of baseline VFs retained}
#'   \item{mean_expr_baseline}{Numeric. Mean expression of baseline VFs (if counts provided)}
#'   \item{mean_expr_filtered}{Numeric. Mean expression of filtered VFs (if counts provided)}
#'   \item{mean_expr_retained}{Numeric. Mean expression of retained VFs (if counts provided)}
#'   \item{mean_expr_lost}{Numeric. Mean expression of lost VFs (if counts provided)}
#'   \item{mean_expr_gained}{Numeric. Mean expression of gained VFs (if counts provided)}
#'
#' @examples
#' \dontrun{
#' # Load baseline VFs from Step 01
#' step01_results <- readRDS("step_01/sample_output.rds")
#' vfs_baseline <- step01_results$vf_genes
#'
#' # Run Step 04 BigSur
#' vfs_filtered <- run_bigsur_standard(counts_filtered, min_cells = 2)
#'
#' # Calculate retention with expression metrics
#' metrics <- calculate_retention_metrics(
#'   vfs_baseline = vfs_baseline,
#'   vfs_filtered = vfs_filtered,
#'   sample_id = "Pat1_P1",
#'   counts_matrix = counts_filtered,
#'   final_passing_cells = final_passing_cells
#' )
#'
#' # Inspect results
#' print(metrics)
#' }
#'
#' @export
calculate_retention_metrics <- function(vfs_baseline, vfs_filtered, sample_id,
                                       counts_matrix = NULL, final_passing_cells = NULL) {

  # ============================================================================
  # Helper function: Compute mean expression for a gene set
  # ============================================================================
  # Schema verified: 2025-10-30
  # Input: gene_set (character vector), counts (dgCMatrix)
  # Output: numeric scalar (mean of means across genes)
  # Data types: gene_set (character), counts (dgCMatrix), returns (numeric)

  compute_mean_expression <- function(gene_set, counts) {
    if (length(gene_set) == 0) {
      return(NA_real_)
    }

    # Find genes present in count matrix
    genes_present <- intersect(gene_set, rownames(counts))

    if (length(genes_present) == 0) {
      return(NA_real_)
    }

    # Calculate mean expression for each gene, then average across genes
    gene_means <- Matrix::rowMeans(counts[genes_present, , drop = FALSE])
    overall_mean <- mean(gene_means)

    return(overall_mean)
  }

  # ============================================================================
  # Calculate set overlaps
  # ============================================================================

  overlap <- intersect(vfs_baseline, vfs_filtered)
  lost <- setdiff(vfs_baseline, vfs_filtered)
  gained <- setdiff(vfs_filtered, vfs_baseline)

  # ============================================================================
  # Calculate expression metrics (if counts provided)
  # ============================================================================

  if (!is.null(counts_matrix) && !is.null(final_passing_cells)) {
    message("\n--- Calculating expression metrics ---")

    # Verify matrix contains final passing cells
    cells_in_matrix <- colnames(counts_matrix)
    if (!all(cells_in_matrix %in% final_passing_cells)) {
      warning("counts_matrix contains cells not in final_passing_cells")
    }

    # Calculate expression for each VF category
    mean_expr_baseline <- compute_mean_expression(vfs_baseline, counts_matrix)
    mean_expr_filtered <- compute_mean_expression(vfs_filtered, counts_matrix)
    mean_expr_retained <- compute_mean_expression(overlap, counts_matrix)
    mean_expr_lost <- compute_mean_expression(lost, counts_matrix)
    mean_expr_gained <- compute_mean_expression(gained, counts_matrix)

    message(sprintf("Mean expression - Baseline VFs: %.3f", mean_expr_baseline))
    message(sprintf("Mean expression - Filtered VFs: %.3f", mean_expr_filtered))
    message(sprintf("Mean expression - Retained VFs: %.3f", mean_expr_retained))
    message(sprintf("Mean expression - Lost VFs: %.3f", mean_expr_lost))
    message(sprintf("Mean expression - Gained VFs: %.3f", mean_expr_gained))

  } else {
    # Expression metrics not available
    mean_expr_baseline <- NA_real_
    mean_expr_filtered <- NA_real_
    mean_expr_retained <- NA_real_
    mean_expr_lost <- NA_real_
    mean_expr_gained <- NA_real_
  }

  # ============================================================================
  # Assemble metrics data.frame
  # Schema verified: 2025-10-30
  # Output: 12 columns (7 retention + 5 expression)
  # Data types: sample_id (character), n_* (integer), *_pct (numeric), mean_expr_* (numeric)
  # ============================================================================

  metrics <- data.frame(
    sample_id = sample_id,
    n_vfs_baseline = length(vfs_baseline),
    n_vfs_filtered = length(vfs_filtered),
    n_vfs_overlap = length(overlap),
    n_vfs_lost = length(lost),
    n_vfs_gained = length(gained),
    retention_pct = (length(overlap) / length(vfs_baseline)) * 100,
    mean_expr_baseline = mean_expr_baseline,
    mean_expr_filtered = mean_expr_filtered,
    mean_expr_retained = mean_expr_retained,
    mean_expr_lost = mean_expr_lost,
    mean_expr_gained = mean_expr_gained,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # ============================================================================
  # Log summary
  # ============================================================================

  message(sprintf("Baseline VFs: %d", metrics$n_vfs_baseline))
  message(sprintf("Filtered VFs: %d", metrics$n_vfs_filtered))
  message(sprintf("Overlap: %d (%.1f%% retention)",
                  metrics$n_vfs_overlap, metrics$retention_pct))
  message(sprintf("Lost: %d", metrics$n_vfs_lost))
  message(sprintf("Gained: %d", metrics$n_vfs_gained))

  return(metrics)
}


#' Main Step 04 post-QC variable feature discovery algorithm
#'
#' @description
#' Coordinates the complete Step 04 workflow: loads raw count data,
#' subsets to final QC-passing cells (from Steps 02+03), re-runs BigSur,
#' and calculates variable feature retention metrics by comparing to
#' Step 01 baseline results.
#'
#' @details
#' Processing workflow:
#' 1. Load raw H5 count matrix from CellRanger output
#' 2. Subset matrix to final passing cells (intersection of Step 02+03)
#' 3. Run BigSur with standard zero-count gene filtering
#' 4. Calculate retention metrics (baseline vs filtered VFs)
#' 5. Return VF list and metrics for downstream analysis
#'
#' Design principles:
#' - Purely descriptive (no QC thresholds applied)
#' - Single comparison: Step 01 baseline vs Step 04 filtered
#' - All parameters passed as function arguments (no hardcoded values)
#' - No file I/O operations (wrapper handles reads/writes)
#' - Memory management with rm() + gc() after large objects
#'
#' Integration with pipeline:
#' - Input: H5 file, Step 01 VFs, Step 02+03 passing cell barcodes
#' - Output: VF list and retention metrics for aggregation
#' - Wrapper: run_step_04_post_qc_vfs_wrapper.R handles file I/O
#'
#' @param h5_path Character. Full path to filtered_feature_bc_matrix.h5 file
#'   from CellRanger output. Must exist and be readable. Same file used in
#'   Step 01.
#'
#' @param sample_id Character. Sample identifier (e.g., "Pat1_P1").
#'   Must match sample_id in raw_data_manifest.yaml and
#'   preprocessing_manifest.yaml.
#'
#' @param vfs_baseline Character vector. Variable feature gene names from
#'   Step 01 baseline analysis. Loaded from Step 01 sample_output.rds
#'   by wrapper function.
#'
#' @param final_passing_cells Character vector. Cell barcodes that passed
#'   both Step 02 (cell-level QC) and Step 03 (doublet detection).
#'   Intersection of passing cells from both steps.
#'
#' @param bigsur_min_cells Integer. Minimum cells expressing a gene for
#'   inclusion in Seurat object (min.cells parameter). Default: 2
#'   (from module_configs.yaml).
#'
#' @return Named list with 2 elements:
#'   \item{vfs_filtered}{Character vector of variable feature gene names
#'     from post-QC BigSur analysis. Length varies by sample.}
#'   \item{retention_metrics}{data.frame with 12 columns (see
#'     calculate_retention_metrics() documentation):
#'     \itemize{
#'       \item sample_id: Character
#'       \item n_vfs_baseline: Integer
#'       \item n_vfs_filtered: Integer
#'       \item n_vfs_overlap: Integer
#'       \item n_vfs_lost: Integer
#'       \item n_vfs_gained: Integer
#'       \item retention_pct: Numeric
#'       \item mean_expr_baseline: Numeric
#'       \item mean_expr_filtered: Numeric
#'       \item mean_expr_retained: Numeric
#'       \item mean_expr_lost: Numeric
#'       \item mean_expr_gained: Numeric
#'     }
#'   }
#'
#' @examples
#' \dontrun{
#' # Load prerequisites from previous steps
#' step01_results <- readRDS("step_01/Pat1_P1/sample_output.rds")
#' step02_results <- readRDS("step_02/Pat1_P1/sample_output.rds")
#' step03_results <- readRDS("step_03/Pat1_P1/sample_output.rds")
#'
#' # Extract inputs
#' vfs_baseline <- step01_results$vf_genes
#' cells_step02 <- step02_results$metrics %>%
#'   filter(qc_pass) %>%
#'   pull(cell_barcode)
#' cells_step03 <- step03_results$metrics %>%
#'   filter(qc_pass) %>%
#'   pull(cell_barcode)
#' final_passing_cells <- intersect(cells_step02, cells_step03)
#'
#' # Run Step 04
#' results <- run_post_qc_vfs(
#'   h5_path = "/path/to/filtered_feature_bc_matrix.h5",
#'   sample_id = "Pat1_P1",
#'   vfs_baseline = vfs_baseline,
#'   final_passing_cells = final_passing_cells,
#'   bigsur_min_cells = 2
#' )
#'
#' # Access results
#' vfs_filtered <- results$vfs_filtered
#' retention_metrics <- results$retention_metrics
#'
#' # Inspect retention
#' print(retention_metrics)
#' }
#'
#' @export
run_post_qc_vfs <- function(
  h5_path,
  sample_id,
  vfs_baseline,
  final_passing_cells,
  bigsur_min_cells = 2
) {

  # ============================================================================
  # VALIDATION
  # ============================================================================

  message(sprintf("\n=== Step 04: Post-QC Variable Features for %s ===", sample_id))

  if (!file.exists(h5_path)) {
    stop(sprintf("H5 file not found: %s", h5_path))
  }

  if (length(final_passing_cells) == 0) {
    stop("final_passing_cells vector is empty - no cells passed QC")
  }

  if (length(vfs_baseline) == 0) {
    stop("vfs_baseline vector is empty - Step 01 produced no VFs")
  }

  # ============================================================================
  # STEP 1: Load raw H5 count matrix
  # ============================================================================

  message("\n--- Step 1: Loading raw H5 count matrix ---")

  counts_raw <- Seurat::Read10X_h5(h5_path)
  message(sprintf("Raw matrix: %d genes × %d cells",
                  nrow(counts_raw), ncol(counts_raw)))

  # ============================================================================
  # STEP 2: Subset to final passing cells
  # ============================================================================

  message(sprintf("\n--- Step 2: Subsetting to %d final passing cells ---",
                  length(final_passing_cells)))

  # Validate that passing cells exist in matrix
  missing_cells <- setdiff(final_passing_cells, colnames(counts_raw))
  if (length(missing_cells) > 0) {
    warning(sprintf("%d passing cells not found in H5 matrix (will be skipped)",
                    length(missing_cells)))
  }

  cells_to_keep <- intersect(final_passing_cells, colnames(counts_raw))

  if (length(cells_to_keep) == 0) {
    stop("No passing cells found in H5 matrix - check cell barcode matching")
  }

  counts_filtered <- counts_raw[, cells_to_keep, drop = FALSE]
  message(sprintf("Filtered matrix: %d genes × %d cells",
                  nrow(counts_filtered), ncol(counts_filtered)))

  # Free memory
  rm(counts_raw)
  gc(verbose = FALSE)

  # ============================================================================
  # STEP 3: Run BigSur with standard zero-count filtering
  # ============================================================================

  message("\n--- Step 3: Running BigSur on filtered cells ---")

  bigsur_results <- run_bigsur_standard(
    counts_matrix = counts_filtered,
    min_cells = bigsur_min_cells
  )

  # Extract VFs and normalized matrix from BigSur results
  vfs_filtered <- bigsur_results$vfs
  normalized_matrix <- bigsur_results$normalized_matrix

  # ============================================================================
  # STEP 4: Calculate retention metrics with expression analysis
  # CRITICAL: Must occur BEFORE freeing counts_filtered matrix
  # ============================================================================

  message("\n--- Step 4: Calculating retention metrics ---")

  retention_metrics <- calculate_retention_metrics(
    vfs_baseline = vfs_baseline,
    vfs_filtered = vfs_filtered,
    sample_id = sample_id,
    counts_matrix = counts_filtered,
    final_passing_cells = cells_to_keep
  )

  # Free memory (AFTER retention metrics calculation)
  rm(counts_filtered)
  gc(verbose = FALSE)

  # ============================================================================
  # SUMMARY
  # ============================================================================

  message("\n=== Step 04 complete ===")
  message(sprintf("Sample: %s", sample_id))
  message(sprintf("Final passing cells: %d", length(cells_to_keep)))
  message(sprintf("Variable features (filtered): %d", length(vfs_filtered)))
  message(sprintf("Retention: %.1f%% (%d / %d VFs)",
                  retention_metrics$retention_pct,
                  retention_metrics$n_vfs_overlap,
                  retention_metrics$n_vfs_baseline))

  # ============================================================================
  # Return structured results (NO FILE I/O)
  # ============================================================================

  return(list(
    vfs_filtered = vfs_filtered,
    retention_metrics = retention_metrics,
    normalized_matrix = normalized_matrix,
    feature_ranks = bigsur_results$feature_ranks
  ))
}
