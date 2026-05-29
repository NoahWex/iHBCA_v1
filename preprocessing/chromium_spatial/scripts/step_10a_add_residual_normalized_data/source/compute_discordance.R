# ==============================================================================
# Compute Discordance Metrics for Debris Detection - VECTORIZED VERSION
# ==============================================================================
# Purpose: Calculate per-cell discordance between LogNorm abundance and BigSur
#          significance using SNN-smoothed signal processing with comprehensive
#          condition testing and feature increment analysis
#
# Algorithm:
#   1. Compute SNN graph from scVI embeddings
#   2. Extract BigSur feature set (2,000 consensus features)
#   3. Extract raw matrices + ScaleData on consensus features
#   4. SNN-weighted smoothing (raw + scaled versions)
#   5. Rank-sum expression scoring (per-cell feature ranking)
#   6. Calculate 4 conditions × 5 increments using VECTORIZED masked moments
#      - raw_unscaled, raw_scaled, smooth_unscaled, smooth_scaled
#      - N ∈ {100, 200, 400, 1000, 2000}
#   7. Write 20 minimal CSVs (5 metrics × 4 conditions)
#   8. Return summary stats
#
# Optimization: MASKED MOMENT EXPANSION (Gemini feedback)
#   - Creates binary mask M[g,c] = 1 if gene g is in top N for cell c
#   - Computes Cov(X,Y) = E[XY] - E[X]E[Y] using element-wise multiplication
#   - Replaces per-cell loops with vectorized matrix operations
#   - ~100x faster than cell loop
#
# Author: Noah Wechter
# Date: 2025-11-19
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(matrixStats)  # For fast colRanks
})

# ==============================================================================
# DEBUG MODE TOGGLE
# ==============================================================================
DEBUG_MODE <- FALSE
DEBUG_N_CELLS <- 1000

# ==============================================================================
# HELPER FUNCTION: Compute Metrics for One Condition (VECTORIZED)
# ==============================================================================
#' @description
#' Calculates 5 metrics (slope, Pearson r/R², Spearman r/R²) for all cells
#' across 5 feature increments using rank-sum adaptive feature selection.
#'
#' VECTORIZED using Masked Moment Expansion:
#' - Creates binary mask M[g,c] = 1 if gene g is in top N for cell c
#' - Computes Cov(X,Y) = E[XY] - E[X]E[Y] using element-wise multiplication
#' - Processes cells in batches to manage RAM
#'
#' @param mat_X Matrix (features × cells) - LogNorm or variant
#' @param mat_Y Matrix (features × cells) - BigSur or variant
#' @param rank_sum Matrix (features × cells) - Combined expression ranking
#' @param precomputed_X_ranked Matrix (features × cells) - Precomputed ranks for Spearman
#' @param precomputed_Y_ranked Matrix (features × cells) - Precomputed ranks for Spearman
#' @param increments Vector of N values (e.g., c(100, 200, 400, 1000, 2000))
#' @param batch_size Cells per batch (default: 5000 for memory management)
#' @param verbose Print progress
#'
#' @return List with matrices (n_increments × n_cells) for each metric
compute_metrics_for_condition <- function(mat_X, mat_Y, rank_sum,
                                          precomputed_X_ranked, precomputed_Y_ranked,
                                          increments, batch_size = 5000, verbose = TRUE) {

  n_cells <- ncol(mat_X)
  n_increments <- length(increments)

  # Initialize result matrices
  slope_mat <- matrix(NA, nrow = n_increments, ncol = n_cells)
  pearson_r_mat <- matrix(NA, nrow = n_increments, ncol = n_cells)
  pearson_r2_mat <- matrix(NA, nrow = n_increments, ncol = n_cells)
  spearman_r_mat <- matrix(NA, nrow = n_increments, ncol = n_cells)
  spearman_r2_mat <- matrix(NA, nrow = n_increments, ncol = n_cells)

  # Use precomputed ranks for Spearman (passed from main function)
  mat_X_ranked <- precomputed_X_ranked
  mat_Y_ranked <- precomputed_Y_ranked

  # Process cells in batches to manage RAM
  starts <- seq(1, n_cells, by = batch_size)
  n_batches <- length(starts)

  if (verbose) {
    cat(sprintf("  Processing %d cells in %d batches (batch size: %d)\n",
                n_cells, n_batches, batch_size))
  }

  # Loop over increments
  for (idx in 1:n_increments) {
    N <- increments[idx]

    if (verbose) {
      cat(sprintf("  Increment N=%d: ", N))
    }

    # Initialize vectors for this increment
    slopes_N <- numeric(n_cells)
    pearson_r_N <- numeric(n_cells)
    spearman_r_N <- numeric(n_cells)

    # Process batches
    for (batch_idx in 1:n_batches) {
      start_idx <- starts[batch_idx]
      end_idx <- min(start_idx + batch_size - 1, n_cells)
      cell_indices <- start_idx:end_idx

      # ========================================================================
      # 1. Subset Data for Batch
      # ========================================================================
      X_batch <- mat_X[, cell_indices, drop = FALSE]
      Y_batch <- mat_Y[, cell_indices, drop = FALSE]
      R_batch <- rank_sum[, cell_indices, drop = FALSE]

      X_ranked_batch <- mat_X_ranked[, cell_indices, drop = FALSE]
      Y_ranked_batch <- mat_Y_ranked[, cell_indices, drop = FALSE]

      # ========================================================================
      # 2. Create Binary Mask (1 = Top N, 0 = Ignore)
      # ========================================================================
      # Use colRanks with descending order (negate rank_sum)
      # Genes with rank <= N are in the mask
      ranks_batch <- colRanks(-R_batch, ties.method = "first", preserveShape = TRUE)
      Mask <- (ranks_batch <= N) * 1  # Convert to numeric 0/1

      # ========================================================================
      # 3. Compute Moments (Vectorized) - Pearson
      # ========================================================================
      sum_X <- colSums(X_batch * Mask)
      sum_Y <- colSums(Y_batch * Mask)
      sum_XY <- colSums(X_batch * Y_batch * Mask)
      sum_XX <- colSums(X_batch^2 * Mask)
      sum_YY <- colSums(Y_batch^2 * Mask)

      # Means and second moments
      mean_X <- sum_X / N
      mean_Y <- sum_Y / N
      mean_XY <- sum_XY / N
      mean_XX <- sum_XX / N
      mean_YY <- sum_YY / N

      # Variance and covariance (population formulas, using N)
      var_X <- mean_XX - mean_X^2
      var_Y <- mean_YY - mean_Y^2
      cov_XY <- mean_XY - (mean_X * mean_Y)

      # ========================================================================
      # 4. Calculate Pearson Metrics
      # ========================================================================
      # Slope = Cov(X,Y) / Var(X)
      batch_slopes <- cov_XY / var_X
      batch_slopes[!is.finite(batch_slopes)] <- 0

      # Pearson r = Cov(X,Y) / sqrt(Var(X) * Var(Y))
      batch_pearson <- cov_XY / sqrt(var_X * var_Y)
      batch_pearson[!is.finite(batch_pearson)] <- 0

      # ========================================================================
      # 5. Compute Spearman (same approach on ranked data)
      # ========================================================================
      sum_XR <- colSums(X_ranked_batch * Mask)
      sum_YR <- colSums(Y_ranked_batch * Mask)
      sum_XRYR <- colSums(X_ranked_batch * Y_ranked_batch * Mask)
      sum_XRXR <- colSums(X_ranked_batch^2 * Mask)
      sum_YRYR <- colSums(Y_ranked_batch^2 * Mask)

      mean_XR <- sum_XR / N
      mean_YR <- sum_YR / N
      mean_XRYR <- sum_XRYR / N
      mean_XRXR <- sum_XRXR / N
      mean_YRYR <- sum_YRYR / N

      var_XR <- mean_XRXR - mean_XR^2
      var_YR <- mean_YRYR - mean_YR^2
      cov_XRYR <- mean_XRYR - (mean_XR * mean_YR)

      batch_spearman <- cov_XRYR / sqrt(var_XR * var_YR)
      batch_spearman[!is.finite(batch_spearman)] <- 0

      # ========================================================================
      # 6. Store Results for This Batch
      # ========================================================================
      slopes_N[cell_indices] <- batch_slopes
      pearson_r_N[cell_indices] <- batch_pearson
      spearman_r_N[cell_indices] <- batch_spearman

      # Progress reporting
      if (verbose && batch_idx %% max(1, n_batches %/% 4) == 0) {
        cat(sprintf("%d%%...", round(100 * batch_idx / n_batches)))
      }
    }

    if (verbose) {
      cat(" Complete\n")
    }

    # Store results for this increment
    slope_mat[idx, ] <- slopes_N
    pearson_r_mat[idx, ] <- pearson_r_N
    pearson_r2_mat[idx, ] <- pearson_r_N^2
    spearman_r_mat[idx, ] <- spearman_r_N
    spearman_r2_mat[idx, ] <- spearman_r_N^2
  }

  return(list(
    slope = slope_mat,
    pearson_r = pearson_r_mat,
    pearson_r2 = pearson_r2_mat,
    spearman_r = spearman_r_mat,
    spearman_r2 = spearman_r2_mat
  ))
}

# ==============================================================================
# MAIN FUNCTION: Compute Discordance Metrics
# ==============================================================================
compute_discordance_metrics <- function(
  seurat_object,
  output_base,
  rna_assay = "RNA",
  bigsur_assay = "BigSur",
  scvi_reduction = "scvi",
  n_scvi_dims = 50,
  k_neighbors = 20,
  verbose = TRUE
) {

  # ============================================================================
  # DEBUG MODE: Subset to small number of cells for testing
  # ============================================================================

  if (DEBUG_MODE) {
    if (verbose) {
      cat("\n")
      cat("*** DEBUG MODE ACTIVE ***\n")
      cat(sprintf("*** Subsetting to %d cells for rapid testing ***\n", DEBUG_N_CELLS))
      cat("\n")
    }

    n_cells_total <- ncol(seurat_object)
    subset_cells <- sample(colnames(seurat_object), min(DEBUG_N_CELLS, n_cells_total))
    seurat_object <- seurat_object[, subset_cells]

    if (verbose) {
      cat(sprintf("Subsetted to %d cells (from %d)\n", ncol(seurat_object), n_cells_total))
    }
  }

  # ============================================================================
  # Step 1: Compute SNN Graph from scVI Embeddings
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 1: Computing SNN graph from scVI embeddings\n")
    cat("============================================================\n")
  }

  # Verify scVI reduction exists
  if (!scvi_reduction %in% names(seurat_object@reductions)) {
    stop(sprintf("scVI reduction '%s' not found in Seurat object!", scvi_reduction))
  }

  # Compute SNN graph
  seurat_object <- FindNeighbors(
    seurat_object,
    reduction = scvi_reduction,
    dims = 1:n_scvi_dims,
    k.param = k_neighbors,
    graph.name = "scvi_snn",
    verbose = FALSE
  )

  # Validate graph creation
  if (!"scvi_snn" %in% names(seurat_object@graphs)) {
    stop("SNN graph creation failed!")
  }

  adj <- seurat_object@graphs$scvi_snn

  if (verbose) {
    cat("SNN graph created:\n")
    cat(sprintf("  Dimensions: %d × %d\n", nrow(adj), ncol(adj)))
    cat(sprintf("  Non-zero edges: %d\n", length(adj@x)))
    cat(sprintf("  Mean degree: %.2f\n", mean(Matrix::rowSums(adj > 0))))
    cat(sprintf("  Sparsity: %.2f%%\n", 100 * (1 - length(adj@x) / (nrow(adj) * ncol(adj)))))
  }

  # ============================================================================
  # Step 2: Extract BigSur Feature Set
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 2: Extracting BigSur feature set\n")
    cat("============================================================\n")
  }

  # Extract BigSur assay data to get consensus feature names
  mat_bigsur <- GetAssayData(seurat_object, assay = bigsur_assay, layer = "counts")
  features <- rownames(mat_bigsur)

  if (verbose) {
    cat(sprintf("  BigSur features: %d genes\n", length(features)))
  }

  # Row-normalize adjacency matrix (convert to transition probabilities)
  row_sums <- Matrix::rowSums(adj)
  adj_norm <- adj / row_sums

  # ============================================================================
  # Step 3: Extract Raw Matrices + ScaleData on Consensus Features
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 3: Extracting and scaling raw matrices\n")
    cat("============================================================\n")
  }

  # Extract LogNorm (data layer), subset to BigSur features
  mat_log <- GetAssayData(seurat_object, assay = rna_assay, layer = "data")
  mat_log <- mat_log[features, ]

  if (verbose) {
    cat(sprintf("  LogNorm matrix: %d genes × %d cells\n", nrow(mat_log), ncol(mat_log)))
    cat(sprintf("  BigSur matrix: %d genes × %d cells\n", nrow(mat_bigsur), ncol(mat_bigsur)))
  }

  # ScaleData on consensus features for both assays
  if (verbose) {
    cat("  Running ScaleData on consensus features...\n")
  }

  # Scale LogNorm (RNA assay)
  DefaultAssay(seurat_object) <- rna_assay
  seurat_object <- ScaleData(seurat_object, features = features, verbose = FALSE)
  mat_log_scaled <- GetAssayData(seurat_object, assay = rna_assay, layer = "scale.data")
  mat_log_scaled <- mat_log_scaled[features, ]

  # Scale BigSur assay (copy counts→data first for ScaleData)
  DefaultAssay(seurat_object) <- bigsur_assay

  # BigSur only has counts layer - copy to data for ScaleData
  counts_data <- GetAssayData(seurat_object, assay = bigsur_assay, layer = "counts")
  seurat_object <- SetAssayData(seurat_object, assay = bigsur_assay, layer = "data", new.data = counts_data)

  seurat_object <- ScaleData(seurat_object, features = features, verbose = FALSE)
  mat_bigsur_scaled <- GetAssayData(seurat_object, assay = bigsur_assay, layer = "scale.data")

  if (verbose) {
    cat(sprintf("  Scaled LogNorm: %d genes × %d cells\n", nrow(mat_log_scaled), ncol(mat_log_scaled)))
    cat(sprintf("  Scaled BigSur: %d genes × %d cells\n", nrow(mat_bigsur_scaled), ncol(mat_bigsur_scaled)))
  }

  # ============================================================================
  # Step 4: SNN-Weighted Smoothing (Raw + Scaled)
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 4: SNN-weighted smoothing (raw + scaled)\n")
    cat("============================================================\n")
  }

  # Smooth raw matrices
  mat_bigsur_smooth <- adj_norm %*% t(mat_bigsur)
  mat_bigsur_smooth <- t(mat_bigsur_smooth)

  mat_log_smooth <- adj_norm %*% t(mat_log)
  mat_log_smooth <- t(mat_log_smooth)

  # Smooth scaled matrices
  mat_bigsur_scaled_smooth <- adj_norm %*% t(mat_bigsur_scaled)
  mat_bigsur_scaled_smooth <- t(mat_bigsur_scaled_smooth)

  mat_log_scaled_smooth <- adj_norm %*% t(mat_log_scaled)
  mat_log_scaled_smooth <- t(mat_log_scaled_smooth)

  if (verbose) {
    cat("  Smoothed 4 matrix pairs (raw + scaled)\n")
  }

  # ============================================================================
  # Step 5: Rank-Sum Expression Scoring (Per-Cell Feature Ranking)
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 5: Rank-sum expression scoring\n")
    cat("============================================================\n")
  }

  # Rank features within each cell (use raw data for ranking)
  # OPTIMIZED: Use matrixStats::colRanks (orders of magnitude faster than apply)
  ranks_log <- colRanks(as.matrix(mat_log), ties.method = "average", preserveShape = TRUE)
  ranks_bigsur <- colRanks(as.matrix(mat_bigsur), ties.method = "average", preserveShape = TRUE)
  rank_sum <- ranks_log + ranks_bigsur

  if (verbose) {
    cat(sprintf("  Rank-sum matrix: %d genes × %d cells\n", nrow(rank_sum), ncol(rank_sum)))
  }

  # Pre-compute ranks for Spearman correlation
  # Key insight: Scaling is monotonic, so Spearman(raw) == Spearman(scaled)
  # Therefore we only need TWO sets of ranks: raw and smooth
  if (verbose) {
    cat("  Pre-computing ranks for Spearman correlation...\n")
  }

  # Raw ranks (already computed above for rank_sum)
  ranks_log_raw <- ranks_log
  ranks_bigsur_raw <- ranks_bigsur

  # Smooth ranks (need to compute from smoothed matrices)
  ranks_log_smooth <- colRanks(as.matrix(mat_log_smooth), ties.method = "average", preserveShape = TRUE)
  ranks_bigsur_smooth <- colRanks(as.matrix(mat_bigsur_smooth), ties.method = "average", preserveShape = TRUE)

  if (verbose) {
    cat("  ✓ Ranks computed (2 sets: raw + smooth)\n")
  }

  # ============================================================================
  # Step 6: Define Increments and Conditions
  # ============================================================================

  increments <- c(100, 200, 400, 1000, 2000)
  conditions <- c("raw_unscaled", "raw_scaled", "smooth_unscaled", "smooth_scaled")
  metrics <- c("slope", "pearson_r", "pearson_r2", "spearman_r", "spearman_r2")

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 6: Calculating metrics for all conditions (VECTORIZED)\n")
    cat("============================================================\n")
    cat(sprintf("  Conditions: %d\n", length(conditions)))
    cat(sprintf("  Increments: %s\n", paste(increments, collapse = ", ")))
    cat(sprintf("  Metrics: %d\n", length(metrics)))
    cat(sprintf("  Algorithm: Masked Moment Expansion\n"))
  }

  # ============================================================================
  # Step 7: Calculate Metrics for All Conditions (VECTORIZED)
  # ============================================================================

  # Condition 1: raw_unscaled (uses raw ranks)
  if (verbose) cat("\n  Condition 1/4: raw_unscaled\n")
  results_raw_unscaled <- compute_metrics_for_condition(
    mat_log, mat_bigsur, rank_sum,
    ranks_log_raw, ranks_bigsur_raw,
    increments, verbose = verbose
  )

  # Condition 2: raw_scaled (uses same raw ranks as condition 1)
  if (verbose) cat("\n  Condition 2/4: raw_scaled\n")
  results_raw_scaled <- compute_metrics_for_condition(
    mat_log_scaled, mat_bigsur_scaled, rank_sum,
    ranks_log_raw, ranks_bigsur_raw,
    increments, verbose = verbose
  )

  # Condition 3: smooth_unscaled (uses smooth ranks)
  if (verbose) cat("\n  Condition 3/4: smooth_unscaled\n")
  results_smooth_unscaled <- compute_metrics_for_condition(
    mat_log_smooth, mat_bigsur_smooth, rank_sum,
    ranks_log_smooth, ranks_bigsur_smooth,
    increments, verbose = verbose
  )

  # Condition 4: smooth_scaled (uses same smooth ranks as condition 3)
  if (verbose) cat("\n  Condition 4/4: smooth_scaled\n")
  results_smooth_scaled <- compute_metrics_for_condition(
    mat_log_scaled_smooth, mat_bigsur_scaled_smooth, rank_sum,
    ranks_log_smooth, ranks_bigsur_smooth,
    increments, verbose = verbose
  )

  # ============================================================================
  # Step 8: Write 20 Minimal CSVs
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 7: Writing 20 minimal CSVs\n")
    cat("============================================================\n")
  }

  # Get cell IDs
  cell_ids <- colnames(seurat_object)

  # Store all results in list for easy access
  all_results <- list(
    raw_unscaled = results_raw_unscaled,
    raw_scaled = results_raw_scaled,
    smooth_unscaled = results_smooth_unscaled,
    smooth_scaled = results_smooth_scaled
  )

  # Nested loop: metric → condition → write CSV
  for (metric in metrics) {
    for (condition in conditions) {

      # Extract matrix from all_results[[condition]][[metric]]
      metric_matrix <- all_results[[condition]][[metric]]  # n_increments × n_cells

      # Build dataframe: cell_id, N100, N200, N400, N1000, N2000
      df <- data.frame(
        cell_id = cell_ids,
        N100 = metric_matrix[1, ],
        N200 = metric_matrix[2, ],
        N400 = metric_matrix[3, ],
        N1000 = metric_matrix[4, ],
        N2000 = metric_matrix[5, ]
      )

      # Write CSV
      csv_path <- file.path(output_base, "discordance/metrics", metric, paste0(condition, ".csv"))
      write.csv(df, csv_path, row.names = FALSE)
    }
  }

  if (verbose) {
    cat(sprintf("  Wrote %d CSVs (5 metrics × 4 conditions)\n", length(metrics) * length(conditions)))
  }

  # ============================================================================
  # Step 9: Compute Summary Statistics
  # ============================================================================

  # Extract representative metrics (smooth_unscaled at N=2000)
  slope_2000 <- results_smooth_unscaled$slope[5, ]
  pearson_r2_2000 <- results_smooth_unscaled$pearson_r2[5, ]
  spearman_r2_2000 <- results_smooth_unscaled$spearman_r2[5, ]

  summary_stats <- list(
    n_cells = ncol(seurat_object),
    n_features = length(features),
    increments = increments,
    conditions = conditions,
    metrics = metrics,
    debug_mode = DEBUG_MODE,
    # Representative stats (smooth_unscaled, N=2000)
    mean_slope = mean(slope_2000, na.rm = TRUE),
    median_slope = median(slope_2000, na.rm = TRUE),
    mean_pearson_r2 = mean(pearson_r2_2000, na.rm = TRUE),
    median_pearson_r2 = median(pearson_r2_2000, na.rm = TRUE),
    mean_spearman_r2 = mean(spearman_r2_2000, na.rm = TRUE),
    median_spearman_r2 = median(spearman_r2_2000, na.rm = TRUE)
  )

  if (verbose) {
    cat("\n============================================================\n")
    cat("Discordance Calculation Complete\n")
    cat("============================================================\n")
    cat(sprintf("  Cells: %d\n", summary_stats$n_cells))
    cat(sprintf("  Features: %d\n", summary_stats$n_features))
    cat(sprintf("  CSVs written: %d (5 metrics × 4 conditions)\n", length(metrics) * length(conditions)))
    cat(sprintf("  Representative stats (smooth_unscaled, N=2000):\n"))
    cat(sprintf("    Mean slope: %.3f (median: %.3f)\n", summary_stats$mean_slope, summary_stats$median_slope))
    cat(sprintf("    Mean Pearson R²: %.3f (median: %.3f)\n", summary_stats$mean_pearson_r2, summary_stats$median_pearson_r2))
    cat(sprintf("    Mean Spearman R²: %.3f (median: %.3f)\n", summary_stats$mean_spearman_r2, summary_stats$median_spearman_r2))
    if (DEBUG_MODE) {
      cat("  *** DEBUG MODE WAS ACTIVE ***\n")
    }
    cat("\n")
  }

  # ============================================================================
  # Return Results
  # ============================================================================

  return(list(
    summary_stats = summary_stats
  ))
}
