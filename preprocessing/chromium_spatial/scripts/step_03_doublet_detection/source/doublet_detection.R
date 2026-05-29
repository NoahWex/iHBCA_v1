# ==============================================================================
# STEP 03: Doublet Detection - Core Algorithm
# ==============================================================================
# Extracted from: project/.archive_pre_refactor_20251029/workflows/qc_refactored/scripts/module_3_doublet_detection.R
# Extraction date: 2025-10-29
# Architecture: Manifest-driven preprocessing pipeline
#
# Purpose:
#   - Detect doublets using scDblFinder with cluster context
#   - Apply score-based thresholding via adaptive elbow detection
#   - Perform cluster enrichment analysis with ascending elbow detection
#   - Generate consensus calls using OR logic (score + cluster methods)
#   - Produce diagnostic visualizations and enrichment metrics
#
# Design Principles:
#   - Pure algorithm: NO file I/O operations
#   - All inputs via function parameters
#   - Returns data structures, NOT file paths
#   - No hardcoded paths or parameter values
#   - Wrapper script handles config loading, file I/O, manifest updates
#
# Key Features:
#   - scDblFinder-based doublet scoring with cluster context
#   - Quadratic doublet rate formula: rate * n^2 / 1000
#   - Adaptive threshold detection using loess smoothing and curvature analysis
#   - Cluster enrichment with ascending elbow detection
#   - Hybrid consensus calling (OR logic: score OR cluster)
#   - Edge case handling (small clusters, NA clusters, MAD=0)
#   - 8 output columns for manifest merge
#   - 8 diagnostic plots
#
# Dependencies (must be loaded by wrapper):
#   - Seurat
#   - scDblFinder
#   - SingleCellExperiment
#   - ggplot2
#   - viridis
#   - dplyr
#   - tidyr
#
# Author: Noah Wechter
# Date: 2025-10-29
# ==============================================================================

#' Detect Elbow in Sorted Values Using Loess Smoothing and Curvature Analysis
#'
#' @description
#' Adaptive threshold detection for doublet scores or enrichment values.
#' Uses loess smoothing to reduce noise, then finds the point of maximum
#' curvature (second derivative) within a specified search window.
#'
#' @details
#' Algorithm:
#' 1. Define search window based on expected index and factors
#' 2. Apply loess smoothing (span = 0.2) to reduce noise
#' 3. Compute first derivative (rate of change)
#' 4. Compute second derivative (curvature)
#' 5. Find maximum absolute curvature (elbow point)
#' 6. Validate elbow is within search bounds
#'
#' Edge Cases Handled:
#' - Too few values: Return expected index
#' - Invalid search window: Return expected index
#' - Elbow outside bounds: Return expected index
#'
#' @param sorted_values Numeric vector. Values sorted in DESCENDING order
#'   (for score-based) or ASCENDING order (for enrichment-based).
#'
#' @param expected_idx Integer. Expected position of threshold based on
#'   prior knowledge (e.g., quadratic formula for doublets).
#'
#' @param search_min_factor Numeric. Minimum search position as fraction
#'   of expected_idx. Default: 0.5 (search from 50% of expected).
#'
#' @param search_max_factor Numeric. Maximum search position as fraction
#'   of expected_idx. Default: 2.0 (search up to 200% of expected).
#'
#' @return Named list with 2 elements:
#'   \item{cutoff_idx}{Integer. Index of elbow point in sorted_values}
#'   \item{cutoff_value}{Numeric. Value at elbow point}
#'
#' @examples
#' # Score-based (descending)
#' scores <- sort(runif(1000), decreasing = TRUE)
#' result <- detect_elbow(scores, expected_idx = 50)
#'
#' # Enrichment-based (ascending)
#' enrichment <- sort(runif(100) * 50)
#' result <- detect_elbow(enrichment, expected_idx = 30)
#'
#' @export
detect_elbow <- function(sorted_values,
                         expected_idx,
                         search_min_factor = 0.5,
                         search_max_factor = 2.0) {

  # Define search window
  min_idx <- max(10, floor(expected_idx * search_min_factor))
  max_idx <- min(length(sorted_values), ceiling(expected_idx * search_max_factor))

  # Edge case: invalid window or too few values
  if (max_idx <= min_idx || length(sorted_values) < 3) {
    return(list(
      cutoff_idx = expected_idx,
      cutoff_value = sorted_values[expected_idx]
    ))
  }

  # Extract search region
  search_indices <- min_idx:max_idx
  search_values <- sorted_values[search_indices]

  # Loess smoothing to reduce noise
  ranks <- 1:length(search_values)
  loess_fit <- stats::loess(search_values ~ ranks, span = 0.2)
  smoothed <- loess_fit$fitted

  # Calculate derivatives
  first_diff <- diff(smoothed)
  second_diff <- diff(first_diff)

  # Find maximum curvature
  if (length(second_diff) == 0) {
    elbow_in_search <- 1
  } else {
    elbow_in_search <- which.max(abs(second_diff)) + 1
  }

  # Map back to original indices
  elbow_idx <- search_indices[min(elbow_in_search, length(search_indices))]

  # Validate elbow is within bounds
  if (elbow_idx < min_idx || elbow_idx > max_idx) {
    elbow_idx <- expected_idx
  }

  return(list(
    cutoff_idx = elbow_idx,
    cutoff_value = sorted_values[elbow_idx]
  ))
}


#' Calculate Per-Cluster Doublet Enrichment with Ascending Elbow Detection
#'
#' @description
#' Computes doublet enrichment percentage for each fine cluster and
#' identifies enriched clusters using ascending elbow detection on
#' sorted enrichment values.
#'
#' @details
#' Algorithm:
#' 1. Calculate doublet % per cluster
#' 2. Sort enrichment values in ASCENDING order
#' 3. Apply loess smoothing to enrichment curve
#' 4. Find ascending elbow (steepest rise point)
#' 5. Apply minimum threshold (10% enrichment)
#' 6. Flag clusters above threshold
#'
#' Edge Cases:
#' - No doublet enrichment: Return empty enriched cluster list
#' - <5 non-zero clusters: Use top cluster approach
#' - All enrichment >50%: Use 75th percentile
#' - Loess fails: Fallback to moving average
#'
#' Search Strategy:
#' - Search up to 50% enrichment (filters extreme outliers)
#' - Adaptive loess span based on cluster count
#' - Minimum threshold of 10% (prevent over-flagging)
#'
#' @param cluster_doublet_df Data frame with columns:
#'   - cluster_fine_doublet: Fine cluster ID
#'   - doublet_call_score: Score-based doublet call ("doublet"/"singlet")
#'
#' @param params List containing enrichment parameters (not currently used,
#'   reserved for future config-based thresholds).
#'
#' @return Named list with 2 elements:
#'   \item{cluster_metrics}{data.frame with columns:
#'     \itemize{
#'       \item cluster_fine_doublet: Cluster ID
#'       \item n_cells: Total cells in cluster
#'       \item n_doublets: Doublets in cluster (score-based)
#'       \item doublet_pct: Enrichment percentage
#'     }
#'     Sorted by descending doublet_pct.
#'   }
#'   \item{enriched_clusters}{Character vector of cluster IDs above threshold}
#'
#' @export
calculate_cluster_enrichment <- function(cluster_doublet_df, params = list()) {

  # Calculate per-cluster metrics
  cluster_metrics <- cluster_doublet_df %>%
    dplyr::group_by(cluster_fine_doublet) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_doublets = sum(doublet_call_score == "doublet"),
      doublet_pct = 100 * n_doublets / n_cells,
      .groups = 'drop'
    ) %>%
    dplyr::arrange(desc(doublet_pct))

  # Extract enrichment values
  enrichment_values <- cluster_metrics$doublet_pct
  nonzero_enrichment <- enrichment_values[enrichment_values > 0]

  # Edge case: no enrichment
  if (length(nonzero_enrichment) == 0) {
    return(list(
      cluster_metrics = cluster_metrics,
      enriched_clusters = character(0),
      enrichment_threshold_pct = Inf
    ))
  }

  # Edge case: too few non-zero clusters
  if (length(nonzero_enrichment) < 5) {
    top_enrichments <- sort(nonzero_enrichment, decreasing = TRUE)
    enrichment_threshold_pct <- ifelse(
      top_enrichments[1] > 10,
      top_enrichments[1] - 0.1,
      10.0
    )
    enriched_clusters <- cluster_metrics$cluster_fine_doublet[
      cluster_metrics$doublet_pct >= enrichment_threshold_pct
    ]

    return(list(
      cluster_metrics = cluster_metrics,
      enriched_clusters = as.character(enriched_clusters),
      enrichment_threshold_pct = enrichment_threshold_pct
    ))
  }

  # ASCENDING ELBOW DETECTION
  sorted_enrichment <- sort(nonzero_enrichment)

  # Search up to 50% enrichment
  max_search_value <- 50.0
  search_end_idx <- max(which(sorted_enrichment <= max_search_value))

  # Edge case: insufficient search range
  if (search_end_idx < 5) {
    enrichment_threshold_pct <- max(10.0, quantile(sorted_enrichment, 0.75))
    enriched_clusters <- cluster_metrics$cluster_fine_doublet[
      cluster_metrics$doublet_pct > enrichment_threshold_pct
    ]

    return(list(
      cluster_metrics = cluster_metrics,
      enriched_clusters = as.character(enriched_clusters),
      enrichment_threshold_pct = enrichment_threshold_pct
    ))
  }

  # Extract search region
  search_values <- sorted_enrichment[1:search_end_idx]

  # Adaptive loess span
  x_vals <- seq_along(search_values)
  span_value <- min(0.5, max(0.2, 10 / length(search_values)))

  # Apply loess smoothing with fallback
  tryCatch({
    loess_fit <- loess(search_values ~ x_vals, span = span_value)
    smoothed <- predict(loess_fit, x_vals)
  }, error = function(e) {
    # Fallback: moving average
    window_size <- max(3, floor(length(search_values) * 0.1))
    smoothed <<- stats::filter(
      search_values,
      rep(1/window_size, window_size),
      sides = 2
    )
    smoothed[is.na(smoothed)] <<- search_values[is.na(smoothed)]
  })

  # Calculate curvature
  first_diff <- diff(smoothed)
  second_diff <- diff(first_diff)

  # Edge case: insufficient derivatives
  if (length(second_diff) < 2) {
    enrichment_threshold_pct <- max(10.0, median(search_values))
  } else {
    # Find ascending elbow (maximum second derivative)
    elbow_idx <- which.max(second_diff) + 2
    elbow_idx <- min(elbow_idx, length(search_values))
    enrichment_threshold_pct <- search_values[elbow_idx]

    # Apply minimum threshold
    if (enrichment_threshold_pct < 10.0) {
      enrichment_threshold_pct <- 10.0
    }
  }

  # Identify enriched clusters
  enriched_clusters <- cluster_metrics$cluster_fine_doublet[
    cluster_metrics$doublet_pct > enrichment_threshold_pct
  ]

  return(list(
    cluster_metrics = cluster_metrics,
    enriched_clusters = as.character(enriched_clusters),
    enrichment_threshold_pct = enrichment_threshold_pct
  ))
}


#' Generate Diagnostic Plots for Doublet Detection
#'
#' @description
#' Creates 8 diagnostic plots to visualize doublet detection results:
#' 1. UMAP colored by scDblFinder score (viridis)
#' 2. UMAP colored by detection method (categorical)
#' 3. UMAP colored by final doublet status (binary)
#' 4. Score histogram with threshold line (red)
#' 5. Cluster enrichment barplot with threshold line (red)
#' 6. UMAP of fine clusters (categorical)
#' 7. UMAP of real vs artificial cells (binary)
#' 8. Provenance summary barplot (counts by category)
#'
#' @param seurat_obj Seurat object with UMAP embeddings and doublet metadata
#' @param seurat_with_doublets Seurat object with real + artificial cells
#' @param cluster_metrics Data frame with cluster enrichment metrics
#' @param enrichment_threshold_pct Numeric. Enrichment threshold percentage
#' @param score_threshold Numeric. Score-based threshold value
#' @param sample_id Character. Sample identifier for plot titles
#'
#' @return Named list of 8 ggplot objects
#'
#' @export
generate_diagnostic_plots <- function(seurat_obj,
                                      seurat_with_doublets,
                                      cluster_metrics,
                                      enrichment_threshold_pct,
                                      score_threshold,
                                      sample_id) {

  plots <- list()

  # Plot 1: UMAP colored by scDblFinder score
  plots$score_umap <- ggplot2::ggplot(
    data.frame(
      UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
      UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
      score = seurat_obj$scDblFinder_score
    ),
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = score)
  ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    viridis::scale_color_viridis(option = "viridis", name = "scDblFinder\nScore") +
    ggplot2::labs(
      title = paste0(sample_id, ": scDblFinder Doublet Score"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_minimal()

  # Plot 2: UMAP colored by detection method
  plots$method_umap <- ggplot2::ggplot(
    data.frame(
      UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
      UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
      method = seurat_obj$doublet_call_method
    ),
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = method)
  ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c(
        "passed" = "#440154",            # Purple (singlet)
        "cluster_enrichment" = "#31688E", # Teal (mild)
        "score_threshold" = "#FDE724",    # Green-yellow (moderate)
        "scDblFinder" = "#B73779"         # Red-orange (highest confidence)
      ),
      name = "Detection\nMethod"
    ) +
    ggplot2::labs(
      title = paste0(sample_id, ": Doublet Detection Method"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_minimal()

  # Plot 3: UMAP colored by final status
  plots$status_umap <- ggplot2::ggplot(
    data.frame(
      UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
      UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
      status = ifelse(seurat_obj$doublet_filter_pass, "Singlet", "Doublet")
    ),
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = status)
  ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c("Singlet" = "#2196F3", "Doublet" = "#F44336"),
      name = "Final\nStatus"
    ) +
    ggplot2::labs(
      title = paste0(sample_id, ": Final Doublet Status"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_minimal()

  # Plot 4: Score histogram with threshold
  plots$score_histogram <- ggplot2::ggplot(
    data.frame(score = seurat_obj$scDblFinder_score),
    ggplot2::aes(x = score)
  ) +
    ggplot2::geom_histogram(bins = 50, fill = "#1976D2", alpha = 0.7) +
    ggplot2::geom_vline(
      xintercept = score_threshold,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::labs(
      title = paste0(sample_id, ": scDblFinder Score Distribution"),
      x = "scDblFinder Score",
      y = "Count"
    ) +
    ggplot2::theme_minimal()

  # Plot 5: Cluster enrichment barplot
  plots$enrichment_barplot <- ggplot2::ggplot(
    cluster_metrics,
    ggplot2::aes(
      x = reorder(cluster_fine_doublet, -doublet_pct),
      y = doublet_pct
    )
  ) +
    ggplot2::geom_bar(stat = "identity", fill = "#1976D2", alpha = 0.7) +
    ggplot2::geom_hline(
      yintercept = enrichment_threshold_pct,
      color = "red",
      linetype = "dashed",
      linewidth = 1
    ) +
    ggplot2::labs(
      title = paste0(sample_id, ": Cluster Doublet Enrichment"),
      x = "Fine Cluster",
      y = "Doublet Enrichment (%)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))

  # Plot 6: UMAP of fine clusters
  plots$cluster_umap <- ggplot2::ggplot(
    data.frame(
      UMAP_1 = seurat_obj@reductions$umap@cell.embeddings[, 1],
      UMAP_2 = seurat_obj@reductions$umap@cell.embeddings[, 2],
      cluster = as.factor(seurat_obj$cluster_fine_doublet)
    ),
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = cluster)
  ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::labs(
      title = paste0(sample_id, ": Fine Clustering"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  # Plot 7: UMAP of real vs artificial cells
  plots$real_artificial_umap <- ggplot2::ggplot(
    data.frame(
      UMAP_1 = seurat_with_doublets@reductions$umap@cell.embeddings[, 1],
      UMAP_2 = seurat_with_doublets@reductions$umap@cell.embeddings[, 2],
      cell_type = seurat_with_doublets$cell_type
    ),
    ggplot2::aes(x = UMAP_1, y = UMAP_2, color = cell_type)
  ) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(
      values = c("Real" = "#4CAF50", "Artificial" = "#FF9800"),
      name = "Cell Type"
    ) +
    ggplot2::labs(
      title = paste0(sample_id, ": Real vs Artificial Cells"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::theme_minimal()

  # Plot 8: Provenance summary barplot
  provenance_counts <- table(seurat_obj$doublet_call_method)
  plots$provenance_barplot <- ggplot2::ggplot(
    data.frame(
      method = names(provenance_counts),
      count = as.numeric(provenance_counts)
    ),
    ggplot2::aes(x = method, y = count, fill = method)
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(
      values = c(
        "passed" = "#2E7D32",
        "high_confidence" = "#D32F2F",
        "low_confidence" = "#FFA726"
      )
    ) +
    ggplot2::labs(
      title = paste0(sample_id, ": Doublet Detection Provenance"),
      x = "Detection Category",
      y = "Cell Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")

  return(plots)
}


#' Run Complete Doublet Detection Pipeline for One Sample
#'
#' @description
#' Main algorithm function that processes one sample through the full
#' doublet detection workflow. Runs scDblFinder with cluster context,
#' applies score-based and cluster-based thresholding, generates
#' consensus calls, and produces diagnostic outputs.
#'
#' @details
#' Processing Phases:
#'
#' Phase 1: Setup and Cluster Context
#' - Add cluster assignments from Step 02 to Seurat object
#' - Identify small clusters (<min_cluster_size cells)
#' - Identify NA clusters (if any)
#' - Mark small/NA clusters for exclusion from scDblFinder
#'
#' Phase 2: Basic Preprocessing
#' - Normalize counts (LogNormalize)
#' - Scale variable features
#' - Run PCA (50 dimensions by default)
#' - Build kNN graph (k=15)
#' - Generate UMAP embeddings
#'
#' Phase 3: scDblFinder Execution
#' - Subset to cells not in small/NA clusters
#' - Convert to SingleCellExperiment
#' - Run scDblFinder with cluster_coarse context
#' - Expected doublet rate from config (0.8% default)
#' - Transfer scores back to full Seurat object
#' - Mark excluded cells as singlets (score = 0.0)
#'
#' Phase 4: Score-Based Thresholding
#' - Calculate expected doublet count: floor(rate * n^2 / 1000)
#' - Sort cells by scDblFinder score (descending)
#' - Apply elbow detection within search window
#' - Assign score-based doublet calls
#'
#' Phase 5: Cluster Enrichment Analysis
#' - Fine clustering on real+artificial cells (resolution = 5.0)
#' - Calculate doublet enrichment % per cluster
#' - Apply ascending elbow detection on enrichment
#' - Flag enriched clusters (>threshold)
#' - Assign cluster-based doublet calls
#'
#' Phase 6: Consensus Calling
#' - Combine score + cluster calls using OR logic
#' - Assign provenance labels:
#'   * passed: All methods agree singlet
#'   * high_confidence: scDblFinder direct call OR both methods
#'   * low_confidence: Only one method flagged
#' - Create final doublet_filter_pass column (TRUE = singlet)
#'
#' Phase 7: Output Generation
#' - Extract 8 metadata columns for manifest merge
#' - Generate enrichment metrics table
#' - Create 8 diagnostic plots
#'
#' Edge Cases Handled:
#' - Small clusters: Marked as singlets, excluded from scDblFinder
#' - NA clusters: Marked as singlets
#' - No enrichment: Skip cluster-based filtering
#' - Too few clusters: Use top cluster approach
#' - MAD=0 situations: Handled in loess fallback
#'
#' Output Column Specifications:
#' 1. doublet_filter_pass: Boolean (TRUE = singlet, FALSE = doublet)
#' 2. scDblFinder_score: Numeric (0-1 scale)
#' 3. scDblFinder_class: Factor (singlet/doublet from scDblFinder)
#' 4. cluster_fine_doublet: Character (fine cluster ID)
#' 5. cluster_doublet_enrichment_pct: Numeric (enrichment % for cell's cluster)
#' 6. doublet_call_method: Factor (passed/high_confidence/low_confidence)
#' 7. doublet_call_score: Character (doublet/singlet from score method)
#' 8. doublet_call_cluster: Character (doublet/singlet from cluster method)
#'
#' @param seurat_obj Seurat object created from raw H5 counts.
#'   Must have RNA assay with raw counts.
#'
#' @param cluster_assignments Character vector of cluster assignments from
#'   Step 02 (cluster_coarse column). Length must equal ncol(seurat_obj).
#'   Order must match colnames(seurat_obj).
#'
#' @param variable_features Character vector of variable feature gene names
#'   from Step 01 BigSur analysis. Used for normalization and PCA.
#'
#' @param sample_id Character. Sample identifier for plot titles and logging.
#'
#' @param params Named list of algorithm parameters from module_configs.yaml:
#'   \itemize{
#'     \item doublet_rate: Numeric. Expected doublet rate (e.g., 0.008 for 0.8%)
#'     \item pca_dims: Integer. Number of PCs for basic Seurat processing (50)
#'     \item pca_dims_scDblFinder: Integer. PCs for scDblFinder algorithm (20)
#'     \item k_neighbors_scDblFinder: Integer. kNN neighbors in scDblFinder (100)
#'     \item min_cluster_size: Integer. Minimum cells per cluster (10)
#'     \item deep_cluster_resolution: Numeric. Fine clustering resolution (5.0)
#'     \item score_elbow_search_min_factor: Numeric. Score search min (0.5)
#'     \item score_elbow_search_max_factor: Numeric. Score search max (2.0)
#'   }
#'
#' @return Named list with 4 elements:
#'   \item{seurat_obj}{Updated Seurat object with doublet calls in metadata}
#'   \item{doublet_metadata}{data.frame with 8 output columns (nrow = ncol(seurat_obj)):
#'     \itemize{
#'       \item cell_barcode: Character. Cell barcode
#'       \item sample_id: Character. Sample identifier
#'       \item doublet_filter_pass: Logical. Final doublet status
#'       \item scDblFinder_score: Numeric. scDblFinder score
#'       \item scDblFinder_class: Character. scDblFinder class
#'       \item cluster_fine_doublet: Character. Fine cluster ID
#'       \item cluster_doublet_enrichment_pct: Numeric. Enrichment %
#'       \item doublet_call_method: Character. Detection provenance
#'     }
#'   }
#'   \item{enrichment_metrics}{data.frame with cluster enrichment statistics}
#'   \item{plots}{Named list of 8 ggplot objects}
#'
#' @examples
#' \dontrun{
#' # Load inputs
#' seurat_obj <- CreateSeuratObject(counts = h5_counts)
#' cluster_assignments <- step02_metadata$cluster_coarse
#' variable_features <- readLines("vfs.txt")
#'
#' # Run doublet detection
#' result <- run_doublet_detection(
#'   seurat_obj = seurat_obj,
#'   cluster_assignments = cluster_assignments,
#'   variable_features = variable_features,
#'   sample_id = "Pat1_P1",
#'   params = config$step_03_doublet$algorithm_params
#' )
#'
#' # Access outputs
#' metadata <- result$doublet_metadata
#' plots <- result$plots
#' }
#'
#' @export
run_doublet_detection <- function(seurat_obj,
                                  cluster_assignments,
                                  variable_features,
                                  sample_id,
                                  params) {

  # ==========================================================================
  # PHASE 1: SETUP AND CLUSTER CONTEXT
  # ==========================================================================

  # Add cluster assignments to Seurat object
  seurat_obj$cluster_coarse <- cluster_assignments

  # Validate cluster assignments
  if (any(is.na(seurat_obj$cluster_coarse))) {
    n_na <- sum(is.na(seurat_obj$cluster_coarse))
    stop(sprintf(
      "FATAL: cluster_coarse has %d NAs! Step 02 must run first.",
      n_na
    ))
  }

  # Identify small clusters
  cluster_sizes <- table(seurat_obj$cluster_coarse)
  small_clusters <- names(cluster_sizes[cluster_sizes < params$min_cluster_size])

  if (length(small_clusters) > 0) {
    seurat_obj$small_cluster <- seurat_obj$cluster_coarse %in% small_clusters
  } else {
    seurat_obj$small_cluster <- FALSE
  }

  n_exclude <- sum(seurat_obj$small_cluster)
  n_process <- sum(!seurat_obj$small_cluster)

  # ==========================================================================
  # PHASE 2: BASIC PREPROCESSING
  # ==========================================================================

  DefaultAssay(seurat_obj) <- "RNA"

  # Set variable features
  VariableFeatures(seurat_obj) <- variable_features

  # Normalize
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

  # Scale
  seurat_obj <- ScaleData(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    verbose = FALSE
  )

  # PCA
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(seurat_obj),
    npcs = params$pca_dims,
    verbose = FALSE
  )

  # UMAP
  seurat_obj <- RunUMAP(
    seurat_obj,
    dims = 1:params$pca_dims,
    verbose = FALSE
  )

  # ==========================================================================
  # PHASE 3: SCDBLFINDER EXECUTION
  # ==========================================================================

  # Handle small cluster exclusion
  if (n_exclude > 0) {
    cells_to_process <- colnames(seurat_obj)[!seurat_obj$small_cluster]
    seurat_obj_filtered <- subset(seurat_obj, cells = cells_to_process)
  } else {
    seurat_obj_filtered <- seurat_obj
  }

  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seurat_obj_filtered, assay = "RNA")

  # Adjust PCA dims for scDblFinder
  pca_dims_use <- min(
    params$pca_dims_scDblFinder,
    ncol(seurat_obj_filtered@reductions$pca)
  )

  # CRITICAL FIX 1: Convert cluster_coarse to factor BEFORE scDblFinder
  # This prevents "subscript out of bounds" errors with numeric cluster IDs
  if (is.numeric(sce$cluster_coarse)) {
    message("  CRITICAL FIX: Converting numeric cluster IDs to factor for scDblFinder compatibility")
    sce$cluster_coarse <- factor(sce$cluster_coarse)
  }

  # Run scDblFinder
  set.seed(42)
  sce_with_doublets <- scDblFinder::scDblFinder(
    sce,
    clusters = sce$cluster_coarse,
    returnType = "full",
    dbr = params$doublet_rate,
    score = "weighted",
    dims = pca_dims_use,
    k = params$k_neighbors_scDblFinder,
    includePCs = pca_dims_use,
    nfeatures = VariableFeatures(seurat_obj_filtered),
    verbose = FALSE
  )

  # Create combined Seurat object (real + artificial)
  counts_data <- CreateAssayObject(counts = sce_with_doublets@assays@data$counts)
  meta_data <- as.data.frame(colData(sce_with_doublets))
  seurat_with_doublets <- CreateSeuratObject(
    counts = counts_data,
    meta.data = meta_data
  )

  # Process combined object
  seurat_with_doublets <- NormalizeData(seurat_with_doublets, verbose = FALSE)
  VariableFeatures(seurat_with_doublets) <- VariableFeatures(seurat_obj)
  seurat_with_doublets <- ScaleData(seurat_with_doublets, verbose = FALSE)
  seurat_with_doublets <- RunPCA(
    seurat_with_doublets,
    verbose = FALSE,
    npcs = pca_dims_use
  )
  seurat_with_doublets <- RunUMAP(
    seurat_with_doublets,
    dims = 1:pca_dims_use,
    verbose = FALSE
  )

  # Identify real vs artificial cells
  real_cell_barcodes <- colnames(seurat_obj_filtered)
  artificial_barcodes <- setdiff(colnames(seurat_with_doublets), real_cell_barcodes)

  seurat_with_doublets$cell_type <- "Artificial"
  seurat_with_doublets$cell_type[real_cell_barcodes] <- "Real"

  # Transfer scores to FULL object
  seurat_obj$scDblFinder_class <- NA_character_
  seurat_obj$scDblFinder_score <- NA_real_

  processed_cells <- colnames(seurat_obj) %in% real_cell_barcodes
  seurat_obj$scDblFinder_class[processed_cells] <-
    seurat_with_doublets[, real_cell_barcodes]$class
  seurat_obj$scDblFinder_score[processed_cells] <-
    seurat_with_doublets[, real_cell_barcodes]$score

  # Fill small cluster cells as singlets
  if (n_exclude > 0) {
    small_cluster_cells <- seurat_obj$small_cluster
    seurat_obj$scDblFinder_class[small_cluster_cells] <- "singlet"
    seurat_obj$scDblFinder_score[small_cluster_cells] <- 0.0
  }

  # ==========================================================================
  # PHASE 4: SCORE-BASED THRESHOLDING
  # ==========================================================================

  # Calculate expected doublets using QUADRATIC formula
  n_cells_total <- ncol(seurat_obj)
  expected_doublets <- floor((params$doublet_rate * (n_process^2) / 1000))

  # Sort by score (descending)
  sorted_indices <- order(seurat_obj$scDblFinder_score, decreasing = TRUE)
  sorted_scores <- seurat_obj$scDblFinder_score[sorted_indices]

  # Elbow detection
  elbow_result_score <- detect_elbow(
    sorted_scores,
    expected_doublets,
    params$score_elbow_search_min_factor,
    params$score_elbow_search_max_factor
  )

  # Assign score-based calls
  doublet_barcodes_score <- colnames(seurat_obj)[
    sorted_indices[1:elbow_result_score$cutoff_idx]
  ]
  seurat_obj$doublet_call_score <- ifelse(
    colnames(seurat_obj) %in% doublet_barcodes_score,
    "doublet",
    "singlet"
  )

  # ==========================================================================
  # PHASE 5: FINE CLUSTERING AND ENRICHMENT ANALYSIS
  # ==========================================================================

  # Fine clustering on combined object
  seurat_with_doublets <- FindNeighbors(
    seurat_with_doublets,
    dims = 1:pca_dims_use,
    verbose = FALSE
  )
  seurat_with_doublets <- FindClusters(
    seurat_with_doublets,
    resolution = params$deep_cluster_resolution,
    verbose = FALSE
  )
  seurat_with_doublets$cluster_fine_doublet <- seurat_with_doublets$seurat_clusters

  # Transfer to original object
  seurat_obj$cluster_fine_doublet <-
    seurat_with_doublets[, real_cell_barcodes]$cluster_fine_doublet

  # Calculate cluster enrichment
  cluster_doublet_df <- data.frame(
    cluster_fine_doublet = seurat_obj$cluster_fine_doublet,
    doublet_call_score = seurat_obj$doublet_call_score,
    stringsAsFactors = FALSE
  )

  enrichment_result <- calculate_cluster_enrichment(cluster_doublet_df, params)
  cluster_metrics <- enrichment_result$cluster_metrics
  enriched_clusters <- enrichment_result$enriched_clusters
  enrichment_threshold_pct <- enrichment_result$enrichment_threshold_pct

  # Assign cluster-based calls
  seurat_obj$doublet_call_cluster <- ifelse(
    seurat_obj$cluster_fine_doublet %in% enriched_clusters,
    "doublet",
    "singlet"
  )

  # Add enrichment % to metadata
  enrichment_map <- setNames(
    cluster_metrics$doublet_pct,
    as.character(cluster_metrics$cluster_fine_doublet)
  )
  seurat_obj$cluster_doublet_enrichment_pct <- as.numeric(
    enrichment_map[as.character(seurat_obj$cluster_fine_doublet)]
  )

  # ==========================================================================
  # PHASE 6: CONSENSUS CALLING AND PROVENANCE
  # ==========================================================================

  # Assign provenance (4 categories: detection point tracking)
  # Priority-based: track WHERE in pipeline each doublet was caught
  seurat_obj$doublet_call_method <- dplyr::case_when(
    # scDblFinder algorithm caught it (highest priority)
    seurat_obj$scDblFinder_class == "doublet" ~ "scDblFinder",

    # Score-based threshold caught it (if not scDblFinder)
    seurat_obj$doublet_call_score == "doublet" ~ "score_threshold",

    # Cluster enrichment caught it (if not other methods)
    seurat_obj$doublet_call_cluster == "doublet" ~ "cluster_enrichment",

    # Passed all filters
    TRUE ~ "passed"
  )

  # Final filter (OR logic: doublet if EITHER method flags)
  seurat_obj$doublet_filter_pass <- (
    seurat_obj$doublet_call_score == "singlet" &
    seurat_obj$doublet_call_cluster == "singlet"
  )

  # ==========================================================================
  # PHASE 7: OUTPUT GENERATION
  # ==========================================================================

  # Extract 8 metadata columns for manifest merge
  doublet_metadata <- data.frame(
    cell_barcode = colnames(seurat_obj),
    sample_id = sample_id,
    doublet_filter_pass = seurat_obj$doublet_filter_pass,
    scDblFinder_score = seurat_obj$scDblFinder_score,
    scDblFinder_class = seurat_obj$scDblFinder_class,
    cluster_fine_doublet = seurat_obj$cluster_fine_doublet,
    cluster_doublet_enrichment_pct = seurat_obj$cluster_doublet_enrichment_pct,
    doublet_call_method = seurat_obj$doublet_call_method,
    stringsAsFactors = FALSE
  )

  # Generate diagnostic plots
  plots <- generate_diagnostic_plots(
    seurat_obj = seurat_obj,
    seurat_with_doublets = seurat_with_doublets,
    cluster_metrics = cluster_metrics,
    enrichment_threshold_pct = enrichment_threshold_pct,
    score_threshold = elbow_result_score$cutoff_value,
    sample_id = sample_id
  )

  # Return all outputs
  return(list(
    seurat_obj = seurat_obj,
    doublet_metadata = doublet_metadata,
    enrichment_metrics = cluster_metrics,
    plots = plots
  ))
}
