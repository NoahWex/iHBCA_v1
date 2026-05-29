# ============================================================================
# Quality Metrics Calculation
# ============================================================================
#
# Purpose: Calculate quality metrics for label transfer results
#
# This file contains pure algorithm functions for computing QC metrics
# on SC/SN transfer predictions and final reconciled labels.
#
# Dependencies: base R
# ============================================================================

#' Calculate Transfer Quality Metrics
#'
#' Computes comprehensive quality metrics for SC/SN label transfer and
#' reconciliation results.
#'
#' @param query_obj Seurat object with SC/SN predictions and final labels.
#'   Required metadata columns:
#'   - predicted.celltype.SC
#'   - prediction.score.SC
#'   - predicted.celltype.SN
#'   - prediction.score.SN
#'   - HBCATransferredLabels.Kumar_2023 (final labels)
#'   - cluster_coarse (optional, for cluster coherence)
#'
#' @return Named list with quality metrics:
#'   - total_cells: Total number of cells
#'   - sc_transfer_success_rate: Proportion of non-NA SC predictions
#'   - sn_transfer_success_rate: Proportion of non-NA SN predictions
#'   - mean_sc_score: Mean SC prediction score
#'   - median_sc_score: Median SC prediction score
#'   - mean_sn_score: Mean SN prediction score
#'   - median_sn_score: Median SN prediction score
#'   - label_distribution: Named vector of cell counts per final label
#'   - reconciliation_success_rate: Proportion of non-Unassigned cells
#'   - unassigned_count: Number of Unassigned cells
#'   - cluster_label_coherence: Data frame with cluster purity metrics (if cluster_coarse present)
#'
#' @export
#'
#' @examples
#' # metrics <- calculate_transfer_metrics(query_obj)
#' # print(metrics$sc_transfer_success_rate)
#' # print(metrics$label_distribution)
calculate_transfer_metrics <- function(query_obj) {

  # ========================================
  # Input Validation
  # ========================================
  if (!inherits(query_obj, "Seurat")) {
    stop("query_obj must be a Seurat object")
  }

  required_cols <- c("predicted.celltype.SC", "prediction.score.SC",
                     "predicted.celltype.SN", "prediction.score.SN",
                     "HBCATransferredLabels.Kumar_2023")

  missing_cols <- setdiff(required_cols, colnames(query_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Query object missing metadata columns: %s",
                 paste(missing_cols, collapse = ", ")))
  }

  message("Calculating transfer quality metrics...")

  # ========================================
  # Extract Metadata
  # ========================================
  metadata <- query_obj@meta.data
  total_cells <- ncol(query_obj)

  # ========================================
  # SC Transfer Metrics
  # ========================================
  sc_predictions <- metadata$predicted.celltype.SC
  sc_scores <- metadata$prediction.score.SC

  sc_non_na <- sum(!is.na(sc_predictions))
  sc_success_rate <- sc_non_na / total_cells

  sc_mean_score <- mean(sc_scores, na.rm = TRUE)
  sc_median_score <- median(sc_scores, na.rm = TRUE)
  sc_min_score <- min(sc_scores, na.rm = TRUE)
  sc_max_score <- max(sc_scores, na.rm = TRUE)

  message(sprintf("SC Transfer: %.1f%% success rate, mean score: %.3f",
                  sc_success_rate * 100, sc_mean_score))

  # ========================================
  # SN Transfer Metrics
  # ========================================
  sn_predictions <- metadata$predicted.celltype.SN
  sn_scores <- metadata$prediction.score.SN

  sn_non_na <- sum(!is.na(sn_predictions))
  sn_success_rate <- sn_non_na / total_cells

  sn_mean_score <- mean(sn_scores, na.rm = TRUE)
  sn_median_score <- median(sn_scores, na.rm = TRUE)
  sn_min_score <- min(sn_scores, na.rm = TRUE)
  sn_max_score <- max(sn_scores, na.rm = TRUE)

  message(sprintf("SN Transfer: %.1f%% success rate, mean score: %.3f",
                  sn_success_rate * 100, sn_mean_score))

  # ========================================
  # Reconciliation Metrics
  # ========================================
  final_labels <- metadata$HBCATransferredLabels.Kumar_2023
  label_distribution <- table(final_labels)

  unassigned_count <- sum(final_labels == "Unassigned", na.rm = TRUE)
  reconciliation_success_rate <- (total_cells - unassigned_count) / total_cells

  message(sprintf("Reconciliation: %.1f%% success rate (%d unassigned)",
                  reconciliation_success_rate * 100, unassigned_count))
  message(sprintf("Unique final labels: %d", length(unique(final_labels))))

  # ========================================
  # Cluster-Label Coherence (optional)
  # ========================================
  cluster_coherence <- NULL
  if ("cluster_coarse" %in% colnames(metadata)) {
    message("Calculating cluster-label coherence...")
    cluster_coherence <- calculate_cluster_coherence(query_obj)
  } else {
    message("Skipping cluster coherence (no cluster_coarse column)")
  }

  # ========================================
  # Compile Results
  # ========================================
  metrics <- list(
    # Overall
    total_cells = total_cells,

    # SC Transfer
    sc_transfer_success_rate = sc_success_rate,
    sc_predictions_non_na = sc_non_na,
    sc_prediction_score_mean = sc_mean_score,
    sc_prediction_score_median = sc_median_score,
    sc_prediction_score_min = sc_min_score,
    sc_prediction_score_max = sc_max_score,

    # SN Transfer
    sn_transfer_success_rate = sn_success_rate,
    sn_predictions_non_na = sn_non_na,
    sn_prediction_score_mean = sn_mean_score,
    sn_prediction_score_median = sn_median_score,
    sn_prediction_score_min = sn_min_score,
    sn_prediction_score_max = sn_max_score,

    # Reconciliation
    reconciliation_success_rate = reconciliation_success_rate,
    unassigned_count = unassigned_count,
    label_distribution = as.list(label_distribution),

    # Cluster coherence (if available)
    cluster_label_coherence = cluster_coherence
  )

  message("Quality metrics calculation complete")
  return(metrics)
}


#' Calculate Cluster-Label Coherence
#'
#' For each cluster, calculates the purity (how homogeneous the cluster is
#' with respect to final cell type labels).
#'
#' @param query_obj Seurat object with cluster_coarse and HBCATransferredLabels.Kumar_2023
#'
#' @return Data frame with columns:
#'   - cluster_id: Cluster identifier
#'   - total_cells: Number of cells in cluster
#'   - dominant_label: Most common label in cluster
#'   - dominant_count: Number of cells with dominant label
#'   - purity_score: Proportion of cells with dominant label (0-1)
#'   - entropy: Shannon entropy of label distribution (lower = more homogeneous)
#'   - n_unique_labels: Number of unique labels in cluster
#'
#' @export
#'
#' @examples
#' # coherence <- calculate_cluster_coherence(query_obj)
#' # mean(coherence$purity_score)  # Average cluster purity
calculate_cluster_coherence <- function(query_obj) {

  # Validate inputs
  if (!inherits(query_obj, "Seurat")) {
    stop("query_obj must be a Seurat object")
  }

  metadata <- query_obj@meta.data

  if (!all(c("cluster_coarse", "HBCATransferredLabels.Kumar_2023") %in% colnames(metadata))) {
    stop("Query object missing cluster_coarse or HBCATransferredLabels.Kumar_2023")
  }

  clusters <- metadata$cluster_coarse
  labels <- metadata$HBCATransferredLabels.Kumar_2023

  # Get unique clusters
  unique_clusters <- sort(unique(clusters))
  n_clusters <- length(unique_clusters)

  # Initialize results
  coherence_results <- data.frame(
    cluster_id = unique_clusters,
    total_cells = integer(n_clusters),
    dominant_label = character(n_clusters),
    dominant_count = integer(n_clusters),
    purity_score = numeric(n_clusters),
    entropy = numeric(n_clusters),
    n_unique_labels = integer(n_clusters),
    stringsAsFactors = FALSE
  )

  # Calculate metrics for each cluster
  for (i in seq_along(unique_clusters)) {
    clust <- unique_clusters[i]

    # Get labels for this cluster
    cluster_labels <- labels[clusters == clust]
    total <- length(cluster_labels)

    # Label distribution
    label_counts <- table(cluster_labels)
    n_unique <- length(label_counts)

    # Dominant label
    dominant_idx <- which.max(label_counts)
    dominant_label <- names(label_counts)[dominant_idx]
    dominant_count <- as.integer(label_counts[dominant_idx])

    # Purity (proportion with dominant label)
    purity <- dominant_count / total

    # Shannon entropy
    proportions <- as.numeric(label_counts) / total
    entropy <- -sum(proportions * log2(proportions + 1e-10))  # Add small value to avoid log(0)

    # Store results
    coherence_results$total_cells[i] <- total
    coherence_results$dominant_label[i] <- dominant_label
    coherence_results$dominant_count[i] <- dominant_count
    coherence_results$purity_score[i] <- purity
    coherence_results$entropy[i] <- entropy
    coherence_results$n_unique_labels[i] <- n_unique
  }

  # Summary statistics
  message("\nCluster-Label Coherence Summary:")
  message(sprintf("  - Total clusters: %d", n_clusters))
  message(sprintf("  - Mean purity: %.3f", mean(coherence_results$purity_score)))
  message(sprintf("  - Median purity: %.3f", median(coherence_results$purity_score)))
  message(sprintf("  - Min purity: %.3f (cluster %s)",
                  min(coherence_results$purity_score),
                  coherence_results$cluster_id[which.min(coherence_results$purity_score)]))
  message(sprintf("  - Max purity: %.3f (cluster %s)",
                  max(coherence_results$purity_score),
                  coherence_results$cluster_id[which.max(coherence_results$purity_score)]))

  return(coherence_results)
}


#' Generate Summary Statistics Table
#'
#' Creates a human-readable summary table from metrics list.
#'
#' @param metrics List returned by calculate_transfer_metrics()
#'
#' @return Data frame with metric names and formatted values
#' @export
format_metrics_summary <- function(metrics) {

  summary_df <- data.frame(
    Metric = character(),
    Value = character(),
    stringsAsFactors = FALSE
  )

  # Add metrics
  summary_df <- rbind(summary_df, data.frame(
    Metric = "Total Cells",
    Value = as.character(metrics$total_cells)
  ))

  summary_df <- rbind(summary_df, data.frame(
    Metric = "SC Transfer Success Rate",
    Value = sprintf("%.2f%%", metrics$sc_transfer_success_rate * 100)
  ))

  summary_df <- rbind(summary_df, data.frame(
    Metric = "SN Transfer Success Rate",
    Value = sprintf("%.2f%%", metrics$sn_transfer_success_rate * 100)
  ))

  summary_df <- rbind(summary_df, data.frame(
    Metric = "Mean SC Score",
    Value = sprintf("%.3f", metrics$sc_prediction_score_mean)
  ))

  summary_df <- rbind(summary_df, data.frame(
    Metric = "Mean SN Score",
    Value = sprintf("%.3f", metrics$sn_prediction_score_mean)
  ))

  summary_df <- rbind(summary_df, data.frame(
    Metric = "Reconciliation Success Rate",
    Value = sprintf("%.2f%%", metrics$reconciliation_success_rate * 100)
  ))

  summary_df <- rbind(summary_df, data.frame(
    Metric = "Unassigned Cells",
    Value = as.character(metrics$unassigned_count)
  ))

  if (!is.null(metrics$cluster_label_coherence)) {
    summary_df <- rbind(summary_df, data.frame(
      Metric = "Mean Cluster Purity",
      Value = sprintf("%.3f", mean(metrics$cluster_label_coherence$purity_score))
    ))
  }

  return(summary_df)
}
