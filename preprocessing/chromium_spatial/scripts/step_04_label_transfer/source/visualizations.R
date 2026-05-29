# Step 07 Label Transfer - Visualization Functions
#
# Part of: Spatial HBCA Preprocessing Pipeline
# Step: 07_LabelTransfer
# Type: Algorithm (Pure R functions, return ggplot objects)
#
# Purpose: Generate diagnostic plots for label transfer QC
#
# Dependencies:
#   - Seurat (>= 5.0) **CRITICAL: Seurat v5**
#   - ggplot2
#   - patchwork (for combining plots)
#
# Author: Agent 3 (Algorithm Developer)
# Date: 2025-11-03
# Pattern: TRACE principles (wrapper + source separation)
# Seurat Version: v5

#' Generate UMAP Colored by Cell Type Labels
#'
#' @param query_obj Seurat object. With UMAP and label metadata
#' @param label_col Character. Metadata column name to visualize
#' @param title Character. Plot title (optional)
#'
#' @return ggplot object (does NOT save to file)
#'
#' @examples
#' p1 <- plot_transfer_umap(query_obj, "predicted.celltype.SC", "SC Predictions")
plot_transfer_umap <- function(query_obj, label_col, title = NULL) {

  # Defensive check
  if (!(label_col %in% colnames(query_obj@meta.data))) {
    stop(sprintf("Label column '%s' not found in metadata", label_col))
  }

  # Generate UMAP plot using Seurat v5 DimPlot
  p <- Seurat::DimPlot(
    query_obj,
    reduction = "umap",
    group.by = label_col,
    label = TRUE,
    label.size = 3,
    repel = TRUE
  ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "right",
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(size = 14, face = "bold")
    )

  # Add custom title if provided
  if (!is.null(title)) {
    p <- p + ggplot2::ggtitle(title)
  }

  return(p)
}


#' Generate Violin Plots of Prediction Scores
#'
#' @param query_obj Seurat object. With prediction.score.SC and prediction.score.SN
#'
#' @return ggplot object with two violin plots (SC and SN scores side-by-side)
#'
#' @examples
#' p2 <- plot_prediction_scores(query_obj)
plot_prediction_scores <- function(query_obj) {

  # Defensive checks
  if (!("prediction.score.SC" %in% colnames(query_obj@meta.data))) {
    stop("prediction.score.SC not found in metadata")
  }
  if (!("prediction.score.SN" %in% colnames(query_obj@meta.data))) {
    stop("prediction.score.SN not found in metadata")
  }

  # SC score violin plot using Seurat v5 VlnPlot
  p_sc <- Seurat::VlnPlot(
    query_obj,
    features = "prediction.score.SC",
    pt.size = 0
  ) +
    ggplot2::ggtitle("SC Prediction Scores") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank()
    ) +
    ggplot2::ylim(0, 1)

  # SN score violin plot
  p_sn <- Seurat::VlnPlot(
    query_obj,
    features = "prediction.score.SN",
    pt.size = 0
  ) +
    ggplot2::ggtitle("SN Prediction Scores") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      axis.title.x = ggplot2::element_blank()
    ) +
    ggplot2::ylim(0, 1)

  # Combine with patchwork
  p_combined <- p_sc | p_sn

  return(p_combined)
}


#' Generate Bar Plot of Label Distribution
#'
#' @param query_obj Seurat object. With HBCATransferredLabels.Kumar_2023
#'
#' @return ggplot object with bar chart of label counts
#'
#' @examples
#' p3 <- plot_label_distribution(query_obj)
plot_label_distribution <- function(query_obj) {

  # Defensive check
  if (!("HBCATransferredLabels.Kumar_2023" %in% colnames(query_obj@meta.data))) {
    stop("HBCATransferredLabels.Kumar_2023 not found in metadata")
  }

  # Count cells per label
  label_counts <- as.data.frame(table(query_obj$HBCATransferredLabels.Kumar_2023))
  colnames(label_counts) <- c("Label", "Count")

  # Calculate percentages
  label_counts$Percentage <- label_counts$Count / sum(label_counts$Count) * 100

  # Sort by count (descending)
  label_counts <- label_counts[order(-label_counts$Count), ]

  # Create bar plot
  p <- ggplot2::ggplot(label_counts, ggplot2::aes(x = reorder(Label, Count), y = Count)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", Percentage)),
      hjust = -0.1,
      size = 3
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Final Label Distribution",
      x = "Cell Type",
      y = "Cell Count"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  return(p)
}


#' Generate All 5 Label Transfer Plots
#'
#' Wrapper function that generates all diagnostic plots for a sample.
#' THIS FUNCTION IS ALLOWED TO WRITE FILES (exception to pure algorithm rule).
#'
#' @param query_obj Seurat object. Fully processed with all labels
#' @param output_dir Character. Directory to save plots
#' @param sample_id Character. Sample identifier for filenames
#'
#' @return List of file paths to saved plots
#'
#' @details
#' Generates 5 PNG files:
#' 1. {sample_id}_umap_sc_labels.png - UMAP colored by SC predictions
#' 2. {sample_id}_umap_sn_labels.png - UMAP colored by SN predictions
#' 3. {sample_id}_umap_final_labels.png - UMAP colored by final reconciled labels
#' 4. {sample_id}_prediction_scores.png - Violin plots of SC and SN scores
#' 5. {sample_id}_label_distribution.png - Bar chart of final label counts
#'
#' @examples
#' plot_paths <- generate_label_plots(
#'   query_obj = query_obj,
#'   output_dir = "/path/to/visualizations",
#'   sample_id = "Pat1_P1"
#' )
generate_label_plots <- function(query_obj, output_dir, sample_id) {

  cat("\n=== Generating Label Transfer Plots ===\n")

  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("  Created output directory: %s\n", output_dir))
  }

  # Initialize file path list
  plot_paths <- list()

  # ============================================================================
  # PLOT 1: UMAP SC Labels
  # ============================================================================

  cat("\n[1/5] Generating UMAP SC labels...\n")
  p1 <- plot_transfer_umap(query_obj, "predicted.celltype.SC", "SC Reference Predictions")
  path1 <- file.path(output_dir, paste0(sample_id, "_umap_sc_labels.png"))
  ggplot2::ggsave(path1, p1, width = 12, height = 8, dpi = 300)
  plot_paths$umap_sc <- path1
  cat(sprintf("  Saved: %s\n", path1))

  # ============================================================================
  # PLOT 2: UMAP SN Labels
  # ============================================================================

  cat("\n[2/5] Generating UMAP SN labels...\n")
  p2 <- plot_transfer_umap(query_obj, "predicted.celltype.SN", "SN Reference Predictions")
  path2 <- file.path(output_dir, paste0(sample_id, "_umap_sn_labels.png"))
  ggplot2::ggsave(path2, p2, width = 12, height = 8, dpi = 300)
  plot_paths$umap_sn <- path2
  cat(sprintf("  Saved: %s\n", path2))

  # ============================================================================
  # PLOT 3: UMAP Final Labels
  # ============================================================================

  cat("\n[3/5] Generating UMAP final labels...\n")
  p3 <- plot_transfer_umap(query_obj, "HBCATransferredLabels.Kumar_2023", "Final Reconciled Labels")
  path3 <- file.path(output_dir, paste0(sample_id, "_umap_final_labels.png"))
  ggplot2::ggsave(path3, p3, width = 12, height = 8, dpi = 300)
  plot_paths$umap_final <- path3
  cat(sprintf("  Saved: %s\n", path3))

  # ============================================================================
  # PLOT 4: Prediction Scores
  # ============================================================================

  cat("\n[4/5] Generating prediction score violin plots...\n")
  p4 <- plot_prediction_scores(query_obj)
  path4 <- file.path(output_dir, paste0(sample_id, "_prediction_scores.png"))
  ggplot2::ggsave(path4, p4, width = 12, height = 8, dpi = 300)
  plot_paths$scores <- path4
  cat(sprintf("  Saved: %s\n", path4))

  # ============================================================================
  # PLOT 5: Label Distribution
  # ============================================================================

  cat("\n[5/5] Generating label distribution bar chart...\n")
  p5 <- plot_label_distribution(query_obj)
  path5 <- file.path(output_dir, paste0(sample_id, "_label_distribution.png"))
  ggplot2::ggsave(path5, p5, width = 12, height = 8, dpi = 300)
  plot_paths$distribution <- path5
  cat(sprintf("  Saved: %s\n", path5))

  cat("\n=== All Plots Generated Successfully ===\n")

  return(plot_paths)
}
