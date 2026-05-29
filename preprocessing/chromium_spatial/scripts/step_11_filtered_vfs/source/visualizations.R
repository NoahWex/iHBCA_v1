# ==============================================================================
# STEP 04: Post-QC Variable Features - Visualizations
# ==============================================================================
# Purpose:
#   - Per-patient scatter plots for VF retention analysis
#   - Plot 1: Baseline vs Filtered VF counts
#   - Plot 2: Cell retention vs VF retention
#
# Design:
#   - Purely descriptive (no QC threshold visualizations)
#   - ggplot2-based scatter plots with sample-level granularity
#
# Author: Noah Wechter
# Date: 2025-10-30
# ==============================================================================

#' Create VF Count Scatter Plot
#'
#' @description
#' Generates a scatter plot comparing baseline VF counts (from Step 01) against
#' filtered VF counts (after Step 04 QC filtering). Includes a 1:1 reference line
#' to indicate where counts would be unchanged. Each point represents a sample,
#' colored by sample_id.
#'
#' @param patient_metrics A data frame containing per-sample metrics with columns:
#'   - `n_vfs_baseline`: Number of VFs at baseline (from Step 01)
#'   - `n_vfs_filtered`: Number of VFs after QC filtering (from Step 04)
#'   - `sample_id`: Sample identifier for coloring points
#' @param patient_id Character string identifying the patient. Used in plot title.
#' @param output_path Character string specifying the full path where the plot
#'   should be saved (e.g., "/path/to/plot.png").
#'
#' @return
#' Invisibly returns the ggplot object. Saves a PNG file at `output_path`
#' with dimensions 8×6 inches at 300 dpi.
#'
#' @details
#' The plot includes:
#' - X-axis: Baseline VF count (Step 01)
#' - Y-axis: Filtered VF count (Post-QC, Step 04)
#' - Color: Sample ID for point distinction
#' - Reference line: Dashed gray 1:1 line (y = x) indicating no change
#' - Points: Size 3, alpha 0.7 for transparency
#' - Theme: Minimal with legend positioned to the right
#'
#' @examples
#' \dontrun{
#' # Sample data
#' patient_metrics <- data.frame(
#'   sample_id = c("Sample_A", "Sample_B", "Sample_C"),
#'   n_vfs_baseline = c(2000, 1800, 2100),
#'   n_vfs_filtered = c(1850, 1600, 1950)
#' )
#'
#' # Create and save plot
#' create_vf_count_scatter(
#'   patient_metrics = patient_metrics,
#'   patient_id = "PATIENT_001",
#'   output_path = "/path/to/output/vf_count_scatter.png"
#' )
#' }
#'
#' @export
create_vf_count_scatter <- function(patient_metrics, patient_id, output_path) {
  library(ggplot2)

  p <- ggplot(patient_metrics, aes(x = n_vfs_baseline, y = n_vfs_filtered, color = sample_id)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = paste0(patient_id, ": VF Counts (Baseline vs Filtered)"),
         x = "Baseline VF Count (Step 01)",
         y = "Filtered VF Count (Post-QC)",
         color = "Sample") +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(output_path, p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

#' Create Cell and VF Retention Scatter Plot
#'
#' @description
#' Generates a scatter plot examining the relationship between the percentage of
#' cells retained after overall QC filtering and the percentage of variable
#' features (VFs) retained after Step 04 filtering. Each point represents a sample,
#' colored by sample_id. This plot is useful for understanding how cell-level QC
#' impacts feature retention.
#'
#' @param patient_metrics A data frame containing per-sample metrics with columns:
#'   - `cells_retained_pct`: Percentage of cells passing overall QC filtering
#'   - `retention_pct`: Percentage of VFs retained after Step 04 filtering
#'   - `sample_id`: Sample identifier for coloring points
#' @param patient_id Character string identifying the patient. Used in plot title.
#' @param output_path Character string specifying the full path where the plot
#'   should be saved (e.g., "/path/to/plot.png").
#'
#' @return
#' Invisibly returns the ggplot object. Saves a PNG file at `output_path`
#' with dimensions 8×6 inches at 300 dpi.
#'
#' @details
#' The plot includes:
#' - X-axis: Percentage of cells retained (overall QC)
#' - Y-axis: Percentage of VFs retained (Step 04)
#' - Color: Sample ID for point distinction
#' - Points: Size 3, alpha 0.7 for transparency
#' - Theme: Minimal with legend positioned to the right
#'
#' @examples
#' \dontrun{
#' # Sample data
#' patient_metrics <- data.frame(
#'   sample_id = c("Sample_A", "Sample_B", "Sample_C"),
#'   cells_retained_pct = c(92.5, 88.3, 95.1),
#'   retention_pct = c(87.5, 81.2, 91.4)
#' )
#'
#' # Create and save plot
#' create_retention_scatter(
#'   patient_metrics = patient_metrics,
#'   patient_id = "PATIENT_001",
#'   output_path = "/path/to/output/retention_scatter.png"
#' )
#' }
#'
#' @export
create_retention_scatter <- function(patient_metrics, patient_id, output_path) {
  library(ggplot2)

  p <- ggplot(patient_metrics, aes(x = cells_retained_pct, y = retention_pct, color = sample_id)) +
    geom_point(size = 3, alpha = 0.7) +
    labs(title = paste0(patient_id, ": Cell vs VF Retention"),
         x = "% Cells Retained (Overall QC)",
         y = "% VFs Retained",
         color = "Sample") +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(output_path, p, width = 8, height = 6, dpi = 300)
  invisible(p)
}

#' Create Expression Boxplots for VF Categories
#'
#' @description
#' Generates boxplots comparing mean expression values across different variable
#' feature (VF) categories: lost, retained, and gained. Each point represents a
#' sample, with boxplots showing the distribution across samples within a patient.
#'
#' @param patient_metrics A data frame containing per-sample metrics with columns:
#'   - `sample_id`: Sample identifier for labeling points
#'   - `mean_expr_lost`: Mean expression of VFs lost during QC
#'   - `mean_expr_retained`: Mean expression of VFs retained after QC
#'   - `mean_expr_gained`: Mean expression of VFs gained during QC
#' @param patient_id Character string identifying the patient. Used in plot title.
#' @param output_path Character string specifying the full path where the plot
#'   should be saved (e.g., "/path/to/plot.png").
#'
#' @return
#' Invisibly returns the ggplot object. Saves a PNG file at `output_path`
#' with dimensions 8×6 inches at 300 dpi.
#'
#' @details
#' The plot includes:
#' - X-axis: VF Category (3 levels: "Lost VFs", "Retained VFs", "Gained VFs")
#' - Y-axis: Mean Expression (linear scale)
#' - Boxplots: Show distribution across samples
#' - Points: Individual sample values (jittered for visibility)
#' - Colors:
#'   - Lost VFs: Red (#E41A1C) - features removed by QC
#'   - Retained VFs: Blue (#377EB8) - stable features across QC
#'   - Gained VFs: Green (#4DAF4A) - features revealed by QC
#'
#' Interpretation:
#' - Low expression in "Lost VFs": QC successfully removed low-signal features
#' - High expression in "Retained VFs": QC preserved robust biological signal
#' - High expression in "Gained VFs": QC revealed previously masked features
#'
#' @examples
#' \dontrun{
#' # Sample data
#' patient_metrics <- data.frame(
#'   sample_id = c("Sample_A", "Sample_B", "Sample_C"),
#'   mean_expr_lost = c(1.85, 1.92, 1.78),
#'   mean_expr_retained = c(2.53, 2.48, 2.61),
#'   mean_expr_gained = c(2.21, 2.18, 2.35)
#' )
#'
#' # Create and save plot
#' create_expression_boxplots(
#'   patient_metrics = patient_metrics,
#'   patient_id = "PATIENT_001",
#'   output_path = "/path/to/output/expression_boxplots.png"
#' )
#' }
#'
#' @export
create_expression_boxplots <- function(patient_metrics, patient_id, output_path) {
  library(ggplot2)
  library(tidyr)

  # Reshape data to long format for ggplot
  expr_data <- patient_metrics %>%
    select(sample_id, mean_expr_lost, mean_expr_retained, mean_expr_gained) %>%
    pivot_longer(
      cols = c(mean_expr_lost, mean_expr_retained, mean_expr_gained),
      names_to = "category",
      values_to = "mean_expr"
    ) %>%
    mutate(
      category = case_when(
        category == "mean_expr_lost" ~ "Lost VFs",
        category == "mean_expr_retained" ~ "Retained VFs",
        category == "mean_expr_gained" ~ "Gained VFs"
      ),
      category = factor(category, levels = c("Lost VFs", "Retained VFs", "Gained VFs"))
    )

  # Create boxplot with jittered points
  p <- ggplot(expr_data, aes(x = category, y = mean_expr, fill = category)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    scale_fill_manual(values = c("Lost VFs" = "#E41A1C",
                                   "Retained VFs" = "#377EB8",
                                   "Gained VFs" = "#4DAF4A")) +
    labs(
      title = paste0(patient_id, ": Mean Expression by VF Category"),
      x = "VF Category",
      y = "Mean Expression",
      fill = "Category"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold")
    )

  # Save plot
  ggsave(output_path, p, width = 8, height = 6, dpi = 300)
  invisible(p)
}
