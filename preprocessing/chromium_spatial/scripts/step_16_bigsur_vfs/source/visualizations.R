# =============================================================================
# Step 16 — BigSur VF Retention Visualizations
# =============================================================================
# Per-patient diagnostic plots for VF retention and expression distribution
# across categories (retained / lost / gained). Purely descriptive —
# invoked optionally by the preview step (step 18) or by ad-hoc QC reports;
# the core pipeline does not depend on these figures.
# =============================================================================

#' Baseline vs filtered VF count scatter (per sample).
#' Reference line is y = x (no change under filtering).
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

#' Cell retention vs VF retention scatter (per sample).
#' Probes whether VF instability tracks with cell removal intensity.
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

#' Boxplots of mean expression per VF category: lost, retained, gained.
#' Interpretation: low mean expression of Lost VFs suggests QC removed
#' low-signal features; high Retained indicates stable biology across the
#' filter; high Gained indicates features unmasked by removing confounding
#' variation.
create_expression_boxplots <- function(patient_metrics, patient_id, output_path) {
    library(ggplot2)
    library(tidyr)

    expr_data <- patient_metrics %>%
        select(sample_id, mean_expr_lost, mean_expr_retained, mean_expr_gained) %>%
        pivot_longer(
            cols      = c(mean_expr_lost, mean_expr_retained, mean_expr_gained),
            names_to  = "category",
            values_to = "mean_expr"
        ) %>%
        mutate(
            category = case_when(
                category == "mean_expr_lost"     ~ "Lost VFs",
                category == "mean_expr_retained" ~ "Retained VFs",
                category == "mean_expr_gained"   ~ "Gained VFs"
            ),
            category = factor(category, levels = c("Lost VFs", "Retained VFs", "Gained VFs"))
        )

    p <- ggplot(expr_data, aes(x = category, y = mean_expr, fill = category)) +
        geom_boxplot(alpha = 0.7, outlier.shape = NA) +
        geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
        scale_fill_manual(values = c("Lost VFs"     = "#E41A1C",
                                     "Retained VFs" = "#377EB8",
                                     "Gained VFs"   = "#4DAF4A")) +
        labs(
            title = paste0(patient_id, ": Mean Expression by VF Category"),
            x     = "VF Category",
            y     = "Mean Expression",
            fill  = "Category"
        ) +
        theme_minimal() +
        theme(
            legend.position = "none",
            axis.text.x     = element_text(size = 11),
            axis.title      = element_text(size = 12, face = "bold"),
            plot.title      = element_text(size = 14, face = "bold")
        )

    ggsave(output_path, p, width = 8, height = 6, dpi = 300)
    invisible(p)
}
