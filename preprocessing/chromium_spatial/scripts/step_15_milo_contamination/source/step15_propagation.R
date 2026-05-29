# =============================================================================
# Step 15 — Cell Flag Propagation from Neighborhoods
# =============================================================================
# Propagate neighborhood-level DA flags to individual cells using weighted
# voting. Each cell votes based on its neighborhood memberships, weighted
# by significance (1/SpatialFDR). A cell is flagged when the weighted
# fraction of flagged-neighborhood memberships exceeds `flag_threshold`
# (default 0.95). The threshold is deliberately conservative so that cells
# belonging to any mixed-signal neighborhood are retained.
#
# Pattern source: hbca_analysis/analyses/AbundanceAnalysis/
#                 position_clustering_20250828/full_object_analysis/
#                 split_01_full_object.Rmd
# =============================================================================

suppressPackageStartupMessages({
    library(miloR)
    library(Matrix)
    library(dplyr)
    library(tibble)
})

propagate_to_cells <- function(milo_obj,
                                da_results,
                                flag_threshold = 0.95,
                                fdr_threshold = 0.05) {

    cat("=== PROPAGATING FLAGS TO CELLS ===\n")
    cat("Flag threshold:     ", flag_threshold, "\n")
    cat("FDR threshold:      ", fdr_threshold, "\n")

    nhood_matrix <- nhoods(milo_obj)

    n_cells  <- nrow(nhood_matrix)
    n_nhoods <- ncol(nhood_matrix)

    cat("Cells:              ", n_cells, "\n")
    cat("Neighborhoods:      ", n_nhoods, "\n")

    cell_ids <- rownames(nhood_matrix)
    if (is.null(cell_ids)) cell_ids <- colnames(milo_obj)

    if (!"is_contamination_enriched" %in% colnames(da_results)) {
        da_results <- da_results %>%
            dplyr::mutate(
                is_contamination_enriched = !is.na(SpatialFDR) &
                                            SpatialFDR < fdr_threshold &
                                            logFC > 0
            )
    }

    if (n_nhoods != nrow(da_results)) {
        warning("Neighborhood count mismatch: matrix has ", n_nhoods,
                " but da_results has ", nrow(da_results), " rows. Using min.")
        n_nhoods <- min(n_nhoods, nrow(da_results))
        nhood_matrix <- nhood_matrix[, seq_len(n_nhoods)]
        da_results   <- da_results[seq_len(n_nhoods), ]
    }

    nhood_weights <- ifelse(
        !is.na(da_results$SpatialFDR) & da_results$SpatialFDR < fdr_threshold,
        1 / pmax(da_results$SpatialFDR, 1e-10),
        0
    )

    nhood_flags <- as.numeric(da_results$is_contamination_enriched)
    nhood_flags[is.na(nhood_flags)] <- 0

    cat("\nNeighborhood flags: ", sum(nhood_flags), " flagged\n")
    cat("Neighborhoods with weight > 0: ", sum(nhood_weights > 0), "\n")

    nhood_matrix_sparse <- as(nhood_matrix, "dgCMatrix")

    cell_n_nhoods <- Matrix::rowSums(nhood_matrix_sparse > 0)

    flag_matrix <- nhood_matrix_sparse %*% Matrix::Diagonal(x = nhood_flags)
    cell_n_flagged <- Matrix::rowSums(flag_matrix > 0)

    weighted_flag_vec <- nhood_weights * nhood_flags
    numerator   <- as.numeric(nhood_matrix_sparse %*% weighted_flag_vec)
    denominator <- as.numeric(nhood_matrix_sparse %*% nhood_weights)

    weighted_flag_score <- ifelse(denominator > 0, numerator / denominator, 0)

    nhood_logfc <- da_results$logFC
    nhood_logfc[is.na(nhood_logfc)] <- 0

    logfc_weighted_sum <- as.numeric(nhood_matrix_sparse %*% (nhood_weights * nhood_logfc))
    weighted_logfc <- ifelse(denominator > 0, logfc_weighted_sum / denominator, NA_real_)

    nhood_neglogFDR <- -log10(pmax(da_results$SpatialFDR, 1e-10))
    nhood_neglogFDR[is.na(nhood_neglogFDR) | is.infinite(nhood_neglogFDR)] <- 0

    neglogFDR_sum <- as.numeric(nhood_matrix_sparse %*% nhood_neglogFDR)
    mean_neglogFDR <- ifelse(cell_n_nhoods > 0, neglogFDR_sum / cell_n_nhoods, NA_real_)

    is_flagged <- weighted_flag_score >= flag_threshold

    flag_reason <- dplyr::case_when(
        cell_n_nhoods == 0     ~ "no_neighborhoods",
        cell_n_flagged == 0    ~ "no_flagged_nhoods",
        is_flagged             ~ paste0("weighted_vote_", round(weighted_flag_score, 2)),
        weighted_flag_score > 0 ~ paste0("minority_flagged_", round(weighted_flag_score, 2)),
        TRUE                   ~ "clean"
    )

    cell_flags <- tibble::tibble(
        cell_id             = cell_ids,
        n_nhoods            = as.integer(cell_n_nhoods),
        n_flagged_nhoods    = as.integer(cell_n_flagged),
        weighted_flag_score = round(weighted_flag_score, 4),
        weighted_logfc      = round(weighted_logfc, 4),
        mean_neglogFDR      = round(mean_neglogFDR, 4),
        is_flagged          = is_flagged,
        flag_reason         = flag_reason
    )

    n_flagged_cells <- sum(cell_flags$is_flagged)
    pct_flagged    <- round(100 * n_flagged_cells / n_cells, 2)

    cat("\n--- Cell Propagation Summary ---\n")
    cat("Total cells:        ", n_cells, "\n")
    cat("Flagged cells:      ", n_flagged_cells, " (", pct_flagged, "%)\n", sep = "")
    cat("Cells with >0 flagged nhoods: ", sum(cell_flags$n_flagged_nhoods > 0), "\n")

    cell_flags
}

summarize_cell_flags <- function(cell_flags) {

    total_cells   <- nrow(cell_flags)
    flagged_cells <- sum(cell_flags$is_flagged)

    flagged_data       <- cell_flags %>% dplyr::filter(is_flagged)
    has_flagged_nhoods <- cell_flags %>% dplyr::filter(n_flagged_nhoods > 0)

    summary_list <- list(
        total_cells                    = total_cells,
        flagged_cells                  = flagged_cells,
        flagged_pct                    = round(100 * flagged_cells / total_cells, 2),

        cells_with_flagged_nhoods     = nrow(has_flagged_nhoods),
        cells_with_flagged_nhoods_pct = round(100 * nrow(has_flagged_nhoods) / total_cells, 2),

        mean_nhoods_per_cell  = round(mean(cell_flags$n_nhoods), 1),
        median_nhoods_per_cell = median(cell_flags$n_nhoods),

        mean_flagged_nhoods_among_flagged = if (flagged_cells > 0) {
            round(mean(flagged_data$n_flagged_nhoods), 1)
        } else NA,

        weighted_score_quantiles = if (nrow(has_flagged_nhoods) > 0) {
            quantile(has_flagged_nhoods$weighted_flag_score,
                     probs = c(0, 0.25, 0.5, 0.75, 1))
        } else NULL
    )

    cat("\n============================================================\n")
    cat("Cell Flag Summary\n")
    cat("============================================================\n")
    cat("Total cells:                    ", summary_list$total_cells, "\n")
    cat("Flagged cells:                  ", summary_list$flagged_cells,
        " (", summary_list$flagged_pct, "%)\n", sep = "")
    cat("\n--- Neighborhood Membership ---\n")
    cat("Mean nhoods/cell:               ", summary_list$mean_nhoods_per_cell, "\n")
    cat("Median nhoods/cell:             ", summary_list$median_nhoods_per_cell, "\n")
    cat("\n--- Cells with Flagged Neighborhoods ---\n")
    cat("Cells with >=1 flagged nhood:   ", summary_list$cells_with_flagged_nhoods,
        " (", summary_list$cells_with_flagged_nhoods_pct, "%)\n", sep = "")

    if (!is.null(summary_list$weighted_score_quantiles)) {
        cat("\n--- Weighted Score Distribution (among cells with flagged nhoods) ---\n")
        cat("Min:    ", summary_list$weighted_score_quantiles[1], "\n")
        cat("Q1:     ", summary_list$weighted_score_quantiles[2], "\n")
        cat("Median: ", summary_list$weighted_score_quantiles[3], "\n")
        cat("Q3:     ", summary_list$weighted_score_quantiles[4], "\n")
        cat("Max:    ", summary_list$weighted_score_quantiles[5], "\n")
    }
    cat("============================================================\n")

    invisible(summary_list)
}
