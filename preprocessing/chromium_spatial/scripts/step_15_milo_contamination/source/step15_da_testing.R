# =============================================================================
# Step 15 — Differential Abundance Testing
# =============================================================================
# MiloR-based differential abundance testing for contamination detection.
# NO FALLBACKS — requires the position column in the Milo object.
#
# Model: ~ position + contamination_status
#   - position:             blocks anatomically distinct positions across patients
#   - contamination_status: binary factor (clean vs contaminated)
#   - contrast:             contamination_statuscontaminated (tests enrichment)
#
# The model blocks on position because contamination is not uniformly
# distributed across anatomical positions — P1 (areola) samples all tend
# to carry epidermal cells and cannot serve as their own clean reference;
# see step15_retention.R for the sensitivity analysis that handles that
# case by excluding position entirely and re-running the DA.
# =============================================================================

suppressPackageStartupMessages({
    library(miloR)
    library(edgeR)
    library(dplyr)
    library(tibble)
})

build_design_matrix <- function(milo_obj, contaminated_samples, exclude_samples = NULL) {
    cat("=== BUILDING DESIGN MATRIX ===\n")

    cell_meta <- as.data.frame(colData(milo_obj))

    required_cols <- c("sample_id", "patient_id", "position")
    missing <- setdiff(required_cols, colnames(cell_meta))
    if (length(missing) > 0) {
        stop("Required columns missing from Milo colData: ", paste(missing, collapse = ", "),
             "\n\nAvailable columns: ", paste(colnames(cell_meta), collapse = ", "),
             "\n\nEnsure Milo object was built with position column.")
    }

    nhood_sample_ids <- colnames(nhoodCounts(milo_obj))
    cat("Samples in nhoodCounts: ", length(nhood_sample_ids), "\n", sep = "")

    design_df <- cell_meta %>%
        dplyr::distinct(sample_id, patient_id, position) %>%
        dplyr::filter(sample_id %in% nhood_sample_ids)

    cat("Unique samples in design: ", nrow(design_df), "\n", sep = "")

    if (!is.null(exclude_samples) && length(exclude_samples) > 0) {
        n_before <- nrow(design_df)
        design_df <- design_df %>%
            dplyr::filter(!sample_id %in% exclude_samples)
        n_excluded <- n_before - nrow(design_df)
        cat("Excluded ", n_excluded, " samples from untestable positions\n", sep = "")
    }

    design_df <- design_df %>%
        dplyr::mutate(
            contamination_status = factor(
                sample_id %in% contaminated_samples,
                levels = c(FALSE, TRUE),
                labels = c("clean", "contaminated")
            ),
            position = factor(position),
            patient_id = factor(patient_id)
        ) %>%
        as.data.frame()

    rownames(design_df) <- design_df$sample_id

    cat("\n--- Design Matrix Summary ---\n")
    cat("Total samples:   ", nrow(design_df), "\n", sep = "")
    cat("Contaminated:    ", sum(design_df$contamination_status == "contaminated"), "\n", sep = "")
    cat("Clean:           ", sum(design_df$contamination_status == "clean"), "\n", sep = "")
    cat("Positions:       ", length(levels(design_df$position)), " unique\n", sep = "")
    cat("Patients:        ", paste(levels(design_df$patient_id), collapse = ", "), "\n", sep = "")

    pos_check <- design_df %>%
        dplyr::group_by(position) %>%
        dplyr::summarise(
            n_clean = sum(contamination_status == "clean"),
            n_contaminated = sum(contamination_status == "contaminated"),
            .groups = "drop"
        ) %>%
        dplyr::filter(n_clean > 0 & n_contaminated > 0)

    cat("Positions with both clean & contaminated: ", nrow(pos_check), "\n", sep = "")

    if (nrow(pos_check) == 0) {
        stop("No positions have both clean and contaminated samples.\n",
             "DA testing requires at least one position with samples from both groups.")
    }

    design_df
}

run_da_testing <- function(milo_obj, design_df) {
    cat("\n=== DIFFERENTIAL ABUNDANCE TESTING ===\n")
    cat("Model: ~ position + contamination_status\n")
    cat("Contrast: contamination_statuscontaminated\n")

    model_mat <- model.matrix(~ position + contamination_status, data = design_df)

    cat("Model matrix dimensions: ", nrow(model_mat), " x ", ncol(model_mat), "\n", sep = "")
    cat("Model terms: ", paste(colnames(model_mat), collapse = ", "), "\n", sep = "")

    cat("\nRunning testNhoods()...\n")

    da_results <- testNhoods(
        milo_obj,
        design = model_mat,
        design.df = design_df,
        model.contrasts = "contamination_statuscontaminated",
        fdr.weighting = "graph-overlap",
        reduced.dim = "SCVI"
    )

    cat("Testing complete.\n")

    da_results <- da_results %>%
        tibble::as_tibble(rownames = "nhood_name") %>%
        dplyr::mutate(nhood_id = dplyr::row_number())

    n_total     <- nrow(da_results)
    n_testable  <- sum(!is.na(da_results$SpatialFDR))
    n_sig       <- sum(da_results$SpatialFDR < 0.05, na.rm = TRUE)
    n_enriched  <- sum(da_results$SpatialFDR < 0.05 & da_results$logFC > 0, na.rm = TRUE)
    n_depleted  <- sum(da_results$SpatialFDR < 0.05 & da_results$logFC < 0, na.rm = TRUE)

    cat("\n--- DA Testing Results ---\n")
    cat("Total neighborhoods:     ", n_total, "\n", sep = "")
    cat("Testable:                ", n_testable, "\n", sep = "")
    cat("Significant (FDR<0.05):  ", n_sig, "\n", sep = "")
    cat("  Enriched (+logFC):     ", n_enriched, "\n", sep = "")
    cat("  Depleted (-logFC):     ", n_depleted, "\n", sep = "")

    da_results
}

classify_contamination <- function(da_results, fdr_threshold = 0.05, logfc_threshold = 0) {
    cat("\n=== CONTAMINATION CLASSIFICATION ===\n")
    cat("FDR threshold:    ", fdr_threshold, "\n", sep = "")
    cat("logFC threshold:  ", logfc_threshold, "\n", sep = "")

    da_results <- da_results %>%
        dplyr::mutate(
            is_contamination_enriched = !is.na(SpatialFDR) &
                                         SpatialFDR < fdr_threshold &
                                         logFC > logfc_threshold,
            confidence = dplyr::case_when(
                is.na(SpatialFDR) ~ "not_testable",
                SpatialFDR < 0.01 & logFC > 1.0 ~ "high",
                SpatialFDR < 0.05 & logFC > 0.5 ~ "moderate",
                SpatialFDR < 0.1  & logFC > 0   ~ "low",
                TRUE ~ "not_significant"
            )
        )

    n_flagged <- sum(da_results$is_contamination_enriched, na.rm = TRUE)
    cat("\nNeighborhoods flagged as contamination-enriched: ", n_flagged, "\n", sep = "")

    conf_summary <- da_results %>%
        dplyr::count(confidence) %>%
        dplyr::mutate(pct = round(100 * n / sum(n), 1))
    cat("\n--- Confidence Distribution ---\n")
    print(as.data.frame(conf_summary), row.names = FALSE)

    da_results
}

summarize_da_results <- function(da_results) {
    list(
        method   = "differential_abundance",
        model    = "~ position + contamination_status",
        contrast = "contamination_statuscontaminated",

        total_neighborhoods    = nrow(da_results),
        testable_neighborhoods = sum(!is.na(da_results$SpatialFDR)),

        sig_enriched_fdr05 = sum(da_results$SpatialFDR < 0.05 & da_results$logFC > 0, na.rm = TRUE),
        sig_depleted_fdr05 = sum(da_results$SpatialFDR < 0.05 & da_results$logFC < 0, na.rm = TRUE),

        flagged_contamination = sum(da_results$is_contamination_enriched, na.rm = TRUE),

        confidence_high     = sum(da_results$confidence == "high",     na.rm = TRUE),
        confidence_moderate = sum(da_results$confidence == "moderate", na.rm = TRUE),
        confidence_low      = sum(da_results$confidence == "low",      na.rm = TRUE),

        logfc_range = list(
            min = round(min(da_results$logFC, na.rm = TRUE), 3),
            max = round(max(da_results$logFC, na.rm = TRUE), 3)
        )
    )
}

print_da_summary <- function(summary) {
    cat("============================================================\n")
    cat("Differential Abundance Testing Summary\n")
    cat("============================================================\n")
    cat("Method:               ", summary$method, "\n", sep = "")
    cat("Model:                ", summary$model, "\n", sep = "")
    cat("Contrast:             ", summary$contrast, "\n", sep = "")
    cat("\n--- Neighborhoods ---\n")
    cat("Total:                ", summary$total_neighborhoods, "\n", sep = "")
    cat("Testable:             ", summary$testable_neighborhoods, "\n", sep = "")
    cat("\n--- Significant (FDR < 0.05) ---\n")
    cat("Enriched (+logFC):    ", summary$sig_enriched_fdr05, "\n", sep = "")
    cat("Depleted (-logFC):    ", summary$sig_depleted_fdr05, "\n", sep = "")
    cat("\n--- Contamination Classification ---\n")
    cat("Flagged:              ", summary$flagged_contamination, "\n", sep = "")
    cat("  High confidence:    ", summary$confidence_high, "\n", sep = "")
    cat("  Moderate:           ", summary$confidence_moderate, "\n", sep = "")
    cat("  Low:                ", summary$confidence_low, "\n", sep = "")
    cat("\n--- Effect Size ---\n")
    cat("logFC range:          [", summary$logfc_range$min, ", ", summary$logfc_range$max, "]\n", sep = "")
    cat("============================================================\n")
}
