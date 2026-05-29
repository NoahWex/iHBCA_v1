# =============================================================================
# Step 15 — Pre-Filter: Identify Untestable Positions
# =============================================================================
# Positions without a clean sample cannot be DA-tested (no counterfactual).
# The pre-filter classifies each position as:
#   NO_CLEAN_REFERENCE : all samples contaminated (untestable via DA)
#   ALL_CLEAN          : no contaminated samples (no contamination signal)
#   MIXED              : has both clean and contaminated (testable via DA)
#
# Cells from NO_CLEAN_REFERENCE positions are flagged directly and removed
# from the DA design matrix; step15_retention.R handles sensitivity analysis
# for whether those cells actually represent global contamination.
# =============================================================================

suppressPackageStartupMessages({
    library(dplyr)
})

identify_untestable_positions <- function(sample_stats) {
    cat("=== PRE-FILTER: IDENTIFYING UNTESTABLE POSITIONS ===\n")

    pos_col <- if ("position" %in% colnames(sample_stats)) "position" else "position_id"

    position_summary <- sample_stats %>%
        dplyr::group_by(.data[[pos_col]]) %>%
        dplyr::summarise(
            n_samples      = dplyr::n(),
            n_contaminated = sum(has_epidermal == TRUE, na.rm = TRUE),
            n_clean        = sum(has_epidermal == FALSE, na.rm = TRUE),
            patients       = paste(unique(patient_id), collapse = ", "),
            .groups        = "drop"
        ) %>%
        dplyr::mutate(
            status = dplyr::case_when(
                n_clean == 0        ~ "NO_CLEAN_REFERENCE",
                n_contaminated == 0 ~ "ALL_CLEAN",
                TRUE                ~ "MIXED"
            )
        ) %>%
        dplyr::rename(position = .data[[pos_col]])

    untestable_positions <- position_summary %>%
        dplyr::filter(status == "NO_CLEAN_REFERENCE") %>%
        dplyr::pull(position)

    testable_positions <- position_summary %>%
        dplyr::filter(status != "NO_CLEAN_REFERENCE") %>%
        dplyr::pull(position)

    untestable_samples <- sample_stats %>%
        dplyr::filter(.data[[pos_col]] %in% untestable_positions) %>%
        dplyr::pull(sample_id)

    cat("\n--- Position Classification Summary ---\n")
    status_counts <- position_summary %>%
        dplyr::count(status) %>%
        dplyr::mutate(pct = round(100 * n / sum(n), 1))
    print(as.data.frame(status_counts), row.names = FALSE)

    cat("\n--- Untestable Positions (NO_CLEAN_REFERENCE) ---\n")
    if (length(untestable_positions) > 0) {
        untestable_df <- position_summary %>%
            dplyr::filter(status == "NO_CLEAN_REFERENCE")
        print(as.data.frame(untestable_df), row.names = FALSE)
        cat("\nSamples excluded from DA testing: ",
            paste(untestable_samples, collapse = ", "), "\n", sep = "")
    } else {
        cat("None — all positions have clean samples for comparison\n")
    }

    list(
        untestable_positions = untestable_positions,
        testable_positions   = testable_positions,
        position_summary     = position_summary,
        untestable_samples   = untestable_samples
    )
}

get_prefilter_cells <- function(milo_obj, untestable_samples) {
    if (length(untestable_samples) == 0) {
        cat("No untestable samples — no pre-filter cells to flag.\n")
        return(data.frame(
            cell_id     = character(),
            sample_id   = character(),
            patient_id  = character(),
            position    = character(),
            is_flagged  = logical(),
            flag_reason = character(),
            stringsAsFactors = FALSE
        ))
    }

    cell_meta <- as.data.frame(SummarizedExperiment::colData(milo_obj))

    if (!"cell_id" %in% colnames(cell_meta)) {
        cell_meta$cell_id <- rownames(cell_meta)
    }

    prefilter_cells <- cell_meta %>%
        dplyr::filter(sample_id %in% untestable_samples) %>%
        dplyr::mutate(
            is_flagged  = TRUE,
            flag_reason = "no_clean_reference"
        ) %>%
        dplyr::select(cell_id, sample_id, patient_id, position, is_flagged, flag_reason)

    cat("\n=== PRE-FILTER FLAGGED CELLS ===\n")
    cat("Samples:  ", paste(unique(prefilter_cells$sample_id), collapse = ", "), "\n", sep = "")
    cat("Cells:    ", nrow(prefilter_cells), "\n", sep = "")

    if (nrow(prefilter_cells) > 0) {
        cat("\nBreakdown by sample:\n")
        sample_counts <- prefilter_cells %>%
            dplyr::count(sample_id) %>%
            dplyr::arrange(dplyr::desc(n))
        print(as.data.frame(sample_counts), row.names = FALSE)
    }

    prefilter_cells
}
