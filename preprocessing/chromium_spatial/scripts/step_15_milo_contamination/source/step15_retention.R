# =============================================================================
# Step 15 — Untestable Position Retention (Sensitivity Analysis via Exclusion)
# =============================================================================
# Positions without clean counterparts (e.g. P1 areola, where every sample
# carries epidermal cells) cannot be DA-tested directly. For each such
# position we run a sensitivity DA that EXCLUDES the position's samples
# entirely and ask: do the previously flagged neighborhoods retain their
# significance without this position's signal?
#
#   Loses significance when excluded  -> the position drove the flag ->
#                                        classify as RETAIN (position-specific)
#   Remains significant when excluded -> signal is shared across positions ->
#                                        classify as FLAG (global contamination)
#
# Exclusion (rather than reclassifying the position as clean) is used to
# avoid the coefficient-absorption problem where blocking the position
# soaks up the contamination signal into the intercept term.
# =============================================================================

suppressPackageStartupMessages({
    library(miloR)
    library(edgeR)
    library(dplyr)
    library(tibble)
})

run_sensitivity_da_for_position <- function(milo_obj, cell_meta, contaminated_samples,
                                             position_to_exclude, exclude_samples = NULL) {
    cat("\n--- Sensitivity DA: EXCLUDING ", position_to_exclude, " ---\n", sep = "")

    position_samples <- cell_meta %>%
        dplyr::filter(position == position_to_exclude) %>%
        dplyr::pull(sample_id) %>%
        unique()

    if (length(position_samples) == 0) {
        cat("No samples found for position ", position_to_exclude, "\n", sep = "")
        return(NULL)
    }

    cat("Samples to EXCLUDE: ", paste(position_samples, collapse = ", "), "\n", sep = "")

    nhood_sample_ids <- colnames(nhoodCounts(milo_obj))

    design_df <- cell_meta %>%
        dplyr::distinct(sample_id, patient_id, position) %>%
        dplyr::filter(sample_id %in% nhood_sample_ids) %>%
        dplyr::filter(!sample_id %in% position_samples)

    if (!is.null(exclude_samples)) {
        other_untestable <- setdiff(exclude_samples, position_samples)
        if (length(other_untestable) > 0) {
            design_df <- design_df %>%
                dplyr::filter(!sample_id %in% other_untestable)
        }
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

    pos_check <- design_df %>%
        dplyr::group_by(position) %>%
        dplyr::summarise(
            n_clean = sum(contamination_status == "clean"),
            n_contaminated = sum(contamination_status == "contaminated"),
            .groups = "drop"
        ) %>%
        dplyr::filter(n_clean > 0 & n_contaminated > 0)

    if (nrow(pos_check) == 0) {
        cat("No testable positions after exclusion — skipping\n")
        return(NULL)
    }

    cat("Samples remaining: ", nrow(design_df), " | Contaminated: ",
        sum(design_df$contamination_status == "contaminated"),
        " | Clean: ", sum(design_df$contamination_status == "clean"), "\n", sep = "")

    remaining_positions <- unique(design_df$position)
    if (length(remaining_positions) > 1) {
        model_mat <- model.matrix(~ position + contamination_status, data = design_df)
    } else {
        model_mat <- model.matrix(~ contamination_status, data = design_df)
    }

    da_sensitivity <- tryCatch({
        testNhoods(
            milo_obj,
            design = model_mat,
            design.df = design_df,
            model.contrasts = "contamination_statuscontaminated",
            fdr.weighting = "graph-overlap",
            reduced.dim = "SCVI"
        )
    }, error = function(e) {
        cat("Error in sensitivity DA: ", e$message, "\n", sep = "")
        return(NULL)
    })

    if (is.null(da_sensitivity)) return(NULL)

    da_sensitivity %>%
        tibble::as_tibble(rownames = "nhood_name") %>%
        dplyr::mutate(
            nhood_id = dplyr::row_number(),
            sensitivity_position = position_to_exclude
        ) %>%
        dplyr::rename_with(
            ~ paste0(.x, "_sens_", position_to_exclude),
            .cols = c("logFC", "logCPM", "F", "PValue", "FDR", "SpatialFDR")
        )
}

get_position_samples <- function(cell_meta, position) {
    cell_meta %>%
        dplyr::filter(position == !!position) %>%
        dplyr::pull(sample_id) %>%
        unique()
}

analyze_untestable_retention <- function(milo_obj, da_original, cell_meta,
                                          contaminated_samples,
                                          untestable_positions,
                                          untestable_samples = NULL,
                                          fdr_threshold = 0.05) {
    cat("\n============================================================\n")
    cat("UNTESTABLE POSITION RETENTION — SENSITIVITY ANALYSIS\n")
    cat("============================================================\n")
    cat("Positions without clean counterparts: ", paste(untestable_positions, collapse = ", "), "\n")
    cat("Method: For each position, compare DA with position excluded\n")
    cat("============================================================\n")

    if (length(untestable_positions) == 0) {
        cat("No untestable positions — all positions have clean samples.\n")
        return(list(
            sensitivity_results = list(),
            comparison = da_original %>%
                dplyr::mutate(
                    retain_action     = NA_character_,
                    retain_reason     = NA_character_,
                    retain_confidence = NA_character_,
                    retain_positions  = NA_character_
                ),
            summary = list(method = "skipped", reason = "no_untestable_positions")
        ))
    }

    sensitivity_results <- list()

    for (pos in untestable_positions) {
        cat("\n")
        da_sens <- run_sensitivity_da_for_position(
            milo_obj             = milo_obj,
            cell_meta            = cell_meta,
            contaminated_samples = contaminated_samples,
            position_to_exclude  = pos,
            exclude_samples      = untestable_samples
        )

        if (!is.null(da_sens)) {
            sensitivity_results[[pos]] <- da_sens
        }
    }

    if (length(sensitivity_results) == 0) {
        cat("\nAll sensitivity DAs failed — using fallback.\n")
        return(list(
            sensitivity_results = list(),
            comparison = da_original %>%
                dplyr::mutate(
                    retain_action     = ifelse(is_contamination_enriched, "flag", NA_character_),
                    retain_reason     = "sensitivity_failed",
                    retain_confidence = "LOW",
                    retain_positions  = NA_character_
                ),
            summary = list(method = "fallback", reason = "all_sensitivity_failed")
        ))
    }

    comparison <- compare_multi_sensitivity(
        da_original         = da_original,
        sensitivity_results = sensitivity_results,
        fdr_threshold       = fdr_threshold
    )

    summary <- summarize_retention_analysis(comparison, untestable_positions)

    list(
        sensitivity_results = sensitivity_results,
        comparison          = comparison,
        summary             = summary
    )
}

compare_multi_sensitivity <- function(da_original, sensitivity_results, fdr_threshold = 0.05) {
    cat("\n=== COMPARING ORIGINAL VS SENSITIVITY DAs ===\n")

    comparison <- da_original
    comparison$positions_lost_sig <- ""
    comparison$min_logfc_change   <- NA_real_
    comparison$any_lost_sig       <- FALSE

    for (pos in names(sensitivity_results)) {
        sens_df <- sensitivity_results[[pos]]

        logfc_col <- paste0("logFC_sens_", pos)
        fdr_col   <- paste0("SpatialFDR_sens_", pos)

        comparison <- comparison %>%
            dplyr::left_join(
                sens_df %>% dplyr::select(nhood_id, !!logfc_col, !!fdr_col),
                by = "nhood_id"
            )

        comparison <- comparison %>%
            dplyr::mutate(
                sig_original = !is.na(SpatialFDR) & SpatialFDR < fdr_threshold & logFC > 0,
                sig_sens     = !is.na(.data[[fdr_col]]) & .data[[fdr_col]] < fdr_threshold & .data[[logfc_col]] > 0,
                lost_sig_this = sig_original & !sig_sens,

                positions_lost_sig = ifelse(
                    lost_sig_this,
                    paste0(positions_lost_sig, ifelse(positions_lost_sig == "", "", ", "), pos),
                    positions_lost_sig
                ),

                logfc_change_this = .data[[logfc_col]] - logFC,
                min_logfc_change  = pmin(min_logfc_change, logfc_change_this, na.rm = TRUE),

                any_lost_sig = any_lost_sig | lost_sig_this
            ) %>%
            dplyr::select(-sig_original, -sig_sens, -lost_sig_this, -logfc_change_this)
    }

    comparison <- comparison %>%
        dplyr::mutate(
            retain_action = dplyr::case_when(
                !is_contamination_enriched                              ~ NA_character_,
                any_lost_sig                                            ~ "retain",
                !is.na(min_logfc_change) & min_logfc_change < -0.5      ~ "retain",
                is_contamination_enriched                               ~ "flag",
                TRUE                                                    ~ NA_character_
            ),

            retain_reason = dplyr::case_when(
                is.na(retain_action)                                    ~ NA_character_,
                any_lost_sig                                            ~ paste0("lost_significance_", positions_lost_sig),
                !is.na(min_logfc_change) & min_logfc_change < -0.5      ~ "large_logfc_drop",
                retain_action == "flag"                                 ~ "retained_significance_shared",
                TRUE                                                    ~ NA_character_
            ),

            retain_confidence = dplyr::case_when(
                is.na(retain_action)                                    ~ NA_character_,
                any_lost_sig                                            ~ "HIGH",
                !is.na(min_logfc_change) & min_logfc_change < -0.5      ~ "MODERATE",
                retain_action == "flag"                                 ~ "HIGH",
                TRUE                                                    ~ "LOW"
            ),

            retain_positions = ifelse(positions_lost_sig == "", NA_character_, positions_lost_sig)
        )

    n_retain <- sum(comparison$retain_action == "retain", na.rm = TRUE)
    n_flag   <- sum(comparison$retain_action == "flag",   na.rm = TRUE)

    cat("Contamination-enriched neighborhoods: ",
        sum(comparison$is_contamination_enriched, na.rm = TRUE), "\n")
    cat("RETAINED (position-specific):         ", n_retain, "\n")
    cat("FLAGGED (shared contamination):       ", n_flag, "\n")

    if (n_retain > 0) {
        cat("\nRetained by position:\n")
        for (pos in names(sensitivity_results)) {
            n_pos <- sum(grepl(pos, comparison$retain_positions, fixed = TRUE), na.rm = TRUE)
            if (n_pos > 0) cat("  ", pos, ": ", n_pos, " neighborhoods\n", sep = "")
        }
    }

    comparison
}

summarize_retention_analysis <- function(comparison, untestable_positions) {
    n_enriched <- sum(comparison$is_contamination_enriched, na.rm = TRUE)
    n_retain   <- sum(comparison$retain_action == "retain", na.rm = TRUE)
    n_flag     <- sum(comparison$retain_action == "flag",   na.rm = TRUE)

    by_position <- list()
    for (pos in untestable_positions) {
        n_pos <- sum(grepl(pos, comparison$retain_positions, fixed = TRUE), na.rm = TRUE)
        by_position[[pos]] <- n_pos
    }

    list(
        method      = "sensitivity_analysis",
        description = "Compare DA with each untestable position excluded",
        untestable_positions = untestable_positions,

        contamination_enriched    = n_enriched,
        retained_position_specific = n_retain,
        flagged_shared             = n_flag,

        retained_by_position = by_position,

        retained_lost_significance = sum(grepl("lost_significance", comparison$retain_reason), na.rm = TRUE),
        retained_logfc_drop        = sum(comparison$retain_reason == "large_logfc_drop", na.rm = TRUE)
    )
}

analyze_retention_sensitivity <- function(milo_obj, da_original, cell_meta,
                                           contaminated_samples,
                                           exclude_samples = NULL,
                                           fdr_threshold = 0.05) {
    all_positions <- unique(cell_meta$position)

    contaminated_positions <- cell_meta %>%
        dplyr::filter(sample_id %in% contaminated_samples) %>%
        dplyr::pull(position) %>%
        unique()

    clean_positions <- cell_meta %>%
        dplyr::filter(!sample_id %in% contaminated_samples) %>%
        dplyr::pull(position) %>%
        unique()

    untestable_positions <- setdiff(contaminated_positions, clean_positions)

    cat("Detected untestable positions: ", paste(untestable_positions, collapse = ", "), "\n")

    analyze_untestable_retention(
        milo_obj             = milo_obj,
        da_original          = da_original,
        cell_meta            = cell_meta,
        contaminated_samples = contaminated_samples,
        untestable_positions = untestable_positions,
        untestable_samples   = exclude_samples,
        fdr_threshold        = fdr_threshold
    )
}
