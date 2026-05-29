# ============================================================================
# Label Reconciliation Algorithm
# ============================================================================
#
# Purpose: Merge SC and SN predictions using rule-based reconciliation
#
# This file contains pure algorithm functions for reconciling SC and SN
# predictions into a final harmonized cell type label.
#
# Dependencies: base R
# ============================================================================

#' Reconcile SC and SN Labels
#'
#' Merges SC and SN predictions using a rule-based mapping file with wildcard
#' support. Rules are applied in order (first match wins).
#'
#' @param sc_pred Character vector. SC predictions for each cell
#' @param sn_pred Character vector. SN predictions for each cell
#' @param mapping_file_path Character. Path to CSV file with reconciliation rules.
#'   Expected columns: SC_Label_Input, SN_Label_Input, Final_Label_Output,
#'   Priority_Rule, Notes
#'
#' @return Character vector of final reconciled labels (same length as inputs)
#'   Returns NULL on error.
#'
#' @details
#' The mapping file supports wildcard patterns:
#'   - "*" matches any value (including NA)
#'   - Exact string matches are case-insensitive
#'   - Rules are evaluated in file order (first match wins)
#'   - Should always include a fallback rule: "*,*,Unassigned,Default_NoMatch"
#'
#' @export
#'
#' @examples
#' # sc_pred <- c("Luminal_HR", "T cells", "Fibroblasts", NA)
#' # sn_pred <- c("lumHR", "Tcells", "fibro", "adipo")
#' # mapping_path <- "project/ReferenceDatasets/HBCA/sc_sn_label_mapping.txt"
#' # final_labels <- reconcile_labels(sc_pred, sn_pred, mapping_path)
reconcile_labels <- function(sc_pred, sn_pred, mapping_file_path) {

  # ========================================
  # Input Validation
  # ========================================
  if (!file.exists(mapping_file_path)) {
    stop(sprintf("Mapping file not found: %s", mapping_file_path))
  }

  if (length(sc_pred) != length(sn_pred)) {
    stop(sprintf("sc_pred and sn_pred must have same length (got %d vs %d)",
                 length(sc_pred), length(sn_pred)))
  }

  n_cells <- length(sc_pred)
  message(sprintf("Starting reconciliation for %d cells", n_cells))

  # ========================================
  # Load Mapping Rules
  # ========================================
  mapping_rules <- tryCatch({
    message(sprintf("Loading mapping rules from: %s", mapping_file_path))
    rules <- read.csv(mapping_file_path, stringsAsFactors = FALSE)

    # Validate columns
    required_cols <- c("SC_Label_Input", "SN_Label_Input", "Final_Label_Output",
                      "Priority_Rule", "Notes")
    missing_cols <- setdiff(required_cols, colnames(rules))
    if (length(missing_cols) > 0) {
      stop(sprintf("Mapping file missing columns: %s",
                   paste(missing_cols, collapse = ", ")))
    }

    message(sprintf("Loaded %d reconciliation rules", nrow(rules)))
    rules

  }, error = function(e) {
    message(sprintf("Error loading mapping file: %s", e$message))
    return(NULL)
  })

  if (is.null(mapping_rules)) {
    return(NULL)
  }

  # Verify fallback rule exists
  has_fallback <- any(mapping_rules$SC_Label_Input == "*" &
                      mapping_rules$SN_Label_Input == "*")
  if (!has_fallback) {
    warning("Mapping file missing fallback rule (*,*,Unassigned) - some cells may not match")
  }

  # ========================================
  # Apply Reconciliation Rules
  # ========================================
  final_labels <- character(n_cells)
  rule_applied <- character(n_cells)  # Track which rule was used

  # Process each cell
  for (i in seq_len(n_cells)) {
    sc_val <- sc_pred[i]
    sn_val <- sn_pred[i]

    # Normalize inputs (lowercase, handle NA)
    sc_lower <- normalize_for_matching(sc_val)
    sn_lower <- normalize_for_matching(sn_val)

    # Find first matching rule
    matched <- FALSE
    for (rule_idx in seq_len(nrow(mapping_rules))) {
      rule <- mapping_rules[rule_idx, ]

      # Check if rule matches
      sc_match <- matches_pattern(sc_lower, rule$SC_Label_Input)
      sn_match <- matches_pattern(sn_lower, rule$SN_Label_Input)

      if (sc_match && sn_match) {
        final_labels[i] <- rule$Final_Label_Output
        rule_applied[i] <- rule$Priority_Rule
        matched <- TRUE
        break  # First match wins
      }
    }

    # If no match found (should not happen with fallback rule)
    if (!matched) {
      final_labels[i] <- "Unassigned"
      rule_applied[i] <- "NoMatch"
      warning(sprintf("Cell %d: No rule matched (SC=%s, SN=%s)",
                     i, sc_val, sn_val))
    }
  }

  # ========================================
  # Report Reconciliation Summary
  # ========================================
  message("\nReconciliation Summary:")
  message(sprintf("  - Total cells: %d", n_cells))
  message(sprintf("  - Unique final labels: %d", length(unique(final_labels))))

  # Count by rule type
  rule_counts <- table(rule_applied)
  message("\nRules Applied:")
  print(rule_counts)

  # Count by final label
  label_counts <- table(final_labels)
  message("\nFinal Label Distribution:")
  print(label_counts)

  # Calculate success rate (non-Unassigned)
  success_rate <- sum(final_labels != "Unassigned") / n_cells
  message(sprintf("\nReconciliation success rate: %.2f%%", success_rate * 100))

  # Return final labels
  return(final_labels)
}


#' Normalize Value for Pattern Matching
#'
#' Converts a value to lowercase for matching. Handles NA values.
#'
#' @param val Character or NA. Input value
#' @return Character. Normalized value ("NA" for NA inputs, lowercase otherwise)
#' @keywords internal
normalize_for_matching <- function(val) {
  if (is.na(val) || is.null(val) || val == "") {
    return("NA")
  }
  return(tolower(as.character(val)))
}


#' Check if Value Matches Pattern
#'
#' Checks if a value matches a pattern (with wildcard support).
#'
#' @param value Character. Normalized value to match
#' @param pattern Character. Pattern from mapping file (supports "*" wildcard)
#' @return Logical. TRUE if matches, FALSE otherwise
#' @keywords internal
matches_pattern <- function(value, pattern) {
  # Wildcard matches everything
  if (pattern == "*") {
    return(TRUE)
  }

  # Exact match (case-insensitive, already normalized)
  pattern_lower <- tolower(as.character(pattern))
  return(value == pattern_lower)
}


#' Reconcile Labels with Rule Tracking
#'
#' Extended version that returns both final labels and rule information.
#' Useful for QC and debugging.
#'
#' @param sc_pred Character vector. SC predictions
#' @param sn_pred Character vector. SN predictions
#' @param mapping_file_path Character. Path to mapping CSV
#'
#' @return Data frame with columns:
#'   - final_label: Final reconciled label
#'   - rule_applied: Rule type that was applied
#'   - rule_priority: Priority value from mapping
#'
#' @export
reconcile_labels_with_tracking <- function(sc_pred, sn_pred, mapping_file_path) {

  # Load mapping rules
  if (!file.exists(mapping_file_path)) {
    stop(sprintf("Mapping file not found: %s", mapping_file_path))
  }

  mapping_rules <- read.csv(mapping_file_path, stringsAsFactors = FALSE)
  n_cells <- length(sc_pred)

  # Initialize results
  results <- data.frame(
    final_label = character(n_cells),
    rule_applied = character(n_cells),
    rule_priority = character(n_cells),
    stringsAsFactors = FALSE
  )

  # Process each cell
  for (i in seq_len(n_cells)) {
    sc_lower <- normalize_for_matching(sc_pred[i])
    sn_lower <- normalize_for_matching(sn_pred[i])

    # Find first matching rule
    for (rule_idx in seq_len(nrow(mapping_rules))) {
      rule <- mapping_rules[rule_idx, ]

      sc_match <- matches_pattern(sc_lower, rule$SC_Label_Input)
      sn_match <- matches_pattern(sn_lower, rule$SN_Label_Input)

      if (sc_match && sn_match) {
        results$final_label[i] <- rule$Final_Label_Output
        results$rule_applied[i] <- rule$Priority_Rule
        results$rule_priority[i] <- as.character(rule_idx)
        break
      }
    }

    # Fallback
    if (results$final_label[i] == "") {
      results$final_label[i] <- "Unassigned"
      results$rule_applied[i] <- "NoMatch"
      results$rule_priority[i] <- "NA"
    }
  }

  return(results)
}
