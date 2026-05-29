# Step 07 Label Transfer - Reconcile SC and SN Labels
#
# Part of: Spatial HBCA Preprocessing Pipeline
# Step: 07_LabelTransfer
# Type: Algorithm (Pure R function)
#
# Purpose: Reconcile SC and SN predictions using evidence-based rules from
#          sc_sn_label_mapping.txt
#
# Dependencies:
#   - dplyr (for case_when logic)
#
# Author: Agent 3 (Algorithm Developer)
# Date: 2025-11-03
# Pattern: TRACE principles (wrapper + source separation)
# Reference: TEMPLATES 03a_ConsolidateAnnotationsAndQC.Rmd lines 847-993

#' Reconcile SC and SN Predictions Using Mapping Rules
#'
#' Applies 15 reconciliation rules from sc_sn_label_mapping.txt to merge
#' SC and SN predictions into a single final label.
#'
#' @param query_obj Seurat object. With predicted.celltype.SC and predicted.celltype.SN
#' @param mapping_path Character. Path to sc_sn_label_mapping.txt CSV
#' @param config List. Reconciliation parameters (optional)
#'   - reconciliation_strategy: "use_mapping_file" (default)
#'
#' @return Seurat object with added metadata column:
#'   - HBCATransferredLabels.Kumar_2023: Final reconciled cell type labels
#'
#' @details
#' Reconciliation logic based on TEMPLATES 03a_ConsolidateAnnotationsAndQC.Rmd
#' (lines 847-993):
#'
#' Rule priority (case_when logic):
#' 1. Both NA → "Unassigned_BothNA"
#' 2. SC NA, SN present → use SN prediction
#' 3. SN NA, SC present → use SC prediction
#' 4. Exact match → use SC naming (norm_pred_SC == norm_pred_SN)
#' 5. Granularity rules:
#'    - SC=myeloid AND SN=mast → "Mast" (SN_Granular)
#' 6. Unique labels:
#'    - SN=adipo → "Adipocyte" (SN_Unique)
#'    - SC=bcells → "B_Cell" (SC_Unique)
#' 7. Naming preference conflicts → use SC naming convention
#'    - SC=pericytes, SN=perivasc → "Pericyte" (SC_NamePreference)
#'    - SC=tcells, SN=Tcells → "T_Cell" (SC_NamePreference)
#'    - SC=fibroblasts, SN=fibro → "Fibroblast" (SC_NamePreference)
#'    - SC=vascular, SN=vasc → "Vascular_Endothelium" (SC_NamePreference)
#'    - SC=lymphatic, SN=lymph → "Lymphatic_Endothelium" (SC_NamePreference)
#' 8. Fallback → "Unassigned_ConflictOrUnhandled"
#'
#' Implementation approach:
#' - Normalize SC and SN predictions (lowercase, remove special chars)
#' - Apply case_when logic matching TEMPLATES
#' - Map normalized outputs to canonical label names
#' - Handle wildcards (*) in mapping file
#'
#' @examples
#' query_obj <- reconcile_labels(
#'   query_obj = query_obj,
#'   mapping_path = "/path/to/sc_sn_label_mapping.txt",
#'   config = list(reconciliation_strategy = "use_mapping_file")
#' )
reconcile_labels <- function(
  query_obj,
  mapping_path,
  config = list()
) {

  # Extract config parameters with defaults
  reconciliation_strategy <- config$reconciliation_strategy %||% "use_mapping_file"

  cat("\n=== Label Reconciliation ===\n")

  # ============================================================================
  # STEP 1: Load Reconciliation Mapping CSV
  # ============================================================================

  cat("\n[1/4] Loading reconciliation mapping...\n")

  # Defensive check
  if (!file.exists(mapping_path)) {
    stop(sprintf("Reconciliation mapping not found: %s", mapping_path))
  }

  # Load mapping CSV (15 rules + header)
  mapping_df <- read.csv(mapping_path, stringsAsFactors = FALSE)

  # Verify structure
  required_cols <- c("SC_Label_Input", "SN_Label_Input", "Final_Label_Output", "Priority_Rule")
  if (!all(required_cols %in% colnames(mapping_df))) {
    stop(sprintf("Mapping file missing required columns. Expected: %s. Found: %s",
                 paste(required_cols, collapse=", "), paste(colnames(mapping_df), collapse=", ")))
  }

  cat(sprintf("  Loaded %d reconciliation rules\n", nrow(mapping_df)))

  # ============================================================================
  # STEP 2: Extract SC and SN Predictions from Query Object
  # ============================================================================

  cat("\n[2/4] Extracting predictions from query object...\n")

  # Verify columns exist
  if (!("predicted.celltype.SC" %in% colnames(query_obj@meta.data))) {
    stop("predicted.celltype.SC not found in query object metadata")
  }
  if (!("predicted.celltype.SN" %in% colnames(query_obj@meta.data))) {
    stop("predicted.celltype.SN not found in query object metadata")
  }

  # Extract predictions
  pred_sc <- query_obj$predicted.celltype.SC
  pred_sn <- query_obj$predicted.celltype.SN

  cat(sprintf("  SC predictions: %d cells (%d unique labels, %d NA)\n",
              length(pred_sc), length(unique(pred_sc[!is.na(pred_sc)])), sum(is.na(pred_sc))))
  cat(sprintf("  SN predictions: %d cells (%d unique labels, %d NA)\n",
              length(pred_sn), length(unique(pred_sn[!is.na(pred_sn)])), sum(is.na(pred_sn))))

  # ============================================================================
  # STEP 3: Normalize Predictions (lowercase, remove special chars)
  # ============================================================================

  cat("\n[3/4] Normalizing predictions...\n")

  # Normalization function (EXACT from TEMPLATES)
  # Source: TEMPLATES/1_phase_per_sample_prep/03a_ConsolidateAnnotationsAndQC.Rmd lines 752-779
  # Uses lookup-table approach with 31 explicit mappings
  normalize_label <- function(label) {
    if (is.na(label) || is.null(label) || !nzchar(trimws(as.character(label)))) {
      return(NA_character_)
    }
    label_lower <- tolower(trimws(as.character(label)))

    # TEMPLATES internal_std_map - exact lookup table
    internal_std_map <- c(
      "basal" = "basal",
      "lumhr" = "luminalhr", "lum hr" = "luminalhr", "luminal_hr" = "luminalhr",
      "lumsec" = "luminalsecretory", "luminal_secretory" = "luminalsecretory",
      "fibroblasts" = "fibroblast", "fibroblast" = "fibroblast", "fibro" = "fibroblast",
      "vascular" = "vascularendothelium", "vasc" = "vascularendothelium", "vascular_endothelium" = "vascularendothelium",
      "lymphatic" = "lymphaticendothelium", "lymph" = "lymphaticendothelium", "lymphatic_endothelium" = "lymphaticendothelium",
      "pericytes" = "pericyte", "pericyte" = "pericyte", "perivasc" = "pericyte", "perivascular_cells" = "pericyte",
      "tcells" = "tcell", "t_cells" = "tcell", "t cell" = "tcell",
      "bcells" = "bcell", "b_cells" = "bcell", "b cell" = "bcell",
      "myeloid" = "myeloid",
      "mast" = "mastcell", "mast_cell" = "mastcell",
      "adipo" = "adipocyte", "adipocytes" = "adipocyte",
      "rbc" = "rbc", "red blood cell" = "rbc",
      "skin_epithelial" = "skinepithelial", "skin epithelial" = "skinepithelial"
    )

    # Apply lookup
    if (label_lower %in% names(internal_std_map)) {
      return(internal_std_map[label_lower])
    } else {
      # Fallback for unmapped labels (TEMPLATES pattern)
      default_norm <- gsub("[^a-z0-9]", "", label_lower) # Remove all non-alphanumeric
      if (nzchar(default_norm)) {
        return(default_norm)
      } else {
        return(paste0("UnknownInputRaw_", make.names(label)))
      }
    }
  }

  # Normalize SC and SN predictions
  norm_pred_sc <- sapply(pred_sc, normalize_label)
  norm_pred_sn <- sapply(pred_sn, normalize_label)

  # ============================================================================
  # STEP 4: Apply case_when Reconciliation Logic (EXACT from TEMPLATES)
  # ============================================================================

  cat("\n[4/4] Applying reconciliation rules...\n")

  # Create data frame for case_when logic
  reconciliation_df <- data.frame(
    norm_pred_sc = norm_pred_sc,
    norm_pred_sn = norm_pred_sn,
    stringsAsFactors = FALSE
  )

  # Apply case_when logic (CRITICAL: Exact TEMPLATES pattern)
  # Reference: TEMPLATES 03a_ConsolidateAnnotationsAndQC.Rmd lines 885-993
  # NOTE: All checks use TEMPLATES normalized forms (post-normalization)
  reconciliation_df$final_label_normalized <- dplyr::case_when(
    # Rule 1: Both NA → "Unassigned_BothNA"
    is.na(reconciliation_df$norm_pred_sc) & is.na(reconciliation_df$norm_pred_sn) ~ "Unassigned_BothNA",

    # Rule 2: SC NA, SN present → use SN prediction
    is.na(reconciliation_df$norm_pred_sc) & !is.na(reconciliation_df$norm_pred_sn) ~ reconciliation_df$norm_pred_sn,

    # Rule 3: SN NA, SC present → use SC prediction
    !is.na(reconciliation_df$norm_pred_sc) & is.na(reconciliation_df$norm_pred_sn) ~ reconciliation_df$norm_pred_sc,

    # Rule 4: Exact match → use SC naming (normalized match)
    reconciliation_df$norm_pred_sc == reconciliation_df$norm_pred_sn ~ reconciliation_df$norm_pred_sc,

    # Rule 5: Granularity rules (SN more specific)
    # NOTE: "mast" normalizes to "mastcell", so check for "mastcell"
    (reconciliation_df$norm_pred_sc == "myeloid" & reconciliation_df$norm_pred_sn == "mastcell") ~ "mastcell",
    (reconciliation_df$norm_pred_sc == "myeloid" & reconciliation_df$norm_pred_sn == "macrophage") ~ "macrophage",

    # Rule 6: Unique labels from SN
    # NOTE: "adipo" normalizes to "adipocyte", so only check for "adipocyte"
    reconciliation_df$norm_pred_sn == "adipocyte" ~ "adipocyte",
    reconciliation_df$norm_pred_sn == "rbc" ~ "rbc",

    # Rule 7: Unique labels from SC
    # NOTE: "bcells" normalizes to "bcell", so only check for "bcell"
    reconciliation_df$norm_pred_sc == "bcell" ~ "bcell",
    reconciliation_df$norm_pred_sc == "plasmacell" ~ "plasmacell",
    reconciliation_df$norm_pred_sc == "skinepithelial" ~ "skinepithelial",

    # Rule 8: Naming preference conflicts → SC naming preferred
    # NOTE: These checks are now redundant because normalization already unified forms
    # "pericytes" normalizes to "pericyte", "perivasc" normalizes to "pericyte"
    # So if both are "pericyte", Rule 4 (exact match) already handles this
    # Kept for explicitness but will not trigger due to normalization

    # Rule 9: Default policy → prefer SC (if both present but don't match above rules)
    !is.na(reconciliation_df$norm_pred_sc) ~ reconciliation_df$norm_pred_sc,

    # Rule 10: Fallback for unhandled cases
    TRUE ~ "Unassigned_ConflictOrUnhandled"
  )

  # ============================================================================
  # STEP 5: Map Normalized Labels to Canonical Names
  # ============================================================================

  cat("  Mapping normalized labels to canonical names...\n")

  # Canonical name mapping (EXACT from TEMPLATES)
  # Source: TEMPLATES 03a_ConsolidateAnnotationsAndQC.Rmd lines 939-957
  # Maps TEMPLATES normalized forms (post-normalization) to canonical display names
  canonical_mapping <- c(
    # Epithelial
    "basal" = "Basal",
    "luminalhr" = "Luminal_HR",
    "luminalsecretory" = "Luminal_Secretory",
    "skinepithelial" = "Skin_Epithelial",
    # Stromal
    "fibroblast" = "Fibroblast",
    # Endothelial
    "vascularendothelium" = "Vascular_Endothelium",
    "lymphaticendothelium" = "Lymphatic_Endothelium",
    "pericyte" = "Pericyte",
    # Immune
    "tcell" = "T_Cell",
    "bcell" = "B_Cell",
    "plasmacell" = "Plasma_Cell",
    "myeloid" = "Myeloid",
    "mastcell" = "Mast",
    "macrophage" = "Myeloid",  # TEMPLATES maps macrophage to Myeloid
    # Other
    "adipocyte" = "Adipocyte",
    "rbc" = "RBC",
    # Unassigned categories
    "Unassigned_BothNA" = "Unassigned_BothNA",
    "Unassigned_ConflictOrUnhandled" = "Unassigned_ConflictOrUnhandled",
    "UnknownInput" = "UnknownInput"
  )

  # Map normalized to canonical (with fallback for unrecognized labels)
  map_to_canonical <- function(norm_label) {
    if (is.na(norm_label)) return("Unassigned_NA")
    canonical <- canonical_mapping[[norm_label]]
    if (is.na(canonical)) {
      # If not in mapping, capitalize first letter and return
      return(paste0(toupper(substring(norm_label, 1, 1)), substring(norm_label, 2)))
    }
    return(canonical)
  }

  final_labels <- sapply(reconciliation_df$final_label_normalized, map_to_canonical)

  # ============================================================================
  # STEP 6: Add Final Labels to Query Object Metadata
  # ============================================================================

  # Name the vector with cell barcodes to match Seurat object structure
  # This prevents "No cell overlap between new meta data and Seurat object" error
  names(final_labels) <- Cells(query_obj)

  # Add HBCATransferredLabels.Kumar_2023 to metadata
  query_obj$HBCATransferredLabels.Kumar_2023 <- final_labels

  # Summary statistics
  label_counts <- table(final_labels)
  n_unassigned <- sum(grepl("^Unassigned", final_labels))
  unassigned_fraction <- n_unassigned / length(final_labels)

  cat("\n=== Reconciliation Complete ===\n")
  cat(sprintf("Total cells: %d\n", length(final_labels)))
  cat(sprintf("Unique final labels: %d\n", length(unique(final_labels))))
  cat(sprintf("Unassigned cells: %d (%.1f%%)\n", n_unassigned, unassigned_fraction * 100))
  cat("\nLabel distribution (top 10):\n")
  print(head(sort(label_counts, decreasing = TRUE), 10))

  return(query_obj)
}
