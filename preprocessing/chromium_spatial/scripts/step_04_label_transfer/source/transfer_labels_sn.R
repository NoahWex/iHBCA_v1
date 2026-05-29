# Step 07 Label Transfer - SN Reference Transfer
#
# Part of: Spatial HBCA Preprocessing Pipeline
# Step: 07_LabelTransfer
# Type: Algorithm (Pure R function, loads reference from RDS)
#
# Purpose: Perform Seurat label transfer using SN (single-nucleus) reference from
#          Kumar 2023 HBCA atlas
#
# Dependencies:
#   - Seurat (>= 5.0) **CRITICAL: Seurat v5**
#
# Author: Agent 3 (Algorithm Developer)
# Date: 2025-11-03
# Pattern: TRACE principles (wrapper + source separation)
# Seurat Version: v5

#' Transfer Cell Type Labels from SN Reference
#'
#' Performs Seurat label transfer using SN (single-nucleus) reference from
#' Kumar 2023 HBCA atlas.
#'
#' @param query_obj Seurat object. Query object (with SC predictions already added)
#' @param sn_ref_path Character. Path to SN reference RDS file
#' @param config List. Transfer parameters from module_configs.yaml
#'   - transfer_dims: PCA dimensions for FindTransferAnchors (default: 30)
#'   - sn_reference_assay: SN reference assay name (default: "integrated")
#'   - sn_label_column: SN cell type column (default: "nuc_celltype")
#'   - query_assay: Query assay name (default: "RNA")
#'   - random_seed: Random seed (default: 42)
#'
#' @return Seurat object with added metadata columns:
#'   - predicted.celltype.SN: SN reference cell type predictions
#'   - prediction.score.SN: SN prediction confidence scores (0-1)
#'
#' @details
#' Uses Seurat FindTransferAnchors + TransferData workflow:
#' 1. Load SN reference from RDS
#' 2. Find transfer anchors between query and SN reference
#' 3. Transfer cell type labels with prediction scores
#' 4. Add predictions to query object metadata
#'
#' SN reference captures:
#' - Adipocytes (nuclear-enriched, often lost in SC)
#' - Finer myeloid resolution (Mast vs Myeloid)
#'
#' @examples
#' query_obj <- transfer_labels_sn(
#'   query_obj = query_obj,
#'   sn_ref_path = "/path/to/sn_reference.rds",
#'   config = list(
#'     transfer_dims = 30,
#'     sn_reference_assay = "integrated",
#'     sn_label_column = "nuc_celltype",
#'     query_assay = "RNA",
#'     random_seed = 42
#'   )
#' )
transfer_labels_sn <- function(
  query_obj,
  sn_ref_path,
  config = list()
) {

  # Extract config parameters with defaults
  transfer_dims <- config$transfer_dims %||% 30
  sn_reference_assay <- config$sn_reference_assay %||% "integrated"
  sn_label_column <- config$sn_label_column %||% "nuc_celltype"
  query_assay <- config$query_assay %||% "RNA"
  random_seed <- config$random_seed %||% 42

  # Set random seed for reproducibility
  set.seed(random_seed)

  cat("\n=== SN Reference Label Transfer ===\n")

  # ============================================================================
  # STEP 1: Load SN Reference from RDS
  # ============================================================================

  cat("\n[1/4] Loading SN reference...\n")

  # Defensive check
  if (!file.exists(sn_ref_path)) {
    stop(sprintf("SN reference RDS not found: %s", sn_ref_path))
  }

  # Load reference
  sn_ref <- readRDS(sn_ref_path)

  # Verify it's a Seurat object
  if (!inherits(sn_ref, "Seurat")) {
    stop(sprintf("SN reference is not a Seurat object: %s", class(sn_ref)))
  }

  # Update to Seurat v5 format if needed (handles v3/v4 objects)
  # Critical for loading legacy RDS files that contain deprecated slots (e.g., images)
  cat("  Updating SN reference to Seurat v5 format...\n")
  sn_ref <- Seurat::UpdateSeuratObject(sn_ref)

  cat(sprintf("  Loaded SN reference: %d cells, %d genes\n",
              ncol(sn_ref), nrow(sn_ref)))

  # ============================================================================
  # STEP 2: Set Default Assay and Verify Label Column
  # ============================================================================

  cat("\n[2/4] Setting reference assay and verifying label column...\n")

  # Set default assay
  if (!(sn_reference_assay %in% Seurat::Assays(sn_ref))) {
    stop(sprintf("SN reference assay '%s' not found in reference. Available: %s",
                 sn_reference_assay, paste(Seurat::Assays(sn_ref), collapse=", ")))
  }
  Seurat::DefaultAssay(sn_ref) <- sn_reference_assay
  cat(sprintf("  Set default assay: %s\n", sn_reference_assay))

  # Verify label column exists
  if (!(sn_label_column %in% colnames(sn_ref@meta.data))) {
    stop(sprintf("SN label column '%s' not found in metadata. Available: %s",
                 sn_label_column, paste(colnames(sn_ref@meta.data), collapse=", ")))
  }
  cat(sprintf("  Verified label column: %s\n", sn_label_column))
  cat(sprintf("  Unique cell types in reference: %d\n", length(unique(sn_ref@meta.data[[sn_label_column]]))))

  # ============================================================================
  # STEP 3: Find Transfer Anchors
  # ============================================================================

  cat("\n[3/4] Finding transfer anchors...\n")

  # Set query assay
  Seurat::DefaultAssay(query_obj) <- query_assay

  # Find anchors using Seurat v5 syntax
  # Use reference.reduction = "pca" for v5 compatibility
  transfer_anchors <- Seurat::FindTransferAnchors(
    reference = sn_ref,
    query = query_obj,
    dims = 1:transfer_dims,
    reference.reduction = "pca",
    verbose = TRUE
  )

  cat(sprintf("  Found %d transfer anchors\n", nrow(transfer_anchors@anchors)))

  # ============================================================================
  # STEP 4: Transfer Labels
  # ============================================================================

  cat("\n[4/4] Transferring labels...\n")

  # Transfer data using Seurat v5 syntax
  predictions <- Seurat::TransferData(
    anchorset = transfer_anchors,
    refdata = sn_ref@meta.data[[sn_label_column]],
    dims = 1:transfer_dims,
    verbose = TRUE
  )

  cat(sprintf("  Transferred predictions: %d cells\n", nrow(predictions)))

  # Add predictions to query object metadata
  # TransferData returns a data frame with predicted.id and prediction.score.max columns
  query_obj$predicted.celltype.SN <- predictions$predicted.id
  query_obj$prediction.score.SN <- predictions$prediction.score.max

  # Summary statistics
  n_assigned <- sum(!is.na(query_obj$predicted.celltype.SN))
  success_rate <- n_assigned / ncol(query_obj)
  mean_score <- mean(query_obj$prediction.score.SN, na.rm = TRUE)

  cat("\n=== SN Transfer Complete ===\n")
  cat(sprintf("Cells with predictions: %d / %d (%.1f%%)\n",
              n_assigned, ncol(query_obj), success_rate * 100))
  cat(sprintf("Mean prediction score: %.3f\n", mean_score))
  cat(sprintf("Unique predicted cell types: %d\n", length(unique(query_obj$predicted.celltype.SN))))

  return(query_obj)
}
