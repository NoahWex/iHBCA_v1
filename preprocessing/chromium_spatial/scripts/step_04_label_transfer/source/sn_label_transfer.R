# ============================================================================
# SN Label Transfer Algorithm
# ============================================================================
#
# Purpose: Transfer cell type labels from SN reference to query object
#
# This file contains pure algorithm functions for SN reference label transfer.
# No file I/O - all paths and data passed as function arguments.
#
# Dependencies: Seurat
# ============================================================================

#' Perform SN Label Transfer
#'
#' Transfers cell type labels from the SN (single-nucleus) reference dataset
#' to a query Seurat object using Seurat's FindTransferAnchors and TransferData
#' workflow.
#'
#' @param query_obj Seurat object (must be processed: normalized, scaled, PCA, UMAP)
#' @param sn_ref_path Character. Path to SN reference RDS file
#' @param config List with parameters:
#'   - reference_assay: Character. Assay name in reference (default: "integrated")
#'   - query_assay: Character. Assay name in query (default: "RNA")
#'   - label_column: Character. Metadata column in reference with cell types (default: "nuc_celltype")
#'   - transfer_dims: Integer. Number of PCA dimensions for transfer (default: 30)
#'
#' @return Modified query_obj with added metadata columns:
#'   - predicted.celltype.SN: Character vector of SN predictions
#'   - prediction.score.SN: Numeric vector of prediction scores (0-1)
#'   Returns NULL on error.
#'
#' @export
#'
#' @examples
#' # config <- list(
#' #   reference_assay = "integrated",
#' #   query_assay = "RNA",
#' #   label_column = "nuc_celltype",
#' #   transfer_dims = 30
#' # )
#' # query_obj <- perform_sn_label_transfer(query_obj, "path/to/sn_ref.rds", config)
perform_sn_label_transfer <- function(query_obj, sn_ref_path, config) {

  # ========================================
  # Input Validation
  # ========================================
  if (!inherits(query_obj, "Seurat")) {
    stop("query_obj must be a Seurat object")
  }

  if (!file.exists(sn_ref_path)) {
    stop(sprintf("SN reference file not found: %s", sn_ref_path))
  }

  # Set defaults (note: SN uses "nuc_celltype" by default)
  reference_assay <- config$reference_assay %||% "integrated"
  query_assay <- config$query_assay %||% "RNA"
  label_column <- config$label_column %||% "nuc_celltype"
  transfer_dims <- config$transfer_dims %||% 30

  message(sprintf("Starting SN label transfer with %d dimensions", transfer_dims))

  # Memory diagnostics helper
  report_memory <- function(label) {
    mem_info <- gc()
    mem_used_mb <- sum(mem_info[, 2])  # Ncells + Vcells in MB
    message(sprintf("  [MEM] %s: %.1f MB used", label, mem_used_mb))
  }

  # ========================================
  # Load SN Reference
  # ========================================
  message(sprintf("Loading SN reference from: %s", sn_ref_path))
  report_memory("Before loading SN reference")

  sn_ref <- tryCatch({
    readRDS(sn_ref_path)
  }, error = function(e) {
    stop(sprintf(
      "Failed to load SN reference RDS file\n  Path: %s\n  Error: %s\n  Hint: Verify file exists and is readable",
      sn_ref_path, e$message
    ))
  })

  # Validate reference structure
  if (!inherits(sn_ref, "Seurat")) {
    stop(sprintf(
      "SN reference is not a Seurat object\n  Path: %s\n  Class: %s\n  Expected: Seurat",
      sn_ref_path, class(sn_ref)[1]
    ))
  }

  if (!reference_assay %in% Assays(sn_ref)) {
    stop(sprintf(
      "SN reference missing required assay\n  Requested: '%s'\n  Available: %s\n  Hint: Check module_configs.yaml sn_reference_assay parameter",
      reference_assay, paste(Assays(sn_ref), collapse = ", ")
    ))
  }

  if (!label_column %in% colnames(sn_ref@meta.data)) {
    stop(sprintf(
      "SN reference missing label column\n  Requested: '%s'\n  Available: %s\n  Hint: Check module_configs.yaml sn_label_column parameter",
      label_column, paste(head(colnames(sn_ref@meta.data), 10), collapse = ", ")
    ))
  }

  message(sprintf("SN reference loaded: %d cells, %d features",
                  ncol(sn_ref), nrow(sn_ref)))

  # Report object sizes
  ref_size_mb <- as.numeric(object.size(sn_ref)) / 1024^2
  message(sprintf("  SN reference object size: %.1f MB", ref_size_mb))
  report_memory("After loading SN reference")

  # Clean up memory before processing
  gc(verbose = FALSE)
  report_memory("After gc() post-load")

  # ========================================
  # Update to Seurat v5 format
  # ========================================
  # UpdateSeuratObject handles deprecated slots (e.g., images) internally
  # Critical for loading legacy RDS files from Seurat v3/v4
  message("  Updating SN reference to Seurat v5 format...")
  update_start <- Sys.time()
  sn_ref <- UpdateSeuratObject(sn_ref)
  update_duration <- round(difftime(Sys.time(), update_start, units = "secs"), 1)
  message(sprintf("  SN reference updated successfully (%.1fs)", update_duration))
  report_memory("After UpdateSeuratObject")

  # ========================================
  # Set Default Assays
  # ========================================
  DefaultAssay(sn_ref) <- reference_assay
  DefaultAssay(query_obj) <- query_assay

  message(sprintf("Using reference assay: %s, query assay: %s",
                  reference_assay, query_assay))

  # Clean up memory before anchor finding (CRITICAL - FindTransferAnchors is memory-intensive)
  gc(verbose = FALSE)
  report_memory("Before FindTransferAnchors (CRITICAL CHECKPOINT)")

  query_size_mb <- as.numeric(object.size(query_obj)) / 1024^2
  message(sprintf("  Query object size: %.1f MB", query_size_mb))
  message(sprintf("  Total objects in memory before anchors: %.1f MB", ref_size_mb + query_size_mb))

  # ========================================
  # Find Transfer Anchors
  # ========================================
  message("Finding transfer anchors between query and SN reference...")
  anchor_start <- Sys.time()

  transfer_anchors_SN <- tryCatch({
    FindTransferAnchors(
      reference = sn_ref,
      query = query_obj,
      normalization.method = "LogNormalize",
      reference.assay = reference_assay,
      query.assay = query_assay,
      dims = 1:transfer_dims,
      verbose = FALSE
    )
  }, error = function(e) {
    stop(sprintf(
      "Failed to find SN transfer anchors\n  Query cells: %d\n  Reference cells: %d\n  Dims: 1:%d\n  Error: %s\n  Hint: Query and reference may have incompatible features or insufficient overlap",
      ncol(query_obj), ncol(sn_ref), transfer_dims, e$message
    ))
  })

  anchor_duration <- round(difftime(Sys.time(), anchor_start, units = "secs"), 1)
  n_anchors <- nrow(transfer_anchors_SN@anchors)
  anchors_size_mb <- as.numeric(object.size(transfer_anchors_SN)) / 1024^2
  message(sprintf("Found %d transfer anchors (%.1fs, %.1f MB)", n_anchors, anchor_duration, anchors_size_mb))
  report_memory("After FindTransferAnchors")

  # Warn if anchor count is low
  if (n_anchors < 100) {
    warning(sprintf(
      "LOW ANCHOR COUNT (%d) - SN predictions may be unreliable\n  Query cells: %d\n  Reference cells: %d",
      n_anchors, ncol(query_obj), ncol(sn_ref)
    ))
  }

  # Clean up memory after anchor finding
  gc(verbose = FALSE)
  report_memory("After gc() post-anchors")

  # ========================================
  # Transfer Data
  # ========================================
  message(sprintf("Transferring cell type labels from column: %s", label_column))
  transfer_start <- Sys.time()

  # Extract reference labels
  refdata_celltype_SN <- sn_ref[[label_column, drop = TRUE]]

  predictions_SN <- tryCatch({
    TransferData(
      anchorset = transfer_anchors_SN,
      refdata = refdata_celltype_SN,
      dims = 1:transfer_dims,
      verbose = FALSE
    )
  }, error = function(e) {
    stop(sprintf(
      "Failed to transfer SN labels\n  Anchors: %d\n  Dims: 1:%d\n  Error: %s",
      nrow(transfer_anchors_SN@anchors), transfer_dims, e$message
    ))
  })

  transfer_duration <- round(difftime(Sys.time(), transfer_start, units = "secs"), 1)
  message(sprintf("SN label transfer complete (%.1fs)", transfer_duration))
  report_memory("After TransferData")

  # ========================================
  # Add Predictions to Query Object
  # ========================================
  query_obj$predicted.celltype.SN <- predictions_SN$predicted.id
  query_obj$prediction.score.SN <- predictions_SN$prediction.score.max

  message(sprintf("Added SN predictions to query object metadata"))
  message(sprintf("  - Unique SN labels transferred: %d",
                  length(unique(query_obj$predicted.celltype.SN))))
  message(sprintf("  - Mean prediction score: %.3f",
                  mean(query_obj$prediction.score.SN, na.rm = TRUE)))
  message(sprintf("  - Median prediction score: %.3f",
                  median(query_obj$prediction.score.SN, na.rm = TRUE)))

  # Show label distribution
  label_counts <- table(query_obj$predicted.celltype.SN)
  message("\nSN Label Distribution:")
  print(label_counts)

  # Clean up SN reference and intermediate objects to free memory
  message("Cleaning up SN reference and intermediate objects...")
  rm(sn_ref, transfer_anchors_SN, predictions_SN)
  gc(verbose = FALSE)
  report_memory("After final cleanup")

  message("SN transfer memory cleaned up")

  # Return modified query object
  return(query_obj)
}


#' %||% Operator (NULL coalescing)
#'
#' Returns left side if not NULL, otherwise returns right side
#' @keywords internal
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}
