# ==============================================================================
# Step 10a: Add BigSur Normalization to Integrated Seurat Object
# ==============================================================================
# Purpose: Aggregate Step 05 per-sample BigSur normalized matrices and add as
#          new assay to Step 09 integrated Seurat object
#
# Input:
#   - integrated_seurat_path: Step 09 integrated Seurat object (RDS)
#   - bigsur_matrix_paths: Named list of Step 05 BigSur matrix paths (sample_id => path)
#   - output_path: Where to save updated object
#
# Output:
#   - Updated Seurat object with BigSur assay added
#   - Summary statistics
#
# Design:
#   - Creates Seurat objects from Step 05 BigSur matrices
#   - Uses Seurat's native merge function to combine samples
#   - Adds merged BigSur assay to Step 09 integrated object
#   - Maintains Seurat compatibility throughout
#
# Author: Noah Wechter
# Date: 2025-11-18
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

#' Compute Consensus Top Features Using Rank Aggregation
#'
#' @description
#' Identifies consensus variable features across all samples using rank
#' aggregation. Uses mean rank (not median) to capture artifact markers
#' that appear strongly in a subset of samples.
#'
#' @param feature_ranks_paths Named list of paths to feature_ranks.rds files
#'   (one per sample, from Step 05)
#' @param n_consensus Integer. Number of top consensus features to select
#' @param verbose Logical. Print progress messages
#'
#' @return Character vector of consensus feature gene names (length = n_consensus)
compute_consensus_features <- function(feature_ranks_paths, n_consensus = 2000,
                                      verbose = TRUE) {

  if (verbose) {
    cat("\n============================================================\n")
    cat(sprintf("Computing Consensus Top %d Features\n", n_consensus))
    cat("============================================================\n")
  }

  # ============================================================================
  # Load all feature ranks
  # ============================================================================

  all_ranks <- list()
  max_ranks <- numeric()

  for (sample_id in names(feature_ranks_paths)) {
    ranks_file <- feature_ranks_paths[[sample_id]]
    if (!file.exists(ranks_file)) {
      warning(sprintf("Feature ranks not found for %s: %s", sample_id, ranks_file))
      next
    }

    ranks_df <- readRDS(ranks_file)
    all_ranks[[sample_id]] <- setNames(ranks_df$rank, ranks_df$gene)
    max_ranks[sample_id] <- max(ranks_df$rank, na.rm = TRUE)
  }

  if (verbose) {
    cat(sprintf("Loaded feature ranks from %d samples\n", length(all_ranks)))
  }

  # ============================================================================
  # Collect all unique genes across all samples
  # ============================================================================

  all_genes <- unique(unlist(lapply(all_ranks, names)))

  if (verbose) {
    cat(sprintf("Total unique genes across all samples: %d\n", length(all_genes)))
  }

  # ============================================================================
  # Compute mean rank for each gene (with penalty for missing)
  # ============================================================================

  gene_avg_ranks <- sapply(all_genes, function(gene) {
    ranks_for_gene <- numeric(length(all_ranks))

    for (i in seq_along(all_ranks)) {
      sample_id <- names(all_ranks)[i]
      if (gene %in% names(all_ranks[[sample_id]])) {
        # Gene is VF in this sample - use its rank
        ranks_for_gene[i] <- all_ranks[[sample_id]][gene]
      } else {
        # Gene is NOT VF in this sample - assign penalty (max_rank + 1)
        ranks_for_gene[i] <- max_ranks[sample_id] + 1
      }
    }

    mean(ranks_for_gene)
  })

  # ============================================================================
  # Select top N consensus features by mean rank
  # ============================================================================

  consensus_genes <- names(sort(gene_avg_ranks))[1:n_consensus]

  if (verbose) {
    cat(sprintf("\nSelected top %d consensus features by mean rank\n", n_consensus))
    cat(sprintf("  Best consensus gene: %s (mean rank: %.1f)\n",
                consensus_genes[1], gene_avg_ranks[consensus_genes[1]]))
    cat(sprintf("  Worst consensus gene: %s (mean rank: %.1f)\n",
                consensus_genes[n_consensus], gene_avg_ranks[consensus_genes[n_consensus]]))
    cat("============================================================\n")
  }

  return(consensus_genes)
}

#' Add BigSur Normalization to Integrated Seurat Object
#'
#' @param integrated_seurat_path Path to Step 09 integrated Seurat object (RDS)
#' @param bigsur_matrix_paths Named list: sample_id => Step 05 BigSur matrix path
#' @param feature_ranks_paths Named list: sample_id => Step 05 feature_ranks.rds path
#' @param output_path Path to save updated Seurat object
#' @param n_consensus_features Number of consensus features to select via rank aggregation (default: 2000)
#' @param verbose Print progress messages
#'
#' @return List with summary_stats and updated_seurat
add_bigsur_to_integrated_object <- function(
  integrated_seurat_path,
  bigsur_matrix_paths,
  feature_ranks_paths,
  output_path,
  n_consensus_features = 2000,
  verbose = TRUE
) {

  # ============================================================================
  # Step 1: Load Integrated Object
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 1: Loading Step 09 integrated Seurat object\n")
    cat("============================================================\n")
  }

  if (!file.exists(integrated_seurat_path)) {
    stop("Integrated Seurat object not found: ", integrated_seurat_path)
  }

  seu_integrated <- readRDS(integrated_seurat_path)

  # Extract metadata
  n_cells_total <- ncol(seu_integrated)
  current_assays <- Assays(seu_integrated)

  if (verbose) {
    cat("Loaded integrated object:\n")
    cat("  - Cells:", n_cells_total, "\n")
    cat("  - Genes (RNA):", nrow(seu_integrated[["RNA"]]), "\n")
    cat("  - Assays:", paste(current_assays, collapse = ", "), "\n")
  }

  # Extract cell order for later matching
  cell_order <- colnames(seu_integrated)

  # ============================================================================
  # Step 2: Compute Consensus Features via Rank Aggregation
  # ============================================================================

  consensus_vfs <- compute_consensus_features(
    feature_ranks_paths = feature_ranks_paths,
    n_consensus = n_consensus_features,
    verbose = verbose
  )

  # ============================================================================
  # Step 3: Load BigSur Matrices with Vertical Slicing
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 3: Loading BigSur matrices (vertically sliced to consensus features)\n")
    cat("============================================================\n")
  }

  # Store subsetted data matrices (not Seurat objects)
  data_matrices <- list()

  for (sample_id in names(bigsur_matrix_paths)) {

    if (verbose) cat("  Processing:", sample_id, "... ")

    matrix_path <- bigsur_matrix_paths[[sample_id]]

    if (!file.exists(matrix_path)) {
      warning("BigSur matrix not found for ", sample_id, ": ", matrix_path)
      next
    }

    # Load BigSur normalized matrix (sparse dgCMatrix from Step 05)
    bigsur_matrix <- readRDS(matrix_path)

    # Add sample prefix to cell barcodes to match Step 09
    # Step 05 matrices have unprefixed barcodes
    # Step 09 uses prefixed: {sample_id}_{barcode}
    cell_barcodes_prefixed <- paste0(sample_id, "_", colnames(bigsur_matrix))
    colnames(bigsur_matrix) <- cell_barcodes_prefixed

    # VERTICAL SLICE: Subset to consensus features IMMEDIATELY
    # This prevents loading unwanted genes and avoids dense matrix memory explosion
    genes_in_sample <- intersect(consensus_vfs, rownames(bigsur_matrix))
    bigsur_matrix_subset <- bigsur_matrix[genes_in_sample, , drop = FALSE]

    # Store subsetted matrix
    data_matrices[[sample_id]] <- bigsur_matrix_subset

    if (verbose) {
      cat(sprintf("%d consensus genes × %d cells\n",
                  length(genes_in_sample), ncol(bigsur_matrix_subset)))
    }

    # Memory cleanup
    rm(bigsur_matrix, bigsur_matrix_subset, cell_barcodes_prefixed)
    gc(verbose = FALSE)
  }

  if (length(data_matrices) == 0) {
    stop("No BigSur matrices loaded successfully")
  }

  if (verbose) {
    cat(sprintf("Loaded %d samples with consensus feature slicing\n", length(data_matrices)))
  }

  # ============================================================================
  # Step 4: Create Seurat Objects with BigSur Data (Following Step 09 Pattern)
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 4: Creating Seurat objects with BigSur data\n")
    cat("============================================================\n")
  }

  # Create list of Seurat objects (one per sample)
  # Pattern from Step 09: integration_preview.Rmd:519-527
  seurat_list <- list()

  for (sample_id in names(data_matrices)) {
    mat <- data_matrices[[sample_id]]

    if (verbose) {
      cat(sprintf("  Creating Seurat object for %s: %d genes × %d cells\n",
                  sample_id, nrow(mat), ncol(mat)))
    }

    # Create Seurat object with BigSur data as normalized data
    # Use counts = mat as placeholder (Seurat requires counts slot)
    # The "data" slot will be set to the same normalized matrix
    seurat_obj <- CreateSeuratObject(
      counts = mat,
      project = sample_id,
      min.cells = 0,
      min.features = 0
    )

    # Store in list
    seurat_list[[sample_id]] <- seurat_obj
  }

  # Clean up data matrices
  rm(data_matrices)
  gc(verbose = FALSE)

  if (verbose) {
    cat(sprintf("\n  Created %d Seurat objects\n", length(seurat_list)))
  }

  # ============================================================================
  # Step 5: Merge Seurat Objects (Following Step 09 Pattern)
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 5: Merging Seurat objects\n")
    cat("============================================================\n")
  }

  # Pattern from Step 09: integration_preview.Rmd:544-578
  if (length(seurat_list) == 0) {
    stop("No Seurat objects to merge!")
  }

  # Get first object
  first_obj <- seurat_list[[1]]

  # Merge remaining objects
  if (length(seurat_list) > 1) {
    remaining_objs <- seurat_list[2:length(seurat_list)]

    seurat_bigsur <- merge(
      x = first_obj,
      y = remaining_objs,
      add.cell.ids = NULL,  # Cell IDs already prefixed from Step 05
      project = "bigsur_integrated"
    )
  } else {
    seurat_bigsur <- first_obj
  }

  # Clean up list
  rm(seurat_list, first_obj)
  if (exists("remaining_objs")) rm(remaining_objs)
  gc(verbose = FALSE)

  if (verbose) {
    cat(sprintf("  Merged: %d cells × %d genes\n",
                ncol(seurat_bigsur), nrow(seurat_bigsur)))
  }

  # ============================================================================
  # Step 6: JoinLayers (Following Step 09 Pattern)
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 6: Joining layers\n")
    cat("============================================================\n")
  }

  # Pattern from Step 09: integration_preview.Rmd:735-748
  # Memory cleanup before JoinLayers
  cat("  Running gc() before JoinLayers...\n")
  gc()

  # Join layers (CRITICAL for Seurat v5)
  cat("  Executing JoinLayers...\n")
  seurat_bigsur <- JoinLayers(seurat_bigsur)
  cat("  JoinLayers complete\n")

  # Memory cleanup after JoinLayers
  cat("  Running gc() after JoinLayers...\n")
  gc()

  if (verbose) {
    cat(sprintf("  Final BigSur object: %d cells × %d genes\n",
                ncol(seurat_bigsur), nrow(seurat_bigsur)))
  }

  # ============================================================================
  # Step 7: Extract BigSur Assay and Add to Integrated Object
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 7: Adding BigSur assay to integrated object\n")
    cat("============================================================\n")
  }

  # Extract the RNA assay (contains BigSur data) and rename it
  bigsur_assay <- seurat_bigsur[["RNA"]]

  # Update key
  Key(bigsur_assay) <- "bigsur_"

  # Clean up merged object
  rm(seurat_bigsur)
  gc(verbose = FALSE)

  # Add to integrated object
  seu_integrated[["BigSur"]] <- bigsur_assay

  # Clean up assay
  rm(bigsur_assay)
  gc(verbose = FALSE)

  # Verify structure
  stopifnot("BigSur" %in% Assays(seu_integrated))
  stopifnot("RNA" %in% Assays(seu_integrated))

  # Keep RNA as default assay
  DefaultAssay(seu_integrated) <- "RNA"

  if (verbose) {
    cat("  Added BigSur assay to integrated object\n")
    cat("  - Assays:", paste(Assays(seu_integrated), collapse = ", "), "\n")
    cat("  - Default assay:", DefaultAssay(seu_integrated), "\n")
    cat("  - BigSur genes:", nrow(seu_integrated[["BigSur"]]), "\n")
    cat("  - BigSur cells:", ncol(seu_integrated[["BigSur"]]), "\n")
  }

  # ============================================================================
  # Step 8: Save Updated Object & Return Summary
  # ============================================================================

  if (verbose) {
    cat("\n============================================================\n")
    cat("Step 8: Saving updated Seurat object\n")
    cat("============================================================\n")
  }

  # Save RDS
  saveRDS(seu_integrated, output_path)

  if (verbose) {
    cat("Saved to:", output_path, "\n")
  }

  # Compute summary statistics
  summary_stats <- list(
    n_cells = ncol(seu_integrated),
    n_genes_bigsur = nrow(seu_integrated[["BigSur"]]),
    n_genes_rna = nrow(seu_integrated[["RNA"]]),
    assays = Assays(seu_integrated),
    object_size_mb = as.numeric(object.size(seu_integrated)) / 1e6,
    n_consensus_features = n_consensus_features,
    n_samples_loaded = length(bigsur_matrix_paths)
  )

  if (verbose) {
    cat("\nSummary:\n")
    cat("  - Cells:", summary_stats$n_cells, "\n")
    cat("  - BigSur genes:", summary_stats$n_genes_bigsur, "\n")
    cat("  - RNA genes:", summary_stats$n_genes_rna, "\n")
    cat("  - Object size:", round(summary_stats$object_size_mb, 1), "MB\n")
    cat("  - Consensus features:", summary_stats$n_consensus_features, "\n")
    cat("  - Samples loaded:", summary_stats$n_samples_loaded, "\n")
  }

  # Return
  return(list(
    summary_stats = summary_stats,
    updated_seurat = seu_integrated
  ))
}
