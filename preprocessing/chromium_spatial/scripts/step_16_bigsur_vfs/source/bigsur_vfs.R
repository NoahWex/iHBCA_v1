# =============================================================================
# Step 16 — BigSur Variable Features: core algorithm
# =============================================================================
# Re-runs BigSur variable-feature discovery on a cell set that has already
# been filtered through the FIVE-WAY gate. Purely descriptive: no QC
# thresholds are applied here, the filter is enforced upstream by the
# wrapper. Returns data structures only — the wrapper handles file I/O,
# retention tracking, and sample manifest generation.
#
# Design principles:
#   - Zero-count gene filter BEFORE CreateSeuratObject (memory efficiency)
#   - Standard BigSur call: correlations = FALSE, variable.features = TRUE
#   - Pass-through of the Pearson-residual normalized data layer for
#     downstream normalization comparisons
# =============================================================================

#' Standard BigSur execution with zero-count gene filtering.
#'
#' Filters zero-count genes, constructs a Seurat object with
#' `CreateSeuratObject(min.cells = min_cells)`, runs BigSur, and returns
#' the variable features, the Pearson-residual normalized matrix, and the
#' per-gene feature ranks for consensus feature selection.
#'
#' @param counts_matrix dgCMatrix; genes x cells.
#' @param min_cells Integer; min cells expressing a gene to include.
#' @return list with elements `vfs`, `normalized_matrix`, `feature_ranks`.
run_bigsur_standard <- function(counts_matrix, min_cells = 2) {

    message("\n--- Filtering zero-count genes ---")
    gene_counts <- Matrix::rowSums(counts_matrix)
    genes_keep  <- names(gene_counts)[gene_counts > 0]
    counts_filtered <- counts_matrix[genes_keep, , drop = FALSE]
    n_genes_removed <- nrow(counts_matrix) - nrow(counts_filtered)

    message(sprintf("Removed %d zero-count genes", n_genes_removed))
    message(sprintf("Filtered matrix: %d genes x %d cells",
                    nrow(counts_filtered), ncol(counts_filtered)))

    message(sprintf("\n--- Creating Seurat object (min.cells = %d) ---", min_cells))
    seu <- Seurat::CreateSeuratObject(counts = counts_filtered, min.cells = min_cells)
    message(sprintf("Seurat object: %d cells x %d features",
                    ncol(seu), nrow(seu)))

    rm(counts_filtered); gc(verbose = FALSE)

    message("\n--- Running BigSur ---")
    message("BigSur parameters: correlations=FALSE, variable.features=TRUE, log.file=FALSE")

    seu <- BigSur::BigSur(
        seu,
        correlations      = FALSE,
        variable.features = TRUE,
        log.file          = FALSE
    )
    message("BigSur complete")

    vfs <- Seurat::VariableFeatures(seu)
    message(sprintf("BigSur identified %d variable features", length(vfs)))

    message("\n--- Extracting normalized data layer ---")
    normalized_matrix <- Seurat::GetAssayData(seu, layer = "data", assay = "RNA")
    message(sprintf("Normalized matrix: %d genes x %d cells",
                    nrow(normalized_matrix), ncol(normalized_matrix)))

    message("\n--- Extracting feature ranks ---")
    feature_metadata <- seu@assays$RNA@meta.data
    feature_ranks <- data.frame(
        gene     = rownames(feature_metadata),
        rank     = feature_metadata$var.features.rank,
        variance = feature_metadata$mcfanos,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    feature_ranks <- feature_ranks[!is.na(feature_ranks$rank), ]
    message(sprintf("Extracted ranks for %d variable features", nrow(feature_ranks)))

    rm(seu); gc(verbose = FALSE)

    list(
        vfs               = vfs,
        normalized_matrix = normalized_matrix,
        feature_ranks     = feature_ranks
    )
}


#' Compute VF retention metrics (baseline step 01 vs filtered step 16).
#'
#' Returns a 1-row, 12-column data.frame with overlap / gain / loss counts,
#' retention percentage, and mean expression of each VF category. Expression
#' metrics use the per-sample count matrix and are reported separately for
#' baseline, filtered, retained (intersection), lost, and gained VF sets.
calculate_retention_metrics <- function(vfs_baseline, vfs_filtered, sample_id,
                                        counts_matrix = NULL, final_passing_cells = NULL) {

    compute_mean_expression <- function(gene_set, counts) {
        if (length(gene_set) == 0) return(NA_real_)
        genes_present <- intersect(gene_set, rownames(counts))
        if (length(genes_present) == 0) return(NA_real_)
        gene_means <- Matrix::rowMeans(counts[genes_present, , drop = FALSE])
        mean(gene_means)
    }

    overlap <- intersect(vfs_baseline, vfs_filtered)
    lost    <- setdiff(vfs_baseline, vfs_filtered)
    gained  <- setdiff(vfs_filtered, vfs_baseline)

    if (!is.null(counts_matrix) && !is.null(final_passing_cells)) {
        message("\n--- Calculating expression metrics ---")
        cells_in_matrix <- colnames(counts_matrix)
        if (!all(cells_in_matrix %in% final_passing_cells)) {
            warning("counts_matrix contains cells not in final_passing_cells")
        }

        mean_expr_baseline <- compute_mean_expression(vfs_baseline, counts_matrix)
        mean_expr_filtered <- compute_mean_expression(vfs_filtered, counts_matrix)
        mean_expr_retained <- compute_mean_expression(overlap,      counts_matrix)
        mean_expr_lost     <- compute_mean_expression(lost,         counts_matrix)
        mean_expr_gained   <- compute_mean_expression(gained,       counts_matrix)

        message(sprintf("Mean expression - Baseline VFs: %.3f", mean_expr_baseline))
        message(sprintf("Mean expression - Filtered VFs: %.3f", mean_expr_filtered))
        message(sprintf("Mean expression - Retained VFs: %.3f", mean_expr_retained))
        message(sprintf("Mean expression - Lost VFs: %.3f",     mean_expr_lost))
        message(sprintf("Mean expression - Gained VFs: %.3f",   mean_expr_gained))
    } else {
        mean_expr_baseline <- NA_real_
        mean_expr_filtered <- NA_real_
        mean_expr_retained <- NA_real_
        mean_expr_lost     <- NA_real_
        mean_expr_gained   <- NA_real_
    }

    metrics <- data.frame(
        sample_id          = sample_id,
        n_vfs_baseline     = length(vfs_baseline),
        n_vfs_filtered     = length(vfs_filtered),
        n_vfs_overlap      = length(overlap),
        n_vfs_lost         = length(lost),
        n_vfs_gained       = length(gained),
        retention_pct      = (length(overlap) / length(vfs_baseline)) * 100,
        mean_expr_baseline = mean_expr_baseline,
        mean_expr_filtered = mean_expr_filtered,
        mean_expr_retained = mean_expr_retained,
        mean_expr_lost     = mean_expr_lost,
        mean_expr_gained   = mean_expr_gained,
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    message(sprintf("Baseline VFs: %d", metrics$n_vfs_baseline))
    message(sprintf("Filtered VFs: %d", metrics$n_vfs_filtered))
    message(sprintf("Overlap: %d (%.1f%% retention)",
                    metrics$n_vfs_overlap, metrics$retention_pct))
    message(sprintf("Lost: %d",   metrics$n_vfs_lost))
    message(sprintf("Gained: %d", metrics$n_vfs_gained))

    metrics
}


#' Main step 16 BigSur VF algorithm (no file I/O).
#'
#' Loads the per-sample raw H5, subsets to FIVE-WAY-filtered cells, runs
#' BigSur, and returns the variable features, retention metrics, and the
#' normalized matrix. Cell barcode matching is case-sensitive: the wrapper
#' strips any sample_id prefix from barcodes before calling this function.
run_post_qc_vfs <- function(h5_path, sample_id, vfs_baseline,
                            final_passing_cells, bigsur_min_cells = 2) {

    message(sprintf("\n=== Step 16: BigSur Variable Features for %s ===", sample_id))

    if (!file.exists(h5_path)) stop(sprintf("H5 file not found: %s", h5_path))
    if (length(final_passing_cells) == 0) stop("final_passing_cells is empty")
    if (length(vfs_baseline) == 0) stop("vfs_baseline is empty")

    message("\n--- Step 1: Loading raw H5 count matrix ---")
    counts_raw <- Seurat::Read10X_h5(h5_path)
    message(sprintf("Raw matrix: %d genes x %d cells",
                    nrow(counts_raw), ncol(counts_raw)))

    message(sprintf("\n--- Step 2: Subsetting to %d final passing cells ---",
                    length(final_passing_cells)))
    missing_cells <- setdiff(final_passing_cells, colnames(counts_raw))
    if (length(missing_cells) > 0) {
        warning(sprintf("%d passing cells not found in H5 matrix (will be skipped)",
                        length(missing_cells)))
    }
    cells_to_keep <- intersect(final_passing_cells, colnames(counts_raw))
    if (length(cells_to_keep) == 0) {
        stop("No passing cells found in H5 matrix - check cell barcode matching")
    }

    counts_filtered <- counts_raw[, cells_to_keep, drop = FALSE]
    message(sprintf("Filtered matrix: %d genes x %d cells",
                    nrow(counts_filtered), ncol(counts_filtered)))

    rm(counts_raw); gc(verbose = FALSE)

    message("\n--- Step 3: Running BigSur on filtered cells ---")
    bigsur_results <- run_bigsur_standard(
        counts_matrix = counts_filtered,
        min_cells     = bigsur_min_cells
    )

    vfs_filtered      <- bigsur_results$vfs
    normalized_matrix <- bigsur_results$normalized_matrix

    message("\n--- Step 4: Calculating retention metrics ---")
    retention_metrics <- calculate_retention_metrics(
        vfs_baseline        = vfs_baseline,
        vfs_filtered        = vfs_filtered,
        sample_id           = sample_id,
        counts_matrix       = counts_filtered,
        final_passing_cells = cells_to_keep
    )

    rm(counts_filtered); gc(verbose = FALSE)

    message(sprintf("\n=== Step 16 complete === %s | %d VFs | retention %.1f%%",
                    sample_id, length(vfs_filtered),
                    retention_metrics$retention_pct))

    list(
        vfs_filtered      = vfs_filtered,
        retention_metrics = retention_metrics,
        normalized_matrix = normalized_matrix,
        feature_ranks     = bigsur_results$feature_ranks
    )
}
