# =============================================================================
# Step 15 — Milo Construction (source algorithm, no I/O)
# =============================================================================
# Pure-algorithm module for building a Milo object from a Seurat object
# carrying scVI embeddings. Called by 1_build_milo_object.R (frozen in the
# canonical migration; kept for reproducibility).
#
# Pattern source: hbca_analysis/analyses/CellFiltering/
#                 milor_contamination_detection_20250808/ (lines 255-337)
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)
    library(SingleCellExperiment)
    library(SummarizedExperiment)
    library(S4Vectors)
    library(miloR)
    library(Matrix)
})

#' Create Milo Object from Seurat
#'
#' Builds neighborhoods on scVI embeddings and counts cells per neighborhood
#' per sample. Returns a Milo object with `calcNhoodDistance` and
#' `buildNhoodGraph` already applied.
create_milo_from_seurat <- function(
    seurat_obj,
    k = 30,
    d = 50,
    prop = 0.1,
    sample_col = "sample_id",
    reduced_dim = "scvi",
    refined = TRUE
) {
    cat("============================================================\n")
    cat("Creating Milo Object from Seurat\n")
    cat("============================================================\n")
    cat("Parameters:\n")
    cat("  k (neighbors):    ", k, "\n")
    cat("  d (dimensions):   ", d, "\n")
    cat("  prop (sampling):  ", prop, "\n")
    cat("  sample_col:       ", sample_col, "\n")
    cat("  reduced_dim:      ", reduced_dim, "\n")
    cat("  refined:          ", refined, "\n")
    cat("============================================================\n\n")

    cat("Step 1/6: Validating Seurat embeddings...\n")
    validate_seurat_embeddings(seurat_obj, reduced_dim, sample_col)

    cat("Step 2/6: Converting Seurat to SingleCellExperiment...\n")
    sce <- seurat_to_sce(seurat_obj, reduced_dim)
    cat("  ", ncol(sce), "cells\n\n")

    cat("Step 3-6: Building Milo neighborhoods...\n")
    milo_obj <- build_neighborhoods(
        sce = sce,
        k = k,
        d = d,
        prop = prop,
        sample_col = sample_col,
        reduced_dim = toupper(reduced_dim),
        refined = refined
    )

    cat("\n============================================================\n")
    cat("Milo Object Construction Complete\n")
    cat("============================================================\n")
    cat("Summary:\n")
    cat("  Cells:         ", ncol(milo_obj), "\n")
    cat("  Neighborhoods: ", ncol(nhoods(milo_obj)), "\n")
    cat("  Samples:       ", ncol(nhoodCounts(milo_obj)), "\n")
    cat("============================================================\n")

    return(milo_obj)
}

validate_seurat_embeddings <- function(
    seurat_obj,
    reduced_dim = "scvi",
    sample_col = "sample_id"
) {
    if (!inherits(seurat_obj, "Seurat")) {
        stop("ERROR: Input must be a Seurat object")
    }

    available_reductions <- names(seurat_obj@reductions)
    matching_reduction <- available_reductions[tolower(available_reductions) == tolower(reduced_dim)]

    if (length(matching_reduction) == 0) {
        stop("ERROR: Reduction '", reduced_dim, "' not found.\n",
             "  Available reductions: ", paste(available_reductions, collapse = ", "))
    }

    actual_reduction_name <- matching_reduction[1]
    cat("  Found reduction: '", actual_reduction_name, "'\n", sep = "")

    embedding_dims <- ncol(Embeddings(seurat_obj, actual_reduction_name))
    cat("  Embedding dimensions: ", embedding_dims, "\n", sep = "")

    umap_variants <- c("umap", "UMAP", "Umap")
    has_umap <- any(umap_variants %in% available_reductions)
    if (!has_umap) {
        warning("UMAP reduction not found — visualization may be limited")
    } else {
        cat("  UMAP found: yes\n")
    }

    required_cols <- c(sample_col, "patient_id")
    missing <- setdiff(required_cols, colnames(seurat_obj@meta.data))
    if (length(missing) > 0) {
        stop("ERROR: Missing required metadata columns: ", paste(missing, collapse = ", "))
    }
    cat("  Required metadata present: ", paste(required_cols, collapse = ", "), "\n", sep = "")

    n_na_sample  <- sum(is.na(seurat_obj@meta.data[[sample_col]]))
    n_na_patient <- sum(is.na(seurat_obj@meta.data[["patient_id"]]))
    if (n_na_sample > 0)  warning("Found ", n_na_sample, " NA values in '", sample_col, "' column")
    if (n_na_patient > 0) warning("Found ", n_na_patient, " NA values in 'patient_id' column")

    cat("  Total cells: ", ncol(seurat_obj), "\n", sep = "")
    return(TRUE)
}

seurat_to_sce <- function(seurat_obj, reduced_dim = "scvi") {
    available_reductions <- names(seurat_obj@reductions)
    actual_reduction <- available_reductions[tolower(available_reductions) == tolower(reduced_dim)][1]

    scvi_matrix <- Embeddings(seurat_obj, actual_reduction)
    cat("  scVI matrix: ", nrow(scvi_matrix), " cells x ", ncol(scvi_matrix), " dimensions\n", sep = "")

    umap_variants <- c("umap", "UMAP", "Umap")
    umap_name <- umap_variants[umap_variants %in% available_reductions][1]

    if (!is.na(umap_name)) {
        umap_matrix <- Embeddings(seurat_obj, umap_name)
        umap_matrix <- umap_matrix[rownames(scvi_matrix), , drop = FALSE]
        cat("  UMAP matrix: ", nrow(umap_matrix), " cells x ", ncol(umap_matrix), " dimensions\n", sep = "")
    } else {
        umap_matrix <- NULL
        cat("  UMAP matrix: not available\n")
    }

    # Milo API requires counts/logcounts assays even though it does not use
    # them for DA testing — supply a sparse placeholder to satisfy the
    # SingleCellExperiment constructor.
    n_cells <- nrow(scvi_matrix)
    n_genes <- 100

    cat("  Creating placeholder counts matrix (", n_genes, " x ", n_cells, ")...\n", sep = "")

    set.seed(42)
    dummy_counts <- Matrix::Matrix(
        data = sample(0:10, n_genes * n_cells, replace = TRUE),
        nrow = n_genes,
        ncol = n_cells,
        dimnames = list(
            paste0("placeholder_gene_", seq_len(n_genes)),
            rownames(scvi_matrix)
        ),
        sparse = TRUE
    )

    metadata <- seurat_obj@meta.data
    metadata <- metadata[rownames(scvi_matrix), , drop = FALSE]

    reduced_dims_list <- list(SCVI = scvi_matrix)
    if (!is.null(umap_matrix)) reduced_dims_list$UMAP <- umap_matrix

    sce <- SingleCellExperiment(
        assays = list(
            counts    = dummy_counts,
            logcounts = dummy_counts
        ),
        colData = metadata,
        reducedDims = SimpleList(reduced_dims_list)
    )

    return(sce)
}

build_neighborhoods <- function(
    sce,
    k = 30,
    d = 50,
    prop = 0.1,
    sample_col = "sample_id",
    reduced_dim = "SCVI",
    refined = TRUE
) {
    available_dims <- names(reducedDims(sce))
    if (!reduced_dim %in% available_dims) {
        match_idx <- which(tolower(available_dims) == tolower(reduced_dim))
        if (length(match_idx) > 0) {
            reduced_dim <- available_dims[match_idx[1]]
        } else {
            stop("ERROR: reducedDim '", reduced_dim, "' not found in SCE.\n",
                 "  Available: ", paste(available_dims, collapse = ", "))
        }
    }

    actual_dims <- ncol(reducedDim(sce, reduced_dim))
    if (d > actual_dims) {
        cat("  Note: Requested d=", d, " exceeds available dims (", actual_dims, "). Using d=", actual_dims, "\n", sep = "")
        d <- actual_dims
    }

    milo_obj <- Milo(sce)

    milo_obj <- buildGraph(milo_obj, k = k, d = d, reduced.dim = reduced_dim)
    milo_obj <- makeNhoods(milo_obj, reduced_dims = reduced_dim, k = k,
                           prop = prop, refined = refined)
    n_nhoods <- ncol(nhoods(milo_obj))
    cat("  ", n_nhoods, " neighborhoods defined\n", sep = "")

    milo_obj <- countCells(
        milo_obj,
        meta.data = as.data.frame(colData(milo_obj)),
        sample = sample_col
    )

    nhood_sizes <- Matrix::colSums(nhoods(milo_obj))
    cat("  Neighborhood sizes: median=", median(nhood_sizes),
        ", min=", min(nhood_sizes),
        ", max=", max(nhood_sizes), "\n", sep = "")

    milo_obj <- calcNhoodDistance(milo_obj, d = d, reduced.dim = reduced_dim)
    milo_obj <- buildNhoodGraph(milo_obj)

    return(milo_obj)
}

get_milo_summary <- function(milo_obj) {
    nhood_sizes <- Matrix::colSums(nhoods(milo_obj))

    list(
        n_cells           = ncol(milo_obj),
        n_neighborhoods   = ncol(nhoods(milo_obj)),
        n_samples         = ncol(nhoodCounts(milo_obj)),
        sample_names      = colnames(nhoodCounts(milo_obj)),
        nhood_size_median = median(nhood_sizes),
        nhood_size_min    = min(nhood_sizes),
        nhood_size_max    = max(nhood_sizes),
        nhood_size_q25    = quantile(nhood_sizes, 0.25),
        nhood_size_q75    = quantile(nhood_sizes, 0.75),
        reduced_dims      = names(reducedDims(milo_obj))
    )
}

print_milo_summary <- function(milo_obj) {
    s <- get_milo_summary(milo_obj)

    cat("============================================================\n")
    cat("Milo Object Summary\n")
    cat("============================================================\n")
    cat("Cells:              ", s$n_cells, "\n")
    cat("Neighborhoods:      ", s$n_neighborhoods, "\n")
    cat("Samples:            ", s$n_samples, "\n")
    cat("Sample names:       ", paste(s$sample_names, collapse = ", "), "\n")
    cat("Neighborhood sizes:\n")
    cat("  Median:           ", s$nhood_size_median, "\n")
    cat("  Range:            ", s$nhood_size_min, " - ", s$nhood_size_max, "\n")
    cat("  IQR:              ", s$nhood_size_q25, " - ", s$nhood_size_q75, "\n")
    cat("Reduced dims:       ", paste(s$reduced_dims, collapse = ", "), "\n")
    cat("============================================================\n")
}
