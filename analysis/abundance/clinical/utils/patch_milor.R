# patch_milor.R — Expose BPPARAM + BNPARAM in miloR::makeNhoods
#
# Why: makeNhoods(refined=TRUE) hardcodes single-threaded exact KNN via
#      BiocNeighbors::findKNN(). At 2.12M cells this takes 5-7 days.
#      findKNN already supports BPPARAM (parallelization) and BNPARAM
#      (approximate methods like Annoy/HNSW) — miloR just doesn't pass
#      them through.
#
#      Additionally, the second KNN call in .refined_sampling() uses
#      k = n_medians + 1 (~212K at 2.12M cells, prop=0.1) to skip past
#      median coordinates and find the nearest real cell. This allocates
#      a 212K x 212K result matrix (~180 GB) and is the primary bottleneck.
#      Replaced with queryKNN(X=real_cells, query=medians, k=1) which
#      returns the nearest real cell directly — same result, 1.7 MB output.
#
# What changes:
#   1. .refined_sampling() gains BPPARAM + BNPARAM args, passes to findKNN/queryKNN
#   2. .refined_sampling() KNN call 2 replaced: findKNN(combined, k=N+1) -> queryKNN(real, query=medians, k=1)
#   3. makeNhoods() gains BPPARAM + BNPARAM args, passes to .refined_sampling()
#   4. makeNhoods() nhood matrix: sequential for-loop sparse assignment -> vectorized sparseMatrix()
#   5. calcNhoodDistance(): O(N_entries × N_nhoods) filter -> O(N_entries) split()
#
# Defaults: SerialParam() + KmknnParam() — identical to original behavior.
# The patch ONLY changes behavior when non-default args are passed.
#
# Usage:
#   source("patch_milor.R")  # before calling makeNhoods()
#   milo <- makeNhoods(milo, ...,
#                      BPPARAM = MulticoreParam(16),
#                      BNPARAM = AnnoyParam(ntrees = 50))
#
# Tested against: miloR 1.10.0 / BiocNeighbors 1.20.2 / R 4.3.3
# Date: 2026-03-12
# Plan: 04_DifferentialAbundance/da_v1_milo_approx
#
# If you're reading this, congratulations on finding the only part of miloR
# that thinks 2 million cells should be a "quick" exact KNN problem.

suppressPackageStartupMessages({
  library(miloR)
  library(BiocNeighbors)
  library(BiocParallel)
})

# --- Version guard ---
.milor_version <- packageVersion("miloR")
if (.milor_version != "1.10.0") {
  warning(
    "patch_milor.R was written for miloR 1.10.0 but found ",
    as.character(.milor_version),
    ". The patch may not work correctly — verify .refined_sampling source."
  )
}

# =============================================================================
# ORIGINAL .refined_sampling (miloR 1.10.0, verbatim)
# =============================================================================
# function (random_vertices, X_reduced_dims, k)
# {
#     message("Running refined sampling with reduced_dim")
#     vertex.knn <- findKNN(X = X_reduced_dims, k = k, subset = as.vector(random_vertices),
#         get.index = TRUE, get.distance = FALSE)
#     nh_reduced_dims <- t(apply(vertex.knn$index, 1, function(x) colMedians(X_reduced_dims[x,
#         ])))
#     if (is.null(rownames(X_reduced_dims))) {
#         warning("Rownames not set on reducedDims - setting to row indices")
#         rownames(X_reduced_dims) <- as.character(seq_len(nrow(X_reduced_dims)))
#     }
#     colnames(nh_reduced_dims) <- colnames(X_reduced_dims)
#     rownames(nh_reduced_dims) <- paste0("nh_", seq_len(nrow(nh_reduced_dims)))
#     all_reduced_dims <- rbind(nh_reduced_dims, X_reduced_dims)
#     nn_mat <- findKNN(all_reduced_dims, k = nrow(nh_reduced_dims) +
#         1, subset = rownames(nh_reduced_dims))[["index"]]
#     nh_ixs <- seq_len(nrow(nh_reduced_dims))
#     i = 1
#     sampled_vertices <- rep(0, nrow(nn_mat))
#     while (any(sampled_vertices <= max(nh_ixs))) {
#         update_ix <- which(sampled_vertices <= max(nh_ixs))
#         sampled_vertices[update_ix] <- nn_mat[update_ix, i]
#         i <- i + 1
#     }
#     sampled_vertices <- sampled_vertices - max(nh_ixs)
#     return(sampled_vertices)
# }

# =============================================================================
# PATCHED .refined_sampling
# =============================================================================
# Changes vs original:
#   1. BPPARAM + BNPARAM passed to findKNN/queryKNN
#   2. KNN call 2 replaced: findKNN(combined, k=N+1) -> queryKNN(real, query=medians, k=1)
#      Original allocates N_medians x (N_medians+1) result matrix (~180 GB at 2.12M cells).
#      queryKNN(k=1) returns N_medians x 1 (~1.7 MB). Same result: nearest real cell.
#   3. Progress reporting on each phase
.refined_sampling_patched <- function(random_vertices, X_reduced_dims, k,
                                      BPPARAM = BiocParallel::SerialParam(),
                                      BNPARAM = BiocNeighbors::KmknnParam()) {
  n_cells <- nrow(X_reduced_dims)
  n_sampled <- length(random_vertices)
  message("Running refined sampling with reduced_dim")
  message(sprintf("  [refine] %s cells, %s sampled vertices, k=%d",
                  formatC(n_cells, format = "d", big.mark = ","),
                  formatC(n_sampled, format = "d", big.mark = ","), k))

  # Step 1: find k nearest neighbors for each sampled vertex
  message("  [refine 1/3] Finding k neighbors for sampled vertices...")
  t1 <- proc.time()
  vertex.knn <- findKNN(
    X = X_reduced_dims, k = k,
    subset = as.vector(random_vertices),
    get.index = TRUE, get.distance = FALSE,
    BPPARAM = BPPARAM, BNPARAM = BNPARAM
  )
  elapsed1 <- (proc.time() - t1)[["elapsed"]]
  message(sprintf("  [refine 1/3] Done in %.1f min", elapsed1 / 60))

  # Step 2: compute median coordinate per neighborhood
  message(sprintf("  [refine 2/3] Computing median coordinates (%s neighborhoods x %d dims)...",
                  formatC(n_sampled, format = "d", big.mark = ","),
                  ncol(X_reduced_dims)))
  t2 <- proc.time()
  nh_reduced_dims <- t(apply(vertex.knn$index, 1, function(x)
    colMedians(X_reduced_dims[x, ])))
  elapsed2 <- (proc.time() - t2)[["elapsed"]]
  message(sprintf("  [refine 2/3] Done in %.1f min", elapsed2 / 60))

  if (is.null(rownames(X_reduced_dims))) {
    warning("Rownames not set on reducedDims - setting to row indices")
    rownames(X_reduced_dims) <- as.character(seq_len(nrow(X_reduced_dims)))
  }
  colnames(nh_reduced_dims) <- colnames(X_reduced_dims)

  # Step 3: find nearest real cell to each median coordinate
  # Original miloR: findKNN(rbind(medians, cells), k=n_medians+1, subset=medians)
  #   -> allocates n_medians x (n_medians+1) matrix, then walks to skip medians
  # Optimized: queryKNN(X=cells, query=medians, k=1) -> nearest real cell directly
  message(sprintf("  [refine 3/3] queryKNN: nearest real cell for each of %s medians (k=1)...",
                  formatC(nrow(nh_reduced_dims), format = "d", big.mark = ",")))
  t3 <- proc.time()
  nn_result <- queryKNN(
    X = X_reduced_dims,
    query = nh_reduced_dims,
    k = 1,
    BPPARAM = BPPARAM, BNPARAM = BNPARAM
  )
  sampled_vertices <- as.integer(nn_result$index[, 1])
  elapsed3 <- (proc.time() - t3)[["elapsed"]]
  message(sprintf("  [refine 3/3] Done in %.1f min", elapsed3 / 60))

  total <- elapsed1 + elapsed2 + elapsed3
  message(sprintf("  [refine] Complete: %.1f min total (KNN1=%.0f%%, median=%.0f%%, queryKNN=%.0f%%)",
                  total / 60,
                  100 * elapsed1 / total, 100 * elapsed2 / total,
                  100 * elapsed3 / total))
  return(sampled_vertices)
}

# =============================================================================
# ORIGINAL makeNhoods (miloR 1.10.0, verbatim) — see dump_milor_source.R output
# for full source. Only the call to .refined_sampling changes.
# =============================================================================

# =============================================================================
# PATCHED makeNhoods
# =============================================================================
# Changes: Added BPPARAM + BNPARAM to signature; passed through to
# .refined_sampling() when refinement_scheme == "reduced_dim".
makeNhoods_patched <- function(x, prop = 0.1, k = 21, d = 30, refined = TRUE,
                               reduced_dims = "PCA",
                               refinement_scheme = "reduced_dim",
                               BPPARAM = BiocParallel::SerialParam(),
                               BNPARAM = BiocNeighbors::KmknnParam()) {
  if (is(x, "Milo")) {
    message("Checking valid object")
    if (!miloR:::.valid_graph(graph(x))) {
      stop("Not a valid Milo object - graph is missing. Please run buildGraph() first.")
    }
    graph <- graph(x)
    if (isTRUE(refined) & refinement_scheme == "reduced_dim") {
      X_reduced_dims <- reducedDim(x, reduced_dims)
      if (d > ncol(X_reduced_dims)) {
        warning("Specified d is higher than the total number of dimensions in reducedDim(x, reduced_dims).\n                        Falling back to using",
                ncol(X_reduced_dims), "dimensions\n")
        d <- ncol(X_reduced_dims)
      }
      X_reduced_dims <- X_reduced_dims[, seq_len(d)]
      mat_cols <- ncol(x)
      match.ids <- all(rownames(X_reduced_dims) == colnames(x))
      if (!match.ids) {
        stop("Rownames of reduced dimensions do not match cell IDs")
      }
    }
  } else if (is(x, "igraph")) {
    if (isTRUE(refined) & refinement_scheme == "reduced_dim" &
        !is.matrix(reduced_dims)) {
      stop("No reduced dimensions matrix provided - required for refined sampling with refinement_scheme = reduced_dim.")
    }
    graph <- x
    if (isTRUE(refined) & refinement_scheme == "reduced_dim") {
      X_reduced_dims <- reduced_dims
      mat_cols <- nrow(X_reduced_dims)
      if (is.null(rownames(X_reduced_dims))) {
        stop("Reduced dim rownames are missing - required to assign cell IDs to neighbourhoods")
      }
    }
    if (isTRUE(refined) & refinement_scheme == "graph" &
        is.matrix(reduced_dims)) {
      warning("Ignoring reduced dimensions matrix because refinement_scheme = graph was selected.")
    }
  } else {
    stop("Data format: ", class(x), " not recognised. Should be Milo or igraph.")
  }

  random_vertices <- miloR:::.sample_vertices(graph, prop, return.vertices = TRUE)

  if (isFALSE(refined)) {
    sampled_vertices <- random_vertices
  } else if (isTRUE(refined)) {
    if (refinement_scheme == "reduced_dim") {
      # PATCHED: pass BPPARAM + BNPARAM through
      sampled_vertices <- .refined_sampling_patched(
        random_vertices, X_reduced_dims, k,
        BPPARAM = BPPARAM, BNPARAM = BNPARAM
      )
    } else if (refinement_scheme == "graph") {
      sampled_vertices <- miloR:::.graph_refined_sampling(random_vertices, graph)
    } else {
      stop("When refined == TRUE, refinement_scheme must be one of \"reduced_dim\" or \"graph\".")
    }
  } else {
    stop("refined must be TRUE or FALSE")
  }

  sampled_vertices <- unique(sampled_vertices)
  n_sv <- length(sampled_vertices)

  # --- Build nhood sparse matrix ---
  # Original: for loop with 212K individual sparse matrix column assignments.
  # Each assignment reallocates CSC storage → O(N_nhoods²) in practice.
  # Fix: collect all (row, col) pairs via vectorized neighborhood(), then
  # build sparse matrix in one shot with sparseMatrix().
  message(sprintf("  [nhoods] Building neighborhood matrix (%s vertices)...",
                  formatC(n_sv, format = "d", big.mark = ",")))
  t_nhood <- proc.time()

  # igraph::neighborhood with vector of nodes returns a list of integer vectors
  nhood_lists <- igraph::neighborhood(graph, order = 1, nodes = sampled_vertices)

  # Build COO triplets: row indices and column indices for all non-zero entries
  col_indices <- rep(seq_len(n_sv), lengths(nhood_lists))
  row_indices <- unlist(nhood_lists)

  if (is(x, "Milo")) {
    n_rows <- ncol(x)
    rnames <- colnames(x)
  } else if (is(x, "igraph")) {
    n_rows <- length(igraph::V(x))
    v.class <- igraph::V(graph)$name
    if (is.null(v.class) & refinement_scheme == "reduced_dim") {
      rnames <- rownames(X_reduced_dims)
    } else if (!is.null(v.class)) {
      rnames <- igraph::V(graph)$name
    } else {
      rnames <- NULL
    }
  }

  nh_mat <- Matrix::sparseMatrix(
    i = row_indices, j = col_indices,
    x = rep(1, length(row_indices)),
    dims = c(n_rows, n_sv),
    dimnames = list(rnames, as.character(sampled_vertices))
  )

  elapsed_nhood <- (proc.time() - t_nhood)[["elapsed"]]
  message(sprintf("  [nhoods] Done in %.1f min (%s entries)",
                  elapsed_nhood / 60,
                  formatC(length(row_indices), format = "d", big.mark = ",")))

  if (is(x, "Milo")) {
    nhoodIndex(x) <- as(sampled_vertices, "list")
    nhoods(x) <- nh_mat
    return(x)
  } else {
    return(nh_mat)
  }
}

# =============================================================================
# PATCHED calcNhoodDistance
# =============================================================================
# Original: sapply over N_nhoods, each iteration scans ALL sparse matrix entries
#   to filter rows for that neighborhood: non.zero.nhoods[non.zero.nhoods[,"col"] == X, "row"]
#   Cost: O(N_entries × N_nhoods) — at 212K nhoods × 6.4M entries = 1.3 trillion comparisons
#
# Fix: pre-split row indices by column with split() — O(N_entries) once.
#   Then iterate over pre-grouped list. Also adds progress reporting.
calcNhoodDistance_patched <- function(x, d, reduced.dim = NULL,
                                      use.assay = "logcounts") {
  if (!is(x, "Milo")) stop("Input is not a valid Milo object")

  if ((length(reducedDimNames(x)) == 0) & is.null(reduced.dim)) {
    message("Computing PCA on input")
    x_pca <- irlba::prcomp_irlba(t(SummarizedExperiment::assay(x, use.assay)),
                                  n = min(d + 1, ncol(x) - 1),
                                  scale. = TRUE, center = TRUE)
    reducedDim(x, "PCA") <- x_pca$x
  } else if ((length(reducedDimNames(x)) == 0) & is.character(reduced.dim)) {
    stop(reduced.dim, " not found in reducedDim slot")
  }

  # Resolve which reduced dim to use
  if (is.character(reduced.dim)) {
    if (!any(names(reducedDims(x)) %in% reduced.dim)) {
      stop(reduced.dim, " not found in the reducedDim slot")
    }
    red_dim_mat <- reducedDim(x, reduced.dim)[, seq_len(d), drop = FALSE]
  } else if (is(reduced.dim, "matrix")) {
    red_dim_mat <- reduced.dim[, seq_len(d), drop = FALSE]
  } else if (is.null(reduced.dim)) {
    if (any(names(reducedDims(x)) %in% "PCA")) {
      red_dim_mat <- reducedDim(x, "PCA")[, seq_len(d), drop = FALSE]
    } else {
      stop("No reduced.dim slot specified")
    }
  }

  n_nhoods <- ncol(nhoods(x))
  message(sprintf("  [nhoodDist] Computing distances for %s neighborhoods...",
                  formatC(n_nhoods, format = "d", big.mark = ",")))
  t_start <- proc.time()

  # Pre-split: O(N_entries) instead of O(N_entries × N_nhoods)
  non.zero.nhoods <- which(nhoods(x) != 0, arr.ind = TRUE)
  message(sprintf("  [nhoodDist] Sparse entries: %s, pre-splitting by column...",
                  formatC(nrow(non.zero.nhoods), format = "d", big.mark = ",")))
  nhood_members <- split(non.zero.nhoods[, "row"], non.zero.nhoods[, "col"])

  # Original .calc_distance: double R loop (lapply × apply) computing pairwise
  # Euclidean distances, then assembles a sparseMatrix per nhood.
  # Replaced with dist() — single C-level call. At 201K nhoods the per-call
  # overhead compounds heavily.
  nhood.dists <- lapply(nhood_members, function(rows) {
    sub <- red_dim_mat[rows, , drop = FALSE]
    d_mat <- as.matrix(dist(sub))
    Matrix::Matrix(d_mat, sparse = TRUE)
  })
  names(nhood.dists) <- nhoodIndex(x)

  elapsed <- (proc.time() - t_start)[["elapsed"]]
  message(sprintf("  [nhoodDist] Done in %.1f min", elapsed / 60))

  nhoodDistances(x) <- nhood.dists
  return(x)
}

# =============================================================================
# PATCHED annotateNhoods
# =============================================================================
# Original: vapply over N_nhoods, each doing which(nhoods(x)[,n]==1) on a
#   2.12M-row sparse column, then table(). O(N_cells × N_nhoods).
# Fix: pre-split sparse entries by column, then table() per group.
annotateNhoods_patched <- function(x, da.res, coldata_col, subset.nhoods = NULL) {
  if (!is(x, "Milo")) stop("Unrecognised input type - must be of class Milo")
  if (!coldata_col %in% names(colData(x))) {
    stop(coldata_col, " is not a column in colData(x)")
  }

  if (is.null(subset.nhoods)) {
    if (ncol(nhoods(x)) != nrow(da.res)) {
      stop("the number of rows in da.res does not match the number of neighbourhoods in nhoods(x)")
    }
    nh_mat <- nhoods(x)
  } else {
    nh_mat <- nhoods(x)[, subset.nhoods, drop = FALSE]
  }

  anno_vec <- colData(x)[[coldata_col]]
  if (!is.factor(anno_vec)) {
    message("Converting ", coldata_col, " to factor...")
    anno_vec <- factor(anno_vec, levels = unique(anno_vec))
  }

  n_nhoods <- ncol(nh_mat)
  message(sprintf("  [annotate] %s neighborhoods, %s levels...",
                  formatC(n_nhoods, format = "d", big.mark = ","),
                  length(levels(anno_vec))))
  t_start <- proc.time()

  # Pre-split: O(N_entries) instead of O(N_entries × N_nhoods)
  nz <- which(nh_mat != 0, arr.ind = TRUE)
  nz_anno <- anno_vec[nz[, "row"]]
  anno_by_nhood <- split(nz_anno, nz[, "col"])

  all_levels <- levels(anno_vec)
  max_val <- character(n_nhoods)
  max_frac <- numeric(n_nhoods)

  for (i in seq_len(n_nhoods)) {
    nh_id <- as.character(i)
    if (nh_id %in% names(anno_by_nhood)) {
      tb <- table(anno_by_nhood[[nh_id]])
      idx <- which.max(tb)
      max_val[i] <- names(tb)[idx]
      max_frac[i] <- tb[idx] / sum(tb)
    } else {
      max_val[i] <- NA_character_
      max_frac[i] <- NA_real_
    }
  }

  elapsed <- (proc.time() - t_start)[["elapsed"]]
  message(sprintf("  [annotate] Done in %.1f min", elapsed / 60))

  da.res[coldata_col] <- max_val
  da.res[paste0(coldata_col, "_fraction")] <- max_frac
  return(da.res)
}

# =============================================================================
# Patched .build_nhood_adjacency — sparse-safe overlap thresholding
# =============================================================================
# Why: original does `nh_intersect_mat[nh_intersect_mat < overlap] <- 0`.
#      The LHS expression generates linear indices via Matrix:::int2i, which
#      coerces to integer. At >= ~46K nhoods, n_nhoods^2 exceeds
#      .Machine$integer.max (2.15e9), the indices overflow to NA, and the
#      `[<-` fails with "missing value where TRUE/FALSE needed".
#
#      Fix: zero @x directly, then drop0 to compact. Skips the mask entirely.
#      At overlap=1 (default) the threshold is a no-op (sparse stores no zeros)
#      so the work is skipped altogether.
.build_nhood_adjacency_patched <- function(nhoods, overlap = 1) {
  nh_intersect_mat <- Matrix::crossprod(nhoods)
  if (overlap > 1) {
    nh_intersect_mat@x[nh_intersect_mat@x < overlap] <- 0
    nh_intersect_mat <- Matrix::drop0(nh_intersect_mat)
  }
  rownames(nh_intersect_mat) <- colnames(nhoods)
  colnames(nh_intersect_mat) <- colnames(nhoods)
  nh_intersect_mat
}

# =============================================================================
# Patched .group_nhoods_from_adjacency — edge-list filtering
# =============================================================================
# Why: original constructs dense N x N matrices for discordance and lfc.delta
#      filters, then uses sparse-mask `[<-` (same overflow as above) and finally
#      densifies via as.matrix() before graph_from_adjacency_matrix. At 60K+
#      nhoods this both crashes and would need ~30 GB of dense matrices.
#
#      Fix: pull the sparse adjacency to an edge list (Matrix::summary), apply
#      all filters as vectorized ops on the edge list, rebuild a sparse
#      adjacency from the kept edges, hand to igraph (which accepts sparse).
#      No dense N x N anywhere.
.group_nhoods_from_adjacency_patched <- function(nhs, nhood.adj, da.res, is.da,
                                                 merge.discord = FALSE,
                                                 max.lfc.delta = NULL,
                                                 overlap = 1,
                                                 subset.nhoods = NULL,
                                                 orig.behave = TRUE) {
  if (is.null(colnames(nhs))) {
    warning("No names attributed to nhoods. Converting indices to names")
    colnames(nhs) <- as.character(seq_len(ncol(nhs)))
  }

  if (!is.null(subset.nhoods)) {
    if (mode(subset.nhoods) %in% c("character", "logical", "numeric")) {
      sub.log <- if (mode(subset.nhoods) == "character") {
        colnames(nhs) %in% subset.nhoods
      } else if (mode(subset.nhoods) == "numeric") {
        colnames(nhs) %in% colnames(nhs)[subset.nhoods]
      } else subset.nhoods
      nhood.adj <- nhood.adj[sub.log, sub.log]
      if (length(is.da) == ncol(nhs)) {
        nhs <- nhs[sub.log]
        is.da <- is.da[sub.log]
        da.res <- da.res[sub.log, ]
      } else stop("Subsetting `is.da` vector length does not equal nhoods length")
    } else stop("Incorrect subsetting vector provided:", class(subset.nhoods))
  } else if (length(is.da) != ncol(nhood.adj)) {
    stop("Subsetting `is.da` vector length is not the same dimension as adjacency")
  }

  n_nh <- ncol(nhood.adj)
  nh_names <- colnames(nhood.adj)

  # Edge list from upper triangle of sparse adjacency
  edges <- Matrix::summary(methods::as(nhood.adj, "TsparseMatrix"))
  edges <- edges[edges$i < edges$j, , drop = FALSE]

  if (overlap > 1) edges <- edges[edges$x >= overlap, , drop = FALSE]

  if (isFALSE(merge.discord) && nrow(edges) > 0) {
    sign_i <- sign(da.res$logFC[edges$i])
    sign_j <- sign(da.res$logFC[edges$j])
    discord <- (sign_i != sign_j) & (sign_i != 0) & (sign_j != 0)
    if (isTRUE(orig.behave)) {
      # Original: only drop discordant edges between two DA nhoods
      drop <- discord & is.da[edges$i] & is.da[edges$j]
    } else {
      drop <- discord
    }
    edges <- edges[!drop, , drop = FALSE]
  }

  if (!is.null(max.lfc.delta) && nrow(edges) > 0) {
    lfc_diff <- abs(da.res$logFC[edges$i] - da.res$logFC[edges$j])
    edges <- edges[lfc_diff <= max.lfc.delta, , drop = FALSE]
  }

  # Symmetric sparse adjacency from kept edges
  if (nrow(edges) > 0) {
    adj_filt <- Matrix::sparseMatrix(
      i = c(edges$i, edges$j),
      j = c(edges$j, edges$i),
      x = 1,
      dims = c(n_nh, n_nh),
      symmetric = FALSE
    )
  } else {
    adj_filt <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0), dims = c(n_nh, n_nh)
    )
  }
  rownames(adj_filt) <- colnames(adj_filt) <- nh_names

  g <- igraph::graph_from_adjacency_matrix(adj_filt, mode = "undirected", diag = FALSE)
  groups <- igraph::cluster_louvain(g)$membership
  names(groups) <- nh_names
  groups
}

# =============================================================================
# Apply patches
# =============================================================================

# Update internal namespace functions
assignInNamespace(".refined_sampling", .refined_sampling_patched, ns = "miloR")
assignInNamespace("makeNhoods", makeNhoods_patched, ns = "miloR")
assignInNamespace("calcNhoodDistance", calcNhoodDistance_patched, ns = "miloR")
assignInNamespace("annotateNhoods", annotateNhoods_patched, ns = "miloR")
assignInNamespace(".build_nhood_adjacency", .build_nhood_adjacency_patched, ns = "miloR")
assignInNamespace(".group_nhoods_from_adjacency", .group_nhoods_from_adjacency_patched, ns = "miloR")

# Also update the exported copies in the package search path.
# assignInNamespace modifies the namespace, but library(miloR) already placed
# the originals in the package:miloR environment on the search path. Without
# this, calls resolve to the pre-loaded unpatched versions.
pkg_env <- as.environment("package:miloR")
for (.fn_name in c("makeNhoods", "calcNhoodDistance", "annotateNhoods")) {
  .fn_patched <- get(paste0(.fn_name, "_patched"))
  unlockBinding(.fn_name, pkg_env)
  assign(.fn_name, .fn_patched, envir = pkg_env)
  lockBinding(.fn_name, pkg_env)
}

message("[patch_milor] Patched makeNhoods + .refined_sampling + calcNhoodDistance + annotateNhoods + .build_nhood_adjacency + .group_nhoods_from_adjacency in miloR ",
        as.character(.milor_version),
        " — BPPARAM/BNPARAM + vectorized nhood distances + vectorized annotation + sparse-safe nhood graph")
