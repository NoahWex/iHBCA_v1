#!/usr/bin/env Rscript
# =============================================================================
# group_nhoods_sparse.R — Sparse replacement for miloR::groupNhoods
# =============================================================================
#
# miloR's .group_nhoods_from_adjacency has two O(N^2) dense allocations:
#   1. as.matrix(nhood.adj > 0) — 201K x 201K dense (323 GB)
#   2. sapply(logFC, "-", logFC) — pairwise LFC diff matrix (323 GB)
# Then it calls graph_from_adjacency_matrix on the dense matrix.
#
# This replacement:
#   - Takes the precomputed sparse overlap matrix (crossprod(nhoods), dgCMatrix)
#   - Filters to DA-significant neighborhoods
#   - Applies overlap, discordance, and LFC-delta thresholds on COO triplets
#   - Builds igraph from surviving edge list (NOT adjacency matrix)
#   - Runs Louvain community detection
#   - Returns da_results with NhoodGroup column
#
# All operations stay sparse. The 201K x 201K matrix is never materialized.
#
# Usage:
#   source("group_nhoods_patched.R")
#   da_results <- group_nhoods_sparse(
#     overlap = "/path/to/nhood_overlap.rds",  # or a dgCMatrix
#     da_results = da_res,
#     da.fdr = 0.1,
#     max.lfc.delta = NULL,
#     min_overlap = 1,
#     merge.discord = FALSE
#   )
#
# Tested against: miloR 1.10.0 / R 4.3.3
# Date: 2026-03-26
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(igraph)
})

#' Sparse neighborhood grouping for miloR DA results
#'
#' @param overlap Either a file path to nhood_overlap.rds (dgCMatrix) or the
#'   matrix itself. Must be N_nhoods x N_nhoods symmetric sparse with
#'   overlap counts as values. Diagonal should be zero.
#' @param da_results Data frame with at least columns: logFC, SpatialFDR.
#'   Must have one row per neighborhood (matching overlap matrix dimensions).
#' @param da.fdr FDR threshold for significance. Only significant neighborhoods
#'   are grouped; others get NhoodGroup = NA. Default 0.1.
#' @param max.lfc.delta Maximum absolute logFC difference between connected
#'   neighborhoods. Edges exceeding this are pruned. NULL = no pruning.
#' @param min_overlap Minimum overlap count to retain an edge. Default 1
#'   (any shared cells). Matches miloR's default overlap parameter.
#' @param merge.discord If FALSE (default), prune edges between neighborhoods
#'   with opposite logFC signs. If TRUE, keep discordant edges.
#' @param seed Random seed for Louvain reproducibility. Default 42.
#'
#' @return da_results with NhoodGroup column populated (integer for significant
#'   neighborhoods, NA for non-significant).
group_nhoods_sparse <- function(overlap,
                                da_results,
                                da.fdr = 0.1,
                                max.lfc.delta = NULL,
                                min_overlap = 1,
                                merge.discord = FALSE,
                                seed = 42L) {

  # ---- Validate inputs ----
  if (is.character(overlap)) {
    cat("[group_nhoods_sparse] Loading overlap matrix:", overlap, "\n")
    t_load <- proc.time()
    overlap <- readRDS(overlap)
    cat(sprintf("[group_nhoods_sparse] Loaded in %.1f s\n",
                (proc.time() - t_load)[["elapsed"]]))
  }

  stopifnot(
    inherits(overlap, "dgCMatrix") || inherits(overlap, "dsCMatrix"),
    nrow(overlap) == ncol(overlap),
    nrow(overlap) == nrow(da_results),
    "logFC" %in% names(da_results),
    "SpatialFDR" %in% names(da_results)
  )

  n_total <- nrow(da_results)
  cat(sprintf("[group_nhoods_sparse] %s total neighborhoods, %s nonzero overlap entries\n",
              formatC(n_total, format = "d", big.mark = ","),
              formatC(nnzero(overlap), format = "d", big.mark = ",")))

  # ---- Step 1: Identify significant neighborhoods ----
  sig_idx <- which(da_results$SpatialFDR < da.fdr)
  n_sig <- length(sig_idx)
  cat(sprintf("[group_nhoods_sparse] %s significant neighborhoods (FDR < %g)\n",
              formatC(n_sig, format = "d", big.mark = ","), da.fdr))

  if (n_sig == 0) {
    cat("[group_nhoods_sparse] No significant neighborhoods, returning with NhoodGroup = NA\n")
    da_results$NhoodGroup <- NA_integer_
    return(da_results)
  }

  if (n_sig == 1) {
    cat("[group_nhoods_sparse] Only 1 significant neighborhood, assigning group 1\n")
    da_results$NhoodGroup <- NA_integer_
    da_results$NhoodGroup[sig_idx] <- 1L
    return(da_results)
  }

  # ---- Step 2: Subset overlap matrix to significant neighborhoods ----
  # This reduces the matrix from N_total x N_total to N_sig x N_sig
  t_sub <- proc.time()
  overlap_sig <- overlap[sig_idx, sig_idx, drop = FALSE]
  cat(sprintf("[group_nhoods_sparse] Subsetted to %s x %s (%s nonzero) in %.1f s\n",
              formatC(n_sig, format = "d", big.mark = ","),
              formatC(n_sig, format = "d", big.mark = ","),
              formatC(nnzero(overlap_sig), format = "d", big.mark = ","),
              (proc.time() - t_sub)[["elapsed"]]))

  # ---- Step 3: Extract COO triplets and apply thresholds ----
  t_coo <- proc.time()

  # summary() on a dgCMatrix returns a data.frame with i, j, x columns (1-indexed)
  # For symmetric matrices, ensure we only process upper triangle to avoid
  # double-counting, then igraph makes it undirected
  triplets <- summary(overlap_sig)
  cat(sprintf("[group_nhoods_sparse] Extracted %s triplets in %.1f s\n",
              formatC(nrow(triplets), format = "d", big.mark = ","),
              (proc.time() - t_coo)[["elapsed"]]))

  # Keep upper triangle only (i < j) to avoid duplicate edges
  triplets <- triplets[triplets$i < triplets$j, , drop = FALSE]
  n_before <- nrow(triplets)
  cat(sprintf("[group_nhoods_sparse] Upper triangle: %s edges\n",
              formatC(n_before, format = "d", big.mark = ",")))

  # 3a. Overlap threshold
  if (min_overlap > 0) {
    triplets <- triplets[triplets$x >= min_overlap, , drop = FALSE]
    cat(sprintf("[group_nhoods_sparse] After overlap >= %d: %s edges (%s pruned)\n",
                min_overlap,
                formatC(nrow(triplets), format = "d", big.mark = ","),
                formatC(n_before - nrow(triplets), format = "d", big.mark = ",")))
  }

  if (nrow(triplets) == 0) {
    cat("[group_nhoods_sparse] No edges survive overlap threshold. Each neighborhood is its own group.\n")
    da_results$NhoodGroup <- NA_integer_
    da_results$NhoodGroup[sig_idx] <- seq_len(n_sig)
    return(da_results)
  }

  # 3b. Map logFC to the sig-indexed neighborhoods
  lfc_sig <- da_results$logFC[sig_idx]
  triplets$lfc_i <- lfc_sig[triplets$i]
  triplets$lfc_j <- lfc_sig[triplets$j]

  # 3c. Discordance pruning
  if (!merge.discord) {
    n_pre <- nrow(triplets)
    # Remove edges where signs differ (opposite direction effects)
    # Keep edges where both are zero (no effect) or same sign
    concordant <- sign(triplets$lfc_i) == sign(triplets$lfc_j)
    triplets <- triplets[concordant, , drop = FALSE]
    cat(sprintf("[group_nhoods_sparse] After discordance pruning: %s edges (%s discordant removed)\n",
                formatC(nrow(triplets), format = "d", big.mark = ","),
                formatC(n_pre - nrow(triplets), format = "d", big.mark = ",")))
  }

  if (nrow(triplets) == 0) {
    cat("[group_nhoods_sparse] No edges survive discordance pruning. Each neighborhood is its own group.\n")
    da_results$NhoodGroup <- NA_integer_
    da_results$NhoodGroup[sig_idx] <- seq_len(n_sig)
    return(da_results)
  }

  # 3d. LFC delta pruning
  if (!is.null(max.lfc.delta)) {
    n_pre <- nrow(triplets)
    lfc_diff <- abs(triplets$lfc_i - triplets$lfc_j)
    triplets <- triplets[lfc_diff <= max.lfc.delta, , drop = FALSE]
    cat(sprintf("[group_nhoods_sparse] After max.lfc.delta <= %g: %s edges (%s pruned)\n",
                max.lfc.delta,
                formatC(nrow(triplets), format = "d", big.mark = ","),
                formatC(n_pre - nrow(triplets), format = "d", big.mark = ",")))
  }

  if (nrow(triplets) == 0) {
    cat("[group_nhoods_sparse] No edges survive LFC delta pruning. Each neighborhood is its own group.\n")
    da_results$NhoodGroup <- NA_integer_
    da_results$NhoodGroup[sig_idx] <- seq_len(n_sig)
    return(da_results)
  }

  # ---- Step 4: Build igraph from edge list and run Louvain ----
  t_graph <- proc.time()

  # Build edge list using sig-local indices (1..n_sig)
  el <- data.frame(from = triplets$i, to = triplets$j)
  g <- graph_from_data_frame(el, directed = FALSE, vertices = data.frame(id = seq_len(n_sig)))

  # Weight edges by overlap count (Louvain uses edge weights)
  E(g)$weight <- triplets$x

  cat(sprintf("[group_nhoods_sparse] Graph: %s vertices, %s edges, built in %.1f s\n",
              formatC(vcount(g), format = "d", big.mark = ","),
              formatC(ecount(g), format = "d", big.mark = ","),
              (proc.time() - t_graph)[["elapsed"]]))

  # Louvain community detection
  set.seed(seed)
  t_louv <- proc.time()
  communities <- cluster_louvain(g)
  n_groups <- max(membership(communities))
  cat(sprintf("[group_nhoods_sparse] Louvain: %d communities detected in %.1f s\n",
              n_groups, (proc.time() - t_louv)[["elapsed"]]))

  # ---- Step 5: Map group assignments back to da_results ----
  da_results$NhoodGroup <- NA_integer_
  memberships <- membership(communities)

  # memberships is keyed by vertex id (1..n_sig) which maps to sig_idx
  da_results$NhoodGroup[sig_idx] <- as.integer(memberships)

  cat(sprintf("[group_nhoods_sparse] Assigned %d groups to %s significant neighborhoods\n",
              n_groups, formatC(n_sig, format = "d", big.mark = ",")))

  # Summary stats
  group_sizes <- table(da_results$NhoodGroup)
  cat("[group_nhoods_sparse] Group sizes:\n")
  cat(sprintf("  min=%d, median=%d, max=%d\n",
              min(group_sizes), median(group_sizes), max(group_sizes)))

  return(da_results)
}


# =============================================================================
# CLI entry point — for run script usage
# =============================================================================
if (sys.nframe() == 0 && !interactive()) {
  library(argparse)

  parser <- ArgumentParser(
    description = "Sparse neighborhood grouping for miloR DA results (201K nhoods)"
  )
  parser$add_argument("--overlap", required = TRUE,
                      help = "Path to nhood_overlap.rds (dgCMatrix, N x N)")
  parser$add_argument("--da-results", required = TRUE,
                      help = "Path to da_results.csv (must have logFC, SpatialFDR columns)")
  parser$add_argument("--output", required = TRUE,
                      help = "Output path for grouped da_results.csv")
  parser$add_argument("--da-fdr", type = "double", default = 0.1,
                      help = "FDR threshold for significance [default: 0.1]")
  parser$add_argument("--max-lfc-delta", type = "double", default = NULL,
                      help = "Max abs logFC difference between connected nhoods [default: NULL]")
  parser$add_argument("--min-overlap", type = "integer", default = 1L,
                      help = "Minimum overlap count to retain edge [default: 1]")
  parser$add_argument("--merge-discord", action = "store_true", default = FALSE,
                      help = "If set, keep edges between opposite-sign neighborhoods")
  parser$add_argument("--seed", type = "integer", default = 42L,
                      help = "Random seed for Louvain [default: 42]")

  args <- parser$parse_args()

  cat("=== group_nhoods_sparse CLI ===\n")
  cat("Overlap:", args$overlap, "\n")
  cat("DA results:", args$da_results, "\n")
  cat("Output:", args$output, "\n")
  cat("FDR threshold:", args$da_fdr, "\n")
  cat("Max LFC delta:", ifelse(is.null(args$max_lfc_delta), "NULL", args$max_lfc_delta), "\n")
  cat("Min overlap:", args$min_overlap, "\n")
  cat("Merge discord:", args$merge_discord, "\n")
  cat("Seed:", args$seed, "\n")
  cat("\n")

  # Load DA results
  cat("Loading DA results...\n")
  da_res <- read.csv(args$da_results, stringsAsFactors = FALSE)
  cat(sprintf("Loaded %s rows x %d columns\n",
              formatC(nrow(da_res), format = "d", big.mark = ","),
              ncol(da_res)))

  # Run grouping
  t_total <- proc.time()
  da_grouped <- group_nhoods_sparse(
    overlap = args$overlap,
    da_results = da_res,
    da.fdr = args$da_fdr,
    max.lfc.delta = args$max_lfc_delta,
    min_overlap = args$min_overlap,
    merge.discord = args$merge_discord,
    seed = args$seed
  )
  elapsed_total <- (proc.time() - t_total)[["elapsed"]]

  # Save
  cat(sprintf("\nSaving to: %s\n", args$output))
  write.csv(da_grouped, args$output, row.names = FALSE)
  cat(sprintf("File size: %.1f MB\n", file.size(args$output) / 1e6))

  # Summary YAML
  yaml_path <- sub("\\.csv$", "_grouping_summary.yaml", args$output)
  n_sig <- sum(!is.na(da_grouped$NhoodGroup))
  n_groups <- length(unique(na.omit(da_grouped$NhoodGroup)))
  yaml_lines <- c(
    "grouping_summary:",
    sprintf("  script: group_nhoods_patched.R"),
    sprintf("  date: '%s'", Sys.Date()),
    sprintf("  total_nhoods: %d", nrow(da_grouped)),
    sprintf("  significant_nhoods: %d", n_sig),
    sprintf("  n_groups: %d", n_groups),
    sprintf("  da_fdr: %g", args$da_fdr),
    sprintf("  max_lfc_delta: %s", ifelse(is.null(args$max_lfc_delta), "null", args$max_lfc_delta)),
    sprintf("  min_overlap: %d", args$min_overlap),
    sprintf("  merge_discord: %s", tolower(as.character(args$merge_discord))),
    sprintf("  seed: %d", args$seed),
    sprintf("  elapsed_seconds: %.1f", elapsed_total),
    sprintf("  input_overlap: '%s'", args$overlap),
    sprintf("  input_da_results: '%s'", args$da_results),
    sprintf("  output: '%s'", args$output)
  )
  writeLines(yaml_lines, yaml_path)
  cat("Summary written to:", yaml_path, "\n")

  cat(sprintf("\nComplete in %.1f s (%.1f min)\n", elapsed_total, elapsed_total / 60))
}
