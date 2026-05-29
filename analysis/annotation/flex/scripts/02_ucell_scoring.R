#!/usr/bin/env Rscript
# 02_ucell_scoring.R — UCell scoring of cells against V1 L2 marker signatures.
#
# Per V1 yaml schema, each non-artifact label has both canonical_markers (broader
# sanity-check set) and identity_markers (tighter identity set). For each cell,
# UCell computes a rank-based score for both signatures per label.
#
# UCell is computed ONCE per cell. Per-cluster aggregation (median) is then
# performed for each resolution column in the clusters CSV. Output is long-format
# with `resolution` as one of the keys.
#
# Inputs:
#   --bundle-dir      scvi_n100 dir with counts.mtx.gz, genes.tsv, cells.tsv
#   --clusters-csv    outputs/clusters/{Compartment}_clusters.csv (with leiden_<res> cols)
#   --resolutions     list of resolution suffixes to aggregate (e.g. 0.3 0.5 1.0 5.0)
#   --yaml            V1 yaml (annotation_v2_{epi,imm,str}.yaml)
#   --out-dir         destination
#   --compartment     Epithelial / Immune / Stromal
#   --ncores          UCell BPPARAM cores (default 4)
#
# Outputs:
#   {Compartment}_ucell_per_cell.csv      cell_id × signature columns
#   {Compartment}_ucell_per_cluster.csv   long: resolution, cluster, label, score_type, median_score, n_cells

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(Matrix)
  library(UCell)
  library(BiocParallel)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser(description = "UCell scoring against V1 L2 markers (multi-resolution aggregation)")
parser$add_argument("--bundle-dir", required = TRUE)
parser$add_argument("--clusters-csv", required = TRUE)
parser$add_argument("--resolutions", nargs = "+", required = TRUE,
                    help = "Resolution suffixes whose leiden columns to aggregate over")
parser$add_argument("--yaml", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--compartment", required = TRUE,
                    choices = c("Epithelial", "Immune", "Stromal"))
parser$add_argument("--ncores", type = "integer", default = 4)
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load V1 yaml signatures ---
cat("Loading V1 yaml:", args$yaml, "\n")
yml <- yaml::read_yaml(args$yaml)
labels <- yml$labels
cat(sprintf("  %d labels in yaml\n", length(labels)))

signatures <- list()
label_meta <- list()
for (label_name in names(labels)) {
  spec <- labels[[label_name]]
  if (isTRUE(spec$is_artifact)) {
    cat(sprintf("  skipping artifact label: %s\n", label_name))
    next
  }
  canonical <- unique(unlist(spec$canonical_markers))
  identity  <- unique(unlist(spec$identity_markers))
  if (length(canonical) > 0) {
    signatures[[paste0(label_name, "__canonical")]] <- canonical
  }
  if (length(identity) > 0) {
    signatures[[paste0(label_name, "__identity")]] <- identity
  }
  label_meta[[label_name]] <- list(
    n_canonical = length(canonical),
    n_identity  = length(identity),
    lineage     = spec$lineage %||% NA_character_
  )
}
cat(sprintf("  %d signatures across %d non-artifact labels\n",
            length(signatures), length(label_meta)))

# --- Load counts (genes × cells) ---
cat("Loading counts from:", args$bundle_dir, "\n")
mtx_path <- file.path(args$bundle_dir, "counts.mtx.gz")
genes <- read.table(file.path(args$bundle_dir, "genes.tsv"),
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
cells <- read.table(file.path(args$bundle_dir, "cells.tsv"),
                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
gene_names <- genes[[1]]
cell_names <- cells[[1]]

X <- readMM(mtx_path)
if (nrow(X) == length(cell_names) && ncol(X) == length(gene_names)) {
  cat("  mtx oriented cells × genes; transposing to genes × cells\n")
  X <- t(X)
} else if (nrow(X) == length(gene_names) && ncol(X) == length(cell_names)) {
  cat("  mtx already oriented genes × cells\n")
} else {
  stop(sprintf("mtx shape (%dx%d) does not match genes (%d) and cells (%d)",
               nrow(X), ncol(X), length(gene_names), length(cell_names)))
}
rownames(X) <- gene_names
colnames(X) <- cell_names
X <- as(X, "CsparseMatrix")
cat(sprintf("  matrix: %d genes × %d cells\n", nrow(X), ncol(X)))

# --- Load clusters ---
clusters <- read.csv(args$clusters_csv, stringsAsFactors = FALSE)
cat(sprintf("  clusters CSV: %d rows × %d cols\n", nrow(clusters), ncol(clusters)))

leiden_cols <- paste0("leiden_", args$resolutions)
missing_cols <- setdiff(leiden_cols, colnames(clusters))
if (length(missing_cols) > 0) {
  stop(sprintf("clusters CSV missing leiden columns: %s",
               paste(missing_cols, collapse = ", ")))
}

# Restrict to cells present in cluster CSV (handles subset test)
common <- intersect(cell_names, clusters$cell_id)
cat(sprintf("  %d common cells between counts and clusters\n", length(common)))
if (length(common) == 0) stop("No overlap between counts and clusters CSV")
X <- X[, common]
clusters <- clusters[match(common, clusters$cell_id), ]
stopifnot(all(clusters$cell_id == colnames(X)))

# --- UCell scoring (single pass over cells) ---
cat(sprintf("Running UCell::ScoreSignatures_UCell on %d cells × %d signatures (cores=%d)...\n",
            ncol(X), length(signatures), args$ncores))
t0 <- Sys.time()
bpp <- if (args$ncores > 1) BiocParallel::MulticoreParam(workers = args$ncores) else BiocParallel::SerialParam()
scores <- UCell::ScoreSignatures_UCell(
  X,
  features = signatures,
  ncores = args$ncores,
  BPPARAM = bpp
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat(sprintf("  UCell complete: %d cells × %d sigs [%.1fs]\n",
            nrow(scores), ncol(scores), elapsed))

# Strip UCell suffix
colnames(scores) <- sub("_UCell$", "", colnames(scores))

# --- Per-cell output (long-lived; downstream render reads this directly) ---
per_cell <- data.frame(
  cell_id = rownames(scores),
  scores,
  check.names = FALSE
)
per_cell_path <- file.path(args$out_dir,
                           paste0(args$compartment, "_ucell_per_cell.csv"))
write.csv(per_cell, per_cell_path, row.names = FALSE)
cat(sprintf("  wrote %s (%d rows × %d cols)\n",
            per_cell_path, nrow(per_cell), ncol(per_cell)))

# --- Per-(resolution × cluster) median aggregation ---
sig_cols <- setdiff(colnames(per_cell), "cell_id")
all_long <- list()

for (i in seq_along(args$resolutions)) {
  res <- args$resolutions[i]
  cluster_col <- leiden_cols[i]
  cat(sprintf("Aggregating medians for resolution %s (col=%s)...\n", res, cluster_col))
  cluster_vals <- as.character(clusters[[cluster_col]])
  medians <- aggregate(per_cell[, sig_cols, drop = FALSE],
                       by = list(cluster = cluster_vals), FUN = median)
  n_per_cluster <- as.data.frame(table(cluster_vals), stringsAsFactors = FALSE)
  colnames(n_per_cluster) <- c("cluster", "n_cells")

  for (sig in sig_cols) {
    parts <- strsplit(sig, "__", fixed = TRUE)[[1]]
    if (length(parts) != 2) next
    label <- parts[1]
    score_type <- parts[2]
    for (k in seq_len(nrow(medians))) {
      cl <- medians$cluster[k]
      all_long[[length(all_long) + 1L]] <- data.frame(
        resolution = res,
        cluster = cl,
        label = label,
        score_type = score_type,
        median_score = medians[k, sig],
        n_cells = n_per_cluster$n_cells[match(cl, n_per_cluster$cluster)],
        stringsAsFactors = FALSE
      )
    }
  }
}

long_df <- do.call(rbind, all_long)
long_path <- file.path(args$out_dir,
                       paste0(args$compartment, "_ucell_per_cluster.csv"))
write.csv(long_df, long_path, row.names = FALSE)
cat(sprintf("  wrote %s (%d rows across %d resolutions)\n",
            long_path, nrow(long_df), length(args$resolutions)))

cat("\nDone.\n")
