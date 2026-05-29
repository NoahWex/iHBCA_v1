#!/usr/bin/env Rscript
# ==============================================================================
# Step 04a: Train SingleR Reference Models from Kumar 2023 HBCA (runs once)
# ==============================================================================
# Loads the Kumar SC source (714K cells, MTX format from CxG) and the Kumar SN
# reference (Seurat RDS), builds trained SingleR models via pseudobulk
# aggregation, and saves two small RDS files consumed by 03b_classify_singler.R.
#
# Usage:
#   Rscript 03a_train_singler_ref.R \
#     --sc_ref_dir /path/to/kumar_converted/ \
#     --sn_ref /path/to/annotated_hbca_nuclei_celltype.rds \
#     [--output_dir /path/to/singler_models/]
#
# Outputs (under output_dir):
#   singler_model_sc.rds  — trained SC model
#   singler_model_sn.rds  — trained SN model
#
# Path conventions:
#   --output_dir defaults to ${CFG_PREPROCESSING_STEP_04}/models/singler_models/
#   The Kumar reference paths are provided as required CLI args because the
#   canonical flow does not use them and they are NOT stored in paths.yaml.
#   See run/run_step_04_train.sh for how SC_REF / SN_REF are passed.
#
# HPC: Tier 5 (128G, 4h) — runs ONCE before the per-sample classify array
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(SingleR)
  library(Matrix)
  library(SummarizedExperiment)
  library(Seurat)
})

parser <- ArgumentParser()
parser$add_argument("--sc_ref_dir", required = TRUE,
                    help = "Kumar SC reference dir (counts.mtx, metadata.csv, genes.txt, barcodes.txt)")
parser$add_argument("--sn_ref", required = TRUE,
                    help = "Kumar SN reference RDS (annotated_hbca_nuclei_celltype.rds)")
parser$add_argument("--output_dir", default = Sys.getenv("CFG_PREPROCESSING_STEP_04", ""),
                    help = "Output directory for trained models. Defaults to CFG_PREPROCESSING_STEP_04/models/singler_models/")
parser$add_argument("--max_cells_per_type", type = "integer", default = 2000L,
                    help = "Max cells per type for training (default: 2000)")
parser$add_argument("--gene_map", default = NULL,
                    help = "Two-column CSV (ensembl,symbol) for Ensembl/symbol namespace conversion")
args <- parser$parse_args()

# Resolve output_dir: if env var gave the step root, append the models subdir
if (args$output_dir == "") {
  stop("--output_dir is required (or set CFG_PREPROCESSING_STEP_04)", call. = FALSE)
}
if (!grepl("singler_model", args$output_dir)) {
  args$output_dir <- file.path(args$output_dir, "models", "singler_models")
}
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Output dir: %s\n", args$output_dir))

# Load gene map if provided (for namespace conversion)
gene_map_ensembl <- NULL
gene_map_symbol <- NULL
if (!is.null(args$gene_map)) {
  gm <- read.csv(args$gene_map, stringsAsFactors = FALSE)
  gene_map_ensembl <- gm$ensembl
  gene_map_symbol <- gm$symbol
  cat(sprintf("Gene map loaded: %d entries\n", nrow(gm)))
}

is_ensembl <- function(ids) any(grepl("^ENSG[0-9]", head(ids, 100)))

# === SC Reference ===
cat("=== Training SC SingleR Reference ===\n")

cat("Loading MTX...\n")
sc_counts <- readMM(file.path(args$sc_ref_dir, "counts.mtx"))
sc_genes   <- readLines(file.path(args$sc_ref_dir, "genes.txt"))
sc_barcodes <- readLines(file.path(args$sc_ref_dir, "barcodes.txt"))
sc_meta    <- read.csv(file.path(args$sc_ref_dir, "metadata.csv"), stringsAsFactors = FALSE)

# Orient: genes x cells
if (nrow(sc_counts) == length(sc_barcodes) && ncol(sc_counts) == length(sc_genes)) {
  sc_counts <- t(sc_counts)
}
rownames(sc_counts) <- sc_genes
colnames(sc_counts) <- sc_barcodes
cat(sprintf("  %d genes x %d cells\n", nrow(sc_counts), ncol(sc_counts)))

sc_labels <- sc_meta$broad_cell_type
cat(sprintf("  %d unique broad types\n", length(unique(sc_labels))))

# Stratified subsample
set.seed(42)
keep_idx <- c()
for (lbl in unique(sc_labels)) {
  type_idx <- which(sc_labels == lbl)
  n <- min(length(type_idx), args$max_cells_per_type)
  keep_idx <- c(keep_idx, sample(type_idx, n))
}
sc_counts <- sc_counts[, keep_idx]
sc_labels <- sc_labels[keep_idx]
cat(sprintf("  Subsampled: %d cells\n", ncol(sc_counts)))

# Log-normalize
sc_lib  <- colSums(sc_counts)
sc_norm <- t(t(sc_counts) / sc_lib * 1e4)
sc_norm@x <- log1p(sc_norm@x)
sc_sce <- SummarizedExperiment(assays = list(logcounts = sc_norm),
                                colData = DataFrame(label = sc_labels))
rm(sc_counts, sc_norm); gc()

cat("  Training...\n")
sc_model <- trainSingleR(ref = sc_sce, labels = sc_labels, de.method = "wilcox")
sc_out <- file.path(args$output_dir, "singler_model_sc.rds")
saveRDS(sc_model, sc_out)
cat(sprintf("  Saved: %s\n", sc_out))
rm(sc_sce, sc_model); gc()

# === SN Reference ===
cat("\n=== Training SN SingleR Reference ===\n")

sn_seu <- readRDS(args$sn_ref)
cat(sprintf("  %d cells\n", ncol(sn_seu)))

tryCatch({ sn_seu <- UpdateSeuratObject(sn_seu) }, error = function(e) NULL)

sn_label_col <- if ("nuc_celltype" %in% colnames(sn_seu@meta.data)) "nuc_celltype" else "celltype"
sn_labels <- sn_seu@meta.data[[sn_label_col]]
cat(sprintf("  Label column: %s (%d types)\n", sn_label_col, length(unique(sn_labels))))

sn_seu <- NormalizeData(sn_seu, verbose = FALSE)
tryCatch({
  sn_sce <- as.SingleCellExperiment(sn_seu)
}, error = function(e) {
  cat("  Rebuilding SCE from counts...\n")
  counts_mat <- GetAssayData(sn_seu, layer = "counts")
  lib <- colSums(counts_mat)
  norm <- t(t(counts_mat) / lib * 1e4)
  norm@x <- log1p(norm@x)
  sn_sce <<- SummarizedExperiment(assays = list(logcounts = norm),
                                   colData = DataFrame(label = sn_labels))
})
rm(sn_seu); gc()

# Convert SN rownames to Ensembl if gene_map provided
if (!is.null(gene_map_ensembl) && !is_ensembl(head(rownames(sn_sce), 100))) {
  sn_genes <- rownames(sn_sce)
  idx <- match(sn_genes, gene_map_symbol)
  has_map <- !is.na(idx)
  cat(sprintf("  SN symbols → Ensembl: mapped %d / %d (dropping %d unmapped)\n",
              sum(has_map), length(sn_genes), sum(!has_map)))
  sn_sce <- sn_sce[has_map, ]
  rownames(sn_sce) <- gene_map_ensembl[idx[has_map]]
}

if (ncol(sn_sce) > 50000) {
  set.seed(42)
  keep <- c()
  for (lbl in unique(sn_labels)) {
    idx <- which(sn_labels == lbl)
    keep <- c(keep, sample(idx, min(length(idx), args$max_cells_per_type)))
  }
  sn_sce <- sn_sce[, keep]
  sn_labels <- sn_labels[keep]
  cat(sprintf("  Subsampled: %d cells\n", ncol(sn_sce)))
}

cat("  Training...\n")
sn_model <- trainSingleR(ref = sn_sce, labels = sn_labels, de.method = "wilcox")
sn_out <- file.path(args$output_dir, "singler_model_sn.rds")
saveRDS(sn_model, sn_out)
cat(sprintf("  Saved: %s\n", sn_out))

cat("\nDone. Models ready for 03b_classify_singler.R.\n")
