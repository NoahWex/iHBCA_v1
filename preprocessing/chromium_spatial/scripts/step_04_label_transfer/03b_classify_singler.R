#!/usr/bin/env Rscript
# ==============================================================================
# Step 04b: Classify FLEX cells against pre-trained SingleR models (per-sample)
# ==============================================================================
# Loads pre-trained SC + SN models (from 03a_train_singler_ref.R) and classifies
# one FLEX sample. Fast and low-memory — models are small, training already done.
#
# Usage:
#   Rscript 03b_classify_singler.R \
#     --sample_id Pat1_P1 \
#     --h5_path /path/to/sample_filtered_feature_bc_matrix.h5 \
#     [--models_dir /path/to/singler_models/] \
#     [--output_dir /path/to/outputs/Pat1_P1/]
#
# Outputs:
#   {output_dir}/singler_labels.csv — cell_id, SC_label, SC_score,
#                                      SN_label, SN_score, consensus_label,
#                                      consensus_source
#
# Path conventions:
#   --models_dir defaults to ${CFG_PREPROCESSING_STEP_04}/models/singler_models/
#   --output_dir defaults to ${CFG_PREPROCESSING_STEP_04}/{sample_id}/
#   run/run_step_04_array.sh resolves both from CFG_PREPROCESSING_STEP_04.
#
# HPC: Tier 2 (16G, 30min per sample) — submitted as --array=0-61
# ==============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(SingleR)
  library(Seurat)
  library(Matrix)
  library(SummarizedExperiment)
})

step04_root <- Sys.getenv("CFG_PREPROCESSING_STEP_04", "")

parser <- ArgumentParser(description = "Step 04b: per-sample SingleR classification")
parser$add_argument("--sample_id", required = TRUE)
parser$add_argument("--h5_path", required = TRUE,
                    help = "CellRanger filtered H5 for this sample")
parser$add_argument("--models_dir", default = "",
                    help = "Dir containing singler_model_sc.rds and singler_model_sn.rds. Defaults to CFG_PREPROCESSING_STEP_04/models/singler_models/")
parser$add_argument("--output_dir", default = "",
                    help = "Per-sample output dir. Defaults to CFG_PREPROCESSING_STEP_04/{sample_id}/")
args <- parser$parse_args()

# Resolve models_dir
if (args$models_dir == "") {
  if (step04_root == "") stop("--models_dir is required (or set CFG_PREPROCESSING_STEP_04)", call. = FALSE)
  args$models_dir <- file.path(step04_root, "models", "singler_models")
}

# Resolve output_dir
if (args$output_dir == "") {
  if (step04_root == "") stop("--output_dir is required (or set CFG_PREPROCESSING_STEP_04)", call. = FALSE)
  args$output_dir <- file.path(step04_root, args$sample_id)
}

cat(sprintf("=== Step 04b: SingleR Classify: %s ===\n", args$sample_id))
cat(sprintf("  Models dir:  %s\n", args$models_dir))
cat(sprintf("  Output dir:  %s\n", args$output_dir))

stopifnot(file.exists(args$h5_path))
stopifnot(dir.exists(args$models_dir))
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# Load query (Ensembl IDs — models trained on Ensembl)
cat("Loading query H5...\n")
counts <- Read10X_h5(args$h5_path, use.names = FALSE)
seu <- CreateSeuratObject(counts = counts, project = args$sample_id,
                           min.cells = 0, min.features = 0)
new_names <- paste0(args$sample_id, "_", colnames(seu))
seu <- RenameCells(seu, new.names = new_names)
seu <- NormalizeData(seu, verbose = FALSE)
cat(sprintf("  %d cells, %d genes (Ensembl)\n", ncol(seu), nrow(seu)))
query_sce <- as.SingleCellExperiment(seu)
rm(seu, counts); gc()

# Load pre-trained models
cat("Loading SC model...\n")
sc_model <- readRDS(file.path(args$models_dir, "singler_model_sc.rds"))
cat("Loading SN model...\n")
sn_model <- readRDS(file.path(args$models_dir, "singler_model_sn.rds"))

# Pad query to reference gene space (SingleR 2.4.1 requires all reference genes)
pad_to_ref <- function(query_sce, trained) {
  ref_genes  <- rownames(trained$ref)
  test_genes <- rownames(query_sce)
  missing    <- setdiff(ref_genes, test_genes)
  cat(sprintf("  Test: %d, Ref: %d, Padding: %d (%.1f%%)\n",
              length(test_genes), length(ref_genes), length(missing),
              100 * length(missing) / length(ref_genes)))
  if (length(missing) == 0) return(query_sce)
  mat <- assay(query_sce, "logcounts")
  zero_rows <- sparseMatrix(i = integer(0), j = integer(0),
                            dims = c(length(missing), ncol(mat)),
                            dimnames = list(missing, colnames(mat)))
  padded <- padded <- rbind(mat, zero_rows)[ref_genes, ]
  SummarizedExperiment(assays = list(logcounts = padded))
}

cat("Classifying (SC)...\n")
sc_query <- pad_to_ref(query_sce, sc_model)
sc_pred  <- classifySingleR(test = sc_query, trained = sc_model)
cat(sprintf("  SC: %d unique labels\n", length(unique(sc_pred$labels))))
rm(sc_query); gc()

cat("Classifying (SN)...\n")
sn_query <- pad_to_ref(query_sce, sn_model)
sn_pred  <- classifySingleR(test = sn_query, trained = sn_model)
cat(sprintf("  SN: %d unique labels\n", length(unique(sn_pred$labels))))
rm(sn_query, query_sce, sc_model, sn_model); gc()

# Confidence: delta.next (gap between best and second-best score)
sc_conf <- as.numeric(if (!is.null(sc_pred$delta.next)) sc_pred$delta.next else rep(NA_real_, nrow(sc_pred)))
sn_conf <- as.numeric(if (!is.null(sn_pred$delta.next)) sn_pred$delta.next else rep(NA_real_, nrow(sn_pred)))
sc_conf[is.na(sc_conf)] <- 0
sn_conf[is.na(sn_conf)] <- 0

# Consensus: use reference with higher delta.next (more confident)
consensus_label  <- ifelse(sc_conf >= sn_conf, sc_pred$labels, sn_pred$labels)
consensus_source <- ifelse(sc_conf >= sn_conf, "SC", "SN")

out_df <- data.frame(
  cell_id          = rownames(sc_pred),
  SC_label         = sc_pred$labels,
  SC_score         = round(sc_conf, 4),
  SN_label         = sn_pred$labels,
  SN_score         = round(sn_conf, 4),
  consensus_label  = consensus_label,
  consensus_source = consensus_source,
  stringsAsFactors = FALSE
)

out_path <- file.path(args$output_dir, "singler_labels.csv")
write.csv(out_df, out_path, row.names = FALSE)

cat(sprintf("\n  %d cells — SC: %d types, SN: %d types\n",
            nrow(out_df), length(unique(out_df$SC_label)), length(unique(out_df$SN_label))))
cat(sprintf("  Consensus: %d SC (%.0f%%), %d SN (%.0f%%)\n",
            sum(out_df$consensus_source == "SC"), 100 * mean(out_df$consensus_source == "SC"),
            sum(out_df$consensus_source == "SN"), 100 * mean(out_df$consensus_source == "SN")))
cat(sprintf("  Saved: %s\n", out_path))
