#!/usr/bin/env Rscript
# 03b_limma_markers.R — pseudobulk limma-voom one-vs-rest markers per cluster.
#
# FLEX is single-study (Spatial HBCA), so the design is:
#   ~ library_id + condition
# with duplicateCorrelation blocked on patient_id when n_patients >= 3.
# library_id serves as the batch covariate (analogous to study covariate in
# multi-study designs); the `condition` factor is the cluster contrast for the
# one-vs-rest comparison.
#
# Usage:
#   Rscript 03b_limma_markers.R \
#     --counts /path/{Compartment}_pseudobulk_counts.csv \
#     --meta   /path/{Compartment}_pseudobulk_meta.csv \
#     --output-dir /path/limma/ \
#     --level-tag {Compartment}_leiden_1.0 \
#     --n-top 50

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(argparse)
})

parser <- ArgumentParser(description = "FLEX pseudobulk limma-voom markers")
parser$add_argument("--counts",     required = TRUE, help = "Pseudobulk counts CSV (rows=samples, cols=genes)")
parser$add_argument("--meta",       required = TRUE, help = "Pseudobulk meta CSV (sample_id, cluster, library_id, patient_id, n_cells)")
parser$add_argument("--output-dir", required = TRUE)
parser$add_argument("--level-tag",  required = TRUE, help = "Output filename tag, e.g. Epithelial_leiden_1.0")
parser$add_argument("--n-top",      type = "integer", default = 50, help = "Top N markers per cluster (default 50)")
parser$add_argument("--hvg-n",      type = "integer", default = 2000, help = "Variance filter HVGs (default 2000)")
parser$add_argument("--min-gene-count", type = "integer", default = 10)
parser$add_argument("--padj-threshold", type = "double",  default = 0.01)
parser$add_argument("--logfc-threshold", type = "double", default = 1.0)
parser$add_argument("--no-dupcor",  action = "store_true",
                    help = "Skip duplicateCorrelation even if n_patients >= 3")
args <- parser$parse_args()

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load pseudobulk counts ---
cat("Loading pseudobulk counts:", args$counts, "\n")
counts_df <- read.csv(args$counts, row.names = 1, check.names = FALSE)
cat(sprintf("  %d samples x %d genes\n", nrow(counts_df), ncol(counts_df)))

cat("Loading meta:", args$meta, "\n")
meta <- read.csv(args$meta, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

required <- c("sample_id", "cluster", "library_id", "patient_id", "n_cells")
missing <- setdiff(required, colnames(meta))
if (length(missing) > 0) {
  stop(sprintf("meta CSV missing required columns: %s",
               paste(missing, collapse = ", ")))
}

common <- intersect(rownames(counts_df), meta$sample_id)
if (length(common) == 0) stop("No overlap between counts and meta sample_ids")
counts_df <- counts_df[common, , drop = FALSE]
meta <- meta[common, , drop = FALSE]

n_clusters  <- length(unique(meta$cluster))
n_libraries <- length(unique(meta$library_id))
n_patients  <- length(unique(meta$patient_id))
cat(sprintf("  %d clusters, %d libraries, %d patients\n",
            n_clusters, n_libraries, n_patients))

use_dupcor <- (!args$no_dupcor) && (n_patients >= 3)
if (use_dupcor) {
  cat("  Design: ~ library_id + condition; dupCor block = patient_id\n")
} else {
  cat(sprintf("  Design: ~ library_id + condition; dupCor disabled (n_patients=%d, no_dupcor=%s)\n",
              n_patients, args$no_dupcor))
}

# Transpose to genes × samples
count_mat <- t(as.matrix(counts_df))

# Low-expression filter
gene_sums <- rowSums(count_mat)
keep <- gene_sums >= args$min_gene_count
count_mat <- count_mat[keep, ]
cat(sprintf("  Genes after low-expr filter (>=%d total counts): %d / %d\n",
            args$min_gene_count, sum(keep), length(keep)))

# Variance filter
if (nrow(count_mat) > args$hvg_n) {
  gene_vars <- apply(count_mat, 1, var)
  hvg_idx <- order(gene_vars, decreasing = TRUE)[seq_len(args$hvg_n)]
  count_mat <- count_mat[hvg_idx, ]
  cat(sprintf("  Genes after variance filter (top %d): %d\n",
              args$hvg_n, nrow(count_mat)))
}

# --- One-vs-rest per cluster ---
clusters <- sort(unique(meta$cluster))
cat(sprintf("Running limma-voom for %d clusters...\n", length(clusters)))

all_results <- list()

for (cl in clusters) {
  t0 <- Sys.time()
  cat(sprintf("  cluster %s ...", cl))

  meta$condition <- factor(ifelse(meta$cluster == cl, "target", "rest"),
                           levels = c("rest", "target"))

  # Drop libraries with only one pseudobulk sample (singleton library factor)
  use_meta <- meta
  lib_tab <- table(use_meta$library_id)
  singleton_libs <- names(lib_tab[lib_tab < 2])
  if (length(singleton_libs) > 0) {
    use_meta <- use_meta[!use_meta$library_id %in% singleton_libs, ]
  }
  use_meta$library_id <- factor(use_meta$library_id)

  n_target <- sum(use_meta$condition == "target", na.rm = TRUE)
  n_rest   <- sum(use_meta$condition == "rest",   na.rm = TRUE)
  if (n_target < 3 || n_rest < 3) {
    cat(" skipped (target/rest < 3 samples)\n")
    next
  }

  tryCatch({
    design <- model.matrix(~ library_id + condition, data = use_meta)
    dge    <- DGEList(counts = count_mat[, rownames(use_meta)])
    dge    <- calcNormFactors(dge)

    consensus_cor <- NA_real_
    if (use_dupcor) {
      v0     <- voom(dge, design, plot = FALSE)
      corfit <- duplicateCorrelation(v0, design, block = use_meta$patient_id)
      consensus_cor <- corfit$consensus.correlation
      v <- voom(dge, design, plot = FALSE,
                block = use_meta$patient_id, correlation = consensus_cor)
      cat(sprintf(" [voom+dupCor cor=%.2f]", consensus_cor))
    } else {
      v <- voom(dge, design, plot = FALSE)
      cat(" [voom]")
    }
    flush.console()

    if (!is.na(consensus_cor)) {
      fit <- lmFit(v, design, block = use_meta$patient_id, correlation = consensus_cor)
    } else {
      fit <- lmFit(v, design)
    }
    fit <- eBayes(fit)

    res <- topTable(fit, coef = "conditiontarget", number = Inf, sort.by = "P")
    res_df <- as.data.frame(res)

    # DESeq2-compatible column names (downstream tools expect these)
    names(res_df)[names(res_df) == "logFC"]     <- "log2FoldChange"
    names(res_df)[names(res_df) == "AveExpr"]   <- "baseMean"
    names(res_df)[names(res_df) == "t"]         <- "stat"
    names(res_df)[names(res_df) == "P.Value"]   <- "pvalue"
    names(res_df)[names(res_df) == "adj.P.Val"] <- "padj"

    res_df$gene    <- rownames(res_df)
    res_df$cluster <- cl

    n_sig <- sum(res_df$padj < args$padj_threshold &
                 res_df$log2FoldChange > args$logfc_threshold, na.rm = TRUE)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf(" %d sig markers (padj<%.2g & lfc>%.1f) [%.1fs]\n",
                n_sig, args$padj_threshold, args$logfc_threshold, elapsed))

    all_results[[as.character(cl)]] <- res_df
  }, error = function(e) {
    cat(sprintf(" ERROR: %s\n", conditionMessage(e)))
  })
  gc()
}

if (length(all_results) == 0) {
  stop("No clusters produced limma results")
}

full_df <- do.call(rbind, all_results)
rownames(full_df) <- NULL
full_path <- file.path(args$output_dir,
                       paste0("limma_markers_", args$level_tag, ".csv"))
write.csv(full_df, full_path, row.names = FALSE)
cat(sprintf("\nFull results: %s (%d rows)\n", full_path, nrow(full_df)))

# Top-N per cluster (filtered by significance + LFC, sorted by padj asc, lfc desc)
top_df <- do.call(rbind, lapply(all_results, function(df) {
  df <- df[!is.na(df$padj) &
           df$padj < args$padj_threshold &
           df$log2FoldChange > args$logfc_threshold, ]
  df <- df[order(df$padj, -df$log2FoldChange), ]
  head(df, args$n_top)
}))
rownames(top_df) <- NULL
top_path <- file.path(args$output_dir,
                      paste0("limma_top_markers_", args$level_tag, ".csv"))
write.csv(top_df, top_path, row.names = FALSE)
cat(sprintf("Top %d markers: %s (%d rows)\n",
            args$n_top, top_path, nrow(top_df)))

cat("\n=== Top 5 markers per cluster ===\n")
for (cl in sort(unique(top_df$cluster))) {
  sub   <- top_df[top_df$cluster == cl, ]
  genes <- head(sub$gene, 5)
  lfcs  <- head(round(sub$log2FoldChange, 1), 5)
  cat(sprintf("  cluster %s: %s\n", cl,
              paste(sprintf("%s(lfc=%s)", genes, lfcs), collapse = ", ")))
}

cat("\nDone.\n")
