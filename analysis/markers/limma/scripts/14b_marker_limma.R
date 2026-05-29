#!/usr/bin/env Rscript
# Pseudobulk limma-voom marker gene analysis for iHBCA cell types
#
# Replaces the DESeq2 implementation with limma-voom for substantially faster
# one-vs-rest DE across 7 Leiden resolutions. Design: ~ study + condition
# accounts for batch effects across the 7 source studies.
#
# Usage (single run):
#   Rscript 14b_marker_limma.R \
#     --counts /path/to/pseudobulk_counts.csv \
#     --meta /path/to/pseudobulk_meta.csv \
#     --output-dir /path/to/output \
#     --level-tag leiden_0.1 \
#     --n-top 20
#
# Usage (chunked SLURM array mode):
#   Rscript 14b_marker_limma.R ... --chunk-idx 0 --chunk-size 8
#   Outputs: limma_markers_{level}_c0.csv  (partial; merge with merge_14b_limma_chunks.R)
#
# Inputs:
#   *_pseudobulk_counts.csv — rows = samples (cell_type__donor), cols = genes
#   *_pseudobulk_meta.csv   — rows = samples, cols = cell_type, donor, study, n_cells
#
# Outputs (unchunked):
#   limma_markers_{level}.csv     — full results (DESeq2-compatible column names)
#   limma_top_markers_{level}.csv — top N per cell type
# Outputs (chunked):
#   limma_markers_{level}_c{N}.csv     — partial results for chunk N

suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(argparse)
})

parser <- ArgumentParser(description = "Pseudobulk limma-voom markers")
parser$add_argument("--counts", required = TRUE, help = "Pseudobulk counts CSV")
parser$add_argument("--meta", default = NULL,
                    help = "Pseudobulk metadata CSV (with study column)")
parser$add_argument("--output-dir", required = TRUE, help = "Output directory")
parser$add_argument("--level-tag", default = "level1",
                    help = "Output file tag (default: level1)")
parser$add_argument("--n-top", type = "integer", default = 20,
                    help = "Top N markers per type")
parser$add_argument("--workers", type = "integer", default = 1,
                    help = "Deprecated; retained for CLI compat.")
parser$add_argument("--n-cells", type = "integer", default = NULL,
                    help = "Subsample to N pseudobulk samples (smoke-test flag)")
parser$add_argument("--chunk-idx", type = "integer", default = NULL,
                    help = "0-based chunk index for SLURM array mode (use with --chunk-size)")
parser$add_argument("--chunk-size", type = "integer", default = 8L,
                    help = "Cell types per chunk in array mode (default: 8)")
args <- parser$parse_args()

output_dir <- args$output_dir
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load pseudobulk counts ---
cat("Loading pseudobulk counts...\n")
counts_df <- read.csv(args$counts, row.names = 1, check.names = FALSE)
cat(sprintf("  %d samples x %d genes\n", nrow(counts_df), ncol(counts_df)))

# Load metadata: prefer external CSV, fall back to parsing row names
if (!is.null(args$meta) && file.exists(args$meta)) {
  cat("Loading metadata from:", args$meta, "\n")
  meta <- read.csv(args$meta, row.names = 1, stringsAsFactors = FALSE)
  common <- intersect(rownames(counts_df), rownames(meta))
  counts_df <- counts_df[common, ]
  meta <- meta[common, ]
} else {
  cat("Parsing metadata from row names (no --meta provided)\n")
  sample_ids <- rownames(counts_df)
  meta <- data.frame(
    cell_type = sub("__.*", "", sample_ids),
    donor     = sub(".*__", "", sample_ids),
    stringsAsFactors = FALSE
  )
  rownames(meta) <- sample_ids
}

has_study <- "study" %in% colnames(meta)
cat(sprintf("  %d cell types, %d donors\n",
            length(unique(meta$cell_type)), length(unique(meta$donor))))
if (has_study) {
  cat(sprintf("  %d studies: %s\n", length(unique(meta$study)),
              paste(sort(unique(meta$study)), collapse = ", ")))
  cat("  Design: ~ study + condition\n")
} else {
  cat("  WARNING: No study column -- design: ~ condition (no batch correction)\n")
}

# Optional: subsample pseudobulk samples for smoke-testing
if (!is.null(args$n_cells) && args$n_cells < nrow(counts_df)) {
  set.seed(42)
  idx <- sample(nrow(counts_df), args$n_cells)
  counts_df <- counts_df[idx, ]
  meta      <- meta[idx, ]
  cat(sprintf("  [--n-cells] Subsampled to %d pseudobulk samples\n", nrow(counts_df)))
}

# Transpose: limma expects genes x samples
count_mat <- t(as.matrix(counts_df))

# Filter low-expressed genes
gene_sums <- rowSums(count_mat)
keep <- gene_sums >= 10
count_mat <- count_mat[keep, ]
cat(sprintf("  Genes after low-expr filter: %d / %d\n", sum(keep), length(keep)))

# Variance filter: top 2K genes matches the scANVI integration gene space and
# keeps dupCor runtime tractable (~5x faster than 10K).
n_hvg <- 2000L
if (nrow(count_mat) > n_hvg) {
  gene_vars <- apply(count_mat, 1, var)
  hvg_idx   <- order(gene_vars, decreasing = TRUE)[seq_len(n_hvg)]
  count_mat <- count_mat[hvg_idx, ]
  cat(sprintf("  Genes after variance filter (top %d): %d\n", n_hvg, nrow(count_mat)))
} else {
  cat(sprintf("  Variance filter not applied (%d genes <= %d threshold)\n",
              nrow(count_mat), n_hvg))
}

# --- Run limma-voom one-vs-rest for each cell type ---
cell_types_all <- sort(unique(meta$cell_type))

# Chunk mode: process a fixed window of cell types for SLURM array parallelism.
# Array task_id decodes to chunk_idx externally; this script just processes its slice.
if (!is.null(args$chunk_idx)) {
  start_1 <- args$chunk_idx * args$chunk_size + 1L  # 1-indexed
  end_1   <- min((args$chunk_idx + 1L) * args$chunk_size, length(cell_types_all))
  if (start_1 > length(cell_types_all)) {
    cat(sprintf("\nChunk %d out of range (%d types total, chunk_size=%d) — nothing to do.\n",
                args$chunk_idx, length(cell_types_all), args$chunk_size))
    quit(status = 0)
  }
  cell_types <- cell_types_all[start_1:end_1]
  cat(sprintf("\nChunk %d/%d: types %d-%d of %d total\n",
              args$chunk_idx,
              ceiling(length(cell_types_all) / args$chunk_size) - 1L,
              start_1, end_1, length(cell_types_all)))
} else {
  cell_types <- cell_types_all
}

cat(sprintf("Running limma-voom for %d cell types...\n", length(cell_types)))

all_results <- list()

for (ct in cell_types) {
  t0 <- Sys.time()
  cat(sprintf("  %s...", ct))

  meta$condition <- factor(ifelse(meta$cell_type == ct, "target", "rest"),
                           levels = c("rest", "target"))
  use_meta <- meta

  if (has_study) {
    study_tab        <- table(use_meta$study)
    singleton_studies <- names(study_tab[study_tab < 2])
    if (length(singleton_studies) > 0) {
      use_meta <- use_meta[!use_meta$study %in% singleton_studies, ]
    }
    use_meta$study <- factor(use_meta$study)
  }

  n_target <- sum(use_meta$condition == "target", na.rm = TRUE)
  n_rest   <- sum(use_meta$condition == "rest",   na.rm = TRUE)
  if (n_target < 3 || n_rest < 3) {
    cat(" skipped (< 3 samples)\n")
    next
  }

  tryCatch({
    design_formula <- if (has_study) ~ study + condition else ~ condition
    design <- model.matrix(design_formula, data = use_meta)

    dge <- DGEList(counts = count_mat[, rownames(use_meta)])
    dge <- calcNormFactors(dge)

    # Iterative voom + duplicateCorrelation: models within-donor correlation
    # across pseudobulk samples (same donor in both target and rest groups).
    # Falls back to simple voom if fewer than 3 unique donors.
    n_donors <- length(unique(use_meta$donor))
    step_t <- Sys.time()
    if (n_donors >= 3) {
      v0          <- voom(dge, design, plot = FALSE)
      corfit      <- duplicateCorrelation(v0, design, block = use_meta$donor)
      consensus_cor <- corfit$consensus.correlation
      v <- voom(dge, design, plot = FALSE,
                block = use_meta$donor, correlation = consensus_cor)
      cat(sprintf(" [voom+dupCor %.1fs cor=%.2f]",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs")),
                  consensus_cor))
    } else {
      v <- voom(dge, design, plot = FALSE)
      consensus_cor <- NA
      cat(sprintf(" [voom %.1fs (no dupCor, %d donors)]",
                  as.numeric(difftime(Sys.time(), step_t, units = "secs")),
                  n_donors))
    }
    flush.console()

    step_t <- Sys.time()
    if (!is.na(consensus_cor)) {
      fit <- lmFit(v, design, block = use_meta$donor, correlation = consensus_cor)
    } else {
      fit <- lmFit(v, design)
    }
    fit <- eBayes(fit)
    cat(sprintf(" [eBayes %.1fs]",
                as.numeric(difftime(Sys.time(), step_t, units = "secs"))))
    flush.console()

    res    <- topTable(fit, coef = "conditiontarget", number = Inf, sort.by = "P")
    res_df <- as.data.frame(res)

    # Rename to DESeq2-compatible column names for downstream consumers
    # (marker_summary.py, branch_report.py, annotation_review.py all expect these)
    names(res_df)[names(res_df) == "logFC"]     <- "log2FoldChange"
    names(res_df)[names(res_df) == "AveExpr"]   <- "baseMean"
    names(res_df)[names(res_df) == "t"]         <- "stat"
    names(res_df)[names(res_df) == "P.Value"]   <- "pvalue"
    names(res_df)[names(res_df) == "adj.P.Val"] <- "padj"

    res_df$gene      <- rownames(res_df)
    res_df$cell_type <- ct

    n_sig   <- sum(res_df$padj < 0.05, na.rm = TRUE)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cat(sprintf(" %d significant (padj<0.05)  [%.1fs]\n", n_sig, elapsed))
    flush.console()

    all_results[[as.character(ct)]] <- res_df
  }, error = function(e) {
    cat(sprintf(" ERROR: %s\n", conditionMessage(e)))
    flush.console()
  })
  gc()
}

# --- Combine and save ---
if (length(all_results) == 0) {
  cat("No results produced.\n")
  quit(status = 1)
}

full_df      <- do.call(rbind, all_results)
rownames(full_df) <- NULL

# In chunk mode, write a partial file; merge_14b_limma_chunks.R combines them.
chunk_suffix <- if (!is.null(args$chunk_idx)) sprintf("_c%d", args$chunk_idx) else ""
full_path <- file.path(output_dir,
  paste0("limma_markers_", args$level_tag, chunk_suffix, ".csv"))
write.csv(full_df, full_path, row.names = FALSE)
cat(sprintf("\nFull results: %s (%d rows)\n", full_path, nrow(full_df)))

# Top markers by padj (written only when not chunked; merge script produces top file)
if (is.null(args$chunk_idx)) {
  top_df <- do.call(rbind, lapply(all_results, function(df) {
    df <- df[!is.na(df$padj) & df$padj < 0.05, ]
    df <- df[order(df$padj, -abs(df$log2FoldChange)), ]
    head(df, args$n_top)
  }))
  rownames(top_df) <- NULL

  top_path <- file.path(output_dir,
    paste0("limma_top_markers_", args$level_tag, ".csv"))
  write.csv(top_df, top_path, row.names = FALSE)
  cat(sprintf("Top %d markers: %s (%d rows)\n", args$n_top, top_path, nrow(top_df)))

  cat("\n=== Top 5 markers per cell type ===\n")
  for (ct in sort(unique(top_df$cell_type))) {
    sub   <- top_df[top_df$cell_type == ct, ]
    genes <- head(sub$gene, 5)
    lfcs  <- head(round(sub$log2FoldChange, 1), 5)
    cat(sprintf("  %s: %s\n", ct,
                paste(sprintf("%s(lfc=%s)", genes, lfcs), collapse = ", ")))
  }
} else {
  cat(sprintf("  Chunk mode: run merge_14b_limma_chunks.R to combine chunks.\n"))
}

cat("\nDone.\n")
