#!/usr/bin/env Rscript
# limma-voom pseudobulk DE for FLEX L2S markers — ANNOTATED vocabulary.
# Uses canonical flex_l2s_labels.csv (l2s_label column) instead of raw cluster IDs.
# Excludes ARTIFACT_* groups. Otherwise identical to prior version.

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(data.table)
})

parser <- ArgumentParser()
parser$add_argument("--seurat-path", required = TRUE)
parser$add_argument("--l2s-csv", required = TRUE,
                    help = "flex_l2s_labels.csv with annotated l2s_label column")
parser$add_argument("--donor-col", default = "orig.ident")
parser$add_argument("--min-cells-per-pb", type = "integer", default = 10)
parser$add_argument("--min-donors-per-L2S", type = "integer", default = 3)
parser$add_argument("--outdir", required = TRUE)
args <- parser$parse_args()

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(args$outdir, "flex_l2s_limma_de_annotated.csv")

cat("=== FLEX L2S limma-voom pseudobulk DE (annotated vocabulary) ===\n")

cat("[1/5] loading annotated L2S labels\n")
l2s <- fread(args$l2s_csv)
cat("      rows:", nrow(l2s), "  unique l2s_label:", uniqueN(l2s$l2s_label), "\n")

# Exclude artifacts
artifact_mask <- grepl("^ARTIFACT_", l2s$l2s_label)
cat("      excluding", sum(artifact_mask), "artifact cells (",
    paste(unique(l2s$l2s_label[artifact_mask]), collapse = ", "), ")\n")
l2s <- l2s[!artifact_mask]
cat("      after artifact exclusion:", nrow(l2s), "cells,",
    uniqueN(l2s$l2s_label), "L2S types\n")

cat("[2/5] loading Seurat object\n")
srt <- readRDS(args$seurat_path)
srt <- UpdateSeuratObject(srt)
cat("      cells x genes:", ncol(srt), "x", nrow(srt), "\n")

if (!(args$donor_col %in% colnames(srt@meta.data))) {
  stop("--donor-col not in metadata. Available: ",
       paste(colnames(srt@meta.data), collapse = ", "))
}

idx <- match(colnames(srt), l2s$cell_id)
srt$l2s_label <- l2s$l2s_label[idx]
srt$donor_id <- as.character(srt@meta.data[[args$donor_col]])
srt <- subset(srt, cells = colnames(srt)[!is.na(srt$l2s_label)])
cat("      after L2S join + artifact exclusion:", ncol(srt), "cells\n")
cat("      donors:", uniqueN(srt$donor_id),
    "  L2S types:", uniqueN(srt$l2s_label), "\n")
cat("      L2S distribution:\n")
print(sort(table(srt$l2s_label), decreasing = TRUE))
cat("\n")

cat("[3/5] building pseudobulk (donor x L2S)\n")
counts <- GetAssayData(srt, layer = "counts")
meta <- data.table(
  cell_id = colnames(srt),
  donor   = srt$donor_id,
  l2s     = srt$l2s_label
)
meta[, pb_key := paste(donor, l2s, sep = "__")]

n_cells_per_pb <- meta[, .N, by = pb_key]
keep_pb <- n_cells_per_pb[N >= args$min_cells_per_pb, pb_key]
cat("      pseudobulks passing", args$min_cells_per_pb, "cell filter:",
    length(keep_pb), "of", uniqueN(meta$pb_key), "\n")

pb_counts <- matrix(0, nrow = nrow(counts), ncol = length(keep_pb),
                    dimnames = list(rownames(counts), keep_pb))
for (k in keep_pb) {
  cell_idx <- which(meta$pb_key == k)
  if (length(cell_idx) == 1) {
    pb_counts[, k] <- as.numeric(counts[, cell_idx])
  } else {
    pb_counts[, k] <- rowSums(counts[, cell_idx, drop = FALSE])
  }
}

pb_meta <- meta[pb_key %in% keep_pb, .(n_cells = .N,
                                        donor = first(donor),
                                        l2s = first(l2s)),
                by = pb_key]
setkey(pb_meta, pb_key)
pb_meta <- pb_meta[colnames(pb_counts)]

donors_per_l2s <- pb_meta[, .(n_donors = uniqueN(donor)), by = l2s]
keep_l2s <- donors_per_l2s[n_donors >= args$min_donors_per_L2S, l2s]
dropped <- setdiff(unique(pb_meta$l2s), keep_l2s)
cat("      L2S passing", args$min_donors_per_L2S, "donor filter:",
    length(keep_l2s), " (dropped:", length(dropped),
    if (length(dropped) > 0) paste0(" -- ", paste(dropped, collapse = ", ")) else "",
    ")\n")

pb_meta_kept <- pb_meta[l2s %in% keep_l2s]
pb_counts <- pb_counts[, pb_meta_kept$pb_key]
cat("      pseudobulk matrix:", nrow(pb_counts), "x", ncol(pb_counts), "\n\n")

cat("[4/5] limma-voom one-vs-rest per L2S\n")
dge <- DGEList(counts = pb_counts)
keep_genes <- filterByExpr(dge, group = pb_meta_kept$l2s)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat("      genes after filterByExpr:", nrow(dge), "\n")

all_rows <- list()
for (l2s_type in keep_l2s) {
  is_target <- factor(ifelse(pb_meta_kept$l2s == l2s_type, "target", "rest"),
                      levels = c("rest", "target"))
  if (any(table(is_target) < 2)) {
    cat("      skip:", l2s_type, "(degenerate)\n")
    next
  }
  design <- model.matrix(~ is_target)
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = "is_targettarget", number = Inf, sort.by = "none")
  tt$gene_symbol <- rownames(tt)
  tt$cell_type <- l2s_type
  tt$n_pseudobulks_in <- sum(is_target == "target")
  tt$n_pseudobulks_out <- sum(is_target == "rest")
  tt$n_cells_in <- sum(pb_meta_kept[l2s == l2s_type, n_cells])
  all_rows[[l2s_type]] <- as.data.table(tt)
}

res <- rbindlist(all_rows, use.names = TRUE, fill = TRUE)
cat("      total DE rows:", nrow(res), "\n")

cat("[5/5] writing output\n")
setnames(res,
         old = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
         new = c("log2FC", "avg_log_cpm", "t_statistic", "p_val", "FDR", "B_statistic"))
res[, label_level := "L2S"]
res <- res[, .(gene_symbol, cell_type, label_level,
               log2FC, avg_log_cpm, t_statistic,
               p_val, FDR, B_statistic,
               n_pseudobulks_in, n_pseudobulks_out, n_cells_in)]
setorder(res, cell_type, FDR, p_val)
fwrite(res, out_csv)
cat("      wrote:", out_csv, "  (", nrow(res), "rows x", ncol(res), "cols)\n")
cat("=== Complete ===\n")
