#!/usr/bin/env Rscript
# AUCell scoring for multi-panel (tumor + preneoplastic) iCAF signatures
suppressPackageStartupMessages({ library(Matrix); library(AUCell) })

args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(name) {
  i <- which(args == paste0("--", name))
  if (length(i) == 0) stop(sprintf("missing --%s", name))
  args[i + 1]
}
mtx_path    <- parse_arg("mtx")
cells_path  <- parse_arg("cells")
genes_path  <- parse_arg("genes")
panels_path <- parse_arg("panels")
out_path    <- parse_arg("out")

cat("[R] reading MTX:", mtx_path, "\n")
mat <- readMM(mtx_path)
mat <- as(mat, "CsparseMatrix")
cell_ids <- readLines(cells_path)
gene_syms <- readLines(genes_path)
rownames(mat) <- gene_syms
colnames(mat) <- cell_ids
cat("[R] mat:", nrow(mat), "genes x", ncol(mat), "cells\n")

cat("[R] AUCell_buildRankings\n")
set.seed(42)
rankings <- AUCell_buildRankings(mat, plotStats = FALSE, verbose = FALSE)

panels_df <- read.table(panels_path, header = TRUE, sep = "\t",
                        stringsAsFactors = FALSE)
panel_names <- unique(panels_df$panel)
gene_sets <- lapply(panel_names, function(nm) {
  unique(panels_df$gene[panels_df$panel == nm])
})
names(gene_sets) <- panel_names

for (nm in panel_names) {
  found <- intersect(gene_sets[[nm]], rownames(mat))
  cat(sprintf("[R] %s: %d/%d genes present\n", nm, length(found),
              length(gene_sets[[nm]])))
}

cat("[R] AUCell_calcAUC\n")
auc <- AUCell_calcAUC(geneSets = gene_sets, rankings = rankings,
                     aucMaxRank = ceiling(0.05 * nrow(rankings)),
                     verbose = FALSE)
auc_mat <- AUCell::getAUC(auc)
cat("[R] auc_mat dim:", dim(auc_mat)[1], "x", dim(auc_mat)[2], "\n")

out <- data.frame(cell_id = colnames(auc_mat), stringsAsFactors = FALSE)
for (nm in rownames(auc_mat)) {
  out[[nm]] <- as.numeric(auc_mat[nm, ])
}
write.csv(out, out_path, row.names = FALSE)
cat("[R] wrote:", out_path, "\n")
