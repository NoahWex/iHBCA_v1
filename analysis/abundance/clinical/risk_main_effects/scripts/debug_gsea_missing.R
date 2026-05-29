#!/usr/bin/env Rscript
# debug_gsea_missing.R
# Run fgsea against the 7 BRCA1 marker files that are missing from
# per_nhoodgroup_nes.csv to diagnose Mode 1 (F.3 skip stub) vs Mode 2
# (fgsea silent failure).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fgsea); library(msigdbr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "debug_gsea_missing.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

f3_dir <- file.path(paths$outputs$stageF3, "parity_x_HR_BRCA1")
missing_groups <- c("imm__CD8_Trm_1", "imm__Th17_20", "str__Fibro_SFRP4_48",
                     "str__Lym_major_80", "str__Pericyte_active_26",
                     "str__Vas_arterial_81", "str__VSMC_40")

cat("Loading hallmark...\n")
hallmark <- split(msigdbr::msigdbr(species = "Homo sapiens", category = "H")$ensembl_gene,
                   msigdbr::msigdbr(species = "Homo sapiens", category = "H")$gs_name)
cat(sprintf("  %d sets\n\n", length(hallmark)))

for (g in missing_groups) {
  csv <- file.path(f3_dir, paste0(g, "_vs_parent.csv"))
  cat(sprintf("====================\nGROUP: %s\n", g))
  d <- tryCatch(read_csv(csv, show_col_types = FALSE),
                 error = function(e) {cat(sprintf("  read error: %s\n", conditionMessage(e))); NULL})
  if (is.null(d)) next
  cat(sprintf("  rows: %d, cols: %s\n", nrow(d), paste(colnames(d), collapse=",")))
  if ("skipped" %in% colnames(d)) {
    cat(sprintf("  >>> MODE 1: F.3 skip stub. reason=%s\n", paste(d$reason, collapse=";")))
    next
  }
  if (!all(c("gene_id", "logFC") %in% colnames(d))) {
    cat("  missing gene_id/logFC\n"); next
  }
  pval_col <- intersect(c("P.Value", "PValue"), colnames(d))[1]
  if (is.na(pval_col)) { cat("  no pval col!\n"); next }
  rank_vec <- sign(d$logFC) * (-log10(pmax(d[[pval_col]], 1e-300)))
  names(rank_vec) <- d$gene_id
  rank_vec <- rank_vec[!is.na(rank_vec) & is.finite(rank_vec)]
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  cat(sprintf("  rank_vec n=%d, head=%s, tail=%s\n", length(rank_vec),
              paste(round(head(rank_vec, 3), 2), collapse=","),
              paste(round(tail(rank_vec, 3), 2), collapse=",")))
  n_dup <- sum(duplicated(names(rank_vec)))
  cat(sprintf("  duplicate gene_ids in rank_vec: %d\n", n_dup))
  if (length(rank_vec) < 100) {
    cat("  >>> MODE A: <100 ranks (script floor)\n"); next
  }
  res <- tryCatch(
    fgsea::fgseaMultilevel(pathways = hallmark, stats = rank_vec,
                            minSize = 10, maxSize = 500,
                            nPermSimple = 1000, eps = 0, nproc = 1),
    error = function(e) {cat(sprintf("  >>> MODE 2A: fgsea error: %s\n", conditionMessage(e))); NULL},
    warning = function(w) {cat(sprintf("  >>> WARN: %s\n", conditionMessage(w))); NULL}
  )
  if (is.null(res)) next
  cat(sprintf("  fgsea returned %d rows\n", nrow(res)))
  if (nrow(res) == 0) {
    cat("  >>> MODE 2B: fgsea returned empty\n"); next
  }
  cat(sprintf("  top 3 by NES: %s\n",
              paste(head(res[order(-abs(res$NES)), c("pathway", "NES", "padj")], 3) %>%
                      apply(1, function(x) sprintf("%s NES=%.2f padj=%.2g",
                                                      x[1], as.numeric(x[2]), as.numeric(x[3]))),
                    collapse=" | ")))
}
cat("\n=== debug done ===\n")
