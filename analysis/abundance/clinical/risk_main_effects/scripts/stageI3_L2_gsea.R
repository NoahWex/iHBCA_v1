#!/usr/bin/env Rscript
# stageI3_L2_gsea.R
# fgsea preranked over per-(contrast × parent_L2_joint) limma-voom marker tables
# from stageF3_L2_markers.R. Same pathway DBs as Stage I.3 NhoodGroup version.
#
# Output: outputs/stageI3_gsea/per_L2_nes.csv

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fgsea); library(msigdbr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI3_L2_gsea.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

out_dir <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea")
ensure_dir(out_dir)

# Pathway DBs
cat("Loading MSigDB Hallmark + Reactome...\n")
hallmark <- {
  df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  split(df$ensembl_gene, df$gs_name)
}
reactome <- {
  df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2",
                          subcategory = "CP:REACTOME")
  split(df$ensembl_gene, df$gs_name)
}
collections <- list(hallmark = hallmark, reactome = reactome)
cat(sprintf("  hallmark=%d, reactome=%d\n", length(hallmark), length(reactome)))

l2_root <- file.path(paths$inquiry_root, "outputs", "stageF3_L2_markers")
contrast_dirs <- list.dirs(l2_root, recursive = FALSE)

run_one <- function(csv, contrast, L2) {
  d <- tryCatch(read_csv(csv, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if ("skipped" %in% colnames(d)) return(NULL)
  if (!all(c("gene_id", "logFC") %in% colnames(d))) return(NULL)
  pval_col <- intersect(c("P.Value", "PValue", "p_value", "pvalue"), colnames(d))[1]
  rank_vec <- if (is.na(pval_col)) d$logFC
                else sign(d$logFC) * (-log10(pmax(d[[pval_col]], 1e-300)))
  names(rank_vec) <- d$gene_id
  rank_vec <- rank_vec[!is.na(rank_vec) & is.finite(rank_vec)]
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  if (length(rank_vec) < 100) return(NULL)
  out <- list()
  for (cn in names(collections)) {
    res <- tryCatch(
      fgsea::fgseaMultilevel(pathways = collections[[cn]], stats = rank_vec,
                              minSize = 10, maxSize = 500, nPermSimple = 1000,
                              eps = 0, nproc = 1),
      error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    res <- as.data.frame(res)
    res$leadingEdge <- vapply(res$leadingEdge,
                                function(x) paste(x, collapse = ";"), character(1))
    res$collection <- cn
    res$contrast <- contrast
    res$parent_L2_joint <- L2
    out[[cn]] <- res
  }
  bind_rows(out)
}

all_results <- list()
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  files <- list.files(cdir, pattern = "(_vs_parity|_markers)\\.csv$", full.names = TRUE)
  cat(sprintf("\n[%s] %d L2 marker files\n", cname, length(files)))
  for (f in files) {
    L2_safe <- sub("(_vs_parity|_markers)\\.csv$", "", basename(f))
    L2 <- gsub("_", "-", gsub("__", "::", L2_safe))   # reverse the safe-name munge
    res <- run_one(f, cname, L2)
    if (!is.null(res) && nrow(res) > 0) {
      all_results[[paste(cname, L2_safe, sep = "::")]] <- res
      cat(sprintf("  [%s/%s] %d rows\n", cname, L2_safe, nrow(res)))
    }
  }
}

if (length(all_results) == 0) {
  cat("No L2 GSEA results. Exit.\n"); quit(status = 0)
}

big <- bind_rows(all_results) %>%
  select(contrast, parent_L2_joint, collection, pathway, NES, padj, pval,
          size, leadingEdge)
out_csv <- file.path(out_dir, "per_L2_nes.csv")
write_csv(big, out_csv)
cat(sprintf("\nWrote: %s (%d rows)\n", out_csv, nrow(big)))

top_sum <- big %>%
  group_by(contrast, parent_L2_joint, collection) %>%
  arrange(desc(NES)) %>% mutate(rank_up = row_number()) %>%
  arrange(NES) %>% mutate(rank_dn = row_number()) %>%
  filter(rank_up <= 5 | rank_dn <= 5) %>% ungroup() %>%
  arrange(contrast, parent_L2_joint, collection, desc(NES))
out_top <- file.path(out_dir, "per_L2_top_summary.csv")
write_csv(top_sum, out_top)
cat(sprintf("Wrote: %s (%d rows)\n", out_top, nrow(top_sum)))

cat("\n=== L2 GSEA done ===\n")
