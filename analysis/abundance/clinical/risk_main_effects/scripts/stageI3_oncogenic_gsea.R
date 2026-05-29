#!/usr/bin/env Rscript
# stageI3_oncogenic_gsea.R
# Extend GSEA to MSigDB C6 (oncogenic signatures) + C4 (cancer-oriented).
# Runs against both F.3 NhoodGroup markers and F.3 L2 markers.
#
# Outputs:
#   outputs/stageI3_gsea/per_nhoodgroup_C6C4_nes.csv
#   outputs/stageI3_gsea/per_L2_C6C4_nes.csv

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fgsea); library(msigdbr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI3_oncogenic_gsea.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

out_dir <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea")
ensure_dir(out_dir)

cat("Loading MSigDB C6 (oncogenic) + C4 (cancer-oriented) collections...\n")

# C6: oncogenic signatures (no subcategory)
c6_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C6")
c6 <- split(c6_df$ensembl_gene, c6_df$gs_name)
cat(sprintf("  C6 oncogenic: %d sets\n", length(c6)))

# C4: cancer-oriented (has multiple subcategories — pull all under C4)
c4_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "C4")
# Tag the subcategory in pathway name so we can disambiguate at output time
c4_df$gs_name_sub <- paste(c4_df$gs_subcat, c4_df$gs_name, sep = "::")
c4 <- split(c4_df$ensembl_gene, c4_df$gs_name_sub)
cat(sprintf("  C4 cancer: %d sets across subcats: %s\n",
            length(c4),
            paste(unique(c4_df$gs_subcat), collapse = ", ")))

collections <- list(C6 = c6, C4 = c4)

# Helper: run fgsea on one marker file
run_one <- function(csv, contrast, key_col_name, key_val) {
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
      error = function(e) {cat(sprintf("  fgsea err: %s\n", conditionMessage(e))); NULL}
    )
    if (is.null(res) || nrow(res) == 0) next
    res <- as.data.frame(res)
    res$leadingEdge <- vapply(res$leadingEdge,
                                function(x) paste(x, collapse = ";"), character(1))
    res$collection <- cn
    res$contrast <- contrast
    res[[key_col_name]] <- key_val
    out[[cn]] <- res
  }
  bind_rows(out)
}

# ---------- NhoodGroup grain ----------
cat("\n=== NhoodGroup grain ===\n")
ng_results <- list()
for (cdir in list.dirs(paths$outputs$stageF3, recursive = FALSE)) {
  cname <- basename(cdir)
  files <- list.files(cdir, pattern = "_vs_parent\\.csv$", full.names = TRUE)
  cat(sprintf("[%s] %d files\n", cname, length(files)))
  for (f in files) {
    g <- sub("_vs_parent\\.csv$", "", basename(f))
    res <- run_one(f, cname, "NhoodGroup_renamed", g)
    if (!is.null(res) && nrow(res) > 0) {
      ng_results[[paste(cname, g, sep = "::")]] <- res
    }
  }
}
if (length(ng_results) > 0) {
  big_ng <- bind_rows(ng_results) %>%
    select(contrast, NhoodGroup_renamed, collection, pathway,
            NES, padj, pval, size, leadingEdge)
  write_csv(big_ng, file.path(out_dir, "per_nhoodgroup_C6C4_nes.csv"))
  cat(sprintf("Wrote per_nhoodgroup_C6C4_nes.csv: %d rows\n", nrow(big_ng)))
}

# ---------- L2 grain ----------
cat("\n=== L2 grain ===\n")
l2_root <- file.path(paths$inquiry_root, "outputs", "stageF3_L2_markers")
l2_results <- list()
for (cdir in list.dirs(l2_root, recursive = FALSE)) {
  cname <- basename(cdir)
  files <- list.files(cdir, pattern = "_vs_parity\\.csv$", full.names = TRUE)
  cat(sprintf("[%s] %d files\n", cname, length(files)))
  for (f in files) {
    L2_safe <- sub("_vs_parity\\.csv$", "", basename(f))
    L2 <- gsub("_", "-", gsub("__", "::", L2_safe))
    res <- run_one(f, cname, "parent_L2_joint", L2)
    if (!is.null(res) && nrow(res) > 0) {
      l2_results[[paste(cname, L2_safe, sep = "::")]] <- res
    }
  }
}
if (length(l2_results) > 0) {
  big_l2 <- bind_rows(l2_results) %>%
    select(contrast, parent_L2_joint, collection, pathway,
            NES, padj, pval, size, leadingEdge)
  write_csv(big_l2, file.path(out_dir, "per_L2_C6C4_nes.csv"))
  cat(sprintf("Wrote per_L2_C6C4_nes.csv: %d rows\n", nrow(big_l2)))
}

cat("\n=== oncogenic GSEA done ===\n")
