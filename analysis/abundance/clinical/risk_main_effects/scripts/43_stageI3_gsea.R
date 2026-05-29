#!/usr/bin/env Rscript
# 43_stageI3_gsea.R
# Stage I.3 — fgsea preranked over per-(contrast × nhoodgroup) marker tables.
#
# Inputs:
#   stageF3_markers/<contrast>/<group>_vs_parent.csv  per-group F.3 markers
#     columns: gene_id, symbol, logFC, FDR, P.Value, ...
#
# For each (contrast, nhoodgroup) marker table:
#   - Build pre-ranked vector: gene_id -> sign(logFC) * -log10(P.Value)
#     (logFC alone is fine for ranking; we use the signed -log10p form so that
#      tied/near-zero logFC genes don't dominate).
#   - Run fgseaMultilevel against:
#       Hallmark (H)
#       Reactome (C2 CP:REACTOME)
#   - Aggregate top results into a single long table.
#
# Outputs:
#   stageI3_gsea/per_nhoodgroup_nes.csv
#     Columns: contrast, NhoodGroup_renamed, parent_L2_joint, parent_compartment,
#              pathway, collection, NES, padj, size, leadingEdge
#   stageI3_gsea/per_nhoodgroup_top_summary.csv  (top 5 enriched + 5 depleted per group)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(fgsea); library(msigdbr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "43_stageI3_gsea.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

out_dir <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea")
ensure_dir(out_dir)

NPROC <- as.integer(Sys.getenv("OMP_NUM_THREADS",
                    Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
cat(sprintf("nproc: %d\n", NPROC))

# Build pathway lists once -------------------------------------------------
cat("Loading MSigDB Hallmark + Reactome...\n")
load_collection <- function(category, subcategory = NULL) {
  if (is.null(subcategory)) {
    df <- msigdbr::msigdbr(species = "Homo sapiens", category = category)
  } else {
    df <- msigdbr::msigdbr(species = "Homo sapiens", category = category,
                            subcategory = subcategory)
  }
  split(df$ensembl_gene, df$gs_name)
}
hallmark <- load_collection("H")
reactome <- load_collection("C2", "CP:REACTOME")
collections <- list(hallmark = hallmark, reactome = reactome)
cat(sprintf("  hallmark: %d sets, reactome: %d sets\n",
            length(hallmark), length(reactome)))

# Iterate F.3 outputs -------------------------------------------------------
f3_root <- paths$outputs$stageF3
contrast_dirs <- list.dirs(f3_root, recursive = FALSE)
cat(sprintf("F.3 contrasts: %d\n", length(contrast_dirs)))

# Helpful index from F.1: NhoodGroup_renamed -> parent_L2_joint, parent_compartment.
build_group_index <- function() {
  rows <- list()
  f1 <- paths$outputs$stageF1
  for (cdir in list.dirs(f1, recursive = FALSE)) {
    cname <- basename(cdir)
    gs_csv <- file.path(cdir, "nhood_groups_summary.csv")
    if (!file.exists(gs_csv)) next
    gs <- read_csv(gs_csv, show_col_types = FALSE) %>%
      transmute(contrast = cname,
                NhoodGroup_renamed = as.character(NhoodGroup_renamed),
                parent_L2_joint = as.character(parent_L2_joint),
                parent_compartment = as.character(parent_compartment),
                n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group)))
    rows[[cname]] <- gs
  }
  bind_rows(rows)
}
group_index <- build_group_index()
cat(sprintf("Group index: %d rows\n", nrow(group_index)))

run_one_marker <- function(csv_path, contrast_name, group_renamed) {
  d <- tryCatch(read_csv(csv_path, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if (!all(c("gene_id", "logFC") %in% colnames(d))) {
    cat(sprintf("    [%s/%s] missing gene_id or logFC, skip\n",
                contrast_name, group_renamed))
    return(NULL)
  }
  pval_col <- intersect(c("P.Value", "PValue", "p_value", "pvalue"),
                         colnames(d))[1]
  if (is.na(pval_col)) {
    rank_vec <- d$logFC
  } else {
    rank_vec <- sign(d$logFC) * (-log10(pmax(d[[pval_col]], 1e-300)))
  }
  names(rank_vec) <- d$gene_id
  rank_vec <- rank_vec[!is.na(rank_vec) & is.finite(rank_vec)]
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  if (length(rank_vec) < 100) {
    cat(sprintf("    [%s/%s] only %d genes ranked, skip\n",
                contrast_name, group_renamed, length(rank_vec)))
    return(NULL)
  }
  out_rows <- list()
  for (coll_name in names(collections)) {
    pw <- collections[[coll_name]]
    res <- tryCatch(
      fgsea::fgseaMultilevel(pathways = pw, stats = rank_vec,
                              minSize = 10, maxSize = 500,
                              nPermSimple = 1000, eps = 0,
                              nproc = NPROC),
      error = function(e) {
        cat(sprintf("    fgsea error %s: %s\n", coll_name, conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res) || nrow(res) == 0) next
    res <- as.data.frame(res)
    res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse=";"),
                                character(1))
    res$collection <- coll_name
    res$contrast <- contrast_name
    res$NhoodGroup_renamed <- group_renamed
    out_rows[[coll_name]] <- res
  }
  bind_rows(out_rows)
}

all_results <- list()
n_total <- 0; n_done <- 0
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  files <- list.files(cdir, pattern = "_vs_parent\\.csv$", full.names = TRUE)
  cat(sprintf("\n[%s] %d marker files\n", cname, length(files)))
  n_total <- n_total + length(files)
  for (f in files) {
    group_renamed <- sub("_vs_parent\\.csv$", "", basename(f))
    res <- run_one_marker(f, cname, group_renamed)
    if (!is.null(res) && nrow(res) > 0) {
      all_results[[paste(cname, group_renamed, sep = "::")]] <- res
      n_done <- n_done + 1
    }
    if (n_done %% 10 == 0) cat(sprintf("  done: %d/%d\n", n_done, n_total))
  }
}
cat(sprintf("\nTotal marker tables processed: %d / %d\n", n_done, n_total))

if (length(all_results) == 0) {
  cat("No GSEA results produced — check F.3 outputs.\n"); quit(status = 1)
}

big <- bind_rows(all_results) %>%
  left_join(group_index %>%
              select(contrast, NhoodGroup_renamed, parent_L2_joint,
                     parent_compartment, n_nhoods_in_group),
            by = c("contrast", "NhoodGroup_renamed")) %>%
  select(contrast, NhoodGroup_renamed, parent_L2_joint, parent_compartment,
         n_nhoods_in_group, collection, pathway, NES, padj = padj, pval,
         size, leadingEdge)

out_full <- file.path(out_dir, "per_nhoodgroup_nes.csv")
write_csv(big, out_full)
cat(sprintf("\nWrote: %s (%d rows)\n", out_full, nrow(big)))

# Top summary: top 5 enriched + 5 depleted per (contrast × group × collection)
top_sum <- big %>%
  group_by(contrast, NhoodGroup_renamed, collection) %>%
  arrange(desc(NES)) %>%
  mutate(rank_up = row_number()) %>%
  arrange(NES) %>%
  mutate(rank_dn = row_number()) %>%
  filter(rank_up <= 5 | rank_dn <= 5) %>%
  ungroup() %>%
  arrange(contrast, NhoodGroup_renamed, collection, desc(NES))
out_top <- file.path(out_dir, "per_nhoodgroup_top_summary.csv")
write_csv(top_sum, out_top)
cat(sprintf("Wrote: %s (%d rows)\n", out_top, nrow(top_sum)))

cat("=== Stage I.3 GSEA done ===\n")
quit(status = 0)
