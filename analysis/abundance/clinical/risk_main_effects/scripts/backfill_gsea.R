#!/usr/bin/env Rscript
# backfill_gsea.R
# Identify F.3 marker files that produced no rows in per_nhoodgroup_nes.csv,
# rerun fgsea (nproc=1 to avoid parallel race), append to a backfill table.
# Then re-annotate the union to produce per_nhoodgroup_nes_annotated.csv v2.
#
# Mode 1 (F.3 skip stubs) are tagged but not retried — they need a different
# reference frame (group-vs-rest-of-compartment), out of scope for this fix.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(fgsea); library(msigdbr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "backfill_gsea.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

out_dir <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea")
existing_path <- file.path(out_dir, "per_nhoodgroup_nes.csv")

# ------------------------------------------------------------------------
# 1. Load existing GSEA — get (contrast, NhoodGroup_renamed) covered set
# ------------------------------------------------------------------------
existing <- read_csv(existing_path, show_col_types = FALSE)
covered <- existing %>% distinct(contrast, NhoodGroup_renamed)
cat(sprintf("Existing GSEA: %d rows, %d (contrast×group) tuples\n",
            nrow(existing), nrow(covered)))

# ------------------------------------------------------------------------
# 2. Enumerate all F.3 marker files; filter to "real markers" (not stubs)
# ------------------------------------------------------------------------
classify_file <- function(csv) {
  d <- tryCatch(read_csv(csv, show_col_types = FALSE, n_max = 5,
                          progress = FALSE),
                 error = function(e) NULL)
  if (is.null(d)) return("read_error")
  if ("skipped" %in% colnames(d)) return("stub")
  if (!all(c("gene_id", "logFC") %in% colnames(d))) return("missing_cols")
  return("real")
}

f3_root <- paths$outputs$stageF3
all_rows <- list()
for (cdir in list.dirs(f3_root, recursive = FALSE)) {
  cname <- basename(cdir)
  files <- list.files(cdir, pattern = "_vs_parent\\.csv$", full.names = TRUE)
  for (f in files) {
    g <- sub("_vs_parent\\.csv$", "", basename(f))
    cls <- classify_file(f)
    all_rows[[length(all_rows) + 1]] <-
      tibble(contrast = cname, NhoodGroup_renamed = g, file = f, status = cls)
  }
}
files_df <- bind_rows(all_rows)
cat(sprintf("\nF.3 marker files enumerated: %d\n", nrow(files_df)))
cat("  by status:\n")
print(files_df %>% count(status))

# ------------------------------------------------------------------------
# 3. Determine which need backfill: real markers AND not in covered set
# ------------------------------------------------------------------------
to_backfill <- files_df %>% filter(status == "real") %>%
  anti_join(covered, by = c("contrast", "NhoodGroup_renamed"))
cat(sprintf("\nTo backfill: %d real-marker groups missing from GSEA\n",
            nrow(to_backfill)))
print(to_backfill %>% count(contrast))

if (nrow(to_backfill) == 0) {
  cat("Nothing to backfill — exiting cleanly.\n")
  quit(status = 0)
}

# ------------------------------------------------------------------------
# 4. Load pathway collections
# ------------------------------------------------------------------------
cat("\nLoading MSigDB Hallmark + Reactome...\n")
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

# ------------------------------------------------------------------------
# 5. Run fgsea per missing group (nproc=1 to avoid parallel race)
# ------------------------------------------------------------------------
run_one <- function(csv_path, contrast_name, group_renamed) {
  d <- tryCatch(read_csv(csv_path, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if (!all(c("gene_id", "logFC") %in% colnames(d))) return(NULL)
  pval_col <- intersect(c("P.Value", "PValue", "p_value", "pvalue"),
                          colnames(d))[1]
  rank_vec <- if (is.na(pval_col)) d$logFC
                else sign(d$logFC) * (-log10(pmax(d[[pval_col]], 1e-300)))
  names(rank_vec) <- d$gene_id
  rank_vec <- rank_vec[!is.na(rank_vec) & is.finite(rank_vec)]
  rank_vec <- sort(rank_vec, decreasing = TRUE)
  if (length(rank_vec) < 100) {
    cat(sprintf("    [%s/%s] only %d ranked genes, skip\n",
                contrast_name, group_renamed, length(rank_vec)))
    return(NULL)
  }
  out_rows <- list()
  for (coll_name in names(collections)) {
    res <- tryCatch(
      fgsea::fgseaMultilevel(pathways = collections[[coll_name]],
                              stats = rank_vec,
                              minSize = 10, maxSize = 500,
                              nPermSimple = 1000, eps = 0,
                              nproc = 1),
      error = function(e) {
        cat(sprintf("    [%s/%s/%s] fgsea error: %s\n",
                    contrast_name, group_renamed, coll_name,
                    conditionMessage(e)))
        NULL
      }
    )
    if (is.null(res) || nrow(res) == 0) next
    res <- as.data.frame(res)
    res$leadingEdge <- vapply(res$leadingEdge,
                                function(x) paste(x, collapse = ";"),
                                character(1))
    res$collection <- coll_name
    res$contrast <- contrast_name
    res$NhoodGroup_renamed <- group_renamed
    out_rows[[coll_name]] <- res
  }
  bind_rows(out_rows)
}

backfill_results <- list()
for (i in seq_len(nrow(to_backfill))) {
  r <- to_backfill[i, ]
  cat(sprintf("[%d/%d] %s / %s\n", i, nrow(to_backfill),
              r$contrast, r$NhoodGroup_renamed))
  res <- run_one(r$file, r$contrast, r$NhoodGroup_renamed)
  if (!is.null(res) && nrow(res) > 0) {
    backfill_results[[paste(r$contrast, r$NhoodGroup_renamed, sep = "::")]] <- res
    cat(sprintf("    -> %d pathways\n", nrow(res)))
  }
}

if (length(backfill_results) == 0) {
  cat("No backfill rows produced — exiting.\n")
  quit(status = 0)
}

# ------------------------------------------------------------------------
# 6. Merge with existing, write updated GSEA table + re-annotate
# ------------------------------------------------------------------------
new_rows <- bind_rows(backfill_results)
cat(sprintf("\nBackfilled: %d rows across %d (contrast×group) tuples\n",
            nrow(new_rows), length(backfill_results)))

# Match existing schema
new_rows <- new_rows %>%
  transmute(contrast, NhoodGroup_renamed,
             parent_L2_joint = NA_character_,
             parent_compartment = NA_character_,
             n_nhoods_in_group = NA_real_,
             collection, pathway, NES, padj, pval,
             size, leadingEdge)

merged <- bind_rows(existing, new_rows)
write_csv(merged, existing_path)
cat(sprintf("Wrote: %s (%d rows total, was %d, +%d)\n",
            existing_path, nrow(merged), nrow(existing), nrow(new_rows)))

# Re-annotate parent_L2 via NhoodGroup_renamed normalization (canonical -> fs)
parent_lookup_rows <- list()
for (cdir in list.dirs(paths$outputs$stageF1, recursive = FALSE)) {
  cname <- basename(cdir)
  gs_csv <- file.path(cdir, "nhood_groups_summary.csv")
  if (!file.exists(gs_csv)) next
  gs <- read_csv(gs_csv, show_col_types = FALSE) %>%
    transmute(NhoodGroup_renamed = as.character(NhoodGroup_renamed),
              parent_L2_joint = as.character(parent_L2_joint),
              parent_compartment = as.character(parent_compartment),
              n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group))) %>%
    mutate(NhoodGroup_renamed = gsub("-", "_",
                                       gsub("::", "__", NhoodGroup_renamed)),
           contrast = cname)
  parent_lookup_rows[[cname]] <- gs
}
parent_lookup <- bind_rows(parent_lookup_rows)

annotated <- merged %>%
  select(-any_of(c("parent_L2_joint", "parent_compartment", "n_nhoods_in_group"))) %>%
  left_join(parent_lookup, by = c("contrast", "NhoodGroup_renamed"))
cat(sprintf("Annotated rows: %d / %d (%.1f%% with parent_L2)\n",
            sum(!is.na(annotated$parent_L2_joint)), nrow(annotated),
            100 * sum(!is.na(annotated$parent_L2_joint)) / nrow(annotated)))

annot_path <- file.path(out_dir, "per_nhoodgroup_nes_annotated.csv")
write_csv(annotated, annot_path)
cat(sprintf("Wrote: %s\n", annot_path))

cat("\n=== gsea backfill done ===\n")
