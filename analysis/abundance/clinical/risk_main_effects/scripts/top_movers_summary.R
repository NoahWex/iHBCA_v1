#!/usr/bin/env Rscript
# top_movers_summary.R
# Per-(contrast × parent_L2_joint) "top mover" summary:
#   - Pick the viable NhoodGroup (n_nhoods >= 10) with max |group_med_lfc|
#   - If no viable NhoodGroup, fall back to L2-grain
#   - Pull top 8 UP + 8 DOWN markers (limma-voom, by FDR)
#   - Pull top 5 UP + 5 DOWN GSEA pathways (Hallmark + Reactome, padj<0.05)
#
# Output:
#   outputs/top_movers_summary.csv    long table
#   outputs/top_movers_summary.md     markdown report (one section per row)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "top_movers_summary.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
out_dir <- file.path(paths$inquiry_root, "outputs")

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1",
                "parity_x_HR_sporadic", "parity_x_HR_BRCA2")
N_MARKERS <- 8
N_PATHWAYS <- 5
VIABLE_MIN <- 10

# Convert canonical NhoodGroup_renamed to filesystem-safe form (matches F.3 + GSEA)
fs_safe <- function(x) gsub("-", "_", gsub("::", "__", x))

# Pre-load per-nhoodgroup GSEA (one big file; index by contrast+ng_safe)
gsea_ng <- read_csv(file.path(out_dir, "stageI3_gsea", "per_nhoodgroup_nes.csv"),
                     show_col_types = FALSE)
gsea_ng$collection <- toupper(gsea_ng$collection)
gsea_l2 <- read_csv(file.path(out_dir, "stageI3_gsea", "per_L2_nes.csv"),
                     show_col_types = FALSE)
gsea_l2$collection <- toupper(gsea_l2$collection)

# Helper: extract top markers (by FDR) for a marker file
read_markers <- function(csv) {
  d <- tryCatch(read_csv(csv, show_col_types = FALSE,
                          col_types = cols(.default = "?")),
                 error = function(e) NULL)
  if (is.null(d) || nrow(d) == 0) return(NULL)
  if ("skipped" %in% colnames(d)) return(NULL)
  if (!all(c("logFC", "gene_id") %in% colnames(d))) return(NULL)
  fdr_col <- intersect(c("FDR", "adj.P.Val"), colnames(d))[1]
  if (is.na(fdr_col)) return(NULL)
  d %>%
    transmute(symbol = if ("symbol" %in% colnames(.)) symbol else NA_character_,
               gene_id, logFC, FDR = .data[[fdr_col]])
}

# Helper: does a NhoodGroup's F.3 marker file have real markers (not stub,
# at least one FDR<0.05)?
has_real_markers <- function(csv) {
  if (!file.exists(csv)) return(FALSE)
  if (file.size(csv) < 200) return(FALSE)   # likely stub
  d <- read_markers(csv)
  if (is.null(d) || nrow(d) == 0) return(FALSE)
  any(!is.na(d$FDR) & d$FDR < 0.05)
}

format_markers <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(list(up = "", dn = ""))
  up <- d %>% filter(logFC > 0, FDR < 0.05) %>%
    arrange(FDR) %>% head(N_MARKERS) %>%
    mutate(label = sprintf("%s (lfc=%+.2f, FDR=%.2g)",
                            ifelse(is.na(symbol), gene_id, symbol),
                            logFC, FDR))
  dn <- d %>% filter(logFC < 0, FDR < 0.05) %>%
    arrange(FDR) %>% head(N_MARKERS) %>%
    mutate(label = sprintf("%s (lfc=%+.2f, FDR=%.2g)",
                            ifelse(is.na(symbol), gene_id, symbol),
                            logFC, FDR))
  list(up = paste(up$label, collapse = "; "),
        dn = paste(dn$label, collapse = "; "))
}

format_gsea <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(list(up = "", dn = ""))
  up <- d %>% filter(NES > 0, padj < 0.05) %>%
    arrange(desc(NES)) %>% head(N_PATHWAYS) %>%
    mutate(label = sprintf("%s [%s] (NES=%+.2f, padj=%.2g)",
                            pathway, collection, NES, padj))
  dn <- d %>% filter(NES < 0, padj < 0.05) %>%
    arrange(NES) %>% head(N_PATHWAYS) %>%
    mutate(label = sprintf("%s [%s] (NES=%+.2f, padj=%.2g)",
                            pathway, collection, NES, padj))
  list(up = paste(up$label, collapse = "; "),
        dn = paste(dn$label, collapse = "; "))
}

# Build long table
long_rows <- list()

for (cname in CONTRASTS) {
  cat(sprintf("\n=== %s ===\n", cname))
  f1_csv <- file.path(paths$outputs$stageF1, cname, "nhood_groups_summary.csv")
  if (!file.exists(f1_csv)) { cat("  no F.1 summary, skip\n"); next }
  ng <- read_csv(f1_csv, show_col_types = FALSE) %>%
    mutate(n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group)),
            group_med_lfc      = suppressWarnings(as.numeric(group_med_lfc)),
            n_sig              = suppressWarnings(as.numeric(n_sig)))

  parent_L2s <- unique(na.omit(ng$parent_L2_joint))
  cat(sprintf("  parent L2s in F.1: %d\n", length(parent_L2s)))

  for (L2 in parent_L2s) {
    ng_sub <- ng %>% filter(parent_L2_joint == L2,
                              !is.na(n_nhoods_in_group),
                              n_nhoods_in_group >= VIABLE_MIN) %>%
      arrange(desc(abs(group_med_lfc)))
    # Filter to NhoodGroups whose F.3 file has real markers (skip Mode 1 stubs)
    if (nrow(ng_sub) > 0) {
      ng_sub$has_markers <- vapply(ng_sub$NhoodGroup_renamed, function(name) {
        has_real_markers(file.path(paths$outputs$stageF3, cname,
                                     paste0(fs_safe(name), "_vs_parent.csv")))
      }, logical(1))
      ng_sub <- ng_sub %>% filter(has_markers)
    }
    use_l2 <- nrow(ng_sub) == 0
    if (use_l2) {
      # Fallback: L2-grain
      L2_safe <- fs_safe(L2)
      mfile <- file.path(paths$inquiry_root, "outputs", "stageF3_L2_markers",
                          cname, paste0(L2_safe, "_vs_parity.csv"))
      if (!file.exists(mfile)) {
        cat(sprintf("    [%s] no viable NG, no L2 markers — skip\n", L2)); next
      }
      m <- read_markers(mfile)
      g <- gsea_l2 %>% filter(contrast == cname,
                                parent_L2_joint == L2 |
                                  parent_L2_joint == sub("-", "_", L2) |
                                  parent_L2_joint == gsub("-", "_", L2))
      mfmt <- format_markers(m); gfmt <- format_gsea(g)
      long_rows[[length(long_rows)+1]] <- tibble(
        contrast = cname, parent_L2 = L2, level = "L2",
        top_mover_id = NA_character_,
        n_nhoods_top = NA_integer_, n_sig_top = NA_integer_,
        lfc_top = NA_real_,
        markers_up = mfmt$up, markers_dn = mfmt$dn,
        gsea_up = gfmt$up, gsea_dn = gfmt$dn)
      cat(sprintf("    [%s] L2-grain (no viable NG)\n", L2))
    } else {
      top <- ng_sub %>% arrange(desc(abs(group_med_lfc))) %>% slice(1)
      ng_safe <- fs_safe(top$NhoodGroup_renamed)
      mfile <- file.path(paths$outputs$stageF3, cname,
                          paste0(ng_safe, "_vs_parent.csv"))
      m <- read_markers(mfile)
      g <- gsea_ng %>% filter(contrast == cname, NhoodGroup_renamed == ng_safe)
      mfmt <- format_markers(m); gfmt <- format_gsea(g)
      long_rows[[length(long_rows)+1]] <- tibble(
        contrast = cname, parent_L2 = L2, level = "NhoodGroup",
        top_mover_id = top$NhoodGroup_renamed,
        n_nhoods_top = as.integer(top$n_nhoods_in_group),
        n_sig_top = as.integer(top$n_sig),
        lfc_top = top$group_med_lfc,
        markers_up = mfmt$up, markers_dn = mfmt$dn,
        gsea_up = gfmt$up, gsea_dn = gfmt$dn)
      cat(sprintf("    [%s] top=%s lfc=%+.2f n=%d\n",
                    L2, top$NhoodGroup_renamed,
                    top$group_med_lfc, as.integer(top$n_nhoods_in_group)))
    }
  }
}

long <- bind_rows(long_rows) %>%
  arrange(contrast, desc(abs(lfc_top %||% 0)), parent_L2)
cat(sprintf("\nTotal rows: %d\n", nrow(long)))

out_csv <- file.path(out_dir, "top_movers_summary.csv")
write_csv(long, out_csv)
cat(sprintf("Wrote: %s\n", out_csv))

# Markdown render
md <- c("# Top movers per (contrast × L2)",
        "",
        sprintf("Generated %s. For each (contrast × parent_L2_joint): the NhoodGroup with max |group_med_lfc| (n_nhoods ≥ %d). Falls back to L2-grain when no viable NhoodGroup exists in that contrast.",
                Sys.Date(), VIABLE_MIN),
        "",
        "Markers: top 8 UP + 8 DOWN by FDR (limma-voom). Pathways: top 5 UP + 5 DOWN by NES (padj<0.05).",
        "",
        "---",
        "")
for (cname in CONTRASTS) {
  rows <- long %>% filter(contrast == cname)
  if (nrow(rows) == 0) next
  md <- c(md, sprintf("## %s", cname), "")
  for (i in seq_len(nrow(rows))) {
    r <- rows[i, ]
    head_line <- if (r$level == "L2")
      sprintf("### %s  →  L2-grain (no viable sub-L2 NhoodGroup in this cohort)",
                r$parent_L2)
    else
      sprintf("### %s  →  top mover: %s  (lfc=%+.2f, n_nhoods=%d, n_sig=%d)",
                r$parent_L2, r$top_mover_id, r$lfc_top,
                r$n_nhoods_top, r$n_sig_top)
    md <- c(md, head_line, "")
    md <- c(md, "**Markers UP:**", "")
    md <- c(md, ifelse(r$markers_up == "", "_(none at FDR<0.05)_",
                        gsub("; ", "  \n", r$markers_up)))
    md <- c(md, "", "**Markers DOWN:**", "")
    md <- c(md, ifelse(r$markers_dn == "", "_(none at FDR<0.05)_",
                        gsub("; ", "  \n", r$markers_dn)))
    md <- c(md, "", "**GSEA UP:**", "")
    md <- c(md, ifelse(r$gsea_up == "", "_(none at padj<0.05)_",
                        gsub("; ", "  \n", r$gsea_up)))
    md <- c(md, "", "**GSEA DOWN:**", "")
    md <- c(md, ifelse(r$gsea_dn == "", "_(none at padj<0.05)_",
                        gsub("; ", "  \n", r$gsea_dn)))
    md <- c(md, "", "---", "")
  }
}
out_md <- file.path(out_dir, "top_movers_summary.md")
writeLines(md, out_md)
cat(sprintf("Wrote: %s (%d lines)\n", out_md, length(md)))
cat("\n=== top movers summary done ===\n")
