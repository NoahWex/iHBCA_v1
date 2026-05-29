#!/usr/bin/env Rscript
# 03_summarize.R — Cross-formula comparison report.
#
# Inputs:
#   outputs/design_audit.csv
#   outputs/de_results/{cohort_id}/{formula_id}/{l2}.csv
#
# Outputs:
#   outputs/comparisons/testability_matrix.csv
#   outputs/comparisons/hit_counts.csv
#   outputs/comparisons/jaccard_top50.csv      (per L2: formula_a, formula_b, jaccard)
#   outputs/comparisons/logfc_concordance.csv  (per L2: pearson_r on union of sig hits)
#   outputs/comparisons/recommended_formula_per_l2.csv (strictest design with signal)
#   reports/CP3_formula_choice.md              (narrative summary)
#
# Rigor ranking (highest to lowest):
#   B3_full > B2_study_facs > B1_study > B4_austin >        # Cohort B (balanced)
#   A3_full > A2_sampletype > A1_bio > A0_bare              # Cohort A (full)
# B-formulas outrank A-formulas because Cohort B can block study (the largest
# technical confound), so the same formula on B is strictly more rigorous.

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(dplyr)
  library(readr)
  library(tidyr)
})

p <- ArgumentParser()
p$add_argument("--project-root", required = TRUE)
p$add_argument("--inquiry-yaml", required = TRUE)
p$add_argument("--min-sig", type = "integer", default = 20,
               help = "Min significant genes (FDR<0.05) to count as 'signal'")
args <- p$parse_args()

project_root <- args$project_root
inquiry      <- yaml::read_yaml(args$inquiry_yaml)
fdr_thr      <- inquiry$significance$fdr_threshold
lfc_thr      <- inquiry$significance$logfc_threshold

out_dir   <- file.path(project_root, "outputs", "comparisons")
report_dir <- file.path(project_root, "reports")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

# Rigor ranking (descending — highest rigor first). 7 formulas after
# B4_austin was dropped (sample_type_coarse is 77% NA on iHBCA donors).
rigor_order <- c("B3_full", "B2_study_facs", "B1_study",
                 "A3_full", "A2_sampletype", "A1_bio", "A0_bare")

# ----- load audit -------------------------------------------------------
audit <- readr::read_csv(file.path(project_root, "outputs", "design_audit.csv"),
                          show_col_types = FALSE)
cat(sprintf("Loaded audit: %d (cohort x formula x L2) rows\n", nrow(audit)))

# ----- testability matrix ----------------------------------------------
testability <- audit %>%
  dplyr::select(L2, formula_id, testable) %>%
  tidyr::pivot_wider(names_from = formula_id, values_from = testable,
                     values_fill = FALSE)
# Reorder columns by rigor
cols_present <- intersect(rigor_order, colnames(testability))
testability <- testability %>%
  dplyr::select(L2, dplyr::all_of(cols_present))
readr::write_csv(testability, file.path(out_dir, "testability_matrix.csv"))
cat(sprintf("Wrote testability_matrix.csv: %d L2 x %d formulas\n",
            nrow(testability), length(cols_present)))

# ----- hit counts -------------------------------------------------------
hit_counts <- audit %>%
  dplyr::filter(testable) %>%
  dplyr::select(L2, formula_id, n_sig_fdr05, n_sig_fdr05_lfc05) %>%
  tidyr::pivot_wider(names_from = formula_id,
                     values_from = c(n_sig_fdr05, n_sig_fdr05_lfc05),
                     values_fill = NA_integer_)
readr::write_csv(hit_counts, file.path(out_dir, "hit_counts.csv"))
cat(sprintf("Wrote hit_counts.csv\n"))

# ----- per-L2 cross-formula hit overlap + logFC concordance -----------
jaccard_top50 <- list()
logfc_concord <- list()

load_de <- function(cohort_id, formula_id, l2_safe) {
  path <- file.path(project_root, "outputs", "de_results",
                    cohort_id, formula_id, paste0(l2_safe, ".csv"))
  if (!file.exists(path)) return(NULL)
  # Skipped-stub CSVs (per stageF3_L2_markers_risk.R:189 convention) have
  # only `skipped, reason` columns. Skip them in summary loading.
  hdr <- readr::read_lines(path, n_max = 1)
  if (grepl("^skipped", hdr)) return(NULL)
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

safe_name <- function(s) gsub("[^A-Za-z0-9_]", "_", s)

testable_pairs_by_l2 <- audit %>%
  dplyr::filter(testable) %>%
  dplyr::select(cohort_id, formula_id, L2) %>%
  dplyr::group_split(L2)

for (l2_grp in testable_pairs_by_l2) {
  if (nrow(l2_grp) < 2) next
  l2 <- l2_grp$L2[1]
  l2_safe <- safe_name(l2)
  # Load all testable DE tables for this L2
  de_tables <- list()
  for (i in seq_len(nrow(l2_grp))) {
    fid <- l2_grp$formula_id[i]
    cid <- l2_grp$cohort_id[i]
    tt <- load_de(cid, fid, l2_safe)
    if (!is.null(tt)) de_tables[[fid]] <- tt
  }
  if (length(de_tables) < 2) next

  formula_ids <- names(de_tables)
  # Pairwise comparisons
  for (i in 1:(length(formula_ids) - 1)) {
    for (j in (i + 1):length(formula_ids)) {
      a <- formula_ids[i]; b <- formula_ids[j]
      tt_a <- de_tables[[a]]; tt_b <- de_tables[[b]]

      # Top-50 Jaccard
      top_a <- tt_a %>% dplyr::arrange(adj.P.Val) %>%
                       dplyr::slice_head(n = 50) %>% dplyr::pull(gene_id)
      top_b <- tt_b %>% dplyr::arrange(adj.P.Val) %>%
                       dplyr::slice_head(n = 50) %>% dplyr::pull(gene_id)
      jacc <- length(intersect(top_a, top_b)) /
              length(union(top_a, top_b))
      jaccard_top50[[length(jaccard_top50) + 1L]] <- data.frame(
        L2 = l2, formula_a = a, formula_b = b, jaccard_top50 = jacc
      )

      # logFC concordance on union of significant genes
      sig_a <- tt_a %>% dplyr::filter(adj.P.Val < fdr_thr) %>% dplyr::pull(gene_id)
      sig_b <- tt_b %>% dplyr::filter(adj.P.Val < fdr_thr) %>% dplyr::pull(gene_id)
      union_sig <- union(sig_a, sig_b)
      if (length(union_sig) >= 10) {
        a_sub <- tt_a %>% dplyr::filter(gene_id %in% union_sig) %>%
                          dplyr::select(gene_id, logFC_a = logFC)
        b_sub <- tt_b %>% dplyr::filter(gene_id %in% union_sig) %>%
                          dplyr::select(gene_id, logFC_b = logFC)
        merged <- dplyr::inner_join(a_sub, b_sub, by = "gene_id")
        pearson_r <- suppressWarnings(stats::cor(merged$logFC_a, merged$logFC_b))
        median_abs_delta <- median(abs(merged$logFC_a - merged$logFC_b))
      } else {
        pearson_r <- NA_real_
        median_abs_delta <- NA_real_
      }
      logfc_concord[[length(logfc_concord) + 1L]] <- data.frame(
        L2 = l2, formula_a = a, formula_b = b,
        n_union_sig = length(union_sig),
        pearson_r = pearson_r,
        median_abs_delta_lfc = median_abs_delta
      )
    }
  }
}

jaccard_df <- dplyr::bind_rows(jaccard_top50)
logfc_df   <- dplyr::bind_rows(logfc_concord)
readr::write_csv(jaccard_df, file.path(out_dir, "jaccard_top50.csv"))
readr::write_csv(logfc_df,   file.path(out_dir, "logfc_concordance.csv"))
cat(sprintf("Wrote jaccard_top50.csv (%d rows) and logfc_concordance.csv (%d rows)\n",
            nrow(jaccard_df), nrow(logfc_df)))

# ----- recommended formula per L2 (strictest with signal) -------------
# For each L2: walk the rigor_order list; the first formula that is
# testable AND has >= min_sig significant genes is the recommendation.
recs <- audit %>%
  dplyr::filter(testable) %>%
  dplyr::mutate(rigor_rank = match(formula_id, rigor_order)) %>%
  dplyr::group_by(L2) %>%
  dplyr::arrange(rigor_rank, .by_group = TRUE) %>%
  dplyr::mutate(has_signal = n_sig_fdr05 >= args$min_sig) %>%
  dplyr::filter(has_signal) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::select(L2, recommended_formula = formula_id,
                cohort_id, n_sig_fdr05, n_sig_fdr05_lfc05,
                n_AR, n_BR1)

# L2s with NO formula meeting min_sig
no_signal <- setdiff(unique(audit$L2), recs$L2)
no_signal_df <- data.frame(L2 = no_signal,
                            recommended_formula = NA_character_,
                            cohort_id = NA_character_,
                            n_sig_fdr05 = 0L,
                            n_sig_fdr05_lfc05 = 0L,
                            n_AR = NA_integer_, n_BR1 = NA_integer_)
recs_full <- dplyr::bind_rows(recs, no_signal_df) %>%
  dplyr::arrange(is.na(recommended_formula), L2)
readr::write_csv(recs_full, file.path(out_dir, "recommended_formula_per_l2.csv"))
cat(sprintf("Wrote recommended_formula_per_l2.csv (%d with signal, %d without)\n",
            nrow(recs), length(no_signal)))

# ----- narrative report -----------------------------------------------
report_path <- file.path(report_dir, "CP3_formula_choice.md")
sink(report_path)
cat("# CP3 — Cross-formula comparison and recommended designs\n\n")
cat(sprintf("Generated: %s\n", Sys.time()))
cat(sprintf("Min significant genes (FDR<%.2f, |logFC|>%.1f) to count as signal: %d\n\n",
            fdr_thr, lfc_thr, args$min_sig))

cat("## Testability summary (n L2 testable per formula)\n\n")
ts <- audit %>%
  dplyr::group_by(formula_id) %>%
  dplyr::summarise(n_testable = sum(testable), n_total = dplyr::n(),
                    .groups = "drop") %>%
  dplyr::mutate(rigor_rank = match(formula_id, rigor_order)) %>%
  dplyr::arrange(rigor_rank) %>%
  dplyr::select(-rigor_rank)
print(as.data.frame(ts))
cat("\n")

cat("## Hit-count summary (median across testable L2s)\n\n")
hs <- audit %>%
  dplyr::filter(testable) %>%
  dplyr::group_by(formula_id) %>%
  dplyr::summarise(median_n_sig = stats::median(n_sig_fdr05),
                    max_n_sig    = max(n_sig_fdr05),
                    n_L2_with_signal = sum(n_sig_fdr05 >= args$min_sig),
                    .groups = "drop") %>%
  dplyr::mutate(rigor_rank = match(formula_id, rigor_order)) %>%
  dplyr::arrange(rigor_rank) %>%
  dplyr::select(-rigor_rank)
print(as.data.frame(hs))
cat("\n")

cat("## Recommended formula per L2 (strictest design with signal)\n\n")
cat(sprintf("L2 types with signal: %d\n", nrow(recs)))
cat(sprintf("L2 types with NO formula passing min_sig=%d: %d\n\n",
            args$min_sig, length(no_signal)))
print(as.data.frame(recs_full))

cat("\n## Caveats\n\n")
for (cav in inquiry$caveats) {
  cat(sprintf("- **%s**: %s\n", cav$id, cav$description))
}
sink()
cat(sprintf("Wrote report: %s\n", report_path))

cat("\nDone.\n")
