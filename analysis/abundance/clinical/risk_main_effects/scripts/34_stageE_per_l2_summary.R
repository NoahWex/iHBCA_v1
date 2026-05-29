#!/usr/bin/env Rscript
# 34_stageE_per_l2_summary.R
# Stage E — per-L2 effect summary across ALL Stage D contrasts.
# Builds one row per (L2 = compartment::label) × contrast, plus a wide
# pivoted summary that pairs main and interaction contrasts to assign
# Pattern α / Pattern β / null quadrants.
#
# Inputs (auto-discovered): all stageD_da_results/<contrast>/da_results.csv
#
# Outputs:
#   stageE_per_l2_summary.csv:
#     columns: L2_joint, L2_compartment, L2_label,
#              <contrast>_n, <contrast>_n_sig, <contrast>_med_lfc,
#              <contrast>_pct_up, <contrast>_<study>_n_sig, ...
#     plus quadrant assignment for paired (main, interaction) contrasts
#     declared via inquiry$contrast_pairs (optional).
#
# A long-form table is also written for easy downstream joins:
#   stageE_per_l2_long.csv

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(yaml); library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "34_stageE_per_l2_summary.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry
in_dir <- paths$outputs$stageD
out_csv <- paths$outputs$stageE
ensure_dir(dirname(out_csv))

cat(sprintf("=== Stage E: per-L2 summary for '%s' ===\n", inquiry$inquiry))

contrast_dirs <- list.dirs(in_dir, recursive = FALSE)
if (length(contrast_dirs) == 0) stop("No Stage D contrast outputs in: ", in_dir)

# -- helpers ------------------------------------------------------------------

summarise_one <- function(da, contrast_name, fdr_thresh = 0.05) {
  # per L2_joint: n_sig, n, median_lfc, pct_up
  if (!"compartment" %in% colnames(da)) {
    if ("L2_compartment" %in% colnames(da))
      da$compartment <- da$L2_compartment
    else stop("DA results missing compartment column for ", contrast_name)
  }
  if (!"label" %in% colnames(da)) {
    if ("L2_label" %in% colnames(da)) da$label <- da$L2_label
    else stop("DA results missing label column for ", contrast_name)
  }
  da$L2_joint <- paste(da$compartment, da$label, sep = "::")
  da$is_sig   <- !is.na(da$SpatialFDR) & da$SpatialFDR < fdr_thresh

  by_l2 <- da %>%
    group_by(L2_joint, compartment, label) %>%
    summarise(
      n         = n(),
      n_sig     = sum(is_sig),
      n_sig_up  = sum(is_sig & logFC > 0),
      n_sig_dn  = sum(is_sig & logFC < 0),
      med_lfc   = median(logFC, na.rm = TRUE),
      pct_up    = round(100 * sum(logFC > 0, na.rm = TRUE) / n(), 1),
      pct_dn    = round(100 * sum(logFC < 0, na.rm = TRUE) / n(), 1),
      .groups = "drop"
    ) %>%
    mutate(contrast = contrast_name)
  by_l2
}

per_study_one <- function(da, contrast_name, fdr_thresh = 0.05) {
  study_col <- intersect(c("index_study", "study"), colnames(da))[1]
  if (is.na(study_col)) return(NULL)
  if (!"L2_joint" %in% colnames(da)) {
    if ("L2_compartment" %in% colnames(da)) da$compartment <- da$L2_compartment
    if ("L2_label" %in% colnames(da)) da$label <- da$L2_label
    da$L2_joint <- paste(da$compartment, da$label, sep = "::")
  }
  da$is_sig <- !is.na(da$SpatialFDR) & da$SpatialFDR < fdr_thresh
  da %>%
    group_by(L2_joint, study = .data[[study_col]]) %>%
    summarise(
      n         = n(),
      n_sig     = sum(is_sig),
      n_sig_up  = sum(is_sig & logFC > 0),
      n_sig_dn  = sum(is_sig & logFC < 0),
      .groups = "drop"
    ) %>%
    mutate(contrast = contrast_name)
}

# -- iterate contrasts --------------------------------------------------------

l2_long       <- list()
study_long    <- list()
contrast_meta <- list()

for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  da_csv <- file.path(cdir, "da_results.csv")
  if (!file.exists(da_csv)) next
  rs <- yaml::read_yaml(file.path(cdir, "run_summary.yaml"))
  fdr_t <- rs$spatial_fdr %||% 0.05

  da <- read_csv(da_csv, show_col_types = FALSE)
  cat(sprintf("  contrast %-25s nhoods=%d sig=%d (FDR<%.2f)\n",
              cname, nrow(da),
              sum(da$SpatialFDR < fdr_t, na.rm = TRUE), fdr_t))

  l2_long[[cname]]    <- summarise_one(da, cname, fdr_t)
  study_long[[cname]] <- per_study_one(da, cname, fdr_t)

  contrast_meta[[cname]] <- list(
    fdr_thresh = fdr_t,
    cohort     = rs$cohort,
    type       = rs$type,
    coef       = rs$coef,
    contrast_str = rs$contrast_str
  )
}

l2_long_df    <- bind_rows(l2_long)
study_long_df <- bind_rows(study_long)

write_csv(l2_long_df,    file.path(dirname(out_csv), "stageE_per_l2_long.csv"))
write_csv(study_long_df, file.path(dirname(out_csv), "stageE_per_l2_per_study_long.csv"))

# Wide form: one row per L2_joint, columns per (contrast × metric)
wide <- l2_long_df %>%
  pivot_wider(
    id_cols     = c(L2_joint, compartment, label),
    names_from  = contrast,
    values_from = c(n, n_sig, n_sig_up, n_sig_dn, med_lfc, pct_up, pct_dn),
    names_glue  = "{contrast}__{.value}"
  )

# Direct correlation per L2 between paired contrasts (if user defines pairs).
# For the parity_x_risk inquiry, the pairs are e.g. (parity_main,
# parity_x_proph_red). Pairs are declared optionally in inquiry.yaml.
pairs <- inquiry$contrast_pairs %||% list()
corr_rows <- list()
for (p in pairs) {
  m <- p$main; i <- p$interaction
  # join Nhood-level lfc columns
  da_m <- read_csv(file.path(in_dir, m, "da_results.csv"), show_col_types = FALSE)
  da_i <- read_csv(file.path(in_dir, i, "da_results.csv"), show_col_types = FALSE)
  if (!"compartment" %in% colnames(da_m) && "L2_compartment" %in% colnames(da_m))
    da_m$compartment <- da_m$L2_compartment
  if (!"label"       %in% colnames(da_m) && "L2_label"      %in% colnames(da_m))
    da_m$label       <- da_m$L2_label
  if (!"compartment" %in% colnames(da_i) && "L2_compartment" %in% colnames(da_i))
    da_i$compartment <- da_i$L2_compartment
  if (!"label"       %in% colnames(da_i) && "L2_label"      %in% colnames(da_i))
    da_i$label       <- da_i$L2_label
  da_m$L2_joint <- paste(da_m$compartment, da_m$label, sep = "::")
  da_i$L2_joint <- paste(da_i$compartment, da_i$label, sep = "::")
  joined <- da_m %>%
    select(Nhood, lfc_main = logFC, L2_joint) %>%
    inner_join(da_i %>% select(Nhood, lfc_int = logFC),
               by = "Nhood")
  by_l2_corr <- joined %>%
    group_by(L2_joint) %>%
    summarise(
      n_nhoods = n(),
      pearson_r_main_vs_int = suppressWarnings(cor(lfc_main, lfc_int)),
      .groups = "drop"
    ) %>%
    mutate(pair = sprintf("%s__vs__%s", m, i))
  corr_rows[[length(corr_rows) + 1]] <- by_l2_corr
}
if (length(corr_rows) > 0) {
  corr_df <- bind_rows(corr_rows) %>%
    pivot_wider(id_cols = L2_joint,
                names_from = pair, values_from = pearson_r_main_vs_int,
                names_prefix = "r_") %>%
    mutate(across(starts_with("r_"), ~ round(.x, 3)))
  wide <- wide %>% left_join(corr_df, by = "L2_joint")
}

write_csv(wide, out_csv)
cat(sprintf("\nWrote: %s  (rows=%d cols=%d)\n", out_csv, nrow(wide), ncol(wide)))
cat("=== Stage E done ===\n")
