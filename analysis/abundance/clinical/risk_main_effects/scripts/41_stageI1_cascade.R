#!/usr/bin/env Rscript
# 41_stageI1_cascade.R
# Stage I.1 — cascade filter + risk-modifier classification.
#
# For each L2 (compartment::label):
#   - parity_altered_in_AR: |median_logFC_AR| > adaptive_lfc_AR AND n_sig_AR ≥ min
#   - direction_AR: sign(median_logFC) restricted to sig nhoods
#   - l2_signal_type: "clean_direction" (single dominant sign) /
#                     "split"           (~50% up / 50% down) /
#                     "ns_in_AR"
#   - For each non-AR contrast: modifier_pattern ∈
#       {concordant_or_amplified, attenuated, attenuated_to_null,
#        reversed, no_baseline, noise}
#
# Reads:
#   - Stage D da_results.csv per contrast (annotated with L2)
#   - Stage D run_summary.yaml per contrast (adaptive_lfc_cutoff, spatial_fdr)
#
# Output:
#   outputs/stageI1_cascade_table.csv

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(yaml); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "41_stageI1_cascade.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
source(file.path(script_dir, "lib", "da_helpers.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

cat(sprintf("=== Stage I.1: cascade filter for '%s' ===\n", inquiry$inquiry))

# Tunables (could be moved to inquiry.yaml later)
MIN_SIG_NHOODS <- 30        # min sig nhoods to call a population "altered"
SPLIT_PCT_BAND <- 15        # |pct_up - 50| < this → split signal
LFC_RATIO_AMP  <- 1.0       # |non-AR| ≥ this × |AR| → concordant_or_amplified

# Load all contrasts ----------------------------------------------------------
contrast_dirs <- list.dirs(paths$outputs$stageD, recursive = FALSE)
if (length(contrast_dirs) == 0) stop("No Stage D output dirs.")

per_contrast <- list()
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  da_csv <- file.path(cdir, "da_results.csv")
  rs_path <- file.path(cdir, "run_summary.yaml")
  if (!file.exists(da_csv)) next
  da <- read_csv(da_csv, show_col_types = FALSE)
  fdr_t <- 0.05
  adaptive_lfc <- 0.5
  if (file.exists(rs_path)) {
    rs <- yaml::read_yaml(rs_path)
    fdr_t <- rs$spatial_fdr %||% 0.05
    al <- suppressWarnings(as.numeric(rs$adaptive_lfc_cutoff))
    if (length(al) && !is.na(al) && al > 0) adaptive_lfc <- al
  } else {
    al <- compute_adaptive_lfc_cutoff(da, fdr_threshold = fdr_t)
    adaptive_lfc <- if (is.finite(al) && al > 0) al else 0.5
  }
  comp_col <- intersect(c("compartment", "L2_compartment"), colnames(da))[1]
  lab_col  <- intersect(c("label", "L2_label"), colnames(da))[1]
  if (is.na(comp_col) || is.na(lab_col)) next

  da <- da %>%
    rename(L2_compartment = all_of(comp_col), L2_label = all_of(lab_col)) %>%
    filter(!is.na(L2_compartment), !is.na(L2_label)) %>%
    mutate(L2_joint = paste(L2_compartment, L2_label, sep = "::"),
           sig = !is.na(SpatialFDR) & SpatialFDR < fdr_t,
           strong = sig & abs(logFC) > adaptive_lfc)

  l2_summary <- da %>%
    group_by(L2_compartment, L2_label, L2_joint) %>%
    summarise(
      n_nhoods       = n(),
      n_sig          = sum(sig, na.rm = TRUE),
      n_sig_up       = sum(sig & logFC > 0, na.rm = TRUE),
      n_sig_dn       = sum(sig & logFC < 0, na.rm = TRUE),
      n_strong       = sum(strong, na.rm = TRUE),
      median_logFC   = median(logFC[sig], na.rm = TRUE),
      pct_up_of_sig  = ifelse(sum(sig) == 0, NA_real_,
                              100 * sum(sig & logFC > 0, na.rm = TRUE) / sum(sig)),
      .groups = "drop"
    ) %>%
    mutate(contrast = cname,
           adaptive_lfc = adaptive_lfc,
           spatial_fdr = fdr_t,
           direction = ifelse(is.na(median_logFC) | n_sig < 5, "ns",
                       ifelse(median_logFC > 0, "up", "dn")))

  per_contrast[[cname]] <- l2_summary
  cat(sprintf("  %s: %d L2s, adaptive_lfc=%.4f, n_sig_total=%d\n",
              cname, nrow(l2_summary), adaptive_lfc, sum(l2_summary$n_sig)))
}

# Pivot: one row per L2, columns per contrast ---------------------------------
all_l2 <- bind_rows(per_contrast)
all_l2_wide <- all_l2 %>%
  select(L2_compartment, L2_label, L2_joint, contrast,
         n_sig, n_sig_up, n_sig_dn, n_strong,
         median_logFC, pct_up_of_sig, direction, adaptive_lfc) %>%
  pivot_wider(names_from = contrast,
              values_from = c(n_sig, n_sig_up, n_sig_dn, n_strong,
                              median_logFC, pct_up_of_sig, direction,
                              adaptive_lfc))

ar <- "parity_in_AR"
if (!sprintf("n_sig_%s", ar) %in% colnames(all_l2_wide)) {
  stop(sprintf("AR baseline contrast '%s' not found in Stage D outputs.", ar))
}

n_sig_AR     <- all_l2_wide[[sprintf("n_sig_%s", ar)]]
median_AR    <- all_l2_wide[[sprintf("median_logFC_%s", ar)]]
pct_up_AR    <- all_l2_wide[[sprintf("pct_up_of_sig_%s", ar)]]
adapt_AR     <- all_l2_wide[[sprintf("adaptive_lfc_%s", ar)]]

parity_altered_AR <- !is.na(median_AR) &
  n_sig_AR >= MIN_SIG_NHOODS &
  abs(median_AR) > adapt_AR
direction_AR <- ifelse(parity_altered_AR & median_AR > 0, "up",
                ifelse(parity_altered_AR & median_AR < 0, "dn", "ns"))
l2_signal_type <- ifelse(!parity_altered_AR, "ns_in_AR",
                  ifelse(!is.na(pct_up_AR) & abs(pct_up_AR - 50) < SPLIT_PCT_BAND,
                         "split", "clean_direction"))

# Per-non-AR-contrast modifier classification --------------------------------
classify_modifier <- function(med_x, n_sig_x, adapt_x, direction_AR, abs_med_AR) {
  if (direction_AR == "ns") return("no_baseline")
  if (is.na(med_x) | is.na(n_sig_x) | n_sig_x < 5) return("attenuated_to_null")
  sign_x <- ifelse(med_x > 0, "up", ifelse(med_x < 0, "dn", "ns"))
  if (sign_x == "ns") return("attenuated_to_null")
  if (sign_x == direction_AR) {
    if (abs(med_x) >= LFC_RATIO_AMP * abs_med_AR) return("concordant_or_amplified")
    return("attenuated")
  } else {
    if (abs(med_x) >= adapt_x & n_sig_x >= MIN_SIG_NHOODS) return("reversed")
    return("noise")
  }
}

cascade <- all_l2_wide %>%
  mutate(parity_altered_in_AR = parity_altered_AR,
         direction_AR = direction_AR,
         n_sig_AR = n_sig_AR,
         median_logFC_AR = median_AR,
         adaptive_lfc_AR = adapt_AR,
         pct_up_AR = pct_up_AR,
         l2_signal_type = l2_signal_type,
         abs_med_AR = abs(median_AR))

other_contrasts <- setdiff(unique(all_l2$contrast), ar)
for (con in other_contrasts) {
  med_col   <- sprintf("median_logFC_%s", con)
  ns_col    <- sprintf("n_sig_%s", con)
  adapt_col <- sprintf("adaptive_lfc_%s", con)
  cascade[[sprintf("modifier_%s", con)]] <-
    mapply(classify_modifier,
           med_x = cascade[[med_col]] %||% NA,
           n_sig_x = cascade[[ns_col]] %||% NA,
           adapt_x = cascade[[adapt_col]] %||% NA,
           direction_AR = cascade$direction_AR,
           abs_med_AR = cascade$abs_med_AR,
           USE.NAMES = FALSE)
}

cascade_out <- cascade %>%
  select(L2_compartment, L2_label, L2_joint,
         parity_altered_in_AR, direction_AR, l2_signal_type,
         n_sig_AR, median_logFC_AR, adaptive_lfc_AR, pct_up_AR,
         starts_with("modifier_"),
         everything()) %>%
  arrange(desc(parity_altered_in_AR), L2_compartment, desc(abs_med_AR))

out_csv <- paths$outputs$stageI1
write_csv(cascade_out, out_csv)
cat(sprintf("\nWrote: %s\n", out_csv))
cat(sprintf("L2s passing AR-altered filter: %d / %d\n",
            sum(parity_altered_AR), length(parity_altered_AR)))
cat(sprintf("  clean_direction: %d, split: %d\n",
            sum(l2_signal_type == "clean_direction" & parity_altered_AR),
            sum(l2_signal_type == "split" & parity_altered_AR)))

# Quick console summary of modifier patterns
for (con in other_contrasts) {
  col <- sprintf("modifier_%s", con)
  patterns <- table(cascade_out[parity_altered_AR, ][[col]], useNA = "ifany")
  cat(sprintf("\n  %s modifier pattern (within AR-altered L2s):\n", con))
  for (pat in names(patterns)) {
    cat(sprintf("    %-25s %d\n", pat, patterns[[pat]]))
  }
}

cat("\n=== Stage I.1 done ===\n")
