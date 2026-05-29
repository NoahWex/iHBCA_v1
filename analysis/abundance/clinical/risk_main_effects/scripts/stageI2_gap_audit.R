#!/usr/bin/env Rscript
# stageI2_gap_audit.R
#
# Audit which L2s have parity_in_AR signal but no viable NhoodGroup in
# non-AR contrasts. These are gap candidates that would need direct
# L2-grain DE (Stage I.2 cross-cohort runner) rather than NhoodGroup-grain
# F.3 markers, because their carving doesn't carry into the modification
# contrasts.
#
# Inputs:
#   outputs/stageE_per_l2_summary.csv      (per L2 x contrast: n_sig, med_lfc)
#   outputs/nhoodgroup_renaming/lookup.csv (viable groups per L2 x contrast)
#
# Output:
#   outputs/stageI2_gap_audit.csv          (one row per (L2 x non-AR contrast)
#                                           with flags + counts)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI2_gap_audit.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

AR_CONTRAST <- "parity_in_AR"
SIG_THRESHOLD <- 30  # n_sig threshold for "L2 has parity_in_AR signal"

# Stage E per-L2 summary (wide format)
stE <- read_csv(file.path(paths$inquiry_root, "outputs",
                           "stageE_per_l2_summary.csv"),
                 show_col_types = FALSE)

# Pivot to long: one row per (L2 x contrast)
metric_cols <- grep("^[A-Za-z_0-9]+__[a-z_]+$", colnames(stE), value = TRUE)
stE_long <- stE %>%
  select(L2_joint, compartment, label, all_of(metric_cols)) %>%
  pivot_longer(all_of(metric_cols),
               names_to = c("contrast", "metric"),
               names_pattern = "(.*)__(.*)") %>%
  pivot_wider(names_from = metric, values_from = value)

# Viable NhoodGroup count per (L2 x contrast)
lookup <- read_csv(file.path(paths$inquiry_root, "outputs",
                              "nhoodgroup_renaming", "lookup.csv"),
                    show_col_types = FALSE) %>%
  filter(viable)

ng_counts <- lookup %>%
  group_by(parent_L2_joint, contrast) %>%
  summarise(n_viable_NG = dplyr::n(),
            n_sig_NG = sum(n_sig >= 5, na.rm = TRUE),
            .groups = "drop") %>%
  rename(L2_joint = parent_L2_joint)

# Join Stage E + viable NG counts
audit <- stE_long %>%
  left_join(ng_counts, by = c("L2_joint", "contrast")) %>%
  mutate(n_viable_NG = coalesce(n_viable_NG, 0L),
         n_sig_NG    = coalesce(n_sig_NG, 0L))

# Identify L2s with parity_in_AR signal
ar_signal <- audit %>%
  filter(contrast == AR_CONTRAST) %>%
  mutate(has_AR_signal = !is.na(n_sig) & n_sig >= SIG_THRESHOLD) %>%
  select(L2_joint, has_AR_signal,
         AR_n_sig = n_sig, AR_med_lfc = med_lfc, AR_n_viable_NG = n_viable_NG)

# Flag gap candidates: AR has signal AND target contrast has 0 viable NG
gap_table <- audit %>%
  filter(contrast != AR_CONTRAST) %>%
  left_join(ar_signal, by = "L2_joint") %>%
  mutate(
    gap_candidate = has_AR_signal & n_viable_NG == 0L,
    gap_severity = case_when(
      !has_AR_signal             ~ "no_AR_baseline",
      n_viable_NG == 0L          ~ "no_viable_NG__needs_direct_DE",
      n_sig_NG == 0L             ~ "viable_NG_no_signal__skipping_DE_ok",
      n_viable_NG == 1L          ~ "single_NG__limited_resolution",
      TRUE                        ~ "carved_well"
    )
  ) %>%
  select(any_of(c("compartment", "L2_joint", "label", "contrast",
                   "n", "n_sig", "n_sig_up", "n_sig_dn",
                   "med_lfc", "pct_up", "pct_dn",
                   "n_viable_NG", "n_sig_NG",
                   "AR_n_sig", "AR_med_lfc", "AR_n_viable_NG",
                   "has_AR_signal", "gap_candidate", "gap_severity"))) %>%
  arrange(desc(gap_candidate), compartment, L2_joint, contrast)

out_csv <- file.path(paths$inquiry_root, "outputs",
                     "stageI2_gap_audit.csv")
write_csv(gap_table, out_csv)
cat(sprintf("Wrote %s (%d rows, %d gap candidates)\n",
            out_csv, nrow(gap_table), sum(gap_table$gap_candidate, na.rm = TRUE)))

# Summary by contrast x severity
cat("\n--- Gap severity distribution by contrast ---\n")
gap_table %>%
  count(contrast, gap_severity) %>%
  pivot_wider(names_from = gap_severity, values_from = n, values_fill = 0L) %>%
  print(n = Inf)

# List the actual gap candidates (most actionable)
cat("\n--- Gap candidates (AR has signal, contrast has 0 viable NG) ---\n")
gap_table %>%
  filter(gap_candidate) %>%
  select(compartment, L2_joint, contrast, AR_n_sig, AR_med_lfc, n_sig, med_lfc) %>%
  print(n = Inf)
