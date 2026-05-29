#!/usr/bin/env Rscript
# stageI2_gap_audit_risk.R
#
# Risk-side mirror of stageI2_gap_audit.R. Audits which (L2 x contrast) cells
# have L2-grain DA signal but no viable NhoodGroup carving — i.e., gap
# candidates where direct L2-grain DE (Stage F.3 L2 markers) is the only
# gene-level evidence path.
#
# Difference from parity version:
#   - parity has an AR-baseline reference (parity_in_AR) against which all
#     non-AR contrasts are compared. risk_main_effects has no such anchor;
#     each contrast (BR1_vs_AR, BR2_vs_AR, HRS_vs_AR) is a main-effect test
#     in its own right. Gap detection is per-contrast, no cross-contrast
#     comparison.
#
# Inputs:
#   outputs/stageE_per_l2_summary.csv      (per L2 x contrast)
#   outputs/nhoodgroup_renaming/lookup.csv (viable groups per L2 x contrast)
#
# Output:
#   outputs/stageI2_gap_audit.csv          (one row per (L2 x contrast))

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI2_gap_audit_risk.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

SIG_THRESHOLD <- 30  # n_sig threshold for "L2 has DA signal in this contrast"

# Stage E per-L2 summary (wide format expected: L2 x contrast metrics)
stE <- read_csv(file.path(paths$inquiry_root, "outputs",
                          "stageE_per_l2_summary.csv"),
                show_col_types = FALSE)

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

audit <- stE_long %>%
  left_join(ng_counts, by = c("L2_joint", "contrast")) %>%
  mutate(n_viable_NG = coalesce(n_viable_NG, 0L),
         n_sig_NG    = coalesce(n_sig_NG, 0L),
         has_signal = !is.na(n_sig) & n_sig >= SIG_THRESHOLD,
         gap_candidate = has_signal & n_viable_NG == 0L,
         gap_severity = case_when(
           !has_signal                ~ "no_L2_signal",
           n_viable_NG == 0L          ~ "no_viable_NG__needs_direct_DE",
           n_sig_NG == 0L             ~ "viable_NG_no_signal__skipping_DE_ok",
           n_viable_NG == 1L          ~ "single_NG__limited_resolution",
           TRUE                       ~ "carved_well"
         )) %>%
  select(any_of(c("compartment", "L2_joint", "label", "contrast",
                  "n", "n_sig", "n_sig_up", "n_sig_dn",
                  "med_lfc", "pct_up", "pct_dn",
                  "n_viable_NG", "n_sig_NG",
                  "has_signal", "gap_candidate", "gap_severity"))) %>%
  arrange(desc(gap_candidate), compartment, L2_joint, contrast)

out_csv <- file.path(paths$inquiry_root, "outputs", "stageI2_gap_audit.csv")
write_csv(audit, out_csv)
cat(sprintf("Wrote %s (%d rows, %d gap candidates)\n",
            out_csv, nrow(audit), sum(audit$gap_candidate, na.rm = TRUE)))

cat("\n--- Gap severity distribution by contrast ---\n")
audit %>%
  count(contrast, gap_severity) %>%
  pivot_wider(names_from = gap_severity, values_from = n, values_fill = 0L) %>%
  print(n = Inf)

cat("\n--- Gap candidates (signal at L2 grain, 0 viable NGs) ---\n")
audit %>%
  filter(gap_candidate) %>%
  select(compartment, L2_joint, contrast, n_sig, med_lfc) %>%
  print(n = Inf)
