#!/usr/bin/env Rscript
# prep_derive_risk_class.R
# One-shot preprocessor for risk_main_effects inquiry. Reads the canonical
# harmonized donor metadata and emits a derived donor metadata CSV with a
# `risk_class` column added. Pooled cohorts in the inquiry then filter on
# this single derived column.
#
# risk_class derivation (mutually exclusive, complete):
#   AR  : risk_status_binary == "AR"
#   BR1 : brca_genotype == "BRCA1"
#   BR2 : brca_genotype == "BRCA2"
#   HRS : risk_status_binary == "HR" AND brca_genotype %in% {NA, "negative", "RAD51C"-excluded}
#   NA  : everything else (untested HR donors with no genotype info)
#
# Usage:
#   Rscript prep_derive_risk_class.R --inquiry-dir <path>
#
# Output:
#   <inquiry>/config/derived_donor_metadata.csv

suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "prep_derive_risk_class.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

# Canonical donor metadata (read from the path the inquiry's paths.yaml points
# to — for risk_main_effects we override this to be the canonical iHBCA
# publication harmonized_donor_metadata, NOT the derived one we are about to
# write).
canonical_path <- paths$inputs$donor_metadata
if (!file.exists(canonical_path))
  stop(sprintf("Canonical donor metadata not found: %s", canonical_path))

donors <- read_csv(canonical_path, show_col_types = FALSE)
cat(sprintf("Loaded %d donors from canonical metadata.\n", nrow(donors)))

donors <- donors %>%
  mutate(
    risk_class = case_when(
      risk_status_binary == "AR"                                    ~ "AR",
      brca_genotype == "BRCA1"                                      ~ "BR1",
      brca_genotype == "BRCA2"                                      ~ "BR2",
      risk_status_binary == "HR" &
        (is.na(brca_genotype) |
           brca_genotype %in% c("negative", "Negative", "neg"))     ~ "HRS",
      TRUE                                                           ~ NA_character_
    )
  )

cat("\nrisk_class distribution:\n")
print(donors %>% count(risk_class, sort = TRUE))

cat("\nrisk_class × parity_binary cross-tab:\n")
print(donors %>% count(risk_class, parity_binary) %>%
        tidyr::pivot_wider(names_from = parity_binary, values_from = n, values_fill = 0L))

out_path <- file.path(paths$inquiry_root, "config", "derived_donor_metadata.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write_csv(donors, out_path)
cat(sprintf("\nWrote: %s (%d rows)\n", out_path, nrow(donors)))
cat("=== prep_derive_risk_class done ===\n")
