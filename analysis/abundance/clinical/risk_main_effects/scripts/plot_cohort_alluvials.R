#!/usr/bin/env Rscript
# plot_cohort_alluvials.R
#
# Pooled donor-flow alluvial showing all donors across the inquiry's testable
# axes. Reveals confounding structure across study, genotype, technical
# covariates, and risk-class assignment. Parity is the rightmost (terminal)
# tier so the reader's eye lands on the testing condition.
#
# Source data preference (in order):
#   1. config/derived_donor_metadata.csv  (preferred — has unified risk_class)
#   2. pooled stageA_cohorts/cohort_*_design.csv (deduplicated by donor id)
#
# Output: outputs/stageB_exploration/alluvials/alluvial_pooled.pdf
#
# Locked axis order (Noah review 2026-05-07):
#   study → genotype → facs_status → cancer_history → risk_class → parity_binary
# Earlier per-cohort alluvials are dropped per the same review.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})
have_alluvial <- requireNamespace("ggalluvial", quietly = TRUE)
if (have_alluvial) suppressPackageStartupMessages(library(ggalluvial))
have_patchwork <- requireNamespace("patchwork", quietly = TRUE)
if (have_patchwork) suppressPackageStartupMessages(library(patchwork))

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_cohort_alluvials.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inq_root <- paths$inquiry_root

out_dir <- file.path(inq_root, "outputs", "stageB_exploration", "alluvials")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Axis order: read from inquiry.yaml's `alluvial.axes` block if present,
# otherwise fall back to the parity-x-risk locked default.
# The rightmost axis is the testing condition (parity in PXR/RM; age in AXM).
inquiry <- cfg$inquiry
alluvial_cfg <- inquiry$alluvial %||% list()
ordered_axes <- alluvial_cfg$axes %||%
  c("study", "facs_status", "brca_genotype",
    "cancer_history", "risk_class", "parity_binary")
exclusion_rules <- alluvial_cfg$exclusion_priority %||% list(
  list(field = "study",          value = "twigger", label = "excl: twigger study"),
  list(field = "parity_binary",  value = "unknown", label = "excl: parity unknown"),
  list(field = "risk_class",     value = "unknown", label = "excl: risk_class unknown")
)
cat(sprintf("Alluvial axes (from %s): %s\n",
            if (is.null(inquiry$alluvial)) "default" else "inquiry.yaml",
            paste(ordered_axes, collapse = " -> ")))

# ---- Pool donors ----
derived_path <- file.path(inq_root, "config", "derived_donor_metadata.csv")
if (file.exists(derived_path)) {
  cat(sprintf("Using derived donor metadata: %s\n", derived_path))
  df <- read_csv(derived_path, show_col_types = FALSE)
} else {
  cat("derived_donor_metadata.csv missing — pooling cohort designs\n")
  stageA_dir <- file.path(inq_root, "outputs", "stageA_cohorts")
  cohort_files <- list.files(stageA_dir, pattern = "^cohort_.*_design\\.csv$",
                              full.names = TRUE)
  pieces <- lapply(cohort_files, function(fp) {
    suppressWarnings(read_csv(fp, show_col_types = FALSE))
  })
  df <- bind_rows(pieces)
  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "donor", "sample_id"),
                          colnames(df))[1]
  if (!is.na(donor_col)) df <- df %>% distinct(.data[[donor_col]], .keep_all = TRUE)
  cat(sprintf("Pooled %d cohort files -> %d unique donors\n",
              length(cohort_files), nrow(df)))
}

# Synthesize risk_class if missing (same logic as prep_derive_risk_class.R)
if (!"risk_class" %in% colnames(df)) {
  if (all(c("risk_status_binary", "brca_genotype") %in% colnames(df))) {
    df <- df %>% mutate(risk_class = case_when(
      risk_status_binary == "AR"          ~ "AR",
      brca_genotype     == "BRCA1"        ~ "BR1",
      brca_genotype     == "BRCA2"        ~ "BR2",
      risk_status_binary == "HR"          ~ "HRS",
      TRUE                                  ~ NA_character_
    ))
    cat("Synthesized risk_class column from risk_status_binary + brca_genotype.\n")
  }
}

# Resolve genotype column name flexibly
if (!"brca_genotype" %in% colnames(df)) {
  geno_alt <- intersect(c("genotype", "germline_genotype", "germline_status"),
                         colnames(df))[1]
  if (!is.na(geno_alt)) {
    df$brca_genotype <- df[[geno_alt]]
    cat(sprintf("Mapped genotype column: %s -> brca_genotype\n", geno_alt))
  }
}

# Final axis set: keep only those present
axes <- intersect(ordered_axes, colnames(df))
cat(sprintf("Axes (left -> right): %s\n", paste(axes, collapse = " -> ")))

# Coerce missing / blank to "unknown" (visible label) for stratum
for (ax in axes) {
  v <- df[[ax]]
  v <- ifelse(is.na(v) | v == "" | v %in% c("NA", "Unknown", "UNKNOWN"),
               "unknown", as.character(v))
  df[[ax]] <- v
}

# Normalize parity_binary explicitly
if ("parity_binary" %in% axes) {
  df$parity_binary <- dplyr::recode(as.character(df$parity_binary),
    "p" = "parous", "P" = "parous", "1" = "parous",
    "n" = "nulliparous", "N" = "nulliparous", "0" = "nulliparous",
    "parous" = "parous", "nulliparous" = "nulliparous",
    .default = "unknown")
}

# NO snipping (revised 2026-05-07 v3): "unknown" at any tier is a legitimate
# stratum, not an exclusion. Genotype-untested Kumar donors still flow through
# to risk_class (AR vs HRS by cancer-history alone). All donors stay in the
# alluvium across all tiers; "unknown" is a visible stratum width.

combo <- df %>% count(across(all_of(axes)), name = "n_donors")
combo$alluvium_id <- seq_len(nrow(combo))
cat(sprintf("Pooled alluvial: %d donors, %d unique node-paths\n",
            sum(combo$n_donors), nrow(combo)))

# Per-donor inclusion status: drives the bottom-row exclusion panel.
# Rules come from inquiry.yaml's alluvial.exclusion_priority (first-match wins;
# rules iterate in reverse so the first rule's hit is the final label).
df$analysis_status <- "in analysis"
for (i in rev(seq_along(exclusion_rules))) {
  rule <- exclusion_rules[[i]]
  fld <- rule$field; val <- rule$value
  if (!fld %in% colnames(df)) next
  hit <- df[[fld]] == val | (val == "unknown" & is.na(df[[fld]]))
  df$analysis_status[hit] <- rule$label
}
status_levels <- c("in analysis",
                    sapply(exclusion_rules, function(r) r$label))
df$analysis_status <- factor(df$analysis_status,
                              levels = intersect(status_levels,
                                                  unique(df$analysis_status)))
status_pal_default <- c("#117733", "#CC6677", "#DDCC77", "#882255",
                         "#88CCEE", "#332288", "#999933", "#AA4499")
status_pal <- setNames(status_pal_default[seq_along(status_levels)],
                        status_levels)
status_summary <- df %>% count(analysis_status, name = "n_donors")
cat("Donor inclusion breakdown:\n"); print(status_summary)

# ---- Render ----
if (have_alluvial) {
  long <- combo %>%
    pivot_longer(cols = all_of(axes), names_to = "axis",
                  values_to = "category") %>%
    mutate(axis = factor(axis, levels = axes))

  # Strata ordering within each axis: exclusion-causing categories sink to
  # the bottom of the stack. ggalluvial's geom_stratum + geom_flow with
  # decreasing=FALSE respects factor-level order (first level = bottom).
  exclusion_cats <- intersect(c("twigger", "unknown"),
                               unique(long$category))
  non_excl       <- setdiff(unique(long$category), exclusion_cats)
  long$category  <- factor(long$category,
                            levels = c(exclusion_cats, sort(non_excl)))
  p <- ggplot(long,
               aes(x = axis, stratum = category, alluvium = alluvium_id,
                   y = n_donors, fill = category, label = category)) +
    ggalluvial::geom_flow(stat = "alluvium", alpha = 0.7,
                            lode.guidance = "frontback",
                            decreasing = FALSE,
                            na.rm = TRUE) +
    ggalluvial::geom_stratum(width = 0.4, color = "grey20", na.rm = TRUE,
                              decreasing = FALSE) +
    geom_text(stat = "stratum", size = 2.4, color = "black", na.rm = TRUE,
              decreasing = FALSE) +
    labs(x = NULL, y = "n donors") +
    theme_bw(base_size = 9) +
    theme(legend.position = "none",
           axis.text.x = element_text(angle = 25, hjust = 1, size = 8),
           panel.grid = element_blank(),
           panel.background = element_blank())
} else {
  p <- combo %>%
    mutate(combo_label = apply(.[, axes], 1,
                                function(r) paste(r, collapse = " | "))) %>%
    ggplot(aes(x = reorder(combo_label, -n_donors), y = n_donors)) +
    geom_col() + coord_flip() + theme_bw(base_size = 8) +
    labs(x = NULL, y = "n donors")
}

pdf_path <- file.path(out_dir, "alluvial_pooled.pdf")
ggsave(pdf_path, p, width = 11.5, height = 5.5,
        device = cairo_pdf, limitsize = FALSE)
cat(sprintf("Wrote %s\n", pdf_path))

# Drop superseded per-cohort alluvials (keep only pooled)
all_pdfs <- list.files(out_dir, pattern = "\\.pdf$", full.names = TRUE)
old <- all_pdfs[basename(all_pdfs) != "alluvial_pooled.pdf"]
if (length(old) > 0) {
  cat(sprintf("Removing %d superseded per-cohort alluvial(s)\n", length(old)))
  file.remove(old)
}

cat("=== plot_cohort_alluvials done (pooled) ===\n")
