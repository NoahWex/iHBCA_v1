#!/usr/bin/env Rscript
# Render Fig 2 supplemental: donor flow alluvial across cohort-design axes.
#
# Plot structure
#   Six axes (left -> right): study, facs_status, tissue_indication,
#     brca_genotype, risk_class, parity_binary. Co-variate axes (technical
#     nuisance) precede definition axes; the rightmost axes (risk_class +
#     parity_binary) are the testing condition for the PXR design.
#   Strata heights = donor counts per category at each axis. Flows = unique
#     donor-paths through the six axes. Each donor contributes weight 1.
#   facs_status: the legacy "no_sort" label is rewritten to "milk-derived"
#     for display (the remaining no_sort entries are milk samples after the
#     twigger study exclusion).
#   Donors with unknown values at any axis flow through "unknown" strata —
#     not snipped — so the eye reads metadata completeness as stratum width.
#
# Inputs (resolved via --project-root + --inquiry-name + --substrate-dir)
#   <substrate>/config/derived_donor_metadata.csv  (preferred, has risk_class)
#     OR <substrate>/outputs/stageA_cohorts/cohort_*_design.csv  (fallback,
#         deduplicated by donor id and pooled)
#   <substrate>/inquiry.yaml  (optional axis override via alluvial.axes)
#
# Output
#   <out_dir>/s_ihbca_study_alluvials.pdf

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(ggalluvial)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--project-root", required = TRUE,
                    help = "Repository root (resolves to publication/config/)")
parser$add_argument("--inquiry-name", default = "parity_findings",
                    help = "Inquiry directory under publication/analysis/abundance/clinical/")
parser$add_argument("--substrate-dir", default = NULL,
                    help = paste("Override substrate directory for cohort-design",
                                 "metadata; defaults to promoted location"))
parser$add_argument("--out-dir", default = NULL,
                    help = "Output directory (default: <root>/publication/figures/output/supplemental/covariate_landscape/)")
args <- parser$parse_args()

root <- normalizePath(args$project_root)

# ----------------------------------------------------------------------------
# Aesthetics framework
# ----------------------------------------------------------------------------
source(file.path(root, "publication", "config", "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(root, "publication", "config"))

# Alluvial strata share a single 10-color Tol pool (nhoodgroup_rank). Flow
# diagrams read by trajectory, not per-axis category identity; a shared pool
# avoids fragmenting the visual into six independent palettes. Cycled across
# the union of strata categories observed in the data.
tol_pool <- get_palette("nhoodgroup_rank", aes_cfg)

# ----------------------------------------------------------------------------
# Substrate resolution
# ----------------------------------------------------------------------------
substrate_dir <- args$substrate_dir
if (is.null(substrate_dir)) {
  substrate_dir <- file.path(root, "publication", "analysis", "abundance",
                              "clinical", args$inquiry_name, "cohort_design")
}
substrate_dir <- normalizePath(substrate_dir, mustWork = FALSE)

derived_path <- file.path(substrate_dir, "config", "derived_donor_metadata.csv")
stageA_dir   <- file.path(substrate_dir, "outputs", "stageA_cohorts")
inquiry_yaml_cfg <- file.path(substrate_dir, "config", "inquiry.yaml")
inquiry_yaml <- if (file.exists(inquiry_yaml_cfg)) inquiry_yaml_cfg else file.path(substrate_dir, "inquiry.yaml")

# ----------------------------------------------------------------------------
# Pool donor metadata
# ----------------------------------------------------------------------------
if (file.exists(derived_path)) {
  df <- read_csv(derived_path, show_col_types = FALSE)
} else if (dir.exists(stageA_dir)) {
  cohort_files <- list.files(stageA_dir, pattern = "^cohort_.*_design\\.csv$",
                              full.names = TRUE)
  if (length(cohort_files) == 0) {
    stop("No cohort_*_design.csv found under: ", stageA_dir)
  }
  pieces <- lapply(cohort_files, function(fp) {
    suppressWarnings(read_csv(fp, show_col_types = FALSE))
  })
  df <- bind_rows(pieces)
  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "donor", "sample_id"),
                          colnames(df))[1]
  if (!is.na(donor_col)) {
    df <- df %>% distinct(.data[[donor_col]], .keep_all = TRUE)
  }
} else {
  stop("Neither derived_donor_metadata.csv nor stageA_cohorts/ exists under: ",
       substrate_dir)
}

# Synthesize risk_class if absent (same logic as the cohort-design pipeline)
if (!"risk_class" %in% colnames(df)) {
  if (all(c("risk_status_binary", "brca_genotype") %in% colnames(df))) {
    df <- df %>% mutate(risk_class = case_when(
      risk_status_binary == "AR"   ~ "AR",
      brca_genotype     == "BRCA1" ~ "BR1",
      brca_genotype     == "BRCA2" ~ "BR2",
      risk_status_binary == "HR"   ~ "HRS",
      TRUE                          ~ NA_character_
    ))
  }
}

# Normalize genotype column name if upstream uses a synonym
if (!"brca_genotype" %in% colnames(df)) {
  geno_alt <- intersect(c("genotype", "germline_genotype", "germline_status"),
                         colnames(df))[1]
  if (!is.na(geno_alt)) df$brca_genotype <- df[[geno_alt]]
}

# ----------------------------------------------------------------------------
# Axis order — locked design (Noah/PI 2026-05-12)
# Co-variate axes first (technical nuisance is controlled), then definition
# axes ending at the comparison columns (risk_class, parity_binary).
# cancer_history dropped per PI 2026-05-12.
# ----------------------------------------------------------------------------
default_axes <- c("study", "facs_status", "tissue_indication",
                  "brca_genotype", "risk_class", "parity_binary")

ordered_axes <- default_axes
if (file.exists(inquiry_yaml) && requireNamespace("yaml", quietly = TRUE)) {
  inq <- yaml::read_yaml(inquiry_yaml)
  alluvial_cfg <- inq$alluvial %||% list()
  ordered_axes <- alluvial_cfg$axes %||% default_axes
}
axes <- intersect(ordered_axes, colnames(df))

# Display-only label rewrite: facs_status no_sort -> milk-derived. The twigger
# study (the only study tagged no_sort that was not from milk) is excluded
# upstream, so the remaining no_sort entries are milk samples.
if ("facs_status" %in% colnames(df)) {
  df$facs_status <- ifelse(df$facs_status == "no_sort", "milk-derived",
                            df$facs_status)
}

# Coerce missing / blank to "unknown" so the alluvium remains connected for
# every donor; "unknown" is a visible stratum, not a snip.
for (ax in axes) {
  v <- df[[ax]]
  v <- ifelse(is.na(v) | v == "" | v %in% c("NA", "Unknown", "UNKNOWN"),
               "unknown", as.character(v))
  df[[ax]] <- v
}

# Normalize parity_binary to {parous, nulliparous, unknown}
if ("parity_binary" %in% axes) {
  df$parity_binary <- dplyr::recode(as.character(df$parity_binary),
    "p" = "parous", "P" = "parous", "1" = "parous",
    "n" = "nulliparous", "N" = "nulliparous", "0" = "nulliparous",
    "parous" = "parous", "nulliparous" = "nulliparous",
    .default = "unknown")
}

# ----------------------------------------------------------------------------
# Build alluvial frame
# ----------------------------------------------------------------------------
combo <- df %>% count(across(all_of(axes)), name = "n_donors")
combo$alluvium_id <- seq_len(nrow(combo))

long <- combo %>%
  pivot_longer(cols = all_of(axes), names_to = "axis",
                values_to = "category") %>%
  mutate(axis = factor(axis, levels = axes))

# Stratum ordering: sink "unknown" + (legacy "twigger") to the bottom so the
# in-cohort strata cluster at the top within each axis.
sink_cats <- intersect(c("twigger", "unknown"), unique(long$category))
top_cats  <- setdiff(unique(long$category), sink_cats)
long$category <- factor(long$category,
                         levels = c(sink_cats, sort(top_cats)))

# Map every unique category to a Tol-pool color (cycled). Single shared pool
# across all axes — no inline hex, no per-axis palette fragmentation.
unique_cats <- levels(long$category)
fill_vec <- setNames(
  tol_pool[((seq_along(unique_cats) - 1L) %% length(tol_pool)) + 1L],
  unique_cats
)

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
p <- ggplot(long,
             aes(x = axis, stratum = category, alluvium = alluvium_id,
                 y = n_donors, fill = category, label = category)) +
  ggalluvial::geom_flow(stat = "alluvium", alpha = 0.7,
                          lode.guidance = "frontback",
                          decreasing = FALSE,
                          na.rm = TRUE) +
  ggalluvial::geom_stratum(width = 0.4, color = "grey20", na.rm = TRUE,
                            decreasing = FALSE) +
  geom_text(stat = "stratum", size = 2.0, color = "white", na.rm = TRUE,
            decreasing = FALSE) +
  scale_fill_manual(values = fill_vec, guide = "none") +
  labs(x = NULL, y = "n donors") +
  get_theme(aes_cfg) +
  theme(legend.position = "none",
         axis.text.x = element_text(angle = 25, hjust = 1),
         panel.grid = element_blank())

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
out_dir <- args$out_dir %||% file.path(root, "publication", "figures",
                                        "output", "supplemental",
                                        "covariate_landscape")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_id   <- "s_ihbca_study_alluvials"
panel_type <- "alluvial_flow"
n_strata   <- length(unique_cats)

validate_panel(p, panel_type = panel_type, n_groups = n_strata,
               config = aes_cfg)
save_panel(p, panel_id = panel_id, panel_type = panel_type,
           output_dir = out_dir, config = aes_cfg)
