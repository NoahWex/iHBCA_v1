#!/usr/bin/env Rscript
# 40_stageI0_per_study.R
# Stage I.0 — per-study breakdown of per-L2 DA metrics. Transparency layer.
# For every (contrast × L2 × study), reports n_donors, n_cells, n_nhoods_in_l2,
# n_sig (FDR<spatial_fdr), n_sig_up, n_sig_dn, median_logFC, dominant_direction,
# study_share_of_l2 (% of L2 cells contributed by this study within the cohort).
#
# Reads:
#   - milo + labels for cell × study × L2 mapping
#   - cohort design CSVs for donor membership per cohort
#   - Stage D da_results.csv per contrast (annotated with compartment+label
#     by Nhood-join in Stage D)
#   - Stage D run_summary.yaml per contrast (for spatial_fdr)
#
# Output:
#   outputs/stageI0_per_study_breakdown.csv

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "40_stageI0_per_study.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

cat(sprintf("=== Stage I.0: per-study breakdown for '%s' ===\n", inquiry$inquiry))

# Load milo + labels for cell×study×L2 ---------------------------------------
source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)
cd <- as.data.frame(colData(milo))
cd$cell_id <- colnames(milo)
study_col <- intersect(c("study", "dataset", "source"), colnames(cd))[1]
if (is.na(study_col)) stop("No study/dataset column in milo coldata.")
donor_col <- intersect(c("ihbca_donor_id", "sample_id", "library_id"),
                       colnames(cd))[1]
if (is.na(donor_col)) stop("No donor identifier column in milo coldata.")
cd <- cd %>% select(cell_id, all_of(c(study_col, donor_col))) %>%
  rename(study = all_of(study_col), donor = all_of(donor_col))

labels <- read_csv(paths$inputs$labels, show_col_types = FALSE)
lab_join <- labels %>%
  select(cell_id,
         L2_compartment = any_of(c("compartment", "L2_compartment")),
         L2_label = any_of(c("label", "L2_label")))
cd <- cd %>% left_join(lab_join, by = "cell_id") %>%
  filter(!is.na(L2_compartment), !is.na(L2_label)) %>%
  mutate(L2_joint = paste(L2_compartment, L2_label, sep = "::"))

cat(sprintf("milo cells with L2 labels: %d\n", nrow(cd)))

# Per-contrast loop ----------------------------------------------------------
contrast_dirs <- list.dirs(paths$outputs$stageD, recursive = FALSE)
if (length(contrast_dirs) == 0) stop("No Stage D output dirs.")

all_rows <- list()
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  da_csv <- file.path(cdir, "da_results.csv")
  rs_path <- file.path(cdir, "run_summary.yaml")
  if (!file.exists(da_csv)) next
  da <- read_csv(da_csv, show_col_types = FALSE)
  fdr_t <- 0.05
  cohort_name <- NA_character_
  if (file.exists(rs_path)) {
    rs <- yaml::read_yaml(rs_path)
    fdr_t <- rs$spatial_fdr %||% 0.05
    cohort_name <- rs$cohort %||% NA_character_
  }

  comp_col <- intersect(c("compartment", "L2_compartment"), colnames(da))[1]
  lab_col  <- intersect(c("label", "L2_label"), colnames(da))[1]
  if (is.na(comp_col) || is.na(lab_col)) {
    cat(sprintf("  %s: missing L2 cols, skipping\n", cname)); next
  }

  da <- da %>%
    rename(L2_compartment = all_of(comp_col), L2_label = all_of(lab_col)) %>%
    filter(!is.na(L2_compartment), !is.na(L2_label)) %>%
    mutate(L2_joint = paste(L2_compartment, L2_label, sep = "::"),
           sig = !is.na(SpatialFDR) & SpatialFDR < fdr_t)

  # Cohort design: which donors in this contrast's cohort
  cohort_csv <- if (!is.na(cohort_name))
    file.path(paths$outputs$stageA,
               sprintf("cohort_%s_design.csv", cohort_name)) else NA
  cohort_donors <- NULL
  if (!is.na(cohort_csv) && file.exists(cohort_csv)) {
    cd_design <- read_csv(cohort_csv, show_col_types = FALSE)
    donor_field <- intersect(c("ihbca_donor_id", "donor", "sample_id"),
                              colnames(cd_design))[1]
    cohort_donors <- if (!is.na(donor_field)) unique(cd_design[[donor_field]]) else NULL
  }

  # Cells in cohort (for per-L2 study composition)
  if (!is.null(cohort_donors)) {
    cd_in_cohort <- cd %>% filter(donor %in% cohort_donors)
  } else {
    cd_in_cohort <- cd
  }

  # Per-L2 cell × study counts (for study_share)
  l2_study_cells <- cd_in_cohort %>%
    count(L2_joint, study, name = "n_cells")
  l2_total_cells <- l2_study_cells %>%
    group_by(L2_joint) %>%
    summarise(total_cells = sum(n_cells), .groups = "drop")

  # Per-L2 cell × study donor counts
  l2_study_donors <- cd_in_cohort %>%
    distinct(L2_joint, study, donor) %>%
    count(L2_joint, study, name = "n_donors")

  # Per-L2 × study DA metrics: requires nhood→study mapping. Each nhood's index
  # cell carries a study label (from milo coldata). Stage D's da_results
  # currently does not carry per-nhood study; recompute from nhood index.
  nh_idx <- nhoodIndex(milo)  # numeric vector: nhood_id → cell_index
  if (length(nh_idx) != nrow(da)) {
    cat(sprintf("  %s: nhoodIndex length=%d != da rows=%d, skipping\n",
                cname, length(nh_idx), nrow(da))); next
  }
  da$index_cell <- colnames(milo)[nh_idx]
  cell_study <- setNames(cd$study, cd$cell_id)
  da$index_study <- cell_study[da$index_cell]

  da_per_l2_study <- da %>%
    filter(!is.na(index_study)) %>%
    group_by(L2_joint, L2_compartment, L2_label, index_study) %>%
    summarise(n_nhoods = n(),
              n_sig = sum(sig, na.rm = TRUE),
              n_sig_up = sum(sig & logFC > 0, na.rm = TRUE),
              n_sig_dn = sum(sig & logFC < 0, na.rm = TRUE),
              median_logFC = median(logFC, na.rm = TRUE),
              .groups = "drop") %>%
    rename(study = index_study) %>%
    mutate(direction = ifelse(n_sig == 0, "ns",
                       ifelse(n_sig_up > n_sig_dn * 1.5, "up",
                       ifelse(n_sig_dn > n_sig_up * 1.5, "dn", "mixed"))))

  # Join study composition
  out <- da_per_l2_study %>%
    left_join(l2_study_cells, by = c("L2_joint", "study")) %>%
    left_join(l2_total_cells, by = "L2_joint") %>%
    left_join(l2_study_donors, by = c("L2_joint", "study")) %>%
    mutate(study_share = ifelse(is.na(total_cells) | total_cells == 0, NA,
                                  round(100 * n_cells / total_cells, 2)),
           contrast = cname,
           cohort = cohort_name,
           spatial_fdr = fdr_t)

  all_rows[[cname]] <- out
  cat(sprintf("  %s: %d (L2 × study) rows\n", cname, nrow(out)))
}

if (length(all_rows) == 0) stop("No per-study rows produced.")

result <- bind_rows(all_rows) %>%
  select(contrast, cohort, spatial_fdr,
         L2_compartment, L2_label, L2_joint, study,
         n_donors, n_cells, study_share,
         n_nhoods, n_sig, n_sig_up, n_sig_dn, median_logFC, direction)

out_csv <- paths$outputs$stageI0
write_csv(result, out_csv)
cat(sprintf("\nWrote: %s (%d rows, %d contrasts × ~%d L2 × ~%d studies)\n",
            out_csv, nrow(result),
            length(unique(result$contrast)),
            length(unique(result$L2_joint)),
            length(unique(result$study))))

cat("=== Stage I.0 done ===\n")
