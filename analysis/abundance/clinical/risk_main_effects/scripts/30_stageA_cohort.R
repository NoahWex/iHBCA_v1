#!/usr/bin/env Rscript
# 30_stageA_cohort.R
# Stage A — cohort design: build per-cohort library-level design tables
# from milo coldata + harmonized donor metadata, applying per-cohort
# filter rules from inquiry.yaml.
#
# Inputs (resolved via lib/load_paths.R):
#   inputs.milo                — Milo .rds (read for coldata only here)
#   inputs.donor_metadata      — harmonized donor table
#   inquiry.cohorts.<name>     — filter rules + formula + joint_factor (optional)
#
# Outputs (per cohort defined in inquiry.yaml):
#   outputs/stageA_cohorts/cohort_<name>_design.csv
#     One row per LIBRARY (sample_id) passing filters, with all design vars
#     plus joint_factor column if defined.
#   outputs/stageA_cohorts/cohort_summary.md
#     Cross-tabs, exclusion counts, per-study counts.
#
# Usage:
#   Rscript 30_stageA_cohort.R --inquiry-dir /path/to/inquiry

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(dplyr); library(readr); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# -- locate framework lib/ -----------------------------------------------------
this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "30_stageA_cohort.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
source(file.path(script_dir, "lib", "cohort_filters.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry
out_dir <- paths$outputs$stageA
ensure_dir(out_dir)

cat(sprintf("=== Stage A: cohort design for inquiry '%s' ===\n", inquiry$inquiry))

# -- 1. Load milo coldata -> library-level table ------------------------------
cat("Loading milo coldata...\n")
milo <- readRDS(paths$inputs$milo)
cd <- as.data.frame(colData(milo))
cat(sprintf("  cells: %d  cols: %s\n", nrow(cd), paste(colnames(cd), collapse = ", ")))

# Decide donor + study columns
donor_col <- intersect(c("ihbca_donor_id", "donor_id", "donor"), colnames(cd))[1]
if (is.na(donor_col)) stop("No donor identifier column in milo coldata.")
study_col <- intersect(c("study", "dataset", "source"), colnames(cd))[1]

# Library/sample identifier: derive from colnames(nhoodCounts(milo)) if the
# milo was already countCells'd. Otherwise infer the column whose value most
# closely matches.
nh_counts <- tryCatch(nhoodCounts(milo), error = function(e) NULL)
if (!is.null(nh_counts) && ncol(nh_counts) > 0) {
  sample_ids <- colnames(nh_counts)
  cat(sprintf("  nhoodCounts already exists: %d sample cols (head: %s)\n",
              length(sample_ids), paste(head(sample_ids, 3), collapse = ", ")))
  # Find the coldata column whose values are a subset of sample_ids
  match_col <- NA_character_
  for (cc in colnames(cd)) {
    vals <- unique(as.character(cd[[cc]]))
    if (length(vals) > 1 && length(vals) <= length(sample_ids) * 2 &&
        mean(vals %in% sample_ids) > 0.9) {
      match_col <- cc
      break
    }
  }
  if (is.na(match_col)) {
    cat("  WARN: no coldata column matches nhoodCounts colnames; falling back to donor_id\n")
    lib_col <- donor_col
  } else {
    lib_col <- match_col
    cat(sprintf("  inferred library/sample column: %s (matches nhoodCounts)\n", lib_col))
  }
} else {
  lib_col <- intersect(c("sample_id", "library_id", "library", "sample",
                         "dataset", "donor_id", "ihbca_donor_id"),
                       colnames(cd))[1]
  cat(sprintf("  nhoodCounts not present; using inferred library_col=%s\n", lib_col))
}
cat(sprintf("  library_col=%s  donor_col=%s  study_col=%s\n",
            lib_col, donor_col, study_col %||% "NA"))

# Aggregate to library-level: take first non-NA value per library for each
# candidate covariate column (libraries are within-donor, so this is safe).
candidate_cols <- intersect(
  c(lib_col, donor_col, study_col, "dataset_block", "facs_status",
    "tissue_indication", "parity_binary", "menopausal_status_binary",
    "age_binary", "carrier_status", "carrier_binary"),
  colnames(cd)
)
lib_tbl <- cd[, candidate_cols, drop = FALSE] |>
  distinct() |>
  group_by(.data[[lib_col]]) |>
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop")
cat(sprintf("  libraries: %d\n", nrow(lib_tbl)))

# -- 2. Join harmonized donor metadata -----------------------------------------
cat("Loading harmonized donor metadata...\n")
hm <- read_csv(paths$inputs$donor_metadata, show_col_types = FALSE)
cat(sprintf("  donor metadata rows: %d  cols: %s\n",
            nrow(hm), paste(head(colnames(hm), 30), collapse = ", ")))

donor_join_col <- intersect(c("ihbca_donor_id", "donor_id", "donor"), colnames(hm))[1]
if (is.na(donor_join_col)) stop("No donor id col in harmonized donor metadata.")

# Avoid duplicate columns: drop hm columns that already exist in lib_tbl
# (except the join key)
shared <- setdiff(intersect(colnames(lib_tbl), colnames(hm)), donor_join_col)
hm_keep <- hm[, setdiff(colnames(hm), shared), drop = FALSE]

lib_full <- lib_tbl |>
  left_join(hm_keep, by = setNames(donor_join_col, donor_col))
cat(sprintf("  libraries after donor metadata join: %d  cols: %d\n",
            nrow(lib_full), ncol(lib_full)))

# -- 2b. Cell-level (L2) composition table -----------------------------------
# For each cohort downstream: per (L2_joint × study), n_cells + n_libraries.
# Output is cohort-specific, used by Stage G to compute composition-weighted
# collinearities ("X% of this L2's cells come from study Y, where parity-dissoc
# was COLLINEAR — caveat applies").
cat("Loading L2 labels...\n")
labels <- read_csv(paths$inputs$labels, show_col_types = FALSE)
cat(sprintf("  labels rows: %d  cols: %s\n",
            nrow(labels), paste(head(colnames(labels), 10), collapse=", ")))
# Drop artifacts (per labels_full schema)
art_col <- intersect(c("is_artifact", "is_artifact_chr"), colnames(labels))[1]
if (!is.na(art_col)) {
  before <- nrow(labels)
  labels <- labels %>% filter(.data[[art_col]] != "artifact" |
                               is.na(.data[[art_col]]))
  cat(sprintf("  dropped %d artifact rows; %d kept\n", before - nrow(labels), nrow(labels)))
}
# Build cell_id × library × study × compartment × label
cell_col <- intersect(c("cell_id", "cellID"), colnames(labels))[1]
if (is.na(cell_col)) stop("No cell_id column in labels_full.csv")
comp_col <- intersect(c("compartment", "L2_compartment"), colnames(labels))[1]
lab_col  <- intersect(c("label", "L2_label", "L2"), colnames(labels))[1]
labels_min <- labels %>%
  select(all_of(c(cell_col, comp_col, lab_col))) %>%
  rename(cell_id = !!cell_col, compartment = !!comp_col, label = !!lab_col)

cd_cell <- cd %>%
  mutate(cell_id = if ("cellID" %in% colnames(cd)) cellID else rownames(cd))
cd_cell <- cd_cell %>%
  select(cell_id, all_of(c(lib_col, study_col))) %>%
  rename(library = !!lib_col, study = !!study_col) %>%
  inner_join(labels_min, by = "cell_id")
cat(sprintf("  cell × L2 rows after labels join: %d\n", nrow(cd_cell)))

# -- 3. Per-cohort: apply filters, build joint factor, write design ----------
summary_lines <- c(
  sprintf("# Cohort summary — inquiry: %s", inquiry$inquiry),
  "",
  sprintf("Total libraries (pre-filter): %d", nrow(lib_full)),
  sprintf("Total donors (pre-filter): %d", n_distinct(lib_full[[donor_col]])),
  ""
)

for (cohort_name in names(inquiry$cohorts)) {
  spec <- inquiry$cohorts[[cohort_name]]
  cat(sprintf("\n--- cohort %s ---\n", cohort_name))

  filters <- spec$filters
  # Validate filter columns exist
  missing_cols <- setdiff(names(filters), colnames(lib_full))
  if (length(missing_cols) > 0) {
    stop(sprintf("Cohort '%s' filter columns not found: %s",
                 cohort_name, paste(missing_cols, collapse = ", ")))
  }

  keep <- apply_cohort_filters(lib_full, filters)
  d <- lib_full[keep, , drop = FALSE]
  cat(sprintf("  filtered: %d libraries, %d donors\n",
              nrow(d), n_distinct(d[[donor_col]])))

  # Derived columns (e.g., stratum factor combining multiple metadata cols)
  if (!is.null(spec$derived_columns)) {
    d <- apply_derived_columns(d, spec$derived_columns)
    for (dc in names(spec$derived_columns)) {
      cat(sprintf("  derived '%s' levels: %s\n",
                  dc, paste(table(d[[dc]]), collapse=", ")))
    }
  }

  # Joint factor
  if (!is.null(spec$joint_factor)) {
    jf_col <- spec$joint_factor$name
    d[[jf_col]] <- as.character(build_joint_factor(d, spec$joint_factor))
    cat(sprintf("  joint factor '%s' levels: %s\n",
                jf_col, paste(sort(unique(d[[jf_col]])), collapse = ", ")))
  }

  out_file <- file.path(out_dir, sprintf("cohort_%s_design.csv", cohort_name))
  write_csv(d, out_file)
  cat(sprintf("  wrote: %s\n", out_file))

  # -- L2 composition per (compartment::label × study) for this cohort --------
  cohort_libs <- d[[lib_col]]
  cd_in <- cd_cell %>% filter(library %in% cohort_libs)
  l2_comp <- cd_in %>%
    mutate(L2_joint = paste(compartment, label, sep = "::")) %>%
    group_by(L2_joint, compartment, label, study) %>%
    summarise(n_cells = n(),
              n_libraries = n_distinct(library),
              .groups = "drop")
  l2_totals <- l2_comp %>%
    group_by(L2_joint) %>%
    summarise(total_cells = sum(n_cells),
              total_libraries = n_distinct(study),  # studies contributing
              .groups = "drop")
  l2_comp_out <- l2_comp %>%
    left_join(l2_totals, by = "L2_joint") %>%
    mutate(pct_of_l2_cells = round(100 * n_cells / total_cells, 2)) %>%
    arrange(L2_joint, desc(pct_of_l2_cells))

  comp_file <- file.path(out_dir,
                          sprintf("cohort_%s_l2_composition.csv", cohort_name))
  write_csv(l2_comp_out, comp_file)
  cat(sprintf("  wrote: %s  (rows=%d L2s=%d)\n",
              comp_file, nrow(l2_comp_out), n_distinct(l2_comp_out$L2_joint)))

  # Summary section
  summary_lines <- c(summary_lines,
    sprintf("## cohort_%s", cohort_name),
    "",
    sprintf("- description: %s", spec$description %||% ""),
    sprintf("- libraries: %d", nrow(d)),
    sprintf("- donors:    %d", n_distinct(d[[donor_col]])),
    "",
    "### per-study counts",
    "",
    "```",
    capture.output(print(table(study = d[[study_col]], useNA = "ifany"))),
    "```",
    ""
  )

  # Cross-tabs across terms_of_interest
  for (col in inquiry$terms_of_interest) {
    if (!col %in% colnames(d)) next
    tab <- table(d[[col]], useNA = "ifany")
    summary_lines <- c(summary_lines,
      sprintf("### %s", col),
      "",
      "```",
      capture.output(print(tab)),
      "```",
      ""
    )
  }

  # Filter trace (which rule dropped what)
  summary_lines <- c(summary_lines,
    "### filter trace",
    "",
    "```"
  )
  cum_keep <- rep(TRUE, nrow(lib_full))
  for (col in names(filters)) {
    rule_keep <- apply_filter_rule(lib_full[[col]], filters[[col]])
    new_cum <- cum_keep & rule_keep
    summary_lines <- c(summary_lines,
      sprintf("  filter %-30s drops %d -> %d libraries",
              col, sum(cum_keep) - sum(new_cum), sum(new_cum)))
    cum_keep <- new_cum
  }
  summary_lines <- c(summary_lines, "```", "")
}

writeLines(summary_lines, file.path(out_dir, "cohort_summary.md"))
cat(sprintf("\nWrote: %s\n", file.path(out_dir, "cohort_summary.md")))
cat("=== Stage A done ===\n")
