#!/usr/bin/env Rscript
# 33_stageD_da_array.R
# Stage D — DA execution. SLURM array task: ONE contrast per task.
# SLURM_ARRAY_TASK_ID indexes into experiment_metadata.yaml's contrasts list.
#
# For each contrast:
#   1. Load milo (read-only) and patch_milor.R
#   2. Load cohort design CSV
#   3. Build sample-level design data frame matching milo's sample_id column
#   4. Apply joint factor encoding if specified
#   5. testNhoods(milo, design = ~..., model.contrasts = | coef =, ...)
#   6. Annotate L2 + study via Nhood-join (NOT annotateNhoods, which has O(N^2)
#      bug per session memory). Annotation source: milo coldata's nhood index +
#      labels_full.csv.
#   7. Write da_results.csv and run_summary.yaml.
#
# Outputs (per contrast):
#   stageD_da_results/<contrast_name>/da_results.csv
#   stageD_da_results/<contrast_name>/run_summary.yaml

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(dplyr); library(readr); library(yaml); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "33_stageD_da_array.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
source(file.path(script_dir, "lib", "cohort_filters.R"))
source(file.path(script_dir, "lib", "da_helpers.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
if (TASK_ID < 1) stop("SLURM_ARRAY_TASK_ID must be >= 1")

meta <- yaml::read_yaml(file.path(paths$outputs$stageC, "experiment_metadata.yaml"))
if (TASK_ID > length(meta$contrasts)) {
  cat(sprintf("TASK_ID %d > #contrasts %d — exiting\n",
              TASK_ID, length(meta$contrasts)))
  quit(status = 0)
}
con <- meta$contrasts[[TASK_ID]]

cat(sprintf("=== Stage D task %d: %s (cohort=%s) ===\n",
            TASK_ID, con$name, con$cohort))

# Output paths
out_sub <- file.path(paths$outputs$stageD, con$name)
ensure_dir(out_sub)
da_raw_path <- file.path(out_sub, "da_raw.rds")  # post-testNhoods cache

# Cache fast-path (DEV-ONLY, opt-in): if --use-cache is passed AND da_raw.rds
# exists, skip the expensive testNhoods step and resume from annotation.
# Default behavior is ALWAYS to re-run testNhoods so upstream metadata changes
# propagate cleanly. Use --use-cache only when iterating on post-test code.
args_all <- commandArgs(trailingOnly = TRUE)
USE_CACHE <- "--use-cache" %in% args_all
cache_hit <- USE_CACHE &&
             file.exists(da_raw_path) && file.size(da_raw_path) > 1000

source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)

if (cache_hit) {
  cat(sprintf("[DEV] --use-cache: loading da_raw.rds (%.1f MB) — skipping testNhoods\n",
              file.size(da_raw_path) / 1e6))
  da <- readRDS(da_raw_path)
  cat(sprintf("DA results: %d nhoods, %d sig (FDR<%.2f)\n",
              nrow(da), sum(da$SpatialFDR < con$spatial_fdr, na.rm = TRUE),
              con$spatial_fdr))
} else {
  cat(sprintf("milo cells=%d nhoods=%d (running testNhoods%s)\n",
              ncol(milo), ncol(nhoods(milo)),
              if (USE_CACHE) "; --use-cache passed but no cache" else ""))
}

if (!cache_hit) {
# Load cohort design (sample-level). Library column may be sample_id /
# library_id / library / sample / donor_id / ihbca_donor_id depending on
# how milo was countCells'd. Resolve by matching against nhoodCounts colnames.
design_lib <- read_csv(con$cohort_design, show_col_types = FALSE)
nh_sample_ids <- colnames(nhoodCounts(milo))
candidates <- intersect(
  c("sample_id", "library_id", "library", "sample",
    "ihbca_donor_id", "donor_id", "donor"),
  colnames(design_lib))
lib_col <- NA_character_
for (cc in candidates) {
  if (mean(as.character(design_lib[[cc]]) %in% nh_sample_ids) > 0.9) {
    lib_col <- cc; break
  }
}
if (is.na(lib_col)) {
  stop(sprintf("Could not find a design column matching nhoodCounts cols. ",
               "Tried: %s. Sample of nhoodCounts: %s",
               paste(candidates, collapse = ", "),
               paste(head(nh_sample_ids, 3), collapse = ", ")))
}
cat(sprintf("design libraries: %d (lib_col=%s)\n", nrow(design_lib), lib_col))

# Use existing nhoodCounts (populated when milo was originally constructed)
nhood_counts <- nhoodCounts(milo)
if (ncol(nhood_counts) == 0) stop("milo has no nhoodCounts; cannot run testNhoods")
in_cohort <- colnames(nhood_counts) %in% design_lib[[lib_col]]
cat(sprintf("nhood-count libs: total=%d in-cohort=%d\n",
            ncol(nhood_counts), sum(in_cohort)))
if (sum(in_cohort) < 5) stop("Too few libraries after cohort filter")

# Filter design rows to in-cohort libraries (do NOT subset nhoodCounts(milo) —
# replacing the matrix breaks miloR's internal slot alignment with nhoods() /
# nhoodIndex, causing testNhoods rownames<- length errors. Working pattern
# from clinical_da_v4 (1_clinical_da_v4.Rmd:483-497): leave milo intact and
# pass a design.df containing only cohort donors. testNhoods handles the
# alignment internally via rownames(design.df).
cohort_libs_in_milo <- colnames(nhood_counts)[in_cohort]
design_lib <- design_lib %>% filter(.data[[lib_col]] %in% cohort_libs_in_milo)
design_lib <- design_lib[match(cohort_libs_in_milo, design_lib[[lib_col]]), ]

# Do NOT subset nhoodCounts(milo) — see comment above. testNhoods uses
# rownames(design.df) to identify in-cohort libraries.
cat(sprintf("nhoodCounts kept full: %s; cohort libs in design.df: %d\n",
            paste(dim(nhoodCounts(milo)), collapse=" x "),
            nrow(design_lib)))

# testNhoods wants design.df rownames to match nhoodCounts colnames
design_lib <- as.data.frame(design_lib)
rownames(design_lib) <- design_lib[[lib_col]]

# Per-contrast factor level override (v4 pattern). For interaction tests we
# refactor the stratum column so the desired interaction term is the LAST
# column of the design matrix; testNhoods (no model.contrasts) defaults to
# testing the last coef via glmQLFTest(coef=ncol(model)). This sidesteps
# limma::makeContrasts which rejects ":" in design-matrix column names.
if (!is.null(con$stratum_levels)) {
  if (!"stratum" %in% colnames(design_lib)) {
    stop("Contrast '", con$name, "' requires stratum_levels override but ",
         "design has no 'stratum' column.")
  }
  desired_levels <- as.character(con$stratum_levels)
  current_levels <- levels(factor(design_lib$stratum))
  if (!setequal(desired_levels, current_levels)) {
    stop("stratum_levels override mismatch for '", con$name, "'. Yaml: ",
         paste(desired_levels, collapse=","), " | data: ",
         paste(current_levels, collapse=","))
  }
  design_lib$stratum <- factor(design_lib$stratum, levels = desired_levels)
  cat(sprintf("Refactored stratum levels for '%s': %s\n",
              con$name, paste(desired_levels, collapse=",")))
}

# Generic per-contrast factor-level overrides (any column). Format in yaml:
#   factor_levels:
#     age_binary: [young, old]              # young is reference; coef = age_binaryold
#     menopausal_status_binary: ["0", "1"]  # 0 reference; coef = menopausal_status_binary1
# Applied AFTER stratum_levels (which is preserved for back-compat).
if (!is.null(con$factor_levels)) {
  for (col_name in names(con$factor_levels)) {
    if (!col_name %in% colnames(design_lib)) {
      stop("Contrast '", con$name, "' factor_levels references missing col: ",
           col_name)
    }
    desired_levels <- as.character(con$factor_levels[[col_name]])
    current_levels <- as.character(unique(design_lib[[col_name]]))
    current_levels <- current_levels[!is.na(current_levels)]
    if (!all(current_levels %in% desired_levels)) {
      stop("factor_levels mismatch for '", con$name, "' col '", col_name,
           "': yaml=[", paste(desired_levels, collapse=","), "] data=[",
           paste(current_levels, collapse=","), "]")
    }
    design_lib[[col_name]] <- factor(design_lib[[col_name]],
                                      levels = desired_levels)
    cat(sprintf("Refactored '%s' levels for '%s': %s\n",
                col_name, con$name, paste(desired_levels, collapse=",")))
  }
}

# Build model matrix to verify the LAST coef matches con$coef (sanity check).
mm <- model.matrix(as.formula(con$formula), data = design_lib)
last_col <- colnames(mm)[ncol(mm)]
expected_coef <- if (con$type == "coef") con$coef else NA_character_
cat(sprintf("Model matrix: %d cols, last='%s' (expected '%s')\n",
            ncol(mm), last_col, expected_coef %||% "<n/a>"))
if (!is.na(expected_coef) && last_col != expected_coef) {
  stop("Last coef mismatch for '", con$name, "': model.matrix last='",
       last_col, "' but yaml expects '", expected_coef, "'. ",
       "Adjust formula term order or stratum_levels override.")
}
cat(sprintf("Rank: %d / %d\n", qr(mm)$rank, ncol(mm)))

# Run testNhoods WITHOUT model.contrasts — default tests last coef.
# This is the v4 pattern (clinical_da_v4 age_x_brca1) that avoids
# limma::makeContrasts and its colon-in-name validation.
cat(sprintf("Running testNhoods (testing last coef = '%s')...\n", last_col))
da <- testNhoods(
  milo,
  design        = as.formula(con$formula),
  design.df     = design_lib,
  fdr.weighting = con$fdr_weighting,
  reduced.dim   = "scVI_100"
)
cat(sprintf("DA results: %d nhoods, %d sig (FDR<%.2f)\n",
            nrow(da), sum(da$SpatialFDR < con$spatial_fdr, na.rm = TRUE),
            con$spatial_fdr))

# Cache raw da result (post-testNhoods, pre-annotation) so dev iterations
# on post-test code can skip testNhoods via --use-cache.
saveRDS(da, da_raw_path)
cat(sprintf("Cached: %s (%.1f MB)\n",
            da_raw_path, file.size(da_raw_path) / 1e6))

}  # end !cache_hit

# Annotate via Nhood-join (avoid O(N^2) annotateNhoods)
# - For each nhood (= an index cell), look up its compartment/label/study
labels <- read_csv(paths$inputs$labels, show_col_types = FALSE)
cd <- as.data.frame(colData(milo))
nhood_index_cells <- as.integer(nhoodIndex(milo))
nh_cells <- rownames(cd)[nhood_index_cells]
nh_anno <- tibble(
  Nhood          = seq_along(nhood_index_cells),
  index_cell_id  = nh_cells
) %>%
  left_join(labels, by = c("index_cell_id" = "cell_id"))

# study column from milo coldata
study_col <- intersect(c("study", "dataset", "source"), colnames(cd))[1]
if (!is.na(study_col)) {
  nh_anno$index_study <- cd[nhood_index_cells, study_col]
}

da <- da %>% left_join(nh_anno, by = "Nhood")

# Write outputs
out_sub <- file.path(paths$outputs$stageD, con$name)
ensure_dir(out_sub)
write_csv(da, file.path(out_sub, "da_results.csv"))

adaptive_lfc <- compute_adaptive_lfc_cutoff(da, fdr_threshold = con$spatial_fdr)
cat(sprintf("adaptive_lfc_cutoff = %.4f\n", adaptive_lfc))

run_summary <- list(
  contrast            = con$name,
  cohort              = con$cohort,
  formula             = con$formula,
  type                = con$type,
  coef                = con$coef,
  tested_coef         = last_col,
  stratum_levels      = if (!is.null(con$stratum_levels)) as.character(con$stratum_levels) else NULL,
  fdr_weighting       = con$fdr_weighting,
  spatial_fdr         = con$spatial_fdr,
  adaptive_lfc_cutoff = round(adaptive_lfc, 4),
  n_nhoods            = nrow(da),
  n_sig               = sum(da$SpatialFDR < con$spatial_fdr, na.rm = TRUE),
  n_sig_up            = sum(da$SpatialFDR < con$spatial_fdr & da$logFC > 0, na.rm = TRUE),
  n_sig_dn            = sum(da$SpatialFDR < con$spatial_fdr & da$logFC < 0, na.rm = TRUE),
  n_libraries         = sum(in_cohort),
  completed_at        = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write_yaml(run_summary, file.path(out_sub, "run_summary.yaml"))
cat(sprintf("Wrote: %s\n", file.path(out_sub, "da_results.csv")))
cat("=== Stage D task done ===\n")
