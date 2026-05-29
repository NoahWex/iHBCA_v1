#!/usr/bin/env Rscript
# 32_stageC_contrasts.R
# Stage C — contrast specification: validate inquiry.yaml's contrasts against
# the cohort designs produced by Stage A, then emit a flattened
# experiment_metadata.yaml that Stage D's array tasks consume directly.
#
# Each row in experiment_metadata.yaml's `contrasts:` list is one DA test:
#   {name, cohort, formula, type, coef|contrast, fdr_weighting,
#    spatial_fdr, prop, design_csv}
#
# Validation:
#   - cohort design CSV exists for each referenced cohort
#   - formula refers only to columns present in design CSV
#   - tested coef is a valid coefficient name when the formula is built
#     against the cohort design (lm(model.matrix(formula, design))) — not run
#     here, since the design factor levels need the joint-factor levels;
#     instead we just confirm referenced terms exist as columns.

suppressPackageStartupMessages({
  library(yaml); library(readr); library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "32_stageC_contrasts.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry
in_dir <- paths$outputs$stageA
out_dir <- paths$outputs$stageC
ensure_dir(out_dir)

cat(sprintf("=== Stage C: contrast specification for '%s' ===\n", inquiry$inquiry))

flat <- list()
errors <- character()

for (cohort_name in names(inquiry$cohorts)) {
  spec <- inquiry$cohorts[[cohort_name]]
  design_csv <- file.path(in_dir, sprintf("cohort_%s_design.csv", cohort_name))
  if (!file.exists(design_csv)) {
    errors <- c(errors, sprintf("cohort %s: design CSV missing (%s)",
                                cohort_name, design_csv))
    next
  }
  d <- read_csv(design_csv, show_col_types = FALSE, n_max = 5)
  cols <- colnames(d)
  formula_str <- spec$formula
  formula_terms <- all.vars(as.formula(formula_str))
  missing_terms <- setdiff(formula_terms, cols)
  if (length(missing_terms) > 0) {
    errors <- c(errors, sprintf("cohort %s: formula references missing columns: %s",
                                cohort_name, paste(missing_terms, collapse = ", ")))
  }

  for (con in spec$contrasts) {
    # Per-contrast formula override allows two contrasts to share a cohort but
    # test different last coefs (different formula term order).
    con_formula <- con$formula %||% formula_str
    if (!identical(con_formula, formula_str)) {
      con_terms <- all.vars(as.formula(con_formula))
      con_missing <- setdiff(con_terms, cols)
      if (length(con_missing) > 0) {
        errors <- c(errors, sprintf("contrast %s: formula references missing columns: %s",
                                    con$name, paste(con_missing, collapse = ", ")))
      }
    }
    flat[[length(flat) + 1]] <- list(
      name           = con$name,
      cohort         = cohort_name,
      cohort_design  = design_csv,
      formula        = con_formula,
      type           = con$type,
      coef           = con$coef %||% NA_character_,
      contrast       = con$contrast %||% NA_character_,
      stratum_levels = con$stratum_levels,
      factor_levels  = con$factor_levels,
      fdr_weighting  = spec$fdr_weighting %||% "graph-overlap",
      spatial_fdr    = spec$spatial_fdr %||% 0.05,
      prop           = spec$prop %||% 0.1,
      joint_factor   = spec$joint_factor %||% NA,
      description    = con$description %||% ""
    )
    cat(sprintf("  contrast %-25s cohort=%s type=%s\n",
                con$name, cohort_name, con$type))
  }
}

if (length(errors) > 0) {
  cat("\nVALIDATION ERRORS:\n")
  for (e in errors) cat(" -", e, "\n")
  stop("Stage C validation failed; fix inquiry.yaml or rerun Stage A.")
}

meta <- list(
  inquiry = inquiry$inquiry,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  contrasts = flat
)
write_yaml(meta, file.path(out_dir, "experiment_metadata.yaml"))
cat(sprintf("\nWrote: %s  (%d contrasts)\n",
            file.path(out_dir, "experiment_metadata.yaml"), length(flat)))
cat("=== Stage C done ===\n")
