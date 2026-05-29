#!/usr/bin/env Rscript
#
# Step 15 — Milo Contamination (canonical port, frozen-handoff mode)
#
# Step 15 is THE biological-artifact removal step that was absent from the
# prior FLEX preprocessing pipeline. It uses Milo differential abundance
# testing on position-matched contrasts to identify neighborhoods enriched
# in contaminated samples (DermEpithelial-bearing positions from step 14's
# epidermal_samples.yaml), then propagates neighborhood-level flags down
# to ~9.9K individual cells — almost all keratinocytes and melanocytes that
# leaked into the breast atlas during dissection.
#
# Dev-tree execution plan (preserved below for provenance):
#
#   Stage 1  1_build_milo_object.R            Milo graph + neighborhoods on
#                                             the step 12 scVI embedding.
#   Stage 2  2_da_testing.Rmd                 edgeR GLM on nhood counts with
#                                             design ~ position + contamination_status
#                                             and sensitivity analysis for
#                                             positions without a clean counterpart
#                                             (step15_retention.R).
#   Stage 3  3_exception_investigation.Rmd    INTERACTIVE. Reviewer inspects
#                                             marker-expression evidence per
#                                             flagged neighborhood cluster and
#                                             writes EXCLUDE/RETAIN decisions
#                                             into contamination_decisions.yaml.
#   Stage 4  4_decision_application.Rmd       INTERACTIVE. Applies the human
#                                             decisions to produce the final
#                                             contamination_specific_cells.csv.
#
# Stages 3-4 cannot be re-executed without human review. Stages 1-2 could in
# principle be re-executed, but the canonical migration strategy freezes the
# full step 15 output as a handoff artifact: downstream steps (16, 17, 18)
# read `contamination_specific_cells.csv` directly and do not trigger a
# re-run. This wrapper therefore ONLY validates the frozen artifact's
# existence, schema, and row count. It never re-executes the Milo graph or
# the DA testing.
#
# To unfreeze (if a reviewer needs to regenerate flags against a new L2
# integration), use the source modules under source/ and the original
# Rmds in this directory. They are preserved for reproducibility.
#
# Exempt from the V1 dry-run matrix until a re-execution decision is made.

suppressPackageStartupMessages({
    library(optparse)
})

# ============================================================
# CLI + canonical path resolution
# ============================================================

option_list <- list(
    make_option("--project-root", type = "character", default = NULL),
    make_option("--step-15-dir", type = "character", default = NULL),
    make_option("--contamination-cells", type = "character", default = NULL,
                help = "Path to contamination_specific_cells.csv (defaults to CFG_PREPROCESSING_CONTAMINATION_CELLS)"),
    make_option("--epidermal-samples", type = "character", default = NULL,
                help = "Path to step 14 epidermal_samples.yaml (defaults to CFG_PREPROCESSING_EPIDERMAL_SAMPLES)"),
    make_option("--dry-run", action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(option_list = option_list))

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

get_env <- function(key, default = NULL) {
    val <- Sys.getenv(key, unset = NA_character_)
    if (is.na(val) || val == "") return(default)
    val
}

PROJECT_ROOT        <- opts$`project-root` %||% get_env("CFG_PROJECT_ROOT")
STEP_15_DIR         <- opts$`step-15-dir` %||% get_env("CFG_PREPROCESSING_STEP_15")
CONTAMINATION_CELLS <- opts$`contamination-cells` %||% get_env("CFG_PREPROCESSING_CONTAMINATION_CELLS")
EPIDERMAL_SAMPLES   <- opts$`epidermal-samples` %||% get_env("CFG_PREPROCESSING_EPIDERMAL_SAMPLES")

if (is.null(PROJECT_ROOT) || !nzchar(PROJECT_ROOT)) {
    stop("FATAL: CFG_PROJECT_ROOT not set. Invoke via run/run_step_15_milo_contamination.sh.",
         call. = FALSE)
}
if (is.null(STEP_15_DIR) || !nzchar(STEP_15_DIR)) {
    stop("FATAL: CFG_PREPROCESSING_STEP_15 not set.", call. = FALSE)
}
if (is.null(CONTAMINATION_CELLS) || !nzchar(CONTAMINATION_CELLS)) {
    stop("FATAL: CFG_PREPROCESSING_CONTAMINATION_CELLS not set.", call. = FALSE)
}

cat("============================================================\n")
cat("Step 15 — Milo Contamination (frozen-handoff validator)\n")
cat("============================================================\n")
cat("PROJECT_ROOT:       ", PROJECT_ROOT, "\n")
cat("STEP_15_DIR:        ", STEP_15_DIR, "\n")
cat("CONTAMINATION_CELLS:", CONTAMINATION_CELLS, "\n")
cat("EPIDERMAL_SAMPLES:  ", EPIDERMAL_SAMPLES, "\n")
cat("Dry-run:            ", opts$`dry-run`, "\n")
cat("============================================================\n\n")

# ============================================================
# Validation of frozen handoff artifact
# ============================================================

cat("--- Frozen artifact validation ---\n")

if (!file.exists(CONTAMINATION_CELLS)) {
    stop(paste0(
        "FATAL: contamination_specific_cells.csv not found.\n",
        "       Expected: ", CONTAMINATION_CELLS, "\n",
        "       Step 15 outputs are consumed as a frozen handoff in the\n",
        "       canonical migration. Stages 3-4 of the original workflow\n",
        "       are interactive and cannot be re-executed from SLURM.\n",
        "       If regenerating, use scripts/step_15_milo_contamination/\n",
        "       source/ modules and the original Rmd stages, then promote\n",
        "       the resulting CSV to this path."),
         call. = FALSE)
}

info <- file.info(CONTAMINATION_CELLS)
cat(sprintf("  OK  contamination_specific_cells.csv  %.1f MB\n", info$size / 1024^2))

# Schema check — required columns for downstream consumers
required_cols <- c("cell_id", "is_flagged", "flag_reason", "weighted_flag_score")
hdr <- readLines(CONTAMINATION_CELLS, n = 1)
hdr_cols <- strsplit(hdr, ",", fixed = TRUE)[[1]]
hdr_cols <- gsub('^"|"$', "", hdr_cols)
missing <- setdiff(required_cols, hdr_cols)
if (length(missing) > 0) {
    stop(paste0("FATAL: frozen contamination CSV missing required columns: ",
                paste(missing, collapse = ", "), "\n",
                "       Found columns: ", paste(hdr_cols, collapse = ", ")),
         call. = FALSE)
}
cat(sprintf("  OK  schema             required cols present (%d cols total)\n",
            length(hdr_cols)))

# Row count (lightweight — count lines minus header)
n_rows <- length(readLines(CONTAMINATION_CELLS)) - 1L
cat(sprintf("  OK  cell count         %d\n", n_rows))

if (n_rows < 200000L) {
    warning(sprintf(
        "contamination_specific_cells.csv has %d rows; spec expects ~276K. Verify upstream L2 cell set.",
        n_rows))
}

# Optional: epidermal_samples.yaml (step 14 handoff) — not strictly needed
# by the downstream path but useful to flag if absent
if (!is.null(EPIDERMAL_SAMPLES) && nzchar(EPIDERMAL_SAMPLES)) {
    if (file.exists(EPIDERMAL_SAMPLES)) {
        cat(sprintf("  OK  epidermal_samples   %s\n", EPIDERMAL_SAMPLES))
    } else {
        cat(sprintf("  WARN epidermal_samples  %s (missing — not consumed downstream)\n",
                    EPIDERMAL_SAMPLES))
    }
}

cat("\n--- Validation passed ---\n\n")

# ============================================================
# Dry-run / frozen exit
# ============================================================

if (isTRUE(opts$`dry-run`)) {
    cat("=== DRY RUN MODE ===\n")
    cat("Step: step_15_milo_contamination\n")
    cat("Mode: frozen-handoff validator (no compute)\n")
    cat("VALIDATION PASSED\n")
    quit(status = 0)
}

cat("============================================================\n")
cat("Step 15 — FROZEN HANDOFF VERIFIED\n")
cat("============================================================\n")
cat("No compute executed. Downstream steps 16/17/18 will read:\n")
cat("  ", CONTAMINATION_CELLS, "\n")
cat("============================================================\n")
quit(status = 0)
