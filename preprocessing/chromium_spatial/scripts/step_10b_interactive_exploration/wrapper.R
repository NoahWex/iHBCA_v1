#!/usr/bin/env Rscript

# ============================================================================
# Step 10b: Interactive Artifact Exploration - Wrapper
# ============================================================================
# SKIP_DRY_RUN: exempt from V1 dry-run matrix — see migration notes below.
#
# Purpose
#   Validate that the frozen step 10b manifest update file exists and is
#   well-formed. Step 10b is the manual artifact-discovery stage of the
#   canonical pipeline: an R Markdown report is rendered interactively in
#   RStudio, a reviewer inspects global/local manifold outliers, and the
#   result is committed as step10b_manifest_update.csv (the four-way
#   filter). Per the revised migration model, this file is FROZEN — the
#   canonical copy lives at $CFG_PREPROCESSING_CELL_RETENTION_LIST and is
#   treated as an input by all downstream wrappers (11, 12, 13). This
#   wrapper therefore performs no batch compute; it only validates the
#   frozen file and emits a dry-run-compatible plan.
#
# Why no batch execution
#   - The original 01_Preprocessing step 10b runs an R Markdown notebook
#     that requires reviewer judgement (cluster inspection, GMM threshold
#     tuning, trash cluster labelling). That is not reproducible in batch
#     mode without pre-committed thresholds.
#   - The current frozen output (step10b_manifest_update.csv, ~22 MB at
#     $CFG_PREPROCESSING_STEP_10B/exports/) is the agreed canonical L2 cell
#     gate. Re-running step 10b would invalidate the downstream step 11/12
#     outputs, which per I3 audit ARE populated and consumed by later steps.
#   - If a reviewer wants to redo step 10b manually, they can render
#     source/artifact_discovery.Rmd in RStudio against the step 10a
#     outputs; this wrapper does not orchestrate that.
#
# Inputs (validated only)
#   - $CFG_PREPROCESSING_CELL_RETENTION_LIST
#       (= ${CFG_PREPROCESSING_STEP_10B}/exports/step10b_manifest_update.csv)
#     Schema: cell_id, step_10b_manifold_pass [, step_10b_removal_reason]
#
# Outputs
#   - (none written by this wrapper; canonical file is preserved in place)
#
# Note on manifest-update merge:
#   This wrapper does not merge step_10b_manifold_pass into
#   central_cell_status.csv. That merge is handled by step 12's wrapper at
#   read time (it supports both "column in central" and "merge from
#   retention list" paths).
# ============================================================================

suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Step 10b: Interactive Exploration (frozen validator)")
parser$add_argument("--project-root", type = "character",
                    default = Sys.getenv("CFG_PROJECT_ROOT", ""))
parser$add_argument("--cell-retention-list", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_CELL_RETENTION_LIST", ""))
parser$add_argument("--output-root", type = "character",
                    default = Sys.getenv("CFG_PREPROCESSING_STEP_10B", ""))
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

log_event <- function(a, v = NULL, l = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  if (is.null(v)) cat(sprintf("[%s] %s step_10b %s\n", ts, l, a))
  else cat(sprintf("[%s] %s step_10b %s=%s\n", ts, l, a, as.character(v)))
}

fail_loud <- function(msg, hint = NULL) {
  log_event("FAIL", msg, l = "ERROR")
  if (!is.null(hint)) cat(sprintf("  hint: %s\n", hint))
  quit(save = "no", status = 2)
}

log_event("start")
log_event("mode", "frozen_validator_only")

if (identical(args$cell_retention_list, "")) {
  fail_loud("cell-retention-list empty",
            "set --cell-retention-list or CFG_PREPROCESSING_CELL_RETENTION_LIST")
}
if (!file.exists(args$cell_retention_list)) {
  fail_loud(sprintf("frozen step 10b retention list not found: %s",
                    args$cell_retention_list),
            paste0("the canonical file is ",
                   "${CFG_PREPROCESSING_STEP_10B}/exports/step10b_manifest_update.csv; ",
                   "if it is missing, restore from the 01_Preprocessing dev repo"))
}
size <- file.info(args$cell_retention_list)$size
log_event("cell_retention_list", args$cell_retention_list)
log_event("bytes", size)
if (size == 0L) {
  fail_loud("frozen step 10b retention list is empty (0 bytes)")
}

# Schema check — peek at the header
header_line <- readLines(args$cell_retention_list, n = 1L)
header_cols <- strsplit(header_line, ",", fixed = TRUE)[[1]]
header_cols <- gsub("^\"|\"$", "", header_cols)
required_cols <- c("cell_id", "step_10b_manifold_pass")
missing <- setdiff(required_cols, header_cols)
if (length(missing) > 0) {
  fail_loud(sprintf("retention list missing required columns: %s",
                    paste(missing, collapse = ", ")),
            "expected schema: cell_id, step_10b_manifold_pass [, step_10b_removal_reason]")
}
log_event("schema_ok", paste(required_cols, collapse = ","))

if (isTRUE(args$dry_run)) {
  cat("=== DRY RUN MODE ===\n")
  cat("Step: step_10b_interactive_exploration (SKIP_DRY_RUN=frozen)\n")
  cat("Container: r_spatial (no batch execution)\n")
  cat("Inputs validated:\n")
  cat(sprintf("  %s: OK (%d bytes)\n", args$cell_retention_list, size))
  cat("Outputs planned:\n")
  cat("  (none — this step is a frozen validator; the retention list is canonical)\n")
  cat("Resources: N/A (validation only)\n")
  cat("VALIDATION PASSED\n")
  quit(save = "no", status = 0)
}

# Non-dry-run path: also validate only, do not rewrite the retention list.
# This keeps the canonical file immutable across pipeline re-runs.
cat("[step_10b] Frozen mode: canonical retention list validated, no compute executed.\n")
cat(sprintf("  File: %s\n", args$cell_retention_list))
cat(sprintf("  Size: %d bytes\n", size))
cat("  See wrapper.R header for the rationale (revised migration model).\n")
log_event("done")
quit(save = "no", status = 0)
