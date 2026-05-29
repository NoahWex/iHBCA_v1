#!/usr/bin/env Rscript
#
# Step 19b — Integration Preview (canonical port, pre-sweep version)
#
# Publication-grade QC report of the final integrated atlas. Reassembles
# a Seurat object from the step 17 scVI sidecars + step 14 compartment
# labels + step 15 contamination flags + raw H5 counts, runs per-cluster
# marker discovery at each Leiden resolution, and renders an HTML preview
# with UMAP composition, retention tracking, and marker tables.
#
# STATUS: Pre-sweep version ported for provenance. Step 17 in the
# canonical pipeline is a sweep pipeline (Phase B), not the single-config
# scVI in the dev tree — so this wrapper points at the canonical sweep
# WINNER outputs (CFG_CANONICAL_SWEEP_WINNER) rather than the old
# 17_ScviIntegration directory. The Rmd notebook, which is the heaviest
# part of step 19b, will be re-ported from the dev tree AFTER the sweep
# winner exists and its sidecar schema is locked; for V1 we keep the
# wrapper in place and mark the step exempt from the dry-run matrix.
#
# Exempt from V1 dry-run matrix (per migration plan): no sweep winner
# sidecars exist yet, so the validation preamble would fail closed by
# design.

suppressPackageStartupMessages({
    library(optparse)
})

option_list <- list(
    make_option("--project-root", type = "character", default = NULL),
    make_option("--step-19b-dir", type = "character", default = NULL),
    make_option("--sweep-winner", type = "character", default = NULL),
    make_option("--step-14-dir", type = "character", default = NULL),
    make_option("--step-15-contamination", type = "character", default = NULL),
    make_option("--notebook", type = "character", default = NULL),
    make_option("--dry-run", action = "store_true", default = FALSE)
)
opts <- parse_args(OptionParser(option_list = option_list))

`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b
get_env <- function(key, default = NULL) {
    val <- Sys.getenv(key, unset = NA_character_)
    if (is.na(val) || val == "") return(default)
    val
}

PROJECT_ROOT     <- opts$`project-root` %||% get_env("CFG_PROJECT_ROOT")
STEP_19B_DIR     <- opts$`step-19b-dir` %||% get_env("CFG_PREPROCESSING_STEP_19B")
SWEEP_WINNER     <- opts$`sweep-winner` %||% get_env("CFG_CANONICAL_SWEEP_WINNER")
STEP_14_DIR      <- opts$`step-14-dir` %||% get_env("CFG_PREPROCESSING_STEP_14")
CONTAMINATION    <- opts$`step-15-contamination` %||% get_env("CFG_PREPROCESSING_CONTAMINATION_CELLS")

if (is.null(PROJECT_ROOT) || !nzchar(PROJECT_ROOT)) {
    stop("FATAL: CFG_PROJECT_ROOT not set. Invoke via run/run_step_19b_integration_preview.sh.",
         call. = FALSE)
}
if (is.null(STEP_19B_DIR) || !nzchar(STEP_19B_DIR)) {
    stop("FATAL: CFG_PREPROCESSING_STEP_19B not set.", call. = FALSE)
}

SCRIPT_DIR <- local({
    cmd_args <- commandArgs(trailingOnly = FALSE)
    script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
    if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path, mustWork = FALSE))
})
NOTEBOOK_PATH <- opts$notebook %||% file.path(SCRIPT_DIR, "integration_preview.Rmd")

cat("============================================================\n")
cat("Step 19b — Integration Preview (canonical port)\n")
cat("============================================================\n")
cat("PROJECT_ROOT:     ", PROJECT_ROOT, "\n")
cat("STEP_19B_DIR:     ", STEP_19B_DIR, "\n")
cat("SWEEP_WINNER:     ", SWEEP_WINNER, "\n")
cat("STEP_14_DIR:      ", STEP_14_DIR, "\n")
cat("CONTAMINATION:    ", CONTAMINATION, "\n")
cat("Notebook:         ", NOTEBOOK_PATH, "\n")
cat("Dry-run:          ", opts$`dry-run`, "\n")
cat("============================================================\n\n")

# ============================================================
# Validation preamble (fail-loud)
# ============================================================

cat("--- Input validation ---\n")

if (is.null(SWEEP_WINNER) || !nzchar(SWEEP_WINNER) || !dir.exists(SWEEP_WINNER)) {
    stop(paste0(
        "FATAL: sweep winner directory not found.\n",
        "       Expected: ", SWEEP_WINNER, "\n",
        "       Step 19b consumes the step 17 sweep winner sidecars:\n",
        "         latent.csv, umap.csv, leiden_multi.csv, metadata_with_umap.csv\n",
        "       EXEMPT FROM V1 DRY-RUN until Phase B produces a winner."),
         call. = FALSE)
}

# Per-compartment sweep: sidecars live under compartments/{Epithelial,Stromal,Immune}/
compartments <- c("Epithelial", "Stromal", "Immune")
required_sidecars <- c("metadata_with_umap.csv")
missing <- character(0)
for (comp in compartments) {
    comp_dir <- file.path(SWEEP_WINNER, "compartments", comp)
    for (f in required_sidecars) {
        p <- file.path(comp_dir, f)
        key <- paste0(comp, "/", f)
        if (!file.exists(p)) {
            missing <- c(missing, key)
            cat(sprintf("  XX  %s\n", key))
        } else {
            info <- file.info(p)
            cat(sprintf("  OK  %-45s %.1f MB\n", key, info$size / 1024^2))
        }
    }
}
if (length(missing) > 0) {
    stop(paste0("FATAL: sweep winner missing required sidecars: ",
                paste(missing, collapse = ", ")), call. = FALSE)
}

# Step 14 compartments (validation overlay)
cell_compartments <- get_env("CFG_PREPROCESSING_CELL_COMPARTMENTS",
                              file.path(STEP_14_DIR, "metadata/cell_compartments.csv"))
if (file.exists(cell_compartments)) {
    cat(sprintf("  OK  cell_compartments.csv  %s\n", cell_compartments))
} else {
    cat(sprintf("  WARN cell_compartments.csv  %s (missing — validation overlay disabled)\n",
                cell_compartments))
}

# Step 15 contamination flags (validation overlay)
if (!is.null(CONTAMINATION) && file.exists(CONTAMINATION)) {
    cat(sprintf("  OK  contamination_cells    %s\n", CONTAMINATION))
} else {
    cat(sprintf("  WARN contamination_cells   %s (missing)\n", CONTAMINATION))
}

# Output dirs
for (d in c(STEP_19B_DIR,
            file.path(STEP_19B_DIR, "reports"),
            file.path(STEP_19B_DIR, "metadata"),
            file.path(STEP_19B_DIR, "markers_by_resolution"),
            file.path(STEP_19B_DIR, "preview_plots"),
            file.path(STEP_19B_DIR, "seurat_objects"),
            file.path(STEP_19B_DIR, "logs"))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
cat(sprintf("  OK  Outputs writable     %s\n", STEP_19B_DIR))

# Rmd notebook — ported from the dev tree during Phase B unfreezing.
if (!file.exists(NOTEBOOK_PATH)) {
    stop(paste0(
        "FATAL: step 19b integration_preview.Rmd not yet ported.\n",
        "       Expected: ", NOTEBOOK_PATH, "\n",
        "       The Rmd is the heaviest part of step 19b and will be re-ported\n",
        "       from the dev tree once the canonical sweep winner schema is\n",
        "       locked. See scripts/step_19b_integration_preview/README.md."),
         call. = FALSE)
}
cat(sprintf("  OK  notebook             %s\n", NOTEBOOK_PATH))

cat("\n--- Validation passed ---\n\n")

# ============================================================
# Dry-run exit
# ============================================================

if (isTRUE(opts$`dry-run`)) {
    cat("=== DRY RUN MODE (EXEMPT FROM V1 MATRIX) ===\n")
    cat("Step: step_19b_integration_preview\n")
    cat("Note: exempt from V1 matrix until Phase B sweep winner exists.\n")
    cat("VALIDATION PASSED\n")
    quit(status = 0)
}

# ============================================================
# Render notebook
# ============================================================

Sys.setenv(PROJECT_ROOT = PROJECT_ROOT)
Sys.setenv(CFG_PREPROCESSING_STEP_19B = STEP_19B_DIR)
Sys.setenv(CFG_CANONICAL_SWEEP_WINNER = SWEEP_WINNER)
setwd(PROJECT_ROOT)

suppressPackageStartupMessages(library(rmarkdown))

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_HTML <- file.path(STEP_19B_DIR, "reports",
                         paste0("integration_preview_", TIMESTAMP, ".html"))
LOG_FILE    <- file.path(STEP_19B_DIR, "logs",
                         paste0("wrapper_", TIMESTAMP, ".log"))

cat(sprintf("--- Rendering notebook ---\n  input:  %s\n  output: %s\n\n",
            NOTEBOOK_PATH, OUTPUT_HTML))

log_conn <- file(LOG_FILE, open = "wt")
sink(log_conn, type = "output", split = TRUE)
sink(log_conn, type = "message")

render_success <- FALSE
tryCatch({
    rmarkdown::render(
        input         = NOTEBOOK_PATH,
        output_file   = basename(OUTPUT_HTML),
        output_dir    = dirname(OUTPUT_HTML),
        output_format = "html_document",
        params        = list(project_root = PROJECT_ROOT),
        envir         = new.env(),
        quiet         = FALSE
    )
    render_success <- TRUE
}, error = function(e) {
    cat("\nERROR during rendering:\n", conditionMessage(e), "\n", sep = "")
})

sink(type = "message")
sink(type = "output")
close(log_conn)

if (!render_success || !file.exists(OUTPUT_HTML)) {
    stop(paste0("FATAL: notebook rendering failed. See log: ", LOG_FILE), call. = FALSE)
}

cat("\n============================================================\n")
cat("Step 19b — COMPLETE\n")
cat("============================================================\n")
cat("Report:  ", OUTPUT_HTML, "\n")
cat("Log:     ", LOG_FILE, "\n")
cat("============================================================\n")
quit(status = 0)
