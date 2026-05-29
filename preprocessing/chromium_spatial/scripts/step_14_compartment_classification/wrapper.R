#!/usr/bin/env Rscript
#
# Step 14 — Compartment Classification (canonical port)
#
# Assigns each cell to a tissue compartment (Epithelial, Stromal, Immune,
# Endothelial, DermEpithelial) from L2 integrated embeddings and identifies
# samples carrying dermal/epidermal contamination. Compartment calls come
# from a two-layer rule:
#
#   (1) Marker-based scoring. Each compartment is defined by a canonical
#       marker panel (module_configs.yaml → step_14_compartment_classification.
#       algorithm_params.compartment_markers): Epithelial {EPCAM, KRT8, KRT18,
#       KRT19}, DermEpithelial {KRT5, KRT14, KRT15}, Stromal {COL1A1, COL1A2,
#       DCN, LUM}, Immune {PTPRC, CD3D, CD79A, LYZ}, Endothelial {PECAM1, VWF,
#       CDH5}. The DermEpithelial panel is KRT5/14/15 on purpose — these
#       identify keratinocyte/basal squamous cells that leak in from skin
#       during breast tissue dissection.
#
#   (2) Leiden majority vote at resolution 1.0. For each cluster, the
#       dominant marker-score compartment assignment becomes the cluster
#       label, which smooths single-cell noise. A sample is declared
#       contaminated (has_epidermal = TRUE) when its DermEpithelial
#       fraction exceeds the 5% epidermal_threshold_pct configured in
#       module_configs.yaml.
#
# Step 04 HBCA label transfer outputs (SingleR-style Kumar 2023 labels) are
# read as a VALIDATION OVERLAY only — they appear in the report alongside
# the marker-based calls but do not drive compartment assignment. Missing
# step 04 outputs produce a warning, not an error, so the pipeline can run
# post-canonical-L1-supersession without needing to regenerate Kumar labels.
#
# Downstream: cell_compartments.csv gates per-cell compartment membership
# for L4 integration; epidermal_samples.yaml defines the contaminated
# sample set that step 15 Milo DA uses as its positive class.
#
# No cells are removed here — this step is annotation-only.

suppressPackageStartupMessages({
    library(yaml)
    library(rmarkdown)
    library(optparse)
})

# ============================================================
# 1. CLI and canonical path resolution
# ============================================================

option_list <- list(
    make_option("--project-root", type = "character", default = NULL,
                help = "FLEX preprocessing project root (defaults to CFG_PROJECT_ROOT)"),
    make_option("--step-14-dir", type = "character", default = NULL,
                help = "Step 14 output root (defaults to CFG_PREPROCESSING_STEP_14)"),
    make_option("--step-04-dir", type = "character", default = NULL,
                help = "Step 04 label transfer root (defaults to CFG_PREPROCESSING_STEP_04, optional)"),
    make_option("--step-12-dir", type = "character", default = NULL,
                help = "Step 12 filtered integration root (defaults to CFG_PREPROCESSING_STEP_12)"),
    make_option("--step-08-dir", type = "character", default = NULL,
                help = "Step 08 integration fallback root (defaults to CFG_PREPROCESSING_STEP_08)"),
    make_option("--module-config", type = "character", default = NULL,
                help = "module_configs.yaml path (defaults to config/module_configs.yaml)"),
    make_option("--notebook", type = "character", default = NULL,
                help = "Rmd notebook path (defaults to sibling compartment_classification.Rmd)"),
    make_option("--dry-run", action = "store_true", default = FALSE,
                help = "Validate inputs and exit without rendering the notebook")
)
opts <- parse_args(OptionParser(option_list = option_list))

# Null-coalesce helper (NA/NULL/empty-string → fallback)
`%||%` <- function(a, b) if (!is.null(a) && !is.na(a) && nzchar(a)) a else b

get_env_required <- function(key) {
    val <- Sys.getenv(key, unset = NA_character_)
    if (is.na(val) || val == "") {
        stop(paste0("FATAL: required environment variable not set: ", key, "\n",
                    "       This wrapper must be invoked via ",
                    "run/run_step_14_compartment_classification.sh which sources run/_common.sh."),
             call. = FALSE)
    }
    val
}

PROJECT_ROOT      <- opts$`project-root` %||% get_env_required("CFG_PROJECT_ROOT")
STEP_14_DIR       <- opts$`step-14-dir` %||% get_env_required("CFG_PREPROCESSING_STEP_14")
STEP_12_DIR       <- opts$`step-12-dir` %||% Sys.getenv("CFG_PREPROCESSING_STEP_12", "")
STEP_08_DIR       <- opts$`step-08-dir` %||% Sys.getenv("CFG_PREPROCESSING_STEP_08", "")
STEP_04_DIR       <- opts$`step-04-dir` %||% Sys.getenv("CFG_PREPROCESSING_STEP_04", "")
MODULE_CONFIG     <- opts$`module-config` %||%
    file.path(PROJECT_ROOT, "config", "module_configs.yaml")

SCRIPT_DIR <- local({
    cmd_args <- commandArgs(trailingOnly = FALSE)
    script_path <- sub("--file=", "", cmd_args[grep("--file=", cmd_args)])
    if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path, mustWork = FALSE))
})
NOTEBOOK_PATH <- opts$notebook %||% file.path(SCRIPT_DIR, "compartment_classification.Rmd")

REPORTS_DIR  <- file.path(STEP_14_DIR, "reports")
METADATA_DIR <- file.path(STEP_14_DIR, "metadata")
VIZ_DIR      <- file.path(STEP_14_DIR, "visualizations")
LOG_DIR      <- file.path(STEP_14_DIR, "logs")

# ============================================================
# 2. Validation preamble (fail-loud)
# ============================================================

cat("============================================================\n")
cat("Step 14 — Compartment Classification\n")
cat("============================================================\n")
cat("PROJECT_ROOT:   ", PROJECT_ROOT, "\n")
cat("STEP_14_DIR:    ", STEP_14_DIR, "\n")
cat("STEP_12_DIR:    ", STEP_12_DIR, "\n")
cat("STEP_08_DIR:    ", STEP_08_DIR, "\n")
cat("STEP_04_DIR:    ", STEP_04_DIR, "\n")
cat("Notebook:       ", NOTEBOOK_PATH, "\n")
cat("Module config:  ", MODULE_CONFIG, "\n")
cat("Dry-run:        ", opts$`dry-run`, "\n")
cat("============================================================\n\n")

assert_path <- function(label, path, must_be_dir = FALSE) {
    if (is.null(path) || is.na(path) || !nzchar(path)) {
        stop(sprintf("FATAL: %s not set\n", label), call. = FALSE)
    }
    if (!file.exists(path)) {
        stop(sprintf("FATAL: %s does not exist: %s\n", label, path), call. = FALSE)
    }
    if (must_be_dir && !dir.exists(path)) {
        stop(sprintf("FATAL: %s exists but is not a directory: %s\n", label, path), call. = FALSE)
    }
    info <- file.info(path)
    cat(sprintf("  OK  %-20s %s (%s)\n", label,
                path,
                if (info$isdir) "dir" else sprintf("%.1f MB", info$size / 1024^2)))
}

cat("--- Input validation ---\n")
assert_path("PROJECT_ROOT",  PROJECT_ROOT, must_be_dir = TRUE)
assert_path("Notebook",      NOTEBOOK_PATH)
assert_path("module_configs", MODULE_CONFIG)

# Primary L2 input: step 12 filtered integration preferred, step 08 fallback.
# We check metadata exports; the full Seurat object is reassembled by the Rmd.
primary_meta <- if (nzchar(STEP_12_DIR))
    file.path(STEP_12_DIR, "metadata", "integrated_metadata.csv") else ""
fallback_meta <- if (nzchar(STEP_08_DIR))
    file.path(STEP_08_DIR, "metadata", "integrated_metadata.csv") else ""

have_primary <- nzchar(primary_meta) && file.exists(primary_meta)
have_fallback <- nzchar(fallback_meta) && file.exists(fallback_meta)

if (have_primary) {
    cat(sprintf("  OK  L2 input          %s (step 12)\n", primary_meta))
} else if (have_fallback) {
    cat(sprintf("  OK  L2 input          %s (step 08 fallback)\n", fallback_meta))
} else {
    stop(paste0(
        "FATAL: no L2 integrated metadata found.\n",
        "       Expected one of:\n",
        "         primary:  ", primary_meta, "\n",
        "         fallback: ", fallback_meta, "\n",
        "       Run step 12 (or step 08) before step 14."),
         call. = FALSE)
}

# Step 04 HBCA labels — SOFT dependency (validation overlay only)
if (nzchar(STEP_04_DIR) && dir.exists(file.path(STEP_04_DIR, "cell_metadata"))) {
    n_label_files <- length(list.files(file.path(STEP_04_DIR, "cell_metadata"),
                                       pattern = "_labels\\.csv$"))
    cat(sprintf("  OK  Step 04 labels    %s (%d per-sample files)\n",
                file.path(STEP_04_DIR, "cell_metadata"), n_label_files))
} else {
    cat("  WARN Step 04 labels   not found; compartment report will skip Kumar 2023 overlay\n")
}

# R_LIBS_USER check (required inside container)
r_libs <- Sys.getenv("R_LIBS_USER", "")
if (!nzchar(r_libs)) {
    warning("R_LIBS_USER not set — proceeding but miloR/presto may be missing")
} else {
    cat(sprintf("  OK  R_LIBS_USER       %s\n", r_libs))
}

# Output directories
for (d in c(REPORTS_DIR, METADATA_DIR, VIZ_DIR, LOG_DIR)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(d)) stop(sprintf("FATAL: could not create output dir: %s", d), call. = FALSE)
}
cat(sprintf("  OK  Outputs writable  %s\n", STEP_14_DIR))

cat("\n--- Validation passed ---\n\n")

# ============================================================
# 3. Dry-run exit
# ============================================================

if (isTRUE(opts$`dry-run`)) {
    cat("=== DRY RUN MODE ===\n")
    cat("Step: step_14_compartment_classification\n")
    cat("Container: r_spatial\n")
    cat("Outputs planned:\n")
    cat("  ", file.path(METADATA_DIR, "cell_compartments.csv"), "\n")
    cat("  ", file.path(METADATA_DIR, "compartment_markers.csv"), "\n")
    cat("  ", file.path(METADATA_DIR, "sample_compartment_stats.csv"), "\n")
    cat("  ", file.path(METADATA_DIR, "epidermal_samples.csv"), "\n")
    cat("  ", file.path(METADATA_DIR, "epidermal_samples.yaml"), "\n")
    cat("Resources: 6 CPUs, 64 GB, 2h\n")
    cat("VALIDATION PASSED\n")
    quit(status = 0)
}

# ============================================================
# 4. Render notebook
# ============================================================

# The Rmd resolves downstream paths from PROJECT_ROOT. Export it so the
# knit environment sees the canonical preprocessing root (not the dev tree).
Sys.setenv(PROJECT_ROOT = PROJECT_ROOT)
Sys.setenv(CFG_PREPROCESSING_STEP_14 = STEP_14_DIR)
setwd(PROJECT_ROOT)

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_HTML <- file.path(REPORTS_DIR, paste0("compartment_classification_", TIMESTAMP, ".html"))
LOG_FILE <- file.path(LOG_DIR, paste0("wrapper_", TIMESTAMP, ".log"))

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
    stop(sprintf("FATAL: notebook rendering failed. See log: %s", LOG_FILE), call. = FALSE)
}

# ============================================================
# 5. Output manifest
# ============================================================

expected <- c(
    file.path(METADATA_DIR, "cell_compartments.csv"),
    file.path(METADATA_DIR, "compartment_markers.csv"),
    file.path(METADATA_DIR, "sample_compartment_stats.csv"),
    file.path(METADATA_DIR, "epidermal_samples.csv"),
    file.path(METADATA_DIR, "epidermal_samples.yaml")
)

cat("\n--- Output validation ---\n")
missing <- character(0)
for (out in expected) {
    if (file.exists(out)) {
        info <- file.info(out)
        cat(sprintf("  OK  %-40s %.1f MB\n", basename(out), info$size / 1024^2))
    } else {
        missing <- c(missing, basename(out))
        cat(sprintf("  XX  %-40s MISSING\n", basename(out)))
    }
}

if (length(missing) > 0) {
    stop(sprintf("FATAL: expected outputs missing: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)
}

cat("\n============================================================\n")
cat("Step 14 — COMPLETE\n")
cat("============================================================\n")
cat("Report:     ", OUTPUT_HTML, "\n")
cat("Log:        ", LOG_FILE, "\n")
cat("Metadata:   ", METADATA_DIR, "\n")
cat("============================================================\n")
quit(status = 0)
