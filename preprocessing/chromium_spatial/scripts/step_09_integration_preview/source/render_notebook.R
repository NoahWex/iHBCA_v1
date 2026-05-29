#!/usr/bin/env Rscript

# ============================================================================
# Step 09: Integration Preview - R Markdown Renderer
# ============================================================================
# Purpose: Render integration_preview.Rmd with PROJECT_ROOT parameter
# Pattern: PROVEN from hbca_analysis render_notebook.R
# Called by: wrapper.R (SLURM worker)
# ============================================================================

suppressPackageStartupMessages({
    library(rmarkdown)
})

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: Rscript render_notebook.R <notebook_path> <output_html> <project_root>")
}

notebook_path <- args[1]
output_html <- args[2]
project_root <- args[3]

cat("============================================================================\n")
cat("Step 09: Rendering Integration Preview\n")
cat("============================================================================\n")
cat("Notebook:      ", notebook_path, "\n")
cat("Output HTML:   ", output_html, "\n")
cat("Project Root:  ", project_root, "\n")
cat("============================================================================\n\n")

# Validate notebook exists
if (!file.exists(notebook_path)) {
    stop("ERROR: Notebook not found: ", notebook_path)
}

# Set working directory to project root
setwd(project_root)
cat("Working directory set to:", getwd(), "\n\n")

# Render notebook with PROJECT_ROOT parameter
cat("Rendering notebook...\n\n")

tryCatch({
    rmarkdown::render(
        input = notebook_path,
        output_file = basename(output_html),
        output_dir = dirname(output_html),
        params = list(
            PROJECT_ROOT = project_root
        ),
        envir = new.env()
    )
    cat("\n============================================================================\n")
    cat("Notebook rendering completed successfully\n")
    cat("============================================================================\n")
    cat("Output HTML:", output_html, "\n")
    cat("============================================================================\n")
}, error = function(e) {
    cat("\n============================================================================\n")
    cat("ERROR: Notebook rendering failed\n")
    cat("============================================================================\n")
    cat(as.character(e), "\n")
    cat("============================================================================\n")
    quit(status = 1)
})
