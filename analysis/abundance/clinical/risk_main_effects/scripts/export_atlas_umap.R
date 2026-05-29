#!/usr/bin/env Rscript
# export_atlas_umap.R - one-shot: write atlas_umap.tsv (cell_id, UMAP_1, UMAP_2)
# alongside cells.tsv / nhoods.mtx in atlas_nhoods/. Idempotent.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(readr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "export_atlas_umap.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
cfg <- parse_args_inquiry()
paths <- cfg$paths

force_flag <- "--force" %in% commandArgs(trailingOnly = TRUE)
out_dir <- file.path(paths$inquiry_root, "outputs", "atlas_nhoods")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_path <- file.path(out_dir, "atlas_umap.tsv")
if (!force_flag && file.exists(out_path)) {
  cat("Already exported; pass --force to redo. Exiting.\n"); quit(status = 0)
}

cat("Loading milo (~3-5 min)...\n")
source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)
umap_key <- intersect(c("UMAP_scVI", "UMAP", "X_umap", "umap"),
                       reducedDimNames(milo))[1]
if (is.na(umap_key)) stop("No UMAP reducedDim. Found: ",
                          paste(reducedDimNames(milo), collapse = ", "))
umap_df <- as.data.frame(reducedDim(milo, umap_key))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$cell_id <- colnames(milo)
umap_df <- umap_df[, c("cell_id", "UMAP_1", "UMAP_2")]
cat(sprintf("Writing %s (%d cells)\n", out_path, nrow(umap_df)))
write_tsv(umap_df, out_path)
cat("Done.\n")
