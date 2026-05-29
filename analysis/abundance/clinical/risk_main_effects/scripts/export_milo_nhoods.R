suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(Matrix); library(readr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "export_milo_nhoods.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

force_flag <- "--force" %in% commandArgs(trailingOnly = TRUE)

out_dir <- file.path(paths$inquiry_root, "outputs", "atlas_nhoods")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
mtx_path  <- file.path(out_dir, "nhoods.mtx")
cell_path <- file.path(out_dir, "cells.tsv")
nh_path   <- file.path(out_dir, "nhood_index.tsv")

if (!force_flag && all(file.exists(c(mtx_path, cell_path, nh_path)))) {
  cat("Already exported; pass --force to redo. Exiting.\n")
  quit(status = 0)
}

cat("Loading milo (~3-5 min)...\n")
source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)
cat(sprintf("milo: %d cells, %d nhoods\n", ncol(milo), ncol(nhoods(milo))))

nh <- nhoods(milo)
cat(sprintf("nhoods sparse: %s, nnz=%d\n",
            paste(dim(nh), collapse = "x"), Matrix::nnzero(nh)))

cat("Writing nhoods.mtx...\n")
Matrix::writeMM(nh, mtx_path)

cat("Writing cells.tsv...\n")
writeLines(colnames(milo), cell_path)

cat("Writing nhood_index.tsv...\n")
# nhood_index = the cell index in milo for each nhood's index cell (the cell
# whose neighborhood defines the nhood). For the panel pipeline we just need
# 1..N here so consumers can map column index <-> Nhood ID consistent with
# da_results.csv / nhood_groups.csv.
writeLines(as.character(seq_len(ncol(nh))), nh_path)

cat("Done.\n")
cat(sprintf("  %s (%.1f MB)\n", mtx_path,
            file.info(mtx_path)$size / 1024^2))
