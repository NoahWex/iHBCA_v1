#!/usr/bin/env Rscript
# fix_L2_marker_symbols.R
# Post-process stageF3_L2_markers/<contrast>/<L2>_vs_parity.csv files —
# the original run left `symbol` as NA because the gene mapping TSV column
# names didn't match the auto-detect logic. Re-join here from the canonical
# gene_symbol_to_ensembl.tsv.

suppressPackageStartupMessages({
  library(dplyr); library(readr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "fix_L2_marker_symbols.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

gm <- read_tsv(paths$inputs$gene_mapping, show_col_types = FALSE)
cat(sprintf("Gene mapping: %d rows, cols=%s\n", nrow(gm),
            paste(colnames(gm), collapse=",")))
# Standardize: the canonical TSV has gene_symbol + ensembl_id
sym_lookup <- setNames(gm$gene_symbol, gm$ensembl_id)
cat(sprintf("Lookup size: %d\n", length(sym_lookup)))

l2_root <- file.path(paths$inquiry_root, "outputs", "stageF3_L2_markers")
contrasts <- list.dirs(l2_root, recursive = FALSE)

n_total <- 0
for (cdir in contrasts) {
  files <- list.files(cdir, pattern = "_vs_parity\\.csv$", full.names = TRUE)
  for (f in files) {
    d <- read_csv(f, show_col_types = FALSE)
    if (!"gene_id" %in% colnames(d) || nrow(d) == 0) next
    n_na_before <- sum(is.na(d$symbol) | d$symbol == "NA")
    d$symbol <- unname(sym_lookup[d$gene_id])
    n_filled <- sum(!is.na(d$symbol))
    write_csv(d, f)
    n_total <- n_total + 1
    cat(sprintf("  %s: %d/%d symbols (was %d NA)\n",
                basename(f), n_filled, nrow(d), n_na_before))
  }
}
cat(sprintf("\nProcessed %d files\n", n_total))
