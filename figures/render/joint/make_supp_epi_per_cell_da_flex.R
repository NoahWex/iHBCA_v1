#!/usr/bin/env Rscript
# Compute per-cell DA logFC for the C5 epi_content contrast on a FLEX Milo.
#
# Each cell belongs to one or more nhoods (via miloR's nhoods() membership
# matrix). For each cell, we take the average logFC across the nhoods it
# belongs to, weighted equally. Also report counts of significant nhoods
# (positive / negative) the cell sits in.
#
# Output (per_cell_logfc.csv): cell_id, mean_logfc, n_nhoods, n_sig_pos,
#                              n_sig_neg, status ("sig+", "sig-", "ns")
#
# This is the substrate for `render_c5_umap.R`. Each panel answers one question:
# "where in the FLEX cell-space does the DA signal land?"

suppressPackageStartupMessages({
  library(argparse)
  library(Matrix)
  library(miloR)
  library(data.table)
})

parser <- ArgumentParser()
parser$add_argument("--milo-rds",     required = TRUE,
                    help = "Path to Milo .rds (the scope's Milo)")
parser$add_argument("--da-csv",       required = TRUE,
                    help = "da_results.csv for this scope+contrast")
parser$add_argument("--out-csv",      required = TRUE)
parser$add_argument("--sig-fdr",      type = "double", default = 0.05)
parser$add_argument("--min-nhoods-for-sig", type = "integer", default = 1,
                    help = "Cell counts as sig+ / sig- only if at least this many of its nhoods are sig in that direction")
args <- parser$parse_args()

cat("[per-cell-da] loading milo: ", args$milo_rds, "\n", sep = "")
milo <- readRDS(args$milo_rds)
N <- nhoods(milo)  # cells × nhoods (sparse)
cat("  membership: ", nrow(N), " cells × ", ncol(N), " nhoods\n", sep = "")

cat("[per-cell-da] loading DA: ", args$da_csv, "\n", sep = "")
da <- fread(args$da_csv)
stopifnot(nrow(da) == ncol(N))

logfc       <- da$logFC
spatial_fdr <- da$SpatialFDR
sig         <- !is.na(spatial_fdr) & spatial_fdr < args$sig_fdr
sig_pos     <- sig & logfc > 0
sig_neg     <- sig & logfc < 0

# Per-cell sums via sparse mat-vec
mem_count   <- rowSums(N)                                   # n nhoods per cell
sum_logfc   <- as.numeric(N %*% logfc)
mean_logfc  <- ifelse(mem_count > 0, sum_logfc / mem_count, NA_real_)
n_sig_pos   <- as.integer(N %*% as.numeric(sig_pos))
n_sig_neg   <- as.integer(N %*% as.numeric(sig_neg))

status <- ifelse(n_sig_pos >= args$min_nhoods_for_sig & n_sig_pos > n_sig_neg, "sig+",
          ifelse(n_sig_neg >= args$min_nhoods_for_sig & n_sig_neg > n_sig_pos, "sig-",
          "ns"))

cell_ids <- rownames(N)
if (is.null(cell_ids)) cell_ids <- colnames(milo)
stopifnot("cannot extract cell_ids from Milo" = !is.null(cell_ids))

out <- data.table(
  cell_id    = cell_ids,
  mean_logfc = round(mean_logfc, 5),
  n_nhoods   = mem_count,
  n_sig_pos  = n_sig_pos,
  n_sig_neg  = n_sig_neg,
  status     = status
)

fwrite(out, args$out_csv)
cat("[per-cell-da] wrote ", args$out_csv, " (", nrow(out), " cells)\n", sep = "")
cat("  status: sig+ ", sum(status == "sig+"),
    "  sig- ", sum(status == "sig-"),
    "  ns ", sum(status == "ns"), "\n", sep = "")
cat("  mean_logfc range: [", round(min(mean_logfc, na.rm = TRUE), 3), ", ",
    round(max(mean_logfc, na.rm = TRUE), 3), "]\n", sep = "")
