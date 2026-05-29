#!/usr/bin/env Rscript
# 09d_concat_l1.R — Concatenate per-compartment L1 harmonized labels.
#
# Reads harmonized_l1_labels_{Epithelial,Stromal,Immune}.csv and writes a
# single harmonized_l1_labels.csv at the top-level L1 outputs dir.
#
# Validation gates:
#   - exactly 270,035 rows in the concatenated output
#   - no duplicate cell_ids
#   - every row has a non-NA L1.0 value

suppressPackageStartupMessages({
  library(argparse)
})

parser <- ArgumentParser(description = "Concatenate per-compartment L1 harmonized labels")
parser$add_argument("--epithelial", required = TRUE,
                    help = "harmonized_l1_labels_Epithelial.csv")
parser$add_argument("--stromal", required = TRUE,
                    help = "harmonized_l1_labels_Stromal.csv")
parser$add_argument("--immune", required = TRUE,
                    help = "harmonized_l1_labels_Immune.csv")
parser$add_argument("--out", required = TRUE,
                    help = "Output concatenated CSV")
parser$add_argument("--expected-rows", type = "integer", default = 270035L,
                    help = "Expected total row count (default: 270035)")
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

cat("=== 09d_concat_l1 ===\n")
cat("Epithelial: ", args$epithelial, "\n")
cat("Stromal:    ", args$stromal, "\n")
cat("Immune:     ", args$immune, "\n")
cat("Out:        ", args$out, "\n")
cat("Expected:   ", args$expected_rows, " rows\n")
cat("Dry run:    ", args$dry_run, "\n\n")

inputs <- c(
  Epithelial = args$epithelial,
  Stromal    = args$stromal,
  Immune     = args$immune
)

missing_items <- character(0)
for (nm in names(inputs)) {
  f <- inputs[[nm]]
  ok <- file.exists(f)
  cat(sprintf("  [%s] %-12s %s\n", if (ok) "OK " else "MISS", nm, f))
  if (!ok) missing_items <- c(missing_items, nm)
}
if (length(missing_items) > 0) {
  stop("Missing required inputs: ", paste(missing_items, collapse = ", "))
}

dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)

if (isTRUE(args$dry_run)) {
  cat("\nVALIDATION PASSED (dry run) — exiting before concat\n")
  quit(status = 0)
}

# Read + concat ---------------------------------------------------------------
dfs <- list()
for (nm in names(inputs)) {
  df <- read.csv(inputs[[nm]], stringsAsFactors = FALSE)
  df$compartment <- nm
  cat(sprintf("  %-12s %d rows, %d cols\n", nm, nrow(df), ncol(df)))
  dfs[[nm]] <- df
}

# Verify column consistency
col_sets <- lapply(dfs, colnames)
if (length(unique(col_sets)) != 1) {
  # Keep the intersection + compartment
  common <- Reduce(intersect, col_sets)
  cat("WARNING: column mismatch across compartments. Using intersection:",
      paste(common, collapse = ", "), "\n")
  dfs <- lapply(dfs, function(df) df[, common, drop = FALSE])
}

concat <- do.call(rbind, dfs)
rownames(concat) <- NULL
cat(sprintf("\nConcatenated: %d rows, %d cols\n", nrow(concat), ncol(concat)))

# Validation gates ------------------------------------------------------------
cat("\n=== Validation gates ===\n")
n_rows <- nrow(concat)
cat(sprintf("  Row count: %d (expected %d)\n", n_rows, args$expected_rows))
if (n_rows != args$expected_rows) {
  stop(sprintf("FAIL: row count %d != expected %d",
               n_rows, args$expected_rows))
}
cat("  Row count: OK\n")

n_dup <- sum(duplicated(concat$cell_id))
cat(sprintf("  Duplicate cell_ids: %d\n", n_dup))
if (n_dup > 0) {
  stop(sprintf("FAIL: %d duplicate cell_ids in concatenated output", n_dup))
}
cat("  Duplicate check: OK\n")

n_na_l1 <- sum(is.na(concat$L1.0) | concat$L1.0 == "")
cat(sprintf("  NA L1.0: %d\n", n_na_l1))
if (n_na_l1 > 0) {
  cat("  WARNING: ", n_na_l1, "cells have missing L1.0\n")
}

# Save ------------------------------------------------------------------------
write.csv(concat, args$out, row.names = FALSE)
cat("\nSaved:", args$out, "\n")
cat("Unique L1.0 types:", length(unique(concat$L1.0)), "\n")
cat("\nL1.0 distribution:\n")
print(sort(table(concat$L1.0), decreasing = TRUE))
cat("\nL1.0_source breakdown:\n")
print(table(concat$L1.0_source))
cat("\nDone.\n")
