#!/usr/bin/env Rscript
# Merge per-chunk limma CSVs produced by 14b_marker_limma.R --chunk-idx into
# the standard limma_markers_{level}.csv and limma_top_markers_{level}.csv.
#
# Run after all array tasks for a given resolution complete:
#   Rscript merge_14b_limma_chunks.R \
#     --output-dir /path/to/limma/leiden_0.3 \
#     --level-tag leiden_0.3 \
#     --n-top 20

suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser(description = "Merge chunked limma outputs")
parser$add_argument("--output-dir", required = TRUE,
                    help = "Directory containing limma_markers_{level}_c*.csv files")
parser$add_argument("--level-tag", required = TRUE, help = "Resolution tag (e.g. leiden_0.3)")
parser$add_argument("--n-top", type = "integer", default = 20L,
                    help = "Top N markers per type for top-markers file")
args <- parser$parse_args()

pattern  <- paste0("limma_markers_", args$level_tag, "_c[0-9]+\\.csv$")
chunk_files <- list.files(args$output_dir, pattern = pattern, full.names = TRUE)

if (length(chunk_files) == 0) {
  cat(sprintf("No chunk files matching '%s' in %s\n", pattern, args$output_dir))
  quit(status = 1)
}

chunk_files <- chunk_files[order(as.integer(
  sub(".*_c([0-9]+)\\.csv$", "\\1", basename(chunk_files))))]
cat(sprintf("Merging %d chunk files:\n", length(chunk_files)))
for (f in chunk_files) cat(sprintf("  %s\n", basename(f)))

all_df <- do.call(rbind, lapply(chunk_files, function(f) {
  df <- read.csv(f, stringsAsFactors = FALSE)
  cat(sprintf("  Read %d rows from %s\n", nrow(df), basename(f)))
  df
}))
rownames(all_df) <- NULL
cat(sprintf("Total rows: %d across %d cell types\n",
            nrow(all_df), length(unique(all_df$cell_type))))

full_path <- file.path(args$output_dir,
  paste0("limma_markers_", args$level_tag, ".csv"))
write.csv(all_df, full_path, row.names = FALSE)
cat(sprintf("Written: %s\n", full_path))

top_df <- do.call(rbind, lapply(split(all_df, all_df$cell_type), function(df) {
  df <- df[!is.na(df$padj) & df$padj < 0.05, ]
  df <- df[order(df$padj, -abs(df$log2FoldChange)), ]
  head(df, args$n_top)
}))
rownames(top_df) <- NULL

top_path <- file.path(args$output_dir,
  paste0("limma_top_markers_", args$level_tag, ".csv"))
write.csv(top_df, top_path, row.names = FALSE)
cat(sprintf("Top %d markers: %s (%d rows)\n", args$n_top, top_path, nrow(top_df)))

cat("\n=== Top 5 markers per cell type ===\n")
for (ct in sort(unique(top_df$cell_type))) {
  sub   <- top_df[top_df$cell_type == ct, ]
  genes <- head(sub$gene, 5)
  lfcs  <- head(round(sub$log2FoldChange, 1), 5)
  cat(sprintf("  %s: %s\n", ct,
              paste(sprintf("%s(lfc=%s)", genes, lfcs), collapse = ", ")))
}

cat("\nDone.\n")
