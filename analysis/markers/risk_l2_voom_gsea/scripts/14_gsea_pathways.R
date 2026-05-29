#!/usr/bin/env Rscript
# fgsea-based pre-ranked GSEA on voom DE outputs.
#
# Per-track design: one invocation processes ALL DE CSVs in --de-dir, ranks
# genes by limma t-statistic, runs fgsea against MSigDB Hallmark + Reactome,
# writes per-cell-type pathway tables + a combined long table for heatmap
# rendering.
#
# Designed to run on each parallel track:
#   A_full × L1   : --de-dir outputs/de_results/A_full/A2_full   --label l1
#   A_tested × L1 : --de-dir outputs/de_results/A_tested/A2_tested --label l1
#   A_full × L2   : --de-dir outputs/de_results/A_full/A2_sampletype --label l2
#   A_tested × L2 : --de-dir outputs/de_results/C_tested/A2_sampletype_tested --label l2
#
# Outputs (per track):
#   <out-dir>/pathways_long.csv          # cell_type × pathway × NES + padj + leading_edge
#   <out-dir>/pathways_per_cell/<L*>.csv # one CSV per cell type
#   <out-dir>/audit.csv                  # which cell types got fgsea results vs skipped
#
# Code pattern provenance: fgsea + msigdbr idiom matches Bioconductor canonical
# workflow (fgsea vignette). Pre-rank by t-statistic per limma authors'
# recommendation (Phipson 2016).

suppressPackageStartupMessages({
  library(argparse)
  library(fgsea)
  library(msigdbr)
  library(dplyr)
  library(readr)
  library(tibble)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--de-dir", required = TRUE,
                    help = "Directory of voom DE CSVs (one per cell type)")
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--collections", default = "H",
                    help = "Comma-separated MSigDB collection specs (CATEGORY or CATEGORY:SUBCATEGORY). Default H (Hallmark, 50 curated pathways).")
parser$add_argument("--species", default = "Homo sapiens")
parser$add_argument("--min-size", type = "integer", default = 15)
parser$add_argument("--max-size", type = "integer", default = 500)
parser$add_argument("--fdr",       type = "double",  default = 0.10)
parser$add_argument("--seed",      type = "integer", default = 42)
args <- parser$parse_args()

set.seed(args$seed)
dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(args$out_dir, "pathways_per_cell"),
           showWarnings = FALSE, recursive = TRUE)

# ----- load MSigDB gene sets --------------------------------------------
cat("[1] load MSigDB gene sets\n")
collections <- strsplit(args$collections, ",", fixed = TRUE)[[1]] |> trimws()
pathways <- list()
for (col in collections) {
  # Parse "H" or "C2:CP:REACTOME" → category + subcategory.
  # MSigDB subcategories use colon-separated notation (e.g., "CP:REACTOME").
  parts <- strsplit(col, ":", fixed = TRUE)[[1]]
  category <- parts[1]
  subcat   <- if (length(parts) > 1) paste(parts[-1], collapse = ":") else NULL
  cat(sprintf("  loading %s (category=%s, subcat=%s)\n",
              col, category, subcat %||% "<none>"))

  ms <- if (!is.null(subcat)) {
    msigdbr::msigdbr(species = args$species, category = category,
                     subcategory = subcat)
  } else {
    msigdbr::msigdbr(species = args$species, category = category)
  }
  pw <- split(ms$gene_symbol, ms$gs_name)
  # Prefix with collection name for traceability
  names(pw) <- paste0(col, "::", names(pw))
  pathways <- c(pathways, pw)
}
cat(sprintf("  total pathways loaded: %d (collections: %s)\n",
            length(pathways), paste(collections, collapse = ", ")))

# ----- list DE CSVs -----------------------------------------------------
de_csvs <- list.files(args$de_dir, pattern = "\\.csv$", full.names = TRUE)
de_csvs <- de_csvs[!grepl("\\.skipped\\.csv$", de_csvs)]
cat(sprintf("[2] found %d non-skipped DE CSVs in %s\n",
            length(de_csvs), args$de_dir))

# ----- per-cell-type fgsea ----------------------------------------------
audit_rows <- list()
all_pathways <- list()

for (csv in de_csvs) {
  cell_type <- sub("\\.csv$", "", basename(csv))
  cat(sprintf("\n=== %s ===\n", cell_type))

  # Detect skipped-stub vs real DE
  hdr <- readr::read_lines(csv, n_max = 1)
  if (grepl("^skipped", hdr)) {
    cat("  skipped stub — no DE; gating out\n")
    audit_rows[[length(audit_rows) + 1L]] <- list(
      cell_type = cell_type, n_genes_ranked = 0L,
      n_pathways_sig = 0L, status = "skipped_de"
    )
    next
  }

  de <- readr::read_csv(csv, show_col_types = FALSE)
  if (!all(c("t", "symbol") %in% colnames(de))) {
    cat("  missing t/symbol columns; gating out\n")
    audit_rows[[length(audit_rows) + 1L]] <- list(
      cell_type = cell_type, n_genes_ranked = 0L,
      n_pathways_sig = 0L, status = "missing_columns"
    )
    next
  }

  # Pre-rank by t-statistic. Drop NA + duplicate symbols (keep highest |t|).
  ranks <- de |>
    dplyr::filter(!is.na(t), !is.na(symbol), symbol != "") |>
    dplyr::group_by(symbol) |>
    dplyr::slice_max(abs(t), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(desc(t))

  rank_vec <- setNames(ranks$t, ranks$symbol)
  cat(sprintf("  ranked %d genes (t range [%.2f, %.2f])\n",
              length(rank_vec), min(rank_vec), max(rank_vec)))

  fgsea_res <- tryCatch(
    fgsea::fgsea(pathways = pathways, stats = rank_vec,
                 minSize = args$min_size, maxSize = args$max_size,
                 eps = 0),
    error = function(e) NULL
  )
  if (is.null(fgsea_res) || nrow(fgsea_res) == 0L) {
    cat("  fgsea returned 0 rows\n")
    audit_rows[[length(audit_rows) + 1L]] <- list(
      cell_type = cell_type, n_genes_ranked = length(rank_vec),
      n_pathways_sig = 0L, status = "fgsea_empty"
    )
    next
  }

  fgsea_res$cell_type <- cell_type
  fgsea_res$collection <- sub("::.*$", "", fgsea_res$pathway)
  fgsea_res$pathway_name <- sub("^[^:]*::", "", fgsea_res$pathway)
  fgsea_res$leadingEdge <- vapply(fgsea_res$leadingEdge,
                                  function(x) paste(x, collapse = ";"),
                                  character(1))
  fgsea_res <- as.data.frame(fgsea_res) |>
    dplyr::select(cell_type, collection, pathway_name, pathway,
                  pval, padj, ES, NES, size, leadingEdge) |>
    dplyr::arrange(padj)

  n_sig <- sum(fgsea_res$padj < args$fdr, na.rm = TRUE)
  cat(sprintf("  pathways tested: %d; sig at FDR<%.2f: %d\n",
              nrow(fgsea_res), args$fdr, n_sig))

  readr::write_csv(fgsea_res,
                   file.path(args$out_dir, "pathways_per_cell",
                             paste0(cell_type, ".csv")))
  all_pathways[[cell_type]] <- fgsea_res
  audit_rows[[length(audit_rows) + 1L]] <- list(
    cell_type = cell_type, n_genes_ranked = length(rank_vec),
    n_pathways_sig = n_sig, status = "ok"
  )
}

# ----- combined long table ----------------------------------------------
if (length(all_pathways) > 0L) {
  combined <- dplyr::bind_rows(all_pathways)
  readr::write_csv(combined,
                   file.path(args$out_dir, "pathways_long.csv"))
  cat(sprintf("\n[3] wrote %s (%d rows)\n",
              file.path(args$out_dir, "pathways_long.csv"), nrow(combined)))
}

# ----- audit ------------------------------------------------------------
audit_df <- dplyr::bind_rows(lapply(audit_rows, as.data.frame))
readr::write_csv(audit_df, file.path(args$out_dir, "audit.csv"))

cat("\n=== audit summary ===\n")
print(audit_df, n = Inf)
cat("\nDone.\n")
