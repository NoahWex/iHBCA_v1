#!/usr/bin/env Rscript
# Build per-cell logFC projection for the BR1_vs_AR_tested Milo DA context,
# in the format expected by the publication UMAP renderer
# (render_fig2_umap_logfc.py).
#
# Method: for each cell, take the mean logFC across the neighborhoods that
# include it, restricted to neighborhoods with SpatialFDR < 0.05. Expressed
# as a sparse matvec:
#     weighted_logfc = nh %*% (logFC * sig)
#     cell_logfc     = weighted_logfc / pmax(rowSums(nh[, sig]), 1)
#
# Output: cell_projection.csv with columns
#   cell_id, UMAP1, UMAP2, study, donor_id, BR1_vs_AR__cell_logfc
# plus an empty pxr_multi_context.csv stub (label channel is disabled at
# render time via --no-ng-labels).
#
# study + donor_id are extracted directly from colData(milo) to support
# per-study faceted UMAP renders (s_ihbca_umap_brca1_by_study_tested) without
# needing cell_id prefix parsing (fragile across studies with different ID
# formats). Same Option-A pattern as build_ihbca_fig2_aucell_scores.py.

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(Matrix); library(dplyr); library(readr)
})

parser <- ArgumentParser()
parser$add_argument("--project-root", required = TRUE,
                    help = "iHBCA_publication root (path resolver lives at publication/config/load_paths.R)")
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--milo-rds",
                    help = "Optional override for the Milo RDS path (default: resolve via paths.yaml spatial_milo)")
parser$add_argument("--patch-r",
                    help = "Optional override for the miloR patch script")
parser$add_argument("--da-csv",
                    help = "Optional override for the per-nhood DA results CSV (default: BR1_vs_AR_tested track of risk_main_effects_20260507)")
parser$add_argument("--context-name", default = "BR1_vs_AR",
                    help = "Output column suffix: writes '<context_name>__cell_logfc'. Default 'BR1_vs_AR' (for both tested and full cohorts; use 'parity_in_AR' for parity).")
parser$add_argument("--fdr-threshold", type = "double", default = 0.05)
args <- parser$parse_args()

# %||% is not exported by argparse; define locally (must be set before first use).
`%||%` <- function(a, b) if (is.null(a)) b else a

source(file.path(args$project_root, "publication", "config", "load_paths.R"))
paths <- load_paths(file.path(args$project_root, "publication", "config"))

milo_rds <- args$milo_rds %||% resolve_path("spatial_milo", paths)
da_csv   <- args$da_csv %||% file.path(
  paths$sources$ihbca_v1_abundance,
  "da_pipeline/inquiries/risk_main_effects_20260507/outputs/stageD_da_results/BR1_vs_AR_tested/da_results.csv"
)
# patch_milor.R lives alongside the Milo build pipeline; default to the
# Spatial HBCA DA scripts location resolved from the user root.
patch_r <- args$patch_r %||% file.path(
  paths$.root, "Spatial_HBCA/project/04_DifferentialAbundance/da_v1/scripts/patch_milor.R"
)

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

cat("[1] load milo + patch ...\n")
source(patch_r)
milo <- readRDS(milo_rds)
nh <- nhoods(milo)
cat(sprintf("  cells=%d nhoods=%d\n", nrow(nh), ncol(nh)))

cat("[2] UMAP coords (UMAP_scVI) + colData study/donor_id\n")
umap_key <- intersect(c("UMAP_scVI", "UMAP", "X_umap", "umap"),
                      reducedDimNames(milo))[1]
umap_df <- as.data.frame(reducedDim(milo, umap_key))
colnames(umap_df)[1:2] <- c("UMAP1", "UMAP2")
umap_df$cell_id <- colnames(milo)

# Extract study + donor_id from colData (Option-A pattern: bake metadata
# into substrate; never parse cell_id at render time).
cd <- colData(milo)
cd_cols <- colnames(cd)
study_col <- intersect(c("study", "Study", "study_id"), cd_cols)[1]
donor_col <- intersect(c("donor_id", "patientID", "ihbca_donor_id",
                         "patient", "Sample", "sample"), cd_cols)[1]
if (is.na(study_col)) stop("No study column in colData(milo); colData has: ",
                            paste(cd_cols, collapse = ", "))
if (is.na(donor_col)) stop("No donor column in colData(milo); colData has: ",
                            paste(cd_cols, collapse = ", "))
cat(sprintf("  colData: study from '%s', donor from '%s'\n", study_col, donor_col))
umap_df$study <- as.character(cd[[study_col]])
umap_df$donor_id <- as.character(cd[[donor_col]])
cat(sprintf("  unique studies: %d (%s)\n", length(unique(umap_df$study)),
            paste(sort(unique(umap_df$study)), collapse = ", ")))
cat(sprintf("  unique donors: %d\n", length(unique(umap_df$donor_id))))

cat("[3] load DA results\n")
da <- read_csv(da_csv, show_col_types = FALSE)
cat(sprintf("  rows=%d cols=%s\n", nrow(da),
            paste(head(colnames(da), 8), collapse = ", ")))

# Sort by Nhood index so order matches nhoods() cols
da <- da[order(da$Nhood), ]
stopifnot(nrow(da) == ncol(nh))

cat("[4] per-cell logFC projection (sig nhoods only)\n")
sig <- !is.na(da$SpatialFDR) & da$SpatialFDR < args$fdr_threshold
cat(sprintf("  sig nhoods (FDR<%.2f): %d / %d\n",
            args$fdr_threshold, sum(sig), length(sig)))

lfc_vec <- da$logFC
lfc_vec[!sig] <- 0
weighted_logfc <- as.numeric(nh %*% lfc_vec)
sig_counts <- as.numeric(nh %*% as.numeric(sig))
cell_logfc <- ifelse(sig_counts > 0, weighted_logfc / sig_counts, NA_real_)
cat(sprintf("  cells with non-NA cell_logfc: %d / %d (%.1f%%)\n",
            sum(!is.na(cell_logfc)), length(cell_logfc),
            100 * sum(!is.na(cell_logfc)) / length(cell_logfc)))

cat("[5] write cell_projection.csv\n")
# Column name matches the contract expected by render_fig2_umap_logfc.py
# (its --context flag selects which column to colorize). The renderer's
# --context BR1_vs_AR selects the column built here even though the
# underlying data is the tested cohort — distinguished by --panel-id at
# render time (fig2_ihbca_umap_brca1_tested vs fig2_ihbca_umap_brca1_full).
out <- umap_df
context_col <- paste0(args$context_name, "__cell_logfc")
out[[context_col]] <- cell_logfc
out_path <- file.path(args$out_dir, "cell_projection.csv")
cat(sprintf("  context column: %s\n", context_col))
write_csv(out, out_path)
cat(sprintf("  wrote: %s (%.1f MB, %d rows)\n",
            out_path, file.info(out_path)$size / 1e6, nrow(out)))

# Empty pxr_multi_context.csv stub required by the renderer's input contract;
# the label channel is disabled at render time via --no-ng-labels.
stub <- data.frame(pxr_ng_id = integer(0), UMAP_1 = numeric(0),
                   UMAP_2 = numeric(0), context = character(0))
write_csv(stub, file.path(args$out_dir, "pxr_multi_context.csv"))

cat("\nDone.\n")
