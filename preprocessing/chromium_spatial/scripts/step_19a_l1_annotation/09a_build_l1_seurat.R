#!/usr/bin/env Rscript
# 09a_build_l1_seurat.R — Build a per-compartment Seurat object for L1 SingleR.
#
# Lightweight reassembly from the flat intermediate directory written by
# 08a_export_h5ad_intermediate.py. Reuses the build pattern from
# scripts/08_compartment_integration_report.Rmd::build_seurat (chunk at
# 08_compartment_integration_report.Rmd:406-536) without the visualization
# chunks — we only need counts + normalized data + obs metadata + the patched
# leiden_n30_r1.0 column + cluster_annotation.
#
# Inputs
#   --intermediate-dir   integration intermediate dir
#                          counts.mtx.gz + cells.tsv + genes.tsv + obs.csv
#   --patches-dir        dir with extras_n30*.csv patches (expects
#                          leiden_n30_r1.0 to be present in one of them)
#   --contamination-csv  cell_compartments.csv (supplies cluster_annotation)
#   --cell-metadata-csv  assembly_step15/cell_metadata.csv
#   --compartment        Epithelial | Stromal | Immune
#   --out                output RDS path
#   --dry-run            validate inputs, exit before reading counts
#
# Output
#   A Seurat object with RNA counts + normalized data + `leiden_n30_r1.0` +
#   `cluster_annotation` + standard sample/patient/position columns.

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(Matrix)
})

# Parse args ------------------------------------------------------------------
parser <- ArgumentParser(description = "Build per-compartment Seurat for L1 SingleR")
parser$add_argument("--intermediate-dir", required = TRUE)
parser$add_argument("--patches-dir", required = TRUE)
parser$add_argument("--contamination-csv", required = TRUE)
parser$add_argument("--cell-metadata-csv", required = TRUE)
parser$add_argument("--compartment", required = TRUE)
parser$add_argument("--out", required = TRUE)
parser$add_argument("--dry-run", action = "store_true", default = FALSE)
args <- parser$parse_args()

cat("=== 09a_build_l1_seurat ===\n")
cat("Compartment:     ", args$compartment, "\n")
cat("Intermediate:    ", args$intermediate_dir, "\n")
cat("Patches dir:     ", args$patches_dir, "\n")
cat("Contamination:   ", args$contamination_csv, "\n")
cat("Cell metadata:   ", args$cell_metadata_csv, "\n")
cat("Output RDS:      ", args$out, "\n")
cat("Dry run:         ", args$dry_run, "\n\n")

COUNTS_MTX  <- file.path(args$intermediate_dir, "counts.mtx.gz")
CELLS_TSV   <- file.path(args$intermediate_dir, "cells.tsv")
GENES_TSV   <- file.path(args$intermediate_dir, "genes.tsv")
OBS_CSV     <- file.path(args$intermediate_dir, "obs.csv")

# Validation preamble ---------------------------------------------------------
required <- c(
  counts_mtx     = COUNTS_MTX,
  cells_tsv      = CELLS_TSV,
  genes_tsv      = GENES_TSV,
  obs_csv        = OBS_CSV,
  contamination  = args$contamination_csv,
  cell_metadata  = args$cell_metadata_csv,
  patches_dir    = args$patches_dir
)

cat("=== Input validation ===\n")
missing_items <- character(0)
for (nm in names(required)) {
  f <- required[[nm]]
  if (nm == "patches_dir") {
    ok <- dir.exists(f)
  } else {
    ok <- file.exists(f)
  }
  cat(sprintf("  [%s] %-15s %s\n", if (ok) "OK " else "MISS", nm, f))
  if (!ok) missing_items <- c(missing_items, nm)
}
if (length(missing_items) > 0) {
  stop("Missing required inputs: ", paste(missing_items, collapse = ", "))
}

dir.create(dirname(args$out), recursive = TRUE, showWarnings = FALSE)

if (isTRUE(args$dry_run)) {
  cat("\nVALIDATION PASSED (dry run) — exiting before heavy loads\n")
  quit(status = 0)
}

# Load intermediate -----------------------------------------------------------
# Pattern source: 08_compartment_integration_report.Rmd:159-195
cat("\n=== Loading intermediate ===\n")
cell_ids <- readLines(CELLS_TSV)
gene_ids <- readLines(GENES_TSV)
cat("Loaded", length(cell_ids), "cell IDs and", length(gene_ids), "gene IDs\n")

counts_mat <- Matrix::readMM(gzfile(COUNTS_MTX))
counts_mat <- as(counts_mat, "CsparseMatrix")
counts_mat <- Matrix::t(counts_mat)  # MTX is cells x genes; Seurat wants genes x cells
rownames(counts_mat) <- gene_ids
colnames(counts_mat) <- cell_ids
cat("Counts matrix:", nrow(counts_mat), "genes x", ncol(counts_mat), "cells\n")

# obs
obs_df <- read.csv(OBS_CSV, stringsAsFactors = FALSE, check.names = FALSE)
rownames(obs_df) <- obs_df$cell_id
obs_df$cell_id <- NULL
stopifnot(all(rownames(obs_df) == cell_ids))
cat("obs columns:", ncol(obs_df), "\n")

# Build Seurat ----------------------------------------------------------------
# Pattern source: 08_compartment_integration_report.Rmd:406-475
cat("\n=== Building Seurat ===\n")
seurat_obj <- CreateSeuratObject(
  counts       = counts_mat,
  project      = paste0("l1_", args$compartment),
  min.cells    = 0,
  min.features = 0
)
cat("Seurat:", ncol(seurat_obj), "cells x", nrow(seurat_obj), "genes\n")

cells <- colnames(seurat_obj)

# Merge h5ad obs columns
obs_aligned <- obs_df[cells, , drop = FALSE]
cols_from_obs <- setdiff(colnames(obs_aligned), colnames(seurat_obj@meta.data))
for (col in cols_from_obs) {
  seurat_obj@meta.data[[col]] <- obs_aligned[[col]]
}
cat("Added", length(cols_from_obs), "columns from h5ad obs\n")

# Assembly cell metadata (sample_id, patient_id, position)
# Pattern source: 08_compartment_integration_report.Rmd:207-218
cat("\n=== Loading assembly cell metadata ===\n")
cell_meta <- read.csv(args$cell_metadata_csv, stringsAsFactors = FALSE)
cell_id_col <- intersect(c("cell_id", "X", "barcode"), colnames(cell_meta))[1]
if (is.na(cell_id_col)) cell_id_col <- colnames(cell_meta)[1]
rownames(cell_meta) <- cell_meta[[cell_id_col]]
cell_meta_comp <- cell_meta[intersect(cells, rownames(cell_meta)), , drop = FALSE]
cat("cell_metadata matched:", nrow(cell_meta_comp), "/", length(cells), "cells\n")
asm_cols <- setdiff(colnames(cell_meta_comp),
                    c(cell_id_col, colnames(seurat_obj@meta.data)))
for (col in asm_cols) {
  seurat_obj@meta.data[[col]] <- cell_meta_comp[cells, col]
}
cat("Added", length(asm_cols), "columns from cell_metadata\n")

# Contamination CSV — source of truth for cluster_annotation + compartment
# Pattern source: 08_compartment_integration_report.Rmd:226-249, 445-455
cat("\n=== Loading contamination CSV ===\n")
contam <- read.csv(args$contamination_csv, stringsAsFactors = FALSE)
contam_id_col <- intersect(c("cell_id", "X", "barcode"), colnames(contam))[1]
if (is.na(contam_id_col)) contam_id_col <- colnames(contam)[1]
rownames(contam) <- contam[[contam_id_col]]
contam_comp <- contam[intersect(cells, rownames(contam)), , drop = FALSE]
cat("contamination matched:", nrow(contam_comp), "/", length(cells), "cells\n")

priority_cols <- intersect(
  c("cluster_annotation", "compartment", "leiden_1.0"),
  colnames(contam_comp)
)
other_cols <- setdiff(colnames(contam_comp),
                      c(contam_id_col, colnames(seurat_obj@meta.data)))
contam_cols <- unique(c(priority_cols, other_cols))
for (col in contam_cols) {
  seurat_obj@meta.data[[col]] <- contam_comp[cells, col]
}
cat("Added", length(contam_cols), "columns from contamination CSV\n")

# Sanity: compartment column should match args$compartment for every cell
if ("compartment" %in% colnames(seurat_obj@meta.data)) {
  n_mismatch <- sum(seurat_obj$compartment != args$compartment, na.rm = TRUE)
  if (n_mismatch > 0) {
    cat("WARNING:", n_mismatch, "cells have compartment != ", args$compartment, "\n")
  }
}

# Patches: leiden_n30_r1.0 comes from patches/extras_n30_r10.csv
# Pattern source: 08_compartment_integration_report.Rmd:367-396, 476-486
cat("\n=== Loading patches ===\n")
patch_files <- list.files(args$patches_dir, pattern = "\\.csv$", full.names = TRUE)
cat("Patch files:", length(patch_files), "\n")

patch_df <- NULL
for (pf in patch_files) {
  pdf <- read.csv(pf, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell_id" %in% colnames(pdf)) {
    cat("  SKIP", basename(pf), "- no cell_id column\n")
    next
  }
  new_cols <- setdiff(colnames(pdf), "cell_id")
  cat("  +", basename(pf), ":", length(new_cols), "cols\n")
  if (is.null(patch_df)) {
    patch_df <- pdf
  } else {
    patch_df <- merge(patch_df, pdf, by = "cell_id", all = TRUE)
  }
}

if (!is.null(patch_df)) {
  m <- match(cells, patch_df$cell_id)
  patch_cols <- setdiff(colnames(patch_df), "cell_id")
  for (col in patch_cols) {
    seurat_obj@meta.data[[col]] <- patch_df[[col]][m]
  }
  cat("Patch columns merged:", length(patch_cols), "\n")
}

# Verify leiden_n30_r1.0 is present — this is the vote column downstream
if (!"leiden_n30_r1.0" %in% colnames(seurat_obj@meta.data)) {
  stop("leiden_n30_r1.0 not found after patch merge — check patches dir contents")
}
n_na_leiden <- sum(is.na(seurat_obj$leiden_n30_r1.0))
cat("leiden_n30_r1.0: ", length(unique(na.omit(seurat_obj$leiden_n30_r1.0))),
    " clusters, ", n_na_leiden, " NAs\n", sep = "")

# Normalize so GetAssayData(layer='data') works downstream
cat("\n=== Normalizing ===\n")
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
seurat_obj <- JoinLayers(seurat_obj)

# Save ------------------------------------------------------------------------
cat("\n=== Saving ===\n")
saveRDS(seurat_obj, args$out)
cat("Saved:", args$out, "\n")
cat("  cells:", ncol(seurat_obj), "\n")
cat("  genes:", nrow(seurat_obj), "\n")
cat("  leiden_n30_r1.0 clusters:",
    length(unique(na.omit(seurat_obj$leiden_n30_r1.0))), "\n")
cat("  cluster_annotation values:",
    length(unique(na.omit(seurat_obj$cluster_annotation))), "\n")

# Sidecar for L2 scaffold worker (10d_aggregate_l2.py). Two-column CSV with
# cell_id + leiden_n30_r1.0. Written alongside the Seurat RDS in the same
# compartment outputs dir (outputs/annotation/l1/{compartment}/).
leiden_sidecar <- file.path(dirname(args$out), "leiden_n30_r1.0.csv")
sidecar_df <- data.frame(
  cell_id         = colnames(seurat_obj),
  leiden_n30_r1.0 = as.character(seurat_obj$leiden_n30_r1.0),
  stringsAsFactors = FALSE,
  check.names      = FALSE
)
write.csv(sidecar_df, leiden_sidecar, row.names = FALSE)
cat("Saved leiden sidecar:", leiden_sidecar, "\n")
cat("  rows:", nrow(sidecar_df), "\n")

# Validation gate: sidecar must exist and match Seurat cell count
if (!file.exists(leiden_sidecar)) {
  stop("Leiden sidecar write failed: ", leiden_sidecar)
}
if (nrow(sidecar_df) != ncol(seurat_obj)) {
  stop(sprintf("Leiden sidecar row count %d != Seurat cells %d",
               nrow(sidecar_df), ncol(seurat_obj)))
}
cat("Sidecar validation: OK\n")

cat("Done.\n")
