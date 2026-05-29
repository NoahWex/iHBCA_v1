#!/usr/bin/env Rscript
# Step 18 — DESeq2 pseudobulk per cluster resolution
#
# Pseudobulk design: ~ patient_id + position + cluster (LRT vs ~ patient_id + position)
# Per-cluster LFC extracted from Wald test coefficients after LRT gene selection.
# (cluster × sample) groups with < --min-cells are excluded from pseudobulk.
# Clusters represented in < 2 patients are excluded.
#
# Inputs:  integration_intermediate/{compartment}/scvi_n100/{counts,genes,cells,obs}
# Outputs: outputs/deseq/{compartment}/deseq_res_{r}.csv per resolution

suppressPackageStartupMessages({
  library(DESeq2)
  library(BiocParallel)
  library(Matrix)
  library(dplyr)
  library(argparse)
})

parser <- ArgumentParser()
parser$add_argument("--compartment",      required = TRUE)
parser$add_argument("--integration-dir",  required = TRUE,
  help = "Path to integration_intermediate/{compartment}/scvi_n100/")
parser$add_argument("--out-dir",          required = TRUE)
parser$add_argument("--resolution",       required = TRUE,
  help = "Single Leiden resolution to process (e.g. 0.3)")
parser$add_argument("--min-cells",        type = "integer", default = 10,
  help = "Min cells per (cluster x sample) for pseudobulk inclusion [default: 10]")
args <- parser$parse_args()

# ---------------------------------------------------------------------------
# Load canonical intermediates
# ---------------------------------------------------------------------------
cat(sprintf("Compartment: %s\n", args$compartment))
cat(sprintf("Integration dir: %s\n", args$integration_dir))

obs    <- read.csv(file.path(args$integration_dir, "obs.csv"))
genes  <- read.table(file.path(args$integration_dir, "genes.tsv"), header = FALSE)$V1
cells  <- read.table(file.path(args$integration_dir, "cells.tsv"), header = FALSE)$V1
counts <- readMM(gzcon(file(file.path(args$integration_dir, "counts.mtx.gz"), "rb")))
# counts.mtx.gz is stored cells x genes; transpose to genes x cells for DESeq2
counts <- t(counts)
rownames(counts) <- genes
colnames(counts) <- cells

stopifnot(all(cells == obs$cell_id))
cat(sprintf("Loaded: %d cells x %d genes\n", ncol(counts), nrow(counts)))

leiden_col <- paste0("leiden_", args$resolution)
if (!leiden_col %in% colnames(obs))
  stop(sprintf("Column '%s' not found in obs.csv", leiden_col))
cat(sprintf("Resolution: %s\n", args$resolution))
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Per-resolution DESeq2
# ---------------------------------------------------------------------------
run_deseq <- function(obs, counts, leiden_col, min_cells, out_dir) {
  res_label <- sub("leiden_", "", leiden_col)
  cat(sprintf("\n=== Resolution %s ===\n", res_label))

  obs$cluster    <- factor(as.character(obs[[leiden_col]]))
  obs$patient_id <- factor(as.character(obs$patient_id))
  obs$position   <- factor(as.character(obs$position))
  obs$sample_id  <- as.character(obs$sample_id)

  groups <- paste(obs$cluster, obs$sample_id, sep = "__")

  # Filter (cluster x sample) groups below min_cells threshold
  keep_groups <- names(table(groups))[table(groups) >= min_cells]
  mask        <- groups %in% keep_groups
  cat(sprintf("  (cluster x sample) groups: %d kept / %d total (min %d cells)\n",
              length(keep_groups), length(unique(groups)), min_cells))

  obs_f    <- obs[mask, ]
  counts_f <- counts[, mask, drop = FALSE]
  groups_f <- groups[mask]

  # Aggregate counts per (cluster x sample) via sparse matrix multiply
  # Builds a cells × groups indicator, then: (genes × cells) %*% (cells × groups)
  group_fac <- factor(groups_f, levels = keep_groups)
  agg_ind   <- fac2sparse(group_fac)          # groups × cells (sparse)
  pb_mat    <- counts_f %*% t(agg_ind)        # genes × groups
  colnames(pb_mat) <- keep_groups

  # Pseudobulk metadata
  pb_meta <- do.call(rbind, lapply(keep_groups, function(g) {
    i <- which(groups_f == g)[1]
    data.frame(group = g, cluster = obs_f$cluster[i],
               patient_id = obs_f$patient_id[i],
               position   = obs_f$position[i],
               n_cells    = sum(groups_f == g),
               stringsAsFactors = FALSE)
  }))
  rownames(pb_meta) <- keep_groups

  # Require clusters present in >= 2 patients
  n_patients <- tapply(pb_meta$patient_id, pb_meta$cluster, function(x) length(unique(x)))
  keep_cl    <- names(n_patients)[n_patients >= 2]
  if (length(keep_cl) < 2) {
    cat(sprintf("  Skipping: only %d cluster(s) with >= 2 patients\n", length(keep_cl)))
    return(invisible(NULL))
  }
  keep_cols <- rownames(pb_meta)[pb_meta$cluster %in% keep_cl]
  pb_meta   <- pb_meta[keep_cols, ]
  pb_mat    <- pb_mat[, keep_cols, drop = FALSE]
  pb_meta$cluster    <- droplevels(factor(pb_meta$cluster))
  pb_meta$patient_id <- droplevels(factor(pb_meta$patient_id))
  pb_meta$position   <- droplevels(factor(pb_meta$position))
  cat(sprintf("  %d clusters (>= 2 patients), %d pseudobulk samples\n",
              nlevels(pb_meta$cluster), nrow(pb_meta)))

  # DESeq2 LRT: identify genes with cluster-specific variation
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_mat),
    colData   = pb_meta,
    design    = ~ patient_id + position + cluster
  )
  dds <- dds[rowSums(counts(dds) >= 5) >= 2, ]  # basic gene filter

  tryCatch({
    dds <- DESeq(dds, test = "LRT", reduced = ~ patient_id + position,
                 fitType = "parametric", quiet = TRUE)

    # Extract per-cluster LFC from Wald coefficients (cluster vs. reference level)
    ref_cl  <- levels(pb_meta$cluster)[1]
    coef_names <- resultsNames(dds)
    cl_coefs <- grep("^cluster_", coef_names, value = TRUE)

    # LRT results (global cluster effect per gene)
    lrt_res <- as.data.frame(results(dds, alpha = 0.05))
    lrt_res$gene <- rownames(lrt_res)

    # Per-cluster LFC (parallel across cluster coefficients)
    n_workers <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
    cl_lfc <- bplapply(cl_coefs, function(cn) {
      cl_name <- sub("^cluster_", "", sub("_vs_.*$", "", cn))
      r <- as.data.frame(lfcShrink(dds, coef = cn, type = "normal", quiet = TRUE))
      data.frame(gene = rownames(r), cluster = cl_name,
                 lfc = r$log2FoldChange, lfc_se = r$lfcSE)
    }, BPPARAM = MulticoreParam(workers = n_workers))
    # Reference cluster: LFC = 0 by definition
    ref_df <- data.frame(gene = rownames(lrt_res), cluster = ref_cl, lfc = 0, lfc_se = NA)
    cl_lfc_df <- rbind(do.call(rbind, cl_lfc), ref_df)

    # Join: LRT significance + per-cluster LFC
    out <- merge(lrt_res[, c("gene", "baseMean", "stat", "pvalue", "padj")],
                 cl_lfc_df, by = "gene", all.x = TRUE)
    out$resolution  <- res_label
    out$compartment <- args$compartment
    out <- out[order(out$padj, na.last = TRUE), ]

    out_file <- file.path(out_dir, sprintf("deseq_res_%s.csv", res_label))
    write.csv(out, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s (%d gene x cluster rows)\n", out_file, nrow(out)))
    rm(dds); gc()

  }, error = function(e) {
    cat(sprintf("  ERROR at res %s: %s\n", res_label, conditionMessage(e)))
  })
}

run_deseq(obs, counts, leiden_col, args$min_cells, args$out_dir)

cat("\nDone.\n")
