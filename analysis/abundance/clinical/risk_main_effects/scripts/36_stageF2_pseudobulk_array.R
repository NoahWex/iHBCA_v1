#!/usr/bin/env Rscript
# 36_stageF2_pseudobulk_array.R
# Stage F.2 — pseudobulk per (nhoodgroup × donor) per contrast.
# SLURM array task: ONE contrast per task. Reads stage F.1's nhood_groups.csv
# and aggregates raw counts per (NhoodGroup_renamed × donor) by summing
# across cells whose nhood membership puts them in that group.
#
# Cell-level group assignment: a cell can belong to multiple nhoods. We
# assign each cell to a NhoodGroup if it is contained in any nhood whose
# group id matches that group. To avoid double-counting genes, each cell
# contributes to at most one group: the group containing the largest
# fraction of nhoods covering that cell. (This matches the existing
# stage 13 pseudobulk pattern.)
#
# Output (per contrast):
#   stageF2_pseudobulk/<contrast>/pseudobulk.rds
#     SummarizedExperiment with assays$counts (gene × group×donor),
#     colData has NhoodGroup_renamed, donor_id, study, parent_L2_joint, etc.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "36_stageF2_pseudobulk_array.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
if (TASK_ID < 1) stop("SLURM_ARRAY_TASK_ID must be >= 1")
meta <- yaml::read_yaml(file.path(paths$outputs$stageC, "experiment_metadata.yaml"))
if (TASK_ID > length(meta$contrasts)) {
  cat(sprintf("TASK_ID %d > #contrasts %d — exiting\n",
              TASK_ID, length(meta$contrasts)))
  quit(status = 0)
}
con <- meta$contrasts[[TASK_ID]]

out_sub <- file.path(paths$outputs$stageF2, con$name)
ensure_dir(out_sub)
out_rds <- file.path(out_sub, "pseudobulk.rds")
if (file.exists(out_rds) && file.size(out_rds) > 1000) {
  cat("Output exists, skipping:", out_rds, "\n")
  quit(status = 0)
}

cat(sprintf("=== Stage F.2: pseudobulk for '%s' ===\n", con$name))

source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)

ng <- read_csv(file.path(paths$outputs$stageF1, con$name, "nhood_groups.csv"),
               show_col_types = FALSE)
cat(sprintf("milo cells=%d nhoods=%d, ng-rows=%d, groups=%d\n",
            ncol(milo), ncol(nhoods(milo)), nrow(ng),
            n_distinct(ng$NhoodGroup_renamed, na.rm = TRUE)))

# Cell × nhood adjacency (sparse, cells in rows, nhoods in columns)
nh_mat <- nhoods(milo)
# Cell -> nhoods it belongs to
# For each cell, pick the dominant group (most nhoods covering it that share
# a group id).
ng_lookup <- ng %>%
  filter(!is.na(NhoodGroup_renamed)) %>%
  select(Nhood, NhoodGroup_renamed)

# nhood index -> group id vector (NA if not assigned)
group_per_nhood <- rep(NA_character_, ncol(nh_mat))
group_per_nhood[ng_lookup$Nhood] <- ng_lookup$NhoodGroup_renamed

# For each cell: tally # nhoods in each group; pick max.
# nh_mat is cells × nhoods; we need cells × groups.
group_levels <- sort(unique(na.omit(group_per_nhood)))
if (length(group_levels) == 0) {
  cat(sprintf("No NhoodGroup assignments for '%s' — F.1 produced no groups for this contrast (typical for stratified within-HR cohorts). Skipping F.2 gracefully.\n",
              con$name))
  # Touch a sentinel file so downstream stages can distinguish "skipped, no
  # groups" from "missing because F.2 didn't run".
  writeLines(sprintf("skipped: no NhoodGroup assignments for contrast '%s'", con$name),
              file.path(out_sub, "SKIPPED_NO_GROUPS.txt"))
  quit(status = 0)
}

cat("Building cell × group counts...\n")
group_idx <- match(group_per_nhood, group_levels)
# columns in nh_mat NOT in any group get NA in group_idx; drop them
keep_nhood <- !is.na(group_idx)
nh_sub <- nh_mat[, keep_nhood, drop = FALSE]
group_idx_sub <- group_idx[keep_nhood]
# Build a sparse aggregation matrix nhoods × groups
ngrp <- length(group_levels)
agg_M <- sparseMatrix(
  i = seq_along(group_idx_sub),
  j = group_idx_sub,
  x = 1,
  dims = c(length(group_idx_sub), ngrp)
)
# cells × groups: number of nhoods of each group covering each cell
cell_group_ct <- nh_sub %*% agg_M
cat(sprintf("cell×group dim: %s\n", paste(dim(cell_group_ct), collapse=" x ")))

# Each cell's dominant group: argmax (ties broken by first)
cell_group_ct_dense <- as.matrix(cell_group_ct)
cell_group_max_idx <- max.col(cell_group_ct_dense, ties.method = "first")
cell_group_max_n   <- cell_group_ct_dense[
  cbind(seq_len(nrow(cell_group_ct_dense)), cell_group_max_idx)
]
cell_assigned_group <- ifelse(cell_group_max_n > 0,
                              group_levels[cell_group_max_idx],
                              NA_character_)
cat(sprintf("cells assigned to a group: %d / %d\n",
            sum(!is.na(cell_assigned_group)), length(cell_assigned_group)))

# donor id
cd <- as.data.frame(colData(milo))
donor_col <- intersect(c("ihbca_donor_id", "donor_id", "donor"), colnames(cd))[1]
study_col <- intersect(c("study", "dataset", "source"), colnames(cd))[1]

# Aggregate counts to (group × donor) by delegating to the annotation pipeline's
# Python pseudobulk aggregator, which iterates the 3 compartment NPZ files
# (epi/imm/str) sequentially. The substrate has a single canonical gene order
# (via gene_data.csv) — no per-study symbol/Ensembl mismatch. Cell-level groupby
# is provided via a per-contrast aligned CSV keyed on global_numeric_id (= atlas
# obs row position). Same pattern that produced V1 leiden pseudobulks.

# 1. Build atlas obs/cellID → global_numeric_id map (positional, 0-indexed).
suppressPackageStartupMessages(library(rhdf5))
if (is.null(paths$inputs$atlas_h5ad)) {
  stop("paths$inputs$atlas_h5ad not configured")
}
cat("Reading atlas obs/cellID for numeric_id bridge...\n")
atlas_cellID <- as.character(rhdf5::h5read(paths$inputs$atlas_h5ad, "obs/cellID"))
n_atlas <- length(atlas_cellID)
cat(sprintf("  atlas n_cells: %d\n", n_atlas))

milo_cells <- colnames(milo)
milo_to_atlas_pos <- match(milo_cells, atlas_cellID) - 1L  # 0-indexed
n_unmatched <- sum(is.na(milo_to_atlas_pos))
if (n_unmatched > 0) {
  cat(sprintf("  WARNING: %d milo cells not found in atlas\n", n_unmatched))
}

# 2. Write per-contrast aligned CSV with numeric_id + NhoodGroup_renamed +
#    donor + study (all columns the Python script's groupby/donor/study keys
#    will use). Cells without a group assignment are excluded.
aligned_dir <- file.path(out_sub, "aligned")
ensure_dir(aligned_dir)
aligned_csv <- file.path(aligned_dir, "cell_to_nhoodgroup.csv")

donor_per_milo <- as.character(cd[[donor_col]])
study_per_milo <- as.character(cd[[study_col]])

aligned_df <- tibble(
  numeric_id = milo_to_atlas_pos,
  NhoodGroup_renamed = cell_assigned_group,
  patientID = donor_per_milo,
  dataset = study_per_milo
) %>% filter(!is.na(numeric_id), !is.na(NhoodGroup_renamed))

cat(sprintf("aligned table: %d rows (cells with group + atlas match)\n",
            nrow(aligned_df)))
write_csv(aligned_df, aligned_csv)

# 3. Shell out to the annotation pipeline's pseudobulk aggregator. It loads
#    the 3 compartment NPZs sequentially, joins our aligned CSV by
#    global_numeric_id, and emits CSV outputs.
components_root <- paths$inputs$components_root
agg_py <- paths$inputs$pseudobulk_aggregate_py
if (is.null(components_root) || !dir.exists(components_root)) {
  stop("paths$inputs$components_root missing: ", components_root %||% "<NULL>")
}
if (is.null(agg_py) || !file.exists(agg_py)) {
  stop("paths$inputs$pseudobulk_aggregate_py missing: ", agg_py %||% "<NULL>")
}

py_outdir <- file.path(out_sub, "py_pseudobulk")
ensure_dir(py_outdir)
py_tag <- sprintf("parity_%s_nhoodgroup", con$name)

py_cmd <- sprintf(
  "python3 -u %s --epi-npz %s --epi-meta %s --imm-npz %s --imm-meta %s --str-npz %s --str-meta %s --gene-data-csv %s --leiden-aligned %s --groupby-col NhoodGroup_renamed --donor-key patientID --study-key dataset --min-cells 10 --tag %s --output-dir %s",
  shQuote(agg_py),
  shQuote(file.path(components_root, "epi_counts.npz")),
  shQuote(file.path(components_root, "epi_metadata_enriched.csv")),
  shQuote(file.path(components_root, "imm_counts.npz")),
  shQuote(file.path(components_root, "imm_metadata_enriched.csv")),
  shQuote(file.path(components_root, "str_counts.npz")),
  shQuote(file.path(components_root, "str_metadata_enriched.csv")),
  shQuote(file.path(components_root, "gene_data.csv")),
  shQuote(aligned_csv),
  shQuote(py_tag),
  shQuote(py_outdir)
)
cat("Calling Python aggregator:\n  ", py_cmd, "\n")
t_py <- proc.time()
py_status <- system(py_cmd)
cat(sprintf("Python aggregator exit=%d, %ds elapsed\n",
            py_status, as.integer((proc.time() - t_py)[["elapsed"]])))
if (py_status != 0) stop("Python pseudobulk aggregator failed (exit=", py_status, ")")

# 4. Read CSV outputs back, build SE so F.3 schema is unchanged.
counts_csv <- file.path(py_outdir, sprintf("%s_pseudobulk_counts.csv", py_tag))
meta_csv   <- file.path(py_outdir, sprintf("%s_pseudobulk_meta.csv", py_tag))
if (!file.exists(counts_csv) || !file.exists(meta_csv)) {
  stop("Expected Python outputs missing: ", counts_csv, " / ", meta_csv)
}
cat(sprintf("Reading Python pseudobulk: %s\n", counts_csv))
pb_df <- read_csv(counts_csv, show_col_types = FALSE)
sample_meta <- read_csv(meta_csv, show_col_types = FALSE)

# Counts CSV from python script is samples × genes (rows = sample_id, cols = genes).
# Convert to genes × samples for SummarizedExperiment.
sample_id_col <- intersect(c("sample_id", "Unnamed: 0", "...1"), colnames(pb_df))[1]
if (is.na(sample_id_col)) sample_id_col <- colnames(pb_df)[1]
sample_ids <- pb_df[[sample_id_col]]
pb_mat <- as.matrix(pb_df[, setdiff(colnames(pb_df), sample_id_col), drop = FALSE])
rownames(pb_mat) <- sample_ids
pb_counts <- t(pb_mat)  # genes × samples
cat(sprintf("pseudobulk dim (genes × samples): %s\n",
            paste(dim(pb_counts), collapse=" x ")))

# Sample id format from Python aggregator: "<NhoodGroup>__<donor>"
gd_levels <- colnames(pb_counts)

# Build colData
cd_pb <- tibble(
  group_donor = gd_levels,
  NhoodGroup_renamed = sub("__.*$", "", gd_levels),
  donor = sub("^[^_]*__|^.+?__", "", gd_levels)
)
# donor extraction: gd_levels are "<group>__<donor>"; donor is everything after first "__"
cd_pb$donor <- vapply(strsplit(cd_pb$group_donor, "__", fixed = TRUE),
                       function(x) paste(x[-1], collapse = "__"), character(1))

# Per-donor study lookup
donor_study <- cd %>%
  distinct(.data[[donor_col]], .data[[study_col]]) %>%
  rename(donor = !!donor_col, study = !!study_col)
cd_pb <- cd_pb %>% left_join(donor_study, by = "donor")

# Attach parent L2 from group summary
gs <- read_csv(file.path(paths$outputs$stageF1, con$name, "nhood_groups_summary.csv"),
               show_col_types = FALSE)
cd_pb <- cd_pb %>%
  left_join(gs %>%
              select(NhoodGroup_renamed, parent_L2_joint, parent_compartment,
                     parent_label, n_nhoods_in_group, group_med_lfc, group_pct_up),
            by = "NhoodGroup_renamed")

pb <- SummarizedExperiment(
  assays = list(counts = pb_counts),
  colData = DataFrame(cd_pb)
)
saveRDS(pb, out_rds)
cat(sprintf("Wrote: %s  (genes=%d, samples=%d)\n",
            out_rds, nrow(pb), ncol(pb)))
cat("=== Stage F.2 done ===\n")

# Explicit clean exit. Without this, some R session-finalizer (likely a dplyr
# promise being forced during garbage collection of cd_pb's lazy join object)
# raises "could not find function 'join'" after all script work has completed,
# corrupting the SLURM exit code. saveRDS has already written valid output by
# this point — the error is a phantom from the cleanup phase.
quit(status = 0)
