#!/usr/bin/env Rscript
# 35_stageF1_groupnhoods_array.R
# Stage F.1 — groupNhoods per contrast.
# SLURM array task: ONE contrast per task. SLURM_ARRAY_TASK_ID indexes
# experiment_metadata.yaml's contrasts list (same indexing as Stage D).
#
# For each contrast:
#   1. Load milo (for nhoodAdjacency)
#   2. Load Stage D da_results.csv
#   3. groupNhoods(milo, da_res, ...) using thresholds from inquiry$stageF$groupNhoods
#   4. Annotate groups with parent L2 (compartment::label) by majority-vote
#      across nhoods.
#   5. Write nhood_groups.csv with columns:
#        Nhood, NhoodGroup, NhoodGroup_renamed, parent_label, parent_compartment,
#        n_nhoods_in_group, group_med_lfc, group_pct_up
#
# Output:
#   stageF1_nhoodgroups/<contrast>/nhood_groups.csv

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(miloR)
  library(dplyr); library(readr); library(yaml); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "35_stageF1_groupnhoods_array.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
source(file.path(script_dir, "lib", "da_helpers.R"))

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

out_sub <- file.path(paths$outputs$stageF1, con$name)
ensure_dir(out_sub)
out_csv <- file.path(out_sub, "nhood_groups.csv")
if (file.exists(out_csv) && file.size(out_csv) > 100) {
  cat("Output exists, skipping:", out_csv, "\n")
  quit(status = 0)
}

cat(sprintf("=== Stage F.1: groupNhoods for contrast '%s' ===\n", con$name))

# Load milo (nhoodAdjacency only) + DA results
source(paths$inputs$patch_milor)
milo <- readRDS(paths$inputs$milo)
da <- read_csv(file.path(paths$outputs$stageD, con$name, "da_results.csv"),
               show_col_types = FALSE)
cat(sprintf("milo cells=%d, DA nhoods=%d (sig FDR<%.2f: %d)\n",
            ncol(milo), nrow(da), con$spatial_fdr,
            sum(da$SpatialFDR < con$spatial_fdr, na.rm = TRUE)))

# groupNhoods params
gn <- inquiry$stageF$groupNhoods %||% list()
da_fdr <- gn$da_fdr %||% 0.10
overlap <- gn$overlap %||% 1
subset_nh <- gn$subset_nhoods

# Adaptive max.lfc.delta. Source priority:
#   1. Stage D run_summary.yaml's adaptive_lfc_cutoff (if Stage D wrote one)
#   2. Compute from da_results via LOESS crossing (script 11 pattern)
#   3. inquiry.yaml override > 0.5 fallback
# Passing NA skips groupNhoods edge filtering → full 201K-nhood dense graph → OOM.
rs_path <- file.path(paths$outputs$stageD, con$name, "run_summary.yaml")
rs <- if (file.exists(rs_path)) yaml::read_yaml(rs_path) else list()
adaptive_lfc <- suppressWarnings(as.numeric(rs$adaptive_lfc_cutoff))
if (length(adaptive_lfc) == 0 || is.na(adaptive_lfc)) {
  cat("Stage D run_summary missing adaptive_lfc_cutoff — computing from da_results\n")
  adaptive_lfc <- compute_adaptive_lfc_cutoff(da, fdr_threshold = con$spatial_fdr,
                                                fallback_lfc = 0.5)
}
if (!is.finite(adaptive_lfc) || adaptive_lfc <= 0) {
  cat("Adaptive LFC computation failed — using inquiry override or 0.5\n")
  adaptive_lfc <- suppressWarnings(as.numeric(gn$max_lfc))
  if (length(adaptive_lfc) == 0 || is.na(adaptive_lfc) || adaptive_lfc <= 0) {
    adaptive_lfc <- 0.5
  }
}
max_lfc_delta <- max(adaptive_lfc / 2, 0.5)
cat(sprintf("adaptive_lfc=%.4f → max.lfc.delta=%.4f\n",
            adaptive_lfc, max_lfc_delta))

# Filter to candidate nhoods
keep <- !is.na(da$SpatialFDR) & da$SpatialFDR < da_fdr
n_candidates <- sum(keep)
cat(sprintf("groupNhoods candidates (FDR<%.2f): %d\n", da_fdr, n_candidates))

# Graceful skip when no DA neighbourhoods — write empty stub + exit 0 so the
# array task doesn't fail and downstream stages see a "no signal" record.
if (n_candidates == 0) {
  cat("No DA neighbourhoods at this FDR cutoff — writing empty stub and exiting.\n")
  write_csv(tibble(Nhood = integer(), NhoodGroup = integer(),
                   NhoodGroup_renamed = character(), parent_label = character(),
                   parent_compartment = character(), n_nhoods_in_group = integer(),
                   group_med_lfc = numeric(), group_pct_up = numeric(),
                   skipped = TRUE,
                   skip_reason = sprintf("no DA nhoods at FDR<%.2f", da_fdr)),
            out_csv)
  summary_csv <- file.path(out_sub, "nhood_groups_summary.csv")
  write_csv(tibble(NhoodGroup = integer(), NhoodGroup_renamed = character(),
                   parent_L2_joint = character(), parent_compartment = character(),
                   parent_label = character(), n_nhoods_in_group = integer(),
                   n_sig = integer(), group_med_lfc = numeric(),
                   group_pct_up = numeric()),
            summary_csv)
  cat(sprintf("Wrote stubs:\n  %s\n  %s\n", out_csv, summary_csv))
  cat("=== Stage F.1 task done (skipped) ===\n")
  quit(status = 0)
}

# SPARSE groupNhoods replacement (v4 pattern: clinical_da_v4/group_nhoods_patched.R).
# miloR::groupNhoods has two O(N^2) dense allocations (as.matrix(adj > 0) and the
# pairwise logFC delta matrix), each ~323 GB at N=201K nhoods → OOMs at 256 GB.
# group_nhoods_sparse stays in COO triplets, builds the igraph from the edge
# list directly, and runs Louvain — never materializes a 201K x 201K dense matrix.
source(file.path(script_dir, "lib", "group_nhoods_sparse.R"))

cat("Computing nhood overlap matrix (crossprod of nhoods)...\n")
t_ov <- proc.time()
nhood_mat <- nhoods(milo)
overlap_mat <- Matrix::crossprod(nhood_mat)  # sparse N x N cell-overlap counts
diag(overlap_mat) <- 0  # ignore self-overlap
cat(sprintf("Overlap matrix: %d x %d, %d nonzero, %.1f s\n",
            nrow(overlap_mat), ncol(overlap_mat),
            Matrix::nnzero(overlap_mat),
            (proc.time() - t_ov)[["elapsed"]]))

ng <- group_nhoods_sparse(
  overlap        = overlap_mat,
  da_results     = da,
  da.fdr         = da_fdr,
  max.lfc.delta  = max_lfc_delta,
  min_overlap    = overlap,
  merge.discord  = FALSE,
  seed           = 42L
)
cat(sprintf("group_nhoods_sparse returned: %d nhoods, %d groups\n",
            nrow(ng), n_distinct(ng$NhoodGroup, na.rm = TRUE)))

# Build per-group summary using each nhood's annotated L2
# ng already carries compartment + label from da_results (Stage D Nhood-joined).
# Normalize column names if Stage D used the L2_-prefixed variant.
if (!"compartment" %in% colnames(ng) && "L2_compartment" %in% colnames(ng)) {
  ng <- ng %>% rename(compartment = L2_compartment)
}
if (!"label" %in% colnames(ng) && "L2_label" %in% colnames(ng)) {
  ng <- ng %>% rename(label = L2_label)
}
if (!all(c("compartment", "label") %in% colnames(ng))) {
  stop("Expected compartment + label columns from Stage D; got: ",
       paste(colnames(ng), collapse=", "))
}

ng$L2_joint <- paste(ng$compartment, ng$label, sep = "::")

# Majority-vote parent label per group
parent_per_group <- ng %>%
  filter(!is.na(NhoodGroup)) %>%
  group_by(NhoodGroup, L2_joint, compartment, label) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(NhoodGroup) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  slice_head(n = 1) %>%
  rename(parent_label = label, parent_compartment = compartment,
         parent_L2_joint = L2_joint, n_in_parent = n) %>%
  ungroup()

# Group size + LFC summary
group_summary <- ng %>%
  filter(!is.na(NhoodGroup)) %>%
  group_by(NhoodGroup) %>%
  summarise(
    n_nhoods_in_group = n(),
    group_med_lfc     = median(logFC, na.rm = TRUE),
    group_pct_up      = round(100 * sum(logFC > 0, na.rm = TRUE) / n(), 1),
    n_sig             = sum(SpatialFDR < con$spatial_fdr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(parent_per_group %>% select(NhoodGroup, parent_L2_joint,
                                        parent_compartment, parent_label,
                                        n_in_parent), by = "NhoodGroup") %>%
  mutate(NhoodGroup_renamed = paste(parent_L2_joint, NhoodGroup, sep = "_")) %>%
  arrange(parent_L2_joint, NhoodGroup)

# Final per-nhood table joined with renamed group + parent labels
ng_out <- ng %>%
  select(Nhood, NhoodGroup, logFC, SpatialFDR,
         compartment, label, L2_joint) %>%
  left_join(group_summary %>%
              select(NhoodGroup, NhoodGroup_renamed, parent_L2_joint,
                     parent_compartment, parent_label,
                     n_nhoods_in_group, group_med_lfc, group_pct_up),
            by = "NhoodGroup")

write_csv(ng_out, out_csv)
write_csv(group_summary,
          file.path(out_sub, "nhood_groups_summary.csv"))
cat(sprintf("Wrote: %s  (%d nhoods, %d groups)\n",
            out_csv, nrow(ng_out), n_distinct(ng_out$NhoodGroup, na.rm = TRUE)))
cat("=== Stage F.1 task done ===\n")
