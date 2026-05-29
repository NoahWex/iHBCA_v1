#!/usr/bin/env Rscript
# Build inquiry-level long-form derivatives for V1 clinical DA inquiries.
#
# For one inquiry under publication/analysis/abundance/clinical/, this script:
#
#   1. Reads the inquiry's nhoodgroups/lookup.csv (already long over contrasts).
#   2. Computes marker_positive + n_up_markers per NG from its F.3 marker tables.
#   3. For paired_diff contrasts where a paired-diff stats table exists (currently
#      only parity_x_HR_BRCA1 via figure_data/fig2_pxr_v4_br1/marker_positive_ngs_v4_br1.csv),
#      populates int_med + int_ci_low/high + n_nh_paired + ci_excludes_zero +
#      interaction_defined.
#   4. Writes nhoodgroups/lookup.csv with the extended schema (additive).
#   5. Reads all per-NG F.3 marker CSVs under markers/nhoodgroup/ and binds them
#      into markers/nhoodgroup_long.csv (one row per contrast, NG, gene_id).
#   6. Reads all widened expression substrates under figure_data/{contrast}/ and
#      binds them into expression_long.csv (one row per contrast, NG, gene_id).
#
# CLI:
#   build_inquiry_derivatives.R --project-root <root> --inquiry {parity_findings|risk_main_effects}

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--project-root", required = TRUE)
parser$add_argument("--inquiry", required = TRUE,
                    choices = c("parity_findings", "risk_main_effects"))
args <- parser$parse_args()

root <- normalizePath(args$project_root)
inq_dir <- file.path(root, "publication", "analysis", "abundance", "clinical",
                     args$inquiry)
stopifnot(dir.exists(inq_dir))

lookup_csv  <- file.path(inq_dir, "nhoodgroups", "lookup.csv")
markers_dir <- file.path(inq_dir, "markers", "nhoodgroup")
figdata_dir <- file.path(inq_dir, "figure_data")

# Marker-positive threshold (mechanical, applied per NG from F.3)
MP_LOGFC_FLOOR <- 0.5
MP_FDR_CEILING <- 0.01
MP_MIN_COUNT   <- 5

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
f3_path_for_ng <- function(contrast, parent_L2_joint, NhoodGroup) {
  basename_ <- gsub("[-:]", "_", parent_L2_joint)
  file.path(markers_dir, contrast,
            sprintf("%s_%d_vs_parent.csv", basename_, NhoodGroup))
}

is_placeholder <- function(fp) {
  if (!file.exists(fp)) return(TRUE)
  hdr <- read_csv(fp, n_max = 0, show_col_types = FALSE)
  !all(c("logFC", "FDR", "gene_id") %in% colnames(hdr))
}

# ----------------------------------------------------------------------------
# Load lookup
# ----------------------------------------------------------------------------
lookup <- read_csv(lookup_csv, show_col_types = FALSE)
cat(sprintf("inquiry: %s\n", args$inquiry))
cat(sprintf("lookup.csv: %d rows, %d distinct contrasts\n",
            nrow(lookup), length(unique(lookup$contrast))))

# ----------------------------------------------------------------------------
# Step 1: compute marker_positive + n_up_markers per (contrast, NG)
# ----------------------------------------------------------------------------
cat("Computing marker_positive per NG...\n")
lookup$n_up_markers <- map_int(seq_len(nrow(lookup)), function(i) {
  fp <- f3_path_for_ng(lookup$contrast[i], lookup$parent_L2_joint[i],
                       lookup$NhoodGroup[i])
  if (is_placeholder(fp)) return(0L)
  df <- read_csv(fp, col_types = cols_only(logFC = col_double(), FDR = col_double()),
                 show_col_types = FALSE)
  as.integer(sum(df$logFC > MP_LOGFC_FLOOR & df$FDR < MP_FDR_CEILING, na.rm = TRUE))
})
lookup$marker_positive <- lookup$n_up_markers >= MP_MIN_COUNT
cat(sprintf("  marker+ NGs: %d of %d\n",
            sum(lookup$marker_positive), nrow(lookup)))

# ----------------------------------------------------------------------------
# Step 2: hydrate paired_diff stats for parity_x_HR_BRCA1 (only source available)
# ----------------------------------------------------------------------------
paired_diff_cols <- c("int_med", "int_ci_low", "int_ci_high",
                      "n_nh_paired", "ci_excludes_zero", "interaction_defined")
for (col in paired_diff_cols) {
  if (!col %in% colnames(lookup)) lookup[[col]] <- NA
}

mp_br1 <- file.path(figdata_dir, "fig2_pxr_v4_br1", "marker_positive_ngs_v4_br1.csv")
if (file.exists(mp_br1)) {
  cat(sprintf("Hydrating paired_diff stats from %s\n", basename(mp_br1)))
  pd <- read_csv(mp_br1, show_col_types = FALSE) |>
    transmute(
      contrast            = "parity_x_HR_BRCA1",
      NhoodGroup          = as.integer(NhoodGroup),
      int_med_new         = int_med,
      int_ci_low_new      = int_low,
      int_ci_high_new     = int_high,
      n_nh_paired_new     = as.integer(n_nh),
      ci_excludes_zero_new = as.logical(ci_excludes_zero),
      interaction_defined_new = as.logical(interaction_defined)
    )
  lookup <- lookup |>
    left_join(pd, by = c("contrast", "NhoodGroup")) |>
    mutate(
      int_med             = coalesce(int_med_new, int_med),
      int_ci_low          = coalesce(int_ci_low_new, int_ci_low),
      int_ci_high         = coalesce(int_ci_high_new, int_ci_high),
      n_nh_paired         = coalesce(n_nh_paired_new, as.integer(n_nh_paired)),
      ci_excludes_zero    = coalesce(ci_excludes_zero_new, as.logical(ci_excludes_zero)),
      interaction_defined = coalesce(interaction_defined_new, as.logical(interaction_defined))
    ) |>
    select(-ends_with("_new"))
  cat(sprintf("  hydrated %d (parity_x_HR_BRCA1) rows with paired_diff stats\n",
              sum(!is.na(lookup$int_med) & lookup$contrast == "parity_x_HR_BRCA1")))
}

# Provenance columns
lookup$da_method <- "milo_edgeR_QLF_spatialFDR_v1"
lookup$substrate_version <- "bundle_1_1_v1_clinical_da_substrate_20260510"

# ----------------------------------------------------------------------------
# Step 3: write extended lookup.csv (overwrites; additive schema)
# ----------------------------------------------------------------------------
write_csv(lookup, lookup_csv)
cat(sprintf("wrote %s (%d rows, %d cols)\n",
            lookup_csv, nrow(lookup), ncol(lookup)))

# ----------------------------------------------------------------------------
# Step 4: bind all F.3 marker files → nhoodgroup_long.csv
# ----------------------------------------------------------------------------
cat("Binding F.3 marker tables...\n")
markers_long <- map_dfr(seq_len(nrow(lookup)), function(i) {
  fp <- f3_path_for_ng(lookup$contrast[i], lookup$parent_L2_joint[i],
                       lookup$NhoodGroup[i])
  if (is_placeholder(fp)) return(NULL)
  read_csv(fp, show_col_types = FALSE) |>
    transmute(
      contrast     = lookup$contrast[i],
      NhoodGroup   = as.integer(lookup$NhoodGroup[i]),
      gene_id,
      gene_symbol  = symbol,
      logFC, AveExpr, t, P.Value, FDR, B,
      is_dissoc_stress = as.logical(is_dissoc_stress),
      marker_method = "limma_voom_within_parent_L2_F3"
    )
})
markers_out <- file.path(inq_dir, "markers", "nhoodgroup_long.csv")
write_csv(markers_long, markers_out)
cat(sprintf("wrote %s (%d rows, %d distinct (contrast, NG))\n",
            markers_out, nrow(markers_long),
            nrow(distinct(markers_long, contrast, NhoodGroup))))

# ----------------------------------------------------------------------------
# Step 5: bind widened expression substrates → expression_long.csv
# ----------------------------------------------------------------------------
cat("Binding expression substrates...\n")
contrasts_with_figdata <- unique(lookup$contrast)
expr_rows <- list()
for (cn in contrasts_with_figdata) {
  for (subdir in c(cn, "fig2_pxr_v4_br1")) {
    expr_dir_candidate <- file.path(figdata_dir, subdir)
    if (!dir.exists(expr_dir_candidate)) next
    for (comp in c("epi", "imm", "str")) {
      fp <- file.path(expr_dir_candidate,
                      sprintf("%s_marker_expression_full.csv", comp))
      if (!file.exists(fp)) next
      e <- read_csv(fp, show_col_types = FALSE) |>
        transmute(
          contrast    = if (subdir == "fig2_pxr_v4_br1") "parity_x_HR_BRCA1" else cn,
          NhoodGroup  = as.integer(NhoodGroup),
          gene_id     = gene_ensembl,
          gene_symbol,
          mean_log_cpm,
          pct_expressing,
          n_cells_in_NG = as.integer(n_cells),
          aggregation_method = "pct_mean_aggregate_v4"
        )
      expr_rows[[paste(subdir, comp, sep = "_")]] <- e
      cat(sprintf("  [%s/%s] %d rows\n", subdir, comp, nrow(e)))
    }
  }
}
expr_long <- bind_rows(expr_rows) |> distinct()
if (nrow(expr_long) == 0) {
  expr_long <- tibble(
    contrast = character(), NhoodGroup = integer(),
    gene_id = character(), gene_symbol = character(),
    mean_log_cpm = double(), pct_expressing = double(),
    n_cells_in_NG = integer(), aggregation_method = character()
  )
}
expr_out <- file.path(inq_dir, "expression_long.csv")
write_csv(expr_long, expr_out)
cat(sprintf("wrote %s (%d rows, %d distinct (contrast, NG))\n",
            expr_out, nrow(expr_long),
            if (nrow(expr_long) > 0) nrow(distinct(expr_long, contrast, NhoodGroup)) else 0))

cat("\n=== Done. Per-inquiry derivatives written:\n")
cat(sprintf("  %s/nhoodgroups/lookup.csv (extended)\n", args$inquiry))
cat(sprintf("  %s/markers/nhoodgroup_long.csv\n", args$inquiry))
cat(sprintf("  %s/expression_long.csv\n", args$inquiry))
