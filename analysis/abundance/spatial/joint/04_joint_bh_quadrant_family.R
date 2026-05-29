#!/usr/bin/env Rscript
# Phase D, step 04 — joint BH adjustment across the 4-quadrant family.
#
# Reads the per-contrast da_results.csv for {C3, C3a, C3b, C3c} at a given
# scope+tier, concatenates raw p-values, applies BH across the combined pool,
# and writes a family_fdr.csv carrying the per-contrast SpatialFDR alongside
# the family-wise adjusted FDR.
#
# Rationale: the 4 within-P3 quadrant tests (UOQ/UIQ/LOQ/LIQ vs rest) form a
# single family of tests on the same Milo. Per-contrast SpatialFDR (graph-
# overlap-weighted BH) controls within-contrast FDR; family BH controls
# across-contrast FDR. We retain both. Joint adjustment is approximate (BH
# assumes independence, but the 4 tests share nhoods, so PValue correlations
# exist) — interpret family_FDR as conservative.
#
# Pattern source: scripts/da/01_run_contrast.R (CLI + YAML conventions)
#
# Usage:
#   Rscript scripts/da/04_joint_bh_quadrant_family.R \
#     --base-dir outputs/da \
#     --scope    joint_full \
#     --tier     M0 \
#     --out-dir  outputs/da/quadrant_family/joint_full

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser()
parser$add_argument("--base-dir", required = TRUE,
                    help = "Root containing per-contrast outputs (e.g., outputs/da)")
parser$add_argument("--scope", required = TRUE,
                    help = "Milo scope (e.g., joint_full, joint_Epithelial)")
parser$add_argument("--tier", required = TRUE, default = "M0",
                    help = "Model tier (M0/M1/M2). Default M0.")
parser$add_argument("--out-dir", required = TRUE,
                    help = "Output directory for family_fdr.csv")
args <- parser$parse_args()

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

CONTRASTS <- c(
  "C3_p3_uoq_vs_rest",
  "C3a_p3_uiq_vs_rest",
  "C3b_p3_loq_vs_rest",
  "C3c_p3_liq_vs_rest"
)

message("=== 04_joint_bh_quadrant_family.R ===")
message("scope: ", args$scope, " | tier: ", args$tier)

read_one <- function(contrast) {
  path <- file.path(args$base_dir, contrast, args$scope, args$tier, "da_results.csv")
  if (!file.exists(path)) {
    stop("Missing input: ", path)
  }
  dt <- fread(path)
  dt[, contrast := contrast]
  # miloR's Nhood column is the nhood index. Preserve as-is.
  setnames(dt, "SpatialFDR", "SpatialFDR_per_contrast", skip_absent = TRUE)
  dt
}

dts <- lapply(CONTRASTS, read_one)
nrows <- sapply(dts, nrow)
message("rows per contrast: ", paste0(CONTRASTS, "=", nrows, collapse = " | "))
if (length(unique(nrows)) > 1L) {
  warning("Per-contrast row counts differ — nhood index alignment is per-contrast")
}

combined <- rbindlist(dts, use.names = TRUE, fill = TRUE)
message("combined rows: ", nrow(combined), " (expected ~ ", sum(nrows), ")")

# Family BH across the pooled PValue vector
combined[, family_FDR_BH := p.adjust(PValue, method = "BH")]

# Reorder columns for readability
preferred <- c("contrast", "Nhood", "logFC", "PValue",
               "SpatialFDR_per_contrast", "family_FDR_BH",
               "nhood_label_l1p5", "nhood_label_l1p5_frac")
keep <- intersect(preferred, names(combined))
extra <- setdiff(names(combined), keep)
setcolorder(combined, c(keep, extra))

out_csv <- file.path(args$out_dir, "family_fdr.csv")
fwrite(combined, out_csv)
message("wrote ", out_csv)

# One-page summary
fam_thr <- c(0.05, 0.1)
sum_tab <- combined[, .(
  n_nhoods = .N,
  n_family_FDR_lt_0.05 = sum(family_FDR_BH < 0.05, na.rm = TRUE),
  n_family_FDR_lt_0.1  = sum(family_FDR_BH < 0.1,  na.rm = TRUE),
  n_perctrst_FDR_lt_0.05 = sum(SpatialFDR_per_contrast < 0.05, na.rm = TRUE),
  n_perctrst_FDR_lt_0.1  = sum(SpatialFDR_per_contrast < 0.1,  na.rm = TRUE)
), by = contrast]

summary_md <- file.path(args$out_dir, "quadrant_family_summary.md")
cat(
  "# Quadrant family FDR summary\n\n",
  "Scope: ", args$scope, " | Tier: ", args$tier, " | Date: ", as.character(Sys.time()), "\n\n",
  "Pooled p-values across 4 contrasts, BH-adjusted. Per-contrast graph-overlap SpatialFDR also retained.\n\n",
  "## Significant nhoods per contrast\n\n",
  "| Contrast | Total nhoods | family_FDR<0.05 | family_FDR<0.1 | per-contrast SpatialFDR<0.05 | per-contrast SpatialFDR<0.1 |\n",
  "|---|---|---|---|---|---|\n",
  sep = "", file = summary_md
)
for (i in seq_len(nrow(sum_tab))) {
  r <- sum_tab[i]
  cat("| ", r$contrast, " | ", r$n_nhoods, " | ",
      r$n_family_FDR_lt_0.05, " | ", r$n_family_FDR_lt_0.1, " | ",
      r$n_perctrst_FDR_lt_0.05, " | ", r$n_perctrst_FDR_lt_0.1, " |\n",
      sep = "", file = summary_md, append = TRUE)
}

cat("\n## Cross-contrast overlap (family_FDR<0.1, nhood × contrast)\n\n",
    sep = "", file = summary_md, append = TRUE)
# Reshape: rows = nhood, cols = contrast, value = significant or not
sig <- combined[family_FDR_BH < 0.1, .(Nhood, contrast)]
sig_wide <- dcast(sig, Nhood ~ contrast, fun.aggregate = length, value.var = "contrast")
sig_wide[, n_contrasts := rowSums(.SD > 0), .SDcols = intersect(CONTRASTS, names(sig_wide))]
overlap_dist <- table(sig_wide$n_contrasts)
cat("Distribution of n_contrasts in which a nhood is family-significant:\n\n", file = summary_md, append = TRUE)
cat("| n_contrasts | n_nhoods |\n|---|---|\n", file = summary_md, append = TRUE)
for (k in names(overlap_dist)) {
  cat("| ", k, " | ", as.integer(overlap_dist[k]), " |\n", sep = "", file = summary_md, append = TRUE)
}

message("wrote ", summary_md)
message("=== done ===")
