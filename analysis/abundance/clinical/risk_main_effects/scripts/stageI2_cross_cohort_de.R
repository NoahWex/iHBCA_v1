#!/usr/bin/env Rscript
# stageI2_cross_cohort_de.R
#
# Parametric within-L2 cross-cohort limma-voom DE on donor-level pseudobulks.
# Replaces the hand-written Fibro-SFRP4 specific scripts. Single code path
# serves any (parent_L2 × comparison × donor_filter) combination.
#
# Usage:
#   Rscript stageI2_cross_cohort_de.R \
#     --inquiry-dir <inq> \
#     --L2 "str::Fibro-SFRP4" \
#     --comparison "BR1vsAR" \
#     --donor-filter "parous"        # or "nulliparous" or "any"
#
# Comparison codes: "BR1vsAR", "HRSvsAR", "BR2vsAR", "BR1vsHRS"
#                   (anything matching <RISK>vs<REF>; reference goes second)
# donor-filter: parity_binary filter applied before the test.
#                  - "parous"      → only parous donors
#                  - "nulliparous" → only nullip donors
#                  - "any"         → both, with parity_binary added as covariate
#
# Outputs:
#   outputs/stageI2_within_l2_de/<L2_safe>/<comparison>_<filter>.csv
#   outputs/stageI2_within_l2_de/<L2_safe>/<comparison>_<filter>_design.tsv
#   outputs/stageI2_within_l2_de/INDEX.tsv  (appended; one row per run)

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(dplyr); library(readr); library(tidyr); library(yaml); library(argparse)
  library(edgeR); library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI2_cross_cohort_de.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

# ---- Parse args ----
parser <- ArgumentParser()
parser$add_argument("--inquiry-dir", required = TRUE)
parser$add_argument("--L2", required = TRUE,
                     help = "parent_L2_joint, e.g. str::Fibro-SFRP4")
parser$add_argument("--comparison", required = TRUE,
                     help = "<RISK>vs<REF>, e.g. BR1vsAR, HRSvsAR, BR2vsAR, BR1vsHRS")
parser$add_argument("--donor-filter", default = "any",
                     choices = c("parous", "nulliparous", "any"))
args <- parser$parse_args()

cfg <- list(paths = parse_args_inquiry(c("--inquiry-dir", args$inquiry_dir))$paths)
paths <- cfg$paths

L2_TARGET <- args$L2
DONOR_FILTER <- args$donor_filter

# Parse comparison: "<RISK>vs<REF>"
m <- regmatches(args$comparison,
                 regexec("^([A-Za-z0-9]+)vs([A-Za-z0-9]+)$", args$comparison))[[1]]
if (length(m) != 3) stop(sprintf("Bad comparison format: %s", args$comparison))
RISK_LBL <- m[2]; REF_LBL <- m[3]

# ---- cohort label -> contrast/design lookups ----
COHORT_REGISTRY <- list(
  AR  = list(contrast = "parity_in_AR",         design_csv = "AR_baseline"),
  BR1 = list(contrast = "parity_x_HR_BRCA1",    design_csv = "HR_BRCA1"),
  HRS = list(contrast = "parity_x_HR_sporadic", design_csv = "HR_sporadic"),
  BR2 = list(contrast = "parity_x_HR_BRCA2",    design_csv = "HR_BRCA2")
)
if (!RISK_LBL %in% names(COHORT_REGISTRY) || !REF_LBL %in% names(COHORT_REGISTRY))
  stop(sprintf("Unknown cohort label(s): %s, %s", RISK_LBL, REF_LBL))

# ---- Settings ----
BASE_COVARS  <- c("study", "facs_status", "age_binary", "menopausal_status_binary")
EXTRA_COVARS <- c("cancer_history")  # parity_binary added below if needed
DROP_PRIORITY <- c("cancer_history", "facs_status",
                    "menopausal_status_binary", "age_binary")
MIN_DONORS_PER_GROUP <- 2

# ---- gene_id -> symbol ----
gm_path_share <- file.path(paths$inquiry_root, "share", "kai_v1", "06_markers",
                            "gene_symbol_to_ensembl.tsv")
gm_path_default <- paths$inputs$gene_mapping
gm_path <- if (file.exists(gm_path_share)) gm_path_share else gm_path_default
gm <- read_tsv(gm_path, show_col_types = FALSE) %>% rename_with(tolower)
ens_col <- intersect(c("ensembl_id", "ensembl", "ensembl_gene_id", "gene_id"),
                      colnames(gm))[1]
sym_col <- intersect(c("gene_symbol", "symbol", "gene_name"),
                      colnames(gm))[1]
gene_id_to_symbol <- setNames(gm[[sym_col]], gm[[ens_col]])

# ---- substrate aggregation ----
agg_donor_pb <- function(cohort_lbl) {
  spec <- COHORT_REGISTRY[[cohort_lbl]]
  pb_path <- file.path(paths$outputs$stageF2, spec$contrast, "pseudobulk.rds")
  if (!file.exists(pb_path)) return(NULL)
  pb <- readRDS(pb_path)
  cd <- as.data.frame(colData(pb))
  if (!"parent_L2_joint" %in% colnames(cd)) {
    cd$parent_L2_joint <- cd$L2_joint %||% cd$parent_L2 %||% NA_character_
  }
  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "patientID", "donor"),
                          colnames(cd))[1]
  cd$donor_id <- cd[[donor_col]]
  cols_in_L2 <- which(cd$parent_L2_joint == L2_TARGET)
  if (length(cols_in_L2) < 2) return(NULL)
  sub <- assay(pb, "counts")[, cols_in_L2, drop = FALSE]
  sub_cd <- cd[cols_in_L2, ]
  donors <- unique(sub_cd$donor_id)
  agg <- matrix(0, nrow = nrow(sub), ncol = length(donors),
                 dimnames = list(rownames(sub), donors))
  for (d in donors) {
    cd_d <- which(sub_cd$donor_id == d)
    agg[, d] <- if (length(cd_d) == 1) as.numeric(sub[, cd_d]) else
                  as.numeric(rowSums(sub[, cd_d, drop = FALSE]))
  }
  list(counts = agg, donors = donors)
}

read_design <- function(cohort_lbl) {
  spec <- COHORT_REGISTRY[[cohort_lbl]]
  fp <- file.path(paths$outputs$stageA,
                   paste0("cohort_", spec$design_csv, "_design.csv"))
  d <- read_csv(fp, show_col_types = FALSE)
  donor_key <- intersect(c("ihbca_donor_id", "donor_id", "patientID"),
                          colnames(d))[1]
  d$donor_id <- d[[donor_key]]
  d
}

cat(sprintf("\n=== %s : %s vs %s (filter=%s) ===\n",
            L2_TARGET, RISK_LBL, REF_LBL, DONOR_FILTER))

risk_pb  <- agg_donor_pb(RISK_LBL)
ref_pb   <- agg_donor_pb(REF_LBL)
if (is.null(risk_pb) || is.null(ref_pb))
  stop("Missing pseudobulk substrate for one or both cohorts.")

risk_design <- read_design(RISK_LBL)
ref_design  <- read_design(REF_LBL)

if (DONOR_FILTER %in% c("parous", "nulliparous")) {
  risk_design <- risk_design %>% filter(parity_binary == DONOR_FILTER)
  ref_design  <- ref_design  %>% filter(parity_binary == DONOR_FILTER)
}

keep_risk <- intersect(colnames(risk_pb$counts), risk_design$donor_id)
keep_ref  <- intersect(colnames(ref_pb$counts),  ref_design$donor_id)
cat(sprintf("  %s donors with pseudobulks (post filter): %d\n", RISK_LBL, length(keep_risk)))
cat(sprintf("  %s donors with pseudobulks (post filter): %d\n", REF_LBL,  length(keep_ref)))

if (length(keep_risk) < MIN_DONORS_PER_GROUP || length(keep_ref) < MIN_DONORS_PER_GROUP)
  stop(sprintf("Too few donors (need >=%d each).", MIN_DONORS_PER_GROUP))

risk_counts <- risk_pb$counts[, keep_risk, drop = FALSE]
ref_counts  <- ref_pb$counts[,  keep_ref,  drop = FALSE]
stopifnot(identical(rownames(risk_counts), rownames(ref_counts)))
counts_mat <- cbind(ref_counts, risk_counts)

ref_meta  <- tibble(donor_id = keep_ref,  cohort = REF_LBL)  %>%
  left_join(ref_design,  by = "donor_id")
risk_meta <- tibble(donor_id = keep_risk, cohort = RISK_LBL) %>%
  left_join(risk_design, by = "donor_id")
meta <- bind_rows(ref_meta, risk_meta)
stopifnot(identical(meta$donor_id, colnames(counts_mat)))
meta$cohort <- factor(meta$cohort, levels = c(REF_LBL, RISK_LBL))

# NA -> "unknown" in covariates
covar_cols <- intersect(c(BASE_COVARS, EXTRA_COVARS, "parity_binary"),
                         colnames(meta))
for (cv in covar_cols) {
  v <- meta[[cv]]
  if (is.character(v) || is.factor(v)) {
    v <- as.character(v); v[is.na(v) | v == ""] <- "unknown"
    meta[[cv]] <- v
  }
}

cat("\nSample distribution (cohort × parity):\n")
print(meta %>% count(cohort, parity_binary))

# ---- Build candidate covariate set ----
# Include parity_binary as covariate only when filter == "any"
candidate_terms <- c(BASE_COVARS, EXTRA_COVARS,
                      if (DONOR_FILTER == "any") "parity_binary" else character(0),
                      "cohort")
use_terms <- candidate_terms[sapply(candidate_terms, function(x) {
  if (!x %in% colnames(meta)) return(FALSE)
  length(unique(meta[[x]])) > 1
})]

# ---- limma-voom with progressive simplification ----
dge <- DGEList(counts = counts_mat)
dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

build_design <- function(terms, dat) {
  fstr <- paste("~", paste(terms, collapse = " + "))
  list(fstr = fstr, design = model.matrix(as.formula(fstr), data = dat))
}
design <- NULL; formula_str <- NULL
for (attempt in seq_len(length(DROP_PRIORITY) + 1)) {
  bd <- build_design(use_terms, meta)
  if (qr(bd$design)$rank == ncol(bd$design)) {
    design <- bd$design; formula_str <- bd$fstr; break
  }
  cat(sprintf("Rank deficient with [%s] — simplifying.\n", bd$fstr))
  to_drop <- intersect(DROP_PRIORITY, use_terms)
  if (length(to_drop) == 0) break
  use_terms <- setdiff(use_terms, to_drop[1])
}
if (is.null(design)) stop("Could not build full-rank design.")
cat(sprintf("Final formula: %s\n", formula_str))

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)
coef_pat <- paste0("^cohort", RISK_LBL, "$")
coef_name <- grep(coef_pat, colnames(design), value = TRUE)
if (length(coef_name) != 1)
  stop(sprintf("Cohort coef not unique: %s", paste(colnames(design), collapse = ",")))

tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(.),
          symbol = gene_id_to_symbol[gene_id],
          FDR = adj.P.Val,
          PValue = P.Value,
          parent_L2 = L2_TARGET,
          comparison = sprintf("%s_vs_%s_%s", RISK_LBL, REF_LBL, DONOR_FILTER))
rownames(tt) <- NULL

# ---- Outputs ----
L2_safe <- gsub("[/:: ]", "_", gsub("-", "_", L2_TARGET))
out_dir <- file.path(paths$inquiry_root, "outputs",
                      "stageI2_within_l2_de", L2_safe)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir,
                      sprintf("%s_vs_%s_%s.csv", RISK_LBL, REF_LBL, DONOR_FILTER))
out_design <- file.path(out_dir,
                         sprintf("%s_vs_%s_%s_design.tsv",
                                  RISK_LBL, REF_LBL, DONOR_FILTER))
write_csv(tt, out_csv)
cat(sprintf("\nWrote: %s\n  rows=%d  sig FDR<0.05: %d\n",
            out_csv, nrow(tt), sum(tt$FDR < 0.05, na.rm = TRUE)))

# Sidecar
info <- tibble(
  L2 = L2_TARGET,
  comparison = sprintf("%s_vs_%s_%s", RISK_LBL, REF_LBL, DONOR_FILTER),
  risk_cohort = RISK_LBL, ref_cohort = REF_LBL,
  donor_filter = DONOR_FILTER,
  n_risk = length(keep_risk), n_ref = length(keep_ref),
  formula = formula_str,
  coef_tested = coef_name,
  n_genes_tested = nrow(tt),
  n_sig_FDR05_up   = sum(tt$FDR < 0.05 & tt$logFC > 0, na.rm = TRUE),
  n_sig_FDR05_down = sum(tt$FDR < 0.05 & tt$logFC < 0, na.rm = TRUE),
  run_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)
write_tsv(info, out_design)

# Append to top-level INDEX.tsv (one row per run)
idx_path <- file.path(paths$inquiry_root, "outputs",
                       "stageI2_within_l2_de", "INDEX.tsv")
if (file.exists(idx_path)) {
  idx_existing <- read_tsv(idx_path, show_col_types = FALSE)
  # de-duplicate same (L2 × comparison) by overwriting the prior row
  key <- paste(info$L2, info$comparison)
  prior_keys <- paste(idx_existing$L2, idx_existing$comparison)
  idx_existing <- idx_existing[prior_keys != key, ]
  idx_new <- bind_rows(idx_existing, info)
} else {
  idx_new <- info
}
write_tsv(idx_new, idx_path)
cat(sprintf("Index updated: %s (%d rows)\n", idx_path, nrow(idx_new)))

# Sanity prints
cat("\n=== Top 15 UP (cohort effect, FDR<0.05) ===\n")
print(tt %>% filter(FDR < 0.05, logFC > 0) %>%
        arrange(desc(t)) %>% head(15) %>%
        select(symbol, gene_id, logFC, AveExpr, FDR))
cat("\n=== Top 15 DOWN (cohort effect, FDR<0.05) ===\n")
print(tt %>% filter(FDR < 0.05, logFC < 0) %>%
        arrange(t) %>% head(15) %>%
        select(symbol, gene_id, logFC, AveExpr, FDR))

cat(sprintf("\n=== %s | %s vs %s (%s) done ===\n",
            L2_TARGET, RISK_LBL, REF_LBL, DONOR_FILTER))
