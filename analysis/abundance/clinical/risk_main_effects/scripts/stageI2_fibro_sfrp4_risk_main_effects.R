#!/usr/bin/env Rscript
# stageI2_fibro_sfrp4_risk_main_effects.R
#
# Three parallel main-effect DE tests for str::Fibro-SFRP4:
#   AR vs BRCA1, AR vs sporadic-HR, AR vs BRCA2.
#
# Each test pools the AR cohort with one risk cohort, includes parity_binary as
# a covariate (so risk effects are separated from parity effects), and tests
# the cohort coefficient.
#
# Design: ~ study + facs + age + menop [+ cancer_history] + parity_binary + cohort
# (covariates dropped if no within-pool variation). AR is reference.
#
# Companion to stageI2_fibro_sfrp4_BR1vsAR_parous.R, which is the parous-only
# stratified version of the AR-vs-BR1 test.
#
# Outputs:
#   outputs/stageI2_within_l2_de/fibro_sfrp4_BR1_vs_AR_main.csv
#   outputs/stageI2_within_l2_de/fibro_sfrp4_HRS_vs_AR_main.csv
#   outputs/stageI2_within_l2_de/fibro_sfrp4_BR2_vs_AR_main.csv
#   outputs/stageI2_within_l2_de/fibro_sfrp4_risk_main_effects_design.tsv

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(dplyr); library(readr); library(tidyr); library(yaml)
  library(edgeR); library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI2_fibro_sfrp4_risk_main_effects.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

L2_TARGET <- "str::Fibro-SFRP4"
BASE_COVARS  <- c("study", "facs_status", "age_binary", "menopausal_status_binary")
EXTRA_COVARS <- c("cancer_history", "parity_binary")
DROP_PRIORITY <- c("cancer_history", "facs_status",
                    "menopausal_status_binary", "age_binary")
MIN_DONORS_PER_GROUP <- 2

# Map contrast -> short cohort label and stageA cohort design filename
PB_BY_COHORT <- list(
  AR  = list(contrast = "parity_in_AR",         design_csv = "AR_baseline"),
  BR1 = list(contrast = "parity_x_HR_BRCA1",    design_csv = "HR_BRCA1"),
  HRS = list(contrast = "parity_x_HR_sporadic", design_csv = "HR_sporadic"),
  BR2 = list(contrast = "parity_x_HR_BRCA2",    design_csv = "HR_BRCA2")
)

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
if (is.na(ens_col) || is.na(sym_col))
  stop(sprintf("gene mapping at %s has unexpected schema: %s",
                gm_path, paste(colnames(gm), collapse = ",")))
gene_id_to_symbol <- setNames(gm[[sym_col]], gm[[ens_col]])
cat(sprintf("Loaded gene mapping: %s (%d entries)\n", gm_path, length(gene_id_to_symbol)))

# ---- per-cohort donor-level aggregation ----
agg_donor_pb <- function(cname) {
  pb_path <- file.path(paths$outputs$stageF2, cname, "pseudobulk.rds")
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

read_design <- function(cohort_csv) {
  fp <- file.path(paths$outputs$stageA,
                   paste0("cohort_", cohort_csv, "_design.csv"))
  d <- read_csv(fp, show_col_types = FALSE)
  donor_key <- intersect(c("ihbca_donor_id", "donor_id", "patientID"),
                          colnames(d))[1]
  d$donor_id <- d[[donor_key]]
  d
}

# Pre-load all four cohorts' pseudobulks + designs
loaded <- list()
for (lbl in names(PB_BY_COHORT)) {
  spec <- PB_BY_COHORT[[lbl]]
  pb <- agg_donor_pb(spec$contrast)
  if (is.null(pb)) {
    cat(sprintf("[%s] no pseudobulk available, skip.\n", lbl)); next
  }
  ds <- read_design(spec$design_csv)
  loaded[[lbl]] <- list(pb = pb, design = ds)
  cat(sprintf("[%s] pseudobulk donors=%d, design rows=%d\n",
              lbl, length(pb$donors), nrow(ds)))
}

# ---- run one (AR vs <risk>) main-effect test ----
run_main_effect <- function(risk_label) {
  if (!"AR" %in% names(loaded) || !risk_label %in% names(loaded)) {
    cat(sprintf("Cannot run AR vs %s: missing substrate\n", risk_label))
    return(NULL)
  }
  ar  <- loaded[["AR"]]
  rsk <- loaded[[risk_label]]

  keep_ar  <- intersect(colnames(ar$pb$counts),  ar$design$donor_id)
  keep_rsk <- intersect(colnames(rsk$pb$counts), rsk$design$donor_id)

  cat(sprintf("\n--- AR vs %s ---\n", risk_label))
  cat(sprintf("  AR donors with pseudobulks: %d\n", length(keep_ar)))
  cat(sprintf("  %s donors with pseudobulks: %d\n", risk_label, length(keep_rsk)))

  if (length(keep_ar) < MIN_DONORS_PER_GROUP ||
      length(keep_rsk) < MIN_DONORS_PER_GROUP) {
    cat("  too few donors, skip\n")
    return(NULL)
  }

  ar_counts  <- ar$pb$counts[,  keep_ar,  drop = FALSE]
  rsk_counts <- rsk$pb$counts[, keep_rsk, drop = FALSE]
  stopifnot(identical(rownames(ar_counts), rownames(rsk_counts)))
  counts_mat <- cbind(ar_counts, rsk_counts)

  ar_meta  <- tibble(donor_id = keep_ar,  cohort = "AR")  %>%
    left_join(ar$design,  by = "donor_id")
  rsk_meta <- tibble(donor_id = keep_rsk, cohort = risk_label) %>%
    left_join(rsk$design, by = "donor_id")
  meta <- bind_rows(ar_meta, rsk_meta)
  stopifnot(identical(meta$donor_id, colnames(counts_mat)))
  meta$cohort <- factor(meta$cohort, levels = c("AR", risk_label))

  # NA -> "unknown" in covariates so model.matrix doesn't drop rows
  covar_cols <- intersect(c(BASE_COVARS, EXTRA_COVARS), colnames(meta))
  for (cv in covar_cols) {
    v <- meta[[cv]]
    if (is.character(v) || is.factor(v)) {
      v <- as.character(v); v[is.na(v) | v == ""] <- "unknown"
      meta[[cv]] <- v
    }
  }

  cat("  sample distribution:\n")
  print(meta %>% count(cohort, parity_binary))

  dge <- DGEList(counts = counts_mat)
  dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  candidate_terms <- c(BASE_COVARS, EXTRA_COVARS, "cohort")
  use_terms <- candidate_terms[sapply(candidate_terms, function(x) {
    if (!x %in% colnames(meta)) return(FALSE)
    length(unique(meta[[x]])) > 1
  })]
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
    cat(sprintf("  rank deficient with [%s] — simplifying.\n", bd$fstr))
    to_drop <- intersect(DROP_PRIORITY, use_terms)
    if (length(to_drop) == 0) break
    use_terms <- setdiff(use_terms, to_drop[1])
  }
  if (is.null(design)) {
    cat("  could not build full-rank design even after simplification.\n")
    return(NULL)
  }
  cat(sprintf("  formula: %s\n", formula_str))

  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)

  coef_pat <- paste0("^cohort", risk_label, "$")
  coef_name <- grep(coef_pat, colnames(design), value = TRUE)
  if (length(coef_name) != 1) {
    cat(sprintf("  cohort coef not found via %s in: %s\n",
                coef_pat, paste(colnames(design), collapse = ",")))
    return(NULL)
  }

  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
    as.data.frame() %>%
    mutate(gene_id = rownames(.),
            symbol = gene_id_to_symbol[gene_id],
            FDR = adj.P.Val,
            PValue = P.Value,
            parent_L2 = L2_TARGET,
            comparison = sprintf("%s_vs_AR_main", risk_label))
  rownames(tt) <- NULL
  list(table = tt,
        info = tibble(comparison = sprintf("%s_vs_AR_main", risk_label),
                       L2 = L2_TARGET,
                       n_AR_donors = length(keep_ar),
                       n_risk_donors = length(keep_rsk),
                       formula = formula_str,
                       coef_tested = coef_name,
                       n_genes_tested = nrow(tt),
                       n_sig_FDR05 = sum(tt$FDR < 0.05, na.rm = TRUE)))
}

out_dir <- file.path(paths$inquiry_root, "outputs", "stageI2_within_l2_de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

all_info <- list()
for (lbl in c("BR1", "HRS", "BR2")) {
  if (lbl == "BR2") {
    cat("\n[BR2] Skipping main-effect DE: BR2 cohort has too few donors with Fibro-SFRP4 pseudobulks (n=2). Underpowered.\n")
    next
  }
  res <- run_main_effect(lbl)
  if (is.null(res)) next
  out_csv <- file.path(out_dir, sprintf("fibro_sfrp4_%s_vs_AR_main.csv", lbl))
  write_csv(res$table, out_csv)
  cat(sprintf("  Wrote: %s (rows=%d, sig FDR<0.05: %d)\n",
              out_csv, nrow(res$table),
              sum(res$table$FDR < 0.05, na.rm = TRUE)))
  all_info[[lbl]] <- res$info

  # quick sanity prints
  cat(sprintf("\n--- top 12 UP in %s vs AR (main, FDR<0.05) ---\n", lbl))
  print(res$table %>% filter(FDR < 0.05, logFC > 0) %>%
          arrange(desc(t)) %>% head(12) %>%
          select(symbol, gene_id, logFC, AveExpr, FDR))
  cat(sprintf("--- top 12 DOWN in %s vs AR (main, FDR<0.05) ---\n", lbl))
  print(res$table %>% filter(FDR < 0.05, logFC < 0) %>%
          arrange(t) %>% head(12) %>%
          select(symbol, gene_id, logFC, AveExpr, FDR))
}

if (length(all_info) > 0) {
  info_df <- bind_rows(all_info)
  write_tsv(info_df,
            file.path(out_dir, "fibro_sfrp4_risk_main_effects_design.tsv"))
}

cat("\n=== stageI2_fibro_sfrp4_risk_main_effects done ===\n")
