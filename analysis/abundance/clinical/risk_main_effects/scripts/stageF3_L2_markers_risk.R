#!/usr/bin/env Rscript
# stageF3_L2_markers_risk.R
# Risk-side companion to stageF3_L2_markers.R (parity-side).
# Per-(contrast x parent_L2_joint) limma-voom on the risk_class coefficient
# (BR1, BR2, or HRS vs AR), pseudobulks aggregated to L2 grain by summing
# across NhoodGroups within parent_L2 per donor.
#
# Output:
#   stageF3_L2_markers/<contrast>/<L2_safe>_markers.csv  per-L2 marker table
#
# Scope: BR1_vs_AR + BR2_vs_AR + HRS_vs_AR (HRS_no_cancerhx excluded —
# definitional collinearity with HR classification, see inquiry.yaml).

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(dplyr); library(readr); library(tidyr); library(yaml)
  library(edgeR); library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageF3_L2_markers_risk.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

# Per-contrast: (cohort_name, coef_name, test_class)
# cohort_name -> cohort_<name>_design.csv
# coef_name   -> column to test in the model matrix
# test_class  -> the non-AR risk_class level (used for min-donor checks)
RISK_CONTRASTS <- list(
  BR1_vs_AR = list(cohort = "BR1_vs_AR", coef = "risk_classBR1", test_class = "BR1"),
  BR2_vs_AR = list(cohort = "BR2_vs_AR", coef = "risk_classBR2", test_class = "BR2"),
  HRS_vs_AR = list(cohort = "HRS_vs_AR", coef = "risk_classHRS", test_class = "HRS")
)

# Base covariates blocked in every risk main-effect formula (parity adjusted).
BASE_COVARS <- c("study", "facs_status", "age_binary", "menopausal_status_binary", "parity_binary")

MIN_DONORS_PER_GROUP <- 2   # >=2 AR + >=2 (BR1|BR2|HRS) per L2 to fit

# ------------------------------------------------------------------------
gm <- read_tsv(paths$inputs$gene_mapping, show_col_types = FALSE) %>%
  rename_with(tolower)
ensg_col <- intersect(c("ensembl_id", "ensembl_gene_id", "gene_id", "ensembl"), colnames(gm))[1]
sym_col  <- intersect(c("gene_symbol", "symbol", "gene_name"), colnames(gm))[1]
if (is.na(ensg_col) || is.na(sym_col)) stop("gene_mapping schema unrecognized: ", paste(colnames(gm), collapse=","))
gene_id_to_symbol <- setNames(gm[[sym_col]], gm[[ensg_col]])

out_root <- file.path(paths$inquiry_root, "outputs", "stageF3_L2_markers")
ensure_dir(out_root)

run_one_l2 <- function(L2, l2_counts_mat, l2_meta, coef_target) {
  # l2_meta must have risk_class with both AR and test class populated
  rc_tab <- table(l2_meta$risk_class)
  if (length(rc_tab) < 2 || any(rc_tab < MIN_DONORS_PER_GROUP)) {
    return(list(skipped = TRUE,
                reason = sprintf("risk_class counts: %s",
                                  paste(names(rc_tab), rc_tab, sep="=", collapse=","))))
  }
  if (sum(rowSums(l2_counts_mat > 0) >= 3) < 200) {
    return(list(skipped = TRUE, reason = "fewer than 200 expressed genes"))
  }
  keep_col <- colSums(l2_counts_mat) > 0
  if (!all(keep_col)) {
    l2_counts_mat <- l2_counts_mat[, keep_col, drop = FALSE]
    l2_meta <- l2_meta[keep_col, ]
  }
  dge <- DGEList(counts = l2_counts_mat)
  dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  # Build design — drop terms with no variation in this subset
  design_terms <- c(BASE_COVARS, "risk_class")
  design_terms <- design_terms[sapply(design_terms, function(x) {
    if (!x %in% colnames(l2_meta)) return(FALSE)
    length(unique(l2_meta[[x]])) > 1
  })]
  if (!"risk_class" %in% design_terms) {
    return(list(skipped = TRUE, reason = "risk_class has no variation"))
  }
  formula_str <- paste("~", paste(design_terms, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = l2_meta)
  if (qr(design)$rank < ncol(design)) {
    return(list(skipped = TRUE, reason = "design rank deficient"))
  }
  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)
  cn <- grep(coef_target, colnames(design), value = TRUE, fixed = TRUE)
  if (length(cn) != 1) {
    return(list(skipped = TRUE,
                reason = sprintf("coef %s not unique in design: %s",
                                  coef_target,
                                  paste(colnames(design), collapse=","))))
  }
  tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none") %>%
    as.data.frame() %>%
    mutate(gene_id = rownames(.),
           symbol = gene_id_to_symbol[gene_id],
           FDR = adj.P.Val,
           PValue = P.Value,
           parent_L2 = L2,
           contrast = l2_meta$contrast[1])
  rownames(tt) <- NULL
  list(skipped = FALSE, table = tt,
       n_AR = sum(l2_meta$risk_class == "AR"),
       n_test = sum(l2_meta$risk_class != "AR"),
       formula = formula_str)
}

for (cname in names(RISK_CONTRASTS)) {
  spec <- RISK_CONTRASTS[[cname]]
  cat(sprintf("\n=== %s (coef=%s) ===\n", cname, spec$coef))
  pb_path <- file.path(paths$outputs$stageF2, cname, "pseudobulk.rds")
  if (!file.exists(pb_path)) {
    cat("  no F.2 pseudobulk, skip\n"); next
  }
  pb <- readRDS(pb_path)
  cat(sprintf("  F.2 pseudobulk: %s (genes=%d, cols=%d)\n",
              class(pb)[1], nrow(pb), ncol(pb)))
  cd <- as.data.frame(colData(pb))
  if (!"parent_L2_joint" %in% colnames(cd)) {
    cd$parent_L2_joint <- cd$L2_joint %||% cd$parent_L2 %||% NA_character_
  }
  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "patientID", "donor"),
                          colnames(cd))[1]
  cd$donor_id <- cd[[donor_col]]

  cohort_csv <- file.path(paths$outputs$stageA,
                          paste0("cohort_", spec$cohort, "_design.csv"))
  if (!file.exists(cohort_csv)) {
    cat(sprintf("  cohort CSV missing: %s, skip\n", cohort_csv)); next
  }
  cohort <- read_csv(cohort_csv, show_col_types = FALSE)
  donor_key <- intersect(c("ihbca_donor_id", "donor_id", "patientID"),
                          colnames(cohort))[1]
  cohort$donor_id <- cohort[[donor_key]]
  if (!"risk_class" %in% colnames(cohort)) {
    cat("  cohort missing risk_class, skip\n"); next
  }
  cohort <- cohort %>% filter(risk_class %in% c("AR", spec$test_class))
  cat(sprintf("  cohort donors usable: %d (AR=%d, %s=%d)\n",
              nrow(cohort),
              sum(cohort$risk_class == "AR"),
              spec$test_class,
              sum(cohort$risk_class == spec$test_class)))

  L2_list <- unique(na.omit(cd$parent_L2_joint))
  cat(sprintf("  L2s to test: %d\n", length(L2_list)))
  out_dir <- file.path(out_root, cname)
  ensure_dir(out_dir)
  written <- 0; skipped <- 0
  for (L2 in L2_list) {
    L2_safe <- gsub("[/:: ]", "_", gsub("-", "_", L2))
    out_csv <- file.path(out_dir, paste0(L2_safe, "_markers.csv"))
    if (file.exists(out_csv) && file.size(out_csv) > 1000) next
    cols_in_L2 <- which(cd$parent_L2_joint == L2)
    if (length(cols_in_L2) < 2) next
    sub <- assay(pb, "counts")[, cols_in_L2, drop = FALSE]
    sub_meta <- cd[cols_in_L2, ]
    donors <- unique(sub_meta$donor_id)
    agg <- matrix(0, nrow = nrow(sub), ncol = length(donors),
                  dimnames = list(rownames(sub), donors))
    for (d in donors) {
      cols_d <- which(sub_meta$donor_id == d)
      if (length(cols_d) == 1) {
        agg[, d] <- as.numeric(sub[, cols_d])
      } else {
        agg[, d] <- as.numeric(rowSums(sub[, cols_d]))
      }
    }
    donor_meta <- tibble(donor_id = donors) %>%
      left_join(cohort, by = "donor_id") %>%
      mutate(contrast = cname)
    keep_d <- !is.na(donor_meta$risk_class)
    if (sum(keep_d) < 2 * MIN_DONORS_PER_GROUP) {
      skipped <- skipped + 1; next
    }
    agg <- agg[, keep_d, drop = FALSE]
    donor_meta <- donor_meta[keep_d, ]
    res <- run_one_l2(L2, agg, donor_meta, spec$coef)
    if (res$skipped) {
      skipped <- skipped + 1
      writeLines(c("skipped,reason", paste0("TRUE,", res$reason)), out_csv)
      next
    }
    write_csv(res$table, out_csv)
    written <- written + 1
    cat(sprintf("    [%s] n_AR=%d n_%s=%d  rows=%d\n",
                L2, res$n_AR, spec$test_class, res$n_test, nrow(res$table)))
  }
  cat(sprintf("  written=%d skipped=%d\n", written, skipped))
}
cat("\n=== stageF3_L2_markers_risk DONE ===\n")
