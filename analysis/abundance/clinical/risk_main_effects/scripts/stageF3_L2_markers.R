#!/usr/bin/env Rscript
# stageF3_L2_markers.R
# Per-(contrast x parent_L2_joint) limma-voom marker tables at L2 grain.
# Aggregates F.2 NhoodGroup pseudobulks to L2 grain (sum across NhoodGroups
# within parent_L2 per donor), then runs DE with cohort covariates.
#
# Output:
#   stageF3_L2_markers/<contrast>/<L2_safe>_vs_<contrast_var>.csv
#
# Inquiry-agnostic. Driven by inquiry.yaml's stageF3_L2 block. Supports any
# 2-level contrast variable (parity_binary, age_binary, menopausal_status_binary,
# etc.) and arbitrary base + extra covariates per contrast.
#
# Required inquiry.yaml structure:
#
#   stageF3_L2:
#     min_donors_per_group: 2          # optional, default 2
#     min_expressed_genes: 200         # optional, default 200
#     contrasts:
#       - name: parity_in_AR           # contrast name (becomes output dir)
#         cohort: AR_baseline          # which Stage A cohort_<X>_design.csv
#         f2_pseudobulk_dir: parity_in_AR  # optional; defaults to <name>
#         coef: parity_binaryparous    # exact design-matrix column to test
#         contrast_var: parity_binary  # column for the 2-level contrast
#         base_covars: [study, facs_status, age_binary, menopausal_status_binary]
#         extra_covars: []             # appended after base_covars
#         output_suffix: vs_parity     # optional; defaults to vs_<contrast_var>
#
# Dry-run mode (`--dry-run`): resolve config, validate cohort + pseudobulk
# existence, print plan. No DE computation; no output files written.

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(dplyr); library(readr); library(tidyr); library(yaml)
  library(edgeR); library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageF3_L2_markers.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

DRY_RUN <- "--dry-run" %in% commandArgs(trailingOnly = TRUE)

# -- Resolve stageF3_L2 config block ------------------------------------------
f3l2 <- inquiry$stageF3_L2 %||% list()
contrasts_cfg <- f3l2$contrasts
if (is.null(contrasts_cfg) || length(contrasts_cfg) == 0) {
  stop("inquiry.yaml is missing stageF3_L2.contrasts block")
}
MIN_DONORS_PER_GROUP <- f3l2$min_donors_per_group %||% 2
MIN_EXPRESSED_GENES  <- f3l2$min_expressed_genes  %||% 200
DEFAULT_BASE_COVARS  <- f3l2$default_base_covars  %||%
  c("study", "facs_status", "age_binary", "menopausal_status_binary")

# -- Load gene mapping (gene_id -> symbol) -------------------------------------
gm <- read_tsv(paths$inputs$gene_mapping, show_col_types = FALSE) %>%
  rename_with(tolower)
ensg_col <- intersect(c("ensembl_id", "ensembl_gene_id", "gene_id", "ensembl"),
                       colnames(gm))[1]
sym_col  <- intersect(c("gene_symbol", "symbol", "gene_name"),
                       colnames(gm))[1]
if (is.na(ensg_col) || is.na(sym_col)) {
  stop("gene_mapping schema unrecognized: ", paste(colnames(gm), collapse=","))
}
gene_id_to_symbol <- setNames(gm[[sym_col]], gm[[ensg_col]])

# -- Stage F.3 L2 main loop ---------------------------------------------------
out_root <- file.path(paths$inquiry_root, "outputs", "stageF3_L2_markers")
if (!DRY_RUN) ensure_dir(out_root)

run_one_l2 <- function(L2, l2_counts_mat, l2_meta, contrast_var, coef_name,
                       base_covars, extra_covars, factor_levels = NULL) {
  if (!is.null(factor_levels)) {
    for (col_name in names(factor_levels)) {
      if (col_name %in% colnames(l2_meta)) {
        desired <- as.character(factor_levels[[col_name]])
        current <- as.character(unique(l2_meta[[col_name]]))
        current <- current[!is.na(current)]
        if (all(current %in% desired)) {
          l2_meta[[col_name]] <- factor(l2_meta[[col_name]], levels = desired)
        }
      }
    }
  }
  level_counts <- table(l2_meta[[contrast_var]])
  if (any(level_counts < MIN_DONORS_PER_GROUP) || length(level_counts) < 2) {
    return(list(skipped = TRUE,
                reason = sprintf("contrast_var=%s level counts=%s",
                                 contrast_var,
                                 paste(names(level_counts), level_counts,
                                       sep = ":", collapse = ","))))
  }
  if (sum(rowSums(l2_counts_mat > 0) >= 3) < MIN_EXPRESSED_GENES) {
    return(list(skipped = TRUE,
                reason = sprintf("fewer than %d expressed genes",
                                 MIN_EXPRESSED_GENES)))
  }
  keep_col <- colSums(l2_counts_mat) > 0
  if (!all(keep_col)) {
    l2_counts_mat <- l2_counts_mat[, keep_col, drop = FALSE]
    l2_meta <- l2_meta[keep_col, ]
  }
  dge <- DGEList(counts = l2_counts_mat)
  dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")

  design_terms <- c(base_covars, extra_covars, contrast_var)
  design_terms <- design_terms[sapply(design_terms, function(x) {
    if (!x %in% colnames(l2_meta)) return(FALSE)
    length(unique(l2_meta[[x]])) > 1
  })]
  if (!contrast_var %in% design_terms) {
    return(list(skipped = TRUE,
                reason = sprintf("contrast_var=%s has no variation after filtering",
                                 contrast_var)))
  }
  formula_str <- paste("~", paste(design_terms, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = l2_meta)
  if (qr(design)$rank < ncol(design)) {
    return(list(skipped = TRUE, reason = "design rank deficient"))
  }

  v <- voom(dge, design)
  fit <- lmFit(v, design)
  fit <- eBayes(fit, robust = TRUE)

  if (coef_name %in% colnames(design)) {
    coef_col <- coef_name
  } else {
    matched <- grep(paste0("^", contrast_var), colnames(design), value = TRUE)
    if (length(matched) == 0) {
      return(list(skipped = TRUE,
                  reason = sprintf("coef %s not in design (cols: %s)",
                                   coef_name,
                                   paste(colnames(design), collapse = ","))))
    }
    coef_col <- tail(matched, 1)
  }

  tt <- topTable(fit, coef = coef_col, number = Inf, sort.by = "none") %>%
    as.data.frame() %>%
    mutate(gene_id = rownames(.),
           symbol = gene_id_to_symbol[gene_id],
           FDR = adj.P.Val,
           PValue = P.Value,
           parent_L2 = L2,
           contrast = l2_meta$contrast[1])
  rownames(tt) <- NULL
  list(skipped = FALSE, table = tt,
       level_counts = level_counts,
       formula = formula_str,
       coef_used = coef_col)
}

if (DRY_RUN) {
  cat("=== DRY RUN: stageF3_L2 plan ===\n\n")
} else {
  cat(sprintf("=== Stage F.3 L2 - %d contrasts to run ===\n\n",
              length(contrasts_cfg)))
}

for (c_cfg in contrasts_cfg) {
  cname        <- c_cfg$name
  cohort_name  <- c_cfg$cohort
  f2_dir       <- c_cfg$f2_pseudobulk_dir %||% cname
  coef_name    <- c_cfg$coef
  contrast_var <- c_cfg$contrast_var
  base_covars  <- c_cfg$base_covars  %||% DEFAULT_BASE_COVARS
  extra_covars <- c_cfg$extra_covars %||% character(0)
  out_suffix   <- c_cfg$output_suffix %||% paste0("vs_", contrast_var)
  factor_levels <- c_cfg$factor_levels  # NULL ok

  cat(sprintf("\n--- %s ---\n", cname))
  cat(sprintf("  cohort=%s  f2_dir=%s  coef=%s  contrast_var=%s\n",
              cohort_name, f2_dir, coef_name, contrast_var))
  cat(sprintf("  covariates: %s\n",
              paste(c(base_covars, extra_covars, contrast_var),
                    collapse = " + ")))

  pb_path <- file.path(paths$outputs$stageF2, f2_dir, "pseudobulk.rds")
  cohort_csv <- file.path(paths$outputs$stageA,
                          paste0("cohort_", cohort_name, "_design.csv"))

  if (!file.exists(pb_path)) {
    cat(sprintf("  [warn] F.2 pseudobulk missing: %s\n", pb_path))
    if (DRY_RUN) next else { cat("  skip\n"); next }
  }
  if (!file.exists(cohort_csv)) {
    cat(sprintf("  [warn] cohort CSV missing: %s\n", cohort_csv))
    if (DRY_RUN) next else { cat("  skip\n"); next }
  }
  if (DRY_RUN) {
    cat(sprintf("  pseudobulk: %s (ok)\n", pb_path))
    cat(sprintf("  cohort:     %s (ok)\n", cohort_csv))
    next
  }

  pb <- readRDS(pb_path)
  cd <- as.data.frame(colData(pb))
  cat(sprintf("  F.2 pseudobulk: %s (genes=%d, cols=%d)\n",
              class(pb)[1], nrow(pb), ncol(pb)))
  if (!"parent_L2_joint" %in% colnames(cd)) {
    cd$parent_L2_joint <- cd$L2_joint %||% cd$parent_L2 %||% NA_character_
  }
  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "patientID", "donor"),
                         colnames(cd))[1]
  cd$donor_id <- cd[[donor_col]]

  cohort <- read_csv(cohort_csv, show_col_types = FALSE)
  donor_key <- intersect(c("ihbca_donor_id", "donor_id", "patientID"),
                         colnames(cohort))[1]
  cohort$donor_id <- cohort[[donor_key]]
  if (!contrast_var %in% colnames(cohort)) {
    cat(sprintf("  cohort missing %s, skip\n", contrast_var)); next
  }
  cohort <- cohort %>% filter(!is.na(.data[[contrast_var]]))
  cat(sprintf("  cohort donors usable (%s known): %d\n",
              contrast_var, nrow(cohort)))

  L2_list <- unique(na.omit(cd$parent_L2_joint))
  cat(sprintf("  L2s to test: %d\n", length(L2_list)))
  out_dir <- file.path(out_root, cname)
  ensure_dir(out_dir)
  written <- 0; skipped <- 0
  for (L2 in L2_list) {
    L2_safe <- gsub("[/:: ]", "_", gsub("-", "_", L2))
    out_csv <- file.path(out_dir,
                         sprintf("%s_%s.csv", L2_safe, out_suffix))
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
    keep_d <- !is.na(donor_meta[[contrast_var]])
    if (sum(keep_d) < 2 * MIN_DONORS_PER_GROUP) {
      skipped <- skipped + 1; next
    }
    agg <- agg[, keep_d, drop = FALSE]
    donor_meta <- donor_meta[keep_d, ]
    res <- run_one_l2(L2, agg, donor_meta, contrast_var, coef_name,
                      base_covars, extra_covars, factor_levels)
    if (res$skipped) {
      skipped <- skipped + 1
      writeLines(c("skipped,reason", paste0("TRUE,", res$reason)), out_csv)
      next
    }
    write_csv(res$table, out_csv)
    written <- written + 1
    lc <- res$level_counts
    cat(sprintf("    [%s] levels=%s  rows=%d coef=%s\n",
                L2,
                paste(names(lc), lc, sep = ":", collapse = ","),
                nrow(res$table),
                res$coef_used))
  }
  cat(sprintf("  written=%d skipped=%d\n", written, skipped))
}

if (DRY_RUN) {
  cat("\n=== DRY RUN complete (no DE computed, no files written) ===\n")
} else {
  cat("\n=== L2 markers done ===\n")
}
