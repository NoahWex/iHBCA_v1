#!/usr/bin/env Rscript
# 02_voom_de_sweep.R — Per-L2 limma voom DE across the 8-design formula ladder.
#
# For each (cohort x formula x L2_type):
#   1. Subset pseudobulk to cohort donors x this L2
#   2. Apply donor-coverage gate (>= min_donors per class with >= min_cells)
#   3. Build design matrix; audit rank (skip if singular)
#   4. filterByExpr -> voom -> lmFit -> eBayes
#   5. topTable on risk_classBR1 coefficient -> save DE CSV
#   6. Append to design_audit.csv: testable, rank, n_sig, etc.
#
# Pattern sources (verified, code-stewardship provenance):
#   - L2-grain voom DE on risk_class: risk_main_effects_20260507/scripts/
#       stageF3_L2_markers_risk.R (closest analog — same unit of analysis,
#       same coefficient target)
#   - Gene filter (rowSums >= 10 + >=200 expressed genes gate): same script lines 65-74
#   - eBayes robust + topTable col schema: same script lines 93, 101-108
#   - Skipped stub CSV pattern: same script line 189
#   - Symbol mapping via gene_symbol_to_ensembl.tsv: same script lines 47-52
#
# Divergences from the L2 marker risk script (intentional):
#   - Fixed formulas across L2s (not adaptive term-dropping) — required by
#     the rigor-sweep design; allows cross-L2 comparison of the same formula
#   - No dupCor (each donor appears once per L2 at L2 grain — no within-unit
#     correlation; dupCor is for the within-parent script's nhoodgroup grain)
#   - Pseudobulk built from source counts (Option A) rather than aggregated
#     from per-nhood pseudobulk RDS (Option B used by L2 markers script)
#
# Outputs:
#   outputs/design_audit.csv
#   outputs/de_results/{cohort_id}/{formula_id}/{l2_safe_name}.csv
#   outputs/de_results/{cohort_id}/{formula_id}/{l2_safe_name}.skipped.csv  (for non-testable cells)

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(edgeR)
  library(limma)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ----- args -------------------------------------------------------------
p <- ArgumentParser()
p$add_argument("--project-root", required = TRUE,
               help = "Inquiry root containing config/ and outputs/")
p$add_argument("--counts",       required = TRUE,
               help = "Pseudobulk counts CSV (rows = sample_id, cols = genes)")
p$add_argument("--pb-meta",      required = TRUE,
               help = "Pseudobulk meta CSV (one row per sample_id: cell_type, donor, n_cells, study)")
p$add_argument("--donor-meta",   required = TRUE,
               help = "Per-donor metadata CSV (includes risk_class, study, parity_binary, age_binary, sample_type, facs_status, sample_type_coarse)")
p$add_argument("--gene-mapping", required = TRUE,
               help = "gene_symbol_to_ensembl.tsv (Ensembl -> symbol; existing-script convention)")
p$add_argument("--inquiry-yaml", required = TRUE,
               help = "config/inquiry.yaml")
args <- p$parse_args()

project_root <- args$project_root
out_dir      <- file.path(project_root, "outputs")
de_dir       <- file.path(out_dir, "de_results")
dir.create(de_dir, showWarnings = FALSE, recursive = TRUE)

inquiry <- yaml::read_yaml(args$inquiry_yaml)
min_cells_l2     <- inquiry$pseudobulk$min_cells_per_donor_l2
min_donors_class <- inquiry$pseudobulk$min_donors_per_class
fdr_thr          <- inquiry$significance$fdr_threshold
lfc_thr          <- inquiry$significance$logfc_threshold

cat(sprintf("=== voom_de_sweep | min_cells_per_L2_donor=%d | min_donors_per_class=%d ===\n",
            min_cells_l2, min_donors_class))

# ----- load substrate ---------------------------------------------------
cat("[step 1] load pseudobulk counts ... ")
counts_df <- readr::read_csv(args$counts, show_col_types = FALSE,
                              progress = FALSE)
counts_mat <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts_mat) <- counts_df[[1]]      # sample_id
counts_mat <- t(counts_mat)                  # genes (rows) x samples (cols)
cat(sprintf("%d genes x %d samples\n", nrow(counts_mat), ncol(counts_mat)))

cat("[step 2] load pseudobulk meta ... ")
pb_meta <- readr::read_csv(args$pb_meta, show_col_types = FALSE,
                            progress = FALSE)
stopifnot(setequal(pb_meta$sample_id, colnames(counts_mat)))
pb_meta <- pb_meta[match(colnames(counts_mat), pb_meta$sample_id), ]
cat(sprintf("%d rows\n", nrow(pb_meta)))

cat("[step 3] load donor meta ... ")
donor_meta <- readr::read_csv(args$donor_meta, show_col_types = FALSE,
                               progress = FALSE)
cat(sprintf("%d donors\n", nrow(donor_meta)))

# Join donor metadata onto pseudobulk meta. The pseudobulk meta has its own
# `study` column (from --study-key in the aggregator); the donor meta also has
# `study`. dplyr auto-suffixes these to study.x/study.y which breaks formulas
# referring to `study`. Drop the pb-meta study before the join — donor-meta
# study is the canonical version.
cat("[step 4] join donor metadata onto pseudobulk meta\n")
if ("study" %in% colnames(pb_meta)) pb_meta$study <- NULL
pb_meta <- pb_meta %>%
  dplyr::left_join(donor_meta, by = c("donor" = "ihbca_donor_id"))

# Sanity: any pb samples without donor metadata?
n_missing <- sum(is.na(pb_meta$risk_class))
if (n_missing > 0) {
  cat(sprintf("  WARNING: %d pseudobulk samples have no donor metadata (donor not in cohort) - dropped per cohort filter\n",
              n_missing))
}

# Gene mapping (Ensembl -> symbol) — convention from stageF3_L2_markers_risk.R:47-52
cat("[step 5] load gene_mapping (ensembl -> symbol)\n")
gm <- readr::read_tsv(args$gene_mapping, show_col_types = FALSE,
                      progress = FALSE) %>%
  dplyr::rename_with(tolower)
ensg_col <- intersect(c("ensembl_id", "ensembl_gene_id", "gene_id", "ensembl"),
                      colnames(gm))[1]
sym_col_in <- intersect(c("gene_symbol", "symbol", "gene_name"),
                         colnames(gm))[1]
if (is.na(ensg_col) || is.na(sym_col_in)) {
  stop("gene_mapping schema unrecognized: ", paste(colnames(gm), collapse = ","))
}
# Build a unique ensg -> symbol vector (first symbol per ensg, matching existing convention)
gm_uniq <- gm %>% dplyr::group_by(.data[[ensg_col]]) %>%
                  dplyr::slice(1) %>% dplyr::ungroup()
gene_id_to_symbol <- setNames(gm_uniq[[sym_col_in]], gm_uniq[[ensg_col]])
cat(sprintf("  loaded %d ensg->symbol mappings\n", length(gene_id_to_symbol)))

# ----- cohort filters ---------------------------------------------------
# yaml::read_yaml returns lists for YAML arrays; unlist to vectors so %in% works.
cohort_donor_ids <- list()
for (cohort_id in names(inquiry$cohorts)) {
  cohort <- inquiry$cohorts[[cohort_id]]
  d <- donor_meta
  rc_vals <- unlist(cohort$donor_filter$risk_class)
  d <- d %>% dplyr::filter(risk_class %in% rc_vals)
  # parity_binary filter is optional (skip if yaml omits the key).
  # Metadata column may be numeric (0/1) or string ("parous"/"nulliparous"),
  # depending on the donor-meta build vintage; yaml values must match the
  # current type. When omitted, all parity statuses are included.
  if (!is.null(cohort$donor_filter$parity_binary)) {
    pb_vals <- unlist(cohort$donor_filter$parity_binary)
    d <- d %>% dplyr::filter(parity_binary %in% pb_vals)
  }
  if (!is.null(cohort$donor_filter$study$exclude)) {
    excl <- unlist(cohort$donor_filter$study$exclude)
    d <- d %>% dplyr::filter(!study %in% excl)
  }
  if (!is.null(cohort$donor_filter$study$include)) {
    incl <- unlist(cohort$donor_filter$study$include)
    d <- d %>% dplyr::filter(study %in% incl)
  }
  # Optional brca_genotype filter (added 2026-05-20 for L1 DE tested-cohort).
  # Backwards-compatible: only applies when the yaml key is set. Mirrors the
  # tested-cohort logic in scripts/07_genotyped_cohort_fibros.R:75.
  if (!is.null(cohort$donor_filter$brca_genotype)) {
    gt_vals <- unlist(cohort$donor_filter$brca_genotype)
    d <- d %>% dplyr::filter(brca_genotype %in% gt_vals)
  }
  cohort_donor_ids[[cohort_id]] <- d$ihbca_donor_id
  cat(sprintf("  cohort %-12s  n_donors = %d (AR=%d, BR1=%d)\n",
              cohort_id, nrow(d),
              sum(d$risk_class == "AR"), sum(d$risk_class == "BR1")))
}

# ----- sweep ------------------------------------------------------------
safe_name <- function(s) gsub("[^A-Za-z0-9_]", "_", s)

# Skipped-stub CSV writer (stageF3_L2_markers_risk.R:189 convention) — drops a
# tiny CSV at the would-be DE output path so consumers can list .csv files and
# distinguish testable from skipped via the `skipped` column.
write_skipped_stub <- function(path, reason) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  writeLines(c("skipped,reason",
               sprintf("TRUE,%s", gsub('"', "'", reason, fixed = TRUE))),
             path)
}

audit_rows <- list()
n_done <- 0

for (formula_id in names(inquiry$formulas)) {
  spec <- inquiry$formulas[[formula_id]]
  cohort_id <- spec$cohort
  formula_str <- spec$formula
  donors_in_cohort <- cohort_donor_ids[[cohort_id]]
  cohort_out_dir <- file.path(de_dir, cohort_id, formula_id)
  dir.create(cohort_out_dir, showWarnings = FALSE, recursive = TRUE)

  cat(sprintf("\n==== formula %s (cohort %s)\n", formula_id, cohort_id))
  cat(sprintf("     %s\n", formula_str))

  # Per-L2 loop
  l2_types <- sort(unique(pb_meta$cell_type))
  for (l2 in l2_types) {
    out_csv <- file.path(cohort_out_dir, paste0(safe_name(l2), ".csv"))

    # Subset to this L2 x cohort donors
    keep <- pb_meta$cell_type == l2 &
            pb_meta$donor %in% donors_in_cohort &
            pb_meta$n_cells >= min_cells_l2 &
            !is.na(pb_meta$risk_class)
    pb_sub <- pb_meta[keep, , drop = FALSE]
    n_AR  <- sum(pb_sub$risk_class == "AR")
    n_BR1 <- sum(pb_sub$risk_class == "BR1")

    audit <- list(
      cohort_id   = cohort_id,
      formula_id  = formula_id,
      formula     = formula_str,
      L2          = l2,
      n_AR        = n_AR,
      n_BR1       = n_BR1,
      n_donors    = n_AR + n_BR1,
      gated_out   = FALSE,
      gate_reason = NA_character_,
      design_rank = NA_integer_,
      design_ncol = NA_integer_,
      testable    = FALSE,
      aliased     = NA_character_,
      n_genes_in  = NA_integer_,
      n_sig_fdr05 = NA_integer_,
      n_sig_fdr05_lfc05 = NA_integer_
    )

    # Donor-coverage gate
    if (n_AR < min_donors_class || n_BR1 < min_donors_class) {
      audit$gated_out <- TRUE
      audit$gate_reason <- sprintf("insufficient donors (AR=%d, BR1=%d, min=%d)",
                                    n_AR, n_BR1, min_donors_class)
      write_skipped_stub(out_csv, audit$gate_reason)
      audit_rows[[length(audit_rows) + 1L]] <- audit
      next
    }

    # Build design matrix
    # Ensure risk_class is a factor with AR as the reference level
    pb_sub$risk_class <- factor(pb_sub$risk_class, levels = c("AR", "BR1"))
    # Other factor variables — coerce explicitly so model.matrix is stable
    for (col in c("study", "facs_status", "sample_type", "sample_type_coarse",
                  "parity_binary", "age_binary", "menopausal_status_binary")) {
      if (col %in% colnames(pb_sub)) {
        pb_sub[[col]] <- factor(pb_sub[[col]])
      }
    }

    f <- as.formula(formula_str)
    mm_err <- NULL
    design <- tryCatch(
      model.matrix(f, data = pb_sub),
      error = function(e) { mm_err <<- conditionMessage(e); NULL }
    )
    if (is.null(design)) {
      audit$gated_out <- TRUE
      audit$gate_reason <- sprintf("model.matrix failed: %s",
                                    ifelse(is.null(mm_err), "unknown", mm_err))
      write_skipped_stub(out_csv, audit$gate_reason)
      audit_rows[[length(audit_rows) + 1L]] <- audit
      next
    }
    # NA-dropping behavior: if model.matrix dropped rows due to NA covariate
    # values, sync pb_sub / counts to the kept rows and re-check the donor-
    # coverage gate. Only gate out if min_donors_class can no longer be met.
    # (Original behavior gated out unconditionally; relaxed 2026-05-20 for the
    # L1 sweep, where the L2 sweep's implicit parity_binary filter no longer
    # masks NA-covariate donors.)
    if (nrow(design) < nrow(pb_sub)) {
      kept_rows <- as.integer(rownames(design))
      n_dropped <- nrow(pb_sub) - nrow(design)
      pb_sub <- pb_sub[kept_rows, , drop = FALSE]
      n_AR  <- sum(pb_sub$risk_class == "AR")
      n_BR1 <- sum(pb_sub$risk_class == "BR1")
      audit$n_AR <- n_AR; audit$n_BR1 <- n_BR1; audit$n_donors <- n_AR + n_BR1
      cat(sprintf("    [NA-drop] %d rows lost to NA covariates; kept %d (AR=%d, BR1=%d)\n",
                  n_dropped, nrow(pb_sub), n_AR, n_BR1))
      if (n_AR < min_donors_class || n_BR1 < min_donors_class) {
        audit$gated_out <- TRUE
        audit$gate_reason <- sprintf("post NA-drop: AR=%d, BR1=%d, min=%d",
                                      n_AR, n_BR1, min_donors_class)
        write_skipped_stub(out_csv, audit$gate_reason)
        audit_rows[[length(audit_rows) + 1L]] <- audit
        next
      }
    }
    audit$design_ncol <- ncol(design)
    audit$design_rank <- qr(design)$rank

    if (audit$design_rank < ncol(design)) {
      # find aliased columns
      qr_obj <- qr(design)
      aliased_cols <- colnames(design)[(qr_obj$pivot[(qr_obj$rank + 1L):ncol(design)])]
      audit$aliased <- paste(aliased_cols, collapse = ";")
      audit$gated_out <- TRUE
      audit$gate_reason <- "rank-deficient design"
      write_skipped_stub(out_csv, audit$gate_reason)
      audit_rows[[length(audit_rows) + 1L]] <- audit
      next
    }

    # Ensure risk_classBR1 coefficient is present
    if (!"risk_classBR1" %in% colnames(design)) {
      audit$gated_out <- TRUE
      audit$gate_reason <- "risk_classBR1 not in design"
      write_skipped_stub(out_csv, audit$gate_reason)
      audit_rows[[length(audit_rows) + 1L]] <- audit
      next
    }

    # Subset counts to these samples; drop zero-sum donors (no cells for this L2
    # despite passing min_cells gate after CSV round-trip)
    cm <- counts_mat[, pb_sub$sample_id, drop = FALSE]
    cm <- round(cm)  # CSV-roundtripped counts may be float
    keep_col <- colSums(cm) > 0
    if (!all(keep_col)) {
      cm <- cm[, keep_col, drop = FALSE]
      pb_sub <- pb_sub[keep_col, , drop = FALSE]
      design <- design[keep_col, , drop = FALSE]
      # Re-audit n after dropping
      n_AR  <- sum(pb_sub$risk_class == "AR")
      n_BR1 <- sum(pb_sub$risk_class == "BR1")
      audit$n_AR <- n_AR; audit$n_BR1 <- n_BR1; audit$n_donors <- n_AR + n_BR1
      if (n_AR < min_donors_class || n_BR1 < min_donors_class) {
        audit$gated_out <- TRUE
        audit$gate_reason <- sprintf("post zero-sum drop: AR=%d, BR1=%d, min=%d",
                                      n_AR, n_BR1, min_donors_class)
        write_skipped_stub(out_csv, audit$gate_reason)
        audit_rows[[length(audit_rows) + 1L]] <- audit
        next
      }
    }

    # Gene-level expression gate (stageF3_L2_markers_risk.R:65-66):
    # skip if fewer than 200 genes expressed in >=3 samples
    n_expressed_genes <- sum(rowSums(cm > 0) >= 3)
    if (n_expressed_genes < 200) {
      audit$gated_out <- TRUE
      audit$gate_reason <- sprintf("only %d genes expressed in >=3 samples",
                                    n_expressed_genes)
      write_skipped_stub(out_csv, audit$gate_reason)
      audit_rows[[length(audit_rows) + 1L]] <- audit
      next
    }

    # voom -> lmFit -> eBayes -> topTable
    # Gene filter: rowSums(counts) >= 10 with keep.lib.sizes=FALSE
    # (stageF3_L2_markers_risk.R:74). NOT filterByExpr.
    # NB: don't use <<- inside tryCatch handlers — handler envs don't find
    # the outer `audit` list. Capture state in local closures instead.
    voom_error <- NULL
    n_genes_after_filter <- NA_integer_
    de_result <- tryCatch({
      dge <- DGEList(counts = cm)
      dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
      dge <- calcNormFactors(dge, method = "TMM")
      v   <- voom(dge, design = design)
      fit <- lmFit(v, design = design)
      fit <- eBayes(fit, robust = TRUE)
      tt  <- topTable(fit, coef = "risk_classBR1", n = Inf, sort.by = "none")
      tt$gene_id <- rownames(tt)
      n_genes_after_filter <<- nrow(tt)
      tt
    }, error = function(e) {
      voom_error <<- conditionMessage(e)
      NULL
    })

    audit$n_genes_in <- n_genes_after_filter

    if (is.null(de_result)) {
      audit$gated_out <- TRUE
      audit$gate_reason <- sprintf("voom/lmFit error: %s",
                                    ifelse(is.null(voom_error), "unknown", voom_error))
      write_skipped_stub(out_csv, audit$gate_reason)
      audit_rows[[length(audit_rows) + 1L]] <- audit
      next
    }

    # Annotate symbol (stageF3_L2_markers_risk.R:104)
    de_result$symbol <- gene_id_to_symbol[de_result$gene_id]
    de_result$FDR    <- de_result$adj.P.Val
    de_result$PValue <- de_result$P.Value
    de_result$cohort_id  <- cohort_id
    de_result$formula_id <- formula_id
    de_result$L2         <- l2
    de_result$n_AR       <- n_AR
    de_result$n_BR1      <- n_BR1

    audit$testable <- TRUE
    audit$n_sig_fdr05 <- sum(de_result$adj.P.Val < fdr_thr, na.rm = TRUE)
    audit$n_sig_fdr05_lfc05 <- sum(
      de_result$adj.P.Val < fdr_thr & abs(de_result$logFC) > lfc_thr,
      na.rm = TRUE
    )

    # Save DE CSV (out_csv defined at top of L2 loop)
    readr::write_csv(de_result, out_csv)
    n_done <- n_done + 1L

    audit_rows[[length(audit_rows) + 1L]] <- audit
  }
}

# ----- write audit ------------------------------------------------------
cat(sprintf("\n[step final] write design_audit.csv (%d cells in ladder x L2; %d DE CSVs written)\n",
            length(audit_rows), n_done))
audit_df <- dplyr::bind_rows(lapply(audit_rows, as.data.frame))
readr::write_csv(audit_df, file.path(out_dir, "design_audit.csv"))
cat(sprintf("  audit: %s\n", file.path(out_dir, "design_audit.csv")))

# Quick summary
cat("\n=== summary ===\n")
summary_table <- audit_df %>%
  dplyr::group_by(cohort_id, formula_id) %>%
  dplyr::summarise(
    n_L2_total    = dplyr::n(),
    n_L2_testable = sum(testable),
    n_L2_gated    = sum(gated_out),
    median_n_sig  = stats::median(n_sig_fdr05[testable], na.rm = TRUE),
    .groups = "drop"
  )
print(summary_table, n = Inf)

cat("\nDone.\n")
