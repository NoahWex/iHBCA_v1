#!/usr/bin/env Rscript
# stageI2_fibro_sfrp4_BR1vsAR_parous.R
#
# Cross-cohort, parous-stratified DE for str::Fibro-SFRP4:
#   AR-parous donors vs BRCA1-parous donors.
#
# This is *not* a parity main-effect test — both groups are parous. The contrast
# asks: given parous status, does germline BRCA1 carrier state alter the
# Fibro-SFRP4 expression program?
#
# Design follows the F.3 L2 markers approach (sum F.2 NhoodGroup pseudobulks
# to L2 grain per donor, then limma-voom on donor-level pseudobulks). Cancer
# history is quasi-aliased with cohort here (AR donors are required to have
# no cancer history) so it is dropped if it has no within-cohort variation
# after subsetting to parous.
#
# Outputs:
#   outputs/stageI2_within_l2_de/fibro_sfrp4_BR1vsAR_parous.csv  topTable
#   outputs/stageI2_within_l2_de/fibro_sfrp4_BR1vsAR_parous_design.tsv  design info

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(SummarizedExperiment); library(Matrix)
  library(dplyr); library(readr); library(tidyr); library(yaml)
  library(edgeR); library(limma)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "stageI2_fibro_sfrp4_BR1vsAR_parous.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

L2_TARGET   <- "str::Fibro-SFRP4"
COHORT_MAP  <- list(
  parity_in_AR      = list(short = "AR",  cohort_csv = "AR_baseline"),
  parity_x_HR_BRCA1 = list(short = "BR1", cohort_csv = "HR_BRCA1")
)
# Cohort design column names (canonical from cohort_*_design.csv)
BASE_COVARS  <- c("study", "facs_status", "age_binary", "menopausal_status_binary")
EXTRA_COVARS <- c("cancer_history")
# Drop priority for progressive simplification when design is rank-deficient
# (drop these in order; "study" + "cohort" are the minimum we keep)
DROP_PRIORITY <- c("cancer_history", "facs_status",
                    "menopausal_status_binary", "age_binary")
MIN_DONORS_PER_GROUP <- 2   # BR1 parous donor pool is small

# ---- gene_id -> symbol ----
# Prefer the canonical share/kai_v1 mapping (already-vetted, used by other scripts).
gm_path_share <- file.path(paths$inquiry_root, "share", "kai_v1", "06_markers",
                            "gene_symbol_to_ensembl.tsv")
gm_path_default <- paths$inputs$gene_mapping
gm_path <- if (file.exists(gm_path_share)) gm_path_share else gm_path_default
gm <- read_tsv(gm_path, show_col_types = FALSE) %>%
  rename_with(tolower)
# Normalize possible column-name variants
ens_col <- intersect(c("ensembl_id", "ensembl", "ensembl_gene_id", "gene_id"),
                      colnames(gm))[1]
sym_col <- intersect(c("gene_symbol", "symbol", "gene_name"),
                      colnames(gm))[1]
if (is.na(ens_col) || is.na(sym_col))
  stop(sprintf("gene mapping at %s has unexpected schema: %s",
                gm_path, paste(colnames(gm), collapse = ",")))
gene_id_to_symbol <- setNames(gm[[sym_col]], gm[[ens_col]])
cat(sprintf("Loaded gene mapping: %s (%d entries)\n", gm_path, length(gene_id_to_symbol)))

# ---- per-cohort pseudobulk aggregation (donor-level, sum across NhoodGroups) ----
agg_donor_pb <- function(cname) {
  pb_path <- file.path(paths$outputs$stageF2, cname, "pseudobulk.rds")
  if (!file.exists(pb_path)) stop(sprintf("Missing F.2 pseudobulk: %s", pb_path))
  pb <- readRDS(pb_path)
  cd <- as.data.frame(colData(pb))
  if (!"parent_L2_joint" %in% colnames(cd)) {
    cd$parent_L2_joint <- cd$L2_joint %||% cd$parent_L2 %||% NA_character_
  }
  donor_col <- intersect(c("ihbca_donor_id", "donor_id", "patientID", "donor"),
                          colnames(cd))[1]
  cd$donor_id <- cd[[donor_col]]
  cols_in_L2 <- which(cd$parent_L2_joint == L2_TARGET)
  if (length(cols_in_L2) < 2)
    stop(sprintf("Too few NhoodGroups for %s in %s", L2_TARGET, cname))
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

# ---- combine AR + BR1 ----
ar_pb  <- agg_donor_pb("parity_in_AR")
br1_pb <- agg_donor_pb("parity_x_HR_BRCA1")

# Same gene order assumed (both come from same Milo object).
stopifnot(identical(rownames(ar_pb$counts), rownames(br1_pb$counts)))

# Cohort designs (parous filter)
read_design <- function(cohort_csv) {
  fp <- file.path(paths$outputs$stageA,
                   paste0("cohort_", cohort_csv, "_design.csv"))
  d <- read_csv(fp, show_col_types = FALSE)
  donor_key <- intersect(c("ihbca_donor_id", "donor_id", "patientID"),
                          colnames(d))[1]
  d$donor_id <- d[[donor_key]]
  d
}
ar_design  <- read_design("AR_baseline")  %>% filter(parity_binary == "parous")
br1_design <- read_design("HR_BRCA1")     %>% filter(parity_binary == "parous")

cat(sprintf("AR  parous donors in cohort design: %d\n", nrow(ar_design)))
cat(sprintf("BR1 parous donors in cohort design: %d\n", nrow(br1_design)))

# Restrict pseudobulk donors to the parous set in each cohort
keep_ar  <- intersect(colnames(ar_pb$counts),  ar_design$donor_id)
keep_br1 <- intersect(colnames(br1_pb$counts), br1_design$donor_id)
cat(sprintf("AR  parous donors with Fibro-SFRP4 pseudobulks: %d\n", length(keep_ar)))
cat(sprintf("BR1 parous donors with Fibro-SFRP4 pseudobulks: %d\n", length(keep_br1)))

if (length(keep_ar) < MIN_DONORS_PER_GROUP || length(keep_br1) < MIN_DONORS_PER_GROUP)
  stop(sprintf("Too few parous donors with pseudobulks: AR=%d BR1=%d",
               length(keep_ar), length(keep_br1)))

ar_counts  <- ar_pb$counts[,  keep_ar,  drop = FALSE]
br1_counts <- br1_pb$counts[, keep_br1, drop = FALSE]
counts_mat <- cbind(ar_counts, br1_counts)

# Build sample meta in the same column order
ar_meta  <- tibble(donor_id = keep_ar,  cohort = "AR")  %>%
  left_join(ar_design, by = "donor_id")
br1_meta <- tibble(donor_id = keep_br1, cohort = "BR1") %>%
  left_join(br1_design, by = "donor_id")
meta <- bind_rows(ar_meta, br1_meta)
stopifnot(identical(meta$donor_id, colnames(counts_mat)))

meta$cohort <- factor(meta$cohort, levels = c("AR", "BR1"))   # AR is reference

# Replace NAs in covariate columns with "unknown" so model.matrix doesn't
# silently drop rows. This is safer than excluding NA-bearing donors when
# n is already small.
covar_cols <- intersect(c(BASE_COVARS, EXTRA_COVARS), colnames(meta))
for (cv in covar_cols) {
  v <- meta[[cv]]
  if (is.character(v) || is.factor(v)) {
    v <- as.character(v)
    v[is.na(v) | v == ""] <- "unknown"
    meta[[cv]] <- v
  }
}
cat("\nSample-level summary (cohort × parity):\n")
print(meta %>% count(cohort, parity_binary))
present_covars <- intersect(c(BASE_COVARS, EXTRA_COVARS), colnames(meta))
if (length(present_covars) > 0) {
  cat("Covariate distribution by cohort:\n")
  for (cv in present_covars) {
    cat(sprintf("  %s:\n", cv))
    print(meta %>% count(cohort, !!sym(cv)))
  }
}

# ---- limma-voom ----
dge <- DGEList(counts = counts_mat)
dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat(sprintf("\nGenes after >=10-total filter: %d\n", nrow(dge)))

# Build initial term set: only covariates present + with variation
candidate_terms <- c(BASE_COVARS, EXTRA_COVARS, "cohort")
use_terms <- candidate_terms[sapply(candidate_terms, function(x) {
  if (!x %in% colnames(meta)) return(FALSE)
  length(unique(meta[[x]])) > 1
})]

# Progressive simplification: build design; if rank-deficient, drop the
# highest-priority covariate from DROP_PRIORITY and retry. Always keep cohort.
build_design <- function(terms, dat) {
  fstr <- paste("~", paste(terms, collapse = " + "))
  list(fstr = fstr, design = model.matrix(as.formula(fstr), data = dat))
}

design <- NULL
formula_str <- NULL
for (attempt in seq_len(length(DROP_PRIORITY) + 1)) {
  bd <- build_design(use_terms, meta)
  if (qr(bd$design)$rank == ncol(bd$design)) {
    design <- bd$design; formula_str <- bd$fstr; break
  }
  cat(sprintf("Rank deficient with [%s] — dropping covariate.\n", bd$fstr))
  to_drop <- intersect(DROP_PRIORITY, use_terms)
  if (length(to_drop) == 0) break
  use_terms <- setdiff(use_terms, to_drop[1])
}
if (is.null(design)) stop("Could not build full-rank design even after simplification.")
cat(sprintf("Final formula: %s\n", formula_str))

v <- voom(dge, design)
fit <- lmFit(v, design)
fit <- eBayes(fit, robust = TRUE)

coef_name <- grep("^cohortBR1$", colnames(design), value = TRUE)
if (length(coef_name) != 1)
  stop(sprintf("cohortBR1 coef not unique in design: %s",
               paste(colnames(design), collapse = ",")))

tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none") %>%
  as.data.frame() %>%
  mutate(gene_id = rownames(.),
          symbol = gene_id_to_symbol[gene_id],
          FDR = adj.P.Val,
          PValue = P.Value,
          parent_L2 = L2_TARGET,
          comparison = "BR1_parous_vs_AR_parous")
rownames(tt) <- NULL

out_dir <- file.path(paths$inquiry_root, "outputs", "stageI2_within_l2_de")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_csv <- file.path(out_dir, "fibro_sfrp4_BR1vsAR_parous.csv")
write_csv(tt, out_csv)
cat(sprintf("\nWrote: %s (rows=%d, sig FDR<0.05: %d)\n",
            out_csv, nrow(tt), sum(tt$FDR < 0.05, na.rm = TRUE)))

# Design-info sidecar
info <- tibble(
  comparison = "BR1_parous_vs_AR_parous",
  L2 = L2_TARGET,
  n_AR_donors = length(keep_ar),
  n_BR1_donors = length(keep_br1),
  formula = formula_str,
  coef_tested = coef_name,
  n_genes_tested = nrow(tt),
  n_sig_FDR05 = sum(tt$FDR < 0.05, na.rm = TRUE)
)
write_tsv(info, file.path(out_dir, "fibro_sfrp4_BR1vsAR_parous_design.tsv"))

# Sanity prints
cat("\n=== Top 20 UP in BR1 (vs AR), parous-only, FDR<0.05 ===\n")
print(tt %>% filter(FDR < 0.05, logFC > 0) %>%
        arrange(desc(t)) %>% head(20) %>%
        select(symbol, gene_id, logFC, AveExpr, t, FDR))

cat("\n=== Top 20 DOWN in BR1 (vs AR), parous-only, FDR<0.05 ===\n")
print(tt %>% filter(FDR < 0.05, logFC < 0) %>%
        arrange(t) %>% head(20) %>%
        select(symbol, gene_id, logFC, AveExpr, t, FDR))

cat("\n=== stageI2_fibro_sfrp4_BR1vsAR_parous done ===\n")
