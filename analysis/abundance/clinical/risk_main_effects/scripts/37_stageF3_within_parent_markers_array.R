#!/usr/bin/env Rscript
# 37_stageF3_within_parent_markers_array.R
# Stage F.3 — within-parent markers per nhoodgroup.
# SLURM array task: ONE (contrast × group) per task. The manifest is built
# deterministically inside this script from the union of all Stage F.1 outputs
# and the coverage thresholds in inquiry$stageF$within_parent_markers.
#
# For each task (contrast × group):
#   - Load that contrast's pseudobulk.rds
#   - Subset to samples whose parent_L2_joint matches the group's parent
#   - Build "this group vs OTHER groups in same parent" contrast
#   - edgeR glmQLF
#   - Annotate symbols + dissoc-stress flag
#   - Write markers/<contrast>/<group>_vs_parent.csv

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(edgeR); library(limma)
  library(dplyr); library(readr); library(yaml); library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "37_stageF3_within_parent_markers_array.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

TASK_ID <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
if (TASK_ID < 1) stop("SLURM_ARRAY_TASK_ID must be >= 1")

# -- Build manifest of (contrast × group) tasks -------------------------------
fb_thresholds <- inquiry$stageF$within_parent_markers %||% list()
min_count        <- fb_thresholds$min_count        %||% 5
min_total_count  <- fb_thresholds$min_total_count  %||% 15
fdr_thresh       <- fb_thresholds$fdr              %||% 0.05
min_donors_grp   <- inquiry$stageF$pseudobulk$min_donors_per_group %||% 10

manifest_rows <- list()
contrast_dirs <- list.dirs(paths$outputs$stageF1, recursive = FALSE)
for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  gs_csv <- file.path(cdir, "nhood_groups_summary.csv")
  if (!file.exists(gs_csv)) next
  gs <- read_csv(gs_csv, show_col_types = FALSE)
  for (k in seq_len(nrow(gs))) {
    r <- gs[k, ]
    manifest_rows[[length(manifest_rows) + 1]] <- tibble(
      contrast       = cname,
      NhoodGroup     = r$NhoodGroup,
      group_renamed  = r$NhoodGroup_renamed,
      parent_L2      = r$parent_L2_joint,
      n_nhoods       = r$n_nhoods_in_group
    )
  }
}
manifest <- bind_rows(manifest_rows) %>% arrange(contrast, parent_L2, NhoodGroup)
n_total <- nrow(manifest)
# Filter to limma-viable groups (n_nhoods >= 10). Smaller groups don't sustain
# (group × donor) pseudobulk DE — limma errors on too-few-samples-in-parent.
manifest <- manifest %>% filter(n_nhoods >= 10)
cat(sprintf("Manifest: %d viable tasks (filtered from %d total via n_nhoods >= 10)\n",
            nrow(manifest), n_total))

if (TASK_ID > nrow(manifest)) {
  cat(sprintf("TASK_ID %d > manifest size %d — exiting\n", TASK_ID, nrow(manifest)))
  quit(status = 0)
}
row <- manifest[TASK_ID, ]
cat(sprintf("\n=== TASK %d: contrast=%s group=%s parent=%s ===\n",
            TASK_ID, row$contrast, row$group_renamed, row$parent_L2))

out_sub <- file.path(paths$outputs$stageF3, row$contrast)
ensure_dir(out_sub)
out_csv <- file.path(out_sub,
            sprintf("%s_vs_parent.csv",
                    gsub("[^A-Za-z0-9]", "_", row$group_renamed)))
if (file.exists(out_csv) && file.size(out_csv) > 1000) {
  cat("Output exists, skipping:", out_csv, "\n")
  quit(status = 0)
}

# -- Load pseudobulk for this contrast ----------------------------------------
pb_path <- file.path(paths$outputs$stageF2, row$contrast, "pseudobulk.rds")
if (!file.exists(pb_path)) stop("pseudobulk missing: ", pb_path)
pb <- readRDS(pb_path)
cd <- as.data.frame(colData(pb))
cat(sprintf("pseudobulk: %d genes × %d samples; cols=%s\n",
            nrow(pb), ncol(pb), paste(colnames(cd), collapse=", ")))

# Restrict to samples whose parent matches this group's parent
keep <- !is.na(cd$parent_L2_joint) & cd$parent_L2_joint == row$parent_L2
pb_p <- pb[, keep]
cd_p <- as.data.frame(colData(pb_p))
cat(sprintf("samples in same parent (%s): %d (group renamed %s, others)\n",
            row$parent_L2, ncol(pb_p), row$group_renamed))
if (ncol(pb_p) < 6) stop("Too few samples in parent")

# Define this-group vs other-groups
cd_p$grp <- factor(ifelse(cd_p$NhoodGroup_renamed == row$group_renamed,
                          "this", "other"),
                   levels = c("other", "this"))
n_this  <- sum(cd_p$grp == "this")
n_other <- sum(cd_p$grp == "other")
cat(sprintf("group×donor counts: this=%d  other=%d\n", n_this, n_other))
if (n_this < min_donors_grp || n_other < min_donors_grp) {
  cat(sprintf("Insufficient donor coverage (this=%d, other=%d, min=%d) — skipping\n",
              n_this, n_other, min_donors_grp))
  # Write a stub so we don't reattempt
  write_csv(tibble(skipped = TRUE,
                   reason = sprintf("n_this=%d n_other=%d min=%d",
                                    n_this, n_other, min_donors_grp)),
            out_csv)
  quit(status = 0)
}

# -- DE: limma-voom (preferred over edgeR-QLF for small-N pseudobulk) --------
# voom estimates mean-variance trend on logCPM and uses standard weighted
# linear model + empirical Bayes shrinkage. More robust FDR behavior at
# this scale than edgeR's NB GLM.
cnt <- assay(pb_p, "counts")
keep_g <- filterByExpr(cnt, group = cd_p$grp,
                       min.count = min_count, min.total.count = min_total_count)
cat(sprintf("genes after filterByExpr: %d / %d\n", sum(keep_g), nrow(cnt)))
cnt <- cnt[keep_g, ]

d <- DGEList(counts = cnt)
d <- calcNormFactors(d, method = "TMM")

# Design: study as fixed batch covariate + grp as test factor; donor blocked
# via duplicateCorrelation to handle the paired-sample structure (one donor
# contributes multiple NG samples within a parent L2). Mirrors the convention
# in publication/analysis/annotation/flex/scripts/03b_limma_markers.R
# (~ library_id + condition + dupCor(patient_id)). study is NA-coalesced to a
# single level so the model is well-defined when study isn't populated.
cd_p$study <- ifelse(is.na(cd_p$study) | cd_p$study == "", "_unknown", cd_p$study)
use_study <- length(unique(cd_p$study)) > 1
if (use_study) {
  des <- model.matrix(~ study + grp, data = cd_p)
} else {
  des <- model.matrix(~ grp, data = cd_p)
}

# Donor block (dupCor) requires >=2 donors with multiple samples
donor_counts <- table(cd_p$donor)
n_repeated_donors <- sum(donor_counts > 1)
use_dupcor <- n_repeated_donors >= 2

if (use_dupcor) {
  v0 <- voom(d, design = des, plot = FALSE)
  corfit <- tryCatch(
    duplicateCorrelation(v0, design = des, block = cd_p$donor),
    error = function(e) {
      cat(sprintf("dupCor failed (%s) — falling back to no-block\n",
                  conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(corfit) && is.finite(corfit$consensus.correlation)) {
    consensus_cor <- corfit$consensus.correlation
    v <- voom(d, design = des, plot = FALSE,
              block = cd_p$donor, correlation = consensus_cor)
    fit <- lmFit(v, design = des,
                 block = cd_p$donor, correlation = consensus_cor)
    cat(sprintf("design: %s + dupCor(donor cor=%.3f)\n",
                paste(colnames(des), collapse="+"), consensus_cor))
  } else {
    v <- voom(d, design = des, plot = FALSE)
    fit <- lmFit(v, design = des)
    cat("design: voom (dupCor unstable, no block)\n")
  }
} else {
  v <- voom(d, design = des, plot = FALSE)
  fit <- lmFit(v, design = des)
  cat(sprintf("design: voom (n_repeated_donors=%d, no dupCor)\n",
              n_repeated_donors))
}
fit <- eBayes(fit, robust = TRUE)
tt <- topTable(fit, coef = "grpthis", number = Inf, sort.by = "none")
# Standardize column names to match prior edgeR schema (logFC, FDR, etc.)
tt$gene_id <- rownames(tt)
if ("adj.P.Val" %in% colnames(tt)) tt$FDR <- tt$adj.P.Val
if ("P.Value"   %in% colnames(tt)) tt$PValue <- tt$P.Value
tt$contrast   <- row$contrast
tt$group      <- row$group_renamed
tt$parent_L2  <- row$parent_L2

# Annotate symbols + dissoc-stress flag
gmap <- read_tsv(paths$inputs$gene_mapping, show_col_types = FALSE)
e2s <- gmap %>% group_by(ensembl_id) %>% slice(1) %>% ungroup()
tt <- tt %>%
  left_join(e2s, by = c("gene_id" = "ensembl_id")) %>%
  rename(symbol = gene_symbol)
dissoc_genes <- inquiry$context$dissoc_genes %||% character(0)
tt$is_dissoc_stress <- tt$symbol %in% dissoc_genes

write_csv(tt, out_csv)
n_sig    <- sum(tt$FDR < fdr_thresh, na.rm = TRUE)
n_sig_ds <- sum(tt$FDR < fdr_thresh & tt$is_dissoc_stress, na.rm = TRUE)
cat(sprintf("Wrote: %s  (sig FDR<%.2f: %d  dissoc-stress sig: %d)\n",
            out_csv, fdr_thresh, n_sig, n_sig_ds))
cat("=== Stage F.3 task done ===\n")
