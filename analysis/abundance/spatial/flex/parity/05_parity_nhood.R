#!/usr/bin/env Rscript
# 05_parity_nhood.R — Track D deliverable A
#
# Nhood-level parity between joint embedding DA and FLEX-only DA, per L1.5 class.
# Both da_results.csv files carry per-nhood logFC, SpatialFDR, nhood_label_l1p5,
# nhood_label_l1p5_frac. The two assemblies build different nhoods on different
# embeddings, so we compare distributions per L1.5 class rather than 1:1 nhoods.
#
# Three views per L1.5 class:
#   1. Composition: KS test on nhood_label_l1p5_frac (purity)
#   2. Significance: counts and fractions of SpatialFDR < threshold by direction
#   3. LogFC: KS on logFC distributions, mean/median, sign agreement at the
#      class-aggregate level (sign of mean logFC).

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
})

parser <- ArgumentParser()
parser$add_argument("--joint-da", required = TRUE,
                    help = "Path to joint da_results.csv")
parser$add_argument("--flex-da", required = TRUE,
                    help = "Path to FLEX da_results.csv")
parser$add_argument("--out-dir", required = TRUE,
                    help = "Output directory (created if missing)")
parser$add_argument("--contrast", required = TRUE,
                    help = "Contrast id, written into outputs (e.g., C1_p1p2_vs_p3)")
parser$add_argument("--scope", default = "flex_full",
                    help = "Scope tag (default flex_full)")
parser$add_argument("--tier", default = "M0",
                    help = "Tier tag (default M0)")
parser$add_argument("--fdr", type = "double", default = 0.1,
                    help = "SpatialFDR threshold for significance (default 0.1)")
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

joint <- fread(args$joint_da)
flex  <- fread(args$flex_da)

stopifnot(all(c("logFC", "SpatialFDR", "Nhood",
                "nhood_label_l1p5", "nhood_label_l1p5_frac") %in% names(joint)))
stopifnot(all(c("logFC", "SpatialFDR", "Nhood",
                "nhood_label_l1p5", "nhood_label_l1p5_frac") %in% names(flex)))

joint[, assembly := "joint"]
flex[,  assembly := "flex"]

cat(sprintf("[parity] joint nhoods: %d   FLEX nhoods: %d\n", nrow(joint), nrow(flex)))
cat(sprintf("[parity] L1.5 classes — joint: %d   FLEX: %d   union: %d\n",
            joint[, uniqueN(nhood_label_l1p5)],
            flex[,  uniqueN(nhood_label_l1p5)],
            length(union(joint$nhood_label_l1p5, flex$nhood_label_l1p5))))

l1p5_classes <- sort(union(unique(joint$nhood_label_l1p5),
                           unique(flex$nhood_label_l1p5)))

# ---- View 1: composition (purity KS per class) ------------------------------

ks_safe <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) {
    return(list(stat = NA_real_, p = NA_real_))
  }
  res <- suppressWarnings(stats::ks.test(x, y))
  list(stat = unname(res$statistic), p = res$p.value)
}

composition_rows <- lapply(l1p5_classes, function(cls) {
  jp <- joint[nhood_label_l1p5 == cls, nhood_label_l1p5_frac]
  fp <- flex [nhood_label_l1p5 == cls, nhood_label_l1p5_frac]
  ks <- ks_safe(jp, fp)
  data.table(
    l1p5_class = cls,
    n_joint = length(jp), n_flex = length(fp),
    mean_purity_joint  = if (length(jp)) mean(jp) else NA_real_,
    mean_purity_flex   = if (length(fp)) mean(fp) else NA_real_,
    median_purity_joint = if (length(jp)) median(jp) else NA_real_,
    median_purity_flex  = if (length(fp)) median(fp) else NA_real_,
    ks_stat = ks$stat, ks_p = ks$p
  )
})
composition <- rbindlist(composition_rows)
composition[, ks_p_adj := p.adjust(ks_p, method = "BH")]

# ---- View 2: significance counts per class ----------------------------------

sig_summary <- function(dt, cls, fdr) {
  sub <- dt[nhood_label_l1p5 == cls]
  n_total <- nrow(sub)
  is_sig <- sub$SpatialFDR < fdr
  n_sig <- sum(is_sig, na.rm = TRUE)
  n_sig_pos <- sum(is_sig & sub$logFC > 0, na.rm = TRUE)
  n_sig_neg <- sum(is_sig & sub$logFC < 0, na.rm = TRUE)
  data.table(n_total = n_total, n_sig = n_sig,
             n_sig_pos = n_sig_pos, n_sig_neg = n_sig_neg,
             sig_frac = if (n_total) n_sig / n_total else NA_real_)
}

significance_rows <- lapply(l1p5_classes, function(cls) {
  j <- sig_summary(joint, cls, args$fdr)
  f <- sig_summary(flex,  cls, args$fdr)
  data.table(
    l1p5_class = cls,
    n_joint           = j$n_total,    n_flex           = f$n_total,
    n_sig_joint       = j$n_sig,      n_sig_flex       = f$n_sig,
    n_sig_pos_joint   = j$n_sig_pos,  n_sig_pos_flex   = f$n_sig_pos,
    n_sig_neg_joint   = j$n_sig_neg,  n_sig_neg_flex   = f$n_sig_neg,
    sig_frac_joint    = j$sig_frac,   sig_frac_flex    = f$sig_frac
  )
})
significance <- rbindlist(significance_rows)

# Class-level dominant-direction agreement: do joint and FLEX agree on which
# direction (sig+, sig-, ns) dominates the class?
significance[, dominant_dir_joint := fcase(
  n_sig_pos_joint == 0 & n_sig_neg_joint == 0, "ns",
  n_sig_pos_joint >= n_sig_neg_joint, "pos",
  default = "neg"
)]
significance[, dominant_dir_flex := fcase(
  n_sig_pos_flex == 0 & n_sig_neg_flex == 0, "ns",
  n_sig_pos_flex >= n_sig_neg_flex, "pos",
  default = "neg"
)]
significance[, dominant_dir_agree := dominant_dir_joint == dominant_dir_flex]

# ---- View 3: logFC distributions per class ----------------------------------

logfc_rows <- lapply(l1p5_classes, function(cls) {
  jl <- joint[nhood_label_l1p5 == cls, logFC]
  fl <- flex [nhood_label_l1p5 == cls, logFC]
  ks <- ks_safe(jl, fl)
  data.table(
    l1p5_class = cls,
    n_joint = length(jl), n_flex = length(fl),
    mean_logfc_joint    = if (length(jl)) mean(jl)   else NA_real_,
    mean_logfc_flex     = if (length(fl)) mean(fl)   else NA_real_,
    median_logfc_joint  = if (length(jl)) median(jl) else NA_real_,
    median_logfc_flex   = if (length(fl)) median(fl) else NA_real_,
    sign_mean_joint     = if (length(jl)) sign(mean(jl))   else NA_real_,
    sign_mean_flex      = if (length(fl)) sign(mean(fl))   else NA_real_,
    sign_median_joint   = if (length(jl)) sign(median(jl)) else NA_real_,
    sign_median_flex    = if (length(fl)) sign(median(fl)) else NA_real_,
    ks_stat = ks$stat, ks_p = ks$p
  )
})
logfc <- rbindlist(logfc_rows)
logfc[, ks_p_adj := p.adjust(ks_p, method = "BH")]
logfc[, sign_mean_agree   := !is.na(sign_mean_joint) & !is.na(sign_mean_flex) &
                              sign_mean_joint == sign_mean_flex]
logfc[, sign_median_agree := !is.na(sign_median_joint) & !is.na(sign_median_flex) &
                              sign_median_joint == sign_median_flex]

# ---- Parity summary ---------------------------------------------------------

summary_dt <- merge(
  composition[, .(l1p5_class, n_joint, n_flex,
                  mean_purity_joint, mean_purity_flex,
                  purity_ks_p = ks_p, purity_ks_p_adj = ks_p_adj)],
  significance[, .(l1p5_class, n_sig_joint, n_sig_flex,
                   sig_frac_joint, sig_frac_flex,
                   dominant_dir_joint, dominant_dir_flex, dominant_dir_agree)],
  by = "l1p5_class", all = TRUE
)
summary_dt <- merge(
  summary_dt,
  logfc[, .(l1p5_class, mean_logfc_joint, mean_logfc_flex,
            sign_mean_agree, sign_median_agree,
            logfc_ks_p = ks_p, logfc_ks_p_adj = ks_p_adj)],
  by = "l1p5_class", all = TRUE
)

# Headline metrics
n_classes <- nrow(summary_dt)
n_dir_agree <- sum(summary_dt$dominant_dir_agree, na.rm = TRUE)
n_sign_mean_agree <- sum(summary_dt$sign_mean_agree, na.rm = TRUE)
cat(sprintf("\n[parity headline] %d L1.5 classes\n", n_classes))
cat(sprintf("  dominant-direction agreement: %d/%d (%.1f%%)\n",
            n_dir_agree, n_classes, 100 * n_dir_agree / n_classes))
cat(sprintf("  sign-of-mean-logFC agreement: %d/%d (%.1f%%)\n",
            n_sign_mean_agree, n_classes, 100 * n_sign_mean_agree / n_classes))

# ---- Write outputs ----------------------------------------------------------

fwrite(composition,  file.path(args$out_dir, "composition_per_l1p5.csv"))
fwrite(significance, file.path(args$out_dir, "significance_per_l1p5.csv"))
fwrite(logfc,        file.path(args$out_dir, "logfc_per_l1p5.csv"))
fwrite(summary_dt,   file.path(args$out_dir, "parity_summary.csv"))

# Run log
log_lines <- c(
  sprintf("contrast: %s", args$contrast),
  sprintf("scope: %s", args$scope),
  sprintf("tier: %s", args$tier),
  sprintf("fdr_threshold: %s", args$fdr),
  sprintf("joint_da: %s", args$joint_da),
  sprintf("flex_da: %s", args$flex_da),
  sprintf("n_joint_nhoods: %d", nrow(joint)),
  sprintf("n_flex_nhoods: %d", nrow(flex)),
  sprintf("n_l1p5_classes: %d", n_classes),
  sprintf("dominant_dir_agree: %d/%d", n_dir_agree, n_classes),
  sprintf("sign_mean_logfc_agree: %d/%d", n_sign_mean_agree, n_classes),
  sprintf("timestamp: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
)
writeLines(log_lines, file.path(args$out_dir, "run_log.yaml"))

cat(sprintf("\n[parity] outputs written to %s\n", args$out_dir))
