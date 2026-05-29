#!/usr/bin/env Rscript
# Pairwise NhoodGroup head-to-head DE (limma-trend on logcounts).
#
# Subset pseudobulk to buckets in {group_A, group_B}, normalize on that subset,
# then a single TestA - TestB contrast. This is the direct "what genes distinguish
# subgroup A from subgroup B within the same parent" test, complementing the
# one-vs-rest in markers.csv.

suppressPackageStartupMessages({
  library(argparse)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(edgeR)
  library(limma)
  library(data.table)
})

parser <- ArgumentParser()
parser$add_argument("--pseudobulk", required = TRUE)
parser$add_argument("--group-a",    required = TRUE)
parser$add_argument("--group-b",    required = TRUE)
parser$add_argument("--out-csv",    required = TRUE)
parser$add_argument("--session",    default = "flex_scvi")
parser$add_argument("--min-cells-per-bucket", type = "integer", default = 5)
parser$add_argument("--prior-count", type = "double", default = 2)
args <- parser$parse_args()

sce <- readRDS(args$pseudobulk)
cd  <- as.data.frame(colData(sce))
keep <- cd$n_cells >= args$min_cells_per_bucket &
        cd$NhoodGroup %in% c(args$group_a, args$group_b)
sce <- sce[, keep]; cd <- cd[keep, , drop = FALSE]
message(sprintf("buckets: %s=%d, %s=%d",
                args$group_a, sum(cd$NhoodGroup == args$group_a),
                args$group_b, sum(cd$NhoodGroup == args$group_b)))
stopifnot(sum(cd$NhoodGroup == args$group_a) >= 2,
          sum(cd$NhoodGroup == args$group_b) >= 2)

counts_mat <- assay(sce, "counts")
lib_size   <- Matrix::colSums(counts_mat)
dge <- DGEList(counts = as.matrix(counts_mat), lib.size = lib_size)
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = args$prior_count)

i_meta <- data.frame(Test = ifelse(cd$NhoodGroup == args$group_a, "A", "B"))
i_model <- model.matrix(~ 0 + Test, data = i_meta)
rownames(i_model) <- colnames(logcpm)

i_fit <- lmFit(logcpm, i_model)
i_contrast <- makeContrasts(contrasts = "TestA - TestB", levels = i_model)
i_fit <- contrasts.fit(i_fit, i_contrast)
i_fit <- eBayes(i_fit, trend = TRUE)
i_res <- as.data.frame(topTreat(i_fit, number = Inf, sort.by = "none"))

out <- data.table(
  group_A   = args$group_a,
  group_B   = args$group_b,
  gene      = rownames(i_res),
  logFC     = i_res$logFC,         # positive = up in A vs B
  AveExpr   = i_res$AveExpr,
  t         = i_res$t,
  P.Value   = i_res$P.Value,
  adj.P.Val = i_res$adj.P.Val,
  B         = i_res$B,
  n_A       = sum(cd$NhoodGroup == args$group_a),
  n_B       = sum(cd$NhoodGroup == args$group_b),
  session   = args$session
)
fwrite(out, args$out_csv)
message(sprintf("wrote: %s (%d rows)", args$out_csv, nrow(out)))
