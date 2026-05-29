#!/usr/bin/env Rscript
#
# Shared GSEA utility — runs fgsea on a standardized DE-markers CSV against
# Hallmark / Reactome / GO pathway databases. Used across V1 NhoodGroup,
# V1 latent-subtype DE, FLEX-DA Track D, and any future per-group DE pipeline.
#
# I/O contract:
#   Input markers CSV must carry columns: group_id, gene, logFC, adj.P.Val
#     (additional columns ignored — these four are the contract)
#   Output CSV: group_id, pathway, NES, pval, padj, size, leadingEdge
#
# Pathway DB selection: Hallmark (default), Reactome, GO_BP, GO_MF, GO_CC.
# Pathway DB files resolved via paths.yaml tokens (gsea_hallmark, gsea_reactome,
# gsea_go_bp, gsea_go_mf, gsea_go_cc) — passed via --pathway-gmt arg.

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(fgsea)
})

parser <- ArgumentParser()
parser$add_argument("--markers-csv", required = TRUE,
                    help = "Path to DE markers CSV (group_id, gene, logFC, adj.P.Val)")
parser$add_argument("--pathway-gmt", required = TRUE,
                    help = "Path to pathway .gmt file (e.g. Hallmark, Reactome)")
parser$add_argument("--out-csv", required = TRUE,
                    help = "Path to write GSEA results CSV")
parser$add_argument("--rank-by", default = "logFC",
                    choices = c("logFC", "signed_neglog10p"),
                    help = "Gene ranking metric (default: logFC)")
parser$add_argument("--min-size", type = "integer", default = 10,
                    help = "Minimum gene set size (default: 10)")
parser$add_argument("--max-size", type = "integer", default = 500,
                    help = "Maximum gene set size (default: 500)")
parser$add_argument("--n-perm", type = "integer", default = 10000,
                    help = "fgsea permutation count (default: 10000)")
parser$add_argument("--seed", type = "integer", default = 42)
args <- parser$parse_args()

set.seed(args$seed)

# ---------------------------------------------------------------------------
# Load markers + pathway DB
# ---------------------------------------------------------------------------
mk <- fread(args$markers_csv)
required_cols <- c("group_id", "gene", "logFC", "adj.P.Val")
missing <- setdiff(required_cols, colnames(mk))
if (length(missing) > 0) {
  stop(sprintf("markers-csv missing required columns: %s",
               paste(missing, collapse = ", ")))
}
cat(sprintf("[gsea] markers: %d rows × %d groups\n", nrow(mk), uniqueN(mk$group_id)))

pathways <- fgsea::gmtPathways(args$pathway_gmt)
cat(sprintf("[gsea] pathways: %d (DB: %s)\n", length(pathways), basename(args$pathway_gmt)))

# ---------------------------------------------------------------------------
# Per-group fgsea: rank genes, run fgsea, collect
# ---------------------------------------------------------------------------
rank_gene <- function(dt, metric) {
  if (metric == "logFC") {
    setNames(dt$logFC, dt$gene)
  } else if (metric == "signed_neglog10p") {
    p <- pmax(dt$adj.P.Val, 1e-300)
    setNames(sign(dt$logFC) * -log10(p), dt$gene)
  }
}

results <- list()
for (g in unique(mk$group_id)) {
  sub <- mk[group_id == g]
  ranks <- rank_gene(sub, args$rank_by)
  ranks <- ranks[!is.na(ranks) & !is.infinite(ranks)]
  ranks <- ranks[!duplicated(names(ranks))]
  if (length(ranks) < args$min_size) {
    cat(sprintf("[gsea] skip %s (only %d genes)\n", g, length(ranks)))
    next
  }
  res <- fgsea::fgsea(pathways = pathways, stats = ranks,
                     minSize = args$min_size, maxSize = args$max_size,
                     nPermSimple = args$n_perm)
  res[, group_id := g]
  res[, leadingEdge := vapply(leadingEdge, paste, FUN.VALUE = character(1),
                              collapse = ";")]
  results[[g]] <- res
  cat(sprintf("[gsea] %s: %d pathways tested, %d sig padj<0.05\n",
              g, nrow(res), sum(res$padj < 0.05, na.rm = TRUE)))
}

out <- rbindlist(results)
out_cols <- c("group_id", "pathway", "NES", "pval", "padj", "size", "leadingEdge")
out <- out[, ..out_cols]

dir.create(dirname(args$out_csv), recursive = TRUE, showWarnings = FALSE)
fwrite(out, args$out_csv)
cat(sprintf("[gsea] wrote %d rows to %s\n", nrow(out), args$out_csv))
