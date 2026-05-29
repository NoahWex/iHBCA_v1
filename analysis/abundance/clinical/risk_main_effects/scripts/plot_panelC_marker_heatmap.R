#!/usr/bin/env Rscript
# plot_panelC_marker_heatmap.R
#
# Donor-pseudobulk × gene heatmap, columns grouped by NhoodGroup, for a
# single parent_L2 in a single contrast. Top column-annotation: NG median
# nhood logFC + NG rank.
#
# Source: Stage F.2 pseudobulk.rds (per-NG × donor counts; ~260 MB) — far
# cheaper than re-loading the full milo (which strips assays anyway). Cell-
# level granularity is approximated by donor pseudobulks: each NG has
# multiple donor columns. Z-scored log-CPM expression.
#
# Args:
#   --L2 <"epi::BMYO-basal" | "str::Fibro-major" | ...>
#   --contrast <"parity_x_HR_BRCA1" (default) | ...>
#   --top-n <N>  (top N markers per NG, default 15)
#   --inquiry-dir <path>
#
# Output:
#   reports/figures/panelC_<L2_safe>_<contrast>_heatmap.pdf
#   reports/figures/panelC_<L2_safe>_<contrast>_genes.csv

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(SingleCellExperiment); library(edgeR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(tidyr)
  library(ggplot2); library(scales)
})
`%||%` <- function(a, b) if (is.null(a)) b else a
have_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
have_circlize <- requireNamespace("circlize", quietly = TRUE)
if (!have_ch || !have_circlize)
  stop("ComplexHeatmap + circlize required.")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == name); if (length(i) == 0) return(default); args[i + 1]
}
L2_target  <- get_arg("--L2",  "epi::BMYO-basal")
contrast   <- get_arg("--contrast", "parity_x_HR_BRCA1")
top_n      <- as.integer(get_arg("--top-n", "15"))

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_panelC_marker_heatmap.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
cfg <- parse_args_inquiry()
paths <- cfg$paths

cat(sprintf("=== Panel C — L2='%s' contrast='%s' top_n=%d ===\n",
            L2_target, contrast, top_n))

L2_compartment <- strsplit(L2_target, "::", fixed = TRUE)[[1]][1]
L2_label       <- strsplit(L2_target, "::", fixed = TRUE)[[1]][2]
L2_safe        <- gsub("::", "__", gsub("[^A-Za-z0-9_]", "_", L2_target))

RANK_PALETTE <- c(
  "1" = "#332288", "2" = "#117733", "3" = "#44AA99", "4" = "#88CCEE",
  "5" = "#DDCC77", "6" = "#CC6677", "7" = "#AA4499", "8" = "#882255",
  "9" = "#999933", "10" = "#661100"
)

# ---- Lookup: viable NGs for this L2 in this contrast ----
lookup <- read_csv(file.path(paths$inquiry_root, "outputs", "nhoodgroup_renaming",
                              "lookup.csv"), show_col_types = FALSE) %>%
  filter(viable, contrast == !!contrast, parent_L2_joint == L2_target) %>%
  mutate(NhoodGroup = as.character(NhoodGroup),
         NG_rendered = paste0(parent_L2_joint, "_", NhoodGroup)) %>%
  arrange(size_rank_within_L2)
if (nrow(lookup) == 0) stop("No viable NGs for ", L2_target, " in ", contrast)
cat(sprintf("Viable NGs: %d\n", nrow(lookup)))
print(lookup %>% select(NhoodGroup, n_nhoods_in_group, group_med_lfc,
                        size_rank_within_L2, name_short))

# ---- Top markers per NG (F.3 limma-voom) ----
marker_dir <- file.path(paths$outputs$stageF3, contrast)
top_markers <- list()
for (i in seq_len(nrow(lookup))) {
  ng <- lookup$NhoodGroup[i]; rk <- lookup$size_rank_within_L2[i]
  fname <- sprintf("%s__%s_%s_vs_parent.csv",
                    L2_compartment,
                    gsub("-", "_", L2_label),
                    ng)
  fp <- file.path(marker_dir, fname)
  if (!file.exists(fp)) next
  m <- read_csv(fp, show_col_types = FALSE) %>%
    filter(!is_dissoc_stress, !is.na(symbol), symbol != "",
           FDR < 0.05, abs(logFC) > log2(1.25)) %>%
    arrange(desc(abs(logFC))) %>%
    head(top_n) %>%
    mutate(rank = rk, ng = ng)
  top_markers[[ng]] <- m
}
markers_df <- bind_rows(top_markers)
gene_union <- unique(markers_df$symbol)
cat(sprintf("Marker union (top %d per NG, FDR<0.05, |lfc|>0.32): %d genes\n",
            top_n, length(gene_union)))

# ---- Load Stage F.2 pseudobulk for this contrast ----
pb_path <- file.path(paths$outputs$stageF2, contrast, "pseudobulk.rds")
cat(sprintf("Loading pseudobulk: %s\n", pb_path))
pb <- readRDS(pb_path)
cat(sprintf("pseudobulk dim: %d genes x %d samples\n", nrow(pb), ncol(pb)))

# Filter to samples in target L2
keep <- colData(pb)$parent_L2_joint == L2_target
pb_sub <- pb[, keep]
cat(sprintf("Samples in %s: %d\n", L2_target, ncol(pb_sub)))
if (ncol(pb_sub) == 0) stop("No pseudobulk samples for ", L2_target)

# Restrict to viable NGs
pb_sub <- pb_sub[, colData(pb_sub)$NhoodGroup_renamed %in% lookup$NG_rendered]
cat(sprintf("Samples after viable-NG filter: %d\n", ncol(pb_sub)))

# log-CPM normalize via edgeR
counts <- assay(pb_sub, "counts")
dge <- DGEList(counts = as.matrix(counts))
dge <- calcNormFactors(dge, method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)
cat(sprintf("logCPM matrix: %d genes x %d samples\n", nrow(logcpm), ncol(logcpm)))

# Map gene_union (symbols) to ENSG via markers_df
mr <- rownames(logcpm)
if (mean(grepl("^ENSG", head(mr, 100))) > 0.5) {
  sym2ens <- markers_df %>% distinct(symbol, gene_id)
  gene_ens <- sym2ens$gene_id[match(gene_union, sym2ens$symbol)]
  gene_lookup <- gene_ens[!is.na(gene_ens) & gene_ens %in% mr]
  display_names <- sym2ens$symbol[match(gene_lookup, sym2ens$gene_id)]
} else {
  gene_lookup <- intersect(gene_union, mr)
  display_names <- gene_lookup
}
cat(sprintf("Genes available in pseudobulk: %d / %d\n",
            length(gene_lookup), length(gene_union)))

mat <- logcpm[gene_lookup, , drop = FALSE]
rownames(mat) <- display_names

# Z-score per gene
mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# ---- Column annotation: NG identity, rank, median nhood logFC ----
cd <- as.data.frame(colData(pb_sub))
cd$NG_short <- sub(paste0("^", L2_target, "_"), "", cd$NhoodGroup_renamed)
cd <- cd %>%
  left_join(lookup %>% select(NhoodGroup, size_rank_within_L2, group_med_lfc,
                                name_short),
             by = c("NG_short" = "NhoodGroup")) %>%
  mutate(rank_str = ifelse(size_rank_within_L2 %in% 1:10,
                            as.character(size_rank_within_L2), "10"))

# Order columns by NG rank, then within-NG by donor
ord <- order(cd$size_rank_within_L2, cd$donor)
mat_z <- mat_z[, ord]
cd <- cd[ord, ]

# Order rows by source-NG rank
gene_to_rank <- markers_df %>% group_by(symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(symbol, rank) %>% distinct()
row_rank <- gene_to_rank$rank[match(rownames(mat_z), gene_to_rank$symbol)]
row_rank[is.na(row_rank)] <- max(lookup$size_rank_within_L2) + 1
mat_z <- mat_z[order(row_rank), ]
row_rank <- sort(row_rank)

# ---- Heatmap ----
library(ComplexHeatmap); library(circlize)

ng_levels <- lookup$NhoodGroup
rank_colors <- setNames(RANK_PALETTE[as.character(lookup$size_rank_within_L2)],
                         lookup$NhoodGroup)
lfc_lim <- max(abs(lookup$group_med_lfc), na.rm = TRUE)
lfc_col_fun <- circlize::colorRamp2(
  breaks = c(-lfc_lim, 0, lfc_lim),
  colors = c("#1f77b4", "#f7f7f7", "#d62728"))

col_ann <- HeatmapAnnotation(
  median_logFC = cd$group_med_lfc,
  NhoodGroup   = cd$NG_short,
  col = list(median_logFC = lfc_col_fun, NhoodGroup = rank_colors),
  border = TRUE)

z_lim <- min(max(abs(mat_z), na.rm = TRUE), 3)
z_col_fun <- circlize::colorRamp2(
  breaks = c(-z_lim, 0, z_lim),
  colors = c("#2166AC", "white", "#B2182B"))

ht <- Heatmap(mat_z,
               name = "z(logCPM)",
               col = z_col_fun,
               cluster_rows = FALSE, cluster_columns = FALSE,
               show_row_names = TRUE, show_column_names = FALSE,
               row_names_gp = gpar(fontsize = 6),
               column_split = cd$NG_short,
               row_split = factor(row_rank, levels = sort(unique(row_rank))),
               row_title_gp = gpar(fontsize = 7, fontface = "bold"),
               column_title_gp = gpar(fontsize = 7, fontface = "bold"),
               top_annotation = col_ann,
               border = TRUE,
               use_raster = TRUE, raster_quality = 4,
               heatmap_legend_param = list(title = "z(logCPM)"))

# ---- Render ----
out_dir <- file.path(paths$outputs$reports, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_pdf <- file.path(out_dir,
                      sprintf("panelC_%s_%s_heatmap.pdf", L2_safe, contrast))

pdf(out_pdf,
    width  = max(7, 0.06 * ncol(mat_z) + 4),
    height = max(5, 0.10 * nrow(mat_z) + 2))
draw(ht, merge_legend = TRUE)
dev.off()
cat(sprintf("Wrote %s\n", out_pdf))

write_csv(markers_df %>%
            select(rank, ng, symbol, logFC, FDR, AveExpr) %>%
            arrange(rank, desc(abs(logFC))),
          file.path(out_dir,
                    sprintf("panelC_%s_%s_genes.csv", L2_safe, contrast)))
cat("Done.\n")
