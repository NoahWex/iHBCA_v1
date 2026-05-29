#!/usr/bin/env Rscript
# plot_panelC_leadingedge_heatmap.R
#
# Donor-pseudobulk × leading-edge-gene heatmap. Per L2 in a single contrast:
#   - For each viable NhoodGroup, pick top sig pathways (Hallmark + Reactome,
#     padj < 0.10) by |NES|.
#   - Pull leading-edge gene unions per NG across those pathways.
#   - Render heatmap of (genes × donor-pseudobulk samples) grouped by NG,
#     row-split by source NG.
# Leading-edge genes are persuasive markers because they cluster on a pathway
# AND drive the enrichment signal.
#
# Args:
#   --L2, --contrast, --inquiry-dir
#   --top-paths-per-ng (default 5)  pathways per NG by |NES|
#   --max-genes-per-ng (default 25) leading-edge gene cap per NG

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(SingleCellExperiment); library(edgeR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(tidyr)
  library(ggplot2); library(scales); library(stringr)
})
`%||%` <- function(a, b) if (is.null(a)) b else a
have_ch <- requireNamespace("ComplexHeatmap", quietly = TRUE)
have_circlize <- requireNamespace("circlize", quietly = TRUE)
if (!have_ch || !have_circlize) stop("ComplexHeatmap + circlize required.")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == name); if (length(i) == 0) return(default); args[i + 1]
}
L2_target  <- get_arg("--L2",  "epi::LHS-major")
contrast   <- get_arg("--contrast", "parity_x_HR_BRCA1")
top_paths  <- as.integer(get_arg("--top-paths-per-ng", "5"))
max_genes  <- as.integer(get_arg("--max-genes-per-ng", "25"))

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_panelC_leadingedge_heatmap.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
cfg <- parse_args_inquiry()
paths <- cfg$paths

cat(sprintf("=== Panel C leading-edge — L2='%s' contrast='%s' top_paths=%d max_genes=%d ===\n",
            L2_target, contrast, top_paths, max_genes))

L2_compartment <- strsplit(L2_target, "::", fixed = TRUE)[[1]][1]
L2_label       <- strsplit(L2_target, "::", fixed = TRUE)[[1]][2]
L2_safe        <- gsub("::", "__", gsub("[^A-Za-z0-9_]", "_", L2_target))

RANK_PALETTE <- c(
  "1" = "#332288", "2" = "#117733", "3" = "#44AA99", "4" = "#88CCEE",
  "5" = "#DDCC77", "6" = "#CC6677", "7" = "#AA4499", "8" = "#882255",
  "9" = "#999933", "10" = "#661100"
)

# ---- Lookup ----
lookup <- read_csv(file.path(paths$inquiry_root, "outputs", "nhoodgroup_renaming",
                              "lookup.csv"), show_col_types = FALSE) %>%
  filter(viable, contrast == !!contrast, parent_L2_joint == L2_target) %>%
  mutate(NhoodGroup = as.character(NhoodGroup),
         NG_rendered = paste0(parent_L2_joint, "_", NhoodGroup),
         NG_underscored = gsub("::|-", "_", NG_rendered)) %>%
  arrange(size_rank_within_L2)
if (nrow(lookup) == 0) stop("No viable NGs for ", L2_target, " in ", contrast)
cat(sprintf("Viable NGs: %d\n", nrow(lookup)))

# ---- GSEA: top pathways per NG, leading-edge gene unions ----
gsea_path <- file.path(paths$inquiry_root, "outputs", "stageI3_gsea",
                        "per_nhoodgroup_nes.csv")
gsea <- read_csv(gsea_path, show_col_types = FALSE) %>%
  filter(contrast == !!contrast,
         NhoodGroup_renamed %in% lookup$NG_underscored,
         padj < 0.10)
cat(sprintf("GSEA sig (padj<0.10) for L2: %d rows across %d pathways\n",
            nrow(gsea), n_distinct(gsea$pathway)))

ng_genes <- list()
ng_paths <- list()
for (i in seq_len(nrow(lookup))) {
  ng_under <- lookup$NG_underscored[i]
  ng_id    <- lookup$NhoodGroup[i]
  rk       <- lookup$size_rank_within_L2[i]
  g <- gsea %>%
    filter(NhoodGroup_renamed == ng_under) %>%
    arrange(desc(abs(NES))) %>%
    head(top_paths)
  if (nrow(g) == 0) next
  ng_paths[[ng_id]] <- g %>%
    transmute(NG = ng_id, rank = rk, pathway, NES, padj,
              direction = ifelse(NES > 0, "up", "dn"))
  le <- unlist(strsplit(g$leadingEdge, ";"))
  le <- unique(le[!is.na(le) & le != ""])
  if (length(le) > max_genes) le <- le[seq_len(max_genes)]
  ng_genes[[ng_id]] <- tibble(NG = ng_id, rank = rk, gene_id = le)
  cat(sprintf("  %s (rank %d): %d sig pathways, %d leading-edge genes\n",
              ng_id, rk, nrow(g), length(le)))
}
genes_df <- bind_rows(ng_genes)
gene_union <- unique(genes_df$gene_id)
cat(sprintf("Leading-edge gene union: %d ENSG ids\n", length(gene_union)))

# ENSG -> symbol via gene_data.csv (from assembly outputs)
gene_data_path <- file.path(paths$inputs$components_root, "gene_data.csv")
sym_map <- if (file.exists(gene_data_path)) {
  gd <- read_csv(gene_data_path, show_col_types = FALSE)
  if ("symbol" %in% colnames(gd)) {
    setNames(gd$symbol, gd[[1]])
  } else NULL
} else NULL

# ---- Pseudobulk + log-CPM ----
pb_path <- file.path(paths$outputs$stageF2, contrast, "pseudobulk.rds")
pb <- readRDS(pb_path)
keep <- colData(pb)$parent_L2_joint == L2_target &
         colData(pb)$NhoodGroup_renamed %in% lookup$NG_rendered
pb_sub <- pb[, keep]
cat(sprintf("Samples in L2: %d (after viable filter)\n", ncol(pb_sub)))
counts <- as.matrix(assay(pb_sub, "counts"))
dge <- DGEList(counts = counts) %>% calcNormFactors(method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

mr <- rownames(logcpm)
gene_lookup <- intersect(gene_union, mr)
display_names <- if (!is.null(sym_map)) {
  s <- sym_map[gene_lookup]; ifelse(is.na(s) | s == "", gene_lookup, s)
} else gene_lookup
mat <- logcpm[gene_lookup, , drop = FALSE]
rownames(mat) <- display_names
cat(sprintf("Genes available in pseudobulk: %d / %d\n",
            length(gene_lookup), length(gene_union)))

mat_z <- t(scale(t(mat)))
mat_z[is.na(mat_z)] <- 0

# ---- Column annotation ----
cd <- as.data.frame(colData(pb_sub))
cd$NG_short <- sub(paste0("^", L2_target, "_"), "", cd$NhoodGroup_renamed)
cd <- cd %>%
  left_join(lookup %>% select(NhoodGroup, size_rank_within_L2, group_med_lfc),
             by = c("NG_short" = "NhoodGroup"))

ord <- order(cd$size_rank_within_L2, cd$donor)
mat_z <- mat_z[, ord]
cd <- cd[ord, ]

# Row order by source-NG rank
genes_df <- genes_df %>% mutate(symbol = display_names[match(gene_id, gene_lookup)])
gene_to_rank <- genes_df %>%
  group_by(symbol) %>% slice_min(rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(symbol, rank) %>% distinct()
row_rank <- gene_to_rank$rank[match(rownames(mat_z), gene_to_rank$symbol)]
row_rank[is.na(row_rank)] <- max(lookup$size_rank_within_L2) + 1
mat_z <- mat_z[order(row_rank), ]
row_rank <- sort(row_rank)

library(ComplexHeatmap); library(circlize)
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
z_col_fun <- circlize::colorRamp2(c(-z_lim, 0, z_lim),
                                   c("#2166AC", "white", "#B2182B"))

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

out_dir <- file.path(paths$outputs$reports, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_pdf <- file.path(out_dir,
                      sprintf("panelC_%s_%s_leadingedge.pdf", L2_safe, contrast))
pdf(out_pdf,
    width  = max(7, 0.06 * ncol(mat_z) + 4),
    height = max(5, 0.10 * nrow(mat_z) + 2))
draw(ht, merge_legend = TRUE)
dev.off()
cat(sprintf("Wrote %s\n", out_pdf))

# Sidecar: pathway + LE gene records
write_csv(bind_rows(ng_paths) %>% arrange(rank, desc(abs(NES))),
          file.path(out_dir,
                    sprintf("panelC_%s_%s_leadingedge_pathways.csv",
                            L2_safe, contrast)))
write_csv(genes_df,
          file.path(out_dir,
                    sprintf("panelC_%s_%s_leadingedge_genes.csv",
                            L2_safe, contrast)))
cat("Done.\n")
