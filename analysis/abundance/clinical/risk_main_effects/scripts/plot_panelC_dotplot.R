#!/usr/bin/env Rscript
# plot_panelC_dotplot.R
#
# Per-L2 dotplot: gene (rows) x NhoodGroup (cols).
#   color = mean log-CPM per NG (across donor pseudobulks)
#   size  = fraction of donor pseudobulks with non-zero counts in that NG
# Fallback / complement to plot_panelC_marker_heatmap.R when ComplexHeatmap is
# heavy or when a Seurat::DotPlot-style aesthetic is preferred. Donor-grain.

suppressPackageStartupMessages({
  library(SummarizedExperiment); library(SingleCellExperiment); library(edgeR)
  library(Matrix); library(dplyr); library(readr); library(yaml); library(tidyr)
  library(ggplot2); library(scales)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- which(args == name); if (length(i) == 0) return(default); args[i + 1]
}
L2_target  <- get_arg("--L2",  "epi::BMYO-basal")
contrast   <- get_arg("--contrast", "parity_x_HR_BRCA1")
top_n      <- as.integer(get_arg("--top-n", "10"))

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_panelC_dotplot.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))
cfg <- parse_args_inquiry()
paths <- cfg$paths

cat(sprintf("=== Panel C dotplot — L2='%s' contrast='%s' top_n=%d ===\n",
            L2_target, contrast, top_n))

L2_compartment <- strsplit(L2_target, "::", fixed = TRUE)[[1]][1]
L2_label       <- strsplit(L2_target, "::", fixed = TRUE)[[1]][2]
L2_safe        <- gsub("::", "__", gsub("[^A-Za-z0-9_]", "_", L2_target))

RANK_PALETTE <- c(
  "1" = "#332288", "2" = "#117733", "3" = "#44AA99", "4" = "#88CCEE",
  "5" = "#DDCC77", "6" = "#CC6677", "7" = "#AA4499", "8" = "#882255",
  "9" = "#999933", "10" = "#661100"
)

lookup <- read_csv(file.path(paths$inquiry_root, "outputs", "nhoodgroup_renaming",
                              "lookup.csv"), show_col_types = FALSE) %>%
  filter(viable, contrast == !!contrast, parent_L2_joint == L2_target) %>%
  mutate(NhoodGroup = as.character(NhoodGroup),
         NG_rendered = paste0(parent_L2_joint, "_", NhoodGroup)) %>%
  arrange(size_rank_within_L2)
if (nrow(lookup) == 0) stop("No viable NGs for ", L2_target)

# Top markers per NG
marker_dir <- file.path(paths$outputs$stageF3, contrast)
top_markers <- list()
for (i in seq_len(nrow(lookup))) {
  ng <- lookup$NhoodGroup[i]; rk <- lookup$size_rank_within_L2[i]
  fname <- sprintf("%s__%s_%s_vs_parent.csv",
                    L2_compartment, gsub("-", "_", L2_label), ng)
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
cat(sprintf("Marker union: %d genes\n", length(gene_union)))

# Pseudobulk
pb_path <- file.path(paths$outputs$stageF2, contrast, "pseudobulk.rds")
pb <- readRDS(pb_path)
keep <- colData(pb)$parent_L2_joint == L2_target &
         colData(pb)$NhoodGroup_renamed %in% lookup$NG_rendered
pb_sub <- pb[, keep]
cat(sprintf("Pseudobulk samples: %d\n", ncol(pb_sub)))

counts <- as.matrix(assay(pb_sub, "counts"))
dge <- DGEList(counts = counts) %>% calcNormFactors(method = "TMM")
logcpm <- cpm(dge, log = TRUE, prior.count = 1)

# Map symbols -> ENSG via markers_df
mr <- rownames(logcpm)
sym2ens <- markers_df %>% distinct(symbol, gene_id)
gene_ens <- sym2ens$gene_id[match(gene_union, sym2ens$symbol)]
keep_mask <- !is.na(gene_ens) & gene_ens %in% mr
gene_lookup <- gene_ens[keep_mask]
display_names <- sym2ens$symbol[match(gene_lookup, sym2ens$gene_id)]
mat <- logcpm[gene_lookup, , drop = FALSE]
counts_mat <- counts[gene_lookup, , drop = FALSE]
rownames(mat) <- display_names
rownames(counts_mat) <- display_names

# Per (gene x NG): mean log-CPM, fraction donors non-zero
cd <- as.data.frame(colData(pb_sub))
cd$NG_short <- sub(paste0("^", L2_target, "_"), "", cd$NhoodGroup_renamed)

dot_long <- list()
for (ng in unique(cd$NG_short)) {
  ng_cols <- which(cd$NG_short == ng)
  if (length(ng_cols) == 0) next
  mean_lcpm <- rowMeans(mat[, ng_cols, drop = FALSE])
  pct_pos   <- rowMeans(counts_mat[, ng_cols, drop = FALSE] > 0)
  dot_long[[ng]] <- tibble(gene = rownames(mat),
                            NG = ng,
                            mean_lcpm = mean_lcpm,
                            pct_pos = pct_pos,
                            n_donors = length(ng_cols))
}
dot <- bind_rows(dot_long) %>%
  left_join(lookup %>% select(NG = NhoodGroup, size_rank_within_L2,
                                group_med_lfc, n_nhoods_in_group),
             by = "NG") %>%
  mutate(rank_str = ifelse(size_rank_within_L2 %in% 1:10,
                            as.character(size_rank_within_L2), "10"))

# Order columns by rank, rows by source-NG rank
dot$NG <- factor(dot$NG, levels = lookup$NhoodGroup)
gene_to_rank <- markers_df %>% group_by(symbol) %>%
  slice_min(adj.P.Val, n = 1, with_ties = FALSE) %>%
  ungroup() %>% select(symbol, rank) %>% distinct()
gene_order <- gene_to_rank %>%
  filter(symbol %in% display_names) %>%
  arrange(rank, symbol) %>% pull(symbol)
dot$gene <- factor(dot$gene, levels = rev(gene_order))

# z-score the mean_lcpm across NGs per gene for color (centers expression)
dot <- dot %>% group_by(gene) %>%
  mutate(z_lcpm = scale(mean_lcpm)[, 1]) %>%
  ungroup() %>%
  mutate(z_lcpm = pmin(pmax(z_lcpm, -2), 2))

p <- ggplot(dot, aes(x = NG, y = gene)) +
  geom_point(aes(color = z_lcpm, size = pct_pos), shape = 16) +
  scale_color_gradient2(low = "#2166AC", mid = "#f7f7f7", high = "#B2182B",
                        midpoint = 0,
                        name = "z(mean logCPM)\nacross NGs") +
  scale_size_continuous("frac donors\nnon-zero",
                         range = c(0.2, 4), limits = c(0, 1)) +
  labs(x = sprintf("NhoodGroup (%s)", L2_target),
        y = NULL,
        caption = sprintf("contrast=%s, %d donor pseudobulks",
                           contrast, ncol(pb_sub))) +
  theme_classic(base_size = 8) +
  theme(panel.grid.major = element_line(color = "grey90", linewidth = 0.15),
         axis.text.y = element_text(size = 6),
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
         plot.caption = element_text(size = 6.5, color = "grey30"),
         legend.key.size = unit(0.3, "cm"))

out_dir <- file.path(paths$outputs$reports, "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
out_pdf <- file.path(out_dir,
                      sprintf("panelC_%s_%s_dotplot.pdf", L2_safe, contrast))
ggsave(out_pdf, p,
        width  = max(4, 0.5 * length(unique(dot$NG)) + 2),
        height = max(4, 0.13 * length(unique(dot$gene)) + 1.5),
        limitsize = FALSE)
cat(sprintf("Wrote %s\n", out_pdf))
