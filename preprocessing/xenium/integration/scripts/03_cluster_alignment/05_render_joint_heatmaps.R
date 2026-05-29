#!/usr/bin/env Rscript
# Render per-source joint-cluster marker heatmaps from the markers tables
# written by 04_joint_markers.R. Produces one heatmap per marker source:
#
#   joint_markers_heatmap_nuclear.{pdf,png}    — Xenium nuclear panel markers
#   joint_markers_heatmap_flex_panel.{pdf,png} — FLEX on the Xenium panel
#   joint_markers_heatmap_flex_full.{pdf,png}  — FLEX full-transcriptome markers
#
# Each heatmap: (top_n markers per cluster) × (all clusters), rows are genes
# ordered by first-cluster-of-maximum-expression, columns are clusters ordered
# numerically. Fill is standardized (z-score) mean expression per cluster.
#
# Usage:
#   Rscript 05_render_joint_heatmaps.R \
#       --joint-leiden  <joint_leiden_assignments.csv> \
#       --joint-obs     <joint_obs.csv> \
#       --leiden-res    1.0 \
#       --markers-dir   <pipeline/outputs/annotations> \
#       --nuc-bundle    <pooled_nuclear/> \
#       --flex-dir      <integration_intermediate/> \
#       --xenium-panel  <xenium_panel_genes.txt> \
#       --out-dir       <pipeline/outputs/previews> \
#       [--top-n 5]

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(ggplot2)
  library(optparse)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--joint-leiden", type = "character"),
  make_option("--joint-obs",    type = "character"),
  make_option("--leiden-res",   type = "character", default = "1.0"),
  make_option("--markers-dir",  type = "character"),
  make_option("--nuc-bundle",   type = "character"),
  make_option("--flex-dir",     type = "character"),
  make_option("--xenium-panel", type = "character"),
  make_option("--out-dir",      type = "character"),
  make_option("--top-n",        type = "integer",   default = 5)
)))

dir.create(opts$`out-dir`, showWarnings = FALSE, recursive = TRUE)
leiden_col <- paste0("leiden_", opts$`leiden-res`)

leiden <- fread(opts$`joint-leiden`)
obs    <- fread(opts$`joint-obs`)
setkey(leiden, cell_id); setkey(obs, cell_id)
jmeta  <- leiden[obs, nomatch = 0][!is.na(get(leiden_col))]

panel_genes <- trimws(readLines(opts$`xenium-panel`))
panel_genes <- panel_genes[nchar(panel_genes) > 0]

read_mtx_bundle <- function(dir, mtx_name, genes_name, cells_name) {
  mtx_path <- file.path(dir, mtx_name)
  if (!file.exists(mtx_path)) mtx_path <- paste0(mtx_path, ".gz")
  list(
    mat   = t(as(readMM(mtx_path), "CsparseMatrix")),
    cells = readLines(file.path(dir, cells_name)),
    genes = readLines(file.path(dir, genes_name))
  )
}

render_heatmap <- function(mat, cells, genes, meta, markers, title, out_stem) {
  # Keep only cells with a cluster assignment
  keep <- cells %in% meta$cell_id
  mat <- mat[, keep, drop = FALSE]; cells <- cells[keep]
  rownames(mat) <- genes; colnames(mat) <- cells
  m2 <- meta[match(cells, cell_id)]
  y  <- as.character(m2[[leiden_col]])

  # Top-N markers per cluster by AUC, unique genes preserved
  top <- markers[order(-auc), head(.SD, opts$`top-n`), by = group]
  top_genes <- unique(top$feature)
  top_genes <- intersect(top_genes, rownames(mat))

  # Mean expression per gene per cluster; z-score per gene across clusters.
  # Order clusters numerically when possible (avoids "1, 10, 11, ..., 2" lex sort).
  clusters <- unique(y)
  n <- suppressWarnings(as.integer(clusters))
  clusters <- if (!any(is.na(n))) clusters[order(n)] else sort(clusters)
  M <- sapply(clusters, function(k) {
    idx <- which(y == k)
    Matrix::rowMeans(mat[top_genes, idx, drop = FALSE])
  })
  colnames(M) <- clusters
  z <- t(scale(t(M)))

  # Order genes by argmax cluster (assigned cluster of peak expression)
  peak <- apply(z, 1, which.max)
  gene_order <- top_genes[order(peak, -apply(z, 1, max))]
  z <- z[gene_order, , drop = FALSE]

  df <- as.data.table(as.table(z))
  setnames(df, c("gene", "cluster", "z"))
  df[, cluster := factor(cluster, levels = clusters)]
  df[, gene    := factor(gene,    levels = gene_order)]

  p <- ggplot(df, aes(cluster, gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-3, 3), oob = scales::squish,
                         name = "z(mean\nper cluster)") +
    labs(title = title, x = leiden_col, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          panel.grid = element_blank(),
          plot.title = element_text(size = 10))

  h <- max(4, 0.14 * length(gene_order))
  w <- max(6, 0.35 * length(clusters) + 2)
  ggsave(paste0(out_stem, ".pdf"), p, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(out_stem, ".png"), p, width = w, height = h, dpi = 200)
  message("  wrote ", out_stem, ".{pdf,png}")
}

# ---------------------------------------------------------------------------
message("=== [1/3] Xenium nuclear heatmap ===")
markers_nuc <- fread(file.path(opts$`markers-dir`, "joint_markers_nuclear.csv"))
nuc <- read_mtx_bundle(
  opts$`nuc-bundle`,
  "xenium_nuclear_counts.mtx.gz", "xenium_genes.tsv", "xenium_cells.tsv"
)
render_heatmap(nuc$mat, nuc$cells, nuc$genes,
               jmeta[platform == "xenium"], markers_nuc,
               "Joint-cluster markers — Xenium nuclear (280-gene panel)",
               file.path(opts$`out-dir`, "joint_markers_heatmap_nuclear"))
rm(nuc); gc()

# ---------------------------------------------------------------------------
message("=== [2/3] FLEX on Xenium panel heatmap ===")
markers_fp <- fread(file.path(opts$`markers-dir`, "joint_markers_flex_panel.csv"))
flex_compartments <- c("Immune", "Epithelial", "Stromal")
flex_parts <- lapply(flex_compartments, function(comp) {
  cdir <- file.path(opts$`flex-dir`, comp, "scvi_n100")
  read_mtx_bundle(cdir, "counts.mtx.gz", "genes.tsv", "cells.tsv")
})
ref_genes <- flex_parts[[1]]$genes
flex <- list(
  mat   = do.call(cbind, lapply(flex_parts, `[[`, "mat")),
  cells = do.call(c,     lapply(flex_parts, `[[`, "cells")),
  genes = ref_genes
)
rm(flex_parts); gc()
panel_idx <- which(flex$genes %in% panel_genes)
render_heatmap(flex$mat[panel_idx, , drop = FALSE],
               flex$cells, flex$genes[panel_idx],
               jmeta[platform == "flex"], markers_fp,
               "Joint-cluster markers — FLEX on Xenium panel",
               file.path(opts$`out-dir`, "joint_markers_heatmap_flex_panel"))

# ---------------------------------------------------------------------------
message("=== [3/3] FLEX full transcriptome heatmap ===")
markers_ff <- fread(file.path(opts$`markers-dir`, "joint_markers_flex_full.csv"))
render_heatmap(flex$mat, flex$cells, flex$genes,
               jmeta[platform == "flex"], markers_ff,
               "Joint-cluster markers — FLEX full transcriptome",
               file.path(opts$`out-dir`, "joint_markers_heatmap_flex_full"))

message("=== done ===")
