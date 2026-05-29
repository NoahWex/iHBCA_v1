#!/usr/bin/env Rscript
# Compute joint-cluster markers for three expression sources using
# presto::wilcoxauc:
#   (1) xenium_nuclear    — Xenium nuclear panel, Xenium cells only (280 genes)
#   (2) flex_panel        — FLEX expression restricted to the Xenium panel (280 genes)
#   (3) flex_full         — FLEX full transcriptome, FLEX cells only (~18K genes)
#
# Groups are joint Leiden clusters at --leiden-res (produced by 03_joint_leiden.py).
# Writes one markers CSV per source, consumed by 05_render_joint_heatmaps.R.
#
# Usage:
#   Rscript 04_joint_markers.R \
#       --joint-leiden   <joint_leiden_assignments.csv> \
#       --joint-obs      <joint_obs.csv> \
#       --leiden-res     1.0 \
#       --nuc-bundle     <pooled_nuclear/>  (xenium_nuclear_counts.mtx.gz + genes/cells tsv) \
#       --flex-dir       <integration_intermediate/>  (FLEX MTX + genes + cells) \
#       --xenium-panel   <xenium_panel_genes.txt> \
#       --out-dir        <pipeline/outputs/annotations/>

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(presto)
  library(optparse)
})

opts <- parse_args(OptionParser(option_list = list(
  make_option("--joint-leiden", type = "character"),
  make_option("--joint-obs",    type = "character"),
  make_option("--leiden-res",   type = "character", default = "1.0"),
  make_option("--nuc-bundle",   type = "character"),
  make_option("--flex-dir",     type = "character"),
  make_option("--xenium-panel", type = "character"),
  make_option("--out-dir",      type = "character")
)))

dir.create(opts$`out-dir`, showWarnings = FALSE, recursive = TRUE)
leiden_col <- paste0("leiden_", opts$`leiden-res`)

message("=== loading joint Leiden + obs ===")
leiden <- fread(opts$`joint-leiden`)
obs    <- fread(opts$`joint-obs`)
setkey(leiden, cell_id); setkey(obs, cell_id)
jmeta  <- leiden[obs, nomatch = 0]
jmeta  <- jmeta[!is.na(get(leiden_col))]
message("cells with leiden assignment: ", nrow(jmeta))
message("clusters at ", leiden_col, ": ", jmeta[, uniqueN(get(leiden_col))])

panel_genes <- readLines(opts$`xenium-panel`)
panel_genes <- trimws(panel_genes[nchar(panel_genes) > 0])
message("Xenium panel genes: ", length(panel_genes))

read_mtx_bundle <- function(dir, mtx_name, genes_name, cells_name) {
  # Returns list(mat = sparse cells x genes, cells = chr, genes = chr).
  mtx_path <- file.path(dir, mtx_name)
  if (!file.exists(mtx_path)) mtx_path <- paste0(mtx_path, ".gz")
  mat <- as(readMM(mtx_path), "CsparseMatrix")
  cells <- readLines(file.path(dir, cells_name))
  genes <- readLines(file.path(dir, genes_name))
  # presto wants genes x cells for wilcoxauc(X, y) with y over columns
  list(mat = t(mat), cells = cells, genes = genes)
}

run_presto <- function(mat, cells, genes, meta, label, out_path) {
  # mat: genes x cells (sparse); meta$cell_id and meta[[leiden_col]] align
  keep <- cells %in% meta$cell_id
  mat <- mat[, keep, drop = FALSE]; cells <- cells[keep]
  rownames(mat) <- genes; colnames(mat) <- cells
  m2 <- meta[match(cells, cell_id)]
  y <- as.character(m2[[leiden_col]])
  message(sprintf("  [%s] cells=%d genes=%d groups=%d",
                  label, ncol(mat), nrow(mat), length(unique(y))))
  res <- wilcoxauc(mat, y)
  fwrite(res, out_path)
  message("  wrote ", out_path)
}

# ---------------------------------------------------------------------------
# Source 1: Xenium nuclear counts (Xenium cells, 280-gene panel)
# ---------------------------------------------------------------------------
message("=== [1/3] Xenium nuclear ===")
nuc <- read_mtx_bundle(
  opts$`nuc-bundle`,
  "xenium_nuclear_counts.mtx.gz", "xenium_genes.tsv", "xenium_cells.tsv"
)
run_presto(nuc$mat, nuc$cells, nuc$genes,
           jmeta[platform == "xenium"],
           "xenium_nuclear",
           file.path(opts$`out-dir`, "joint_markers_nuclear.csv"))
rm(nuc); gc()

# ---------------------------------------------------------------------------
# Source 2 + 3: FLEX (on panel, and full transcriptome)
#   FLEX is stored as three compartment bundles under --flex-dir:
#     {flex-dir}/{Immune,Epithelial,Stromal}/scvi_n100/{counts.mtx.gz,genes.tsv,cells.tsv}
#   Genes are identical across compartments (full transcriptome from one
#   upstream object), so we cbind on cells.
# ---------------------------------------------------------------------------
message("=== [2+3/3] FLEX ===")
flex_compartments <- c("Immune", "Epithelial", "Stromal")
flex_parts <- lapply(flex_compartments, function(comp) {
  cdir <- file.path(opts$`flex-dir`, comp, "scvi_n100")
  message("  loading FLEX/", comp)
  read_mtx_bundle(cdir, "counts.mtx.gz", "genes.tsv", "cells.tsv")
})
ref_genes <- flex_parts[[1]]$genes
for (i in seq_along(flex_parts)[-1]) {
  stopifnot(identical(flex_parts[[i]]$genes, ref_genes))
}
flex <- list(
  mat   = do.call(cbind, lapply(flex_parts, `[[`, "mat")),
  cells = do.call(c,     lapply(flex_parts, `[[`, "cells")),
  genes = ref_genes
)
rm(flex_parts); gc()
message("  FLEX combined: cells=", length(flex$cells), " genes=", length(flex$genes))
# Source 2: restrict to Xenium panel
panel_idx <- which(flex$genes %in% panel_genes)
message("FLEX genes matching Xenium panel: ", length(panel_idx), " / ", length(panel_genes))
run_presto(flex$mat[panel_idx, , drop = FALSE],
           flex$cells,
           flex$genes[panel_idx],
           jmeta[platform == "flex"],
           "flex_panel",
           file.path(opts$`out-dir`, "joint_markers_flex_panel.csv"))

# Source 3: full FLEX transcriptome
run_presto(flex$mat, flex$cells, flex$genes,
           jmeta[platform == "flex"],
           "flex_full",
           file.path(opts$`out-dir`, "joint_markers_flex_full.csv"))

message("=== done ===")
