#!/usr/bin/env Rscript
# render_supp_flex_dotplot_l1p5.R
#
# Supplemental Fig 3: FLEX-cohort canonical-marker dotplot at L1.5 resolution.
# Cross-modal companion to fig3e (Xenium nuclear, L1.5) and the cytoplasmic
# variant — same row substrate (L1.5 canonicals + LBridge brackets), same
# column order (l15_order), only the count source differs.
#
# Forked from publication/figures/render/flex/render_fig3d_flex_dotplot_l2s.R:
#   - drops the L2S → LBridge mapping (features CSV already encodes families)
#   - swaps --l2s-labels for L1.5 columns sourced from joint_l1p5.csv filtered
#     to platform == "flex" AND is_artifact == FALSE
#   - swaps --features-csv to fig3_dotplot_features.csv (L1.5 markers)
#   - restores fig3e's 0.42cm row height (L1.5 has ~25 labels, not L2S's ~180)
#
# Drafting status: candidate for submission/figures/fig3/.drafting/panels/
#                  s_flex_dotplot_markers_l1p5/v1/panel.pdf

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(Matrix)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# ---- CLI -------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--features-csv",     required = TRUE,
                    help = "fig3_dotplot_features.csv (L1.5 features; feature, lbridge_family, ...)")
parser$add_argument("--flex-bundle-root", required = TRUE,
                    help = "FLEX integration_intermediate root containing per-compartment scvi_n100/")
parser$add_argument("--joint-l1p5",       required = TRUE,
                    help = "joint_l1p5.csv (cell_id, platform, l1p5_label, is_artifact)")
parser$add_argument("--out-dir",          required = TRUE)
parser$add_argument("--panel-id",         default = "s_flex_dotplot_markers_l1p5")
parser$add_argument("--test", action = "store_true",
                    help = "Subset to 5000 cells per compartment for smoke test")
args <- parser$parse_args()

PANEL_ID <- args$panel_id

# ---- Feature curation ------------------------------------------------------
features <- fread(args$features_csv)
message(sprintf("[%s] %d features × %d LBridge families",
                PANEL_ID, nrow(features),
                length(unique(features$lbridge_family))))

# ---- Load FLEX per-compartment bundles + concatenate -----------------------
load_compartment <- function(comp) {
  bundle <- file.path(args$flex_bundle_root, comp, "scvi_n100")
  cells  <- fread(file.path(bundle, "cells.tsv"), header = FALSE)$V1
  genes  <- fread(file.path(bundle, "genes.tsv"), header = FALSE)$V1
  message(sprintf("[%s]   %s: %d cells × %d genes",
                  PANEL_ID, comp, length(cells), length(genes)))
  mat <- t(as(readMM(file.path(bundle, "counts.mtx.gz")), "CsparseMatrix"))
  dimnames(mat) <- list(genes, cells)
  if (args$test) {
    keep <- seq_len(min(5000L, ncol(mat)))
    mat  <- mat[, keep, drop = FALSE]
    cells <- cells[keep]
    message(sprintf("[%s]   --test: kept %d cells", PANEL_ID, length(cells)))
  }
  list(mat = mat, genes = genes, cells = cells)
}

bundles <- lapply(c("Epithelial", "Stromal", "Immune"), load_compartment)
names(bundles) <- c("Epithelial", "Stromal", "Immune")

all_genes <- Reduce(intersect, lapply(bundles, `[[`, "genes"))
message(sprintf("[%s] %d genes shared across all 3 compartments",
                PANEL_ID, length(all_genes)))

features_present <- as.character(features$feature[features$feature %in% all_genes])
features_missing <- setdiff(as.character(features$feature), all_genes)
if (length(features_missing) > 0) {
  message(sprintf("[%s] %d features missing from FLEX: %s",
                  PANEL_ID, length(features_missing),
                  paste(features_missing, collapse = ", ")))
}
features_kept <- features[feature %in% features_present]
features_kept <- features_kept[!duplicated(feature)]
features_kept[, feature := factor(feature, levels = features_kept$feature)]
features_kept[, lbridge_family := factor(lbridge_family,
                                         levels = unique(features_kept$lbridge_family))]

mats_kept <- lapply(bundles, function(b) {
  b$mat[as.character(features_kept$feature), , drop = FALSE]
})
mat <- do.call(cbind, mats_kept)
message(sprintf("[%s] Concatenated FLEX matrix: %d features × %d cells",
                PANEL_ID, nrow(mat), ncol(mat)))

# ---- L1.5 labels (platform == "flex", artifact-filtered) -------------------
labs <- fread(args$joint_l1p5)[, .(cell_id, platform, l1p5_label, is_artifact)]
labs <- labs[platform == "flex" & is_artifact == FALSE & cell_id %in% colnames(mat)]
mat  <- mat[, labs$cell_id, drop = FALSE]
stopifnot(ncol(mat) == nrow(labs))
message(sprintf("[%s] %d cells × %d features after FLEX/artifact filter",
                PANEL_ID, ncol(mat), nrow(mat)))

# ---- L1.5 column order (epi → str → imm), matching fig3e -------------------
l15_order <- c(
  "BMYO-myo", "LASP", "LASP-basal", "LHS",
  "Fibroblast", "Fibroblast_activated", "Fibroblast_SFRP4",
  "Endothelial", "Vas-capillary", "Lymphatic Endothelial", "Pericyte", "Adipocyte",
  "Macrophage", "cDC", "cDC1", "cDC2", "pDC",
  "Mast cell", "Neutrophil",
  "CD4 T cell", "CD4 Treg", "CD8 T cell", "T-NK",
  "B cell", "Plasma cell"
)
labs[, l1p5_label := factor(l1p5_label,
                            levels = intersect(l15_order, unique(l1p5_label)))]
labs <- labs[!is.na(l1p5_label)]
mat  <- mat[, labs$cell_id, drop = FALSE]
stopifnot(ncol(mat) == nrow(labs))

# ---- Aggregate per (feature × L1.5 label) ----------------------------------
groups <- split(seq_len(ncol(mat)), labs$l1p5_label)
groups <- groups[lengths(groups) > 0]

mat_pct  <- matrix(0, nrow = nrow(mat), ncol = length(groups),
                   dimnames = list(rownames(mat), names(groups)))
mat_mean <- matrix(0, nrow = nrow(mat), ncol = length(groups),
                   dimnames = list(rownames(mat), names(groups)))
for (lab in names(groups)) {
  ix <- groups[[lab]]
  sub <- mat[, ix, drop = FALSE]
  mat_pct[, lab]  <- rowMeans(sub > 0)
  mat_mean[, lab] <- rowMeans(log1p(sub))
}

mat_z <- t(scale(t(mat_mean)))
mat_z[is.na(mat_z)] <- 0
cap <- 3
mat_z[mat_z >  cap] <-  cap
mat_z[mat_z < -cap] <- -cap

# ---- Horizontal orientation: features on x (columns), L1.5 on y (rows) ----
# Transposed from fig3d/fig3e vertical pattern. column_split holds LBridge
# families with a top bracket; row_split is absent (single block of L1.5 labels).
mat_z_h   <- t(mat_z)
mat_pct_h <- t(mat_pct)

col_split <- features_kept$lbridge_family
names(col_split) <- as.character(features_kept$feature)

lbridge_anno <- HeatmapAnnotation(
  lbridge = anno_block(
    gp     = gpar(fill = NA, col = NA),
    labels = NULL
  ),
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  height               = unit(2.5, "cm"),
  which                = "column"
)

# ---- Color + cell function -------------------------------------------------
pal_z <- colorRamp2(c(-cap, 0, cap), c("#3B4992", "#FFFFFF", "#EE0000"))

cell_fun <- function(j, i, x, y, width, height, fill) {
  # Heatmap now has rows = L1.5 labels, cols = features.
  pct <- mat_pct_h[i, j]
  zv  <- mat_z_h[i, j]
  if (is.finite(pct) && pct > 0.02) {
    r_mm <- 2.1 * sqrt(pct)
    grid.circle(
      x  = x, y = y,
      r  = unit(r_mm, "mm"),
      gp = gpar(fill = pal_z(zv), col = NA)
    )
  }
}

ht <- Heatmap(
  mat_z_h,
  name             = "z(log1p)",
  col              = pal_z,
  rect_gp          = gpar(type = "none"),
  cell_fun         = cell_fun,
  cluster_rows     = FALSE, cluster_columns = FALSE,
  show_row_names   = TRUE,  show_column_names = TRUE,
  row_names_gp     = gpar(fontsize = 5),
  column_names_gp  = gpar(fontsize = 5),
  column_names_rot  = 45,
  column_names_side = "bottom",
  column_split      = col_split,
  column_gap       = unit(1.0, "mm"),
  column_title     = NULL,
  top_annotation   = lbridge_anno,
  width            = unit(ncol(mat_z_h) * 0.42, "cm"),
  height           = unit(nrow(mat_z_h) * 0.55, "cm"),
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 5),
    labels_gp     = gpar(fontsize = 4),
    legend_height = unit(2, "cm")
  )
)

decorate_lbridge <- function() {
  fams <- levels(col_split)
  for (i in seq_along(fams)) {
    decorate_annotation("lbridge", slice = i, {
      grid.lines(
        x = c(0.05, 0.95), y = c(0.05, 0.05),
        gp = gpar(col = "black", lwd = 1.0)
      )
      grid.text(
        fams[i],
        x = 0.05, y = 0.10,
        rot = 45, just = c("left", "bottom"),
        gp = gpar(fontsize = 5, fontface = "bold")
      )
    })
  }
}

# ---- Save ------------------------------------------------------------------
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(args$out_dir, paste0(PANEL_ID, ".pdf"))
png_path <- file.path(args$out_dir, paste0(PANEL_ID, ".png"))

# Dims chosen for vertical stacking with the 2 Xenium L1.5 dotplots (composite
# ~16 x 18 in across all 3).
pdf(pdf_path, width = 16, height = 7)
draw(ht)
decorate_lbridge()
dev.off()

png(png_path, width = 16, height = 7, units = "in", res = 600)
draw(ht)
decorate_lbridge()
dev.off()

message(sprintf("[%s] Saved PDF=%s, PNG=%s (16x6 in, horizontal for stacking)",
                PANEL_ID, pdf_path, png_path))
