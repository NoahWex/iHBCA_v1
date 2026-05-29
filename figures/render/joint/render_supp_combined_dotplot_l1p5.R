#!/usr/bin/env Rscript
# render_supp_combined_dotplot_l1p5.R
#
# Combined FLEX + Xenium nuclear L1.5 dotplot: vertical orientation (features
# on y-axis), two panels side by side sharing row axis. Feature set restricted
# to genes present in the Xenium 280-gene panel (the limiting factor).
#
# Pattern source: render_supp_flex_dotplot_l1p5.R + render_supp_xenium_dotplot_nuclear.R
# (same staging directory). Combines both into one ComplexHeatmap ht1 + ht2.

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
parser$add_argument("--features-csv",     required = TRUE)
parser$add_argument("--flex-bundle-root", required = TRUE)
parser$add_argument("--xenium-bundle",    required = TRUE)
parser$add_argument("--joint-l1p5",       required = TRUE)
parser$add_argument("--out-dir",          required = TRUE)
parser$add_argument("--panel-id",         default = "s_combined_dotplot_l1p5")
parser$add_argument("--test", action = "store_true")
args <- parser$parse_args()

PANEL_ID <- args$panel_id

# ---- Feature curation ------------------------------------------------------
features <- fread(args$features_csv)

# Determine which features exist in the Xenium 280-gene panel
xen_genes <- fread(file.path(args$xenium_bundle, "xenium_genes.tsv"),
                   header = FALSE)$V1
message(sprintf("[%s] Xenium panel: %d genes", PANEL_ID, length(xen_genes)))

features_in_panel <- features[feature %in% xen_genes]
features_missing  <- features[!feature %in% xen_genes]
if (nrow(features_missing) > 0) {
  message(sprintf("[%s] %d features not in Xenium panel (dropped): %s",
                  PANEL_ID, nrow(features_missing),
                  paste(features_missing$feature, collapse = ", ")))
}
features_kept <- features_in_panel[!duplicated(feature)]
features_kept[, feature := factor(feature, levels = features_kept$feature)]
features_kept[, lbridge_family := factor(lbridge_family,
                                         levels = unique(features_kept$lbridge_family))]
message(sprintf("[%s] %d shared features × %d LBridge families",
                PANEL_ID, nrow(features_kept),
                length(unique(features_kept$lbridge_family))))

# ---- L1.5 labels (shared) ---------------------------------------------------
labs_all <- fread(args$joint_l1p5)[, .(cell_id, platform, l1p5_label, is_artifact)]
labs_all <- labs_all[is_artifact == FALSE]

l15_order <- c(
  "BMYO-myo", "LASP", "LASP-basal", "LHS",
  "Fibroblast", "Fibroblast_activated", "Fibroblast_SFRP4",
  "Endothelial", "Vas-capillary", "Lymphatic Endothelial", "Pericyte", "Adipocyte",
  "Macrophage", "cDC", "cDC1", "cDC2", "pDC",
  "Mast cell", "Neutrophil",
  "CD4 T cell", "CD4 Treg", "CD8 T cell", "T-NK",
  "B cell", "Plasma cell"
)

# ---- Load FLEX ---------------------------------------------------------------
load_flex_compartment <- function(comp) {
  bundle <- file.path(args$flex_bundle_root, comp, "scvi_n100")
  cells  <- fread(file.path(bundle, "cells.tsv"), header = FALSE)$V1
  genes  <- fread(file.path(bundle, "genes.tsv"), header = FALSE)$V1
  mat <- t(as(readMM(file.path(bundle, "counts.mtx.gz")), "CsparseMatrix"))
  dimnames(mat) <- list(genes, cells)
  if (args$test) {
    keep <- seq_len(min(5000L, ncol(mat)))
    mat  <- mat[, keep, drop = FALSE]
  }
  message(sprintf("[%s] FLEX %s: %d cells × %d genes",
                  PANEL_ID, comp, ncol(mat), nrow(mat)))
  mat
}

flex_mats <- lapply(c("Epithelial", "Stromal", "Immune"), load_flex_compartment)
flex_genes <- Reduce(intersect, lapply(flex_mats, rownames))
feat_vec   <- as.character(features_kept$feature)
feat_flex  <- feat_vec[feat_vec %in% flex_genes]

flex_mats <- lapply(flex_mats, function(m) m[feat_flex, , drop = FALSE])
flex_mat  <- do.call(cbind, flex_mats)

labs_flex <- labs_all[platform == "flex" & cell_id %in% colnames(flex_mat)]
labs_flex[, l1p5_label := factor(l1p5_label,
                                 levels = intersect(l15_order, unique(l1p5_label)))]
labs_flex <- labs_flex[!is.na(l1p5_label)]
flex_mat  <- flex_mat[, labs_flex$cell_id, drop = FALSE]
message(sprintf("[%s] FLEX final: %d cells × %d features",
                PANEL_ID, ncol(flex_mat), nrow(flex_mat)))

# ---- Load Xenium nuclear -----------------------------------------------------
xen_cells <- fread(file.path(args$xenium_bundle, "xenium_cells.tsv"),
                   header = FALSE)$V1
xen_mat   <- t(as(readMM(file.path(args$xenium_bundle,
                                    "xenium_nuclear_counts.mtx.gz")),
                  "CsparseMatrix"))
dimnames(xen_mat) <- list(xen_genes, xen_cells)
feat_xen <- feat_vec[feat_vec %in% xen_genes]
xen_mat  <- xen_mat[feat_xen, , drop = FALSE]

labs_xen <- labs_all[platform == "xenium" & cell_id %in% xen_cells]
labs_xen[, l1p5_label := factor(l1p5_label,
                                levels = intersect(l15_order, unique(l1p5_label)))]
labs_xen <- labs_xen[!is.na(l1p5_label)]
xen_mat  <- xen_mat[, labs_xen$cell_id, drop = FALSE]
message(sprintf("[%s] Xenium nuclear final: %d cells × %d features",
                PANEL_ID, ncol(xen_mat), nrow(xen_mat)))

# ---- Shared feature order (intersection, preserving CSV order) ---------------
shared_features <- intersect(feat_flex, feat_xen)
features_final  <- features_kept[feature %in% shared_features]
features_final[, feature := factor(feature, levels = shared_features)]
message(sprintf("[%s] Shared features: %d", PANEL_ID, length(shared_features)))

flex_mat <- flex_mat[shared_features, , drop = FALSE]
xen_mat  <- xen_mat[shared_features, , drop = FALSE]

# ---- Aggregate per (feature × L1.5 label) ------------------------------------
aggregate_dotplot <- function(mat, labs_dt) {
  groups <- split(seq_len(ncol(mat)), labs_dt$l1p5_label)
  groups <- groups[lengths(groups) > 0]
  pct  <- matrix(0, nrow = nrow(mat), ncol = length(groups),
                 dimnames = list(rownames(mat), names(groups)))
  mn   <- matrix(0, nrow = nrow(mat), ncol = length(groups),
                 dimnames = list(rownames(mat), names(groups)))
  for (lab in names(groups)) {
    ix <- groups[[lab]]
    sub <- mat[, ix, drop = FALSE]
    pct[, lab] <- rowMeans(sub > 0)
    mn[, lab]  <- rowMeans(log1p(sub))
  }
  z <- t(scale(t(mn)))
  z[is.na(z)] <- 0
  cap <- 3
  z[z >  cap] <-  cap
  z[z < -cap] <- -cap
  list(z = z, pct = pct)
}

flex_agg <- aggregate_dotplot(flex_mat, labs_flex)
xen_agg  <- aggregate_dotplot(xen_mat,  labs_xen)

# ---- Heatmap construction (vertical: features = rows, L1.5 = columns) -------
row_split <- features_final$lbridge_family
names(row_split) <- as.character(features_final$feature)

pal_z <- colorRamp2(c(-3, 0, 3), c("#3B4992", "#FFFFFF", "#EE0000"))

make_cell_fun <- function(pct_mat, z_mat) {
  force(pct_mat); force(z_mat)
  function(j, i, x, y, width, height, fill) {
    pct <- pct_mat[i, j]
    zv  <- z_mat[i, j]
    if (is.finite(pct) && pct > 0.02) {
      r_mm <- 2.1 * sqrt(pct)
      grid.circle(x = x, y = y, r = unit(r_mm, "mm"),
                  gp = gpar(fill = pal_z(zv), col = NA))
    }
  }
}

lbridge_anno <- rowAnnotation(
  lbridge = anno_block(
    gp     = gpar(fill = NA, col = NA),
    labels = NULL
  ),
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  width                = unit(2.5, "cm")
)

ht_flex <- Heatmap(
  flex_agg$z,
  name              = "z(log1p)",
  col               = pal_z,
  rect_gp           = gpar(type = "none"),
  cell_fun          = make_cell_fun(flex_agg$pct, flex_agg$z),
  cluster_rows      = FALSE, cluster_columns = FALSE,
  show_row_names    = FALSE,
  show_column_names = TRUE,
  column_names_gp   = gpar(fontsize = 5),
  column_names_rot  = 45,
  column_names_side = "bottom",
  column_title      = "FLEX",
  column_title_side = "top",
  column_title_gp   = gpar(fontsize = 7, fontface = "bold"),
  row_split         = row_split,
  row_gap           = unit(1.0, "mm"),
  row_title         = NULL,
  left_annotation   = lbridge_anno,
  width             = unit(ncol(flex_agg$z) * 0.42, "cm"),
  height            = unit(nrow(flex_agg$z) * 0.42, "cm"),
  show_heatmap_legend = FALSE
)

ht_xen <- Heatmap(
  xen_agg$z,
  name              = "z(log1p)_xen",
  col               = pal_z,
  rect_gp           = gpar(type = "none"),
  cell_fun          = make_cell_fun(xen_agg$pct, xen_agg$z),
  cluster_rows      = FALSE, cluster_columns = FALSE,
  show_row_names    = TRUE,
  row_names_side    = "right",
  row_names_gp      = gpar(fontsize = 5),
  show_column_names = TRUE,
  column_names_gp   = gpar(fontsize = 5),
  column_names_rot  = 45,
  column_names_side = "bottom",
  column_title      = "Xenium (nuclear)",
  column_title_side = "top",
  column_title_gp   = gpar(fontsize = 7, fontface = "bold"),
  row_split         = row_split,
  row_gap           = unit(1.0, "mm"),
  row_title         = NULL,
  width             = unit(ncol(xen_agg$z) * 0.42, "cm"),
  height            = unit(nrow(xen_agg$z) * 0.42, "cm"),
  heatmap_legend_param = list(
    title_gp      = gpar(fontsize = 5),
    labels_gp     = gpar(fontsize = 4),
    legend_height = unit(2, "cm")
  )
)

decorate_lbridge <- function() {
  fams <- levels(row_split)
  for (i in seq_along(fams)) {
    decorate_annotation("lbridge", slice = i, {
      grid.lines(
        x = c(0.95, 0.95), y = c(0.05, 0.95),
        gp = gpar(col = "black", lwd = 1.0)
      )
      grid.text(
        fams[i],
        x = 0.90, y = 0.5,
        rot = 0, just = c("right", "center"),
        gp = gpar(fontsize = 5, fontface = "bold")
      )
    })
  }
}

# ---- Save --------------------------------------------------------------------
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(args$out_dir, paste0(PANEL_ID, ".pdf"))
png_path <- file.path(args$out_dir, paste0(PANEL_ID, ".png"))

ht_list <- ht_flex + ht_xen

pdf(pdf_path, width = 14, height = 16)
draw(ht_list, gap = unit(4, "mm"))
decorate_lbridge()
dev.off()

png(png_path, width = 14, height = 16, units = "in", res = 600)
draw(ht_list, gap = unit(4, "mm"))
decorate_lbridge()
dev.off()

message(sprintf("[%s] Saved: %s, %s", PANEL_ID, pdf_path, png_path))
