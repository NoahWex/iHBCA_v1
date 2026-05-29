#!/usr/bin/env Rscript
# render_supp_xenium_dotplot_cytoplasmic.R
#
# Supplemental Fig 3: Xenium cytoplasmic-counts dotplot at L1.5 resolution.
# Cross-segmentation companion to fig3e (Xenium nuclear, L1.5). Cell-ID lists
# are byte-identical between pooled_nuclear/ and pooled_cytoplasmic/ bundles,
# so labels are joined on cell_id via joint_l1p5.csv (platform == "xenium",
# is_artifact == FALSE). The matrix is the only thing that differs — same
# cells, whole-cell segmentation vs nuclear segmentation.
#
# Hard-fork of publication/figures/render/xenium/render_fig3e_xenium_dotplot_nuclear.R
# with the bundle filename swapped. At promotion the two variants will fold
# into one parameterized render_xenium_dotplot_l1p5.R with --segmentation.
#
# Drafting status: candidate for submission/figures/fig3/.drafting/panels/
#                  s_xenium_cytoplasmic_dotplot_markers_l1p5/v1/panel.pdf

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
parser$add_argument("--features-csv",   required = TRUE,
                    help = "fig3_dotplot_features.csv (feature, lbridge_family, ...)")
parser$add_argument("--xenium-bundle",  required = TRUE,
                    help = "Xenium cytoplasmic bundle dir (xenium_genes.tsv, xenium_cells.tsv, xenium_cytoplasmic_counts.mtx.gz)")
parser$add_argument("--joint-l1p5",     required = TRUE,
                    help = "joint_l1p5.csv (cell_id, l1p5_label, is_artifact)")
parser$add_argument("--out-dir",        required = TRUE)
parser$add_argument("--panel-id",       default = "s_xenium_cytoplasmic_dotplot_markers_l1p5")
args <- parser$parse_args()

PANEL_ID <- args$panel_id

# ---- Feature curation ------------------------------------------------------
features <- fread(args$features_csv)
message(sprintf("[%s] %d features × %d LBridge families",
                PANEL_ID, nrow(features),
                length(unique(features$lbridge_family))))

# ---- Xenium cytoplasmic bundle ---------------------------------------------
genes <- fread(file.path(args$xenium_bundle, "xenium_genes.tsv"), header = FALSE)$V1
cells <- fread(file.path(args$xenium_bundle, "xenium_cells.tsv"), header = FALSE)$V1
message(sprintf("[%s] Loading MTX (cytoplasmic): %d cells × %d genes …",
                PANEL_ID, length(cells), length(genes)))
mat <- t(as(readMM(file.path(args$xenium_bundle, "xenium_cytoplasmic_counts.mtx.gz")),
            "CsparseMatrix"))
dimnames(mat) <- list(genes, cells)

# Subset to features in panel
features_present <- as.character(features$feature[features$feature %in% genes])
features_missing <- setdiff(as.character(features$feature), genes)
if (length(features_missing) > 0) {
  message(sprintf("[%s] %d features not in 280-gene panel: %s",
                  PANEL_ID, length(features_missing),
                  paste(features_missing, collapse = ", ")))
}
features_kept <- features[feature %in% features_present]
features_kept[, feature := factor(feature, levels = features_kept$feature)]
features_kept[, lbridge_family := factor(lbridge_family,
                                         levels = unique(features_kept$lbridge_family))]
mat <- mat[as.character(features_kept$feature), , drop = FALSE]

# ---- L1.5 labels (xenium platform; cell-IDs match nuclear bundle) ----------
labs <- fread(args$joint_l1p5)[, .(cell_id, platform, l1p5_label, is_artifact)]
labs <- labs[platform == "xenium" & is_artifact == FALSE & cell_id %in% cells]
mat  <- mat[, labs$cell_id, drop = FALSE]
stopifnot(ncol(mat) == nrow(labs))
message(sprintf("[%s] %d cells × %d features after filtering (expect ~701K)",
                PANEL_ID, ncol(mat), nrow(mat)))

# ---- L1.5 column order (epi → str → imm) -----------------------------------
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

# ---- Aggregate per (feature × L1.5 label) ----------------------------------
groups <- split(seq_len(ncol(mat)), labs$l1p5_label)

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

# ---- Horizontal orientation for vertical-stacked composite ----------------
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

# Dims chosen for vertical stacking with FLEX-L1.5 + Xenium-nuclear-L1.5
# (composite ~16 x 18 in across all 3).
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
