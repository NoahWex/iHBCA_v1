#!/usr/bin/env Rscript
# render_fig3d_flex_dotplot_l2s.R
#
# Fig 3d: FLEX-cohort canonical-marker dotplot. Rows = LBridge-bracketed
# canonical features. Columns = FLEX L2S labels (Track C annotation_v2s).
# Companion to fig3e (Xenium nuclear, L1.5 columns) for cross-modal marker
# validation. Xenium can't resolve L2S without the FLEX transcriptome — that's
# why FLEX gets the higher-resolution column tier.
#
# L2S labels not in the LBridge vocab need explicit handling:
#   - ARTIFACT_* → dropped (treated as is_artifact)
#   - Adipocyte → LBridge family "Adipocyte" (FLEX/Xenium-only family)
#   - LHS-apocrine → LBridge family "LHS" (apocrine sub-state of LHS)

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
parser$add_argument("--features-csv",       required = TRUE,
                    help = "fig3_dotplot_features.csv")
parser$add_argument("--flex-bundle-root",   required = TRUE,
                    help = "FLEX integration_intermediate root containing per-compartment scvi_n100/")
parser$add_argument("--joint-l1p5",         required = TRUE,
                    help = "joint_l1p5.csv (cell_id, is_artifact)")
parser$add_argument("--l2s-labels",         required = TRUE,
                    help = "flex_l2s_labels.csv (cell_id, l2s_label)")
parser$add_argument("--lbridge-vocab",      required = TRUE,
                    help = "lbridge_vocabulary.csv")
parser$add_argument("--out-dir",            required = TRUE)
parser$add_argument("--panel-id",           default = "fig3d_flex_dotplot_markers_l2s")
args <- parser$parse_args()

PANEL_ID <- args$panel_id

# ---- Feature curation ------------------------------------------------------
features <- fread(args$features_csv)
message(sprintf("[%s] %d features × %d LBridge families",
                PANEL_ID, nrow(features),
                length(unique(features$lbridge_family))))

# ---- L2S → LBridge family map ----------------------------------------------
# Three-source precedence: features CSV > manual_map > lbridge_voc.
#
# Features CSV is authoritative for any L2S label that appears as a row block
# (e.g., BMYO-myo, BMYO-basal, Fibro-SFRP4, Adipo-* — these are treated as
# their own LBridge family for row-block resolution, even though lbridge_voc
# may collapse them to a parent family like BMYO or Adipocyte).
#
# manual_map covers L2S labels that should route to a known family but don't
# have a feature row (e.g., LHS-apocrine → LHS).
#
# lbridge_voc is the fallback for everything else.
features_map <- unique(fread(args$features_csv)[, .(l2s_label = l2s_label,
                                                     lbridge_family,
                                                     lbridge_compartment)])
manual_map <- data.table(
  l2s_label = c("LHS-apocrine"),
  lbridge_family = c("LHS"),
  lbridge_compartment = c("Epithelial")
)
lbridge_voc <- fread(args$lbridge_vocab)
voc_map <- unique(lbridge_voc[!is.na(ihbca_l20_label) & ihbca_l20_label != "NA",
                              .(l2s_label = ihbca_l20_label,
                                lbridge_family,
                                lbridge_compartment)])
l2s_to_lb <- rbind(features_map, manual_map, voc_map)
l2s_to_lb <- l2s_to_lb[, .SD[1], by = l2s_label]
message(sprintf("[%s] L2S→LBridge map: %d L2S labels mapped",
                PANEL_ID, nrow(l2s_to_lb)))

# ---- Load FLEX per-compartment bundles + concatenate ----------------------
load_compartment <- function(comp) {
  bundle <- file.path(args$flex_bundle_root, comp, "scvi_n100")
  cells  <- fread(file.path(bundle, "cells.tsv"), header = FALSE)$V1
  genes  <- fread(file.path(bundle, "genes.tsv"), header = FALSE)$V1
  message(sprintf("[%s]   %s: %d cells × %d genes",
                  PANEL_ID, comp, length(cells), length(genes)))
  mat <- t(as(readMM(file.path(bundle, "counts.mtx.gz")), "CsparseMatrix"))
  dimnames(mat) <- list(genes, cells)
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
# Same gene can map to multiple L2S labels in the curated feature list (e.g.,
# RAMP3 marks both Vas-vein and Vas-capillary). One dotplot row per gene —
# keep the first occurrence, which determines the row block and y-axis order.
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

# ---- L2S labels via flex_l2s_labels.csv -----------------------------------
l2s   <- fread(args$l2s_labels)[, .(cell_id, l2s_label)]
labs  <- fread(args$joint_l1p5)[, .(cell_id, is_artifact)]
labs  <- merge(labs, l2s, by = "cell_id", all.x = FALSE, all.y = FALSE)
labs  <- merge(labs, l2s_to_lb, by = "l2s_label", all.x = TRUE)
labs  <- labs[cell_id %in% colnames(mat) &
              is_artifact == FALSE &
              !startsWith(l2s_label, "ARTIFACT_") &
              !is.na(lbridge_family)]
mat   <- mat[, labs$cell_id, drop = FALSE]
stopifnot(ncol(mat) == nrow(labs))
message(sprintf("[%s] %d cells × %d features after L2S filter (dropped artifacts + unmapped)",
                PANEL_ID, ncol(mat), nrow(mat)))

# ---- L2S column order: by LBridge family (matching feature row order),
#      within family by cell count (descending). LBridge families with no
#      curated features (e.g., Xenium-only) are dropped.
lbridge_order <- levels(features_kept$lbridge_family)
labs[, lbridge_family := factor(lbridge_family, levels = lbridge_order)]
labs_counts <- labs[, .N, by = .(lbridge_family, l2s_label)]
labs_counts <- labs_counts[order(lbridge_family, -N)]
l2s_order   <- as.character(labs_counts$l2s_label)

labs <- labs[lbridge_family %in% lbridge_order]
labs[, l2s_label := factor(l2s_label, levels = l2s_order)]
mat  <- mat[, labs$cell_id, drop = FALSE]
message(sprintf("[%s] %d L2S labels across %d LBridge families (column dimension)",
                PANEL_ID, length(l2s_order), length(unique(labs$lbridge_family))))

# ---- Aggregate per (feature × L2S label) ----------------------------------
groups <- split(seq_len(ncol(mat)), labs$l2s_label)
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

# ---- LBridge row-split factor + bracket -----------------------------------
row_split <- features_kept$lbridge_family
names(row_split) <- as.character(features_kept$feature)

lbridge_anno <- HeatmapAnnotation(
  lbridge = anno_block(
    gp     = gpar(fill = NA, col = NA),
    labels = NULL
  ),
  show_annotation_name = FALSE,
  show_legend          = FALSE,
  width                = unit(2.5, "cm"),
  which                = "row"
)

# ---- Color + cell function -------------------------------------------------
pal_z <- colorRamp2(c(-cap, 0, cap), c("#3B4992", "#FFFFFF", "#EE0000"))

cell_fun <- function(j, i, x, y, width, height, fill) {
  pct <- mat_pct[i, j]
  zv  <- mat_z[i, j]
  if (is.finite(pct) && pct > 0.02) {
    # Dot radius capped to row-height envelope. With L2S subtype expansion
    # row count climbed to ~180; row height was compressed to 0.20cm (2.0mm),
    # so max radius = 0.95mm keeps dots inside the row.
    r_mm <- 0.95 * sqrt(pct)
    grid.circle(
      x  = x, y = y,
      r  = unit(r_mm, "mm"),
      gp = gpar(fill = pal_z(zv), col = NA)
    )
  }
}

ht <- Heatmap(
  mat_z,
  name             = "z(log1p)",
  col              = pal_z,
  rect_gp          = gpar(type = "none"),
  cell_fun         = cell_fun,
  cluster_rows     = FALSE, cluster_columns = FALSE,
  show_row_names   = TRUE,  show_column_names = TRUE,
  row_names_gp     = gpar(fontsize = 5),
  column_names_gp  = gpar(fontsize = 5),
  column_names_rot = 60,
  row_split        = row_split,
  row_gap          = unit(1.0, "mm"),
  row_title        = NULL,
  left_annotation  = lbridge_anno,
  width            = unit(ncol(mat_z) * 0.45, "cm"),
  height           = unit(nrow(mat_z) * 0.20, "cm"),
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
        x = 0.85, y = 0.5,
        rot = 0, just = c("right", "center"),
        gp = gpar(fontsize = 5, fontface = "bold")
      )
    })
  }
}

# ---- Save ------------------------------------------------------------------
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(args$out_dir, paste0(PANEL_ID, ".pdf"))
png_path <- file.path(args$out_dir, paste0(PANEL_ID, ".png"))

pdf(pdf_path, width = 13, height = 16)
draw(ht)
decorate_lbridge()
dev.off()

png(png_path, width = 13, height = 16, units = "in", res = 600)
draw(ht)
decorate_lbridge()
dev.off()

message(sprintf("[%s] Saved PDF=%s, PNG=%s (13x16 in)",
                PANEL_ID, pdf_path, png_path))
