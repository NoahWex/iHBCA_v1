# =============================================================================
# render_fig1c_marker_heatmap.R — all-compartment canonical marker heatmap
# =============================================================================
# Reads:
#   - mean_expression_per_label.csv  (long: compartment, label, gene, mean_log1p)
#   - canonical_markers.csv          (long: compartment, label, gene, marker_rank)
# Produces:
#   - fig1_ihbca_marker_heatmap_allcomp.pdf + .png
#
# Design:
#   - Rows: 43 non-artifact (compartment, label) pairs. Order = compartment
#     (epi -> str -> imm) > LBridge family > label. Row labels show only the
#     L2 label (compartment dropped — column annotation conveys block).
#   - Columns: ~246 unique canonical genes. Each gene assigned to the
#     (compartment, label) where its mean_log1p is HIGHEST. Columns ordered
#     by (assigned label's compartment, LBridge family, source label, rank).
#     This produces the cleanest block-diagonal possible without ambiguity.
#   - Column annotation: LBridge family bars — each family draws a single
#     horizontal segment above the heatmap with the family name. Replaces
#     the prior "source" annotation.
#   - No row "compartment" annotation strip; row labels speak for themselves
#     and column LBridge bars provide grouping.
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(yaml)
})

parser <- ArgumentParser()
parser$add_argument("--mean-expression", required = TRUE)
parser$add_argument("--canonical-markers", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--aesthetics-r", required = TRUE)
parser$add_argument("--project-root", required = TRUE)
args <- parser$parse_args()

Sys.setenv(PROJECT_ROOT = args$project_root)
source(args$aesthetics_r)

# -----------------------------------------------------------------------------
# LBridge family mapping (per ~/.claude/plans/...wombat.md, 30-family draft)
# -----------------------------------------------------------------------------
# Used as column annotation bars + row ordering. iHBCA L2 -> LBridge family.
LBRIDGE <- list(
  # Epithelial (6)
  "BMYO-myo"        = "BMYO",
  "BMYO-basal"      = "BMYO",
  "LASP-major"      = "LASP",
  "LASP-KIT"        = "LASP",
  "LASP-HLA"        = "LASP",
  "LASP-basal"      = "LASP-basal",
  "LHS-major"       = "LHS",
  "Lactocyte-LC1"   = "Lactocyte",
  "Lactocyte-LC2"   = "Lactocyte",
  # Stromal (5; Adipocyte not in iHBCA)
  "Fibro-major"     = "Fibroblast",
  "Fibro-prematrix" = "Fibroblast",
  "Fibro-IGF1"      = "Fibroblast",
  "Fibro-SFRP4"     = "Fibro-SFRP4",
  "Vas-capillary"   = "Endothelial",
  "Vas-vein"        = "Endothelial",
  "Vas-arterial"    = "Endothelial",
  "Lym-major"       = "Lymphatic_Endothelial",
  "Pericyte"        = "Perivascular",
  "Pericyte_active" = "Perivascular",
  "VSMC"            = "Perivascular",
  "VSMC_stress"     = "Perivascular",
  # Immune (15)
  "Macro_C1Q"       = "Macrophage",
  "Macro_FOLR2"     = "Macrophage",
  "Macro_LAM"       = "Macrophage",
  "Macro_inflam"    = "Macrophage",
  "Mono_NC"         = "Monocyte",
  "cDC1"            = "cDC1",
  "cDC2"            = "cDC2",
  "pDC"             = "pDC",
  "Mast"            = "Mast",
  "Neutrophil"      = "Neutrophil",
  "B_cell"          = "B",
  "Plasma"          = "Plasma",
  "Treg"            = "Treg",
  "CD4_Th_like"     = "CD4 T",
  "Th17"            = "CD4 T",
  "CD8_Trm"         = "CD8 T",
  "CD8_Resting"     = "CD8 T",
  "CD8_Tem"         = "CD8 T",
  "NK"              = "NK",
  "NK_ILC_prolif"   = "NK",
  "IFNg_T"          = "CD8 T",
  "ZNF683_T"        = "CD8 T"
)

# Display order across compartments (epi -> str -> imm)
LBRIDGE_ORDER <- c(
  # Epi
  "BMYO", "LASP", "LASP-basal", "LHS", "Lactocyte",
  # Str
  "Fibroblast", "Fibro-SFRP4", "Endothelial", "Lymphatic_Endothelial", "Perivascular",
  # Imm  (T_ambiguous folded into CD8 T)
  "Macrophage", "Monocyte", "cDC1", "cDC2", "pDC", "Mast", "Neutrophil",
  "B", "Plasma", "Treg", "CD4 T", "CD8 T", "NK"
)
COMPARTMENT_ORDER <- c("epi", "str", "imm")

# -----------------------------------------------------------------------------
# Load + filter
# -----------------------------------------------------------------------------

cat("[fig1c] loading mean expression:", args$mean_expression, "\n")
me <- fread(args$mean_expression)
cat("[fig1c]   rows:", nrow(me), "\n")

cat("[fig1c] loading canonical markers:", args$canonical_markers, "\n")
mk <- fread(args$canonical_markers)
mk <- mk[is_artifact == FALSE]   # drop doublets
cat("[fig1c]   non-artifact marker rows:", nrow(mk), "\n")
cat("[fig1c]   unique non-artifact labels:", uniqueN(mk[, .(compartment, label)]), "\n")

me <- merge(me, unique(mk[, .(compartment, label)]),
            by = c("compartment", "label"))

# Attach LBridge to labels (rows + gene-source)
me[, lbridge := unlist(LBRIDGE[label])]
mk[, lbridge := unlist(LBRIDGE[label])]
unmapped <- me[is.na(lbridge), unique(label)]
if (length(unmapped) > 0) {
  cat("[fig1c]   WARNING — labels missing from LBRIDGE:", paste(unmapped, collapse = ", "), "\n")
}

# -----------------------------------------------------------------------------
# Row order: compartment > LBridge family > label
# -----------------------------------------------------------------------------
me[, comp_order := factor(compartment, levels = COMPARTMENT_ORDER)]
me[, lbridge_order := factor(lbridge, levels = LBRIDGE_ORDER)]
row_keys <- unique(me[, .(comp_order, compartment, lbridge_order, lbridge, label)])
setorder(row_keys, comp_order, lbridge_order, label)
row_keys[, row_id := paste(compartment, label, sep = "/")]
cat("[fig1c]   rows:", nrow(row_keys), "\n")

# -----------------------------------------------------------------------------
# Column order: assign each gene to its CANONICAL label (where it has the
# lowest marker_rank — i.e., where the author placed it earliest in their
# canonical list). Honors author intent over signal-strength heuristic.
# Tie-break (gene canonical at the same rank in multiple labels): prefer the
# label with higher mean_log1p for the gene.
# -----------------------------------------------------------------------------
mk_with_signal <- merge(mk, me[, .(compartment, label, gene, mean_log1p)],
                        by = c("compartment", "label", "gene"), all.x = TRUE)
mk_with_signal[, mean_log1p := nafill(mean_log1p, fill = 0)]
setorder(mk_with_signal, gene, marker_rank, -mean_log1p)
gene_to_assigned <- mk_with_signal[, .SD[1L], by = gene][
  , .(gene, assigned_compartment = compartment, assigned_label = label,
      assigned_lbridge = lbridge, rank = marker_rank)
]
gene_to_assigned[, comp_order := factor(assigned_compartment, levels = COMPARTMENT_ORDER)]
gene_to_assigned[, lbridge_order := factor(assigned_lbridge, levels = LBRIDGE_ORDER)]
setorder(gene_to_assigned, comp_order, lbridge_order, assigned_label, rank)
gene_order <- gene_to_assigned$gene
cat("[fig1c]   genes in display order:", length(gene_order), "\n")

# -----------------------------------------------------------------------------
# Wide matrix
# -----------------------------------------------------------------------------
me[, row_id := paste(compartment, label, sep = "/")]
mat <- dcast(me, row_id ~ gene, value.var = "mean_log1p", fill = 0)
mat_rownames <- mat$row_id
mat <- as.matrix(mat[, -1])
rownames(mat) <- mat_rownames

mat <- mat[row_keys$row_id, intersect(gene_order, colnames(mat))]
cat("[fig1c]   matrix shape:", paste(dim(mat), collapse = " x "), "\n")

# Display rownames: just the L2 label (drop compartment prefix)
display_rownames <- row_keys$label
rownames(mat) <- display_rownames

# Column display: gene name only (no annotation suffix)
# colnames already gene symbols.

# Z-score per column, cap +/-3
mat_z <- scale(mat, center = TRUE, scale = TRUE)
mat_z[is.na(mat_z)] <- 0
cap <- 3
mat_z[mat_z > cap] <- cap
mat_z[mat_z < -cap] <- -cap

# -----------------------------------------------------------------------------
# Render with pheatmap. LBridge bars drawn manually as a top decoration grob
# rather than via annotation_col (which only supports per-cell color squares).
# Strategy: render heatmap with no col annotation, then on the saved PDF/PNG
# we'd add manual segments — but pheatmap doesn't expose that easily.
#
# Solution: emit a per-LBridge-family segment table alongside the figure
# (start_col, end_col, family_name) so the user (or a follow-up render) can
# overlay them in Illustrator OR re-render with ComplexHeatmap if we want
# baked-in bars. For now, ship the heatmap clean with column-name labels;
# emit segments for downstream composition.
# -----------------------------------------------------------------------------

# Compute LBridge family segments: start/end column indices in the displayed matrix
gene_order_in_matrix <- colnames(mat_z)
lbridge_per_col <- gene_to_assigned[match(gene_order_in_matrix, gene), assigned_lbridge]
seg_runs <- rle(lbridge_per_col)
seg_starts <- cumsum(c(1, head(seg_runs$lengths, -1)))
seg_ends   <- cumsum(seg_runs$lengths)
seg_table <- data.table(
  lbridge_family = seg_runs$values,
  col_start = seg_starts,
  col_end   = seg_ends,
  n_genes   = seg_runs$lengths
)

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
seg_path <- file.path(args$out_dir, "fig1c_lbridge_segments.csv")
fwrite(seg_table, seg_path)
cat("[fig1c]   wrote LBridge column segments:", seg_path, "\n")

pal_diverging <- colorRamp2(c(-cap, 0, cap), c("#3B4992", "#FFFFFF", "#EE0000"))
pdf_path <- file.path(args$out_dir, "fig1_ihbca_marker_heatmap_allcomp.pdf")
png_path <- file.path(args$out_dir, "fig1_ihbca_marker_heatmap_allcomp.png")

# ComplexHeatmap top annotation: per-LBridge horizontal LINE SEGMENT (no fill)
# spanning the family's column range, with the LBridge name rotated 45 degrees
# above the segment.
seg_centers <- (seg_starts + seg_ends) / 2

lbridge_anno <- HeatmapAnnotation(
  lbridge = anno_block(
    gp = gpar(fill = NA, col = NA),               # invisible rect (we draw a line below)
    labels = NULL,
    labels_gp = gpar(fontsize = 5, fontface = "bold")
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE,
  height = unit(1.6, "cm"),
  which = "column"
)

column_split_factor <- factor(lbridge_per_col, levels = unique(lbridge_per_col))

ht <- Heatmap(
  mat_z,
  name = "z(log1p)",
  col = pal_diverging,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 5),
  column_names_gp = gpar(fontsize = 4),
  column_names_rot = 90,
  top_annotation = lbridge_anno,
  rect_gp = gpar(col = NA),
  width = unit(ncol(mat_z) * 0.16, "cm"),
  height = unit(nrow(mat_z) * 0.24, "cm"),
  column_split = column_split_factor,
  column_gap = unit(1.5, "mm"),
  column_title = NULL,
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 5),
    labels_gp = gpar(fontsize = 4),
    legend_height = unit(2, "cm")
  )
)

decorate_lbridge <- function() {
  # Per-slice decoration: each LBridge column-block gets a thin horizontal
  # line at the bottom + 45-deg rotated label above. Column-split gaps
  # already break continuity between families.
  fams <- levels(column_split_factor)
  for (i in seq_along(fams)) {
    decorate_annotation("lbridge", slice = i, {
      grid.lines(
        x = c(0.05, 0.95), y = c(0.05, 0.05),
        gp = gpar(col = "black", lwd = 1.0)
      )
      grid.text(
        fams[i],
        x = 0.5, y = 0.18,
        rot = 45, just = c("left", "bottom"),
        gp = gpar(fontsize = 5, fontface = "bold")
      )
    })
  }
}

pdf(pdf_path, width = 20, height = 9)
draw(ht)
decorate_lbridge()
dev.off()
cat("[fig1c] wrote PDF:", pdf_path, "\n")

png(png_path, width = 20, height = 9, units = "in", res = 600)
draw(ht)
decorate_lbridge()
dev.off()
cat("[fig1c] wrote PNG:", png_path, "\n")
cat("[fig1c] segments table:", seg_path, "\n")
cat("[fig1c] done\n")
