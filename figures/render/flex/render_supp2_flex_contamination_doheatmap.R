#!/usr/bin/env Rscript
# render_supp2_flex_contamination_doheatmap.R — Fig 3 Supp 2, Panel 2b.
#
# Replicates the DoHeatmap panel from
# Spatial_HBCA/project/01_Preprocessing/scripts/step_15_milo_contamination/
# 4_decision_application.Rmd (lines 852-897), promoted to publication
# render with framework conformance.
#
# Inputs:
#   --seurat-obj     path to integrated_seurat.rds (step 09 or 10a)
#   --cell-comp      path to cell_compartments_step15.csv (cluster_annotation)
#   --out-dir        directory for PDF + companion CSV
#
# Outputs:
#   supp2_flex_contamination_doheatmap.pdf
#   supp2_flex_contamination_doheatmap_data.csv  (top markers per group)
#
# Method (verbatim from 4_decision_application.Rmd):
#   1. Load seurat_obj, attach cluster_annotation from cell_compartments
#   2. Idents <- cluster_annotation, drop NA cells
#   3. presto::wilcoxauc on data slot
#   4. Top 10 markers per group, filter auc > 0.6 & pct_in > 25
#   5. ceiling(70 / n_groups) per group, deduped
#   6. Downsample 1500 cells per group (max)
#   7. DoHeatmap, group.by = cluster_annotation, slot = "data"

suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(dplyr)
  library(presto)
  library(ggplot2)
})

parser <- ArgumentParser()
parser$add_argument("--seurat-obj", required = TRUE)
parser$add_argument("--cell-comp", required = TRUE)
parser$add_argument("--project-root", required = TRUE)
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--aesthetics-r",
                    default = NULL,
                    help = "load_aesthetics.R path (optional; falls back to default)")
parser$add_argument("--auc-threshold", type = "double", default = 0.6)
parser$add_argument("--pct-in-threshold", type = "double", default = 25)
parser$add_argument("--downsample-per-group", type = "integer", default = 1500)
args <- parser$parse_args()

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

# Aesthetics
if (!is.null(args$aesthetics_r) && file.exists(args$aesthetics_r)) {
  source(args$aesthetics_r)
}

cat(sprintf("[%s] Loading seurat_obj: %s\n", Sys.time(), args$seurat_obj))
seurat_obj <- readRDS(args$seurat_obj)
cat(sprintf("[%s] Seurat dims: %d x %d\n",
            Sys.time(), nrow(seurat_obj), ncol(seurat_obj)))

cat(sprintf("[%s] Loading cell_compartments: %s\n", Sys.time(), args$cell_comp))
cc <- read.csv(args$cell_comp, stringsAsFactors = FALSE)
cat(sprintf("[%s] cell_compartments: %d rows\n", Sys.time(), nrow(cc)))

# Attach cluster_annotation by cell_id
cell_ids <- colnames(seurat_obj)
match_ix <- match(cell_ids, cc$cell_id)
cluster_ann <- cc$cluster_annotation[match_ix]
seurat_obj$cluster_annotation <- cluster_ann

n_with_ann <- sum(!is.na(cluster_ann))
cat(sprintf("[%s] Cells with annotation: %d / %d\n",
            Sys.time(), n_with_ann, length(cell_ids)))

seurat_obj <- subset(seurat_obj,
                     cells = cell_ids[!is.na(cluster_ann)])
Idents(seurat_obj) <- "cluster_annotation"

cat(sprintf("[%s] Annotation counts:\n", Sys.time()))
print(table(seurat_obj$cluster_annotation))

# Marker computation (presto::wilcoxauc — 1-vs-rest per group)
cat(sprintf("[%s] Running presto::wilcoxauc\n", Sys.time()))
celltype_markers <- presto::wilcoxauc(seurat_obj, assay = "data")
cat(sprintf("[%s] Markers: %d rows\n", Sys.time(), nrow(celltype_markers)))

# Top markers per group, replicating 4_decision_application.Rmd
top_celltype_markers <- celltype_markers %>%
  filter(auc > args$auc_threshold, pct_in > args$pct_in_threshold) %>%
  group_by(group) %>%
  slice_max(order_by = auc, n = 15) %>%
  arrange(group, desc(auc))

# Markers for heatmap: ceiling(70 / n_groups) per group, deduped
n_groups <- length(unique(seurat_obj$cluster_annotation))
per_group <- ceiling(70 / n_groups)

markers_for_heatmap <- top_celltype_markers %>%
  group_by(group) %>%
  slice_head(n = per_group) %>%
  pull(feature) %>%
  unique()

cat(sprintf("[%s] Plotting %d markers across %d groups (top %d per group)\n",
            Sys.time(), length(markers_for_heatmap), n_groups, per_group))

# Downsample for the heatmap
cells_downsample_prop <- subset(seurat_obj,
                                downsample = args$downsample_per_group) |>
  colnames()

cat(sprintf("[%s] Downsampled to %d cells\n",
            Sys.time(), length(cells_downsample_prop)))

# Order annotations: contamination first, then by compartment
annotation_order <- c(
  "Keratinocyte", "Eccrine_Gland", "Melanocyte",
  "Basal", "Luminal_HR", "Luminal_Secretory",
  "Contractile_Myoepithelial", "PS_Epithelial",
  "Fibroblast", "Pericyte", "VSMC", "Adipocyte",
  "Vascular_Endothelium", "Lymphatic_Endothelium", "Schwann_Cell",
  "T_Cell", "B_Cell", "Plasma_cell", "Myeloid", "Mast",
  "Neutrophil", "pDC", "Langerhans_Cell", "RBC"
)
present_order <- annotation_order[annotation_order %in% unique(seurat_obj$cluster_annotation)]
seurat_obj$cluster_annotation <- factor(seurat_obj$cluster_annotation,
                                        levels = present_order)
Idents(seurat_obj) <- "cluster_annotation"

# Render
p <- DoHeatmap(seurat_obj[, cells_downsample_prop],
               features = markers_for_heatmap,
               group.by = "cluster_annotation",
               slot = "data",
               raster = FALSE,
               size = 2.0) +
  theme(axis.text.y = element_text(size = 4),
        plot.title = element_blank())

out_pdf <- file.path(args$out_dir, "supp2_flex_contamination_doheatmap.pdf")
ggsave(out_pdf, plot = p, width = 9.0, height = 8.5, units = "in",
       device = "pdf", dpi = 600)
cat(sprintf("[%s] Wrote %s\n", Sys.time(), out_pdf))

# Companion CSV: top markers per group
out_csv <- file.path(args$out_dir, "supp2_flex_contamination_doheatmap_data.csv")
write.csv(top_celltype_markers %>% group_by(group) %>% slice_head(n = per_group),
          out_csv, row.names = FALSE)
cat(sprintf("[%s] Wrote %s\n", Sys.time(), out_csv))

cat(sprintf("[%s] === Complete ===\n", Sys.time()))
