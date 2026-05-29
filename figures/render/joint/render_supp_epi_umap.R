#!/usr/bin/env Rscript
# UMAP scatter colored by per-cell DA logFC for the C5 epi_content contrast.
#
# Question per panel: where in the cell-space (FLEX scvi_n100 or joint Concord)
# does the DA signal land?
#
# Plot structure:
#   X/Y:    UMAP1 / UMAP2 from the canonical embedding for this substrate
#   Color:  per-cell mean_logfc on sig cells (status sig+/sig-) via RdBu diverging;
#           non-sig cells in grey background
#   Order:  ns drawn first (background), sig drawn on top
#
# Single-question discipline:
#   - no panel title; no methods text
#   - no facet by patient or platform (those are separate panels if needed)
#   - rasterize the scatter at 600 DPI within a vector PDF wrapper (>100K points)

suppressPackageStartupMessages({
  library(argparse)
  library(data.table)
  library(ggplot2)
  library(ggrastr)
})

# load_aesthetics() uses %||% internally; define here so it's in scope
`%||%` <- function(a, b) if (is.null(a)) b else a

parser <- ArgumentParser()
parser$add_argument("--project-root",   required = TRUE)
parser$add_argument("--per-cell-csv",   required = TRUE,
                    help = "per_cell_logfc.csv from make_per_cell_da_{flex,joint}.{R,py}")
parser$add_argument("--umap-csv",       required = TRUE,
                    help = "UMAP coords CSV with columns cell_id, UMAP_1, UMAP_2 (or umap_1/umap_2)")
parser$add_argument("--out-pdf",        required = TRUE)
parser$add_argument("--width-in",       type = "double", default = 4.5)
parser$add_argument("--height-in",      type = "double", default = 4.5)
parser$add_argument("--point-size",     type = "double", default = 0.05)
parser$add_argument("--lfc-cap-quantile", type = "double", default = 0.98,
                    help = "Symmetric color limit from this quantile of |sig logFC|")
args <- parser$parse_args()

# Aesthetics framework
root <- normalizePath(args$project_root)
source(file.path(root, "publication", "config", "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(root, "publication", "config"))
div <- aes_cfg$scales$expression_diverging

# Load substrate
pc   <- fread(args$per_cell_csv)
umap <- fread(args$umap_csv)
# normalize UMAP column names
nm <- tolower(colnames(umap))
xcol <- colnames(umap)[grep("^umap_?1$", nm)[1]]
ycol <- colnames(umap)[grep("^umap_?2$", nm)[1]]
stopifnot("UMAP CSV missing umap_1/umap_2 columns" = !is.na(xcol) && !is.na(ycol))
setnames(umap, c(xcol, ycol), c("x", "y"))

# Find a cell_id col
cid_col <- if ("cell_id" %in% colnames(umap)) "cell_id" else colnames(umap)[1]
setnames(umap, cid_col, "cell_id")

d <- merge(umap[, .(cell_id, x, y)], pc, by = "cell_id", all.x = FALSE)
cat("[umap-da] cells matched: ", nrow(d), " (umap=", nrow(umap),
    "  per-cell=", nrow(pc), ")\n", sep = "")

# Symmetric color cap on sig cells
sig_lfc <- d$mean_logfc[d$status %in% c("sig+", "sig-")]
lfc_cap <- if (length(sig_lfc) > 0) {
  quantile(abs(sig_lfc), args$lfc_cap_quantile, na.rm = TRUE)
} else {
  0.5
}
if (!is.finite(lfc_cap) || lfc_cap == 0) lfc_cap <- 0.5

# Layered: ns underneath, sig on top
d_ns  <- d[status == "ns"]
d_sig <- d[status != "ns"]
cat("  ns=", nrow(d_ns), " sig=", nrow(d_sig),
    " lfc_cap=", round(lfc_cap, 3), "\n", sep = "")

p <- ggplot() +
  rasterise(geom_point(data = d_ns, aes(x = x, y = y),
                       color = "grey85", size = args$point_size, stroke = 0),
            dpi = 600) +
  rasterise(geom_point(data = d_sig, aes(x = x, y = y, color = mean_logfc),
                       size = args$point_size * 1.4, stroke = 0),
            dpi = 600) +
  scale_color_gradient2(low = div$low, mid = div$mid, high = div$high,
                        midpoint = 0, limits = c(-lfc_cap, lfc_cap),
                        oob = scales::squish,
                        guide = guide_colorbar(title = "logFC",
                                               barwidth = 0.4, barheight = 4)) +
  coord_fixed() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  get_theme(aes_cfg) +
  theme(panel.grid       = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.title       = element_blank(),
        plot.subtitle    = element_blank(),
        legend.position  = "right",
        axis.text        = element_blank(),
        axis.ticks       = element_blank())

dir.create(dirname(args$out_pdf), recursive = TRUE, showWarnings = FALSE)
ggsave(args$out_pdf, plot = p, width = args$width_in, height = args$height_in,
       units = "in", device = cairo_pdf)
cat("[umap-da] wrote ", args$out_pdf, "\n", sep = "")
