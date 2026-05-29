#!/usr/bin/env Rscript
# 06_volcano_fibros_caf_labels.R
# Fibroblast volcano plots labeling a curated CAF marker panel only
# (iCAF + myCAF + matrix-CAF + pan-CAF + related). All 4 fibros under A1_bio.
#
# Output: outputs/plots/fibro_volcanos_caf_labeled/

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

p <- ArgumentParser()
p$add_argument("--project-root", required = TRUE)
p$add_argument("--fdr", type = "double", default = 0.05)
p$add_argument("--lfc", type = "double", default = 0.5)
args <- p$parse_args()

project_root <- args$project_root
out_dir <- file.path(project_root, "outputs", "plots", "fibro_volcanos_caf_labeled")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

FORMULA <- "A2_sampletype"
COHORT  <- "A_full"

# Pretty-display L2 name (drop "str__" prefix, swap underscore for hyphen)
pretty_l2 <- function(l2) {
  s <- sub("^[a-z]+__", "", l2)
  gsub("_", "-", s)
}

# ---- curated CAF marker panel ------------------------------------------
# iCAF (inflammatory): cytokines + chemokines + mediators + surface
ICAF <- c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","CXCL12","CXCL14",
          "CCL2","CCL7","LIF","PTGS2","C3",
          "HAS1","HAS2","MMP1","MMP3","TIMP1","SERPINE1",
          "PDPN","DPT","DPP4","ICAM1")

# myCAF (myofibroblast / contractile)
MYCAF <- c("ACTA2","TAGLN","MYH11","MCAM","RGS5","CNN1")

# Matrix CAF / desmoplastic
MATCAF <- c("POSTN","COMP","COL1A1","COL3A1","COL10A1","COL11A1","COL12A1",
            "FN1","ASPN","VCAN","LRRC15")

# Pan-CAF / activation
PANCAF <- c("FAP","THY1","S100A4","PDGFRA","PDGFRB")

# Related (matrix remodeling, signaling)
OTHER <- c("LOX","LOXL2","MMP2","MMP11","CCN2","CTGF",
            "SFRP2","SFRP4","GREM1","INHBA")

marker_class <- c(
  setNames(rep("iCAF",       length(ICAF)),   ICAF),
  setNames(rep("myCAF",      length(MYCAF)),  MYCAF),
  setNames(rep("matrix",     length(MATCAF)), MATCAF),
  setNames(rep("pan-CAF",    length(PANCAF)), PANCAF),
  setNames(rep("other",      length(OTHER)),  OTHER)
)
all_markers <- names(marker_class)

CLASS_COLORS <- c(
  "iCAF"    = "#E69F00",   # orange — inflammatory
  "myCAF"   = "#0072B2",   # blue — myofibroblast
  "matrix"  = "#009E73",   # green — matrix
  "pan-CAF" = "#CC79A7",   # pink — pan-activation
  "other"   = "#56B4E9"    # light blue — related
)

PALETTE_BG <- c("BR1 up"   = "#FBE4C2",   # light orange
                "BR1 down" = "#C6DDF0",   # light blue
                "ns"       = "#E5E5E5")   # light grey

safe_name <- function(s) gsub("[^A-Za-z0-9_]", "_", s)

panel_theme <- function() {
  theme_classic(base_size = 8) +
    theme(
      panel.grid       = element_blank(),
      axis.line        = element_line(linewidth = 0.3, color = "black"),
      axis.ticks       = element_line(linewidth = 0.25, color = "black"),
      axis.text        = element_text(size = 7, color = "black"),
      axis.title       = element_text(size = 7.5, color = "black"),
      plot.margin      = margin(4, 6, 4, 6),
      legend.position  = "none"
    )
}

build_volcano <- function(L2) {
  csv_path <- file.path(project_root, "outputs", "de_results",
                        COHORT, FORMULA, paste0(safe_name(L2), ".csv"))
  if (!file.exists(csv_path)) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0, y = 0, label = paste0(L2, ": missing")))
  }
  hdr <- readr::read_lines(csv_path, n_max = 1)
  if (grepl("^skipped", hdr)) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0, y = 0, label = paste0(L2, ": gated")))
  }

  de <- readr::read_csv(csv_path, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    dplyr::mutate(
      minus_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
      sig = dplyr::case_when(
        adj.P.Val < args$fdr & logFC >  args$lfc ~ "BR1 up",
        adj.P.Val < args$fdr & logFC < -args$lfc ~ "BR1 down",
        TRUE                                      ~ "ns"
      )
    )

  # Background points (de-emphasized — light palette)
  de_bg <- de %>% dplyr::filter(!(symbol %in% all_markers))
  # Marker points (highlighted by class)
  de_mk <- de %>% dplyr::filter(symbol %in% all_markers) %>%
    dplyr::mutate(marker_class = marker_class[symbol])

  ggplot() +
    geom_point(data = de_bg,
               aes(x = logFC, y = minus_log10_fdr, color = sig),
               size = 0.35, alpha = 0.7, shape = 16) +
    scale_color_manual(values = PALETTE_BG) +
    ggnewscale::new_scale_color() +
    geom_point(data = de_mk,
               aes(x = logFC, y = minus_log10_fdr, fill = marker_class),
               size = 1.8, alpha = 1, shape = 21, color = "black", stroke = 0.25) +
    scale_fill_manual(values = CLASS_COLORS) +
    geom_hline(yintercept = -log10(args$fdr), linetype = "dashed",
               color = "grey50", linewidth = 0.25) +
    geom_vline(xintercept = c(-args$lfc, args$lfc), linetype = "dashed",
               color = "grey50", linewidth = 0.25) +
    ggrepel::geom_text_repel(data = de_mk %>% dplyr::filter(sig != "ns"),
                              aes(x = logFC, y = minus_log10_fdr, label = symbol),
                              size = 2.2, color = "black",
                              max.overlaps = Inf,
                              segment.size = 0.15,
                              segment.color = "grey40",
                              min.segment.length = 0,
                              box.padding = 0.3,
                              point.padding = 0.2,
                              force = 2,
                              seed = 42) +
    labs(title = pretty_l2(L2),
         x = expression(log[2]~FC~(BR1 / AR)),
         y = expression(-log[10]~FDR)) +
    panel_theme() +
    theme(plot.title = element_text(size = 9, face = "bold", hjust = 0,
                                     margin = margin(b = 2)))
}

# ---- legend (built once, attached to the grid) -------------------------
legend_plot <- function() {
  # Dummy data to render legend
  legend_df <- data.frame(
    x = seq_along(CLASS_COLORS),
    y = 0,
    marker_class = factor(names(CLASS_COLORS), levels = names(CLASS_COLORS))
  )
  ggplot(legend_df, aes(x = x, y = y, fill = marker_class)) +
    geom_point(size = 3, shape = 21, color = "black", stroke = 0.3) +
    scale_fill_manual(values = CLASS_COLORS, name = NULL) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7),
          legend.key.size = grid::unit(0.4, "cm"),
          legend.spacing.x = grid::unit(0.3, "cm"))
}

# ---- run ---------------------------------------------------------------
fibros <- c("str__Fibro_major", "str__Fibro_IGF1",
            "str__Fibro_SFRP4", "str__Fibro_prematrix")

# Try loading ggnewscale; install or fail gracefully
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  stop("ggnewscale not installed — please install in container R library")
}

plots <- list()
for (L2 in fibros) {
  cat(sprintf("%s\n", L2))
  p <- build_volcano(L2)
  plots[[L2]] <- p
  ggsave(file.path(out_dir, paste0(safe_name(L2), "_caf.pdf")),
         p, width = 100, height = 100, units = "mm", dpi = 600)
}

# Combined 2x2 + legend strip
combined <- patchwork::wrap_plots(plots, ncol = 2)
final <- combined / legend_plot() + patchwork::plot_layout(heights = c(1, 0.04))
ggsave(file.path(out_dir, "fibros_caf_grid.pdf"),
       final, width = 200, height = 210, units = "mm", dpi = 600,
       limitsize = FALSE)

cat("\nDone.\n")
