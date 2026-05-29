#!/usr/bin/env Rscript
# Fig 3 / supp — Xenium DA beeswarm grid WITH per-contrast position diagrams.
#
# Composes a 2-row figure:
#   Row 1: per-contrast position diagrams (PNG insets, one per column)
#   Row 2: beeswarm grid (motif M0..M8 OR L1.5 cell types on Y, faceted across
#          contrasts on X)
# Same publication_v2 ggbeeswarm aesthetic. Reads same contrasts CSV as
# render_xenium_motif_beeswarm_grid.R, with an added `diagram_path` column.
#
# --level motif:  group by motif via per_nhood_motif_annot join (M0..M8 Y axis)
# --level l1p5:   group by nhood_label_l1p5 from the DA CSV directly
#
# Substrate:
#   --contrasts-csv  columns: contrast_label, da_csv_path, diagram_path
#   --nhood-motif    per_nhood_motif_annot (only needed if --level motif)

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(png)
  library(grid)
  library(cowplot)
})
have_arrow    <- requireNamespace("arrow",      quietly = TRUE)
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)
have_rastr    <- requireNamespace("ggrastr",    quietly = TRUE)

parser <- ArgumentParser()
parser$add_argument("--contrasts-csv", required = TRUE)
parser$add_argument("--nhood-motif", default = "",
                    help = "Required if --level motif")
parser$add_argument("--publication-root", required = TRUE)
parser$add_argument("--root", default = "")
parser$add_argument("--out-dir", required = TRUE)
parser$add_argument("--panel-id",
                    default = "s_xenium_beeswarm_grid_with_diagrams")
parser$add_argument("--level", choices = c("motif", "l1p5", "l2s"), required = TRUE)
parser$add_argument("--label-csv", default = "",
                    help = "Required if --level l2s: nhood_id → dominant_label (e.g., nhood_l2s_full.csv)")
parser$add_argument("--platform-csv", default = "",
                    help = "Optional (l1p5 only): per-nhood dominant_platform; triggers platform-dodge")
parser$add_argument("--motifs", default = "M0,M1,M2,M3,M4,M5,M6,M7,M8")
parser$add_argument("--add-epi-rows", action = "store_true",
                    help = "When --level motif, also include epi L1.5 rows (LASP-basal/LASP/LHS/BMYO-myo) sourced from nhood_label_l1p5")
parser$add_argument("--l1p5-frac-min", type = "double", default = 0.7)
parser$add_argument("--fdr-floor", type = "double", default = 1e-15)
parser$add_argument("--width-per-contrast-in", type = "double", default = 1.1)
parser$add_argument("--width-legend-in", type = "double", default = 1.2)
parser$add_argument("--per-row-height-in", type = "double", default = 0.10)
parser$add_argument("--diagram-row-height-in", type = "double", default = 1.0,
                    help = "Height of the diagram header strip")
parser$add_argument("--no-diagrams", action = "store_true",
                    help = "Skip the diagram header strip (main-panel mode)")
parser$add_argument("--height-overhead-in", type = "double", default = 1.4)
args <- parser$parse_args()

source(file.path(args$publication_root, "publication", "config",
                 "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(args$publication_root, "publication",
                                     "config"))
panel_theme <- get_theme(aes_cfg)
div_scale   <- aes_cfg$scales$expression_diverging

dir.create(args$out_dir, showWarnings = FALSE, recursive = TRUE)

motif_order <- strsplit(args$motifs, ",", fixed = TRUE)[[1]] |> trimws()

# --- Load contrasts CSV ---
contrasts <- read_csv(args$contrasts_csv, show_col_types = FALSE)
cat(sprintf("[1] %d contrasts to render (level=%s)\n",
            nrow(contrasts), args$level))

resolve_path <- function(p) {
  if (substr(p, 1, 1) == "/") p else file.path(args$root, p)
}

# --- Build long-form DA depending on level ---
if (args$level == "motif") {
  if (args$nhood_motif == "") stop("--nhood-motif required for --level motif")
  nhood_motif <- if (grepl("\\.parquet$", args$nhood_motif)) {
    if (!have_arrow) stop("arrow required for parquet")
    arrow::read_parquet(args$nhood_motif)
  } else {
    read_csv(args$nhood_motif, show_col_types = FALSE)
  } |> select(nhood_idx, dominant_motif)

  EPI_L1P5 <- c("LASP-basal", "LASP", "LHS", "BMYO-myo")
  da_all <- contrasts |>
    rowwise() |>
    do({
      p <- resolve_path(.$da_csv_path)
      if (!file.exists(p)) stop("missing DA: ", p)
      da_raw <- read_csv(p, show_col_types = FALSE) |>
        filter(!is.na(SpatialFDR), !is.na(logFC)) |>
        mutate(nhood_label_l1p5 = dplyr::recode(nhood_label_l1p5,
                                                  "BMYO" = "BMYO-myo"))
      # Non-epi motif rows
      non_epi <- da_raw |>
        inner_join(nhood_motif, by = c("Nhood" = "nhood_idx")) |>
        filter(dominant_motif %in% motif_order) |>
        rename(group = dominant_motif) |>
        mutate(contrast = .$contrast_label,
               row_class = "motif")
      out <- non_epi
      if (args$add_epi_rows) {
        epi <- da_raw |>
          filter(!is.na(nhood_label_l1p5),
                 nhood_label_l1p5 %in% EPI_L1P5,
                 nhood_label_l1p5_frac >= args$l1p5_frac_min) |>
          rename(group = nhood_label_l1p5) |>
          mutate(contrast = .$contrast_label,
                 row_class = "epi") |>
          select(any_of(names(non_epi)))
        out <- bind_rows(non_epi, epi)
      }
      out
    }) |>
    ungroup() |>
    as_tibble()
  if (args$add_epi_rows) {
    # Epi rows at top (drawn first because ggplot reverses y-axis factor),
    # motifs below — combined order top→bottom = epi → motifs
    group_levels <- c(EPI_L1P5, motif_order)
  } else {
    group_levels <- motif_order
  }
  y_label <- ""
  per_row_h <- args$per_row_height_in
  axis_text_y_size <- 6
} else if (args$level == "l1p5") {
  da_all <- contrasts |>
    rowwise() |>
    do({
      p <- resolve_path(.$da_csv_path)
      if (!file.exists(p)) stop("missing DA: ", p)
      read_csv(p, show_col_types = FALSE) |>
        filter(!is.na(SpatialFDR), !is.na(logFC),
               !is.na(nhood_label_l1p5),
               nhood_label_l1p5_frac >= args$l1p5_frac_min) |>
        rename(group = nhood_label_l1p5) |>
        mutate(contrast = .$contrast_label)
    }) |>
    ungroup() |>
    as_tibble()
  # Optional platform-dodge: join per-nhood dominant_platform
  if (args$platform_csv != "") {
    nhood_plat <- read_csv(resolve_path(args$platform_csv),
                             show_col_types = FALSE) |>
      select(Nhood, dominant_platform)
    da_all <- da_all |>
      left_join(nhood_plat, by = "Nhood") |>
      filter(!is.na(dominant_platform))
    cat(sprintf("[1b] platform join: %d nhoods (%d xenium-dom, %d flex-dom)\n",
                nrow(da_all),
                sum(da_all$dominant_platform == "xenium"),
                sum(da_all$dominant_platform == "flex")))
  }
  # Match dotplot order (render_fig3e_xenium_dotplot_nuclear.R:74-82).
  # Remap C1's older vocabulary (BMYO → BMYO-myo) so labels align across contrasts.
  da_all$group <- dplyr::recode(da_all$group, "BMYO" = "BMYO-myo")
  L1P5_ORDER <- c(
    "BMYO-myo", "LASP", "LASP-basal", "LHS",
    "Fibroblast", "Fibroblast_activated", "Fibroblast_SFRP4",
    "Endothelial", "Vas-capillary", "Lymphatic Endothelial", "Pericyte", "Adipocyte",
    "Macrophage", "cDC", "cDC1", "cDC2", "pDC",
    "Mast cell", "Neutrophil",
    "CD4 T cell", "CD4 Treg", "CD8 T cell", "T-NK",
    "B cell", "Plasma cell")
  present <- unique(da_all$group)
  group_levels <- c(intersect(L1P5_ORDER, present),
                    setdiff(present, L1P5_ORDER))
  y_label <- ""
  per_row_h <- 0.08
  axis_text_y_size <- 5.5
} else {
  # l2s: join nhood_l2s_full.csv on Nhood == nhood_id
  if (args$label_csv == "") stop("--label-csv required for --level l2s")
  nhood_l2s <- read_csv(resolve_path(args$label_csv),
                         show_col_types = FALSE) |>
    select(nhood_id, dominant_label, label_purity)
  da_all <- contrasts |>
    rowwise() |>
    do({
      p <- resolve_path(.$da_csv_path)
      if (!file.exists(p)) stop("missing DA: ", p)
      read_csv(p, show_col_types = FALSE) |>
        filter(!is.na(SpatialFDR), !is.na(logFC)) |>
        inner_join(nhood_l2s, by = c("Nhood" = "nhood_id")) |>
        filter(!is.na(dominant_label),
               !grepl("^ARTIFACT_", dominant_label),
               label_purity >= args$l1p5_frac_min) |>
        rename(group = dominant_label) |>
        mutate(contrast = .$contrast_label)
    }) |>
    ungroup() |>
    as_tibble()
  # L2.0s ordering: epi → stromal → endothelial → adipo → immune
  L2S_ORDER <- c(
    # Epi luminal
    "LHS-major", "LHS-apocrine", "LASP-major", "LASP-basal",
    # Epi basal
    "BMYO-myo", "BMYO-basal",
    # Stromal fibroblasts
    "Fibro-major", "Fibro-myo", "Fibro-perivascular",
    "Fibro-prematrix", "Fibro-SFRP4",
    # Adipo
    "Adipo-storage", "Adipo-lipogenic", "Adipo-lipolysis",
    # Endothelial / vascular
    "EC", "EC-cap", "EC-art", "EC-vein", "LEC", "PV",
    # Immune T / NK
    "CD4_Th_like", "CD8_Trm", "CD8_Tem", "CD8_Resting", "NK", "Treg",
    # Dendritic
    "cDC1", "cDC2", "pDC",
    # B / Plasma
    "B_cell", "Plasma_cell",
    # Macrophage / mast / neutrophil
    "Macrophage", "Mac_other", "Mast", "Neutrophil")
  present <- unique(da_all$group)
  group_levels <- c(intersect(L2S_ORDER, present),
                    setdiff(present, L2S_ORDER))
  y_label <- ""
  per_row_h <- 0.08
  axis_text_y_size <- 5.5
}
cat(sprintf("[2] loaded %d nhoods\n", nrow(da_all)))

# --- Color: signed −log10(SpatialFDR), shared scale ---
da_all <- da_all |>
  mutate(neglog_fdr = -log10(pmax(SpatialFDR, args$fdr_floor)),
         signed_sig = sign(logFC) * neglog_fdr)
sig_cap <- max(quantile(abs(da_all$signed_sig), 0.99, na.rm = TRUE), 1.0)

da_all$group <- factor(da_all$group, levels = rev(group_levels))
da_all$contrast <- factor(da_all$contrast, levels = contrasts$contrast_label)

# --- Subsample ns ---
NS_SIG_THRESHOLD <- 0.10
NS_SAMPLE_CAP <- if (args$level == "motif") 1200 else 600
set.seed(42)
sub_sig <- da_all |> filter(SpatialFDR <  NS_SIG_THRESHOLD)
sub_ns  <- da_all |> filter(SpatialFDR >= NS_SIG_THRESHOLD) |>
  group_by(contrast, group) |> slice_sample(n = NS_SAMPLE_CAP) |> ungroup()
da_plot <- bind_rows(sub_sig, sub_ns) |>
  arrange(contrast, group, abs(signed_sig))
cat(sprintf("[3] subsample kept %d / %d\n", nrow(da_plot), nrow(da_all)))

x_lim <- max(max(abs(da_plot$logFC), na.rm = TRUE) * 1.05, 1)

use_platform_dodge <- args$platform_csv != "" && args$level == "l1p5"
swarm_geom <- if (have_beeswarm) {
  if (use_platform_dodge) {
    ggbeeswarm::geom_quasirandom(
      aes(x = logFC, y = group, color = signed_sig,
          group = dominant_platform),
      groupOnX = FALSE, size = 0.10, alpha = 0.75,
      bandwidth = 0.35, width = 0.135, dodge.width = 0.30)
  } else {
    ggbeeswarm::geom_quasirandom(
      aes(x = logFC, y = group, color = signed_sig),
      groupOnX = FALSE, size = 0.10, alpha = 0.75,
      bandwidth = 0.35, width = 0.27)
  }
} else {
  geom_jitter(aes(x = logFC, y = group, color = signed_sig),
              width = 0.0, height = 0.4, size = 0.10, alpha = 0.75)
}
swarm_layer <- if (have_rastr) {
  ggrastr::rasterise(swarm_geom, dpi = 600, dev = "cairo")
} else {
  swarm_geom
}

violin_layer <- if (use_platform_dodge) {
  list(
    geom_violin(
      aes(x = logFC, y = group, linetype = dominant_platform,
          group = interaction(group, dominant_platform)),
      orientation = "y",
      position = position_dodge(width = 0.30),
      fill = NA, color = "black", linewidth = 0.30,
      scale = "width", width = 0.27),
    scale_linetype_manual(values = c(xenium = "solid", flex = "dashed"),
                            name = "platform")
  )
} else {
  geom_violin(
    aes(x = logFC, y = group, group = group),
    orientation = "y",
    fill = NA, color = "black", linewidth = 0.30,
    scale = "width", width = 0.55, linetype = "solid")
}

# --- Beeswarm panel (no strip text; diagrams will sit above) ---
beeswarm <- ggplot(da_plot) +
  geom_vline(xintercept = 0, color = "grey60",
             linewidth = 0.3, linetype = "dashed") +
  swarm_layer +
  violin_layer +
  scale_color_gradient2(low = div_scale$low, mid = "grey90",
                        high = div_scale$high,
                        midpoint = 0, limits = c(-sig_cap, sig_cap),
                        oob = scales::squish,
                        name = "signed-log10(FDR)") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(xlim = c(-x_lim, x_lim), clip = "off") +
  facet_grid(. ~ contrast, switch = "y") +
  labs(x = "logFC (per nhood)", y = y_label) +
  panel_theme +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),  # diagrams replace text strip
        legend.position  = "right",
        legend.direction = "vertical",
        legend.title     = element_text(size = 6, angle = 90, hjust = 0.5),
        legend.title.position = "left",
        legend.key.height = unit(0.9, "cm"),
        legend.key.width  = unit(0.30, "cm"),
        panel.grid       = element_blank(),
        panel.background = element_blank(),
        axis.text.y      = element_text(size = axis_text_y_size,
                                          face = "bold"),
        axis.text.x      = element_text(size = 5),
        axis.title.x     = element_text(size = 7),
        axis.ticks.y     = element_blank(),
        panel.spacing.x  = unit(0.4, "lines"),
        plot.margin      = margin(2, 8, 4, 4))

if (args$no_diagrams) {
  # Main-panel mode: show strip text labels above each contrast (no diagrams)
  beeswarm <- beeswarm + theme(strip.text.x =
                                element_text(face = "bold", size = 7))
  full <- beeswarm
  n_groups <- length(group_levels)
  n_contrasts <- nrow(contrasts)
  pdf_width  <- n_contrasts * args$width_per_contrast_in + args$width_legend_in
  pdf_height <- n_groups * per_row_h + args$height_overhead_in
  pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
  ggsave(pdf_path, plot = full, width = pdf_width, height = pdf_height,
         units = "in")
  cat(sprintf("Wrote %s (%.2f x %.2f in; %d contrasts × %d groups; main mode)\n",
              pdf_path, pdf_width, pdf_height, n_contrasts, n_groups))
  quit(save = "no")
}

# --- Diagram strip: one cowplot::ggdraw per contrast, then plot_grid ---
diagram_panels <- lapply(contrasts$diagram_path, function(dp) {
  p <- resolve_path(dp)
  if (!file.exists(p)) {
    cat(sprintf("WARN missing diagram: %s\n", p))
    return(ggdraw())
  }
  img <- png::readPNG(p, native = TRUE)
  ggdraw() +
    draw_image(img, x = 0, y = 0, width = 1, height = 1) +
    theme(plot.margin = margin(0, 0, 0, 0))
})
diagram_strip <- do.call(plot_grid,
                          c(diagram_panels, list(nrow = 1,
                                                   align = "h",
                                                   rel_widths =
                                                     rep(1, nrow(contrasts)))))
# Pad diagram_strip on the right to leave room for the beeswarm legend
diagram_strip_padded <- plot_grid(diagram_strip, NULL, nrow = 1,
                                   rel_widths = c(
                                     nrow(contrasts) * args$width_per_contrast_in,
                                     args$width_legend_in))

# --- Composite ---
full <- plot_grid(diagram_strip_padded, beeswarm, ncol = 1,
                   rel_heights = c(args$diagram_row_height_in,
                                   length(group_levels) * per_row_h +
                                     args$height_overhead_in))

n_groups <- length(group_levels)
n_contrasts <- nrow(contrasts)
pdf_width  <- n_contrasts * args$width_per_contrast_in + args$width_legend_in
pdf_height <- args$diagram_row_height_in + n_groups * per_row_h + args$height_overhead_in
pdf_path <- file.path(args$out_dir, paste0(args$panel_id, ".pdf"))
ggsave(pdf_path, plot = full, width = pdf_width, height = pdf_height,
       units = "in")
cat(sprintf("Wrote %s (%.2f x %.2f in; %d contrasts × %d groups; level=%s)\n",
            pdf_path, pdf_width, pdf_height, n_contrasts, n_groups,
            args$level))
