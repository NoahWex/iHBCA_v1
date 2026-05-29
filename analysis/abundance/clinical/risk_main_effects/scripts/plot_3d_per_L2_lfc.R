#!/usr/bin/env Rscript
# plot_3d_per_L2_lfc.R
# 3D scatter: x=AR_med_lfc, y=BRCA1_med_lfc, z=sporadic_med_lfc per L2.
# Colored by compartment (L0). Static PDF (scatterplot3d) + interactive HTML (plotly).

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})
have_s3d <- requireNamespace("scatterplot3d", quietly = TRUE)
have_plotly <- requireNamespace("plotly", quietly = TRUE)
have_htmlwidgets <- requireNamespace("htmlwidgets", quietly = TRUE)

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_3d_per_L2_lfc.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

e <- read_csv(paths$outputs$stageE, show_col_types = FALSE)

df <- e %>% transmute(
  L2_joint, compartment, label,
  AR        = as.numeric(parity_in_AR__med_lfc),
  BRCA1     = as.numeric(parity_x_HR_BRCA1__med_lfc),
  sporadic  = as.numeric(parity_x_HR_sporadic__med_lfc),
  n_AR      = as.integer(parity_in_AR__n),
  n_sig_AR  = as.integer(parity_in_AR__n_sig),
  n_sig_BRCA1 = as.integer(parity_x_HR_BRCA1__n_sig),
  n_sig_sporadic = as.integer(parity_x_HR_sporadic__n_sig)) %>%
  filter(!is.na(AR), !is.na(BRCA1), !is.na(sporadic))

cat(sprintf("L2 entities: %d\n", nrow(df)))
cat("By compartment:\n")
print(df %>% count(compartment))

comp_colors <- c(epi = "#E41A1C", imm = "#377EB8", str = "#4DAF4A")

# ---- Static PDF via scatterplot3d ----
if (have_s3d) {
  pdf_path <- file.path(fig_dir, "headline_3d_lfc_AR_BRCA1_sporadic.pdf")
  pdf(pdf_path, width = 8, height = 8)
  s3d <- scatterplot3d::scatterplot3d(
    x = df$AR, y = df$BRCA1, z = df$sporadic,
    pch = 19, cex.symbols = 1.4,
    color = comp_colors[df$compartment],
    xlab = "parity_in_AR  median lfc",
    ylab = "parity_x_HR_BRCA1  median lfc",
    zlab = "parity_x_HR_sporadic  median lfc",
    main = "Per-L2 parity median LFC across 3 strata",
    angle = 55, box = FALSE, grid = TRUE,
    xlim = c(-1, 2), ylim = c(-1, 2), zlim = c(-1, 2))
  # Reference plane at zero on each axis (origin)
  s3d$points3d(c(0,0), c(0,0), c(-2,2), type = "l", lty = 2, col = "grey60")
  s3d$points3d(c(-2,2), c(0,0), c(0,0), type = "l", lty = 2, col = "grey60")
  s3d$points3d(c(0,0), c(-2,2), c(0,0), type = "l", lty = 2, col = "grey60")
  legend("topright",
          legend = names(comp_colors),
          col    = comp_colors,
          pch = 19, bty = "n", cex = 1.1)
  # Label highlight points
  hl <- df %>% filter(L2_joint %in% c("epi::BMYO-basal", "str::Fibro-SFRP4",
                                        "epi::Lactocyte-LC2", "epi::LASP-basal",
                                        "imm::Plasma", "imm::CD8_Trm",
                                        "str::Fibro-IGF1", "str::Fibro-prematrix",
                                        "imm::NK", "imm::Treg"))
  if (nrow(hl) > 0) {
    proj <- s3d$xyz.convert(hl$AR, hl$BRCA1, hl$sporadic)
    text(proj$x, proj$y, labels = hl$label, pos = 4, cex = 0.55,
          col = "black")
  }
  dev.off()
  cat(sprintf("Wrote: %s\n", pdf_path))
}

# ---- Interactive HTML via plotly ----
if (have_plotly && have_htmlwidgets) {
  hover_text <- sprintf(
    "%s<br>AR=%+.2f BRCA1=%+.2f spor=%+.2f<br>n_sig: AR=%d  BRCA1=%d  spor=%d",
    df$L2_joint, df$AR, df$BRCA1, df$sporadic,
    df$n_sig_AR %||% 0, df$n_sig_BRCA1 %||% 0, df$n_sig_sporadic %||% 0)
  p <- plotly::plot_ly(
    x = df$AR, y = df$BRCA1, z = df$sporadic,
    text = hover_text, hoverinfo = "text",
    type = "scatter3d", mode = "markers+text",
    marker = list(size = 6, color = comp_colors[df$compartment],
                    line = list(color = "black", width = 0.5)),
    textfont = list(size = 9)) %>%
    plotly::layout(
      title = "Per-L2 parity median LFC across 3 strata",
      scene = list(
        xaxis = list(title = "AR baseline"),
        yaxis = list(title = "BRCA1 modifier"),
        zaxis = list(title = "Sporadic modifier")))
  html_path <- file.path(fig_dir, "headline_3d_lfc_AR_BRCA1_sporadic.html")
  htmlwidgets::saveWidget(p, html_path, selfcontained = TRUE)
  cat(sprintf("Wrote: %s\n", html_path))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Output the table for reference
write_csv(df, file.path(paths$inquiry_root, "outputs", "headline_3d_lfc_table.csv"))
cat("Wrote headline_3d_lfc_table.csv\n")

cat("\n=== 3D plot done ===\n")
