#!/usr/bin/env Rscript
# plot_AR_hero_beeswarm.R
# AR baseline parity main-effects hero panel for Slide 2.
# Beeswarm of per-nhood logFC for parity_in_AR (and parity_in_AR_with_twigger
# when available) faceted by parent_L2, ordered by AR median lfc within
# compartment. Highlight populations called out in the slide narrative.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_AR_hero_beeswarm.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

# Headline populations for Slide 2 narrative — highlighted in plot.
HIGHLIGHTED_L2 <- c(
  # AR parity main effects (per Kai narrative)
  "epi::Lactocyte-LC2",   "epi::LASP-basal",
  "imm::Plasma",          "imm::CD8_Trm",
  "imm::NK",              "imm::Macro_FOLR2",
  "str::Fibro-IGF1",      "str::Fibro-prematrix",
  "str::Fibro-SFRP4"
)

contrast_dirs <- list.dirs(paths$outputs$stageD, recursive = FALSE)
ar_contrasts <- grep("^parity_in_AR", basename(contrast_dirs), value = TRUE)
cat(sprintf("Found AR contrasts: %s\n", paste(ar_contrasts, collapse = ", ")))

all_da <- list()
for (cn in ar_contrasts) {
  da_csv <- file.path(paths$outputs$stageD, cn, "da_results.csv")
  if (!file.exists(da_csv)) next
  da <- read_csv(da_csv, show_col_types = FALSE)
  if (!"compartment" %in% colnames(da)) {
    if ("L2_compartment" %in% colnames(da)) da$compartment <- da$L2_compartment
  }
  if (!"label" %in% colnames(da)) {
    if ("L2_label" %in% colnames(da)) da$label <- da$L2_label
  }
  da$L2_joint <- paste(da$compartment, da$label, sep = "::")
  da$contrast <- cn
  da$is_sig <- !is.na(da$SpatialFDR) & da$SpatialFDR < 0.05
  da$direction <- case_when(
    da$is_sig & da$logFC > 0 ~ "sig_up",
    da$is_sig & da$logFC < 0 ~ "sig_dn",
    TRUE ~ "ns"
  )
  all_da[[cn]] <- da
}
df <- bind_rows(all_da)
cat(sprintf("Total nhoods loaded: %d across %d contrasts\n",
            nrow(df), length(unique(df$contrast))))

# Order L2s within compartment by AR median lfc
order_l2 <- df %>% filter(contrast == "parity_in_AR") %>%
  group_by(compartment, L2_joint) %>%
  summarise(med_lfc = median(logFC, na.rm = TRUE), .groups = "drop") %>%
  arrange(compartment, desc(med_lfc))
df$L2_joint <- factor(df$L2_joint, levels = order_l2$L2_joint)
df$is_highlight <- df$L2_joint %in% HIGHLIGHTED_L2

palette_dir <- c(sig_up = "#d62728", sig_dn = "#1f77b4", ns = "grey75")

# ---- Render: AR-only hero ----
ar_only <- df %>% filter(contrast == "parity_in_AR")

p_hero <- ggplot(ar_only,
                  aes(x = logFC, y = L2_joint, color = direction)) +
  geom_vline(xintercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.3)
p_hero <- if (have_beeswarm) {
  p_hero + ggbeeswarm::geom_quasirandom(size = 0.4, alpha = 0.6,
                                          groupOnX = FALSE, dodge.width = 0.6)
} else {
  p_hero + geom_jitter(size = 0.4, alpha = 0.6, height = 0.25)
}
p_hero <- p_hero +
  scale_color_manual(values = palette_dir) +
  facet_grid(compartment ~ ., scales = "free_y", space = "free_y") +
  labs(x = "logFC (parity, AR baseline)", y = NULL,
        color = "DA direction (FDR<0.05)") +
  theme_bw(base_size = 8) +
  theme(strip.text = element_text(size = 8.5, face = "bold"),
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey95"),
        axis.text.y = element_text(size = 7))

# Highlight overlay: black border swarm on highlighted L2s
hl_df <- ar_only %>% filter(is_highlight)
if (nrow(hl_df) > 0 && have_beeswarm) {
  p_hero <- p_hero +
    ggbeeswarm::geom_quasirandom(data = hl_df,
                                   size = 0.5, alpha = 0.85,
                                   color = "black", shape = 1,
                                   groupOnX = FALSE, dodge.width = 0.6,
                                   stroke = 0.3,
                                   inherit.aes = FALSE,
                                   aes(x = logFC, y = L2_joint))
}

out_dir <- file.path(paths$inquiry_root, "share", "kai_v1", "04_panels")
if (!dir.exists(out_dir))
  out_dir <- file.path(paths$outputs$reports, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_l2 <- nlevels(droplevels(ar_only$L2_joint))
height_in <- max(5, 0.16 * n_l2 + 1)
ggsave(file.path(out_dir, "slide2_AR_hero_beeswarm.pdf"),
        p_hero, width = 7, height = height_in, device = cairo_pdf)
cat(sprintf("Wrote: %s/slide2_AR_hero_beeswarm.pdf (%d L2s, h=%.1f in)\n",
            out_dir, n_l2, height_in))

# ---- If twigger comparison contrast exists, side-by-side render ----
if ("parity_in_AR_with_twigger" %in% ar_contrasts) {
  pair_df <- df %>% filter(contrast %in% c("parity_in_AR",
                                              "parity_in_AR_with_twigger"))
  pair_df$contrast <- factor(pair_df$contrast,
                              levels = c("parity_in_AR",
                                          "parity_in_AR_with_twigger"))
  p_pair <- ggplot(pair_df,
                     aes(x = logFC, y = L2_joint, color = direction))
  p_pair <- p_pair + geom_vline(xintercept = 0, linetype = "dashed",
                                  color = "grey40", linewidth = 0.3)
  p_pair <- if (have_beeswarm) {
    p_pair + ggbeeswarm::geom_quasirandom(size = 0.35, alpha = 0.55,
                                            groupOnX = FALSE)
  } else {
    p_pair + geom_jitter(size = 0.35, alpha = 0.55, height = 0.25)
  }
  p_pair <- p_pair +
    scale_color_manual(values = palette_dir) +
    facet_grid(compartment ~ contrast, scales = "free_y", space = "free_y") +
    labs(x = "logFC (parity)", y = NULL, color = "DA direction (FDR<0.05)") +
    theme_bw(base_size = 8) +
    theme(strip.text = element_text(size = 8.5, face = "bold"),
          legend.position = "bottom",
          axis.text.y = element_text(size = 7))
  ggsave(file.path(out_dir, "slide2_AR_vs_AR_twigger_beeswarm.pdf"),
          p_pair, width = 11, height = height_in, device = cairo_pdf)
  cat(sprintf("Wrote: %s/slide2_AR_vs_AR_twigger_beeswarm.pdf\n", out_dir))
}

cat("\n=== plot_AR_hero_beeswarm done ===\n")
