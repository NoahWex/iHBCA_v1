#!/usr/bin/env Rscript
# plot_headline_panels.R
# Curated panels for the two headline stories (BRCA1 BMYO-basal expansion +
# AR/BRCA1 Fibro-SFRP4 stratum-flip).
#
# Outputs (per-headline focused panels for share package):
#   reports/figures/headline_bmyo_basal_3contrast.pdf
#   reports/figures/headline_fibro_sfrp4_3contrast.pdf
#   reports/figures/headline_population_directionality.pdf

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr); library(stringr)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_headline_panels.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1", "parity_x_HR_sporadic")
CONTRAST_LABELS <- c(
  parity_in_AR        = "AR baseline",
  parity_x_HR_BRCA1   = "BRCA1 carriers",
  parity_x_HR_sporadic = "Sporadic HR")

# Aesthetic — restrained, presentation-clean
COLORS <- c(up = "#B2182B", dn = "#2166AC", ns = "grey75")
HIGHLIGHT_COLOR <- "#1A5276"   # frame for highlighted niches
ALPHA_NS <- 0.20
ALPHA_SIG <- 0.85

VIABLE_MIN <- 10

# Headline highlights — niches we call out by name in panel
HIGHLIGHTS <- list(
  `epi::BMYO-basal` = list(
    parity_x_HR_BRCA1 = c("epi::BMYO-basal_11", "epi::BMYO-basal_33",
                            "epi::BMYO-basal_75"),
    parity_in_AR        = c(),
    parity_x_HR_sporadic = c()),
  `str::Fibro-SFRP4` = list(
    parity_in_AR        = c("str::Fibro-SFRP4_64"),
    parity_x_HR_BRCA1   = c("str::Fibro-SFRP4_48"),
    parity_x_HR_sporadic = c())
)

# ------------------------------------------------------------------------
# Load Stage D + F.1 nhoodgroup membership for the headline L2s
# ------------------------------------------------------------------------
load_per_contrast <- function(L2_target) {
  out <- list()
  for (cn in CONTRASTS) {
    da_csv <- file.path(paths$outputs$stageD, cn, "da_results.csv")
    ng_csv <- file.path(paths$outputs$stageF1, cn, "nhood_groups.csv")
    if (!file.exists(da_csv) || !file.exists(ng_csv)) next
    d  <- read_csv(da_csv, show_col_types = FALSE,
                   col_types = cols(.default = "?"))
    nm <- read_csv(ng_csv, show_col_types = FALSE,
                   col_types = cols(.default = "c")) %>%
      mutate(Nhood = as.integer(Nhood),
              n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group)))
    sub <- d %>% filter(label == sub("^[a-z]+::", "", L2_target)) %>%
      left_join(nm %>% select(Nhood, NhoodGroup_renamed,
                                  n_nhoods_in_group),
                by = "Nhood") %>%
      mutate(contrast = cn,
              sig = !is.na(SpatialFDR) & SpatialFDR < 0.05,
              direction = ifelse(sig & logFC > 0, "up",
                                 ifelse(sig & logFC < 0, "dn", "ns")),
              viable = !is.na(n_nhoods_in_group) & n_nhoods_in_group >= VIABLE_MIN,
              ng_label = ifelse(viable, NhoodGroup_renamed, NA_character_))
    out[[cn]] <- sub
  }
  bind_rows(out)
}

bmyo <- load_per_contrast("epi::BMYO-basal")
sfrp4 <- load_per_contrast("str::Fibro-SFRP4")

cat(sprintf("BMYO-basal: %d nhoods loaded across %d contrasts\n",
            nrow(bmyo), length(unique(bmyo$contrast))))
cat(sprintf("Fibro-SFRP4: %d nhoods loaded across %d contrasts\n",
            nrow(sfrp4), length(unique(sfrp4$contrast))))

# ------------------------------------------------------------------------
# Per-L2 beeswarm panel: each contrast as a row, NhoodGroups as column-bins
# X-axis: NhoodGroup (ordered by median LFC); Y-axis: per-nhood logFC
# Each viable group gets a column with all its nhoods as dots.
# Non-viable nhoods pooled into 'other' column.
# ------------------------------------------------------------------------
build_l2_panel <- function(d, L2_target, hlist, out_pdf, height = 6.5) {
  # Order NhoodGroups within each contrast by median group LFC
  d2 <- d %>%
    mutate(contrast_label = factor(CONTRAST_LABELS[contrast],
                                    levels = CONTRAST_LABELS[CONTRASTS]))
  # Build column factor: viable groups labeled, non-viable as 'other'
  d2 <- d2 %>%
    mutate(col_label = ifelse(viable,
                               sub(paste0("^", L2_target, "_"), "_",
                                    NhoodGroup_renamed),
                               "other"))

  # Within each contrast, order columns by median LFC
  col_order <- d2 %>% group_by(contrast, col_label) %>%
    summarise(med = median(logFC, na.rm = TRUE),
              n = n(), .groups = "drop") %>%
    arrange(contrast, ifelse(col_label == "other", 0, 1), desc(med))

  d2 <- d2 %>% mutate(
    col_label = factor(col_label, levels = unique(col_order$col_label)))

  # Highlight flag: niches called out in headlines
  d2 <- d2 %>% rowwise() %>%
    mutate(highlight = !is.na(NhoodGroup_renamed) &
                       NhoodGroup_renamed %in%
                         (hlist[[as.character(contrast)]] %||% character(0))) %>%
    ungroup()

  # Row counts per contrast for subtitle
  ncont_summary <- d2 %>%
    group_by(contrast_label) %>%
    summarise(n_nhoods = n(),
              n_viable = sum(viable),
              pct_pos = sprintf("%.0f%%", 100 * mean(logFC > 0, na.rm = TRUE)),
              n_sig_up = sum(sig & logFC > 0, na.rm = TRUE),
              n_sig_dn = sum(sig & logFC < 0, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(subtitle = sprintf("n=%d nhoods | %s pos | sig_up:dn = %d:%d | %d viable groups",
                                n_nhoods, pct_pos, n_sig_up, n_sig_dn, n_viable))

  # Median tick per column
  med_per_col <- d2 %>% group_by(contrast_label, col_label) %>%
    summarise(med = median(logFC, na.rm = TRUE), .groups = "drop")

  swarm <- if (have_beeswarm) {
    ggbeeswarm::geom_quasirandom(
      aes(x = col_label, y = logFC, colour = direction, alpha = sig),
      size = 0.7, bandwidth = 0.4, width = 0.4)
  } else {
    geom_jitter(aes(x = col_label, y = logFC, colour = direction, alpha = sig),
                width = 0.25, height = 0, size = 0.7)
  }

  p <- ggplot(d2) +
    swarm +
    geom_point(data = med_per_col,
                aes(x = col_label, y = med),
                shape = "|", size = 6, colour = "black", inherit.aes = FALSE) +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
    facet_wrap(~ contrast_label, ncol = 1, scales = "free_x") +
    scale_colour_manual(values = COLORS, name = "direction") +
    scale_alpha_manual(values = c(`TRUE` = ALPHA_SIG, `FALSE` = ALPHA_NS),
                        guide = "none") +
    coord_cartesian(ylim = c(-4, 4)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          strip.text = element_text(face = "bold", size = 10),
          legend.position = "top") +
    labs(x = NULL, y = "logFC (per nhood, parous vs nullip)",
          title = sprintf("Headline: %s — per-NhoodGroup parity response across germline strata",
                            L2_target),
          subtitle = sprintf("Each column = one viable NhoodGroup (n>=%d) plus 'other' (non-viable nhoods pooled). Black tick = column median LFC. Highlighted niches outlined in caption.",
                              VIABLE_MIN))

  # Annotate highlighted columns with rectangle + label
  highlight_rows <- d2 %>% filter(highlight) %>%
    distinct(contrast_label, col_label, NhoodGroup_renamed) %>%
    group_by(contrast_label, col_label) %>%
    summarise(label = first(NhoodGroup_renamed), .groups = "drop")

  if (nrow(highlight_rows) > 0) {
    p <- p +
      geom_text(data = highlight_rows %>%
                  left_join(med_per_col, by = c("contrast_label", "col_label")),
                aes(x = col_label, y = pmax(med, 2.5) + 0.5,
                    label = "*"),
                size = 5, colour = HIGHLIGHT_COLOR, inherit.aes = FALSE)
  }

  ggsave(out_pdf, p, width = 11, height = height)
  cat(sprintf("Wrote: %s\n", out_pdf))
  list(plot = p, summary = ncont_summary)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

bmyo_out <- build_l2_panel(bmyo, "epi::BMYO-basal",
                            HIGHLIGHTS$`epi::BMYO-basal`,
                            file.path(fig_dir,
                                       "headline_bmyo_basal_3contrast.pdf"),
                            height = 7)
sfrp4_out <- build_l2_panel(sfrp4, "str::Fibro-SFRP4",
                             HIGHLIGHTS$`str::Fibro-SFRP4`,
                             file.path(fig_dir,
                                        "headline_fibro_sfrp4_3contrast.pdf"),
                             height = 7)

cat("\nBMYO-basal summary:\n"); print(bmyo_out$summary)
cat("\nFibro-SFRP4 summary:\n"); print(sfrp4_out$summary)

# ------------------------------------------------------------------------
# Headline directionality summary panel: 2 populations × 3 contrasts
# Tile: % pos as fill + n_sig_up vs n_sig_dn as text annotation
# ------------------------------------------------------------------------
dir_data <- bind_rows(
  bmyo  %>% mutate(L2 = "epi::BMYO-basal"),
  sfrp4 %>% mutate(L2 = "str::Fibro-SFRP4")) %>%
  group_by(L2, contrast) %>%
  summarise(n = n(),
             pct_pos = 100 * mean(logFC > 0, na.rm = TRUE),
             n_sig_up = sum(sig & logFC > 0, na.rm = TRUE),
             n_sig_dn = sum(sig & logFC < 0, na.rm = TRUE),
             med_lfc = median(logFC, na.rm = TRUE),
             .groups = "drop") %>%
  mutate(contrast_label = factor(CONTRAST_LABELS[contrast],
                                   levels = CONTRAST_LABELS[CONTRASTS]),
          # Centered diverging fill: 50% = neutral
          fill_val = pct_pos - 50)

p_dir <- ggplot(dir_data,
                 aes(x = contrast_label, y = L2, fill = fill_val)) +
  geom_tile(colour = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("%.0f%% pos\nsig %d : %d\nmed %+.2f",
                                  pct_pos, n_sig_up, n_sig_dn, med_lfc)),
             size = 3.2, colour = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-25, 25),
                        oob = scales::squish,
                        name = "% pos − 50",
                        breaks = c(-20, 0, 20),
                        labels = c("80% neg", "balanced", "80% pos")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = "right",
        plot.title = element_text(face = "bold")) +
  labs(x = NULL, y = NULL,
        title = "Headline populations: directional shift by stratum",
        subtitle = "Each cell shows the population's per-nhood LFC distribution: % positive (fill), n_sig_up:n_sig_dn, median LFC")

out_dir_pdf <- file.path(fig_dir, "headline_population_directionality.pdf")
ggsave(out_dir_pdf, p_dir, width = 7.5, height = 4)
cat(sprintf("Wrote: %s\n", out_dir_pdf))

cat("\n=== headline panels done ===\n")
