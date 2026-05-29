#!/usr/bin/env Rscript
# plot_nhoodgroup_renaming.R
# Visualize the NhoodGroup renaming lookup so the size-rank assignment is legible.
# One panel per parent_L2 (filtered to those with viable groups in >= 2 contrasts
# OR a headline L2). Within a panel, viable NhoodGroups are points:
#   x = group_med_lfc, y = stratum (ordered AR / HRS / BR1 / BR2 / BR12)
#   point size = n_nhoods_in_group, color = stratum
#   text label = name_short (e.g. "BMYO-basal_1") with original NhoodGroup id below
# Vertical dashed line at lfc = 0.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel); library(tidyr)
})

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_nhoodgroup_renaming.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inq <- paths$inquiry_root

lookup_path <- file.path(inq, "outputs", "nhoodgroup_renaming", "lookup.csv")
df_all <- read_csv(lookup_path, show_col_types = FALSE)

df <- df_all %>% filter(viable)

# Pick parent_L2s to plot: any L2 with >= 2 viable groups summed across contrasts,
# OR the four headline L2s (BMYO-basal, Fibro-SFRP4, LASP-major, LHS-major).
keep_l2 <- df %>%
  count(parent_L2_joint, parent_label, name = "n_viable_total") %>%
  filter(n_viable_total >= 4) %>%
  pull(parent_L2_joint)
keep_l2 <- union(keep_l2, c("epi::BMYO-basal", "str::Fibro-SFRP4",
                              "epi::LASP-major", "epi::LHS-major"))

df <- df %>%
  filter(parent_L2_joint %in% keep_l2) %>%
  mutate(stratum = factor(stratum,
                          levels = c("AR", "HRS", "BR1", "BR2", "BR12",
                                      "BR1m", "BR2m", "HRSm")))

stratum_pal <- c(
  AR     = "#1f77b4",
  HRS    = "#2ca02c",
  BR1    = "#d62728",
  BR2    = "#9467bd",
  BR12   = "#8c564b",
  BR1m   = "#7f7f7f",
  BR2m   = "#7f7f7f",
  HRSm   = "#7f7f7f"
)

# Order parent_L2 facets: epi first, then str, then imm; within each compartment by
# total viable count desc so most-content facets come first.
l2_order <- df %>%
  count(parent_compartment, parent_L2_joint, name = "n_viable_total") %>%
  arrange(parent_compartment, desc(n_viable_total)) %>%
  pull(parent_L2_joint)
df$parent_L2_joint <- factor(df$parent_L2_joint, levels = l2_order)

p <- ggplot(df, aes(x = group_med_lfc, y = stratum,
                     color = stratum, size = n_nhoods_in_group)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  geom_point(alpha = 0.85) +
  geom_text_repel(aes(label = name_short),
                   size = 2.4, color = "black", max.overlaps = Inf,
                   box.padding = 0.25, point.padding = 0.2,
                   segment.size = 0.2, segment.color = "grey50",
                   show.legend = FALSE,
                   min.segment.length = 0) +
  scale_color_manual(values = stratum_pal, drop = FALSE,
                      breaks = c("AR", "HRS", "BR1", "BR2", "BR12")) +
  scale_size_continuous(range = c(1.5, 7), trans = "sqrt",
                         breaks = c(20, 50, 100, 250, 500),
                         name = "n_nhoods") +
  facet_wrap(~ parent_L2_joint, ncol = 3, scales = "free_y") +
  labs(x = "group_med_lfc",
        y = NULL,
        color = "stratum") +
  theme_bw(base_size = 9) +
  theme(strip.text = element_text(size = 8.5, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal")

n_facets <- length(unique(df$parent_L2_joint))
n_rows <- ceiling(n_facets / 3)

out_dir <- file.path(inq, "..", "..", "share", "kai_v1", "04_panels")
out_dir <- normalizePath(out_dir, mustWork = FALSE)
# fall back to inquiry-root share if the relative join doesn't resolve.
fallback_share <- file.path(inq, "share", "kai_v1", "04_panels")
if (dir.exists(fallback_share)) out_dir <- fallback_share
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pdf_path <- file.path(out_dir, "nhoodgroup_naming_overview.pdf")
ggsave(pdf_path, p,
        width = 10, height = max(4.5, 2.4 * n_rows),
        device = cairo_pdf)
cat(sprintf("Wrote: %s\n", pdf_path))

# Companion CSV: just the rows that were plotted, ordered as displayed.
display_csv <- df %>%
  arrange(parent_L2_joint, stratum, size_rank_within_L2) %>%
  select(parent_L2_joint, parent_compartment, parent_label,
         contrast, stratum, NhoodGroup, NhoodGroup_renamed,
         n_nhoods_in_group, n_sig, group_med_lfc, group_pct_up,
         size_rank_within_L2, name_short, name_full)
csv_out <- file.path(out_dir, "nhoodgroup_naming_overview.csv")
write_csv(display_csv, csv_out)
cat(sprintf("Wrote: %s (%d rows)\n", csv_out, nrow(display_csv)))

cat(sprintf("\nFacet count: %d  L2s: %s\n",
            n_facets,
            paste(levels(df$parent_L2_joint), collapse = ", ")))
cat("=== plot_nhoodgroup_renaming done ===\n")
