#!/usr/bin/env Rscript
# plot_nhoodgroup_beeswarm.R
# Per-NhoodGroup beeswarm of nhood-level logFC, one panel per contrast.
#
# Inputs (per F.1 contrast directory):
#   nhood_groups.csv         per-nhood: Nhood, NhoodGroup_renamed, logFC, SpatialFDR,
#                            parent_L2_joint, parent_compartment, n_nhoods_in_group,
#                            group_med_lfc, group_pct_up
#   nhood_groups_summary.csv per-group aggregates (used for sort + viability)
#
# Layout: y = NhoodGroup_renamed (sorted by group_med_lfc, descending — strongest
# positive groups at top, strongest negative at bottom). Faceted by parent
# compartment. x = per-nhood logFC. Color by significance.
#
# Two panels per contrast:
#   nhoodgroup_beeswarm_<contrast>_all.pdf            unfiltered groups
#   nhoodgroup_beeswarm_<contrast>_limma_viable.pdf   groups with n_nhoods >= 10
#                                                      (proxy for F.3 viability)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(yaml)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_nhoodgroup_beeswarm.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

f1_dir <- paths$outputs$stageF1
contrast_dirs <- list.dirs(f1_dir, recursive = FALSE)
if (length(contrast_dirs) == 0) {
  cat("No F.1 outputs at", f1_dir, "\n"); quit(status = 0)
}

# Viability proxy: limma-voom on (group × donor) pseudobulk needs ~10 cells/donor
# per parity class. Without per-cell donor counts here, we use n_nhoods >= 10 as
# the proxy — groups smaller than that rarely sustain viable per-donor pseudobulk.
VIABLE_MIN_NHOODS <- 10
MIN_GROUPS_TO_PLOT <- 10  # skip contrasts with too few groups (stratified)

build_panel <- function(d, fig_path, title, subtitle,
                         facet_mode = c("by_compartment", "by_L2_compartment")) {
  facet_mode <- match.arg(facet_mode)
  if (nrow(d) == 0) {
    cat("  empty data, skipping:", fig_path, "\n"); return(invisible(NULL))
  }
  fdr_t <- 0.05
  d <- d %>%
    mutate(sig = !is.na(SpatialFDR) & SpatialFDR < fdr_t,
           dir = ifelse(sig & logFC > 0, "up",
                  ifelse(sig & logFC < 0, "dn", "ns")))

  # Sort axes:
  #   by_compartment       — sort groups by group_med_lfc within compartment
  #   by_L2_compartment    — group within parent L2, sort groups by med within L2,
  #                          then sort L2 blocks by L2 median group_med_lfc
  if (facet_mode == "by_compartment") {
    group_order <- d %>%
      distinct(NhoodGroup_renamed, parent_compartment, group_med_lfc) %>%
      arrange(parent_compartment, desc(group_med_lfc))
  } else {
    l2_order <- d %>%
      distinct(parent_compartment, parent_L2_joint, NhoodGroup_renamed,
               group_med_lfc) %>%
      group_by(parent_compartment, parent_L2_joint) %>%
      summarise(l2_med = median(group_med_lfc, na.rm = TRUE), .groups = "drop") %>%
      arrange(parent_compartment, desc(l2_med)) %>%
      mutate(l2_rank = row_number())
    group_order <- d %>%
      distinct(parent_compartment, parent_L2_joint, NhoodGroup_renamed,
               group_med_lfc) %>%
      left_join(l2_order, by = c("parent_compartment", "parent_L2_joint")) %>%
      arrange(parent_compartment, l2_rank, desc(group_med_lfc))
    d <- d %>%
      mutate(parent_L2_joint = factor(parent_L2_joint,
                                       levels = l2_order$parent_L2_joint))
  }
  d$NhoodGroup_renamed <- factor(d$NhoodGroup_renamed,
                                  levels = rev(group_order$NhoodGroup_renamed))

  swarm_layer <- if (have_beeswarm) {
    ggbeeswarm::geom_quasirandom(
      aes(x = logFC, y = NhoodGroup_renamed, colour = dir),
      groupOnX = FALSE, size = 0.30, alpha = 0.55, bandwidth = 0.5
    )
  } else {
    geom_jitter(aes(x = logFC, y = NhoodGroup_renamed, colour = dir),
                width = 0.0, height = 0.18, size = 0.30, alpha = 0.55)
  }

  p <- ggplot(d) +
    swarm_layer +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.3) +
    scale_colour_manual(values = c(ns = "grey75", up = "#B2182B", dn = "#2166AC"),
                        guide = guide_legend(override.aes = list(size = 2.5,
                                                                  alpha = 1)))

  if (facet_mode == "by_compartment") {
    p <- p + facet_grid(parent_compartment ~ ., scales = "free_y", space = "free_y")
  } else {
    p <- p + facet_grid(parent_compartment + parent_L2_joint ~ .,
                         scales = "free_y", space = "free_y",
                         labeller = labeller(parent_L2_joint = label_wrap_gen(20)))
  }

  p <- p +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_minimal(base_size = 8) +
    theme(axis.text.y = element_text(size = 5.5),
          axis.text.x = element_text(size = 8),
          panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          strip.text.y = element_text(face = "bold", size = 7, angle = 0),
          legend.position = "right") +
    labs(title = title, subtitle = subtitle,
         x = "logFC (per nhood)", y = NULL, colour = "FDR<0.05")

  n_groups <- nlevels(d$NhoodGroup_renamed)
  fig_h <- max(8, min(32, 0.10 * n_groups + 4))
  fig_w <- if (facet_mode == "by_L2_compartment") 12 else 10
  ggsave(fig_path, p, width = fig_w, height = fig_h, limitsize = FALSE)
  cat(sprintf("  wrote: %s (%d groups, %d nhoods, mode=%s)\n",
              fig_path, n_groups, nrow(d), facet_mode))
}

for (cdir in contrast_dirs) {
  cname <- basename(cdir)
  ng_csv <- file.path(cdir, "nhood_groups.csv")
  if (!file.exists(ng_csv) || file.size(ng_csv) < 200) {
    cat(sprintf("  skip %s: no nhood_groups.csv\n", cname)); next
  }
  ng <- read_csv(ng_csv, show_col_types = FALSE) %>%
    filter(!is.na(NhoodGroup_renamed))
  n_groups <- length(unique(ng$NhoodGroup_renamed))
  cat(sprintf("\n--- %s: %d groups, %d nhoods ---\n",
              cname, n_groups, nrow(ng)))
  if (n_groups < MIN_GROUPS_TO_PLOT) {
    cat("  too few groups, skipping\n"); next
  }

  # Unfiltered — by compartment
  build_panel(
    ng,
    file.path(fig_dir, sprintf("nhoodgroup_beeswarm_%s_all.pdf", cname)),
    title = sprintf("%s — NhoodGroup beeswarm (all groups)", cname),
    subtitle = sprintf("Y sorted by group_med_lfc (desc) within compartment. n=%d groups.",
                        n_groups),
    facet_mode = "by_compartment"
  )

  # Unfiltered — grouped by parent L2, faceted by compartment
  build_panel(
    ng,
    file.path(fig_dir, sprintf("nhoodgroup_beeswarm_%s_all_by_L2.pdf", cname)),
    title = sprintf("%s — NhoodGroup beeswarm grouped by parent L2 (all groups)", cname),
    subtitle = sprintf("Rows nested as compartment :: parent L2; groups within L2 sorted by group_med_lfc desc. n=%d groups.",
                        n_groups),
    facet_mode = "by_L2_compartment"
  )

  # Limma-viable (n_nhoods >= 10) — by compartment
  ng_v <- ng %>% filter(n_nhoods_in_group >= VIABLE_MIN_NHOODS)
  n_v <- length(unique(ng_v$NhoodGroup_renamed))
  cat(sprintf("  viable subset: %d groups (n_nhoods>=%d)\n", n_v, VIABLE_MIN_NHOODS))
  if (n_v >= MIN_GROUPS_TO_PLOT) {
    build_panel(
      ng_v,
      file.path(fig_dir, sprintf("nhoodgroup_beeswarm_%s_limma_viable.pdf", cname)),
      title = sprintf("%s — NhoodGroup beeswarm (limma-viable)", cname),
      subtitle = sprintf("Subset to groups with n_nhoods >= %d (F.3 pseudobulk viability proxy). n=%d groups.",
                          VIABLE_MIN_NHOODS, n_v),
      facet_mode = "by_compartment"
    )
    # Limma-viable — grouped by parent L2
    build_panel(
      ng_v,
      file.path(fig_dir, sprintf("nhoodgroup_beeswarm_%s_limma_viable_by_L2.pdf", cname)),
      title = sprintf("%s — NhoodGroup beeswarm grouped by parent L2 (limma-viable)", cname),
      subtitle = sprintf("Rows nested as compartment :: parent L2. Subset n_nhoods >= %d. n=%d groups.",
                          VIABLE_MIN_NHOODS, n_v),
      facet_mode = "by_L2_compartment"
    )
  }
}

cat("\n=== nhoodgroup beeswarm done ===\n")
