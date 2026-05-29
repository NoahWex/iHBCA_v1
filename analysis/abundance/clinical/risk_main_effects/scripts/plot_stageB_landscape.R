#!/usr/bin/env Rscript
# plot_stageB_landscape.R
# Render Stage B collinearity heatmap for the supplemental.
#
# Inputs (auto-resolved via lib/load_paths.R):
#   stageB_exploration/collinearity_matrix.csv
#
# Outputs to: reports/figures/
#   collinearity_heatmap_main.pdf  (term1 × term2 by within_study; main cohort only)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(scales)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_stageB_landscape.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
inquiry <- cfg$inquiry

fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

cat(sprintf("=== Plot Stage B landscape for '%s' ===\n", inquiry$inquiry))

collin <- read_csv(file.path(paths$outputs$stageB, "collinearity_matrix.csv"),
                    show_col_types = FALSE)
cat(sprintf("collinearity rows: %d\n", nrow(collin)))

# Single collinearity heatmap: main cohort only (interaction cohort is a superset
# and the redundant facet panel adds no decision value).
cohort_name <- "main"
d <- collin %>%
  filter(cohort == cohort_name, !is.na(association)) %>%
  mutate(within_study = ifelse(is.na(within_study), "(pooled)", within_study),
         term1 = factor(term1), term2 = factor(term2))

if (nrow(d) == 0) {
  stop(sprintf("No collinearity rows for cohort '%s'", cohort_name))
}

d_sym <- bind_rows(d, d %>% rename(term1 = term2, term2 = term1))
term_order <- sort(unique(c(as.character(d$term1), as.character(d$term2))))
d_sym$term1 <- factor(d_sym$term1, levels = term_order)
d_sym$term2 <- factor(d_sym$term2, levels = rev(term_order))

p <- ggplot(d_sym, aes(x = term1, y = term2, fill = association)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is.na(association), "",
                                sprintf("%.2f", association))),
            size = 1.8, colour = "grey20") +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                        midpoint = 0.5, limits = c(0, 1),
                        na.value = "grey90", name = "|assoc|") +
  facet_wrap(~ within_study, ncol = 4) +
  theme_minimal(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", colour = NA)) +
  labs(title = sprintf("Pairwise collinearity — cohort %s", cohort_name),
       subtitle = "Cramér's V (categorical) / |Spearman| (numeric) / sqrt(η²) (mixed); pooled + per-study",
       x = NULL, y = NULL)

out <- file.path(fig_dir, sprintf("collinearity_heatmap_%s.pdf", cohort_name))
ggsave(out, p, width = 16, height = 12)
cat(sprintf("  wrote: %s\n", out))

cat("=== plotting done ===\n")
