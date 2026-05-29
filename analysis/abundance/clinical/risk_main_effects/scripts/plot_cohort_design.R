#!/usr/bin/env Rscript
# plot_cohort_design.R
# Cohort design overview for parity_x_risk inquiry:
#   - Alluvial: study → stratum → parity_binary (biological)
#   - Alluvial: study → facs_status → sample_type (technical)
#   - Heatmap table: per-variable × per-study donor counts (bio + technical)
#
# Reads from outputs/stageA_cohorts/cohort_pooled_design.csv (most informative
# single CSV; carries the derived `stratum` column).
# Outputs to: reports/figures/cohort_design_overview.pdf

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2)
  library(stringr); library(patchwork)
})
have_alluvial <- requireNamespace("ggalluvial", quietly = TRUE)
if (have_alluvial) {
  # Need to attach the package so `stat = "stratum"` resolves via search path
  suppressPackageStartupMessages(library(ggalluvial))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_cohort_design.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths

fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

pooled_csv <- file.path(paths$outputs$stageA, "cohort_pooled_design.csv")
if (!file.exists(pooled_csv)) stop("Missing pooled cohort design: ", pooled_csv)

d <- read_csv(pooled_csv, show_col_types = FALSE)
cat(sprintf("Loaded %d donors from pooled cohort\n", nrow(d)))

# Normalize / coerce key factors ---------------------------------------------
d$stratum <- factor(d$stratum, levels = c("AR", "HR_sporadic", "HR_BRCA1", "HR_BRCA2"))
d$parity_binary <- factor(d$parity_binary, levels = c("nulliparous", "parous"))

# Bin dissociation_minutes for table readability
if ("dissociation_minutes" %in% colnames(d)) {
  d$dissoc_bin <- cut(d$dissociation_minutes,
                       breaks = c(-Inf, 10, 30, 120, 240, Inf),
                       labels = c("≤10", "11-30", "31-120", "121-240", ">240"),
                       right = TRUE)
  d$dissoc_bin <- as.character(d$dissoc_bin)
  d$dissoc_bin[is.na(d$dissoc_bin)] <- "NA"
}

study_levels <- d %>% count(study, sort = TRUE) %>% pull(study)
d$study <- factor(d$study, levels = study_levels)

stratum_pal <- c(AR = "#4D7C8A", HR_sporadic = "#E07B39",
                  HR_BRCA1 = "#9B2226", HR_BRCA2 = "#3A86FF")

# 1) Biological alluvial: study → stratum → parity ----------------------------
build_alluvial <- function(df, axes, fill_col, title) {
  flow <- df %>% count(across(all_of(axes)))
  flow$.fill <- flow[[fill_col]]
  if (have_alluvial) {
    aes_list <- lapply(seq_along(axes), function(i) sym(axes[i]))
    names(aes_list) <- paste0("axis", seq_along(axes))
    p <- ggplot(flow, aes(y = n, !!!aes_list)) +
      ggalluvial::geom_alluvium(aes(fill = .fill), alpha = 0.55,
                                  curve_type = "sigmoid",
                                  knot.pos = 0.5, width = 1/5) +
      ggalluvial::geom_stratum(width = 1/4, fill = "white", color = "grey30",
                                 linewidth = 0.4) +
      geom_text(stat = "stratum",
                aes(label = after_stat(stratum)),
                size = 2.6) +
      scale_x_discrete(limits = axes, expand = c(0.05, 0.05)) +
      theme_void(base_size = 9) +
      theme(legend.position = "right",
            plot.title = element_text(face = "bold", size = 10),
            axis.text.x = element_text(size = 8)) +
      labs(title = title, fill = fill_col,
           subtitle = sprintf("n=%d donors", sum(flow$n)))
  } else {
    # Fallback: faceted bar chart across the 3 axes
    long <- flow %>%
      pivot_longer(all_of(axes), names_to = "axis", values_to = "level")
    p <- ggplot(long, aes(x = level, y = n, fill = .fill)) +
      geom_col() + facet_wrap(~ axis, scales = "free_x") +
      theme_minimal(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = title, fill = fill_col)
  }
  p
}

p_combined <- build_alluvial(
  d %>% mutate(facs_status = ifelse(is.na(facs_status), "NA", facs_status),
                sample_type = ifelse(is.na(sample_type), "NA", sample_type)),
  axes = c("study", "stratum", "parity_binary", "facs_status", "sample_type"),
  fill_col = "stratum",
  title = "Cohort design — biological × technical flow") +
  scale_fill_manual(values = stratum_pal, drop = FALSE)

# 2) Heatmap table: var × level × study ---------------------------------------
bio_vars <- c("parity_binary", "stratum", "age_binary", "menopausal_status_binary",
              "cancer_history")
tech_vars <- c("facs_status", "sample_type", "dissoc_bin")
all_vars <- c(bio_vars, tech_vars)
present <- all_vars[all_vars %in% colnames(d)]

table_data <- lapply(present, function(v) {
  d %>%
    mutate(level = ifelse(is.na(.data[[v]]), "NA", as.character(.data[[v]]))) %>%
    count(study, level) %>%
    mutate(variable = v, axis_class = ifelse(v %in% bio_vars, "biological", "technical"))
}) %>% bind_rows()

table_data$variable <- factor(table_data$variable, levels = present)
table_data$axis_class <- factor(table_data$axis_class, levels = c("biological","technical"))

p_table <- ggplot(table_data, aes(x = study, y = level, fill = n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n,
                 color = ifelse(n > max(n) * 0.5, "white", "grey20")),
             size = 2.8, show.legend = FALSE) +
  facet_grid(variable ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_identity() +
  scale_fill_gradient(low = "#F2F2F2", high = "#3A86FF",
                       name = "n donors") +
  scale_x_discrete(position = "top") +
  theme_minimal(base_size = 9) +
  theme(strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
        strip.placement = "outside",
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 10)) +
  labs(title = "Donor counts per (variable × study)",
       subtitle = sprintf("Pooled cohort, n=%d donors. Variables grouped: biological (top) / technical (bottom)",
                           nrow(d)),
       x = NULL, y = NULL)

# 3) Compose ------------------------------------------------------------------
combined <- p_combined / p_table +
  plot_layout(heights = c(1.2, 2.2)) +
  plot_annotation(
    title = "parity × risk DA inquiry — cohort design overview",
    subtitle = "Pooled cohort (n=208). Modifier interaction tests draw from this design.",
    theme = theme(plot.title = element_text(face = "bold", size = 12))
  )

out <- file.path(fig_dir, "cohort_design_overview.pdf")
ggsave(out, combined, width = 14, height = 14, limitsize = FALSE)
cat(sprintf("Wrote: %s\n", out))

# Also write the same composition as PNG for non-vector consumers
out_png <- sub("\\.pdf$", ".png", out)
ggsave(out_png, combined, width = 14, height = 14, dpi = 200, limitsize = FALSE)
cat(sprintf("Wrote: %s\n", out_png))
