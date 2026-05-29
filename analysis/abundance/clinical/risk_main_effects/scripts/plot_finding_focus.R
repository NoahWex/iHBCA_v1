#!/usr/bin/env Rscript
# plot_finding_focus.R
# One focused beeswarm per headline finding. For each (L2, finding) tuple:
#   - Filter Stage D nhoods to that L2 across all 4 contrasts
#   - Beeswarm: x = logFC, y = contrast (4 rows)
#   - Color by FDR significance (sig=red/blue, ns=grey)
#   - Title: L2 + interpretation; subtitle: β_AR + Δβ for each modifier
#
# Output: reports/figures/finding_<L2>.pdf  (one PDF per finding)

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(yaml)
})
have_beeswarm <- requireNamespace("ggbeeswarm", quietly = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "plot_finding_focus.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1",
                "parity_x_HR_sporadic", "parity_x_HR_BRCA2")

# Headline findings: L2 (compartment::label match) + caption
FINDINGS <- tribble(
  ~slug,                   ~compartment, ~label,        ~caption,
  "lasp_major_brca1_amp",  "epi",        "LASP-major",  "LASP-major × BRCA1: AMPLIFY parity-driven enrichment (Bach 2021)",
  "lasp_major_brca2_rev",  "epi",        "LASP-major",  "LASP-major × BRCA2: REVERSE — novel finding",
  "bmyo_basal_brca1_emt",  "epi",        "BMYO-basal",  "BMYO-basal × BRCA1: introduces EMT-hybrid signal (Wang 2023)",
  "cd8_trm_brca1_amp",     "imm",        "CD8_Trm",     "CD8_Trm × BRCA1: AMPLIFY exhaustion (Reed 2024)",
  "fibro_sfrp4_brca1_rev", "str",        "Fibro-SFRP4", "Fibro-SFRP4 × BRCA1: REVERSE"
)

# Load Stage D contrasts (filter to needed L2s for memory).
all_da <- list()
for (cname in CONTRASTS) {
  d_csv <- file.path(paths$outputs$stageD, cname, "da_results.csv")
  if (!file.exists(d_csv)) next
  d <- read_csv(d_csv, show_col_types = FALSE) %>%
    mutate(sig = !is.na(SpatialFDR) & SpatialFDR < 0.05,
           contrast = cname)
  all_da[[cname]] <- d
}
da <- bind_rows(all_da)
cat(sprintf("Stage D loaded: %d rows total\n", nrow(da)))

# Stage E for summary stats
e_long <- read_csv(paths$outputs$stageE, show_col_types = FALSE) %>%
  tidyr::pivot_longer(-c(L2_joint, compartment, label),
                       names_to = "key", values_to = "value") %>%
  tidyr::separate(key, into = c("contrast", "metric"),
                   sep = "__", extra = "merge") %>%
  tidyr::pivot_wider(names_from = metric, values_from = value)

CONTRAST_LABELS <- c(
  parity_in_AR = "AR baseline",
  parity_x_HR_BRCA1 = "BRCA1 modifier",
  parity_x_HR_sporadic = "sporadic modifier",
  parity_x_HR_BRCA2 = "BRCA2 modifier"
)

render_finding <- function(slug, comp, lab, caption) {
  d <- da %>% filter(compartment == comp, label == lab)
  if (nrow(d) == 0) {
    cat(sprintf("  [%s] no nhoods, skipping\n", slug)); return(invisible(NULL))
  }
  d$contrast <- factor(d$contrast, levels = CONTRASTS,
                        labels = CONTRAST_LABELS[CONTRASTS])
  d$dir <- with(d, ifelse(sig & logFC > 0, "up",
                    ifelse(sig & logFC < 0, "dn", "ns")))

  # Stats per contrast
  stats <- e_long %>%
    filter(compartment == comp, label == lab,
           contrast %in% CONTRASTS) %>%
    transmute(contrast_label = CONTRAST_LABELS[contrast],
              med_lfc, n_sig, n)
  stats$contrast_label <- factor(stats$contrast_label,
                                   levels = CONTRAST_LABELS[CONTRASTS])

  # Subtitle text
  ar_row <- stats %>% filter(contrast_label == CONTRAST_LABELS["parity_in_AR"])
  brca1_row <- stats %>% filter(contrast_label == CONTRAST_LABELS["parity_x_HR_BRCA1"])
  brca2_row <- stats %>% filter(contrast_label == CONTRAST_LABELS["parity_x_HR_BRCA2"])
  spor_row <- stats %>% filter(contrast_label == CONTRAST_LABELS["parity_x_HR_sporadic"])
  fmt_row <- function(row, prefix) {
    if (nrow(row) == 0) return(sprintf("%s: —", prefix))
    sprintf("%s: lfc=%+.2f, n_sig=%d/%d", prefix,
             row$med_lfc[1] %||% NA_real_,
             as.integer(row$n_sig[1] %||% 0),
             as.integer(row$n[1] %||% 0))
  }
  subt <- paste(c(
    fmt_row(ar_row, "β_AR"),
    fmt_row(brca1_row, "Δβ_BRCA1"),
    fmt_row(spor_row, "Δβ_sporadic"),
    fmt_row(brca2_row, "Δβ_BRCA2")
  ), collapse = "  |  ")

  # Sample down ns to keep plot legible
  d_sig <- d %>% filter(sig)
  d_ns <- d %>% filter(!sig)
  if (nrow(d_ns) > 8000) d_ns <- d_ns[sample.int(nrow(d_ns), 8000), ]
  d_plot <- bind_rows(d_sig, d_ns)

  swarm <- if (have_beeswarm) {
    ggbeeswarm::geom_quasirandom(aes(x = logFC, y = contrast, colour = dir),
                                  groupOnX = FALSE, size = 0.4, alpha = 0.7,
                                  bandwidth = 0.5)
  } else {
    geom_jitter(aes(x = logFC, y = contrast, colour = dir),
                height = 0.25, width = 0, size = 0.4, alpha = 0.7)
  }

  p <- ggplot(d_plot) +
    swarm +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40",
               linewidth = 0.4) +
    # Median LFC (across sig + ns) per contrast as a black tick
    stat_summary(aes(x = logFC, y = contrast),
                  fun = median, geom = "point",
                  shape = "|", size = 6, colour = "black", inherit.aes = FALSE,
                  data = d) +
    scale_colour_manual(values = c(ns = "grey75", up = "#B2182B", dn = "#2166AC"),
                         guide = guide_legend(override.aes = list(size = 2.5,
                                                                    alpha = 1))) +
    coord_cartesian(xlim = c(-4, 4)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(),
          axis.text.y = element_text(face = "bold")) +
    labs(title = caption,
         subtitle = subt,
         x = "logFC (per nhood)",
         y = NULL,
         colour = "FDR<0.05")

  out <- file.path(fig_dir, sprintf("finding_%s.pdf", slug))
  ggsave(out, p, width = 10, height = 4.5)
  cat(sprintf("  wrote: %s (%d nhoods, %d sig)\n",
              out, nrow(d), nrow(d_sig)))
}

for (i in seq_len(nrow(FINDINGS))) {
  r <- FINDINGS[i, ]
  render_finding(r$slug, r$compartment, r$label, r$caption)
}

cat("\n=== finding focus plots done ===\n")
