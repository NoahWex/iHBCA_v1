#!/usr/bin/env Rscript
# Render Fig 2 supplemental: blocking-schema covariate aliasing heatmap.
#
# Plot structure
#   Rows + columns = covariates from the parity-x-risk inquiry's term universe,
#     sorted by formula role (tested first, then blocked, then context, then
#     excluded). Tile fill encodes the pairwise verdict from Stage B
#     diagnostics:
#       collinear     -> red          (V >= 0.4)
#       orthogonal    -> white        (uncorrelated; viable for blocking)
#       indeterminate -> light gray   (insufficient power to call)
#   Cell text = numeric association value (Cramer's V for categorical,
#     sqrt(eta^2) for mixed, |Spearman| for ordinal).
#   Axis-label color encodes formula role: tested / blocked /
#     blocked_sensitivity / excluded / context / input. This communicates
#     "what the design controls for" alongside the aliasing structure.
#
# Inputs (resolved via --project-root + --inquiry-name + --substrate-dir)
#   <substrate>/outputs/stageB_exploration/collinearity_matrix.csv
#   <substrate>/inquiry.yaml  (term universe + stageB.formula_role map)
#
# Output
#   <out_dir>/s_ihbca_power_matrix.pdf

suppressPackageStartupMessages({
  library(argparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--project-root", required = TRUE,
                    help = "Repository root (resolves to publication/config/)")
parser$add_argument("--inquiry-name", default = "parity_findings",
                    help = "Inquiry directory under publication/analysis/abundance/clinical/")
parser$add_argument("--substrate-dir", default = NULL,
                    help = paste("Override substrate directory for stageB",
                                 "exploration outputs; defaults to promoted location"))
parser$add_argument("--out-dir", default = NULL,
                    help = "Output directory (default: <root>/publication/figures/output/supplemental/covariate_landscape/)")
args <- parser$parse_args()

root <- normalizePath(args$project_root)

# ----------------------------------------------------------------------------
# Aesthetics framework
# ----------------------------------------------------------------------------
source(file.path(root, "publication", "config", "load_aesthetics.R"))
aes_cfg <- load_aesthetics(file.path(root, "publication", "config"))

# Verdict palette: borrow registered scale endpoints (no inline hex).
#   collinear     <- expression_diverging.high (red)
#   orthogonal    <- expression_diverging.mid  (white)
#   indeterminate <- qc_status.score_ramp_low  (light gray)
div_scale <- aes_cfg$scales$expression_diverging
qc_pal    <- get_palette("qc_status", aes_cfg)
verdict_pal <- c(
  collinear     = div_scale$high,
  orthogonal    = div_scale$mid,
  indeterminate = qc_pal[["score_ramp_low"]]
)

# Role palette for axis-label coloring: reuse six Tol slots from
# nhoodgroup_rank. Slot choices map by sort priority (tested + blocked dominant).
rank_pal <- get_palette("nhoodgroup_rank", aes_cfg)
role_pal <- c(
  tested              = rank_pal[["2"]],   # green — direct test variable
  blocked             = rank_pal[["1"]],   # dark blue — in formula
  blocked_sensitivity = rank_pal[["4"]],   # light blue — in formula w/ caveat
  excluded            = rank_pal[["6"]],   # rose — filter-restricted to 1 level
  context             = rank_pal[["5"]],   # sand — present, not in formula
  input               = rank_pal[["9"]]    # olive — defines derived terms
)

# Default role map (overridden by inquiry.yaml stageB.formula_role if present)
default_role <- list(
  study = "blocked",
  facs_status = "blocked",
  age_binary = "blocked_sensitivity",
  menopausal_status_binary = "blocked_sensitivity",
  cancer_history = "blocked",
  stratum = "tested",
  risk_class = "tested",
  parity_binary = "tested",
  ethnicity_grouped = "excluded",
  bmi_continuous = "excluded",
  bmi_category = "excluded",
  tissue_indication = "excluded",
  dissociation_minutes = "context",
  brca_genotype = "input"
)

# ----------------------------------------------------------------------------
# Substrate resolution
# ----------------------------------------------------------------------------
substrate_dir <- args$substrate_dir
if (is.null(substrate_dir)) {
  substrate_dir <- file.path(root, "publication", "analysis", "abundance",
                              "clinical", args$inquiry_name, "cohort_design")
}
substrate_dir <- normalizePath(substrate_dir, mustWork = FALSE)

collin_path  <- file.path(substrate_dir, "outputs", "stageB_exploration",
                           "collinearity_matrix.csv")
inquiry_yaml_cfg <- file.path(substrate_dir, "config", "inquiry.yaml")
inquiry_yaml <- if (file.exists(inquiry_yaml_cfg)) inquiry_yaml_cfg else file.path(substrate_dir, "inquiry.yaml")

if (!file.exists(collin_path)) {
  stop("collinearity_matrix.csv not found at: ", collin_path)
}

inquiry <- list()
if (file.exists(inquiry_yaml) && requireNamespace("yaml", quietly = TRUE)) {
  inquiry <- yaml::read_yaml(inquiry_yaml)
}
fr_cfg <- (inquiry$stageB %||% list())$formula_role %||% default_role
term_role <- function(term) fr_cfg[[term]] %||% "context"

# ----------------------------------------------------------------------------
# Load V matrix; restrict to pooled cohort
# ----------------------------------------------------------------------------
collin_raw <- read_csv(collin_path, show_col_types = FALSE) %>%
  filter(is.na(within_study) | within_study == "NA")
available_cohorts <- unique(collin_raw$cohort)
pooled_cohort_name <- (inquiry$stageB %||% list())$pooled_cohort %||% {
  priority <- c("pooled", "pooled_main_effects",
                 "pooled_interaction", "BR1_vs_AR", "AR_baseline")
  match <- intersect(priority, available_cohorts)
  if (length(match) > 0) match[1] else sort(available_cohorts)[1]
}
if (!pooled_cohort_name %in% available_cohorts) {
  stop(sprintf("pooled cohort '%s' not in collinearity_matrix; available: %s",
               pooled_cohort_name,
               paste(available_cohorts, collapse = ", ")))
}
collin <- collin_raw %>% filter(cohort == pooled_cohort_name)

# Term universe: from inquiry.yaml; restrict to terms present in the matrix.
term_universe <- inquiry$terms_of_interest %||%
  unique(c(collin$term1, collin$term2))
present_terms <- intersect(term_universe,
                            unique(c(collin$term1, collin$term2)))
for (extra in c("stratum", "risk_class", "brca_genotype")) {
  if (extra %in% c(collin$term1, collin$term2)) {
    present_terms <- union(present_terms, extra)
  }
}

# ----------------------------------------------------------------------------
# Build symmetric V matrix + verdict matrix
# ----------------------------------------------------------------------------
mat         <- matrix(NA_real_, nrow = length(present_terms),
                       ncol = length(present_terms),
                       dimnames = list(present_terms, present_terms))
verdict_mat <- matrix(NA_character_, nrow = length(present_terms),
                       ncol = length(present_terms),
                       dimnames = list(present_terms, present_terms))
for (i in seq_len(nrow(collin))) {
  t1 <- collin$term1[i]; t2 <- collin$term2[i]
  v  <- as.numeric(collin$association[i])
  vd <- collin$verdict[i]
  if (t1 %in% present_terms && t2 %in% present_terms) {
    mat[t1, t2] <- v; mat[t2, t1] <- v
    verdict_mat[t1, t2] <- vd; verdict_mat[t2, t1] <- vd
  }
}

long <- expand.grid(term1 = present_terms, term2 = present_terms,
                     stringsAsFactors = FALSE) %>%
  mutate(V       = mapply(function(t1, t2) mat[t1, t2], term1, term2),
         verdict = mapply(function(t1, t2) verdict_mat[t1, t2], term1, term2))

long$verdict_simple <- case_when(
  is.na(long$verdict) | grepl("INDETERMINATE", long$verdict) ~ "indeterminate",
  grepl("^COLLINEAR", long$verdict)  ~ "collinear",
  grepl("^ORTHOGONAL", long$verdict) ~ "orthogonal",
  TRUE                                ~ "indeterminate"
)

# ----------------------------------------------------------------------------
# Order terms by formula role
# ----------------------------------------------------------------------------
term_order_priority <- c(tested = 1, blocked = 2, blocked_sensitivity = 2,
                          input = 3, context = 4, excluded = 5)
present_terms_sorted <- present_terms[order(
  vapply(present_terms,
         function(t) term_order_priority[[term_role(t)]] %||% 6,
         numeric(1)),
  present_terms
)]
long$term1 <- factor(long$term1, levels = present_terms_sorted)
long$term2 <- factor(long$term2, levels = rev(present_terms_sorted))

axis_label_color <- vapply(present_terms_sorted,
                            function(t) role_pal[[term_role(t)]] %||% "black",
                            character(1))

# ----------------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------------
p <- ggplot(long, aes(x = term1, y = term2)) +
  geom_tile(aes(fill = verdict_simple), color = "grey80", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is.na(V) | term1 == term2, "",
                                sprintf("%.2f", V))),
            size = 1.9) +
  scale_fill_manual(values = verdict_pal, name = "Stage B verdict",
                     na.value = qc_pal[["score_ramp_low"]]) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  get_theme(aes_cfg) +
  theme(axis.text.x.top = element_text(angle = 45, hjust = 0,
                                         color = axis_label_color),
         axis.text.y     = element_text(color = rev(axis_label_color)),
         panel.grid      = element_blank(),
         legend.position = "bottom")

# ----------------------------------------------------------------------------
# Save
# ----------------------------------------------------------------------------
out_dir <- args$out_dir %||% file.path(root, "publication", "figures",
                                        "output", "supplemental",
                                        "covariate_landscape")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

panel_id   <- "s_ihbca_power_matrix"
panel_type <- "heatmap_tile"
n_terms    <- length(present_terms)

validate_panel(p, panel_type = panel_type, n_groups = n_terms,
               config = aes_cfg)
save_panel(p, panel_id = panel_id, panel_type = panel_type,
           output_dir = out_dir, config = aes_cfg)
