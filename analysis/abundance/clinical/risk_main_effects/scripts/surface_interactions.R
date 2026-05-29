#!/usr/bin/env Rscript
# surface_interactions.R
# Operationalize the deep-research-brief predictions against pipeline outputs.
#
# CORRECTED INTERPRETATION OF INTERACTION TERMS:
#   - parity_in_AR's median_logFC = β_AR (parity main effect within AR cohort)
#   - parity_x_HR_*'s median_logFC = Δβ_modifier (the INTERACTION; relative to
#     AR, NOT a main effect on its own).
#   - AMPLIFY: sign(Δβ) == sign(β_AR) AND |Δβ| meaningful (modifier adds in
#     the same direction as AR's parity effect → reinforces it).
#   - ATTENUATE: sign(Δβ) == -sign(β_AR) AND |Δβ| < |β_AR| (modifier opposes
#     but doesn't overwhelm AR's effect; net AR-direction with smaller magnitude).
#   - REVERSE: sign(Δβ) == -sign(β_AR) AND |Δβ| > |β_AR| (modifier flips the
#     sign of the net effect in the modified stratum).
#
# Two independent evaluations per prediction:
#   1. DIRECTION (Stage E): does the per-L2 interaction match prediction?
#   2. BIOLOGY (F.3 markers): are the predicted signature genes enriched in
#      the modifier-contrast group markers within the relevant parent L2?
#
# Outputs:
#   outputs/predictions_evidence.csv            per-prediction structured table
#   outputs/predictions_evidence.md             narrative report
#   reports/figures/predictions_signature_heatmap.pdf

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(ggplot2); library(yaml)
  library(purrr); library(stringr)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "surface_interactions.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
fig_dir <- file.path(paths$outputs$reports, "figures")
ensure_dir(fig_dir)
out_dir <- paths$inquiry_root

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1",
                "parity_x_HR_sporadic", "parity_x_HR_BRCA2")
AR <- "parity_in_AR"
MODIFIERS <- setdiff(CONTRASTS, AR)

# Thresholds for direction classification at L2 level.
MIN_NSIG_AR <- 30           # AR baseline must reach this n_sig to count
MIN_NSIG_MOD <- 10          # modifier interaction must reach this n_sig
MIN_DELTA_MEANINGFUL <- 0.15
MIN_AR_BASELINE <- 0.15

# ------------------------------------------------------------------------
# Signatures from the deep research brief — gene symbols.
# ------------------------------------------------------------------------
SIGNATURES <- list(
  milk_biosynthesis = c("CSN2", "CSN3", "LALBA", "WAP", "ELF5",
                         "PIP", "MUCL1", "CSN1S1"),
  lasp_subcluster_priming = c("ALDH1A3", "SLPI", "KRT15", "ELF5", "KIT", "CCL28"),
  proliferation = c("MKI67", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
                     "TOP2A", "CDK1", "PCNA", "BIRC5"),
  basal_classical = c("KRT5", "KRT14", "TP63", "ACTA2", "MYLK"),
  emt_hybrid = c("SNAI2", "VIM", "ZEB1", "ZEB2", "FN1"),
  fibroblast_active = c("TGFB1", "COL1A1", "COL3A1", "LOX", "CXCL12",
                         "PDGFRA", "MMP3", "MMP14"),
  precaf_brca1 = c("MMP3", "MMP14", "MMP1", "POSTN", "ACTA2"),
  cd8_exhaustion = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX",
                      "ENTPD1", "CTLA4", "EOMES"),
  treg_signature = c("FOXP3", "IL2RA", "CTLA4", "IL10", "IKZF2", "TIGIT"),
  m2_lam = c("TREM2", "LIPA", "LPL", "FABP4", "FABP5", "APOE", "APOC1",
              "CD163", "MRC1", "SPP1", "C1QC", "C1QA", "C1QB"),
  er_pr_signaling = c("ESR1", "PGR", "GREB1", "TFF1", "TFF3", "AGR2", "FOXA1"),
  mtorc1_proteostasis = c("RPS6", "EIF4E", "EIF4G1", "EIF4EBP1", "MTOR",
                           "RPS3", "RPL5", "RPL10")
)

# ------------------------------------------------------------------------
# Predictions: cell_type pattern × signature × stratum × predicted dir/class.
# Brief-derived. AR baseline rows establish β_AR direction; modifier rows
# evaluate AMPLIFY/ATTENUATE/REVERSE relative to AR.
# ------------------------------------------------------------------------
PREDICTIONS <- tribble(
  ~cell_type_pattern,    ~signature,                ~stratum,                ~predicted_class,
  # AR baseline
  "LASP",                "milk_biosynthesis",        "parity_in_AR",          "ENRICHED_PARITY",
  "LASP",                "lasp_subcluster_priming",  "parity_in_AR",          "GRADED_SHIFT",
  "BMYO|Basal",          "basal_classical",           "parity_in_AR",          "DEPLETED_PARITY",
  "Fibro",               "fibroblast_active",         "parity_in_AR",          "ENRICHED_PARITY",
  "Macro",               "m2_lam",                    "parity_in_AR",          "ENRICHED_PARITY",
  "Treg|CD4|CD8",        "treg_signature",            "parity_in_AR",          "ENRICHED_PARITY",
  # BRCA1 modifier
  "LASP",                "milk_biosynthesis",         "parity_x_HR_BRCA1",     "AMPLIFY",
  "LASP",                "proliferation",             "parity_x_HR_BRCA1",     "AMPLIFY",
  "CD8|T_",              "cd8_exhaustion",            "parity_x_HR_BRCA1",     "AMPLIFY",
  "BMYO|Basal",          "basal_classical",           "parity_x_HR_BRCA1",     "ATTENUATE",
  "BMYO|Basal",          "emt_hybrid",                "parity_x_HR_BRCA1",     "AMPLIFY",
  "Fibro",               "precaf_brca1",              "parity_x_HR_BRCA1",     "AMPLIFY",
  "Macro",               "m2_lam",                    "parity_x_HR_BRCA1",     "AMPLIFY",
  # BRCA2 modifier
  "LASP",                "milk_biosynthesis",         "parity_x_HR_BRCA2",     "AR_LIKE_SMALLER",
  "LASP",                "mtorc1_proteostasis",       "parity_x_HR_BRCA2",     "ENRICHED_BRCA2",
  # HR_sporadic — cluster check
  "LASP",                "milk_biosynthesis",         "parity_x_HR_sporadic",  "MIXED_DEPENDS_OCCULT",
  "BMYO|Basal",          "basal_classical",           "parity_x_HR_sporadic",  "MIXED_DEPENDS_OCCULT",
  "Macro",               "m2_lam",                    "parity_x_HR_sporadic",  "TREATMENT_RESIDUAL"
)

# ------------------------------------------------------------------------
# Stage E: load + reshape wide → long.
# Wide schema: L2_joint, compartment, label, <contrast>__<metric>
# ------------------------------------------------------------------------
e_path <- paths$outputs$stageE
if (!file.exists(e_path)) stop("Stage E missing: ", e_path)
e_wide <- read_csv(e_path, show_col_types = FALSE)
metric_cols <- setdiff(colnames(e_wide), c("L2_joint", "compartment", "label"))

e_long <- e_wide %>%
  pivot_longer(all_of(metric_cols), names_to = "key", values_to = "value") %>%
  separate(key, into = c("contrast", "metric"), sep = "__", extra = "merge") %>%
  pivot_wider(names_from = metric, values_from = value)

cat(sprintf("Stage E loaded: %d L2s × %d contrasts (long form: %d rows)\n",
            nrow(e_wide), length(unique(e_long$contrast)), nrow(e_long)))

# ------------------------------------------------------------------------
# F.3 markers — for biology check.
# ------------------------------------------------------------------------
load_markers <- function() {
  rows <- list()
  for (cname in CONTRASTS) {
    cdir <- file.path(paths$outputs$stageF3, cname)
    if (!dir.exists(cdir)) next
    files <- list.files(cdir, pattern = "_vs_parent\\.csv$", full.names = TRUE)
    for (f in files) {
      tt <- tryCatch(read_csv(f, show_col_types = FALSE),
                     error = function(e) NULL)
      if (is.null(tt) || nrow(tt) == 0) next
      group_renamed <- sub("_vs_parent\\.csv$", "", basename(f))
      tt$contrast <- cname
      tt$NhoodGroup_renamed <- group_renamed
      rows[[length(rows) + 1]] <- tt
    }
  }
  bind_rows(rows)
}
markers <- load_markers()

# Parent L2 mapping from F.1 with NhoodGroup_renamed normalization to match
# F.3's filesystem-safe form.
parent_map <- list()
for (cname in CONTRASTS) {
  gs_csv <- file.path(paths$outputs$stageF1, cname, "nhood_groups_summary.csv")
  if (!file.exists(gs_csv)) next
  gs <- read_csv(gs_csv, show_col_types = FALSE) %>%
    transmute(NhoodGroup_renamed = as.character(NhoodGroup_renamed),
              parent_L2_joint = as.character(parent_L2_joint),
              parent_compartment = as.character(parent_compartment),
              parent_label = as.character(parent_label),
              n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group))) %>%
    mutate(NhoodGroup_renamed = gsub("-", "_", gsub("::", "__",
                                                      NhoodGroup_renamed)))
  parent_map[[cname]] <- gs
}
parent_df <- bind_rows(parent_map, .id = "contrast")
markers <- markers %>%
  left_join(parent_df, by = c("contrast", "NhoodGroup_renamed"))
sym_col <- intersect(c("symbol", "Symbol", "gene_symbol"), colnames(markers))[1]
if (is.na(sym_col)) sym_col <- "gene_id"
fdr_col <- intersect(c("FDR", "adj.P.Val", "padj"), colnames(markers))[1]
if (is.na(fdr_col)) fdr_col <- "P.Value"

cat(sprintf("F.3 markers loaded: %d rows; sym_col=%s fdr_col=%s\n",
            nrow(markers), sym_col, fdr_col))

# ------------------------------------------------------------------------
# Direction classification helpers
# ------------------------------------------------------------------------
classify_modifier_direction <- function(beta_ar, n_sig_ar,
                                          delta_beta, n_sig_mod) {
  # NO_BASELINE only if AR truly has no signal (few sig nhoods OR β_AR sign is
  # unstable). With many sig nhoods even a small magnitude is informative.
  if (is.na(beta_ar) || is.na(n_sig_ar) || n_sig_ar < MIN_NSIG_AR) {
    return("NO_BASELINE")
  }
  if (is.na(delta_beta) || is.na(n_sig_mod) || n_sig_mod < MIN_NSIG_MOD) {
    return("AR_LIKE_OR_UNDERPOWERED")
  }
  if (abs(delta_beta) < MIN_DELTA_MEANINGFUL) return("AR_LIKE")
  same_sign <- sign(delta_beta) == sign(beta_ar)
  if (same_sign) return("AMPLIFY")
  # Opposite sign — buffer for clear reversal vs attenuation
  if (abs(delta_beta) > abs(beta_ar) + MIN_DELTA_MEANINGFUL) return("REVERSE")
  return("ATTENUATE")
}

classify_baseline_direction <- function(beta_ar, n_sig_ar) {
  # Use sign(β_AR) as long as n_sig is sufficient — small magnitudes with many
  # sig nhoods are still meaningful directional signal.
  if (is.na(beta_ar) || is.na(n_sig_ar) || n_sig_ar < MIN_NSIG_AR) {
    return("UNDERPOWERED")
  }
  if (abs(beta_ar) < 0.02) return("NULL_AT_BASELINE")  # near-zero only
  if (beta_ar > 0) return("ENRICHED")
  return("DEPLETED")
}

# ------------------------------------------------------------------------
# Per-prediction evaluation
# ------------------------------------------------------------------------
evaluate_prediction <- function(cell_type_pattern, signature_name,
                                  stratum, predicted_class) {
  sig_genes <- SIGNATURES[[signature_name]]

  # === 1. DIRECTION (Stage E, per L2 → aggregate to cell_type pattern) ===
  matching_l2 <- e_long %>%
    filter(grepl(cell_type_pattern, label, ignore.case = TRUE) |
           grepl(cell_type_pattern, L2_joint, ignore.case = TRUE))

  ar_rows <- matching_l2 %>% filter(contrast == AR)
  mod_rows <- matching_l2 %>% filter(contrast == stratum)

  # Aggregate L2s by weighted (n_sig-weighted) mean of med_lfc
  agg_lfc <- function(df) {
    d <- df %>% filter(!is.na(med_lfc), !is.na(n_sig), n_sig > 0)
    if (nrow(d) == 0) return(list(mean_lfc = NA_real_, total_n_sig = 0,
                                    n_l2 = 0))
    list(mean_lfc = sum(d$med_lfc * d$n_sig) / sum(d$n_sig),
         total_n_sig = sum(d$n_sig),
         n_l2 = nrow(d))
  }
  ar_agg <- agg_lfc(ar_rows)
  mod_agg <- agg_lfc(mod_rows)

  beta_ar <- ar_agg$mean_lfc
  n_sig_ar <- ar_agg$total_n_sig
  delta_beta <- mod_agg$mean_lfc
  n_sig_mod <- mod_agg$total_n_sig

  # Direction class
  direction_status <- if (stratum == AR) {
    classify_baseline_direction(beta_ar, n_sig_ar)
  } else {
    classify_modifier_direction(beta_ar, n_sig_ar, delta_beta, n_sig_mod)
  }

  # === 2. BIOLOGY (F.3 markers in modifier contrast, this cell-type pattern) ===
  bio_match <- markers %>%
    filter(contrast == stratum,
           grepl(cell_type_pattern, parent_label, ignore.case = TRUE) |
           grepl(cell_type_pattern, parent_L2_joint, ignore.case = TRUE)) %>%
    filter(.data[[sym_col]] %in% sig_genes | gene_id %in% sig_genes)

  bio_mean_lfc <- if (nrow(bio_match) > 0) mean(bio_match$logFC, na.rm = TRUE) else NA_real_
  bio_n_genes <- length(unique(c(bio_match[[sym_col]], bio_match$gene_id)))
  bio_top_genes <- if (nrow(bio_match) > 0) {
    bio_match %>% arrange(desc(abs(logFC))) %>% head(5) %>%
      pull(.data[[sym_col]]) %>% paste(collapse = ",")
  } else NA_character_
  bio_status <- if (is.na(bio_mean_lfc)) "GENES_NOT_IN_MARKERS" else
                  if (abs(bio_mean_lfc) > 0.3) "STRONG" else
                  if (abs(bio_mean_lfc) > 0.15) "PRESENT" else "WEAK"

  # === Combined evidence vs predicted_class ===
  combined <- "INCONCLUSIVE"
  if (predicted_class %in% c("ENRICHED_PARITY", "DEPLETED_PARITY", "GRADED_SHIFT")) {
    expected <- if (predicted_class == "ENRICHED_PARITY") "ENRICHED" else
                if (predicted_class == "DEPLETED_PARITY") "DEPLETED" else "GRADED"
    combined <- if (expected == "GRADED" && direction_status %in% c("ENRICHED", "DEPLETED")) "PRESENT" else
                if (direction_status == expected) "CONFIRMS" else
                if (direction_status == "UNDERPOWERED") "UNDERPOWERED" else "CONTRADICTS_OR_NULL"
  } else if (predicted_class == "AMPLIFY") {
    combined <- if (direction_status == "AMPLIFY") "CONFIRMS" else
                if (direction_status == "ATTENUATE") "ATTENUATES" else
                if (direction_status == "REVERSE") "REVERSES" else
                direction_status
  } else if (predicted_class == "ATTENUATE") {
    combined <- if (direction_status == "ATTENUATE") "CONFIRMS" else
                if (direction_status == "AMPLIFY") "AMPLIFIES_NOT_ATTENUATE" else
                if (direction_status == "REVERSE") "REVERSES" else
                direction_status
  } else if (predicted_class == "AR_LIKE_SMALLER") {
    combined <- if (direction_status %in% c("AR_LIKE", "ATTENUATE")) "CONFIRMS" else
                direction_status
  } else if (predicted_class == "ENRICHED_BRCA2") {
    combined <- if (direction_status == "AMPLIFY") "CONFIRMS" else direction_status
  } else if (predicted_class %in% c("MIXED_DEPENDS_OCCULT", "TREATMENT_RESIDUAL")) {
    # Compare modifier delta sign+magnitude vs BRCA1's modifier
    brca1_mod_rows <- e_long %>%
      filter(contrast == "parity_x_HR_BRCA1") %>%
      filter(grepl(cell_type_pattern, label, ignore.case = TRUE) |
             grepl(cell_type_pattern, L2_joint, ignore.case = TRUE))
    brca1_agg <- agg_lfc(brca1_mod_rows)
    if (!is.na(delta_beta) && !is.na(brca1_agg$mean_lfc) && !is.na(beta_ar)) {
      d_to_brca1 <- abs(delta_beta - brca1_agg$mean_lfc)
      d_to_ar <- abs(delta_beta - 0)  # AR's delta is by definition 0
      combined <- if (d_to_brca1 < d_to_ar) "CLUSTERS_WITH_BRCA1" else "CLUSTERS_WITH_AR"
    } else combined <- "INCONCLUSIVE"
  }

  tibble(
    prediction = sprintf("%s × %s × %s",
                          cell_type_pattern, signature_name, stratum),
    cell_type_pattern, signature = signature_name, stratum,
    predicted_class,
    n_l2_matched = if (stratum == AR) ar_agg$n_l2 else mod_agg$n_l2,
    beta_ar = beta_ar, n_sig_ar = n_sig_ar,
    delta_beta = delta_beta, n_sig_mod = n_sig_mod,
    direction_status, bio_mean_lfc,
    bio_n_genes_present = bio_n_genes, bio_status,
    bio_top_genes = bio_top_genes,
    evidence_status = combined
  )
}

cat("\nEvaluating predictions...\n")
results <- pmap_dfr(PREDICTIONS, function(cell_type_pattern, signature, stratum, predicted_class) {
  evaluate_prediction(cell_type_pattern, signature, stratum, predicted_class)
})
cat(sprintf("Predictions evaluated: %d\n", nrow(results)))

# Output structured CSV
out_csv <- file.path(out_dir, "predictions_evidence.csv")
write_csv(results, out_csv)
cat(sprintf("Wrote: %s\n", out_csv))

# ------------------------------------------------------------------------
# Heatmap visual: per-prediction direction × biology
# ------------------------------------------------------------------------
heat_df <- results %>%
  mutate(label = paste(cell_type_pattern, signature, sep = " — "),
         strat_short = sub("parity_x_HR_|parity_in_", "", stratum),
         display_lfc = ifelse(stratum == AR, beta_ar, delta_beta))
heat_df$display_lfc_capped <- pmin(pmax(heat_df$display_lfc, -1.5), 1.5)

ord <- heat_df %>% filter(strat_short == "AR") %>%
  arrange(beta_ar) %>% pull(label)
heat_df$label <- factor(heat_df$label, levels = unique(c(ord, heat_df$label)))
heat_df$strat_short <- factor(heat_df$strat_short,
                                levels = c("AR", "BRCA1", "sporadic", "BRCA2"))

p <- ggplot(heat_df, aes(x = strat_short, y = label, fill = display_lfc_capped)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = evidence_status), size = 2.3, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "grey92", high = "#B2182B",
                       midpoint = 0, na.value = "grey95",
                       limits = c(-1.5, 1.5), oob = scales::squish,
                       name = "AR: β_AR | mod: Δβ\n(weighted by n_sig)") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank()) +
  labs(title = "Predictions vs evidence — DIRECTION (Stage E β_AR + interaction Δβ)",
       subtitle = paste("Fill = β_AR for AR column, interaction Δβ for modifier",
                         "columns. Text = combined-evidence vs prediction.",
                         sep = "\n"),
       x = NULL, y = NULL)

out_pdf <- file.path(fig_dir, "predictions_signature_heatmap.pdf")
ggsave(out_pdf, p, width = 9,
       height = max(8, min(18, 0.30 * nlevels(heat_df$label) + 3)),
       limitsize = FALSE)
cat(sprintf("Wrote: %s\n", out_pdf))

# ------------------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------------------
md_lines <- c("# Predictions Evidence — corrected interaction interpretation",
              "",
              "## Framework",
              "",
              "- `parity_in_AR.med_lfc` = β_AR (parity main effect within AR cohort)",
              "- `parity_x_HR_*.med_lfc` = Δβ_modifier (interaction, relative to AR)",
              "- **AMPLIFY**: sign(Δβ) == sign(β_AR), |Δβ| > 0.15 → modifier reinforces AR effect",
              "- **ATTENUATE**: sign(Δβ) == −sign(β_AR), |Δβ| < |β_AR|",
              "- **REVERSE**: sign(Δβ) == −sign(β_AR), |Δβ| > |β_AR|",
              "- **AR_LIKE**: |Δβ| < 0.15 (no meaningful modifier effect)",
              "",
              "Aggregated across L2s matching the cell-type pattern, weighted by n_sig.",
              "")
for (st in c(AR, "parity_x_HR_BRCA1", "parity_x_HR_sporadic", "parity_x_HR_BRCA2")) {
  sub <- results %>% filter(stratum == st)
  if (nrow(sub) == 0) next
  md_lines <- c(md_lines, sprintf("### %s", st), "",
                "| Cell type | Signature | Predicted | β_AR | Δβ_mod | n_sig_mod | Direction | Bio status (signature in markers) | Evidence |",
                "|---|---|---|---:|---:|---:|---|---|---|")
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    md_lines <- c(md_lines, sprintf(
      "| %s | %s | %s | %s | %s | %s | %s | %s%s | **%s** |",
      r$cell_type_pattern, r$signature, r$predicted_class,
      ifelse(is.na(r$beta_ar), "—", sprintf("%+.2f", r$beta_ar)),
      # Δβ_mod is the same as β_AR for AR-baseline rows (cosmetic — show "—")
      ifelse(st == AR || is.na(r$delta_beta), "—",
             sprintf("%+.2f", r$delta_beta)),
      ifelse(is.na(r$n_sig_mod), "—", as.character(r$n_sig_mod)),
      r$direction_status,
      r$bio_status,
      ifelse(is.na(r$bio_top_genes), "",
             sprintf(" (%s)", r$bio_top_genes)),
      r$evidence_status))
  }
  md_lines <- c(md_lines, "")
}

md_path <- file.path(out_dir, "predictions_evidence.md")
writeLines(md_lines, md_path)
cat(sprintf("Wrote: %s\n", md_path))

# ------------------------------------------------------------------------
# Per-L2 breakdown for focused findings (avoid pattern aggregation)
# ------------------------------------------------------------------------
focus_patterns <- c("LASP", "BMYO|Basal", "Macro", "Fibro", "CD8|T_")
md2 <- c("# Per-L2 Breakdown — finding-level granularity",
         "",
         "Stage E β_AR + interaction Δβ shown per L2, not aggregated.",
         "Helps disambiguate which sub-L2 (LASP-major vs LASP-basal vs LASP-KIT)",
         "drives a pattern-level direction.",
         "",
         "Columns: β_AR (parity_in_AR med_lfc) | n_sig_AR | Δβ for each modifier contrast | n_sig_mod.",
         "")

per_l2_rows <- list()
for (pat in focus_patterns) {
  pat_rows <- e_long %>%
    filter(grepl(pat, label, ignore.case = TRUE) |
           grepl(pat, L2_joint, ignore.case = TRUE))
  if (nrow(pat_rows) == 0) next
  l2s <- unique(pat_rows$L2_joint)
  md2 <- c(md2, sprintf("## Pattern: %s (%d L2s)", pat, length(l2s)), "")
  md2 <- c(md2,
            "| L2 | β_AR | n_sig_AR | Δβ_BRCA1 | n_sig_BRCA1 | Δβ_sporadic | n_sig_spor | Δβ_BRCA2 | n_sig_BRCA2 | direction (BRCA1) |",
            "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|")
  for (l2 in sort(l2s)) {
    rows <- pat_rows %>% filter(L2_joint == l2)
    get_v <- function(con, col) {
      x <- rows %>% filter(contrast == con)
      if (nrow(x) == 0) NA else x[[col]][1]
    }
    bAR <- get_v(AR, "med_lfc"); nAR <- get_v(AR, "n_sig")
    db1 <- get_v("parity_x_HR_BRCA1", "med_lfc")
    n1  <- get_v("parity_x_HR_BRCA1", "n_sig")
    dsp <- get_v("parity_x_HR_sporadic", "med_lfc")
    nsp <- get_v("parity_x_HR_sporadic", "n_sig")
    db2 <- get_v("parity_x_HR_BRCA2", "med_lfc")
    n2  <- get_v("parity_x_HR_BRCA2", "n_sig")
    dir1 <- classify_modifier_direction(bAR, nAR, db1, n1)

    fmt <- function(x) ifelse(is.na(x), "—", sprintf("%+.2f", x))
    fmt_n <- function(x) ifelse(is.na(x), "—", as.character(as.integer(x)))
    md2 <- c(md2, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                            l2, fmt(bAR), fmt_n(nAR),
                            fmt(db1), fmt_n(n1),
                            fmt(dsp), fmt_n(nsp),
                            fmt(db2), fmt_n(n2),
                            dir1))
    per_l2_rows[[length(per_l2_rows) + 1]] <- tibble(
      pattern = pat, L2_joint = l2,
      beta_AR = bAR, n_sig_AR = nAR,
      delta_BRCA1 = db1, n_sig_BRCA1 = n1, dir_BRCA1 = dir1,
      delta_sporadic = dsp, n_sig_sporadic = nsp,
      delta_BRCA2 = db2, n_sig_BRCA2 = n2)
  }
  md2 <- c(md2, "")
}

md2_path <- file.path(out_dir, "predictions_per_L2_breakdown.md")
writeLines(md2, md2_path)
cat(sprintf("Wrote: %s\n", md2_path))

per_l2_df <- bind_rows(per_l2_rows)
per_l2_csv <- file.path(out_dir, "predictions_per_L2_breakdown.csv")
write_csv(per_l2_df, per_l2_csv)
cat(sprintf("Wrote: %s\n", per_l2_csv))

cat("\n=== surface_interactions done ===\n")
