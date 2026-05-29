#!/usr/bin/env Rscript
# top_markers_browser.R
# Browseable top-markers view per (parent L2 × contrast × NhoodGroup), with
# brief-signature genes flagged for "drawing the connection" to predicted biology.
#
# Outputs:
#   outputs/top_markers_per_group.csv         long table: top N markers per group
#   outputs/top_markers_by_L2.md              markdown report grouped by parent L2,
#                                              all contrasts visible side by side
#   outputs/marker_overlap_jaccard.csv        cross-contrast Jaccard per parent L2

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(stringr)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

this_script <- sub("--file=", "",
                   grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
if (length(this_script) == 0) this_script <- "top_markers_browser.R"
script_dir <- dirname(normalizePath(this_script, mustWork = FALSE))
source(file.path(script_dir, "lib", "load_paths.R"))

cfg <- parse_args_inquiry()
paths <- cfg$paths
out_dir <- paths$inquiry_root

CONTRASTS <- c("parity_in_AR", "parity_x_HR_BRCA1",
                "parity_x_HR_sporadic", "parity_x_HR_BRCA2")

TOP_N <- 30
JACCARD_TOP <- 50
FDR_T <- 0.1

# Brief signatures — flag these in markers
BRIEF_SIGS <- list(
  milk_bio = c("CSN2", "CSN3", "LALBA", "WAP", "ELF5", "PIP", "MUCL1", "CSN1S1"),
  prolif = c("MKI67", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
              "TOP2A", "CDK1", "PCNA", "BIRC5"),
  basal = c("KRT5", "KRT14", "TP63", "ACTA2", "MYLK"),
  emt_hybrid = c("SNAI2", "VIM", "ZEB1", "ZEB2"),
  cd8_exhaust = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX", "ENTPD1", "CTLA4"),
  treg = c("FOXP3", "IL2RA", "IL10", "IKZF2"),
  m2_lam = c("TREM2", "LIPA", "LPL", "FABP4", "FABP5", "APOE", "APOC1",
              "CD163", "MRC1", "SPP1", "C1QC", "C1QA"),
  precaf = c("MMP3", "MMP14", "MMP1", "POSTN", "TGFB1", "LOX", "CXCL12"),
  er_pr = c("ESR1", "PGR", "GREB1", "TFF1", "TFF3", "AGR2", "FOXA1"),
  mtorc1 = c("RPS6", "EIF4E", "EIF4G1", "EIF4EBP1", "MTOR", "RPL5", "RPL10"),
  brca_ddr = c("BRCA1", "BARD1", "RAD51", "FANCD2", "ATM", "CHEK2", "TP53BP1", "PALB2")
)

# Load all F.3 marker tables ----------------------------------------------
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
cat(sprintf("Markers: %d rows × %d cols across %d contrasts\n",
            nrow(markers), ncol(markers), length(unique(markers$contrast))))

sym_col <- intersect(c("symbol", "Symbol", "gene_symbol"), colnames(markers))[1]
if (is.na(sym_col)) sym_col <- "gene_id"
fdr_col <- intersect(c("FDR", "adj.P.Val", "padj"), colnames(markers))[1]
if (is.na(fdr_col)) fdr_col <- intersect(c("P.Value", "PValue"),
                                            colnames(markers))[1]

# Pull parent L2 mapping from F.1 summaries
parent_map <- list()
for (cname in CONTRASTS) {
  gs_csv <- file.path(paths$outputs$stageF1, cname, "nhood_groups_summary.csv")
  if (!file.exists(gs_csv)) next
  gs <- read_csv(gs_csv, show_col_types = FALSE) %>%
    transmute(NhoodGroup_renamed = as.character(NhoodGroup_renamed),
              parent_L2_joint = as.character(parent_L2_joint),
              parent_compartment = as.character(parent_compartment),
              n_nhoods_in_group = suppressWarnings(as.numeric(n_nhoods_in_group)),
              group_med_lfc = suppressWarnings(as.numeric(group_med_lfc)),
              n_sig = suppressWarnings(as.numeric(n_sig))) %>%
    # Normalize NhoodGroup_renamed to F.3's filesystem-safe form
    # (epi::BMYO-basal_23 -> epi__BMYO_basal_23) so the join matches.
    mutate(NhoodGroup_renamed = gsub("-", "_", gsub("::", "__",
                                                      NhoodGroup_renamed)))
  parent_map[[cname]] <- gs
}
parent_df <- bind_rows(parent_map, .id = "contrast")
markers <- markers %>%
  left_join(parent_df, by = c("contrast", "NhoodGroup_renamed"))

# Tag brief-signature hits
flat_sig <- unlist(BRIEF_SIGS)
sig_lookup <- setNames(rep(names(BRIEF_SIGS), lengths(BRIEF_SIGS)), flat_sig)
markers$brief_sig_hit <- sig_lookup[markers[[sym_col]]]
markers$brief_sig_hit[is.na(markers$brief_sig_hit)] <- ""

# Top N UP and TOP N DOWN per (contrast × group), split by direction so the
# UP-marker signal (which characterizes what the group IS) doesn't get mixed
# with DOWN-marker signal (what the group lacks vs siblings).
markers_sig <- markers %>% filter(.data[[fdr_col]] < FDR_T)
markers_sig$direction <- ifelse(markers_sig$logFC > 0, "UP", "DOWN")

top_up <- markers_sig %>%
  filter(logFC > 0) %>%
  group_by(contrast, NhoodGroup_renamed) %>%
  arrange(desc(logFC)) %>%
  slice_head(n = TOP_N) %>%
  ungroup()
top_down <- markers_sig %>%
  filter(logFC < 0) %>%
  group_by(contrast, NhoodGroup_renamed) %>%
  arrange(logFC) %>%
  slice_head(n = TOP_N) %>%
  ungroup()
top_per_group <- bind_rows(top_up, top_down) %>%
  arrange(parent_compartment, parent_L2_joint, contrast,
          NhoodGroup_renamed, desc(logFC))

out_csv <- file.path(out_dir, "top_markers_per_group.csv")
top_per_group %>%
  select(parent_compartment, parent_L2_joint, contrast, NhoodGroup_renamed,
         group_med_lfc, n_nhoods_in_group, n_sig,
         gene_id, all_of(sym_col), logFC, all_of(fdr_col), brief_sig_hit) %>%
  write_csv(out_csv)
cat(sprintf("Wrote: %s (%d rows)\n", out_csv, nrow(top_per_group)))

# Per-L2 markdown — show all groups across all contrasts side by side ------
md_lines <- c("# Top Markers Browser",
              "",
              sprintf("Top %d markers (by |logFC| at FDR<%.2f) per viable NhoodGroup,",
                       TOP_N, FDR_T),
              "organized by parent L2 then contrast. Brief-signature genes are",
              "**bold** with their signature class in (parens).", "")

l2_groups <- top_per_group %>%
  distinct(parent_compartment, parent_L2_joint) %>%
  arrange(parent_compartment, parent_L2_joint)

for (i in seq_len(nrow(l2_groups))) {
  comp <- l2_groups$parent_compartment[i]
  l2 <- l2_groups$parent_L2_joint[i]
  if (is.na(l2) || is.na(comp)) next
  l2_subset <- top_per_group %>%
    filter(parent_compartment == comp, parent_L2_joint == l2)
  if (nrow(l2_subset) == 0) next
  md_lines <- c(md_lines, sprintf("## %s :: %s", comp, l2), "")

  for (cname in CONTRASTS) {
    contrast_subset <- l2_subset %>% filter(contrast == cname)
    if (nrow(contrast_subset) == 0) next
    groups <- unique(contrast_subset$NhoodGroup_renamed)
    md_lines <- c(md_lines, sprintf("### `%s` (%d viable group(s))", cname, length(groups)), "")

    for (g in groups) {
      g_rows <- contrast_subset %>% filter(NhoodGroup_renamed == g)
      header_info <- g_rows %>%
        summarise(n_nhoods = first(n_nhoods_in_group),
                  med_lfc = first(group_med_lfc),
                  n_sig_nhoods = first(n_sig),
                  .groups = "drop")
      direction_arrow <- if (!is.na(header_info$med_lfc)) {
        if (header_info$med_lfc > 0) "↑" else if (header_info$med_lfc < 0) "↓" else "—"
      } else "—"
      md_lines <- c(md_lines,
                     sprintf("**Group `%s`** %s n_nhoods=%s | group_med_lfc=%s | n_sig_nhoods=%s",
                              g,
                              direction_arrow,
                              header_info$n_nhoods,
                              ifelse(is.na(header_info$med_lfc),
                                     "—", sprintf("%.2f", header_info$med_lfc)),
                              header_info$n_sig_nhoods))
      md_lines <- c(md_lines, "")

      render_section <- function(rows, title) {
        if (nrow(rows) == 0) return(character())
        out <- c(sprintf("_%s (n=%d)_", title, nrow(rows)),
                  "",
                  "| Gene | logFC | FDR | Signature |",
                  "|---|---:|---:|---|")
        for (k in seq_len(nrow(rows))) {
          rr <- rows[k, ]
          sym <- rr[[sym_col]]
          flag <- rr$brief_sig_hit %||% ""
          gene_disp <- if (nchar(flag) > 0) {
            sprintf("**%s** _(%s)_", sym, flag)
          } else sym
          out <- c(out, sprintf("| %s | %+.2f | %.2g | %s |",
                                  gene_disp, rr$logFC,
                                  rr[[fdr_col]], flag))
        }
        c(out, "")
      }

      up_rows <- g_rows %>% filter(logFC > 0) %>% arrange(desc(logFC))
      dn_rows <- g_rows %>% filter(logFC < 0) %>% arrange(logFC)
      md_lines <- c(md_lines,
                     render_section(up_rows, "Top UP markers (group enriched vs L2 siblings)"))
      md_lines <- c(md_lines,
                     render_section(dn_rows, "Top DOWN markers (group depleted vs L2 siblings)"))
    }
  }
}

md_path <- file.path(out_dir, "top_markers_by_L2.md")
writeLines(md_lines, md_path)
cat(sprintf("Wrote: %s\n", md_path))

# Cross-contrast marker overlap (Jaccard top-50) per parent L2 -------------
overlap_rows <- list()
for (i in seq_len(nrow(l2_groups))) {
  comp <- l2_groups$parent_compartment[i]
  l2 <- l2_groups$parent_L2_joint[i]
  if (is.na(l2) || is.na(comp)) next
  l2_markers <- markers_sig %>%
    filter(parent_compartment == comp, parent_L2_joint == l2)
  if (nrow(l2_markers) == 0) next
  group_top <- l2_markers %>%
    group_by(contrast, NhoodGroup_renamed) %>%
    arrange(desc(abs(logFC))) %>%
    slice_head(n = JACCARD_TOP) %>%
    ungroup()

  groups_per_contrast <- group_top %>%
    group_by(contrast, NhoodGroup_renamed) %>%
    summarise(genes = list(unique(.data[[sym_col]])), .groups = "drop")
  if (nrow(groups_per_contrast) < 2) next
  for (a in seq_len(nrow(groups_per_contrast) - 1)) {
    for (b in (a + 1):nrow(groups_per_contrast)) {
      ga <- groups_per_contrast$genes[[a]]
      gb <- groups_per_contrast$genes[[b]]
      jacc <- length(intersect(ga, gb)) / length(union(ga, gb))
      overlap_rows[[length(overlap_rows) + 1]] <- tibble(
        parent_L2 = l2,
        contrast_a = groups_per_contrast$contrast[a],
        group_a = groups_per_contrast$NhoodGroup_renamed[a],
        contrast_b = groups_per_contrast$contrast[b],
        group_b = groups_per_contrast$NhoodGroup_renamed[b],
        jaccard_top50 = jacc,
        n_shared = length(intersect(ga, gb))
      )
    }
  }
}
overlap_df <- bind_rows(overlap_rows)
out_jaccard <- file.path(out_dir, "marker_overlap_jaccard.csv")
write_csv(overlap_df, out_jaccard)
cat(sprintf("Wrote: %s (%d pairs)\n", out_jaccard, nrow(overlap_df)))

cat("\n=== top markers browser done ===\n")
