#!/usr/bin/env Rscript
# 07_genotyped_cohort_fibros.R
# BR1 vs AR DE on the genotype-tested cohort: AR donors restricted to those
# with brca_genotype == "negative" (n=21) + BR1 donors (n=34). Removes the
# untested-AR ascertainment confound — all donors went through clinical
# BRCA testing.
#
# Formula: ~ sample_type + parity_binary + age_binary + risk_class
# (same as A2_sampletype on the full cohort)
#
# Runs DE for the 4 fibroblast L2s and renders CAF-labeled volcanos.

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(dplyr)
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  stop("ggnewscale not installed")
}

p <- ArgumentParser()
p$add_argument("--project-root", required = TRUE)
p$add_argument("--counts",       required = TRUE)
p$add_argument("--pb-meta",      required = TRUE)
p$add_argument("--donor-meta",   required = TRUE)
p$add_argument("--gene-mapping", required = TRUE)
p$add_argument("--fdr", type = "double", default = 0.05)
p$add_argument("--lfc", type = "double", default = 0.5)
args <- p$parse_args()

project_root <- args$project_root
de_dir   <- file.path(project_root, "outputs", "de_results",
                       "C_tested", "A2_sampletype_tested")
plot_dir <- file.path(project_root, "outputs", "plots",
                       "fibro_volcanos_genotyped")
dir.create(de_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

FORMULA_STR <- "~ sample_type + parity_binary + age_binary + risk_class"
fibros <- c("str__Fibro_major", "str__Fibro_IGF1",
            "str__Fibro_SFRP4", "str__Fibro_prematrix")

# ---- load substrate ----------------------------------------------------
cat("[1] load pseudobulk counts ... ")
counts_df <- readr::read_csv(args$counts, show_col_types = FALSE, progress = FALSE)
counts_mat <- as.matrix(counts_df[, -1, drop = FALSE])
rownames(counts_mat) <- counts_df[[1]]
counts_mat <- t(counts_mat)
cat(sprintf("%d genes x %d samples\n", nrow(counts_mat), ncol(counts_mat)))

cat("[2] load pseudobulk meta ... ")
pb_meta <- readr::read_csv(args$pb_meta, show_col_types = FALSE, progress = FALSE)
pb_meta <- pb_meta[match(colnames(counts_mat), pb_meta$sample_id), ]
cat(sprintf("%d rows\n", nrow(pb_meta)))

cat("[3] load donor meta ... ")
donor_meta <- readr::read_csv(args$donor_meta, show_col_types = FALSE,
                               progress = FALSE)
cat(sprintf("%d donors\n", nrow(donor_meta)))

# ---- gene symbol mapping -----------------------------------------------
gm <- readr::read_tsv(args$gene_mapping, show_col_types = FALSE,
                      progress = FALSE) %>% dplyr::rename_with(tolower)
ensg_col   <- intersect(c("ensembl_id","ensembl_gene_id","gene_id","ensembl"), colnames(gm))[1]
sym_col_in <- intersect(c("gene_symbol","symbol","gene_name"), colnames(gm))[1]
gm_uniq <- gm %>% dplyr::group_by(.data[[ensg_col]]) %>%
                  dplyr::slice(1) %>% dplyr::ungroup()
gene_id_to_symbol <- setNames(gm_uniq[[sym_col_in]], gm_uniq[[ensg_col]])

# ---- tested-only cohort filter -----------------------------------------
cat("\n[4] build genotype-tested cohort\n")
tested <- donor_meta %>%
  dplyr::filter(risk_class %in% c("AR", "BR1"),
                 parity_binary %in% c("parous", "nulliparous"),
                 brca_genotype %in% c("negative", "BRCA1"))
n_AR  <- sum(tested$risk_class == "AR")
n_BR1 <- sum(tested$risk_class == "BR1")
cat(sprintf("  n_donors = %d (AR-tested-negative=%d, BR1=%d)\n",
            nrow(tested), n_AR, n_BR1))

# Join donor meta to pb meta
if ("study" %in% colnames(pb_meta)) pb_meta$study <- NULL
pb_meta <- pb_meta %>%
  dplyr::left_join(donor_meta, by = c("donor" = "ihbca_donor_id"))

# ---- per-fibro voom DE -------------------------------------------------
results <- list()
for (L2 in fibros) {
  cat(sprintf("\n=== %s ===\n", L2))
  keep <- pb_meta$cell_type == L2 &
          pb_meta$donor %in% tested$ihbca_donor_id &
          pb_meta$n_cells >= 5 &
          !is.na(pb_meta$risk_class)
  sub <- pb_meta[keep, , drop = FALSE]
  nA <- sum(sub$risk_class == "AR"); nB <- sum(sub$risk_class == "BR1")
  cat(sprintf("  AR=%d  BR1=%d  total=%d\n", nA, nB, nA + nB))
  if (nA < 5 || nB < 5) {
    cat("  insufficient donors; skipping\n")
    next
  }
  sub$risk_class <- factor(sub$risk_class, levels = c("AR", "BR1"))
  for (col in c("sample_type", "parity_binary", "age_binary")) {
    sub[[col]] <- factor(sub[[col]])
  }
  design <- model.matrix(as.formula(FORMULA_STR), data = sub)
  if (qr(design)$rank < ncol(design)) {
    cat("  design rank-deficient; skipping\n")
    next
  }
  cm <- round(counts_mat[, sub$sample_id, drop = FALSE])
  keep_col <- colSums(cm) > 0
  cm <- cm[, keep_col, drop = FALSE]
  sub <- sub[keep_col, , drop = FALSE]
  design <- design[keep_col, , drop = FALSE]
  if (sum(rowSums(cm > 0) >= 3) < 200) {
    cat("  <200 expressed genes; skipping\n"); next
  }
  dge <- DGEList(counts = cm)
  dge <- dge[rowSums(dge$counts) >= 10, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  v   <- voom(dge, design = design)
  fit <- lmFit(v, design = design); fit <- eBayes(fit, robust = TRUE)
  tt  <- topTable(fit, coef = "risk_classBR1", n = Inf, sort.by = "none")
  tt$gene_id <- rownames(tt)
  tt$symbol  <- gene_id_to_symbol[tt$gene_id]
  tt$FDR     <- tt$adj.P.Val; tt$PValue <- tt$P.Value
  tt$L2 <- L2; tt$n_AR <- nA; tt$n_BR1 <- nB
  out_csv <- file.path(de_dir, paste0(gsub("[^A-Za-z0-9_]", "_", L2), ".csv"))
  readr::write_csv(tt, out_csv)
  cat(sprintf("  wrote %s  (n_sig FDR<0.05: %d)\n",
              out_csv, sum(tt$adj.P.Val < 0.05, na.rm = TRUE)))
  results[[L2]] <- tt
}

# ---- volcanos (same aesthetic as 06_volcano_fibros_caf_labels.R) -------
ICAF   <- c("IL6","CXCL1","CXCL2","CXCL3","CXCL8","CXCL12","CXCL14",
            "CCL2","CCL7","LIF","PTGS2","C3","HAS1","HAS2",
            "MMP1","MMP3","TIMP1","SERPINE1","PDPN","DPT","DPP4","ICAM1")
MYCAF  <- c("ACTA2","TAGLN","MYH11","MCAM","RGS5","CNN1")
MATCAF <- c("POSTN","COMP","COL1A1","COL3A1","COL10A1","COL11A1","COL12A1",
            "FN1","ASPN","VCAN","LRRC15")
PANCAF <- c("FAP","THY1","S100A4","PDGFRA","PDGFRB")
OTHER  <- c("LOX","LOXL2","MMP2","MMP11","CCN2","CTGF",
            "SFRP2","SFRP4","GREM1","INHBA")

marker_class <- c(
  setNames(rep("iCAF",    length(ICAF)),   ICAF),
  setNames(rep("myCAF",   length(MYCAF)),  MYCAF),
  setNames(rep("matrix",  length(MATCAF)), MATCAF),
  setNames(rep("pan-CAF", length(PANCAF)), PANCAF),
  setNames(rep("other",   length(OTHER)),  OTHER))
all_markers <- names(marker_class)

CLASS_COLORS <- c("iCAF"="#E69F00","myCAF"="#0072B2","matrix"="#009E73",
                  "pan-CAF"="#CC79A7","other"="#56B4E9")
PALETTE_BG <- c("BR1 up"="#FBE4C2","BR1 down"="#C6DDF0","ns"="#E5E5E5")

pretty_l2 <- function(l2) gsub("_", "-", sub("^[a-z]+__", "", l2))

panel_theme <- function() {
  theme_classic(base_size = 8) +
    theme(panel.grid = element_blank(),
          axis.line = element_line(linewidth = 0.3, color = "black"),
          axis.ticks = element_line(linewidth = 0.25, color = "black"),
          axis.text = element_text(size = 7, color = "black"),
          axis.title = element_text(size = 7.5, color = "black"),
          plot.margin = margin(4, 6, 4, 6),
          legend.position = "none",
          plot.title = element_text(size = 9, face = "bold", hjust = 0,
                                     margin = margin(b = 2)))
}

build_volcano <- function(L2, de) {
  if (is.null(de)) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0, y = 0, label = paste0(L2, ": gated")))
  }
  de <- de %>% dplyr::filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
    dplyr::mutate(
      minus_log10_fdr = -log10(pmax(adj.P.Val, .Machine$double.xmin)),
      sig = dplyr::case_when(
        adj.P.Val < args$fdr & logFC >  args$lfc ~ "BR1 up",
        adj.P.Val < args$fdr & logFC < -args$lfc ~ "BR1 down",
        TRUE                                      ~ "ns"))
  de_bg <- de %>% dplyr::filter(!(symbol %in% all_markers))
  de_mk <- de %>% dplyr::filter(symbol %in% all_markers) %>%
    dplyr::mutate(marker_class = marker_class[symbol])
  ggplot() +
    geom_point(data = de_bg, aes(x = logFC, y = minus_log10_fdr, color = sig),
               size = 0.35, alpha = 0.7, shape = 16) +
    scale_color_manual(values = PALETTE_BG) +
    ggnewscale::new_scale_color() +
    geom_point(data = de_mk, aes(x = logFC, y = minus_log10_fdr, fill = marker_class),
               size = 1.8, alpha = 1, shape = 21, color = "black", stroke = 0.25) +
    scale_fill_manual(values = CLASS_COLORS) +
    geom_hline(yintercept = -log10(args$fdr), linetype = "dashed",
               color = "grey50", linewidth = 0.25) +
    geom_vline(xintercept = c(-args$lfc, args$lfc), linetype = "dashed",
               color = "grey50", linewidth = 0.25) +
    ggrepel::geom_text_repel(data = de_mk %>% dplyr::filter(sig != "ns"),
                              aes(x = logFC, y = minus_log10_fdr, label = symbol),
                              size = 2.2, color = "black", max.overlaps = Inf,
                              segment.size = 0.15, segment.color = "grey40",
                              min.segment.length = 0, box.padding = 0.3,
                              point.padding = 0.2, force = 2, seed = 42) +
    labs(title = pretty_l2(L2),
         x = expression(log[2]~FC~(BR1 / AR)),
         y = expression(-log[10]~FDR)) +
    panel_theme()
}

legend_plot <- function() {
  legend_df <- data.frame(x = seq_along(CLASS_COLORS), y = 0,
    marker_class = factor(names(CLASS_COLORS), levels = names(CLASS_COLORS)))
  ggplot(legend_df, aes(x = x, y = y, fill = marker_class)) +
    geom_point(size = 3, shape = 21, color = "black", stroke = 0.3) +
    scale_fill_manual(values = CLASS_COLORS, name = NULL) +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 7),
          legend.key.size = grid::unit(0.4, "cm"),
          legend.spacing.x = grid::unit(0.3, "cm"))
}

cat("\n[5] render volcanos\n")
plots <- list()
for (L2 in fibros) {
  p <- build_volcano(L2, results[[L2]])
  plots[[L2]] <- p
  out_pdf <- file.path(plot_dir,
                       paste0(gsub("[^A-Za-z0-9_]", "_", L2), "_caf_tested.pdf"))
  ggsave(out_pdf, p, width = 100, height = 100, units = "mm", dpi = 600)
}
combined <- patchwork::wrap_plots(plots, ncol = 2)
final <- combined / legend_plot() + patchwork::plot_layout(heights = c(1, 0.04))
ggsave(file.path(plot_dir, "fibros_caf_tested_grid.pdf"),
       final, width = 200, height = 210, units = "mm", dpi = 600,
       limitsize = FALSE)

cat("\nDone.\n")
