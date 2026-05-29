# lib/da_helpers.R
# Shared DA helpers used across Stage D + Stage F.

# Adaptive LFC cutoff: |logFC| where LOESS-smoothed -log10(SpatialFDR) crosses
# -log10(fdr_threshold). Used to set max.lfc.delta in groupNhoods (script 12 pattern).
# Source: clinical_cohort_summary_20260426/scripts/11_run_da.Rmd:70-100
compute_adaptive_lfc_cutoff <- function(da_results, fdr_threshold = 0.05,
                                         fallback_lfc = 0.5, loess_span = 0.35) {
  eps <- 1e-300
  target_y <- -log10(pmin(pmax(fdr_threshold, eps), 1))
  d <- da_results %>%
    dplyr::filter(is.finite(logFC), is.finite(SpatialFDR)) %>%
    dplyr::mutate(SpatialFDR_clamped = pmin(pmax(SpatialFDR, eps), 1),
                   abs_logFC = abs(logFC),
                   neglogFDR = -log10(SpatialFDR_clamped)) %>%
    dplyr::filter(is.finite(abs_logFC), is.finite(neglogFDR))
  if (nrow(d) < 10) return(fallback_lfc)
  out <- tryCatch({
    fit <- loess(neglogFDR ~ abs_logFC, data = d, span = loess_span,
                  degree = 1, family = "symmetric")
    grid_x <- seq(min(d$abs_logFC), max(d$abs_logFC), length.out = 800)
    pred_y <- as.numeric(predict(fit, newdata = data.frame(abs_logFC = grid_x)))
    fit_df <- data.frame(x = grid_x, y = pred_y)
    fit_df <- fit_df[is.finite(fit_df$y), ]
    cross_idx <- which(fit_df$y >= target_y)[1]
    if (!is.na(cross_idx) && cross_idx > 1) {
      x0 <- fit_df$x[cross_idx - 1]; y0 <- fit_df$y[cross_idx - 1]
      x1 <- fit_df$x[cross_idx];     y1 <- fit_df$y[cross_idx]
      w <- min(max((target_y - y0) / (y1 - y0 + 1e-12), 0), 1)
      x0 + w * (x1 - x0)
    } else if (!is.na(cross_idx)) {
      fit_df$x[cross_idx]
    } else {
      max(d$abs_logFC)
    }
  }, error = function(e) fallback_lfc)
  if (!is.finite(out) || out <= 0) fallback_lfc else out
}
