# lib/patch_makeContrasts.R
# Monkey-patch limma::makeContrasts to handle interaction-term column names.
#
# Problem: limma::makeContrasts(levels = model_matrix) validates that all level
# names pass make.names(). R's `:` interaction operator produces design-matrix
# column names like "stratumHR_BRCA1:parity_binaryparous" which fail that check.
# miloR's testNhoods builds the model matrix internally and pipes it straight
# to makeContrasts without sanitizing names, so any contrast referencing an
# interaction term fails.
#
# Fix: intercept makeContrasts, translate non-syntactic colnames via make.names()
# (typically `:` -> `.`), apply the same translation (longest-first to avoid
# partial overlaps) to the user's contrast strings, then call the original.
#
# Sourced after patch_milor.R in scripts that call miloR::testNhoods.

local({
  patched_makeContrasts <- function(..., contrasts = NULL, levels = NULL) {
    args <- list(...)
    if (is.null(contrasts) && length(args) >= 1) contrasts <- args[[1]]
    if (is.null(levels)    && length(args) >= 2) levels    <- args[[2]]

    if (is.matrix(levels) && !is.null(colnames(levels))) {
      orig <- colnames(levels)
      safe <- make.names(orig)
      if (!all(safe == orig)) {
        # Sort by orig-name length descending so longer names (interaction
        # terms like "stratumHR_BRCA1:parity_binaryparous") get matched
        # before any substring (e.g. "stratumHR_BRCA1").
        ord <- order(-nchar(orig))
        for (i in ord) {
          if (orig[i] != safe[i]) {
            contrasts <- gsub(orig[i], safe[i], contrasts, fixed = TRUE)
          }
        }
        colnames(levels) <- safe
      }
    }
    # Use the original (unpatched) function to avoid recursion.
    .orig_makeContrasts(contrasts = contrasts, levels = levels)
  }

  # Stash the original once, then install the patched version.
  if (!exists(".orig_makeContrasts", envir = .GlobalEnv, inherits = FALSE)) {
    assign(".orig_makeContrasts",
           getFromNamespace("makeContrasts", "limma"),
           envir = .GlobalEnv)
  }
  assignInNamespace("makeContrasts", patched_makeContrasts, ns = "limma")

  # Also update the exported binding on the search path so callers that
  # resolve via package:limma pick up the patch.
  pkg_env <- as.environment("package:limma")
  unlockBinding("makeContrasts", pkg_env)
  assign("makeContrasts", patched_makeContrasts, envir = pkg_env)
  lockBinding("makeContrasts", pkg_env)

  message("[patch_makeContrasts] limma::makeContrasts wrapped to translate ",
          "non-syntactic column names (e.g. interaction `:` -> `.`)")
})
