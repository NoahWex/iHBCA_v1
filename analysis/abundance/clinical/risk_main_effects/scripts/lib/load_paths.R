# lib/load_paths.R
# Inquiry-agnostic config loader. All pipeline scripts source this.
#
# Usage:
#   args <- commandArgs(trailingOnly = TRUE)
#   inquiry_dir <- args[which(args == "--inquiry-dir") + 1]
#   cfg <- load_inquiry(inquiry_dir)
#   # cfg$paths, cfg$inquiry, cfg$inquiry_dir

suppressPackageStartupMessages({
  library(yaml)
})

load_inquiry <- function(inquiry_dir) {
  inquiry_dir <- normalizePath(inquiry_dir, mustWork = TRUE)
  paths_file <- file.path(inquiry_dir, "config", "paths.yaml")
  inquiry_file <- file.path(inquiry_dir, "config", "inquiry.yaml")
  if (!file.exists(paths_file)) stop("paths.yaml missing at: ", paths_file)
  if (!file.exists(inquiry_file)) stop("inquiry.yaml missing at: ", inquiry_file)

  paths <- yaml::read_yaml(paths_file)
  inquiry <- yaml::read_yaml(inquiry_file)

  # Resolve dissoc_table relative to config/ if present and not absolute
  dissoc <- paths$inputs$dissoc_table
  if (!is.null(dissoc) && !startsWith(dissoc, "/")) {
    paths$inputs$dissoc_table <- file.path(inquiry_dir, "config", dissoc)
  }

  # Resolve output paths relative to inquiry_root
  inq_root <- paths$inquiry_root
  for (k in names(paths$outputs)) {
    p <- paths$outputs[[k]]
    if (!startsWith(p, "/")) paths$outputs[[k]] <- file.path(inq_root, p)
  }

  # Sanity check inputs
  for (k in setdiff(names(paths$inputs), "dissoc_table")) {
    p <- paths$inputs[[k]]
    if (!file.exists(p)) {
      warning(sprintf("Input '%s' not found at: %s", k, p))
    }
  }

  list(
    inquiry_dir = inquiry_dir,
    paths = paths,
    inquiry = inquiry
  )
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

parse_args_inquiry <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  i <- which(args == "--inquiry-dir")
  if (length(i) == 0) stop("Required argument: --inquiry-dir <path>")
  inquiry_dir <- args[i + 1]
  load_inquiry(inquiry_dir)
}
