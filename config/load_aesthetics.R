# =============================================================================
# load_aesthetics.R — R loader for Submission v1 aesthetic rules
# =============================================================================
# Reads aesthetics.yaml and provides helper functions for consistent rendering.
# Every panel rendering script MUST source this file.
#
# Usage:
#   source("publication/config/load_aesthetics.R")
#   aes <- load_aesthetics()
#   pal <- get_palette("l1")
#   theme <- get_theme()
#   dims <- get_dimensions("umap")
# =============================================================================

suppressPackageStartupMessages({
  library(yaml)
  library(ggplot2)
})

# Force UTF-8 locale so yaml::read_yaml correctly parses non-ASCII characters
# (em-dashes etc.) in aesthetics.yaml. Without this, R on POSIX/C-locale
# containers warns "invalid input found on input connection" and silently
# returns an empty parse, breaking every render-dev script that calls get_palette().
Sys.setlocale("LC_ALL", "C.UTF-8")

# --- Core loader ---

#' Load aesthetics config from YAML
#' @param config_dir Path to config directory (default: auto-detect)
#' @return Parsed YAML as nested list
load_aesthetics <- function(config_dir = NULL) {
  if (is.null(config_dir)) {
    # Auto-detect: look for aesthetics.yaml relative to common locations
    candidates <- c(
      file.path(Sys.getenv("PROJECT_ROOT", ""), "publication/config/aesthetics.yaml"),
      "publication/config/aesthetics.yaml",
      file.path(Sys.getenv("PROJECT_ROOT", ""), "Analysis/Submission_v1/config/aesthetics.yaml"),
      "Analysis/Submission_v1/config/aesthetics.yaml",
      file.path(dirname(sys.frame(1)$ofile %||% "."), "aesthetics.yaml")
    )
    yaml_path <- Find(file.exists, candidates)
    if (is.null(yaml_path)) {
      stop("Cannot find aesthetics.yaml. Set PROJECT_ROOT or pass config_dir explicitly.")
    }
  } else {
    yaml_path <- file.path(config_dir, "aesthetics.yaml")
  }

  config <- yaml::read_yaml(yaml_path)
  config$.source_path <- yaml_path
  config$.config_dir <- dirname(yaml_path)
  config
}

# --- Palette helpers ---

#' Get a named color vector for a label level
#' @param level One of "compartment", "l1", "l2", "anatomic_axes", "patient"
#' @param config Optional pre-loaded config (avoids re-reading YAML)
#' @return Named character vector of hex colors
get_palette <- function(level, config = NULL) {
  if (is.null(config)) config <- load_aesthetics()
  pal <- config$palettes[[level]]
  if (is.null(pal)) {
    stop(sprintf("Unknown palette level: '%s'. Available: %s",
                 level, paste(names(config$palettes), collapse = ", ")))
  }
  unlist(pal)
}

#' Get display labels for L2 types from summary_idents.csv
#' @param label_style "short" or "descriptive" (default from config)
#' @param config Optional pre-loaded config
#' @return Named vector: canonical key -> display label
get_display_labels <- function(label_style = NULL, config = NULL) {
  if (is.null(config)) config <- load_aesthetics()

  # Find summary_idents.csv in promoted stages
  stages_dir <- file.path(dirname(config$.config_dir), "stages", "S2_behavioral_groups")
  idents_path <- file.path(stages_dir, config$labels$source)

  if (!file.exists(idents_path)) {
    warning("summary_idents.csv not found at ", idents_path,
            ". Returning canonical keys as labels.")
    pal <- get_palette("l2", config)
    return(setNames(names(pal), names(pal)))
  }

  idents <- read.csv(idents_path, stringsAsFactors = FALSE)
  if (is.null(label_style)) {
    label_style <- config$typography$l2_labels  # "short" or "descriptive"
  }

  col <- if (label_style == "short") {
    config$labels$short_column
  } else {
    config$labels$descriptive_column
  }

  setNames(idents[[col]], idents[[config$labels$key_column]])
}

# --- Dimension helpers ---

#' Get dimensions for a panel type
#' @param panel_type One of the keys in dimensions.panel_types
#' @param config Optional pre-loaded config
#' @return List with width, height, and type-specific params
get_dimensions <- function(panel_type, config = NULL) {
  if (is.null(config)) config <- load_aesthetics()
  dims <- config$dimensions$panel_types[[panel_type]]
  if (is.null(dims)) {
    stop(sprintf("Unknown panel type: '%s'. Available: %s",
                 panel_type, paste(names(config$dimensions$panel_types), collapse = ", ")))
  }
  dims
}

# --- Theme ---

#' Get a ggplot2 theme conforming to Nature typography rules
#' @param config Optional pre-loaded config
#' @return ggplot2 theme object
get_theme <- function(config = NULL) {
  if (is.null(config)) config <- load_aesthetics()
  typo <- config$typography

  theme_classic(base_size = typo$base_size, base_family = typo$family) +
    theme(
      plot.title = element_text(size = typo$title_size, face = "bold"),
      axis.text = element_text(size = typo$axis_text_size),
      axis.title = element_text(size = typo$axis_title_size),
      legend.text = element_text(size = typo$legend_size),
      legend.title = element_text(size = typo$legend_title_size),
      strip.text = element_text(size = typo$strip_text_size),
      legend.key.size = unit(0.3, "cm"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# --- Validation ---

#' Validate a panel before saving
#' Checks density limits, font sizes, and rendering rules.
#' @param plot A ggplot2 object (or NULL for non-ggplot panels)
#' @param panel_type Panel type key
#' @param n_groups Number of groups/categories in the plot
#' @param n_cells Number of cells (for UMAPs)
#' @param config Optional pre-loaded config
#' @return Invisible TRUE if valid, warns on violations
validate_panel <- function(plot = NULL, panel_type, n_groups = NULL,
                           n_cells = NULL, config = NULL) {
  if (is.null(config)) config <- load_aesthetics()
  density <- config$density
  issues <- character(0)

  # Density checks
  if (!is.null(n_cells) && grepl("umap", panel_type)) {
    if (n_cells > density$umap_max_cells) {
      issues <- c(issues, sprintf(
        "UMAP has %d cells (limit %d). Consider stratified downsampling.",
        n_cells, density$umap_max_cells))
    }
  }

  if (!is.null(n_groups)) {
    if (panel_type == "dotplot" && n_groups > density$dotplot_max_groups) {
      issues <- c(issues, sprintf(
        "Dotplot has %d groups (limit %d). Split into sub-panels.",
        n_groups, density$dotplot_max_groups))
    }
    if (grepl("heatmap", panel_type) && n_groups > density$heatmap_max_groups) {
      issues <- c(issues, sprintf(
        "Heatmap has %d groups (limit %d). Split into sub-panels.",
        n_groups, density$heatmap_max_groups))
    }
    if (panel_type == "stacked_bar" && n_groups > density$stacked_bar_max_categories) {
      issues <- c(issues, sprintf(
        "Stacked bar has %d categories (limit %d). Combine rare into Other.",
        n_groups, density$stacked_bar_max_categories))
    }
  }

  # Rendering checks
  if (!isFALSE(config$rendering$raster)) {
    issues <- c(issues, "rendering.raster must be false for publication")
  }

  if (length(issues) > 0) {
    for (iss in issues) warning("[validate_panel] ", iss)
  } else {
    message("[validate_panel] All checks passed for panel type: ", panel_type)
  }

  invisible(length(issues) == 0)
}

# --- Save helper ---

#' Save a panel as PDF + PNG with correct dimensions
#' @param plot A ggplot2 object
#' @param panel_id Panel identifier (e.g., "l1_umap")
#' @param panel_type Panel type key (for dimensions lookup)
#' @param output_dir Directory for output files
#' @param config Optional pre-loaded config
#' @return List of output file paths
save_panel <- function(plot, panel_id, panel_type, output_dir, config = NULL) {
  if (is.null(config)) config <- load_aesthetics()
  dims <- get_dimensions(panel_type, config)
  dpi <- config$rendering$dpi

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pdf_path <- file.path(output_dir, paste0(panel_id, ".pdf"))
  png_path <- file.path(output_dir, paste0(panel_id, ".png"))

  ggsave(pdf_path, plot = plot, width = dims$width, height = dims$height,
         units = "in", device = cairo_pdf)
  ggsave(png_path, plot = plot, width = dims$width, height = dims$height,
         units = "in", dpi = dpi, bg = "white")

  message(sprintf("[save_panel] Saved %s: PDF=%s, PNG=%s (%gx%g in, %d dpi)",
                  panel_id, pdf_path, png_path, dims$width, dims$height, dpi))

  invisible(list(pdf = pdf_path, png = png_path,
                 width = dims$width, height = dims$height, dpi = dpi))
}
