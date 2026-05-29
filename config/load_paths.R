# =============================================================================
# load_paths.R — R loader for publication path resolution
# =============================================================================
# Reads paths.yaml and resolves all paths for the current environment (HPC or
# local). Every analysis, rendering, and assembly script sources this file
# instead of hardcoding absolute paths.
#
# Usage:
#   source("publication/config/load_paths.R")
#   paths <- load_paths()
#   h5ad  <- resolve_path("ihbca_integrated", paths)
#   sif   <- get_container("r_spatial_4.3.3", paths)
#   binds <- get_bind_mounts(paths)
# =============================================================================

suppressPackageStartupMessages(library(yaml))

# --- Environment detection ---

#' Detect whether we are running on HPC or locally.
#' @return "hpc" or "local"
detect_environment <- function() {
  env_override <- Sys.getenv("IHBCA_ENV", "")
  if (nzchar(env_override)) {
    return(tolower(env_override))
  }
  node <- Sys.info()[["nodename"]]
  if (grepl("hpc3|^login-|^compute-", node, ignore.case = TRUE)) "hpc" else "local"
}

# --- Config location ---

find_paths_yaml <- function(config_dir = NULL) {
  if (!is.null(config_dir)) {
    yaml_path <- file.path(config_dir, "paths.yaml")
    if (file.exists(yaml_path)) return(yaml_path)
    stop("paths.yaml not found at: ", yaml_path)
  }

  candidates <- character(0)

  # From environment variable
  project_root <- Sys.getenv("IHBCA_PROJECT_ROOT", "")
  if (nzchar(project_root)) {
    candidates <- c(candidates,
      file.path(project_root, "publication", "config", "paths.yaml"))
  }

  # Relative to calling script
  calling_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(calling_file)) {
    candidates <- c(candidates,
      file.path(dirname(calling_file), "paths.yaml"))
  }

  # Common working directory candidates
  candidates <- c(candidates,
    "publication/config/paths.yaml",
    file.path("..", "config", "paths.yaml"),
    file.path("..", "..", "config", "paths.yaml"),
    file.path("..", "..", "..", "config", "paths.yaml"))

  yaml_path <- Find(file.exists, candidates)
  if (is.null(yaml_path)) {
    stop("Cannot find paths.yaml. Set IHBCA_PROJECT_ROOT or pass config_dir explicitly.")
  }
  normalizePath(yaml_path, mustWork = TRUE)
}

# --- Core loader ---

#' Load and resolve all paths from paths.yaml
#' @param config_dir Path to config directory (default: auto-detect)
#' @return Nested list with all paths resolved for current environment
load_paths <- function(config_dir = NULL) {
  yaml_path <- find_paths_yaml(config_dir)
  config <- yaml::read_yaml(yaml_path)

  env <- detect_environment()
  root <- config$roots[[env]]
  if (is.null(root)) {
    stop("No root defined for environment '", env, "' in paths.yaml")
  }

  # Resolve sources: absolute paths stay absolute, relative paths get root prepended
  if (!is.null(config$sources)) {
    config$sources <- lapply(config$sources, function(p) {
      if (startsWith(p, "/")) p else file.path(root, p)
    })
  }

  # Resolve projects
  if (!is.null(config$projects)) {
    config$projects <- lapply(config$projects, function(p) {
      if (startsWith(p, "/")) p else file.path(root, p)
    })
  }

  config$.environment <- env
  config$.root <- root
  config$.yaml_path <- yaml_path
  config$.config_dir <- dirname(yaml_path)
  config
}

# --- Convenience accessors ---

#' Resolve a single source path by key
#' @param key Source key from paths.yaml (e.g., "ihbca_integrated")
#' @param paths Pre-loaded config (from load_paths())
#' @return Resolved absolute path string
resolve_path <- function(key, paths = NULL) {
  if (is.null(paths)) paths <- load_paths()
  p <- paths$sources[[key]]
  if (is.null(p)) {
    available <- paste(names(paths$sources), collapse = ", ")
    stop("Unknown source key: '", key, "'. Available: ", available)
  }
  p
}

#' Get container path by name
#' @param name Container key (e.g., "r_spatial_4.3.3")
#' @param paths Pre-loaded config
#' @return Absolute container path
get_container <- function(name, paths = NULL) {
  if (is.null(paths)) paths <- load_paths()
  ctr <- paths$containers[[name]]
  if (is.null(ctr)) {
    available <- paste(names(paths$containers), collapse = ", ")
    stop("Unknown container: '", name, "'. Available: ", available)
  }
  ctr
}

#' Get bind mount strings for singularity --bind
#' @param paths Pre-loaded config
#' @return Character vector of bind mount strings
get_bind_mounts <- function(paths = NULL) {
  if (is.null(paths)) paths <- load_paths()
  unlist(paths$bind_mounts)
}

#' Get R library path for a container
#' @param name R libs key (e.g., "r_spatial_4.3.3")
#' @param paths Pre-loaded config
#' @return Absolute R_LIBS_USER path
get_r_libs <- function(name, paths = NULL) {
  if (is.null(paths)) paths <- load_paths()
  rlib <- paths$r_libs[[name]]
  if (is.null(rlib)) {
    # Fall back to top-level r_libs_user for single-entry configs
    rlib <- paths$r_libs_user
  }
  if (is.null(rlib)) {
    stop("No R library path found for: '", name, "'")
  }
  rlib
}

#' Get SLURM defaults
#' @param paths Pre-loaded config
#' @return Named list with account, partition, default_cpus, default_mem, default_time
get_slurm_defaults <- function(paths = NULL) {
  if (is.null(paths)) paths <- load_paths()
  paths$slurm
}
