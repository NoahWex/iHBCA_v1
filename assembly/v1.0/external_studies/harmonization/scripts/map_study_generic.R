#!/usr/bin/env Rscript
# =============================================================================
# Generic Cell ID Mapping Script for Phase A (R version)
# =============================================================================
# Maps any component study to iHBCA reference based on study_manifest.yaml.
# Handles RDS (Seurat), h5ad (via anndata), identity and batch_aware mapping.
#
# Usage:
#   Rscript map_study_generic.R --study kumar
#   Rscript map_study_generic.R --study gray
#
# Outputs per study:
#   - {study}_cells.csv: All mapped cells with {study}_* prefixed columns
#   - mapping_findings.yaml: Comprehensive mapping report
# =============================================================================

suppressPackageStartupMessages({
  library(argparse)
  library(yaml)
  library(data.table)
})

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

load_manifest <- function(path) {
  yaml::read_yaml(path)
}

load_ihbca_inventory <- function(path, dataset_filter) {
  cat("Loading iHBCA inventory...\n")
  df <- fread(path)
  filtered <- df[dataset == dataset_filter]
  cat(sprintf("  Total iHBCA cells: %s\n", format(nrow(df), big.mark = ",")))
  cat(sprintf("  %s cells: %s\n", dataset_filter, format(nrow(filtered), big.mark = ",")))
  return(filtered)
}

load_component_study <- function(path, file_format) {
  cat(sprintf("Loading component study: %s\n", path))

  if (file_format == "rds") {
    obj <- readRDS(path)

    # Handle Seurat v3 vs v5
    if (inherits(obj, "Seurat")) {
      tryCatch({
        obj <- UpdateSeuratObject(obj)
      }, error = function(e) {
        # Already updated or compatible
      })
    }

    cell_ids <- colnames(obj)
    metadata <- as.data.frame(obj@meta.data)
    metadata$`_cell_id` <- rownames(metadata)

    cat(sprintf("  Cells: %s\n", format(length(cell_ids), big.mark = ",")))
    cat(sprintf("  Metadata columns: %d\n", ncol(metadata)))

    return(list(cell_ids = cell_ids, metadata = metadata))

  } else if (file_format == "h5ad") {
    # Use anndata via reticulate
    if (!requireNamespace("reticulate", quietly = TRUE)) {
      stop("reticulate required for h5ad files")
    }
    reticulate::use_condaenv("base", required = FALSE)
    sc <- reticulate::import("scanpy")

    adata <- sc$read_h5ad(path)
    cell_ids <- as.character(adata$obs_names$values)
    metadata <- as.data.frame(reticulate::py_to_r(adata$obs))
    metadata$`_cell_id` <- cell_ids

    cat(sprintf("  Cells: %s\n", format(length(cell_ids), big.mark = ",")))
    cat(sprintf("  Metadata columns: %d\n", ncol(metadata)))

    return(list(cell_ids = cell_ids, metadata = metadata))

  } else {
    stop(sprintf("Unknown format: %s", file_format))
  }
}

# =============================================================================
# MAPPING FUNCTIONS
# =============================================================================

map_identity <- function(component_ids, ihbca_df) {
  # Identity mapping: cell IDs should match exactly
  # FULLY VECTORIZED for performance (714K cells in seconds, not hours)

  cat("  Running vectorized identity mapping...\n")

  ihbca_ids <- ihbca_df$cell_id

  # Vectorized lookup - O(n) using hash table internally
  is_mapped <- component_ids %in% ihbca_ids

  mapped_ids <- component_ids[is_mapped]
  unmapped_ids <- component_ids[!is_mapped]

  n_mapped <- length(mapped_ids)
  n_unmapped <- length(unmapped_ids)

  cat(sprintf("  Mapped: %s, Unmapped: %s\n",
              format(n_mapped, big.mark = ","),
              format(n_unmapped, big.mark = ",")))

  # Return data.table directly for fast CSV building
  mapped_dt <- data.table(
    ihbca_cell_id = mapped_ids,
    component_cell_id = mapped_ids,
    mapping_method = "identity"
  )

  # For unmapped, only need sample for report (max 20)
  # Create full list only if small, otherwise just sample
  if (n_unmapped <= 100) {
    unmapped <- lapply(unmapped_ids, function(cid) {
      list(component_cell_id = cid, reason = "not_in_ihbca")
    })
  } else {
    # Only sample for report
    sample_ids <- unmapped_ids[1:min(20, n_unmapped)]
    unmapped <- lapply(sample_ids, function(cid) {
      list(component_cell_id = cid, reason = "not_in_ihbca")
    })
    attr(unmapped, "total_count") <- n_unmapped
  }

  # Use a simple list structure that report can handle
  # but store the data.table for fast CSV output
  mapped <- list()
  attr(mapped, "dt") <- mapped_dt
  attr(mapped, "count") <- n_mapped

  return(list(mapped = mapped, unmapped = unmapped))
}

map_custom_pal <- function(component_ids, ihbca_df, config) {
  # Custom Pal mapping using direct (patient_id, type, barcode) lookup
  #
  # Component format: {prefix}_{donor}_{type}_{barcode-1}
  #   e.g., N_1105_epi_AAACCTGAGTATGACA-1
  #   prefix: N (normal) or B1 (BRCA1)
  #   donor: numeric code (1105, 0019, 0023)
  #   type: epi or total
  #
  # iHBCA format: {batch}_{barcode-1}
  #   batch types: *_epi, *_mix, *_mix_BR1
  #   patient_id: N_1105, B1_0023, etc.
  #
  # Key insight: Same patient can have same barcode in different batches (epi vs mix)
  # So we need (patient_id, type_category, barcode) for unique lookup
  #
  # Type mapping:
  #   component "epi"   -> iHBCA batch ending in "_epi"
  #   component "total" -> iHBCA batch ending in "_mix" or "_mix_BR1"

  cat("  Running custom Pal mapping (resolving all ambiguity)...\n")

  # Extract barcode and type_category from iHBCA
  ihbca_df$barcode <- sapply(strsplit(ihbca_df$cell_id, "_"), function(x) {
    bc <- tail(x, 1)
    gsub("-1$", "", bc)
  })

  # Determine type_category from batch
  ihbca_df$type_category <- ifelse(
    grepl("_epi$", ihbca_df$batch), "epi",
    ifelse(grepl("_mix", ihbca_df$batch), "mix", "unknown")
  )

  cat(sprintf("  iHBCA type distribution: epi=%d, mix=%d\n",
              sum(ihbca_df$type_category == "epi"),
              sum(ihbca_df$type_category == "mix")))

  # Build (patient_id, type_category, barcode) -> ihbca_cell_id lookup
  ihbca_df$lookup_key <- paste(ihbca_df$patient_id, ihbca_df$type_category,
                                ihbca_df$barcode, sep = "|")

  # Check for remaining duplicates
  dup_keys <- ihbca_df$lookup_key[duplicated(ihbca_df$lookup_key)]
  if (length(dup_keys) > 0) {
    cat(sprintf("  WARNING: %d duplicate (patient, type, barcode) keys\n", length(dup_keys)))
  }

  patient_type_barcode_lookup <- setNames(ihbca_df$cell_id, ihbca_df$lookup_key)

  ihbca_patients <- unique(ihbca_df$patient_id)
  cat(sprintf("  iHBCA patients (%d): %s...\n",
              length(ihbca_patients),
              paste(head(sort(ihbca_patients), 8), collapse = ", ")))

  # Parse component cell IDs
  # Format 1 (N cells): {prefix}_{donor}_{type}_{barcode-1} e.g., N_1105_epi_AAACCTGAG-1
  # Format 2 (B1 cells): {prefix}_{donor}_{barcode-1} e.g., B1_0023_AAACCTGAG-1
  # B1 cells without explicit type map to "mix" in iHBCA
  cat("  Parsing component cell IDs...\n")

  parsed <- lapply(component_ids, function(cid) {
    parts <- strsplit(cid, "_")[[1]]

    if (length(parts) >= 4) {
      # Format 1: 4+ parts (N cells with explicit type)
      prefix <- parts[1]  # N or B1
      donor <- parts[2]   # numeric code
      comp_type <- parts[3]    # epi or total
      barcode_full <- paste(parts[4:length(parts)], collapse = "_")
      barcode <- gsub("-1$", "", barcode_full)

      # Map component type to iHBCA type_category
      type_category <- if (comp_type == "epi") "epi" else "mix"

      # Build expected patient_id: {prefix}_{donor}
      expected_patient <- paste(prefix, donor, sep = "_")

      list(
        cell_id = cid,
        prefix = prefix,
        donor = donor,
        comp_type = comp_type,
        type_category = type_category,
        barcode = barcode,
        expected_patient = expected_patient
      )
    } else if (length(parts) == 3 && parts[1] == "B1") {
      # Format 2: 3 parts (B1 cells without explicit type)
      # B1 cells map to "mix" type in iHBCA
      prefix <- parts[1]  # B1
      donor <- parts[2]   # numeric code
      barcode <- gsub("-1$", "", parts[3])

      expected_patient <- paste(prefix, donor, sep = "_")

      list(
        cell_id = cid,
        prefix = prefix,
        donor = donor,
        comp_type = "total",  # implicit - B1 cells are "total" samples
        type_category = "mix",  # B1 cells map to mix batches in iHBCA
        barcode = barcode,
        expected_patient = expected_patient
      )
    } else {
      NULL
    }
  })

  parsed <- Filter(Negate(is.null), parsed)
  cat(sprintf("  Parsed %s of %s cell IDs\n",
              format(length(parsed), big.mark = ","),
              format(length(component_ids), big.mark = ",")))

  # Summarize by donor
  donor_counts <- table(sapply(parsed, function(x) x$expected_patient))
  cat(sprintf("  Unique donors: %d\n", length(donor_counts)))

  # Direct mapping using (patient_id, type_category, barcode)
  cat("\n  Mapping cells directly (no voting needed)...\n")

  # Vectorized extraction
  cell_ids <- sapply(parsed, function(x) x$cell_id)
  expected_patients <- sapply(parsed, function(x) x$expected_patient)
  type_categories <- sapply(parsed, function(x) x$type_category)
  barcodes <- sapply(parsed, function(x) x$barcode)
  comp_types <- sapply(parsed, function(x) x$comp_type)

  # Build lookup keys
  lookup_keys <- paste(expected_patients, type_categories, barcodes, sep = "|")

  # Look up iHBCA cell IDs
  ihbca_cell_ids <- patient_type_barcode_lookup[lookup_keys]

  # Identify mapped vs unmapped
  is_mapped <- !is.na(ihbca_cell_ids)

  n_mapped <- sum(is_mapped)
  n_unmapped <- sum(!is_mapped)

  cat(sprintf("  Mapped: %s (%.1f%%)\n",
              format(n_mapped, big.mark = ","),
              100 * n_mapped / length(parsed)))
  cat(sprintf("  Unmapped: %s (%.1f%%)\n",
              format(n_unmapped, big.mark = ","),
              100 * n_unmapped / length(parsed)))

  # Analyze unmapped by donor
  if (n_unmapped > 0) {
    unmapped_donors <- table(expected_patients[!is_mapped])
    cat("\n  Unmapped by donor:\n")
    for (donor in names(sort(unmapped_donors, decreasing = TRUE))[1:min(10, length(unmapped_donors))]) {
      total_for_donor <- sum(expected_patients == donor)
      unmapped_for_donor <- unmapped_donors[[donor]]
      cat(sprintf("    %s: %d / %d unmapped (%.1f%%)\n",
                  donor, unmapped_for_donor, total_for_donor,
                  100 * unmapped_for_donor / total_for_donor))
    }
  }

  # Build mapped data.table
  mapped_dt <- data.table(
    ihbca_cell_id = ihbca_cell_ids[is_mapped],
    component_cell_id = cell_ids[is_mapped],
    mapping_method = "custom_pal_direct",
    patient = expected_patients[is_mapped],
    comp_type = comp_types[is_mapped]
  )

  # Sample unmapped for report
  unmapped_idx <- which(!is_mapped)
  if (length(unmapped_idx) > 0) {
    sample_idx <- head(unmapped_idx, 50)

    # Check why unmapped: patient not in iHBCA, or barcode not found?
    unmapped <- lapply(sample_idx, function(i) {
      patient <- expected_patients[i]
      reason <- if (!(patient %in% ihbca_patients)) {
        "patient_not_in_ihbca"
      } else {
        "barcode_type_not_found"
      }
      list(
        component_cell_id = cell_ids[i],
        expected_patient = patient,
        type_category = type_categories[i],
        barcode = barcodes[i],
        reason = reason
      )
    })
    attr(unmapped, "total_count") <- n_unmapped

    # Summarize reasons
    all_reasons <- sapply(unmapped_idx, function(i) {
      patient <- expected_patients[i]
      if (!(patient %in% ihbca_patients)) "patient_not_in_ihbca" else "barcode_type_not_found"
    })
    reason_summary <- table(all_reasons)
    cat("\n  Unmapped reasons:\n")
    for (r in names(reason_summary)) {
      cat(sprintf("    %s: %d\n", r, reason_summary[[r]]))
    }
  } else {
    unmapped <- list()
  }

  # Build return structure
  mapped <- list()
  attr(mapped, "dt") <- mapped_dt
  attr(mapped, "count") <- n_mapped

  # Build donor mapping summary (for report)
  donor_mapping <- setNames(
    as.list(names(donor_counts)),
    names(donor_counts)
  )

  return(list(
    mapped = mapped,
    unmapped = unmapped,
    donor_mapping = donor_mapping
  ))
}

map_batch_aware <- function(component_ids, ihbca_df, cell_id_regex) {
  # Two-pass batch-aware mapping (as implemented for Gray)

  # Build (patient, barcode) lookup
  ihbca_df$barcode <- sapply(strsplit(ihbca_df$cell_id, "_"), function(x) {
    bc <- tail(x, 1)
    gsub("-1$", "", bc)  # Remove -1 suffix
  })

  patient_barcode_lookup <- setNames(ihbca_df$cell_id,
                                      paste(ihbca_df$patient_id, ihbca_df$barcode, sep = "|"))
  ihbca_patients <- unique(ihbca_df$patient_id)

  cat(sprintf("  iHBCA patients: %s\n", paste(sort(ihbca_patients), collapse = ", ")))

  # Parse component cell IDs
  parsed <- lapply(component_ids, function(cid) {
    m <- regmatches(cid, regexec(cell_id_regex, cid))[[1]]
    if (length(m) == 4) {  # Full match + 3 groups
      list(
        cell_id = cid,
        status = m[2],
        batch = m[3],
        barcode = m[4],
        batch_key = paste(m[2], m[3], sep = "-")
      )
    } else {
      NULL
    }
  })
  parsed <- Filter(Negate(is.null), parsed)

  # Pass 1: Discover batch→patient mapping
  batch_to_patient_votes <- list()

  for (cell in parsed) {
    batch_key <- cell$batch_key
    barcode <- cell$barcode

    for (patient in ihbca_patients) {
      key <- paste(patient, barcode, sep = "|")
      if (key %in% names(patient_barcode_lookup)) {
        if (is.null(batch_to_patient_votes[[batch_key]])) {
          batch_to_patient_votes[[batch_key]] <- list()
        }
        batch_to_patient_votes[[batch_key]][[patient]] <-
          (batch_to_patient_votes[[batch_key]][[patient]] %||% 0) + 1
      }
    }
  }

  # Get majority for each batch
  batch_to_patient <- list()
  cat("\n  Batch -> Patient mapping:\n")
  for (batch_key in sort(names(batch_to_patient_votes))) {
    votes <- batch_to_patient_votes[[batch_key]]
    if (length(votes) > 0) {
      majority <- names(which.max(unlist(votes)))
      count <- votes[[majority]]
      batch_to_patient[[batch_key]] <- majority
      cat(sprintf("    %s -> %s (%d cells)\n", batch_key, majority, count))
    }
  }

  # Pass 2: Map all cells
  mapped <- list()
  unmapped <- list()

  for (cell in parsed) {
    batch_key <- cell$batch_key
    barcode <- cell$barcode
    cell_id <- cell$cell_id

    if (is.null(batch_to_patient[[batch_key]])) {
      unmapped[[length(unmapped) + 1]] <- list(
        component_cell_id = cell_id,
        reason = "batch_not_mapped"
      )
      next
    }

    expected_patient <- batch_to_patient[[batch_key]]
    key <- paste(expected_patient, barcode, sep = "|")

    if (key %in% names(patient_barcode_lookup)) {
      ihbca_cell_id <- patient_barcode_lookup[[key]]
      mapped[[length(mapped) + 1]] <- list(
        ihbca_cell_id = ihbca_cell_id,
        component_cell_id = cell_id,
        mapping_method = "batch_aware",
        batch_key = batch_key,
        patient = expected_patient
      )
    } else {
      unmapped[[length(unmapped) + 1]] <- list(
        component_cell_id = cell_id,
        batch_key = batch_key,
        expected_patient = expected_patient,
        reason = "not_in_ihbca_for_patient"
      )
    }
  }

  return(list(mapped = mapped, unmapped = unmapped, batch_mapping = batch_to_patient))
}

# =============================================================================
# REPORTING
# =============================================================================

generate_report <- function(study, config, mapped, unmapped, metadata_cols, extra = list()) {
  # Handle both list and attribute-based counts
  n_mapped <- if (!is.null(attr(mapped, "count"))) attr(mapped, "count") else length(mapped)
  n_unmapped <- if (!is.null(attr(unmapped, "total_count"))) attr(unmapped, "total_count") else length(unmapped)

  total <- n_mapped + n_unmapped
  match_rate <- if (total > 0) 100 * n_mapped / total else 0
  expected_rate <- config$expected_match_rate %||% 100

  # Analyze unmapped reasons (handle empty case)
  if (length(unmapped) > 0) {
    reasons <- sapply(unmapped, function(x) x$reason %||% "unknown")
    reason_counts <- as.list(table(reasons))
  } else {
    reason_counts <- list()
  }

  report <- list(
    study = study,
    extraction_date = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    source_file = config$source_path,
    mapping_type = config$mapping_type,
    summary = list(
      total_component_cells = total,
      mapped = n_mapped,
      unmapped = n_unmapped,
      match_rate_pct = round(match_rate, 2),
      expected_match_rate_pct = expected_rate,
      rate_within_tolerance = abs(match_rate - expected_rate) < 0.5
    ),
    validation = list(
      expected_issues = config$expected_issues %||% list(),
      verified = list(),
      unexpected_gaps = list()
    ),
    unmapped_by_reason = reason_counts,
    metadata_columns = metadata_cols
  )

  # Merge extra info
  for (n in names(extra)) {
    report[[n]] <- extra[[n]]
  }

  # Sample unmapped
  if (length(unmapped) > 0) {
    report$unmapped_sample <- head(unmapped, 20)
  }

  # Check for unexpected gaps
  if (match_rate < expected_rate - 0.5) {
    report$validation$unexpected_gaps <- list(list(
      issue = "match_rate_below_expected",
      expected = expected_rate,
      actual = round(match_rate, 2),
      gap = round(expected_rate - match_rate, 2)
    ))
  }

  return(report)
}

print_report <- function(report) {
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat(sprintf("%s CELL ID MAPPING REPORT\n", toupper(report$study)))
  cat(strrep("=", 70), "\n")

  s <- report$summary
  cat("\nSUMMARY\n")
  cat(sprintf("  Total component cells: %s\n", format(s$total_component_cells, big.mark = ",")))
  cat(sprintf("  Mapped:                %s\n", format(s$mapped, big.mark = ",")))
  cat(sprintf("  Unmapped:              %s\n", format(s$unmapped, big.mark = ",")))
  cat(sprintf("  Match rate:            %s%%\n", s$match_rate_pct))
  cat(sprintf("  Expected rate:         %s%%\n", s$expected_match_rate_pct))

  if (s$rate_within_tolerance) {
    cat("  Status:                ✓ WITHIN TOLERANCE\n")
  } else {
    cat("  Status:                ✗ BELOW EXPECTED\n")
  }

  if (length(report$unmapped_by_reason) > 0) {
    cat("\nUNMAPPED BY REASON\n")
    for (reason in names(report$unmapped_by_reason)) {
      cat(sprintf("  %s: %s\n", reason, format(report$unmapped_by_reason[[reason]], big.mark = ",")))
    }
  }

  if (!is.null(report$batch_to_patient_mapping)) {
    cat("\nBATCH -> PATIENT MAPPING\n")
    for (batch in sort(names(report$batch_to_patient_mapping))) {
      cat(sprintf("  %s -> %s\n", batch, report$batch_to_patient_mapping[[batch]]))
    }
  }

  v <- report$validation
  if (length(v$expected_issues) > 0) {
    cat("\nEXPECTED ISSUES\n")
    for (issue in v$expected_issues) {
      status <- if (issue %in% v$verified) "✓ verified" else "? pending"
      cat(sprintf("  %s: %s\n", issue, status))
    }
  }

  if (length(v$unexpected_gaps) > 0) {
    cat("\n⚠️  UNEXPECTED GAPS\n")
    for (gap in v$unexpected_gaps) {
      cat(sprintf("  %s: expected %s%%, got %s%%\n",
                  gap$issue, gap$expected, gap$actual))
    }
  }

  cat("\n")
  cat(strrep("=", 70), "\n")
}

# =============================================================================
# MAIN
# =============================================================================

main <- function() {
  parser <- ArgumentParser(description = "Generic cell ID mapping (R version)")
  parser$add_argument("--study", required = TRUE, help = "Study name from manifest")
  parser$add_argument("--manifest",
    default = "${SOURCE_IHBCAV1_HARMONIZATION%/}/config/study_manifest.yaml",
    help = "Path to study manifest")
  parser$add_argument("--dry-run", action = "store_true", help = "Print config only")

  args <- parser$parse_args()

  # Load manifest
  manifest <- load_manifest(args$manifest)
  shared <- manifest$shared

  if (!(args$study %in% names(manifest$studies))) {
    cat(sprintf("ERROR: Study '%s' not in manifest\n", args$study))
    cat(sprintf("Available: %s\n", paste(names(manifest$studies), collapse = ", ")))
    quit(status = 1)
  }

  config <- manifest$studies[[args$study]]

  cat(sprintf("=== Mapping %s ===\n", args$study))
  cat(sprintf("Source: %s\n", config$source_path))
  cat(sprintf("Format: %s\n", config$format))
  cat(sprintf("Mapping type: %s\n", config$mapping_type))
  cat(sprintf("Expected cells: %s\n", format(config$expected_cells, big.mark = ",")))
  cat(sprintf("Expected match rate: %s%%\n", config$expected_match_rate))
  cat("\n")

  if (args$dry_run) {
    cat("DRY RUN - exiting without mapping\n")
    return()
  }

  # Load Seurat for RDS files
  if (config$format == "rds") {
    suppressPackageStartupMessages(library(Seurat))
  }

  # Load iHBCA inventory
  ihbca_df <- load_ihbca_inventory(shared$ihbca_inventory, config$ihbca_dataset)

  # Load component study
  result <- load_component_study(config$source_path, config$format)
  component_ids <- result$cell_ids
  metadata <- result$metadata

  # Map based on type
  mapping_type <- config$mapping_type
  extra_info <- list()

  if (mapping_type == "identity") {
    map_result <- map_identity(component_ids, ihbca_df)
    mapped <- map_result$mapped
    unmapped <- map_result$unmapped

  } else if (mapping_type == "batch_aware") {
    regex <- config$cell_id_regex
    if (is.null(regex)) {
      cat("ERROR: batch_aware mapping requires cell_id_regex\n")
      quit(status = 1)
    }
    map_result <- map_batch_aware(component_ids, ihbca_df, regex)
    mapped <- map_result$mapped
    unmapped <- map_result$unmapped
    extra_info$batch_to_patient_mapping <- map_result$batch_mapping

  } else if (mapping_type == "custom_pal") {
    map_result <- map_custom_pal(component_ids, ihbca_df, config)
    mapped <- map_result$mapped
    unmapped <- map_result$unmapped
    extra_info$donor_mapping <- map_result$donor_mapping

  } else {
    cat(sprintf("ERROR: Mapping type '%s' not yet implemented\n", mapping_type))
    quit(status = 1)
  }

  # Get counts (handle both list and attribute-based)
  n_mapped <- if (!is.null(attr(mapped, "count"))) attr(mapped, "count") else length(mapped)
  n_unmapped <- if (!is.null(attr(unmapped, "total_count"))) attr(unmapped, "total_count") else length(unmapped)

  cat(sprintf("\nMapping results:\n"))
  cat(sprintf("  Mapped: %s\n", format(n_mapped, big.mark = ",")))
  cat(sprintf("  Unmapped: %s\n", format(n_unmapped, big.mark = ",")))

  # Build output dataframe
  cat("\nBuilding output dataframe...\n")

  # Use pre-built data.table if available (from vectorized identity mapping)
  if (!is.null(attr(mapped, "dt"))) {
    mapped_df <- attr(mapped, "dt")
    cat("  Using pre-built data.table (fast path)\n")
  } else {
    # Fallback for batch_aware and other methods
    mapped_df <- rbindlist(lapply(mapped, as.data.frame), fill = TRUE)
  }

  # Merge with metadata
  merged <- merge(mapped_df, metadata,
                  by.x = "component_cell_id", by.y = "_cell_id",
                  all.x = TRUE)

  # ---------------------------------------------------------------------------
  # Pal cell type mapping (from Austin's 10e_pal_data_preparation.R lines 89-104, 222-236)
  # Why use backticks? Because R thinks "0" is a string but list indexing wants it quoted.
  # ---------------------------------------------------------------------------
  PAL_CELLTYPE_MAPS <- list(
    pal_norm_epi = list(
      `0` = list(cell_type = "normEpi_LP", level0 = "Epithelial", level1 = "Luminal Progenitor"),
      `1` = list(cell_type = "normEpi_ML", level0 = "Epithelial", level1 = "Mature Luminal"),
      `2` = list(cell_type = "normEpi_Basal", level0 = "Epithelial", level1 = "Basal"),
      `3` = list(cell_type = "normEpi_Fb", level0 = "Stroma", level1 = "Fibroblast")
    ),
    pal_norm_total = list(
      `0` = list(cell_type = "normTotal_LP", level0 = "Epithelial", level1 = "Luminal Progenitor"),
      `1` = list(cell_type = "normTotal_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `2` = list(cell_type = "normTotal_Basal", level0 = "Epithelial", level1 = "Basal"),
      `3` = list(cell_type = "normTotal_ML", level0 = "Epithelial", level1 = "Mature Luminal"),
      `4` = list(cell_type = "normTotal_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `5` = list(cell_type = "normTotal_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `6` = list(cell_type = "normTotal_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `7` = list(cell_type = "normTotal_Stroma", level0 = "Stroma", level1 = "Stroma")
    ),
    pal_norm_b1 = list(
      `0` = list(cell_type = "normBr1_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `1` = list(cell_type = "normBr1_Epithelial", level0 = "Epithelial", level1 = "Epithelial"),
      `2` = list(cell_type = "normBr1_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `3` = list(cell_type = "normBr1_Epithelial", level0 = "Epithelial", level1 = "Epithelial"),
      `4` = list(cell_type = "normBr1_Epithelial", level0 = "Epithelial", level1 = "Epithelial"),
      `5` = list(cell_type = "normBr1_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `6` = list(cell_type = "normBr1_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `7` = list(cell_type = "normBr1_Epithelial", level0 = "Epithelial", level1 = "Epithelial"),
      `8` = list(cell_type = "normBr1_Stroma", level0 = "Stroma", level1 = "Stroma"),
      `9` = list(cell_type = "normBr1_Stroma", level0 = "Stroma", level1 = "Stroma")
    )
  )

  # Apply Pal cell type mapping if this is a Pal study
  if (args$study %in% names(PAL_CELLTYPE_MAPS)) {
    cat("Applying Austin's cell type mapping for Pal...\n")
    ctype_map <- PAL_CELLTYPE_MAPS[[args$study]]
    clusters <- as.character(merged$seurat_clusters)
    merged$cell_type <- sapply(clusters, function(c) {
      if (c %in% names(ctype_map)) ctype_map[[c]]$cell_type else NA_character_
    })
    merged$level0 <- sapply(clusters, function(c) {
      if (c %in% names(ctype_map)) ctype_map[[c]]$level0 else NA_character_
    })
    merged$level1 <- sapply(clusters, function(c) {
      if (c %in% names(ctype_map)) ctype_map[[c]]$level1 else NA_character_
    })
    cat(sprintf("  Added cell_type, level0, level1 for %d cells\n", nrow(merged)))
  }

  # Prefix columns with study name
  col_renames <- setdiff(names(merged), "ihbca_cell_id")
  new_names <- ifelse(
    grepl(paste0("^", args$study, "_"), col_renames),
    col_renames,
    paste0(args$study, "_", col_renames)
  )
  setnames(merged, col_renames, new_names)

  # Reorder: ihbca_cell_id first
  setcolorder(merged, c("ihbca_cell_id", setdiff(names(merged), "ihbca_cell_id")))

  # Write outputs
  output_dir <- file.path(shared$output_base, args$study, "outputs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  cells_path <- file.path(output_dir, paste0(args$study, "_cells.csv"))
  cat(sprintf("\nWriting cells to: %s\n", cells_path))
  fwrite(merged, cells_path)

  # Generate and write report
  report <- generate_report(
    args$study, config, mapped, unmapped,
    names(metadata), extra_info
  )

  findings_path <- file.path(output_dir, "mapping_findings.yaml")
  cat(sprintf("Writing findings to: %s\n", findings_path))
  write_yaml(report, findings_path)

  # Print report
  print_report(report)

  # Exit with error if unexpected gaps
  if (length(report$validation$unexpected_gaps) > 0) {
    cat("\n⚠️  FAILED: Unexpected gaps found. Investigation required.\n")
    quit(status = 1)
  }

  cat(sprintf("\n✓ %s mapping complete\n", args$study))
}

# Null coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a

# Run
main()
