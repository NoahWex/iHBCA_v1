# lib/cohort_filters.R
# Generic filter primitives for inquiry.yaml cohort definitions.
# A filter spec is a named list whose values are either:
#   - a vector of allowed values (positive filter)
#   - a list with `exclude` key giving values to drop
#   - the literal string "not_na" to require non-NA
# Returns a logical vector indicating rows to KEEP.

apply_filter_rule <- function(values, rule) {
  if (is.character(rule) && length(rule) == 1 && rule == "not_na") {
    return(!is.na(values))
  }
  if (is.list(rule) && !is.null(rule$exclude)) {
    return(!values %in% rule$exclude)
  }
  if (is.list(rule) && !is.null(rule$include)) {
    return(values %in% rule$include)
  }
  if (is.atomic(rule)) {
    return(values %in% rule)
  }
  stop("Unrecognized filter rule: ", paste(deparse(rule), collapse = " "))
}

apply_cohort_filters <- function(df, filters) {
  if (is.null(filters) || length(filters) == 0) return(rep(TRUE, nrow(df)))
  keep <- rep(TRUE, nrow(df))
  for (col in names(filters)) {
    if (!col %in% colnames(df)) {
      stop(sprintf("Filter column '%s' not in data frame columns: %s",
                   col, paste(colnames(df), collapse = ", ")))
    }
    keep <- keep & apply_filter_rule(df[[col]], filters[[col]])
  }
  keep
}

# Apply derived-column rules from inquiry.yaml `derived_columns` spec.
# Each derived column is a named entry with:
#   rules:    ordered list of {condition: {col: val, ...}, value: <new value>}
#             First-match-wins; condition matches when ALL listed (col, val)
#             pairs equal in the row.
#   default:  value for rows that match no rule (defaults to NA).
#   levels:   optional factor levels (first becomes reference).
# Returns df with added columns.
apply_derived_columns <- function(df, derived_spec) {
  if (is.null(derived_spec) || length(derived_spec) == 0) return(df)
  for (col_name in names(derived_spec)) {
    spec <- derived_spec[[col_name]]
    default <- spec$default %||% NA_character_
    new_col <- rep(default, nrow(df))
    if (!is.null(spec$rules)) {
      assigned <- rep(FALSE, nrow(df))
      for (rule in spec$rules) {
        cond <- rule$condition
        m <- rep(TRUE, nrow(df))
        for (cc in names(cond)) {
          if (!cc %in% colnames(df)) {
            stop(sprintf("Derived column '%s' rule references missing col: %s",
                         col_name, cc))
          }
          m <- m & !is.na(df[[cc]]) & as.character(df[[cc]]) == as.character(cond[[cc]])
        }
        new_col[m & !assigned] <- rule$value
        assigned <- assigned | m
      }
    }
    if (!is.null(spec$levels)) {
      df[[col_name]] <- factor(new_col, levels = spec$levels)
    } else {
      df[[col_name]] <- new_col
    }
  }
  df
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Build joint factor from inquiry.yaml `joint_factor` spec.
# Returns a factor of length nrow(df) with reference level set per spec.
build_joint_factor <- function(df, spec) {
  if (is.null(spec)) return(NULL)
  cols <- spec$from
  sep <- if (!is.null(spec$sep)) spec$sep else "."
  ref <- spec$reference
  jf <- apply(df[, cols, drop = FALSE], 1, function(r) paste(r, collapse = sep))
  if (!is.null(ref) && ref %in% jf) {
    levels_all <- unique(c(ref, sort(unique(jf))))
    factor(jf, levels = levels_all)
  } else {
    factor(jf)
  }
}
