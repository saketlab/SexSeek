#' @rdname EstimateSex
#' @details A named list of matrices is scored one element per unit; elements
#'   may have different gene lists and are never zero-padded. `group` is not
#'   accepted for a list. With `per_cell = TRUE` each unit is
#'   `<element>:<column>`.
#' @export
EstimateSex.list <- function(x, ..., group = NULL, per_cell = FALSE) {
  .check_unit_list(x, group)
  .per_unit(x, EstimateSex, ..., per_cell = per_cell, prefix = isTRUE(per_cell))
}

#' f on each element, unit set to the element name (or prefixed by it)
#' @noRd
.per_unit <- function(x, f, ..., prefix = FALSE) {
  out <- do.call(rbind, lapply(names(x), function(u) {
    r <- f(x[[u]], ...)
    r$unit <- if (prefix) paste0(u, ":", r$unit) else u
    r
  }))
  rownames(out) <- NULL
  out
}

#' @noRd
.check_unit_list <- function(x, group = NULL) {
  nm <- names(x)
  if (!length(x) || is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) {
    stop("A list needs unique, non-empty names, one per unit.", call. = FALSE)
  }
  if (!is.null(group)) {
    stop("Each list element is one unit; do not pass 'group'.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Compare calls with recorded sex
#'
#' Joins [EstimateSex()] calls to the sex recorded for each unit. Only
#' `female`/`male` (or `F`/`M`, any case) are scored; other labels are kept as
#' recorded.
#'
#' @param calls Output of [EstimateSex()].
#' @param recorded Recorded sex as a character vector named by unit, matching
#'   `calls$unit`.
#' @return `calls` with three columns added: `recorded` (as given),
#'   `recorded_sex` (`"female"`, `"male"` or `NA`) and `agrees`: `TRUE` when the
#'   verdict is the recorded sex, `FALSE` when it is the other sex, and `NA`
#'   when the call abstains (`uncertain`, `possible_mixed`, `unknown`) or no
#'   usable sex is recorded.
#' @examples
#' ex <- ExampleData("mouse_bulk")
#' calls <- EstimateSex(ex$counts, species = ex$species, group = colnames(ex$counts))
#' cmp <- CompareSex(calls, ex$recorded)
#' table(recorded = cmp$recorded_sex, call = cmp$verdict)
#' @export
CompareSex <- function(calls, recorded) {
  if (!is.data.frame(calls) || !all(c("unit", "verdict") %in% names(calls))) {
    stop("'calls' must be the output of EstimateSex().", call. = FALSE)
  }
  calls$recorded <- .recorded_for(calls$unit, recorded)
  calls$recorded_sex <- .recorded_sex(calls$recorded)
  called <- calls$verdict %in% c("female", "male")
  calls$agrees <- ifelse(called & !is.na(calls$recorded_sex),
    calls$verdict == calls$recorded_sex, NA
  )
  calls
}

#' @noRd
.recorded_for <- function(units, recorded) {
  if (is.null(names(recorded)) || anyNA(names(recorded)) || anyDuplicated(names(recorded))) {
    stop("'recorded' must be named by unit, with unique names.", call. = FALSE)
  }
  unname(as.character(recorded[match(units, names(recorded))]))
}

#' anything but female/male is unscorable
#' @noRd
.recorded_sex <- function(x) {
  x <- tolower(trimws(x))
  out <- rep(NA_character_, length(x))
  out[x %in% c("f", "female")] <- "female"
  out[x %in% c("m", "male")] <- "male"
  out
}

#' Count agreements with recorded sex
#'
#' @param x Output of [CompareSex()], or a named list of them.
#' @return One row per comparison: `samples`, and how many calls `agree`,
#'   `disagree` or `abstain` (including samples with no usable recorded sex).
#' @examples
#' ex <- ExampleData("rat_bulk")
#' calls <- EstimateSex(ex$counts, species = ex$species, group = colnames(ex$counts))
#' SummariseComparison(CompareSex(calls, ex$recorded))
#' @export
SummariseComparison <- function(x) {
  one <- function(d) {
    if (!is.data.frame(d) || !"agrees" %in% names(d)) {
      stop("Give the output of CompareSex(), or a named list of them.", call. = FALSE)
    }
    data.frame(
      samples = nrow(d), agree = sum(d$agrees %in% TRUE),
      disagree = sum(d$agrees %in% FALSE), abstain = sum(is.na(d$agrees))
    )
  }
  if (is.data.frame(x)) {
    return(one(x))
  }
  .check_unit_list(x)
  out <- do.call(rbind, lapply(x, one))
  cbind(comparison = names(x), out, row.names = NULL)
}

#' Fraction of cells carrying any sex-specific marker count
#'
#' For each unit, the share of columns (cells) with at least one count across
#' the species' scoring-tier sex-specific markers (Y in XY species, W in ZW
#' species).
#'
#' @param x Counts with cells in columns (matrix, sparse matrix, data frame or
#'   file path), or a named list of them, one per unit.
#' @param species Species name, scientific or common.
#' @param group For a matrix, the unit of each column. `NULL` treats the whole
#'   matrix as one unit.
#' @return A data frame with `unit`, `n_cells` and `fraction`.
#' @examples
#' ex <- ExampleData("mouse_sc")
#' MarkerCellFraction(ex$counts, ex$species)
#' @export
MarkerCellFraction <- function(x, species, group = NULL) {
  UseMethod("MarkerCellFraction")
}

#' @rdname MarkerCellFraction
#' @export
MarkerCellFraction.list <- function(x, species, group = NULL) {
  .check_unit_list(x, group)
  .per_unit(x, MarkerCellFraction, species = species)
}

#' @rdname MarkerCellFraction
#' @export
MarkerCellFraction.default <- function(x, species, group = NULL) {
  x <- .as_counts(x)
  s <- .scoring_setup(rownames(x), species)
  spec <- if (is.null(s$reason)) .specific_panel(s$panel, s$reg$system)
  if (!NROW(spec)) {
    stop("No sex-specific markers in the scoring panel for '", s$species %||% s$reg$scientific_name,
      "'.", if (!is.null(s$note)) paste0(" ", s$note),
      call. = FALSE
    )
  }
  hit <- Matrix::colSums(x[.match_rows(rownames(x), spec), , drop = FALSE]) > 0
  idx <- .unit_index(x, group, per_cell = FALSE)
  data.frame(
    unit = names(idx), n_cells = lengths(idx),
    fraction = vapply(idx, function(i) mean(hit[i]), 1), row.names = NULL
  )
}
