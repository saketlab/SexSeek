#' Cached package data, loaded once per session.
#'
#' Reading the panel is the one unavoidable I/O cost, so it happens at most once
#' and lives here rather than being passed through every call.
#' @noRd
.cache <- new.env(parent = emptyenv())

#' @noRd
.extdata <- function(file) {
  path <- system.file("extdata", file, package = "SexSeek")
  if (!nzchar(path)) {
    stop("Package data '", file, "' is missing. Reinstall SexSeek.",
      call. = FALSE
    )
  }
  path
}

#' Read a shipped TSV with base R only.
#'
#' `read.delim` with quoting off: gene symbols legitimately contain apostrophes
#' and quotes would silently swallow rows.
#' @noRd
.read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  utils::read.delim(con,
    stringsAsFactors = FALSE, quote = "", comment.char = "",
    colClasses = "character", na.strings = character(0)
  )
}

#' Sex-linked gene panels for every supported species
#'
#' One row per sex-linked gene. `role` is one of `Y`, `X`, `W`, `Z`,
#' `inactivation` (XIST, RSX) or `male_specific` (roX1, roX2).
#'
#' @return A data frame.
#' @export
SexPanels <- function() {
  if (is.null(.cache$panels)) {
    .cache$panels <- .read_tsv(.extdata("sex_panels.tsv.gz"))
  }
  .cache$panels
}

#' Supported species and their sex-determination systems
#'
#' @return A data frame with one row per species.
#' @export
SexSpecies <- function() {
  if (is.null(.cache$species)) {
    .cache$species <- .read_tsv(.extdata("species.tsv"))
  }
  .cache$species
}

#' @noRd
.aliases <- function() {
  if (is.null(.cache$aliases)) {
    .cache$aliases <- .read_tsv(.extdata("species_aliases.tsv"))
  }
  .cache$aliases
}

#' The panel rows for one species
#'
#' Defaults to the `core` tier, which is what scoring uses: a curated set of
#' uniquely mappable, somatic-tissue markers. The full inventory carries every
#' gene on the sex chromosomes and is far too permissive to call from —
#' ampliconic and testis-restricted Y families in particular fire on tissue, not
#' on genotype.
#'
#' @param species A scientific name, as it appears in [SexSpecies()].
#' @param tier Which tiers to keep. `NULL` returns every row.
#' @return A data frame of panel rows.
#' @export
PanelFor <- function(species, tier = "core") {
  p <- SexPanels()
  out <- p[p$scientific_name == species, , drop = FALSE]
  if (nrow(out) == 0) {
    stop("No panel for species '", species, "'.", call. = FALSE)
  }
  if (!is.null(tier)) {
    out <- out[out$tier %in% tier, , drop = FALSE]
  }
  # `exclude` never contributes, whatever the caller asked for.
  out[out$tier != "exclude", , drop = FALSE]
}

#' The registry row for one species
#' @noRd
.species_row <- function(species) {
  s <- SexSpecies()
  row <- s[s$scientific_name == species, , drop = FALSE]
  if (nrow(row) == 0) {
    return(NULL)
  }
  row[1, , drop = FALSE]
}

#' X/Y gametolog pairs
#'
#' Each row pairs a Y gene with its X-linked homolog. Scoring uses complete
#' pairs from the selected marker tier. A common library-size multiplier
#' cancels in Y/(Y+X), but gene-specific capture, length, mapping and tissue
#' effects remain. The resulting fraction requires empirical calibration.
#'
#' @return A data frame of pairs.
#' @export
GametologPairs <- function() {
  if (is.null(.cache$pairs)) {
    .cache$pairs <- .read_tsv(.extdata("gametolog_pairs.tsv"))
  }
  .cache$pairs
}
