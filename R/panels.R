#' Package data cache. SexSeek reads each file once per session.
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

#' Read a TSV
#'
#' Keep `quote = ""`. Some gene symbols contain apostrophes. With quoting on,
#' R drops those rows without a warning.
#' @noRd
.read_tsv <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  utils::read.delim(con,
    stringsAsFactors = FALSE, quote = "", comment.char = "",
    colClasses = "character", na.strings = character(0)
  )
}

#' Sex-linked gene panels for all supported species
#'
#' One row for each sex-linked gene. `role` is `Y`, `X`, `W`, `Z`,
#' `inactivation` (XIST, RSX) or `male_specific` (roX1, roX2).
#'
#' @return A data frame with one row for each gene. The columns include
#'   `scientific_name`, `gene_id`, `gene_name`, `chr`, `role` and `tier`.
#' @examples
#' panels <- SexPanels()
#' table(panels$role)
#' @export
SexPanels <- function() {
  if (is.null(.cache$panels)) {
    .cache$panels <- .read_tsv(.extdata("sex_panels.tsv.gz"))
  }
  .cache$panels
}

#' Supported species and their sex-determination systems
#'
#' @return A data frame with one row for each species. The columns include the
#'   sex-determination `system`, the `heterogametic` sex and the panel `status`.
#' @examples
#' head(SexSpecies()[, c("scientific_name", "common_name", "system", "status")])
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

#' Panel rows for one species
#'
#' The default `core` tier holds the curated somatic markers that the score uses.
#' Ampliconic and testis Y genes in the full inventory follow tissue type.
#'
#' @param species A scientific or common name that
#'   [ResolveSpecies()] accepts.
#' @param tier The tiers to keep. `NULL` keeps all rows.
#' @return A data frame of panel rows.
#' @examples
#' PanelFor("Homo sapiens")[, c("gene_name", "role", "tier")]
#' PanelFor("Gallus gallus")[, c("gene_name", "role", "tier")]
#' @export
PanelFor <- function(species, tier = "core") {
  p <- SexPanels()
  out <- p[p$scientific_name == (ResolveSpecies(species) %||% species), , drop = FALSE]
  if (nrow(out) == 0) {
    stop("No panel for species '", species, "'.", call. = FALSE)
  }
  if (!is.null(tier)) {
    out <- out[out$tier %in% tier, , drop = FALSE]
  }
  # Drop `exclude` rows even when `tier` asks for them.
  out[out$tier != "exclude", , drop = FALSE]
}

#' Registry row for one species
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
#' Each row pairs a Y gene with its X homolog. Scoring uses complete pairs.
#' Y/(Y+X) cancels library size. Gene effects remain, so calibrate it on real data.
#'
#' @return A data frame with one row for each pair: `scientific_name`, the Y and
#'   X gene IDs and symbols, and the panel tier of the Y gene.
#' @examples
#' gp <- GametologPairs()
#' gp[gp$scientific_name == "Homo sapiens", c("y_gene_name", "x_gene_name")]
#' @export
GametologPairs <- function() {
  if (is.null(.cache$pairs)) {
    .cache$pairs <- .read_tsv(.extdata("gametolog_pairs.tsv"))
  }
  .cache$pairs
}
