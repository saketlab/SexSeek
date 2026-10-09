#' Strip an Ensembl version suffix.
#'
#' ENSG00000229807.13 and ENSG00000229807 are the same gene; matrices disagree
#' about which form they carry.
#' @noRd
.strip_version <- function(x) {
  # only Ensembl IDs; dotted clone symbols like AL627309.1 are distinct genes
  sub("^(ENS[A-Z]*[GTP][0-9]+)\\.[0-9]+$", "\\1", x, ignore.case = TRUE)
}

#' @noRd
.norm_symbol <- function(x) toupper(trimws(x))

#' Resolve a user-supplied species name to a scientific name
#'
#' Accepts a scientific name, a common name, or any alias in the shipped alias
#' table, in any case.
#'
#' @param name A species name.
#' @return The scientific name, or `NULL` if unrecognised.
#' @examples
#' ResolveSpecies("mouse")
#' ResolveSpecies("CHICKEN")
#' ResolveSpecies("Danio rerio")
#' @export
ResolveSpecies <- function(name) {
  if (is.null(name) || !nzchar(name)) {
    return(NULL)
  }
  key <- tolower(trimws(name))

  s <- SexSpecies()
  hit <- match(key, tolower(s$scientific_name))
  if (!is.na(hit)) {
    return(s$scientific_name[hit])
  }

  a <- .aliases()
  hit <- match(key, tolower(a$alias))
  if (!is.na(hit)) {
    return(a$scientific_name[hit])
  }
  NULL
}

#' Detect species from gene identifiers or symbols
#'
#' Tries stable-ID prefixes first; they are unambiguous and cost one regex per
#' species. Falls back to scoring the input symbols against each species' panel
#' symbols, which is weaker: panels overlap heavily between close relatives, so
#' the returned confidence is what tells you whether to trust it.
#'
#' @param features A character vector of gene identifiers or symbols, typically
#'   `rownames(counts)`.
#' @return A list with `species`, `confidence` (`high`, `medium`, `low`) and
#'   `method` (`id_prefix`, `symbol`, or `none`). `species` is `NULL` when
#'   nothing matched.
#' @examples
#' DetectSpecies(c("ENSMUSG00000086503", "ENSMUSG00000069045"))
#' DetectSpecies(c("ENSG00000229807.13", "ENSG00000012817"))
#' @export
DetectSpecies <- function(features) {
  features <- features[!is.na(features) & nzchar(features)]
  if (length(features) == 0) {
    return(list(species = NULL, confidence = "low", method = "none"))
  }

  ids <- .strip_version(features)
  s <- SexSpecies()
  # literal "NA" and short prefixes match symbols like NACA and call honey bee
  s <- s[!is.na(s$id_prefix) & nzchar(s$id_prefix) & s$id_prefix != "NA" &
    nchar(s$id_prefix) >= 4, , drop = FALSE]

  # longest first so ENSMUSG beats any shorter prefix it contains
  s <- s[order(-nchar(s$id_prefix)), , drop = FALSE]
  for (i in seq_len(nrow(s))) {
    n <- sum(startsWith(ids, s$id_prefix[i]))
    # a count alone lets human symbols like AGAP call mosquito
    if (n >= 0.2 * length(ids) && n >= min(10L, length(ids))) {
      return(list(
        species = s$scientific_name[i], confidence = "high",
        method = "id_prefix"
      ))
    }
  }

  .detect_by_symbol(features)
}

#' Symbol-overlap fallback.
#'
#' Scores against panel symbols only. Close relatives share most symbols, so a
#' win by a thin margin is reported as low confidence rather than resolved.
#' @noRd
.detect_by_symbol <- function(features) {
  syms <- unique(.norm_symbol(features))
  p <- SexPanels()
  p <- p[nzchar(p$gene_name), , drop = FALSE]
  if (nrow(p) == 0) {
    return(list(species = NULL, confidence = "low", method = "none"))
  }

  by_sp <- split(.norm_symbol(p$gene_name), p$scientific_name)
  sizes <- vapply(by_sp, function(x) length(unique(x)), integer(1))

  # tiny panels score high by coincidence
  MIN_PANEL <- 50L
  usable <- sizes >= MIN_PANEL
  if (!any(usable)) {
    return(list(species = NULL, confidence = "low", method = "none"))
  }
  by_sp <- by_sp[usable]
  sizes <- sizes[usable]

  hits <- vapply(by_sp, function(x) sum(unique(x) %in% syms), integer(1))
  # ZFX, KDM5C and similar are shared across most vertebrates
  if (max(hits) < 10L) {
    return(list(species = NULL, confidence = "low", method = "none"))
  }

  frac <- hits / sizes
  ord <- order(-frac)
  best <- names(frac)[ord[1]]
  margin <- if (length(frac) > 1) frac[ord[1]] - frac[ord[2]] else frac[ord[1]]

  # mammal X symbols are near-identical across species; refuse
  if (!(frac[ord[1]] >= 0.5 && margin >= 0.15)) {
    return(list(species = NULL, confidence = "low", method = "ambiguous"))
  }
  conf <- "medium"

  # a close call across XY and ZW inverts every result
  if (length(ord) > 1 && margin < 0.15) {
    sys <- .system_of(names(frac)[ord[1:2]])
    if (!anyNA(sys) && sys[1] != sys[2]) {
      return(list(species = NULL, confidence = "low", method = "ambiguous"))
    }
  }
  list(species = best, confidence = conf, method = "symbol")
}

#' @noRd
.system_of <- function(species) {
  s <- SexSpecies()
  s$system[match(species, s$scientific_name)]
}
