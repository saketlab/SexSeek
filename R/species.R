#' Strip an Ensembl version suffix.
#'
#' ENSG00000229807.13 and ENSG00000229807 are the same gene. Matrices can use
#' either form.
#' @noRd
.strip_version <- function(x) {
  # Only Ensembl IDs. Clone symbols such as AL627309.1 are different genes.
  sub("^(ENS[A-Z]*[GTP][0-9]+)\\.[0-9]+$", "\\1", x, ignore.case = TRUE)
}

#' @noRd
.norm_symbol <- function(x) toupper(trimws(x))

#' Change a species name to its scientific name
#'
#' Give a scientific name, a common name, or an alias from the package alias
#' table. Matching ignores case.
#'
#' @param name A species name.
#' @return The scientific name, or `NULL` if the name is not known.
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

#' Detect species from gene IDs or symbols
#'
#' The function tries stable ID prefixes first. If no prefix matches, it compares
#' symbols with each species panel. Close species share symbols, so check `confidence`.
#'
#' @param features A character vector of gene IDs or symbols, such as
#'   `rownames(counts)`.
#' @return A list with `species`, `confidence` (`high`, `medium`, `low`) and
#'   `method` (`id_prefix`, `symbol`, `ambiguous` or `none`). `species` is
#'   `NULL` when no species matches. A `low` symbol match is a near tie between
#'   species with the same sex system. `candidates` gives these species.
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
  # "NA" and short prefixes match symbols like NACA and call honey bee.
  s <- s[!is.na(s$id_prefix) & nzchar(s$id_prefix) & s$id_prefix != "NA" &
    nchar(s$id_prefix) >= 4, , drop = FALSE]

  # Longest prefix first, so that ENSMUSG wins over a shorter prefix in it.
  s <- s[order(-nchar(s$id_prefix)), , drop = FALSE]
  for (i in seq_len(nrow(s))) {
    n <- sum(startsWith(ids, s$id_prefix[i]))
    # Without the fraction limit, human AGAP symbols call mosquito.
    if (n >= 0.2 * length(ids) && n >= min(10L, length(ids))) {
      return(list(
        species = s$scientific_name[i], confidence = "high",
        method = "id_prefix"
      ))
    }
  }

  .detect_by_symbol(features)
}

#' Symbol overlap fallback.
#'
#' The score is the fraction of matched input symbols in each species panel.
#' Conserved X symbols match every mammal panel, so coverage cannot rank species.
#' A near tie across sex systems gives no species.
#' @noRd
.detect_by_symbol <- function(features) {
  none <- function(method) list(species = NULL, confidence = "low", method = method)
  syms <- unique(.norm_symbol(features))
  p <- SexPanels()
  p <- p[nzchar(p$gene_name), , drop = FALSE]
  p$sym <- .norm_symbol(p$gene_name)

  by_sp <- lapply(split(p$sym, p$scientific_name), unique)
  # Small panels can get a high score by chance.
  by_sp <- by_sp[lengths(by_sp) >= 50L]
  if (!length(by_sp)) {
    return(none("none"))
  }

  seen <- intersect(syms, unlist(by_sp, use.names = FALSE))
  # Most vertebrates have ZFX, KDM5C and similar genes.
  if (length(seen) < 10L) {
    return(none("none"))
  }

  score <- sort(vapply(by_sp, function(x) mean(seen %in% x), numeric(1)),
    decreasing = TRUE
  )
  best <- names(score)[1]
  margin <- score[[1]] - if (length(score) > 1) score[[2]] else 0
  if (margin >= 0.1) {
    return(list(species = best, confidence = "medium", method = "symbol"))
  }

  close <- names(score)[score >= score[[1]] - 0.1]
  s <- SexSpecies()
  sex_of <- function(sp) paste(s$system, s$heterogametic)[match(sp, s$scientific_name)]
  sex <- sex_of(close)
  if (length(unique(sex)) > 1L) {
    # Z genes are autosomal in mammals. Y/W names appear whatever the sample sex.
    scored <- p[p$tier %in% c("core", "secondary"), , drop = FALSE]
    n_spec <- vapply(close, function(sp) {
      spec <- .specific_panel(
        scored[scored$scientific_name == sp, , drop = FALSE],
        s$system[match(sp, s$scientific_name)]
      )
      sum(unique(spec$sym) %in% syms)
    }, integer(1))
    close <- close[n_spec >= 2L]
    sex <- sex_of(close)
  }
  if (!length(close) || anyNA(sex) || length(unique(sex)) > 1L) {
    return(none("ambiguous"))
  }
  list(
    species = close[1], confidence = "low", method = "symbol",
    candidates = close
  )
}
