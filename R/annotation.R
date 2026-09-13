#' Retrieve sex-chromosome annotations with genevintage
#'
#' Returns a chromosome inventory, not a diagnostic marker panel. Pass the
#' result to `EstimateSex(annotation = ...)` to reconcile curated markers
#' against this annotation. Fetch once and reuse across datasets.
#'
#' @param species Scientific or common species name.
#' @param release,assembly,source Annotation selection passed to
#'   [genevintage::sex_genes()]. Pin these for reproducible analyses.
#' @param mapping Optional local annotation with `id`, `name`, `chr`, `biotype`
#'   columns. Supplying this avoids downloads.
#' @return A data frame with `id`, `name`, `chr`, and `biotype` columns.
#' @export
SexChromosomeGenes <- function(species, release = NULL, assembly = NULL,
                               source = NULL, mapping = NULL) {
  .need("genevintage", "retrieve sex-chromosome annotations")
  out <- genevintage::sex_genes(
    species = species, release = release, assembly = assembly,
    source = source, mapping = mapping, quiet = TRUE
  )
  attr(out, "SexSeek_species") <- ResolveSpecies(species) %||% species
  out
}

#' @noRd
.prepare_annotation <- function(annotation, species) {
  if (is.null(annotation)) {
    return(NULL)
  }
  if (identical(annotation, "genevintage")) {
    annotation <- SexChromosomeGenes(species)
  }
  if (!is.data.frame(annotation) ||
    !all(c("id", "name", "chr", "biotype") %in% names(annotation))) {
    stop("`annotation` must be NULL, 'genevintage', or a data frame with ",
      "id, name, chr and biotype columns.",
      call. = FALSE
    )
  }
  sp <- attr(annotation, "SexSeek_species")
  if (!is.null(sp) && !identical(sp, species)) {
    stop("Annotation species does not match the input species.", call. = FALSE)
  }
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  annotation[] <- lapply(annotation, as.character)
  annotation$chr <- sub("^chr", "", annotation$chr)
  annotation$name[is.na(annotation$name)] <- ""
  if (anyNA(annotation[c("id", "chr")]) || any(!nzchar(annotation$id))) {
    stop("Annotation identifiers and chromosomes must not be missing.",
      call. = FALSE
    )
  }
  annotation
}

#' Reconcile curated markers on their expected chromosome. Ambiguous symbols
#' are excluded instead of selecting an arbitrary gene.
#' @noRd
.annotate_markers <- function(panel, annotation, chromosome) {
  if (is.null(annotation)) {
    return(panel)
  }
  hit <- vapply(seq_len(nrow(panel)), function(i) {
    eligible <- which(annotation$chr == chromosome[i])
    by_id <- eligible[.strip_version(annotation$id[eligible]) ==
      .strip_version(panel$gene_id[i])]
    if (length(by_id) == 1L) {
      return(by_id)
    }
    by_name <- if (nzchar(panel$gene_name[i])) {
      eligible[.norm_symbol(annotation$name[eligible]) ==
        .norm_symbol(panel$gene_name[i])]
    } else {
      integer()
    }
    if (length(by_name) == 1L) by_name else NA_integer_
  }, integer(1))
  panel$gene_id <- annotation$id[hit]
  panel$gene_name <- annotation$name[hit]
  panel
}
