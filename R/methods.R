#' @rdname EstimateSex
#' @param assay Assay to read counts from. Defaults to the object's active or
#'   first assay.
#' @export
EstimateSex.Seurat <- function(x, species = NULL, group = NULL,
                               per_cell = FALSE, assay = NULL,
                               model = c("ratio", "logistic"),
                               annotation = NULL,
                               input_scale = c("auto", "counts", "normalised"), ...) {
  .need("SeuratObject", "read a Seurat object")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(x)

  a <- x[[assay]]
  if (inherits(a, "Assay5")) {
    layers <- SeuratObject::Layers(a, search = "^counts($|\\.)")
    if (!length(layers)) stop("Assay has no counts layers.", call. = FALSE)
    matrices <- lapply(layers, function(layer) {
      .check_matrix(SeuratObject::LayerData(a, layer = layer))
    })
    genes <- rownames(matrices[[1]])
    if (!all(vapply(matrices, function(m) setequal(rownames(m), genes), logical(1)))) {
      stop("Counts layers have different gene sets. Supply a common, unfiltered ",
        "gene set; missing genes cannot be filled with zero safely.",
        call. = FALSE
      )
    }
    cells <- unlist(lapply(matrices, colnames), use.names = FALSE)
    if (anyDuplicated(cells)) {
      stop("Counts layers contain overlapping cells; select an unambiguous assay.",
        call. = FALSE
      )
    }
    m <- do.call(cbind, lapply(matrices, function(m) m[genes, , drop = FALSE]))
  } else {
    m <- SeuratObject::GetAssayData(x, assay = assay, slot = "counts")
  }
  if (prod(dim(m)) == 0) {
    stop("Assay '", assay, "' holds no counts. Sex estimation needs raw counts,",
      " not scaled or normalised data.",
      call. = FALSE
    )
  }
  .estimate_on_matrix(
    .check_matrix(m), species,
    .group_from(group, x[[]], colnames(m)), per_cell,
    model = match.arg(model), annotation = annotation,
    input_scale = match.arg(input_scale)
  )
}

#' @rdname EstimateSex
#' @export
EstimateSex.SingleCellExperiment <- function(x, species = NULL, group = NULL,
                                             per_cell = FALSE, assay = "counts",
                                             model = c("ratio", "logistic"),
                                             annotation = NULL,
                                             input_scale = c("auto", "counts", "normalised"), ...) {
  .need("SummarizedExperiment", "read a SingleCellExperiment")
  m <- SummarizedExperiment::assay(x, assay)
  meta <- as.data.frame(SummarizedExperiment::colData(x))
  .estimate_on_matrix(
    .check_matrix(m), species,
    .group_from(group, meta, colnames(m)), per_cell,
    model = match.arg(model), annotation = annotation,
    input_scale = match.arg(input_scale)
  )
}

#' @rdname EstimateSex
#' @export
EstimateSex.SummarizedExperiment <- EstimateSex.SingleCellExperiment

#' Resolve `group` against object metadata.
#'
#' A length-one `group` naming a metadata column is the common case; users say
#' `group = "donor"`, not `group = obj$donor`. A longer vector is taken as the
#' grouping itself.
#' @noRd
.group_from <- function(group, meta, cells) {
  if (is.null(group)) {
    return(NULL)
  }
  if (length(group) == 1 && is.character(group) && !is.null(meta) &&
    group %in% colnames(meta)) {
    if (is.null(cells) || anyNA(cells) || anyDuplicated(cells) ||
      !all(cells %in% rownames(meta))) {
      stop("Counts columns cannot be aligned uniquely to metadata rows.", call. = FALSE)
    }
    return(as.character(meta[[group]][match(cells, rownames(meta))]))
  }
  if (length(group) == 1 && length(cells) != 1) {
    stop("'", group, "' is not a column in the object's metadata.",
      call. = FALSE
    )
  }
  group
}
