#' Read a 10x MatrixMarket triplet from a directory
#'
#' Expects `matrix.mtx`, a barcode file and a feature file, each optionally
#' gzipped, as CellRanger writes them. `Matrix::readMM` does the parsing, so
#' this costs no dependency beyond what the package already needs.
#'
#' @param dir Directory holding the triplet.
#' @param feature_column Which column of the feature file to use as rownames.
#'   CellRanger v3 writes id, symbol, type; column 1 is the Ensembl identifier,
#'   which is what species detection and panel matching want.
#' @return A sparse matrix with gene rownames and barcode colnames.
#' @export
Read10xDir <- function(dir, feature_column = 1L) {
  if (!dir.exists(dir)) {
    stop("Not a directory: ", dir, call. = FALSE)
  }
  # GEO deposits prefix the triplet with the accession and a library id --
  # `GSM9323718_6854-37_matrix.mtx.gz` -- so match on the suffix, not the whole
  # name. Bare CellRanger output still matches.
  mtx <- .pick_file(dir, "matrix\\.mtx(\\.gz)?$", "matrix.mtx")
  bar <- .pick_file(dir, "barcodes\\.tsv(\\.gz)?$", "barcodes.tsv")
  fea <- .pick_file(
    dir, "(features|genes)\\.tsv(\\.gz)?$", "features.tsv or genes.tsv"
  )

  .read_10x_triplet(mtx, bar, fea, feature_column)
}

#' @noRd
.read_10x_triplet <- function(mtx, bar, fea, feature_column = 1L) {
  m <- Matrix::readMM(mtx)
  features <- .read_col(fea, feature_column)
  barcodes <- .read_col(bar, 1L)

  # A dimension mismatch means the triplet is inconsistent — a partial download,
  # or files from different runs mixed in one directory. Silently recycling
  # names here would attach the wrong gene to every row.
  if (nrow(m) != length(features) || ncol(m) != length(barcodes)) {
    stop(
      "Triplet is inconsistent in ", dirname(mtx), ": matrix is ", nrow(m), "x", ncol(m),
      " but there are ", length(features), " features and ",
      length(barcodes), " barcodes.",
      call. = FALSE
    )
  }
  # CellRanger v2 wrote genes as rows; some deposits transpose. Trust the file
  # counts rather than guessing from shape.
  rownames(m) <- features
  colnames(m) <- barcodes
  .check_matrix(m)
}

#' @noRd
.pick_file <- function(dir, pattern, what) {
  hit <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(hit) == 0) {
    stop("No ", what, " in ", dir, call. = FALSE)
  }
  if (length(hit) > 1) {
    stop("More than one ", what, " in ", dir,
      " \u2014 split the samples into one directory each.",
      call. = FALSE
    )
  }
  hit[1]
}

#' @noRd
.read_col <- function(path, column) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  vapply(
    parts, function(p) if (length(p) >= column) p[column] else p[1],
    character(1)
  )
}

#' Read a 10x HDF5 counts file
#'
#' Handles both the CellRanger v2 layout (a single named group holding the
#' genome) and the v3 layout (a `matrix` group with a `features` subgroup).
#' Needs the `hdf5r` package, which is a suggested dependency.
#'
#' @param path Path to a `.h5` file.
#' @param feature_column For v3 files, `"id"` for Ensembl identifiers or
#'   `"name"` for symbols. Identifiers are preferred — symbols cannot
#'   distinguish species.
#' @return A sparse matrix with gene rownames and barcode colnames.
#' @export
Read10xH5 <- function(path, feature_column = c("id", "name")) {
  .need("hdf5r", "read a 10x HDF5 file")
  feature_column <- match.arg(feature_column)

  h5 <- hdf5r::H5File$new(path, mode = "r")
  on.exit(h5$close_all())

  root <- if ("matrix" %in% names(h5)) "matrix" else names(h5)[1]
  g <- h5[[root]]

  # v3 nests feature metadata under `features`; v2 keeps `genes`/`gene_names`
  # beside the data. Both carry the same CSC arrays.
  ids <- if ("features" %in% names(g)) {
    key <- if (feature_column == "id") "id" else "name"
    g[["features"]][[key]]$read()
  } else {
    key <- if (feature_column == "id") "genes" else "gene_names"
    g[[key]]$read()
  }

  m <- Matrix::sparseMatrix(
    i = g[["indices"]]$read() + 1L,
    p = g[["indptr"]]$read(),
    x = as.numeric(g[["data"]]$read()),
    dims = g[["shape"]]$read(),
    repr = "C"
  )
  rownames(m) <- ids
  colnames(m) <- g[["barcodes"]]$read()
  .check_matrix(m)
}
