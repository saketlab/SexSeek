#' Pull a counts matrix out of whatever the user handed us.
#'
#' Every input path converges here, so the scoring code only ever sees a matrix
#' with gene rownames. Optional containers are reached through
#' `requireNamespace()` so the package installs and works without them.
#' @noRd
.as_counts <- function(x, assay = NULL, ...) UseMethod(".as_counts")

#' @noRd
.as_counts.default <- function(x, assay = NULL, ...) {
  if (is.matrix(x) || methods::is(x, "Matrix")) {
    return(.check_matrix(x))
  }
  if (is.data.frame(x)) {
    return(.check_matrix(.df_to_matrix(x)))
  }
  stop("Cannot get a counts matrix from an object of class '",
    paste(class(x), collapse = "/"), "'.",
    call. = FALSE
  )
}

#' @noRd
.as_counts.character <- function(x, assay = NULL, ...) {
  if (length(x) != 1) {
    stop("Give a single file path.", call. = FALSE)
  }
  if (dir.exists(x)) {
    return(Read10xDir(x))
  }
  if (!file.exists(x)) {
    stop("File or directory not found: ", x, call. = FALSE)
  }
  ext <- tolower(sub(".*\\.", "", sub("\\.gz$", "", x)))
  obj <- switch(ext,
    rds = readRDS(x),
    csv = ,
    tsv = ,
    txt = .read_table_matrix(x, ext),
    h5 = Read10xH5(x),
    mtx = .mtx_sibling(x),
    h5ad = stop(
      "h5ad (AnnData) input is not supported. Convert to a 10x directory, or ",
      "read it in Python and pass the counts matrix.",
      call. = FALSE
    ),
    stop("Unsupported file extension '.", ext, "'.", call. = FALSE)
  )
  # an .rds can hold any supported object
  if (ext == "rds") .as_counts(obj, assay = assay, ...) else obj
}

#' Read a delimited counts table.
#'
#' First column is gene identifiers; everything else must be numeric. Reading
#' the header separately lets us set colClasses, which is the difference between
#' seconds and minutes on a wide matrix.
#' @noRd
.read_table_matrix <- function(path, ext) {
  sep <- if (ext == "csv") "," else .sniff_sep(path)
  con <- if (grepl("\\.gz$", path)) gzfile(path) else path
  df <- utils::read.table(con,
    sep = sep, header = TRUE, row.names = 1,
    check.names = FALSE, stringsAsFactors = FALSE, quote = "\""
  )
  .check_matrix(.df_to_matrix(df))
}

#' Guess the delimiter from the first line.
#'
#' A `.txt` deposit is as likely to be space-delimited as tab-delimited, and
#' reading a space-delimited matrix as TSV yields one giant character column
#' that then fails as "no numeric columns". `sep = ""` is read.table's
#' any-whitespace mode.
#' @noRd
.sniff_sep <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  line <- readLines(con, n = 1L, warn = FALSE)
  if (length(line) == 0) {
    return("\t")
  }
  if (grepl("\t", line)) "\t" else ""
}

#' @noRd
.df_to_matrix <- function(df) {
  num <- vapply(df, is.numeric, logical(1))
  if (!any(num)) {
    stop("No numeric columns found \u2014 is this a counts table?", call. = FALSE)
  }
  m <- as.matrix(df[, num, drop = FALSE])
  rownames(m) <- rownames(df)
  m
}

#' @noRd
.check_matrix <- function(m) {
  if (is.null(rownames(m))) {
    stop("The counts matrix has no rownames, so genes cannot be identified.",
      call. = FALSE
    )
  }
  if (nrow(m) == 0 || ncol(m) == 0) {
    stop("The counts matrix is empty.", call. = FALSE)
  }
  features <- rownames(m)
  if (anyNA(features) || any(!nzchar(trimws(features)))) {
    stop("Gene identifiers must not be missing or empty.", call. = FALSE)
  }
  # reject collisions under the same normalisation
  keys <- .norm_symbol(.strip_version(features))
  if (anyDuplicated(keys)) {
    stop("Duplicate or ambiguous gene identifiers. Supply unique gene IDs; ",
      "aggregate duplicate measurements explicitly only when justified.",
      call. = FALSE
    )
  }
  if (methods::is(m, "sparseMatrix") && !"x" %in% methods::slotNames(m)) {
    stop("Counts must be numeric; pattern matrices are not supported.", call. = FALSE)
  }
  values <- if (methods::is(m, "sparseMatrix")) m@x else as.vector(m)
  if (!is.numeric(values) || any(values < 0, na.rm = TRUE) ||
    any(is.infinite(values))) {
    stop("Counts must be numeric, non-negative and finite (NA is allowed).",
      call. = FALSE
    )
  }
  m
}

#' A bare .mtx path is only meaningful with its barcode and feature files, so
#' treat it as a pointer to the directory that holds all three.
#' @noRd
.mtx_sibling <- function(path) {
  prefix <- sub("matrix\\.mtx(\\.gz)?$", "", basename(path))
  if (identical(prefix, basename(path))) {
    stop("A 10x matrix path must end in matrix.mtx or matrix.mtx.gz.", call. = FALSE)
  }
  sibling <- function(suffixes) {
    candidates <- file.path(dirname(path), paste0(prefix, suffixes))
    found <- candidates[file.exists(candidates)]
    if (length(found) != 1L) {
      stop("Expected exactly one matching 10x sibling for ", basename(path),
        ": ", paste(suffixes, collapse = ", "),
        call. = FALSE
      )
    }
    found
  }
  .read_10x_triplet(
    path,
    sibling(c("barcodes.tsv", "barcodes.tsv.gz")),
    sibling(c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))
  )
}

#' @noRd
.need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required to ", what,
      ". Install it, or pass a plain counts matrix instead.",
      call. = FALSE
    )
  }
}
