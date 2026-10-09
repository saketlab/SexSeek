.example_data <- data.frame(
  dataset = c(
    "human_bulk", "mouse_bulk", "mouse_tx", "human_sc", "mouse_sc",
    "rat_bulk", "zebrafish_bulk"
  ),
  study = c(
    "GSE199208", "GSE85136", "GSE232600", "GSE196554", "GSE132420",
    "GSE100122", "GSE123439"
  ),
  species = c(
    "Homo sapiens", "Mus musculus", "Mus musculus", "Homo sapiens",
    "Mus musculus", "Rattus norvegicus", "Danio rerio"
  ),
  assay = c(
    "bulk", "bulk", "bulk, transcript level", "single-cell", "single-cell",
    "bulk", "bulk"
  ),
  tissue = c(
    "submandibular gland", "nucleus accumbens", "parathyroid gland",
    "skeletal muscle stem cells", "exorbital lacrimal gland",
    "dorsal root ganglion", "gonad"
  ),
  file = c(
    "human_bulk_GSE199208.tsv.gz", "mouse_bulk_GSE85136.tsv.gz",
    "mouse_tx_GSE232600.tsv.gz", "human_sc_GSE196554.rds",
    "mouse_sc_GSE132420.rds", "rat_bulk_GSE100122.tsv.gz",
    "zebrafish_bulk_GSE123439.tsv.gz"
  ),
  tx2gene = c(NA, NA, "mouse_tx_GSE232600_tx2gene.tsv", NA, NA, NA, NA),
  stringsAsFactors = FALSE
)

#' Example datasets with recorded sex
#'
#' Excerpts of public GEO datasets whose sample records state each donor's
#' sex. Each keeps the rows [EstimateSex()] scores for the species plus an
#' `OTHER` row with the rest of each library, so calls match the full matrix.
#' Excerpts are too small for species detection; pass `species`.
#'
#' @param dataset A name from `ExampleData()$dataset`. `NULL` lists them.
#' @return With `dataset = NULL`, a data frame of the available datasets.
#'   Otherwise a list with `dataset`, `study`, `species`, `assay`, `tissue`,
#'   `counts` (a matrix; a named list of sparse matrices, one per donor, for
#'   single-cell; transcript rows for `"mouse_tx"`), `recorded` (the recorded
#'   sex named by sample, as written in the record), `samples` (sample
#'   accession, title and recorded sex) and `tx2gene` (the transcript-to-gene
#'   table for `"mouse_tx"`, otherwise `NULL`).
#' @examples
#' ExampleData()
#' ex <- ExampleData("rat_bulk")
#' dim(ex$counts)
#' ex$recorded
#' @export
ExampleData <- function(dataset = NULL) {
  if (is.null(dataset)) {
    return(.example_data[, c("dataset", "study", "species", "assay", "tissue")])
  }
  i <- match(dataset, .example_data$dataset)
  if (length(dataset) != 1 || is.na(i)) {
    stop("Unknown dataset '", paste(dataset, collapse = ", "),
      "'. See ExampleData().",
      call. = FALSE
    )
  }
  row <- .example_data[i, ]
  path <- .example_file(row$file)
  counts <- if (grepl("\\.rds$", path)) readRDS(path) else .read_table_matrix(path, "tsv")

  if (is.null(.cache$example_samples)) {
    .cache$example_samples <- .read_tsv(.example_file("samples.tsv"))
  }
  samples <- .cache$example_samples
  samples <- samples[samples$dataset == dataset, c("sample", "title", "recorded_as"), drop = FALSE]
  units <- if (is.list(counts)) names(counts) else colnames(counts)
  hit <- match(units, samples$sample)
  if (anyNA(hit)) stop("Shipped sample table is missing units of '", dataset, "'.", call. = FALSE)
  samples <- samples[hit, , drop = FALSE]
  rownames(samples) <- NULL

  tx2gene <- if (is.na(row$tx2gene)) NULL else .read_tsv(.example_file(row$tx2gene))
  list(
    dataset = dataset, study = row$study, species = row$species,
    assay = row$assay, tissue = row$tissue, counts = counts,
    recorded = stats::setNames(samples$recorded_as, samples$sample),
    samples = samples, tx2gene = tx2gene
  )
}

#' @noRd
.example_file <- function(file) .extdata(file.path("examples", file))
