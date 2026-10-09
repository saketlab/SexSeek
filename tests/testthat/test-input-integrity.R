test_that("duplicate matching keys are rejected regardless of representation or order", {
  m <- paired_counts()
  for (names in list(
    c("ENSG00000129824.1", "ENSG00000129824.2"),
    c("RPS4Y1", "rps4y1"), c("DDX3Y", "DDX3Y")
  )) {
    z <- rbind(m, m[1:2, , drop = FALSE])
    rownames(z)[(nrow(z) - 1):nrow(z)] <- names
    for (x in list(z, z[nrow(z):1, , drop = FALSE], Matrix::Matrix(z, sparse = TRUE))) {
      expect_error(EstimateSex(x, species = "human"), "Duplicate or ambiguous")
    }
  }
})

test_that("all values are checked for fractional input", {
  m <- matrix(c(rep(1, 2000), 0.5),
    ncol = 1,
    dimnames = list(paste0("G", 1:2001), "sample")
  )
  for (x in list(m, Matrix::Matrix(m, sparse = TRUE))) {
    expect_true(SexSeek:::.looks_normalised(x, chunk_size = 100L))
    expect_error(
      EstimateSex(x, species = "human", input_scale = "counts"),
      "fractional values"
    )
  }
  expect_false(SexSeek:::.looks_normalised(m * 2, chunk_size = 100L))
  for (scale in c("auto", "counts", "normalised")) {
    r <- EstimateSex(paired_counts(), species = "human", input_scale = scale)
    expect_equal(r$verdict, "male")
    expect_equal(grepl("normalised_input", r$flags), scale == "normalised")
    if (scale == "normalised") expect_equal(r$confidence, "medium")
  }
})

test_that("invalid feature names and duplicate unit names fail explicitly", {
  m <- paired_counts()
  for (name in c(NA_character_, "", " ")) {
    rownames(m)[1] <- name
    expect_error(EstimateSex(m, species = "human"), "missing or empty")
  }
  m <- cbind(paired_counts(), paired_counts())
  expect_error(EstimateSex(m, species = "human", per_cell = TRUE), "unique")
})

test_that("RDS paths redispatch to matrices and containers with all arguments", {
  m <- paired_counts()
  path <- tempfile(fileext = ".RDS")
  on.exit(unlink(path), add = TRUE)
  for (x in list(m, as.data.frame(m), Matrix::Matrix(m, sparse = TRUE))) {
    saveRDS(x, path)
    expect_equal(
      EstimateSex(path,
        species = "human",
        annotation = annotation_fixture(), input_scale = "normalised"
      ),
      EstimateSex(x,
        species = "human",
        annotation = annotation_fixture(), input_scale = "normalised"
      )
    )
  }
  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = m),
    colData = data.frame(donor = "donor1", row.names = colnames(m))
  )
  saveRDS(se, path)
  expect_equal(
    EstimateSex(path,
      species = "human", assay = "raw", group = "donor",
      input_scale = "normalised"
    ),
    EstimateSex(se,
      species = "human", assay = "raw", group = "donor",
      input_scale = "normalised"
    )
  )
})

test_that("an exact matrix path selects its own triplet in a shared directory", {
  d <- tempfile()
  dir.create(d)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- paired_counts()
  for (prefix in c("A_", "B_")) {
    Matrix::writeMM(Matrix::Matrix(m, sparse = TRUE), file.path(d, paste0(prefix, "matrix.mtx")))
    writeLines(rownames(m), file.path(d, paste0(prefix, "features.tsv")))
    writeLines(colnames(m), file.path(d, paste0(prefix, "barcodes.tsv")))
  }
  expect_error(Read10xDir(d), "More than one")
  expect_equal(EstimateSex(file.path(d, "A_matrix.mtx"), species = "human")$verdict, "male")
})


test_that("dotted gene symbols remain distinct from Ensembl versions", {
  ids <- c("AL627309.1", "AL627309.5", "ENSG00000129824.2")
  expect_equal(
    SexSeek:::.strip_version(ids),
    c("AL627309.1", "AL627309.5", "ENSG00000129824")
  )
  m <- paired_counts()
  extra <- matrix(c(100, 200),
    ncol = 1,
    dimnames = list(ids[1:2], colnames(m))
  )
  expect_equal(EstimateSex(rbind(m, extra), species = "human")$verdict, "male")
})

test_that("core estimation works with explicit species and no genevintage annotation", {
  m <- matrix(0,
    nrow = 10, ncol = 1,
    dimnames = list(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "TMSB4Y", "XIST", "TSIX", "JPX", "ACTB"), "sample")
  )
  m[c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "TMSB4Y"), 1] <- 100
  m[c("ACTB", "JPX"), 1] <- 1000
  out <- EstimateSex(m, species = "Homo sapiens")
  expect_s3_class(out, "data.frame")
  expect_identical(out$verdict[[1]], "male")
})
