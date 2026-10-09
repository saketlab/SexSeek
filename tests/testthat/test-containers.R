test_that("SingleCellExperiment RDS preserves donor aggregation", {
  skip_if_not_installed("SingleCellExperiment")
  m <- cbind(paired_counts(), paired_counts())
  colnames(m) <- c("cell1", "cell2")
  x <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = m),
    colData = data.frame(donor = c("A", "B"), row.names = colnames(m))
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path)
  expect_equal(
    EstimateSex(path, species = "human", group = "donor"),
    EstimateSex(m, species = "human", group = c("A", "B"))
  )
})

test_that("Seurat split counts layers and RDS agree with the matrix path", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  old <- options(Seurat.object.assay.version = "v5")
  on.exit(options(old), add = TRUE)
  male <- Matrix::Matrix(paired_counts(), sparse = TRUE)
  female <- male
  a <- annotation_fixture()
  female[a$id[a$chr == "Y"], ] <- 0
  female[a$id[a$name == "XIST"], ] <- 1000
  colnames(male) <- "cellM"
  colnames(female) <- "cellF"
  x <- SeuratObject::CreateSeuratObject(
    counts = list(b = female, a = male),
    meta.data = data.frame(donor = c("M", "F"), row.names = c("cellM", "cellF"))
  )
  expected <- EstimateSex(cbind(male, female),
    species = "human", group = c("M", "F"),
    model = "logistic", input_scale = "normalised"
  )
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path)
  for (object in list(x, path)) {
    expect_equal(EstimateSex(object,
      species = "human", group = "donor",
      model = "logistic", input_scale = "normalised"
    ), expected)
  }
})

test_that("split layers cannot silently fill missing genes or count cells twice", {
  skip_if_not_installed("SeuratObject", "5.0.0")
  old <- options(Seurat.object.assay.version = "v5")
  on.exit(options(old), add = TRUE)
  a <- Matrix::Matrix(paired_counts(), sparse = TRUE)
  b <- a
  colnames(a) <- "A"
  colnames(b) <- "B"
  x <- SeuratObject::CreateSeuratObject(counts = list(a = a, b = b[-1, , drop = FALSE]))
  expect_error(EstimateSex(x, species = "human"), "different gene sets")
  x <- SeuratObject::CreateSeuratObject(counts = list(a = a, b = a))
  expect_error(EstimateSex(x, species = "human"), "overlapping cells")
})

test_that("metadata grouping follows assay columns rather than metadata order", {
  cells <- c("B", "A")
  meta <- data.frame(donor = c("donorA", "donorB", "donorC"), row.names = c("A", "B", "C"))
  expect_equal(SexSeek:::.group_from("donor", meta, cells), c("donorB", "donorA"))
  expect_error(SexSeek:::.group_from("donor", meta, c("A", "unknown")), "aligned")
})
