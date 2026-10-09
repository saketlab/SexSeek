test_that("ratios cannot bypass the expression floor", {
  m <- paired_counts()
  m[, ] <- 0
  m["AUTOSOME", ] <- 1e8
  gp <- subset(GametologPairs(), scientific_name == "Homo sapiens" & y_tier == "core")
  m[gp$y_gene_id[1:2], ] <- 1
  expect_equal(EstimateSex(m, species = "human")$verdict, "uncertain")
  m[gp$y_gene_id, ] <- 0
  p <- PanelFor("Homo sapiens")
  m[p$gene_id[p$role == "inactivation"], ] <- 1
  expect_equal(EstimateSex(m, species = "human")$verdict, "uncertain")
})

test_that("missing diagnostic values cause abstention rather than false negatives", {
  m <- paired_counts()
  p <- PanelFor("Homo sapiens")
  m[p$gene_id[p$role == "Y"], ] <- NA_real_
  m[p$gene_id[p$role == "inactivation"], ] <- 1000
  r <- EstimateSex(m, species = "human")
  expect_equal(r$verdict, "unknown")
  expect_match(r$flags, "missing_marker_values")
  m <- paired_counts()
  m[p$gene_id[p$role == "inactivation"], ] <- NA_real_
  expect_equal(EstimateSex(m, species = "human")$verdict, "unknown")
})

test_that("a fitted model requires positive evidence for either sex", {
  old <- SexSeek:::.cache$models
  on.exit(assign("models", old, envir = SexSeek:::.cache), add = TRUE)
  m <- paired_counts()
  m[, ] <- 0
  m["AUTOSOME", ] <- 1e5
  for (intercept in c(-20, 20)) {
    fake <- data.frame(
      scientific_name = "Homo sapiens",
      feature = SexSeek:::.model_features,
      estimate = c(intercept, 0, 0, 0, 0, 0)
    )
    assign("models", fake, envir = SexSeek:::.cache)
    r <- EstimateSex(m, species = "human", model = "logistic")
    expect_equal(r$verdict, "uncertain")
    expect_equal(r$flags, "model_without_positive_evidence")
  }
})

test_that("invalid expression values and missing donor labels are rejected", {
  m <- paired_counts()
  for (bad in c(-1, Inf, -Inf)) {
    m["AUTOSOME", ] <- bad
    expect_error(EstimateSex(m, species = "human"), "non-negative and finite")
    expect_error(
      EstimateSex(Matrix::Matrix(m, sparse = TRUE), species = "human"),
      "non-negative and finite"
    )
  }
  m <- cbind(paired_counts(), paired_counts())
  expect_error(EstimateSex(m, species = "human", group = c("donor", NA)), "missing")
})

test_that("per-column evidence is explicitly unvalidated", {
  r <- EstimateSex(paired_counts(), species = "human", per_cell = TRUE)
  expect_equal(r$verdict, "male")
  expect_equal(r$confidence, "medium")
  expect_match(r$flags, "per_cell_unvalidated")
})
