# colnames() of a list is NULL, so lists are scored per element
score_ex <- function(ex) {
  EstimateSex(ex$counts, species = ex$species, group = colnames(ex$counts))
}

# r-universe served a 0.1.3 build before aggregate_transcripts landed
has_tx <- function() {
  requireNamespace("genevintage", quietly = TRUE) &&
    "aggregate_transcripts" %in% getNamespaceExports("genevintage")
}

score_tx <- function(ex) {
  genes <- genevintage::aggregate_transcripts(ex$counts,
    mapping = ex$tx2gene, label = "id", unmapped = "keep", quiet = TRUE
  )
  EstimateSex(genes, species = ex$species, group = colnames(genes))
}

exs <- lapply(stats::setNames(nm = ExampleData()$dataset), ExampleData)
calls <- lapply(Filter(function(ex) is.null(ex$tx2gene), exs), score_ex)
tx_calls <- if (has_tx()) score_tx(exs$mouse_tx)

test_that("ExampleData lists and loads every shipped dataset", {
  for (ex in exs) {
    units <- if (is.list(ex$counts)) names(ex$counts) else colnames(ex$counts)
    expect_identical(names(ex$recorded), units)
    expect_false(anyNA(ex$recorded))
  }
  expect_false(is.null(exs$mouse_tx$tx2gene))
  expect_error(ExampleData("nope"), "Unknown dataset")
})

test_that("shipped excerpts reproduce the recorded-sex results", {
  expected <- list(
    human_bulk = c(15, 0, 0), mouse_bulk = c(24, 0, 0), human_sc = c(10, 0, 2),
    mouse_sc = c(12, 0, 0), rat_bulk = c(4, 0, 4), zebrafish_bulk = c(0, 0, 16)
  )
  for (d in names(expected)) {
    got <- SummariseComparison(CompareSex(calls[[d]], exs[[d]]$recorded))
    expect_equal(unlist(got[c("agree", "disagree", "abstain")]), expected[[d]],
      ignore_attr = TRUE, label = d
    )
  }
})

test_that("rat abstains on females and zebrafish has no panel", {
  rat <- CompareSex(calls$rat_bulk, exs$rat_bulk$recorded)
  expect_true(all(rat$verdict[rat$recorded_sex == "female"] == "uncertain"))
  expect_true(all(calls$zebrafish_bulk$verdict == "unknown"))
})

test_that("transcript-level excerpt aggregates with genevintage", {
  skip_if_not(has_tx())
  cmp <- CompareSex(tx_calls, exs$mouse_tx$recorded)
  expect_identical(
    SummariseComparison(cmp),
    data.frame(samples = 8L, agree = 4L, disagree = 0L, abstain = 4L)
  )
  females <- cmp$recorded_sex == "female"
  expect_true(all(grepl("inactivation_marker_missing_from_input", cmp$flags[females])))
})

test_that("inact_measured is FALSE when the marker cannot be measured", {
  expect_true(all(calls$mouse_bulk$inact_measured))
  expect_false(any(calls$rat_bulk$inact_measured))
  expect_false(any(calls$zebrafish_bulk$inact_measured))
  skip_if_not(has_tx())
  expect_false(any(tx_calls$inact_measured))
})

test_that("CompareSex scores only female/male records and real calls", {
  x <- data.frame(
    unit = c("a", "b", "c", "d", "e"),
    verdict = c("female", "male", "uncertain", "male", "female")
  )
  rec <- c(a = "F", b = "female", c = "Male", d = "unknown", e = "M")
  cmp <- CompareSex(x, rec)
  expect_identical(cmp$recorded, c("F", "female", "Male", "unknown", "M"))
  expect_identical(cmp$recorded_sex, c("female", "female", "male", NA, "male"))
  expect_identical(cmp$agrees, c(TRUE, FALSE, NA, NA, FALSE))
  expect_error(CompareSex(x, c("F", "M")), "named")
  expect_error(CompareSex(data.frame(x = 1), rec), "EstimateSex")
})

test_that("SummariseComparison counts agree, disagree and abstain", {
  cmp <- CompareSex(
    data.frame(unit = c("a", "b", "c"), verdict = c("female", "male", "uncertain")),
    c(a = "female", b = "female", c = "male")
  )
  expect_identical(
    SummariseComparison(cmp),
    data.frame(samples = 3L, agree = 1L, disagree = 1L, abstain = 1L)
  )
  expect_identical(SummariseComparison(list(x = cmp, y = cmp))$comparison, c("x", "y"))
  expect_error(SummariseComparison(list(cmp)), "names")
  expect_error(SummariseComparison(data.frame(a = 1)), "CompareSex")
})

test_that("a list is scored one unit per element", {
  libs <- exs$mouse_sc$counts[1:2]
  res <- EstimateSex(libs, species = "Mus musculus")
  expect_identical(res$unit, names(libs))
  one <- EstimateSex(libs[[1]], species = "Mus musculus")
  one$unit <- names(libs)[1]
  expect_equal(res[1, ], one, ignore_attr = TRUE)
  expect_error(EstimateSex(libs, group = "x"), "group")
  expect_error(EstimateSex(unname(libs)), "names")
})

test_that("MarkerCellFraction agrees between list and grouped matrix", {
  libs <- exs$mouse_sc$counts[1:2]
  a <- MarkerCellFraction(libs, "mouse")
  m <- do.call(cbind, libs)
  b <- MarkerCellFraction(m, "Mus musculus", group = rep(names(libs), vapply(libs, ncol, 1L)))
  expect_equal(a, b)
  expect_true(all(a$fraction >= 0 & a$fraction <= 1))
  expect_error(MarkerCellFraction(m, "Danio rerio"), "No sex-specific markers")
  expect_error(MarkerCellFraction(m, "not a species"), "Unknown species")
  expect_error(MarkerCellFraction(unname(libs), "mouse"), "names")
})

test_that("plots are ggplots with circle/triangle shapes and refuse a species with no scores", {
  p <- PlotSexEvidence(calls$mouse_bulk, exs$mouse_bulk$recorded)
  expect_s3_class(p, "ggplot")
  expect_true("recorded_sex" %in% names(p$data))
  b <- ggplot2::ggplot_build(p)$data[[1]]
  sex <- p$data$recorded_sex
  expect_true(all(b$shape[sex == "female"] == 16) && all(b$shape[sex == "male"] == 17))
  expect_s3_class(ggplot2::ggplot_build(PlotSexEvidence(calls$mouse_bulk)), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(PlotSexEvidence(calls$rat_bulk, exs$rat_bulk$recorded)), "ggplot_built")
  f <- PlotCellFraction(exs$mouse_sc$counts, "Mus musculus", exs$mouse_sc$recorded)
  expect_setequal(f$data$unit, names(exs$mouse_sc$counts))
  expect_s3_class(ggplot2::ggplot_build(f), "ggplot_built")
  expect_s3_class(ggplot2::ggplot_build(PlotCellFraction(exs$mouse_sc$counts, "Mus musculus")), "ggplot_built")
  expect_error(PlotSexEvidence(calls$zebrafish_bulk), "No marker scores")
})
