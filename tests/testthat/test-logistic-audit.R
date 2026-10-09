testthat::test_that("logistic training audit is explicit and conservative", {
  path <- testthat::test_path("..", "..", "data-raw", "logistic-audit.tsv")
  testthat::skip_if_not(file.exists(path))
  audit <- utils::read.delim(path, stringsAsFactors = FALSE)
  testthat::expect_true(all(c(
    "scientific_name", "n_homogametic", "n_heterogametic",
    "n_studies", "eligible", "reason"
  ) %in% names(audit)))
  testthat::expect_true(all(audit$n_homogametic >= 0 & audit$n_heterogametic >= 0))
  testthat::expect_true(all(audit$eligible == FALSE))
  testthat::expect_true(all(grepl("insufficient_", audit$reason)))
})
