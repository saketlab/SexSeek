test_that("real seqout fixtures reproduce recorded calls", {
  path <- testthat::test_path("fixtures", "validation.rds")
  skip_if_not(file.exists(path), "extended validation fixture not generated")
  fixture <- readRDS(path)
  expect_gt(length(fixture$samples), 100L)
  checked <- 0L
  for (x in fixture$samples) {
    for (mode in names(x$expected)) {
      ann <- if (mode == "genevintage") fixture$annotations[[x$species]] else NULL
      got <- SexSeek::EstimateSex(x$counts,
        species = x$species, annotation = ann,
        input_scale = x$input_scale
      )
      expect_equal(got$verdict, x$expected[[mode]]$verdict,
        info = paste(x$dataset, x$sample, mode)
      )
      expect_equal(got$confidence, x$expected[[mode]]$confidence,
        info = paste(x$dataset, x$sample, mode)
      )
      checked <- checked + 1L
    }
  }
  expect_gt(checked, 150L)
})
