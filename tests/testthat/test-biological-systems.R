test_that("mouse core markers and Xist have the expected polarity", {
  p <- PanelFor("Mus musculus")
  gp <- subset(GametologPairs(), scientific_name == "Mus musculus" & y_tier == "core")
  ids <- unique(c(p$gene_id, gp$x_gene_id, "AUTOSOME"))
  m <- matrix(0, length(ids), 1, dimnames = list(ids, "sample"))
  m["AUTOSOME", ] <- 1e5
  m[gp$x_gene_id, ] <- 100
  m[p$gene_id[p$role == "Y"], ] <- 100
  expect_equal(EstimateSex(m, species = "mouse")$verdict, "male")
  m[p$gene_id[p$role == "Y"], ] <- 0
  expect_equal(EstimateSex(m, species = "mouse")$verdict, "uncertain")
  m[p$gene_id[p$role == "inactivation"], ] <- 1000
  expect_equal(EstimateSex(m, species = "mouse")$verdict, "female")
})

test_that("a bird with only one usable W marker abstains", {
  p <- SexSeek:::.scoring_panel("Gallus gallus")
  ids <- c(p$gene_id, "AUTOSOME")
  m <- matrix(0, length(ids), 1, dimnames = list(ids, "sample"))
  m["AUTOSOME", ] <- 1e5
  m[p$gene_id[p$role == "W"], ] <- 100
  r <- EstimateSex(m, species = "chicken")
  expect_equal(r$verdict, "uncertain")
  expect_true(r$confidence %in% c("medium", "low"))
  m[p$gene_id[p$role == "W"], ] <- 0
  expect_equal(EstimateSex(m, species = "chicken")$verdict, "uncertain")
})

test_that("annotation coverage is distinguished from observed zero expression", {
  m <- paired_counts()
  a <- annotation_fixture()
  m[a$id[a$chr == "Y"], ] <- 0
  m[a$id[a$name == "XIST"], ] <- 1000
  complete <- EstimateSex(m, species = "human", annotation = a)
  incomplete <- EstimateSex(m, species = "human", annotation = a[a$chr != "Y", ])
  expect_false(grepl("annotation_markers_missing", complete$flags))
  expect_match(incomplete$flags, "annotation_markers_missing")
  expect_equal(incomplete$confidence, "medium")
  expect_match(incomplete$notes, "supplied annotation")
})


test_that("two independent W markers have female polarity in the decision table", {
  f <- list(
    total = 1e5, gfrac = NA_real_, ifrac = NA_real_, n_spec = 2L,
    spec_cpm = 1000, inact_cpm = 0, can_call_homo = FALSE,
    het = "female", homo = "male", species_has_inact = FALSE
  )
  expect_equal(SexSeek:::.decide(f)$verdict, "female")
  f$n_spec <- 0L
  expect_equal(SexSeek:::.decide(f)$verdict, "uncertain")
})
