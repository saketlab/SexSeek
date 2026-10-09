test_that("genevintage accepts offline annotations and respects chromosomes", {
  skip_if_not_installed("genevintage")
  a <- annotation_fixture()
  a <- rbind(a, data.frame(
    id = "AUTOSOME", name = "OTHER", chr = "1",
    biotype = "protein_coding"
  ))
  got <- SexChromosomeGenes("human", mapping = a)
  expect_setequal(got$id, setdiff(a$id, "AUTOSOME"))
  m <- paired_counts()
  expect_equal(EstimateSex(m, species = "human", annotation = got)$verdict, "male")
  expect_error(EstimateSex(m, species = "mouse", annotation = got), "species")
})

test_that("annotation resolves changed IDs without promoting arbitrary Y genes", {
  a <- annotation_fixture()
  m <- paired_counts()
  a$id <- paste0("CUSTOM_", seq_len(nrow(a)))
  rownames(m)[seq_len(nrow(a))] <- a$id
  expect_equal(EstimateSex(m, species = "human", annotation = a)$verdict, "male")
  a$chr[a$chr == "Y"] <- "1"
  expect_equal(EstimateSex(m, species = "human", annotation = a)$verdict, "uncertain")
  expect_error(EstimateSex(m, species = "human", annotation = a["id"]), "data frame")
})

test_that("only complete core gametolog pairs contribute to the ratio", {
  m <- paired_counts()
  r <- EstimateSex(m, species = "human")
  expect_equal(r$n_gametolog_pairs, 6L)
  expect_equal(r$gametolog_frac, 1 / 3)
  gp <- subset(GametologPairs(), scientific_name == "Homo sapiens" & y_tier == "core")
  m <- m[setdiff(rownames(m), gp$x_gene_id[1]), , drop = FALSE]
  m[gp$y_gene_id[1], ] <- 100000
  r <- EstimateSex(m, species = "human")
  expect_equal(r$n_gametolog_pairs, 5L)
  expect_equal(r$gametolog_frac, 1 / 3)
  m <- m[setdiff(rownames(m), gp$x_gene_id), , drop = FALSE]
  r <- EstimateSex(m, species = "human")
  expect_equal(r$n_gametolog_pairs, 0L)
  expect_true(is.na(r$gametolog_frac))
})

test_that("SummarizedExperiment forwards annotation", {
  skip_if_not_installed("SummarizedExperiment")
  m <- paired_counts()
  a <- annotation_fixture()
  a$id <- paste0("CUSTOM_", seq_len(nrow(a)))
  rownames(m)[seq_len(nrow(a))] <- a$id
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = m))
  r <- EstimateSex(se, species = "human", annotation = a)
  expect_equal(r$verdict, "male")
})
