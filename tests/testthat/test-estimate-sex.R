core <- function(role) {
  p <- PanelFor("Homo sapiens", tier = "core")
  p$gene_id[p$role == role]
}

# Filler carries the library size so that marker CPM is realistic rather than
# every marker dominating a tiny matrix.
synth <- function(y = 0, xist = 0, filler = 1e5, n_y = 6L) {
  yg <- core("Y")
  xg <- core("inactivation")
  # Autosomal filler carries Ensembl-style identifiers because that is what a
  # real matrix looks like, and species detection reads the whole rowname set.
  fill <- sprintf("ENSG9%08d", seq_len(50))
  ids <- c(yg, xg, fill)
  m <- matrix(0, nrow = length(ids), ncol = 1, dimnames = list(ids, "s1"))
  if (n_y > 0) m[yg[seq_len(n_y)], 1] <- y
  m[xg, 1] <- xist
  m[fill, 1] <- filler / 50
  m
}

test_that("Y signal without XIST calls male", {
  r <- EstimateSex(synth(y = 50, xist = 0))
  expect_equal(r$verdict, "male")
  expect_equal(r$n_y_core_detected, 6)
})

test_that("XIST without Y calls female", {
  r <- EstimateSex(synth(y = 0, xist = 200))
  expect_equal(r$verdict, "female")
})

test_that("both signals report a mixture rather than picking one", {
  r <- EstimateSex(synth(y = 50, xist = 200))
  expect_equal(r$verdict, "possible_mixed")
})

test_that("neither signal is uncertain, never female", {
  r <- EstimateSex(synth(y = 0, xist = 0))
  expect_equal(r$verdict, "uncertain")
  expect_false(r$verdict == "female")
})

test_that("a single Y marker is not enough", {
  # One marker fires on ambient RNA or a mismap; the floor is two.
  r <- EstimateSex(synth(y = 50, xist = 0, n_y = 1L))
  expect_equal(r$verdict, "uncertain")
})

test_that("shallow libraries return unknown", {
  r <- EstimateSex(synth(y = 5, xist = 0, filler = 100))
  expect_equal(r$verdict, "unknown")
  expect_equal(r$flags, "low_depth")
})

test_that("species is detected from Ensembl identifiers", {
  d <- DetectSpecies(core("Y"))
  expect_equal(d$species, "Homo sapiens")
  expect_equal(d$method, "id_prefix")
})

test_that("unsupported species returns unknown, not a guess", {
  s <- SexSpecies()
  un <- s$scientific_name[s$status == "unsupported"]
  skip_if(length(un) == 0, "no unsupported species in the registry")
  m <- synth(y = 50)
  rownames(m) <- paste0("G", seq_len(nrow(m)))
  r <- EstimateSex(m, species = un[1])
  expect_equal(r$verdict, "unknown")
})

test_that("grouping splits columns into units", {
  m <- cbind(synth(y = 50), synth(y = 0, xist = 200))
  colnames(m) <- c("a", "b")
  r <- EstimateSex(m, group = c("g1", "g2"))
  expect_equal(nrow(r), 2)
  expect_equal(sort(r$verdict), c("female", "male"))
})

test_that("sparse and dense matrices agree", {
  m <- synth(y = 50)
  sp <- Matrix::Matrix(m, sparse = TRUE)
  expect_equal(EstimateSex(m)$verdict, EstimateSex(sp)$verdict)
})

test_that("a species with Y markers but no XIST never calls female", {
  # Most non-human eutherians in the panel are male-detectors: Ensembl has no
  # XIST ortholog for them, so absence of Y must not become a female call.
  p <- PanelFor("Rattus norvegicus", tier = NULL)
  expect_equal(sum(p$role == "inactivation"), 0)

  fill <- sprintf("ENSRNOG9%08d", seq_len(50))
  m <- matrix(0,
    nrow = 50, ncol = 1,
    dimnames = list(fill, "s1")
  )
  m[, 1] <- 2000
  r <- EstimateSex(m)
  expect_equal(r$species, "Rattus norvegicus")
  expect_equal(r$verdict, "uncertain")
  expect_match(r$flags, "no_inactivation_marker")
  # Rat is experimental, so the confidence must be capped too.
  expect_true(r$confidence %in% c("medium", "low"))
})

test_that("a species with XIST but no Y never calls male", {
  # Cat is the mirror case: an inactive-X lncRNA and no assembled Y.
  p <- PanelFor("Felis catus", tier = NULL)
  expect_equal(sum(p$role == "Y"), 0)
  expect_gt(sum(p$role == "inactivation"), 0)

  xg <- p$gene_id[p$role == "inactivation"]
  ids <- c(xg, sprintf("ENSFCTG9%08d", seq_len(50)))
  m <- matrix(0, nrow = length(ids), ncol = 1, dimnames = list(ids, "s1"))
  m[xg, 1] <- 400
  m[grep("^ENSFCTG9", ids), 1] <- 2000
  expect_equal(EstimateSex(m)$verdict, "female")
})

test_that("ZW species invert: the sex-specific chromosome means female", {
  reg <- SexSpecies()
  expect_equal(reg$heterogametic[reg$scientific_name == "Gallus gallus"], "female")
  p <- PanelFor("Gallus gallus", tier = NULL)
  w <- p[p$role == "W" & p$tier == "secondary", ]
  expect_gt(nrow(w), 0)
})

test_that("both signals present is possible_mixed", {
  r <- EstimateSex(synth(y = 50, xist = 200))
  expect_equal(r$verdict, "possible_mixed")
})

test_that("a registry species with no panel returns unknown, not an error", {
  # The registry carries species the panel build could not serve (scaffold-level
  # assemblies, unsupported rows kept for explanation). EstimateSex must explain
  # itself rather than raise.
  s <- SexSpecies()
  p <- SexPanels()
  orphan <- setdiff(s$scientific_name, p$scientific_name)
  skip_if(length(orphan) == 0, "every registry species has panel rows")

  m <- matrix(2000,
    nrow = 50, ncol = 1,
    dimnames = list(sprintf("GENE%03d", 1:50), "s1")
  )
  r <- expect_no_error(EstimateSex(m, species = orphan[1]))
  expect_equal(r$verdict, "unknown")
})

test_that("NA values in the matrix do not poison the call", {
  # Published FPKM tables are full of NA. A missing value means the gene was
  # not quantified, so it contributes no evidence — it must not turn the whole
  # unit's total into NA.
  m <- synth(y = 50)
  fill <- grep("^ENSG9", rownames(m))
  m[fill[1:20], 1] <- NA_real_
  r <- EstimateSex(m)
  expect_equal(r$verdict, "male")
  expect_false(is.na(r$qc_score))
})

test_that("symbol detection refuses rather than guessing between mammals", {
  # Regression: quail's 5-gene W panel once let a few coincidental matches beat
  # human. The deeper problem is that X-linked symbols are conserved — human
  # symbols score 0.44 against naked mole-rat — so no mammal is separable this
  # way. Refusing is correct; a wrong species invalidates every call after it.
  p <- SexPanels()
  human_symbols <- p$gene_name[
    p$scientific_name == "Homo sapiens" & nzchar(p$gene_name)
  ]
  d <- DetectSpecies(head(human_symbols, 300))
  expect_null(d$species)
  expect_equal(d$method, "ambiguous")
})

test_that("identifiers still detect cleanly where symbols cannot", {
  # The point of refusing on symbols is that the ID path is unambiguous.
  p <- SexPanels()
  ids <- p$gene_id[p$scientific_name == "Sus scrofa"]
  d <- DetectSpecies(head(ids, 300))
  expect_equal(d$species, "Sus scrofa")
  expect_equal(d$method, "id_prefix")
})

test_that("a few coincidental symbols detect nothing at all", {
  d <- DetectSpecies(c("ACTB", "GAPDH", "ZFX", "KDM5C", "TUBB"))
  expect_null(d$species)
})

test_that("a gene-symbol matrix is not mistaken for an insect", {
  # Regression: eleven human genes start with AGAP and many start with NA, so a
  # count-only prefix rule called Anopheles or Apis on human symbol matrices.
  syms <- c(
    "A1BG", "A2M", "NACA", "NAMPT", "NAA10", "NAT1", "NANOG", "NAB1", "NARS1",
    "NADK", "NAPA", "NAGK", "AGAP1", "AGAP2", "AGAP3", "AGAP4", "AGAP5",
    "AGAP6", "AGAP9", "AGAP11", "GAPDH", "ACTB", "TP53", "EGFR", "MYC"
  )
  syms <- c(syms, sprintf("GENE%04d", seq_len(500)))
  d <- DetectSpecies(syms)
  expect_false(isTRUE(d$species %in% c("Anopheles gambiae", "Apis mellifera")))
})

test_that("registry id prefixes are all usable", {
  s <- SexSpecies()
  p <- s$id_prefix[!is.na(s$id_prefix) & nzchar(s$id_prefix) & s$id_prefix != "NA"]
  expect_true(all(nchar(p) >= 4))
})

test_that("a 10x MatrixMarket directory loads and scores", {
  p <- PanelFor("Homo sapiens", tier = "core")
  ids <- c(p$gene_id, sprintf("ENSG9%08d", seq_len(50)))
  v <- rep(0, length(ids))
  v[p$role == "Y"] <- 500
  v[grep("^ENSG9", ids)] <- 2000
  m <- Matrix::Matrix(matrix(v, ncol = 1), sparse = TRUE)

  d <- file.path(tempdir(), "tenx")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  Matrix::writeMM(m, file.path(d, "matrix.mtx"))
  writeLines(ids, file.path(d, "features.tsv"))
  writeLines("AAACCCAAGAAACACT-1", file.path(d, "barcodes.tsv"))

  r <- EstimateSex(d)
  expect_equal(r$species, "Homo sapiens")
  expect_equal(r$verdict, "male")
})

test_that("an inconsistent triplet errors instead of mislabelling genes", {
  d <- file.path(tempdir(), "tenx_bad")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- Matrix::Matrix(matrix(1:6, nrow = 3), sparse = TRUE)
  Matrix::writeMM(m, file.path(d, "matrix.mtx"))
  writeLines(c("G1", "G2"), file.path(d, "features.tsv")) # one too few
  writeLines(c("B1", "B2"), file.path(d, "barcodes.tsv"))
  expect_error(EstimateSex(d), "inconsistent")
})

test_that("confidence is capped on an experimental panel", {
  # Pig is experimental: orthology-derived secondary markers, never calibrated.
  # A correct verdict there still must not be reported as high confidence.
  p <- PanelFor("Sus scrofa", tier = NULL)
  y <- p$gene_id[p$role == "Y" & p$tier == "secondary"]
  skip_if(length(y) < 2, "pig has too few secondary Y markers")

  ids <- c(y, sprintf("ENSSSCG9%08d", seq_len(50)))
  m <- matrix(0, nrow = length(ids), ncol = 1, dimnames = list(ids, "s1"))
  m[y, 1] <- 500
  m[grep("^ENSSSCG9", ids), 1] <- 2000
  r <- EstimateSex(m)
  expect_equal(r$verdict, "male")
  expect_true(r$confidence %in% c("medium", "low"))
  expect_true(grepl("experimental_panel", r$flags))
})

test_that("a normalised matrix is flagged and not called with high confidence", {
  m <- synth(y = 50)
  m <- m + 0.5 # FPKM-like: no longer integer counts
  r <- EstimateSex(m)
  expect_true(grepl("normalised_input", r$flags))
  expect_true(r$confidence %in% c("medium", "low"))
})

test_that("a missing marker row is distinguished from a species without one", {
  # Human has XIST; dropping the row from the input is a recoverable input
  # problem, not the permanent limit that rat has.
  m <- synth(y = 0, xist = 0)
  xg <- PanelFor("Homo sapiens", tier = "core")
  xg <- xg$gene_id[xg$role == "inactivation"]
  m <- m[setdiff(rownames(m), xg), , drop = FALSE]
  r <- EstimateSex(m)
  expect_equal(r$flags, "inactivation_marker_missing_from_input")
  expect_match(r$notes, "unfiltered gene set")
})

test_that("duplicate matrix rows are rejected rather than silently selected", {
  m <- synth(y = 50)
  dup <- rbind(m, m[1, , drop = FALSE])
  expect_error(EstimateSex(dup), "Duplicate or ambiguous")
})

test_that("GEO-prefixed 10x triplets are found", {
  # Real deposits ship GSM9323718_6854-37_matrix.mtx.gz, not matrix.mtx.gz.
  # Anchoring the pattern at the start of the name missed every one of them.
  d <- file.path(tempdir(), "prefixed")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)

  m <- Matrix::Matrix(matrix(c(5, 0, 3), ncol = 1), sparse = TRUE)
  Matrix::writeMM(m, file.path(d, "GSM9323718_6854-37_matrix.mtx"))
  writeLines(c("G1", "G2", "G3"), file.path(d, "GSM9323718_6854-37_features.tsv"))
  writeLines("AAAC-1", file.path(d, "GSM9323718_6854-37_barcodes.tsv"))

  got <- Read10xDir(d)
  expect_equal(dim(got), c(3L, 1L))
  expect_equal(rownames(got), c("G1", "G2", "G3"))
})

test_that("two samples in one directory error rather than silently picking one", {
  d <- file.path(tempdir(), "twosamples")
  dir.create(d, showWarnings = FALSE)
  on.exit(unlink(d, recursive = TRUE), add = TRUE)
  m <- Matrix::Matrix(matrix(1, ncol = 1), sparse = TRUE)
  for (p in c("A", "B")) {
    Matrix::writeMM(m, file.path(d, paste0(p, "_matrix.mtx")))
    writeLines("G1", file.path(d, paste0(p, "_features.tsv")))
    writeLines("BC1", file.path(d, paste0(p, "_barcodes.tsv")))
  }
  expect_error(Read10xDir(d), "More than one")
})

test_that("whitespace-delimited tables are read, not mangled into one column", {
  # A .txt deposit is as likely to be space-delimited as tab-delimited. Read as
  # TSV it becomes a single character column and fails as "no numeric columns".
  f <- file.path(tempdir(), "spaced.txt")
  on.exit(unlink(f), add = TRUE)
  writeLines(c('"c1" "c2"', '"ENSG1" 5 7', '"ENSG2" 0 3'), f)
  got <- SexSeek:::.as_counts.character(f)
  expect_equal(dim(got), c(2L, 2L))
  expect_equal(rownames(got), c("ENSG1", "ENSG2"))
  expect_equal(unname(got["ENSG1", ]), c(5, 7))
})
