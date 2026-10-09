annotation_fixture <- function() {
  p <- PanelFor("Homo sapiens", tier = "core")
  gp <- GametologPairs()
  gp <- gp[gp$scientific_name == "Homo sapiens" & gp$y_tier == "core", ]
  unique(rbind(
    data.frame(
      id = p$gene_id, name = p$gene_name,
      chr = ifelse(p$role == "inactivation", "X", p$role),
      biotype = "protein_coding"
    ),
    data.frame(
      id = gp$x_gene_id, name = gp$x_gene_name,
      chr = "X", biotype = "protein_coding"
    )
  ))
}

paired_counts <- function() {
  a <- annotation_fixture()
  m <- matrix(0,
    nrow = nrow(a) + 1L, ncol = 1,
    dimnames = list(c(a$id, "AUTOSOME"), "sample")
  )
  m[a$id[a$chr == "Y"], ] <- 50
  m[a$id[a$chr == "X" & a$name != "XIST"], ] <- 100
  m["AUTOSOME", ] <- 100000
  m
}
