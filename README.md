# SexSeek <img src="man/figures/logo.png" align="right" height="139" alt="SexSeek logo" />

`SexSeek` estimates genetic sex from bulk or single-cell expression counts, using
sex-chromosome marker expression and X-inactivation markers. It predicts the matrix
 to be coming from `male`, `female`, `uncertain`, `possible_mixed` or `unknown`, each with 
confidennce metrics.

## Install

```r
install.packages("SexSeek")

# development version
remotes::install_github("saketlab/SexSeek")
```

## Use

```r
library(SexSeek)

path <- system.file("extdata", "example_counts.tsv.gz", package = "SexSeek")
counts <- as.matrix(read.delim(path, row.names = 1))
EstimateSex(counts, group = colnames(counts))

EstimateSex(counts)                          # matrix, sparse matrix, data frame
EstimateSex("counts.tsv.gz")                 # delimited file
EstimateSex("filtered_feature_bc_matrix/")   # 10x MatrixMarket directory
EstimateSex("sample.h5")                     # 10x HDF5
EstimateSex(seurat_obj, group = "donor")     # aggregate by a metadata column
```

Species is detected from the gene identifiers; pass `species = "mouse"` (common
or scientific name) to set it. 

## Genevintage annotations

You can provide your own set of sex genes or fetch it from [genevintage](https://github.com/saketlab/genevintage):

```r
install.packages("genevintage", repos = "https://saketlab.r-universe.dev")
sex_genes <- SexChromosomeGenes("human", release = 116, source = "ensembl")
EstimateSex(counts, species = "human", annotation = sex_genes)
EstimateSex(counts, species = "human", annotation = "genevintage")
```
