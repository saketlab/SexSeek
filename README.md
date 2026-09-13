# SexSeek

Estimates genetic sex from bulk or single-cell expression counts, using
sex-chromosome marker expression and X-inactivation markers.

## Install

```r
remotes::install_github("saketlab/SexSeek")
```

## Use

```r
library(SexSeek)

EstimateSex(counts)                          # matrix, sparse matrix, data frame
EstimateSex("counts.tsv.gz")                 # delimited file
EstimateSex("filtered_feature_bc_matrix/")   # 10x MatrixMarket directory
EstimateSex("sample.h5")                     # 10x HDF5
EstimateSex(seurat_obj, group = "donor")     # aggregate by a metadata column
```

## Genevintage annotations

You can provide your own set of sex genes or fetch it from [genevintage](https://github.com/saketlab/genevintage):
```r
sex_genes <- SexChromosomeGenes("human", release = 116, source = "ensembl")
EstimateSex(counts, species = "human", annotation = sex_genes)
EstimateSex(counts, species = "human", annotation = "genevintage")
```
