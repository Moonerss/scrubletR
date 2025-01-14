
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scrubletR

<!-- badges: start -->
<!-- badges: end -->

The goal of scrubletR is to run the scrublet in R.

## Installation

You can install the released version of scrubletR from
[Github](https://github.com/Moonerss/scrubletR) with:

``` r
remotes::install_github("Moonerss/scrubletR")
```

## Prerequisites

-   R package: Matrix, reticulate, Seurat

-   Python module: scrublet (v0.2.3)

## Using

### One step run scrublet

``` r
load("~/samples.RData")

res <- scrublet_R(seurat_obj = samples)

## use specific threshold
res <- scrublet_R(seurat_obj = samples, threshold = 0.25)
# Preprocessing...
# Simulating doublets...
# Embedding transcriptomes using PCA...
# Calculating doublet scores...
# Automatically set threshold at doublet score = 0.57
# Detected doublet rate = 0.1%
# Estimated detectable doublet fraction = 0.6%
# Overall doublet rate:
#   Expected   = 6.0%
#   Estimated  = 12.1%
# Elapsed time: 32.5 seconds
```

``` r
head(res@meta.data)
#                       orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACCCAAGAGTCACG-1 SeuratProject       5985         1844   7.552214
# AAACCCAAGCCGTCGT-1 SeuratProject      11923         2532  11.725237
# AAACCCAAGGACACTG-1 SeuratProject       1728          783  22.858796
# AAACCCACAAAGGTTA-1 SeuratProject      10676         2649  22.798801
# AAACCCACAGACGCTC-1 SeuratProject       2235          690  55.078300
# AAACCCACAGGAACCA-1 SeuratProject      14519         3402   4.979682
#                    doublet_scores predicted_doublets
# AAACCCAAGAGTCACG-1     0.01636424              FALSE
# AAACCCAAGCCGTCGT-1     0.03400332              FALSE
# AAACCCAAGGACACTG-1     0.22724604              FALSE
# AAACCCACAAAGGTTA-1     0.02451680              FALSE
# AAACCCACAGACGCTC-1     0.01955671              FALSE
# AAACCCACAGGAACCA-1     0.01955671              FALSE
```

### Step by step run scrublet  

```r
## init scrublet
scrublet_obj = get_init_scrublet(seurat_obj = samples)

## plot histogram
plot_histogram(scrublet_obj)

## update threshold
scrublet_obj = call_doublets(scrublet_obj, threshold = 0.25)

## plot to check again
plot_histogram(scrublet_obj)

## add info to seurat obj
samples[["doublet_scores"]] <- scrublet_obj$doublet_scores_obs_
samples[["predicted_doublets"]] <- scrublet_obj$predicted_doublets_
```
