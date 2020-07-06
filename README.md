
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MMRFBiolinks

## An R package that extends TCGABiolink package for integrative analysis with MMRF-COMMPASS data

<!-- badges: start -->

<!-- badges: end -->

MMRFBiolinks extends TCGABiolink package for searching, downloading and
analyzing MMRF-COMMPASS data available at the NCI’s Genomic Data Commons
(GDC) Data Portal.

## Installation

Once R (version “4.0”) has been started, you can install the released
version of MMRFBiolinks from GitHub with:

``` r
devtools::install_github("MarzyUnicz/MMRFBiolinks", build_vignettes = TRUE)
library(MMRFBiolinks)
```

## Required libraries

``` r
library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)
library(DT)
library(ggplot2)
```

## Vignettes

A list of all currently integrated vignettes can be obtained through:

``` r
vignette(package="MMRFBiolinks")
```

The best way to view vignettes is in your web browser:

``` r
browseVignettes("MMRFBiolinks")
```

Get the list of the example data sets

``` r
data(package = "MMRFBiolinks")
```

## Citation
