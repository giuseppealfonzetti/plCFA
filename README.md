
<!-- README.md is generated from README.Rmd. Please edit that file -->

# plCFA

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/plCFA)](https://CRAN.R-project.org/package=plCFA)
<!-- badges: end -->

The plCFA deals with the estimation of confirmatory factor models for
ordinal data.

## Installation

You can install the development version of plCFA like so:

``` r
devtools::install_github("giuseppealfonzetti/stIsing")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(plCFA)

# Setup the model
seed <- 1
set.seed(seed)
p <- 50; q <- 5; n <- 5000
constrMat <- build_constrMat(p,q,'simple')
true_load <- gen_loadings(FIXED = NULL, CONSTRAINT_MAT = constrMat, SEED = seed)
true_tau <- c(-1.2, 0, 1.2)
true_latent <- matrix(.7,q,q); diag(true_latent) <- 1
true_theta <- get_theta(rep(true_tau, p), true_load, true_latent, cat, constrMat, 0)

# generate the data
manifest <- gen_URV_data(n, true_load, true_tau, true_latent)

# fit the model
fit <- fit_plCFA(
    DATA_LIST = list('DATA' = manifest, 'CONSTRMAT' = constrMat, 'CORRFLAG'=1),
    METHOD = 'ucminf',
    INIT = NULL
)
#> 1. Initialising at default values
#> 2. Optimising with ucminf...
#> 3. Done! (92.08 secs)

mean((fit$theta_init-true_theta)^2)
#> [1] 0.07241336
mean((fit$theta-true_theta)^2)
#> [1] 0.0005276809
```
