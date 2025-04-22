
<!-- README.md is generated from README.Rmd. Please edit that file -->

# svycdiff <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/salernos/svycdiff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/salernos/svycdiff/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Propensity score methods are broadly employed with observational data as
a tool to achieve covariate balance, but how to implement them in
complex surveys is less studied – in particular, when the survey weights
depend on the group variable under comparison.

In this package, we focus on the specific case when sample selection
depends the comparison groups of interest. We implement identification
formulas to properly estimate the *average controlled difference (ACD)*,
or under stronger assumptions, the *population average treatment effect
(ATE)* in outcomes between groups, with appropriate weighting for both
covariate imbalance and generalizability.

This packages also contains the code necessary to reproduce the
motivating data analysis in *“What’s the weight? Estimating controlled
outcome differences in complex surveys for health disparities
research.”* This analysis focuses on data from the National Health and
Nutrition Examination Survey (NHANES), investigating the interplay of
race and social determinants of health when our interest lies in
estimating racial differences in mean telomere length.

## Installation

You can install the development version of svycdiff like so:

``` r
#--- CRAN Version
install.packages("svycdiff")
#--- Development Version
# install.packages("devtools")
devtools::install_github("salernos/svycdiff")
```

## Example

This is a basic example usage via simulated data:

``` r
library(svycdiff)

N <- 1000

dat <- simdat(N)

S <- rbinom(N, 1, dat$pS)

samp <- dat[S == 1,]

y_mod <- Y ~ A * X1

a_mod <- A ~ X1

s_mod <- pS ~ A + X1

fit <- svycdiff(samp, "DR", a_mod, s_mod, y_mod, "gaussian")

fit
```

## Vignette

Once you have `svycdiff` installed, you can type

``` r
vignette("svycdiff")
```

in `R` to bring up a tutorial on `svycdiff` and how to use it. To access
the vignettes in the developer version, please install the package with

``` r
devtools::install_github("salernos/svycdiff", build_vignettes = TRUE)
```

## Questions

For technical details on the method, see please refer to Salerno et
al. (2024+) *“What’s the weight? Estimating controlled outcome
differences in complex surveys for health disparities research.”* To
reproduce the analysis results for the main paper, see
`inst/nhanes.Rmd`. For questions and comments, please contact Stephen
Salerno (<ssalerno@fredhutch.org>).

``` r
sessionInfo()
#> R version 4.4.1 (2024-06-14 ucrt)
#> Platform: x86_64-w64-mingw32/x64
#> Running under: Windows 11 x64 (build 22631)
#> 
#> Matrix products: default
#> 
#> 
#> locale:
#> [1] LC_COLLATE=English_United States.utf8 
#> [2] LC_CTYPE=English_United States.utf8   
#> [3] LC_MONETARY=English_United States.utf8
#> [4] LC_NUMERIC=C                          
#> [5] LC_TIME=English_United States.utf8    
#> 
#> time zone: America/New_York
#> tzcode source: internal
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.4.1    fastmap_1.2.0     cli_3.6.3         tools_4.4.1      
#>  [5] htmltools_0.5.8.1 rstudioapi_0.17.1 yaml_2.3.10       rmarkdown_2.29   
#>  [9] knitr_1.50        xfun_0.52         digest_0.6.37     rlang_1.1.4      
#> [13] evaluate_1.0.3
```
