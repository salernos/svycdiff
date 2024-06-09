
<!-- README.md is generated from README.Rmd. Please edit that file -->

# svycdiff <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->

[![R-CMD-check](https://github.com/salernos/svycdiff/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/salernos/svycdiff/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A basic descriptive question in statistics asks whether there are
differences in mean outcomes between groups based on levels of a
discrete covariate (e.g., racial disparities in health outcomes,
differences in opinion based on political party identification,
heterogeneity in educational outcomes for students in urban vs. rural
school districts, etc). When this categorical covariate of interest is
correlated with other factors related to the outcome, however, direct
comparisons may lead to erroneous estimates and invalid inferential
conclusions without appropriate adjustment.

In this package, we focus on the specific case when sample selection
depends the comparison groups of interest. We implement identification
formulas to properly estimate the *average controlled difference (ACD)*,
or under stronger assumptions, the *population average treatment effect
(ATE)* in outcomes between groups, with appropriate weighting for
covariate imbalance and generalizability.

This packages also contains the code necessary to reproduce the
motivation data analysis in *“What’s the weight? Estimating controlled
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

S <- rbinom(N, 1, dat$P_S_cond_AX)

samp <- dat[S == 1,]

y_mod <- Y ~ A * X

a_mod <- A ~ X

s_mod <- P_S_cond_AX ~ A + X

fit <- svycdiff(samp, "OM", a_mod, s_mod, y_mod, "gaussian")

fit
```

## Vignette

Once you have `svycdiff` installed, you can type

``` r
vignette("svycdiff")
```

or

``` r
vignette("nhanes")
```

in `R` to bring up tutorials on `svycdiff` and how to use it.

## Questions

For questions and comments, please contact Stephen Salerno
(<ssalerno@fredhutch.org>).
