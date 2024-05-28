#===============================================================================
#
#  PROGRAM: sim.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: To simulate data based on the assumed relationship presented
#           in the conceptual diagrams for our simulation studies.
#
#  UPDATED: 2024-05-25
#
#  NOTES:
#
#    - Future Updates:
#
#        1. Allow Different Families for Y | A, X
#        2. Better Way to Allow for Known S | A, X
#
#===============================================================================

#=== SIMULATION FUNCTION =======================================================

#' Simulate data with varying degrees of selection and confounding bias
#'
#' @description
#'
#' Function to simulate data based on specified relationships between the
#' generated (continuous) outcome, variable of interest, confounder, and
#' selection mechanism.
#'
#' @details
#'
#' The data are generated as follows. For a user-given number, \code{N},
#' observations in our so-called super population, we first generate a
#' confounding variable, \code{X}, which relates to our outcome, \code{Y}, our
#' variable of interest, \code{A}, and our selection indicator, \code{S}.
#' We generate population-level data with \code{X ~ N(1,1)} or
#' \code{X ~ Bern(0.5)} depending on whether distribution of \code{X} is
#' chosen to be \code{X_dist = "continous"} or \code{X_dist = "binary"},
#' respectively.
#'
#' We then generate the remaining data from three models:
#'
#' \describe{
#'    \item{1. Propensity Model}{}
#'    \item{2. Selection Model}{}
#'    \item{3. Outcome Model}{}
#' }
#'
#' @param N int - Number of observations to be generated
#'
#' @param X_dist string - Distribution of the confounding variable, X. Defaults
#' to "continuous" for a N(1, 1) variable, or "binary" for a Bernoulli(0.5)
#' variable
#'
#' @param S_known boolean - Logical for whether the selection mechanism should
#' be treated as known (deterministic) or needs to be estimated (simulated with
#' Gaussian error; defaults to FALSE)
#'
#' @param tau_0 double - Intercept for propensity model (defaults to 0)
#'
#' @param tau_X double - Coefficient for X in propensity model (defaults to 1)
#'
#' @param beta_0 double - Intercept for selection model (defaults to 0)
#'
#' @param beta_A double - Coefficient for A in selection model (defaults to 1)
#'
#' @param beta_X double - Coefficient for X in selection model (defaults to 1)
#'
#' @param hetero boolean - Logical for heterogeneous treatment effect in
#' the outcome model (defaults to TRUE)
#'
#' @param alpha_0 double - Intercept for outcome model (defaults to 0)
#'
#' @param alpha_X double - Coefficient for X in outcome model (defaults to 1)
#'
#' @param alpha_A double - Coefficient for A in outcome model (defaults to 1)
#'
#' @param alpha_AX double - Coefficient for interaction between A and X in
#' outcome model (only used if \code{hetero == TRUE}; defaults to 0.1)
#'
#' @returns
#' A \code{data.frame} with \code{N} observations of 7 variables:
#' \describe{
#'    \item{Y}{Observed outcome (continuous)}
#'    \item{A}{Comparison group variable of interest (binary)}
#'    \item{X}{Confounding variable (continuous or binary)}
#'    \item{P_A_cond_X}{True probability of A = 1 conditional on X (continuous)}
#'    \item{P_S_cond_AX}{True probability of selection (S = 1) conditional on A and X (continuous)}
#'    \item{P_S_cond_A1X}{True probability of selection (S = 1) conditional on A = 1 and X (continuous)}
#'    \item{P_S_cond_A0X}{True probability of selection (S = 1) conditional on A = 0 and X (continuous)}
#'    \item{CDIFF}{True controlled difference in outcomes by comparison group (double)}
#' }
#'
#' @importFrom stats rbinom rnorm
#'
#' @examples
#'
#' N <- 100000
#'
#' dat <- simdat(N)
#'
#' head(dat)
#'
#' @export
simdat <- function(N, X_dist = "continuous", S_known = F,

  tau_0 = 0, tau_X = 1, beta_0 = 0, beta_A = 1, beta_X = 1,

  hetero = T, alpha_0 = 0, alpha_X = 1, alpha_A = 1, alpha_AX = 0.1) {

  #-- Covariate

  if (X_dist == "continuous") {

    X <- rnorm(N, 1, 1)

  } else if (X_dist == "binary") {

    X <- rbinom(N, 1, 0.5)

  } else {

    stop("`X_dist` must be either 'continuous' or 'binary'")
  }

  #-- Treatment Model

  P_A_cond_X <- plogis(tau_0 + tau_X*X)

  A <- rbinom(N, 1, P_A_cond_X)

  #-- Selection Model

  if (S_known) {

    P_S_cond_AX  <- plogis(beta_0 + beta_A*A + beta_X*X)

    P_S_cond_A1X <- plogis(beta_0 + beta_A*1 + beta_X*X)

    P_S_cond_A0X <- plogis(beta_0 + beta_A*0 + beta_X*X)

  } else {

    e_S <- rnorm(N, 0, 0.1)

    P_S_cond_AX  <- plogis(beta_0 + beta_A*A + beta_X*X + e_S)

    P_S_cond_A1X <- plogis(beta_0 + beta_A*1 + beta_X*X + e_S)

    P_S_cond_A0X <- plogis(beta_0 + beta_A*0 + beta_X*X + e_S)
  }

  #-- Outcome Model

  if (hetero) {

    aa <- 1; Y1 <- alpha_0 + alpha_A*aa + alpha_X*X + alpha_AX*aa*X + rnorm(N)

    aa <- 0; Y0 <- alpha_0 + alpha_A*aa + alpha_X*X + alpha_AX*aa*X + rnorm(N)

  } else {

    aa <- 1; Y1 <- alpha_0 + alpha_A*aa + alpha_X*X + rnorm(N)

    aa <- 0; Y0 <- alpha_0 + alpha_A*aa + alpha_X*X + rnorm(N)
  }

  Y <- A*Y1 + (1 - A)*Y0

  CDIFF <- mean(Y1 - Y0)

  #-- Final Data

  dat <- data.frame(Y, A, X,

    P_A_cond_X, P_S_cond_AX, P_S_cond_A1X, P_S_cond_A0X, CDIFF)

  return(dat)
}

#=== END =======================================================================
