#===============================================================================
#
#  PROGRAM: sim.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: To simulate data based on the assumed relationship presented
#           in the conceptual diagrams for our simulation studies.
#
#  UPDATED: 2025-05-01
#
#===============================================================================

#=== SIMULATION FUNCTION =======================================================

#' Simulate data with varying degrees of selection and confounding bias
#'
#' @description
#' Function to simulate data based on specified relationships between the
#' generated outcome, group variable, confounder(s), and selection mechanism.
#'
#' @details
#' The function generates data in a hierarchical structure with stratified
#' clusters. The data generation process follows these steps:
#'
#' 1. \strong{Stratum and Cluster Means:} For each of the \code{n_strat}
#'    strata, a matrix of stratum-level means for \code{p} covariates is
#'    generated from a normal distribution with standard deviation
#'    \code{sigma_strat}. Similarly, for each of the \code{n_clust} clusters
#'    within each stratum, cluster-level means are generated from a normal
#'    distribution with standard deviation \code{sigma_clust}.
#'
#' 2. \strong{Covariate Generation:} Within each cluster, covariates,
#'    \code{X}, for \code{N / (n_strat * n_clust)} individuals are generated
#'    from a multivariate normal distribution with mean equal to the sum of
#'    the cluster and stratum means, and an identity covariance matrix.
#'
#' 3. \strong{Covariate Transformation:} If \code{X_fam} is \code{"binary"},
#'    each covariate is discretized at its median, otherwise it remains
#'    continuous.
#'
#' 4. \strong{Propensity Model:} The group variable, \code{A}, is generated
#'    using a logistic regression model with intercept \code{tau_0}, covariate
#'    effects \code{tau_X}, and an interaction effect between the first two
#'    covariates with coefficient \code{tau_X12}. The group membership
#'    probability, \code{pA}, is defined by the logistic model.
#'
#' 5. \strong{Selection Model:} The probability of selection, \code{pS}, is
#'    generated using a logistic regression model with intercept \code{beta_0},
#'    group effect \code{beta_A}, and covariate effects \code{beta_X}. Gaussian
#'    noise is added to the linear predictor.
#'
#' 6. \strong{Outcome Model:} The outcome, \code{Y}, is generated based on a
#'    chosen outcome distribution, \code{Y_fam}. The linear predictor includes
#'    an intercept, \code{alpha_0}, group effect, \code{alpha_A}, covariate
#'    effects, \code{alpha_X}, and an optional interaction effect,
#'    \code{alpha_AX}, between the group variable and covariates.
#'
#' 7. \strong{Controlled Difference:} The true controlled difference in the
#'    outcome between groups is calculated as \code{CDIFF}.
#'
#' The output is a data frame containing the generated outcome, group variable,
#' covariates, and selection probabilities.
#'
#' @param N int - Number of observations to be generated. Defaults to 1000000.
#'
#' @param p int - Number of covariates to be generated. Defaults to 1.
#'
#' @param q int - Number of additional covariates that affect selection to be
#' generated. Defaults to 0.
#'
#' @param n_strat int - Number of strata in the population to be generated.
#' Defaults to 1.
#'
#' @param n_clust int - Number of clusters within each stratum in the
#' population to be generated. Defaults to 1.
#'
#' @param sigma_strat double - Standard deviation of covariate means across
#' strata. Defaults to 1.
#'
#' @param sigma_clust double - Standard deviation of covariate means across
#' clusters. Defaults to 1.
#'
#' @param X_fam string - Distribution of the covariates, \code{X}. Defaults
#' to a multivariate normal distribution with mean equal to the sum of the
#' cluster and stratum means, and an identity covariance matrix. If "binary",
#' continuous covariates are discretized at their median values.
#'
#' @param tau_0 double - Intercept for propensity model. Defaults to 0.
#'
#' @param tau_A double - Scaling factor for group assignment. Defaults to 1.
#'
#' @param tau_X double - Coefficients for \code{X} in propensity model.
#' Defaults to a 1 vector of length \code{p}.
#'
#' @param tau_X12 double - Interaction term coefficient for \code{X1*X2} if
#' p > 1. Defaults to 0.
#'
#' @param beta_0 double - Intercept for selection model. Defaults to 0.
#'
#' @param beta_A double - Coefficient for \code{A} in selection model.
#' Defaults to 1.
#'
#' @param beta_X double - Coefficients for \code{X} in selection model.
#' Defaults to a 1 vector of length \code{p}.
#'
#' @param beta_U double - Coefficients for \code{U} (additional covariates
#' affection only selection) in selection model. Defaults to a 1 vector of
#' length \code{q}.
#'
#' @param Y_fam string - Distribution of the outcome variable, \code{Y}.
#' Defaults to "gaussian" for a normally distributed outcome. Other options
#' include "binary" for a Bernoulli-distributed outcome and "poisson" for a
#' Poisson-distributed outcome.
#'
#' @param alpha_0 double - Intercept for outcome model. Defaults to 0.
#'
#' @param alpha_A double - Coefficient for \code{A} in outcome model.
#' Defaults to 1.
#'
#' @param alpha_X double - Coefficients for \code{X} in outcome model.
#' Defaults to a 1 vector of length \code{p}.
#'
#' @param alpha_AX double - Coefficient for interaction between \code{A} and
#' \code{X} in outcome model. Defaults to 0.
#'
#' @returns
#' A \code{data.frame} with \code{N} observations and the following variables:
#' \describe{
#'    \item{Strata}{Stratum index (integer)}
#'    \item{Cluster}{Cluster index (integer)}
#'    \item{X1, X2, ..., Xp}{Confounding covariates (continuous or binary,
#'    depending on \code{X_fam})}
#'    \item{pA}{True probability of A = 1 conditional on X (continuous)}
#'    \item{A}{Group assignment (binary)}
#'    \item{pS}{True probability of selection conditional on A and X
#'    (continuous)}
#'    \item{Y0}{Potential outcome under A = 0 (continuous, binary, or count
#'    depending on \code{Y_fam})}
#'    \item{Y1}{Potential outcome under A = 1 (continuous, binary, or count
#'    depending on \code{Y_fam})}
#'    \item{Y}{Observed outcome, based on treatment assignment (continuous,
#'    binary, or count depending on \code{Y_fam})}
#'    \item{CDIFF}{True controlled difference in outcomes by comparison group
#'    (double, computed as mean(Y1 - Y0))}
#' }
#'
#' @importFrom stats median rbinom rnorm rpois
#'
#' @importFrom MASS mvrnorm
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
simdat <- function(N = 1000000,
                   p = 1,
                   q = 0,
                   #- Design Variables
                   n_strat = 1,
                   n_clust = 1,
                   sigma_strat = 1,
                   sigma_clust = 1,
                   X_fam = c("gaussian", "binary"),
                   #- Propensity Model
                   tau_0 = 0,
                   tau_A = 1,
                   tau_X = rep(1, p),
                   tau_X12 = 0,
                   #- Selection Model
                   beta_0 = 0,
                   beta_A = 1,
                   beta_X = rep(1, p),
                   beta_U = rep(1, q),
                   #- Outcome Model
                   Y_fam = c("gaussian", "binary", "poisson"),
                   alpha_0 = 0,
                   alpha_A = 1,
                   alpha_X = rep(1, p),
                   alpha_AX = 0) {

  X_fam <- match.arg(X_fam)
  Y_fam <- match.arg(Y_fam)

  n_grp <- N / (n_strat * n_clust)

  # Generate Strata and Cluster Means

  mu_strat <- matrix(rnorm(n_strat * (p + q), 0, sigma_strat), n_strat, (p + q))
  mu_clust <- matrix(rnorm(n_clust * (p + q), 0, sigma_clust), n_clust, (p + q))

  # Generate Covariates

  dat <- do.call(rbind, lapply(1:n_strat, function(s) {

    do.call(rbind, lapply(1:n_clust, function(c) {

      mu <- mu_strat[s, ] + mu_clust[c, ]

      cbind(
        Cluster = (s - 1) * n_clust + c,
        Strata = s,
        mvrnorm(n_grp, mu, diag(rep(1, p + q)))
      )
    }))
  }))

  if (X_fam == "binary") {

    dat[, 3:(2 + p)] <- apply(

      dat[, 3:(2 + p)], 2,

      function(x) {

        as.integer(x > median(x))
      }
    )
  }

  dat <- data.frame(dat)

  if (q > 0) {

    colnames(dat) <- c("Cluster", "Strata", paste0("X", 1:p), paste0("U", 1:q))

    X <- as.matrix(dat[, 3:(2 + p)])

    U <- as.matrix(dat[, (3 + p):(2 + p + q)])

  } else {

    colnames(dat) <- c("Cluster", "Strata", paste0("X", 1:p))

    X <- as.matrix(dat[, 3:(2 + p)])
  }

  # Propensity Model

  if (p > 1) {

    dat$pA <- plogis(tau_0 + tau_A * (X %*% tau_X + tau_X12 * X[, 1] * X[, 2]))

  } else {

    dat$pA <- plogis(tau_0 + tau_A * (X %*% tau_X))
  }

  dat$A <- rbinom(N, 1, dat$pA)

  # Selection Model

  eS <- rnorm(N, 0, 0.1)

  dat$pS <- plogis(beta_0 + beta_A * dat$A + X %*% beta_X + U %*% beta_U  + eS)

  # Outcome Model

  mu_Y0 <- as.vector(alpha_0 + X %*% alpha_X)
  mu_Y1 <- as.vector(alpha_0 + alpha_A + alpha_AX * (X %*% alpha_X))

  if (Y_fam == "gaussian") {

    dat$Y0 <- mu_Y0 + rnorm(N)
    dat$Y1 <- mu_Y1 + rnorm(N)

  } else if (Y_fam == "binary") {

    dat$Y0 <- rbinom(N, 1, plogis(mu_Y0))
    dat$Y1 <- rbinom(N, 1, plogis(mu_Y1))

  } else if (Y_fam == "poisson") {

    dat$Y0 <- rpois(N, lambda = exp(mu_Y0))
    dat$Y1 <- rpois(N, lambda = exp(mu_Y1))

  } else {

    stop("`Y_fam` must be 'gaussian', 'binary', or 'poisson'")
  }

  dat$Y <- dat$A * dat$Y1 + (1 - dat$A) * dat$Y0

  # Calculate Controlled Difference

  dat$CDIFF <- mean(dat$Y1 - dat$Y0)

  return(dat)
}

#=== END =======================================================================
