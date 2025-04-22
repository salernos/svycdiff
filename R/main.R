#===============================================================================
#
#  PROGRAM: main.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Implements the main function (svycdiff) for estimating the
#           population average controlled difference (ACD) or population
#           average treatment effect (ATE) in complex surveys where
#           selection depends on the comparison group of interest.
#
#  UPDATED: 2025-04-22
#
#===============================================================================

#=== MAIN FUNCTION =============================================================

#' Controlled Difference Estimation for Complex Surveys
#'
#' @description
#'
#' This is the main function to estimate population average controlled
#' difference (ACD), or under stronger assumptions, the population average
#' treatment effect (PATE), for a given outcome between levels of a binary
#' treatment, exposure, or other group membership variable of interest for
#' clustered, stratified survey samples where sample selection depends on
#' the comparison group.
#'
#' @details
#'
#' The argument \code{id_form} takes possible values \code{"OM"},
#' \code{"IPW1"}, \code{"IPW2"}, or \code{"DR"}, corresponding to the four
#' formulas presented in Salerno et al. \code{"OM"} refers to the method that
#' uses outcome modeling and direct standardization to estimate the controlled
#' difference, while \code{"IPW1"} and \code{"IPW2"} are inverse probability
#' weighted methods. \code{"IPW1"} and \code{"IPW2"} differ with respect to how
#' the joint propensity and selection mechanisms are factored (see Salerno et
#' al. for additional details). \code{"DR"} refers to the doubly robust form of
#' estimator, which essentially combines \code{"OM"} and \code{"IPW2"}.
#'
#' For \code{id_form = "IPW1"} or \code{id_form = "IPW2"}, \code{y_form} should
#' be of the form \code{Y ~ 1}.
#'
#' For known selection mechanisms, \code{s_form} should be of the form
#' \code{pS ~ 1}, where \code{pS} is the variable corresponding to the
#' probability of selection (e.g., inverse of the selection weight), and there
#' should be two additional variables in the dataset: \code{P_S_cond_A1X} and
#' \code{P_S_cond_A0X}, corresponding to the known probability of selection
#' conditional on \eqn{A = 1} or \eqn{0} and \eqn{X = x}, respectively. If
#' these quantities are not known, \code{s_form} should contain the variables
#' which affect sample selection on the right hand side of the equation,
#' including the comparison group variable of interest.
#'
#' @param df a `data.frame` or `tibble` containing the variables in the models.
#'
#' @param id_form a `string` indicating which identification formula to be used.
#' Options include \code{"OM"}, \code{"IPW1"}, \code{"IPW2"}, or \code{"DR"}.
#' See 'Details' for information.
#'
#' @param a_form an object of class `formula` which describes the propensity
#' score model to be fit.
#'
#' @param s_form an object of class `formula` which describes the selection
#' model to be fit.
#'
#' @param y_form an object of class `formula` which describes the outcome model
#' to be fit. Only used if \code{id_form = "OM"}, else \code{y_form = y ~ 1}.
#'
#' @param y_fam a `family` function. Only used if \code{id_form = "OM"}, else
#' \code{y_fam = NULL}. Current options include \code{gaussian},
#' \code{binomial}, or \code{poisson}.
#'
#' @param strata a `string` indicating strata, else \code{strata = NULL}
#' for no strata.
#'
#' @param cluster a `string` indicating cluster IDs, else \code{cluster = NULL}
#' for no clusters.
#'
#' @returns `svycdiff` returns an object of class "svycdiff" which contains:
#'
#' \describe{
#'    \item{id_form}{A string denoting Which method was selected for estimation}
#'    \item{cdiff}{A named vector containing the point estimate (est),
#'                 standard error (err), lower confidence limit (lcl),
#'                 upper confidence limit (ucl), and p-value (pval) for
#'                 the estimated controlled difference}
#'    \item{fit_y}{An object of class inheriting from "glm" corresponding to
#'                 the outcome model fit, or NULL for IPW1 and IPW2}
#'    \item{fit_a}{An object of class inheriting from "glm" corresponding to
#'                 the propensity model fit}
#'    \item{wtd_fit_a}{An object of class inheriting from "glm" corresponding
#'                     to the weighted propensity model fit}
#'    \item{fit_s}{An object of class "betareg" corresponding to the selection
#'                 model fit, or NULL if the selection mechanism is known}
#' }
#'
#' @examples
#'
#' N <- 1000
#'
#' dat <- simdat(N)
#'
#' S <- rbinom(N, 1, dat$pS)
#'
#' samp <- dat[S == 1,]
#'
#' y_mod <- Y ~ A * X1
#'
#' a_mod <- A ~ X1
#'
#' s_mod <- pS ~ A + X1
#'
#' fit <- svycdiff(samp, "DR", a_mod, s_mod, y_mod, "gaussian")
#'
#' fit
#'
#' summary(fit)
#'
#' @importFrom betareg betareg
#'
#' @importFrom numDeriv jacobian
#'
#' @importFrom survey svydesign svyglm
#'
#' @importFrom stats coef glm lm model.matrix qlogis plogis pnorm predict qnorm
#'
#' @export
svycdiff <- function(df,
                     id_form,
                     a_form,
                     s_form,
                     y_form,
                     y_fam = NULL,
                     strata = NULL,
                     cluster = NULL) {

  #--- SETUP -------------------------------------------------------------------

  #-- DATA

  yy <- df[, all.vars(y_form)[1]]
  aa <- df[, all.vars(a_form)[1]]
  ss <- df[, all.vars(s_form)[1]]

  df$ss <- ss

  df1 <- df; df1[[a_form[[2]]]] <- 1
  df0 <- df; df0[[a_form[[2]]]] <- 0

  S_known <- ifelse(s_form[[3]] == 1, TRUE, FALSE)

  #-- ASSERTIONS

  if (id_form == "OM" || id_form == "DR") {

    if (is.null(y_fam)) stop("Must supply `y_fam` for OM or DR.")
  }

  if (y_form[[3]] == 1) {

    fit_y <- NULL
  }

  if (!all(aa %in% c(0,1))) stop("`A` must be coded 0/1.")

  #--- FIT COMMON MODELS -------------------------------------------------------

  #-- Pr(A = 1 | X, S = 1): Within-Sample Propensity Score

  fit_a <- glm(a_form, family = "quasibinomial", data = df)

  P_A1_cond_X <- predict(fit_a, newdata = df, type = "response")

  #-- Pr(A = 1 | X): Weighted Propensity Score

  des_wtd_a <- svydesign(ids = ~ 1, weights = 1 / df$ss, data = df)

  wtd_fit_a <- svyglm(a_form, des_wtd_a, family = "quasibinomial")

  wtd_P_A1_cond_X <- predict(wtd_fit_a, newdata = df, type = "response")

  #-- Pr(S = 1 | A, X): Probability of Selection

  if (S_known) {

    fit_s <- NULL

    #- Pr(S = 1 | A = 1, X): Probability of Selection Among A = 1

    P_S1_cond_A1X <- df$P_S_cond_A1X

    #- Pr(S = 1 | A = 0, X): Probability of Selection Among A = 0

    P_S1_cond_A0X <- df$P_S_cond_A0X

  } else {

    fit_s <- betareg(s_form, data = df)

    #- Pr(S = 1 | A = 1, X): Probability of Selection Among A = 1

    P_S1_cond_A1X <- predict(fit_s, newdata = df1, type = "response")

    #- Pr(S = 1 | A = 0, X): Probability of Selection Among A = 0

    P_S1_cond_A0X <- predict(fit_s, newdata = df0, type = "response")
  }

  #-- Pr(S = 1 | X): Probability of Selection, Marginalized Over A

  P_S1_cond_X <- P_S1_cond_A1X*wtd_P_A1_cond_X +

    P_S1_cond_A0X*(1 - wtd_P_A1_cond_X)

  #- Pr(S = 1): Marginal Sampling Fraction

  P_S1 <- length(ss) / sum(1/ss)

  #--- METHOD-SPECIFIC CALCULATIONS --------------------------------------------

  if (id_form == "OM"){

    fit_y <- glm(y_form, family = y_fam, data = df)

    m0 <- predict(fit_y, newdata = df0, type = "response")

    m1 <- predict(fit_y, newdata = df1, type = "response")

    est <- mean((m1 - m0)/P_S1_cond_X)*P_S1

    if (S_known) {

      ests <- c(coef(wtd_fit_a), coef(fit_y), cdiff = est)

      ests_dims <- c(length(coef(wtd_fit_a)), length(coef(fit_y)), 1)

    } else {

      ests <- c(coef(wtd_fit_a), coef(fit_s), coef(fit_y), cdiff = est)

      ests_dims <- c(

        length(coef(wtd_fit_a)), length(coef(fit_s)), length(coef(fit_y)), 1)
    }

  } else if (id_form == "IPW1") {

    est <- mean((aa*yy/wtd_P_A1_cond_X/P_S1_cond_A1X) -

      ((1 - aa)*yy/(1 - wtd_P_A1_cond_X)/P_S1_cond_A0X))*P_S1

    if (S_known) {

      ests <- c(coef(wtd_fit_a), cdiff = est)

      ests_dims <- c(length(coef(wtd_fit_a)), 1)

    } else {

      ests <- c(coef(wtd_fit_a), coef(fit_s), cdiff = est)

      ests_dims <- c(length(coef(wtd_fit_a)), length(coef(fit_s)), 1)
    }

  } else if (id_form == "IPW2") {

    est <- mean(((aa*yy/P_A1_cond_X) - ((1 - aa)*yy/(1 - P_A1_cond_X))) /

      P_S1_cond_X)*P_S1

    if (S_known) {

      ests <- c(coef(fit_a), coef(wtd_fit_a), cdiff = est)

      ests_dims <- c(length(coef(fit_a)), length(coef(wtd_fit_a)), 1)

    } else {

      ests <- c(coef(fit_a), coef(wtd_fit_a), coef(fit_s), cdiff = est)

      ests_dims <- c(

        length(coef(fit_a)), length(coef(wtd_fit_a)), length(coef(fit_s)), 1)
    }

  } else if (id_form == "DR") {

    fit_y <- glm(y_form, family = y_fam, data = df)

    m1 <- predict(fit_y, newdata = df1, type = "response")

    m0 <- predict(fit_y, newdata = df0, type = "response")

    mu1 <- (m1 + aa * (yy - m1) / P_A1_cond_X)

    mu0 <- (m0 + (1 - aa) * (yy - m0) / (1 - P_A1_cond_X))

    est <- mean((mu1 - mu0)/P_S1_cond_X)*P_S1

    if (S_known) {

      ests <- c(coef(fit_a), coef(wtd_fit_a), coef(fit_y), cdiff = est)

      ests_dims <- c(length(coef(wtd_fit_a)), length(coef(fit_y)), 1)

    } else {

      ests <- c(coef(fit_a), coef(wtd_fit_a), coef(fit_s),

        coef(fit_y), cdiff = est)

      ests_dims <- c(length(coef(fit_a)), length(coef(wtd_fit_a)),

        length(coef(fit_s)), length(coef(fit_y)), 1)
    }

  } else {

    stop("id_form must be one of 'OM', 'IPW1', 'IPW2', or 'DR'")
  }

  #--- INFERENCE ---------------------------------------------------------------

  halfmeat <- .U(

    ests, ests_dims, df, id_form, a_form, s_form, y_form, y_fam, S_known)

  bread <- jacobian(.G, ests,

    theta_dims = ests_dims, df = df, id_form = id_form, a_form = a_form,

    s_form = s_form, y_form = y_form, y_fam = y_fam, S_known = S_known)

  IF <- tcrossprod(halfmeat, solve(-bread))

  if (!is.null(strata)) {

    err <- matrix(0, nrow = ncol(IF), ncol = ncol(IF))

    stratas <- df[, strata]

    strata_levels <- unique(stratas)

    if (!is.null(cluster)) {

      for (k in 1:length(strata_levels)) {

        df_k <- df[which(stratas == strata_levels[k]), ]

        IF_k <- IF[which(stratas == strata_levels[k]), ]

        clusters_k <- df_k[, cluster]

        clusters_k_levels <- unique(clusters_k)

        err_gk <- matrix(0, nrow = ncol(IF), ncol = ncol(IF))

        tot_IF_g <- vector("list", length(clusters_k_levels))

        for(g in 1:length(clusters_k_levels)) {

          IF_g <- IF_k[which(clusters_k == clusters_k_levels[g]), ]

          tot_IF_g[[g]] <- colSums(IF_g)
        }

        avg_IF_k <- Reduce('+', tot_IF_g) / length(clusters_k_levels)

        err_gk <- Reduce('+', lapply(tot_IF_g, function(g) {

          tcrossprod(g - avg_IF_k)

        }))

        err <- err +

          (length(clusters_k_levels) / (length(clusters_k_levels) - 1))*err_gk
      }

      err <- sqrt(err[nrow(err), ncol(err)])

    } else {

      # NEED
    }

  } else {

    err <- sqrt(sum(IF[, ncol(IF)]^2))

  }

  lcl <- est - qnorm(0.975)*err

  ucl <- est + qnorm(0.975)*err

  pval <- 2*(1 - pnorm(abs(est)/err))

  out <- c(est = est, err = err, lcl = lcl, ucl = ucl, pval = pval)

  return(

    structure(

      list(id_form = id_form, cdiff = out, fit_y = fit_y, fit_a = fit_a,

        wtd_fit_a = wtd_fit_a, fit_s = fit_s),

      class = "svycdiff")
  )
}

#=== END =======================================================================
