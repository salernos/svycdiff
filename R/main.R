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
#' The argument \code{id_form} takes possible values "OM", "IPW1", or "IPW2",
#' corresponding to the three formulas presented in Salerno et al. "OM" refers
#' to the method that uses outcome modeling and direct standardization to
#' estimate the controlled difference, while IPW1 and IPW2 are inverse
#' probability weighted methods. IPW1 and IPW2 differ with respect to how the
#' joint propensity and selection mechanisms are factored (see Salerno et al.
#' for additional details). "OM", "IPW1", or "IPW2" are useful in different
#' settings, which warrants some brief discussion here. "OM" requires you to
#' specify an outcome regression model, whereas "IPW1" and "IPW2" do not
#' require estimation, nor do they assume additivity or interactivity. However,
#' while OM, IPW1, and IPW2 are consistent, OM is most efficient if correctly
#' specified.
#'
#' For IPW1/IPW2, \code{y_form} should be of the form \code{Y ~ 1}.
#'
#' For known S, \code{s_form} should be of the form \code{S ~ 1}, where
#' \code{S} is the variable corresponding to the probability of selection.
#' There should be two additional variables in the dataset: P_S_cond_A1X and
#' P_S_cond_A0X, corresponding to the known probability of selection
#' conditional on A = 1 or 0 and X = x, respectively. If these quantities are
#' not known, s_form should contain the variables which affect sample selection
#' on the right hand side of the equation, including the comparison group
#' variable of interest.
#'
#' @param df a data frame or tibble containing the variables in the models.
#'
#' @param id_form a string indicating which identification formula to be used.
#' Options include "OM", "IPW1", or "IPW2". See 'Details' for more information.
#'
#' @param a_form an object of class `formula` which describes the propensity
#' score model to be fit.
#'
#' @param s_form an object of class `formula` which describes the selection
#' model to be fit.
#'
#' @param y_form an object of class `formula` which describes the outcome model
#' to be fit. Only used if `id_form` = "OM", else `y_form = y ~ 1`.
#'
#' @param y_fam a family function. Only used if `id_form` = "OM", else
#' `y_fam = NULL`.
#'
#' @param strata a string indicating strata, else `strata = NULL` for no strata.
#'
#' @param cluster a string indicating cluster IDs, else `cluster = NULL` for no
#' clusters.
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
#' S <- rbinom(N, 1, dat$P_S_cond_AX)
#'
#' samp <- dat[S == 1,]
#'
#' y_mod <- Y ~ A * X
#'
#' a_mod <- A ~ X
#'
#' s_mod <- P_S_cond_AX ~ A + X
#'
#' fit <- svycdiff(samp, "OM", a_mod, s_mod, y_mod, "gaussian")
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
svycdiff <- function(df, id_form, a_form, s_form, y_form, y_fam = NULL,

  strata = NULL, cluster = NULL) {

  #--- SETUP -------------------------------------------------------------------

  #-- DATA

  yy <- df[, all.vars(y_form)[1]]
  aa <- df[, all.vars(a_form)[1]]
  ss <- df[, all.vars(s_form)[1]]

  df$ss <- ss

  df1 <- df; df1[[a_form[[2]]]] <- 1
  df0 <- df; df0[[a_form[[2]]]] <- 0

  S_known <- ifelse(s_form[[3]] == 1, T, F)

  #-- ASSERTIONS

  if (is.null(y_fam) & id_form == "OM") {

    stop("Must supply `y_fam` when `id_form` is 'OM'.")
  }

  if (y_form[[3]] == 1) {

    fit_y <- NULL
  }

  if (!is.numeric(aa)) {

    stop("`A` must be numeric.")
  }

  if (length(setdiff(0:1, unique(aa))) != 0) {

    stop("`A` must be coded 0/1.")
  }

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

    m0 <- predict(fit_y, newdata = df0)

    m1 <- predict(fit_y, newdata = df1)

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

  } else {

    stop("id_form must be one of `OM`, `IPW1`, or `IPW2`")
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
