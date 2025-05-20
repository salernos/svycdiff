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
#' @importFrom stats coef glm lm model.matrix qlogis plogis pnorm predict qt
#'
#' @export

svycdiff <- function(

  df,
  id_form,
  a_form,
  s_form,
  y_form,
  y_fam = NULL,
  strata = NULL,
  cluster = NULL) {

  #-- ASSERTIONS

  id_form <- match.arg(id_form, c("OM", "IPW1", "IPW2", "DR"))

  if (id_form %in% c("OM", "DR")) stopifnot(!is.null(y_fam))

  #-- HELPER FUNCTIONS

  prep <- prep_data(df, a_form, s_form, y_form)

  models <- fit_models(prep, id_form, a_form, s_form, y_form, y_fam)

  ests <- compute_targets(models, prep, id_form)

  u_func <- make_uclosure(prep, models, id_form, y_fam)

  vcov <- compute_vcov(u_func, ests, prep, strata, cluster)

  #-- FINAL OUTPUT

  est <- ests[length(ests)]

  se <- sqrt(vcov[nrow(vcov), ncol(vcov)])

  ci <- est + c(-1, 1) * qt(0.975, df = nrow(df) - 1) * se

  pval <- 2 * (1 - pnorm(abs(est) / se))

  structure(
    list(
      id_form   = id_form,
      cdiff     = c(est = est, err = se, lcl = ci[1], ucl = ci[2], pval = pval),
      fit_a     = models$fit_a,
      wtd_fit_a = models$wtd_fit_a,
      fit_s     = models$fit_s,
      fit_y     = models$fit_y
    ),

    class = "svycdiff"
  )
}
