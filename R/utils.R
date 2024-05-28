#===============================================================================
#
#  PROGRAM: utils.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Implements helper functions and methods for svycdiff (see main.R):
#
#             .U      - Helper function for calculating score equations (hidden)
#             .G      - Helper function for calculating Jacobian (hidden)
#             print   - Print method
#             summary - Summary method
#
#  UPDATED: 2024-05-24
#
#  NOTES:
#
#    - Future Updates:
#
#        1. Allow Different Families for Y | A, X
#        2. Better Way to Allow for Known S | A, X
#
#===============================================================================

#=== HELPER FUNCTIONS ==========================================================

#--- U FUNCTION ----------------------------------------------------------------

.U <- function(

  theta, theta_dims, df, id_form, a_form, s_form, y_form, y_fam, S_known) {

  #--- DATA

  yy <- df[, all.vars(y_form)[1]]
  aa <- df[, all.vars(a_form)[1]]
  ss <- df[, all.vars(s_form)[1]]

  df1 <- df; df1[[a_form[[2]]]] <- 1
  df0 <- df; df0[[a_form[[2]]]] <- 0

  #--- DESIGN MATRICES FOR INFERENCE

  A_X   <- model.matrix(df,  object = a_form)

  S_AX  <- model.matrix(df,  object = s_form)

  S_A1X <- model.matrix(df1, object = s_form)
  S_A0X <- model.matrix(df0, object = s_form)

  Y_AX  <- model.matrix(df,  object = y_form)

  Y_A1X <- model.matrix(df1, object = y_form)
  Y_A0X <- model.matrix(df0, object = y_form)

  #--- OM

  if (id_form == "OM") {

    #-- PARAMETERS

    tau_w <- theta[1:theta_dims[1]]

    if (S_known) {

      gamma <- theta[theta_dims[1] + 1:theta_dims[2]]

      CDIFF   <- theta[theta_dims[1] + theta_dims[2] + 1]

    } else {

      beta  <- theta[theta_dims[1] + 1:(theta_dims[2] - 1)]

      phi   <- theta[theta_dims[1] + theta_dims[2]]

      gamma <- theta[theta_dims[1] + theta_dims[2] + 1:theta_dims[3]]

      CDIFF   <- theta[theta_dims[1] + theta_dims[2] + theta_dims[3] + 1]
    }

    #-- SCORE EQUATIONS

    #- Propensity Model

    wtd_P_A1_cond_X <- plogis(A_X%*%tau_w)

    wtd_U_A <- A_X*c(aa - wtd_P_A1_cond_X)/ss

    #- Selection Model

    if (S_known) {

      P_S1_cond_A1X <- df$P_S_cond_A1X

      P_S1_cond_A0X <- df$P_S_cond_A0X

    } else {

      P_S1_cond_AX  <- plogis(S_AX%*%beta)

      P_S1_cond_A1X <- plogis(S_A1X%*%beta)

      P_S1_cond_A0X <- plogis(S_A0X%*%beta)

      ss_star <- qlogis(ss)

      mu_star <- digamma(P_S1_cond_AX*phi) - digamma((1 - P_S1_cond_AX)*phi)

      W <- diag(c(P_S1_cond_AX*(1 - P_S1_cond_AX)))

      U_S <- phi*(W%*%S_AX)*c(ss_star - mu_star)

      U_Phi <- (P_S1_cond_AX * (ss_star - mu_star) + log(1 - ss) -

        digamma((1 - P_S1_cond_AX)*phi) + digamma(phi))
    }

    P_S1_cond_X <- P_S1_cond_A1X*wtd_P_A1_cond_X +

      P_S1_cond_A0X*(1 - wtd_P_A1_cond_X)

    P_S1 <- length(ss) / sum(1/ss)

    #- Outcome Model

    U_Y <- Y_AX*c(yy - Y_AX%*%gamma)

    #- CDIFF

    m1 <- Y_A1X%*%gamma
    m0 <- Y_A0X%*%gamma

    U_CDIFF <- CDIFF - ((m1 - m0)/P_S1_cond_X)*P_S1

    #-- STACK ESTIMATING EQUATIONS

    if (S_known) {

      return(cbind(wtd_U_A, U_Y, U_CDIFF))

    } else {

      return(cbind(wtd_U_A, U_S, U_Phi, U_Y, U_CDIFF))
    }

  #--- IPW1

  } else if (id_form == "IPW1") {

    #-- PARAMETERS

    tau_w <- theta[1:theta_dims[1]]

    if (S_known) {

      CDIFF <- theta[theta_dims[1] + 1]

    } else {

      beta <- theta[theta_dims[1] + 1:(theta_dims[2] - 1)]

      phi  <- theta[theta_dims[1] + theta_dims[2]]

      CDIFF  <- theta[theta_dims[1] + theta_dims[2] + 1]
    }

    #-- SCORE EQUATIONS

    #- Propensity Model

    wtd_P_A1_cond_X <- plogis(A_X%*%tau_w)

    wtd_U_A <- A_X*c(aa - wtd_P_A1_cond_X)/ss

    #- Selection Model

    if (S_known) {

      P_S1_cond_A1X <- df$P_S_cond_A1X

      P_S1_cond_A0X <- df$P_S_cond_A0X

    } else {

      P_S1_cond_AX  <- plogis(S_AX%*%beta)

      P_S1_cond_A1X <- plogis(S_A1X%*%beta)

      P_S1_cond_A0X <- plogis(S_A0X%*%beta)

      ss_star <- qlogis(ss)

      mu_star <- digamma(P_S1_cond_AX*phi) - digamma((1 - P_S1_cond_AX)*phi)

      W <- diag(c(P_S1_cond_AX*(1 - P_S1_cond_AX)))

      U_S <- phi*(W%*%S_AX)*c(ss_star - mu_star)

      U_Phi <- (P_S1_cond_AX * (ss_star - mu_star) + log(1 - ss) -

        digamma((1 - P_S1_cond_AX)*phi) + digamma(phi))
    }

    P_S1 <- length(ss) / sum(1/ss)

    #- CDIFF

    U_CDIFF <- CDIFF - ((aa*yy/wtd_P_A1_cond_X/P_S1_cond_A1X) -

      ((1 - aa)*yy/(1 - wtd_P_A1_cond_X)/P_S1_cond_A0X))*P_S1

    #-- STACK ESTIMATING EQUATIONS

    if (S_known) {

      return(cbind(wtd_U_A, U_CDIFF))

    } else {

      return(cbind(wtd_U_A, U_S, U_Phi, U_CDIFF))
    }

  #--- IPW2

  } else if (id_form == "IPW2") {

    #-- PARAMETERS

    tau   <- theta[1:theta_dims[1]]

    tau_w <- theta[theta_dims[1] + 1:theta_dims[2]]

    if (S_known) {

      CDIFF  <- theta[theta_dims[1] + theta_dims[2] + 1]

    } else {

      beta <- theta[theta_dims[1] + theta_dims[2] + 1:(theta_dims[3] - 1)]

      phi  <- theta[theta_dims[1] + theta_dims[2] + theta_dims[3]]

      CDIFF  <- theta[theta_dims[1] + theta_dims[2] + theta_dims[3] + 1]
    }

    #-- SCORE EQUATIONS

    #- Propensity Models

    P_A1_cond_X <- plogis(A_X%*%tau)

    U_A <- A_X*c(aa - P_A1_cond_X)

    wtd_P_A1_cond_X <- plogis(A_X%*%tau_w)

    wtd_U_A <- A_X*c(aa - wtd_P_A1_cond_X)/ss

    #- Selection Model

    if (S_known) {

      P_S1_cond_A1X <- df$P_S_cond_A1X

      P_S1_cond_A0X <- df$P_S_cond_A0X

    } else {

      P_S1_cond_AX  <- plogis(S_AX%*%beta)

      P_S1_cond_A1X <- plogis(S_A1X%*%beta)

      P_S1_cond_A0X <- plogis(S_A0X%*%beta)

      ss_star <- qlogis(ss)

      mu_star <- digamma(P_S1_cond_AX*phi) - digamma((1 - P_S1_cond_AX)*phi)

      W <- diag(c(P_S1_cond_AX*(1 - P_S1_cond_AX)))

      U_S <- phi*(W%*%S_AX)*c(ss_star - mu_star)

      U_Phi <- (P_S1_cond_AX * (ss_star - mu_star) + log(1 - ss) -

        digamma((1 - P_S1_cond_AX)*phi) + digamma(phi))
    }

    P_S1_cond_X <- P_S1_cond_A1X*wtd_P_A1_cond_X +

      P_S1_cond_A0X*(1 - wtd_P_A1_cond_X)

    P_S1 <- length(ss) / sum(1/ss)

    #- CDIFF

    U_CDIFF <- CDIFF - ((aa*yy/P_A1_cond_X) -

      ((1 - aa)*yy/(1 - P_A1_cond_X))) / P_S1_cond_X*P_S1

    #-- STACK ESTIMATING EQUATIONS

    if (S_known) {

      return(cbind(U_A, wtd_U_A, U_CDIFF))

    } else {

      return(cbind(U_A, wtd_U_A, U_S, U_Phi, U_CDIFF))
    }

  } else {

    stop("'id_form' must be one of 'OM', 'IPW1', or 'IPW2'")
  }
}

#--- G FUNCTION ----------------------------------------------------------------

.G <- function(

  theta, theta_dims, df, id_form, a_form, s_form, y_form, y_fam, S_known){

  return(apply(

    .U(theta, theta_dims, df, id_form, a_form, s_form, y_form, y_fam, S_known),

    2, sum))
}

#=== METHODS ===================================================================

#--- PRINT ---------------------------------------------------------------------

#' @export
print.svycdiff <- function(x, digits = 4L, ...) {

  cat(

    "\nOutcome Model:  \n",

    paste(deparse(x$fit_y$call), sep = "\n", collapse = "\n"), "\n",

    "\nTreatment Model:  \n",

    paste(deparse(x$fit_a$call), sep = "\n", collapse = "\n"), "\n",

    "\nSelection Model:  \n",

    paste(deparse(x$fit_s$call), sep = "\n", collapse = "\n"), "\n", sep = "")

  cat("\nCDIFF:  \n")

  cdiff <- x$cdiff

  names(cdiff) <- c("CDIFF", "SE", "LCL", "UCL", "P-Value")

  print(round(cdiff, digits))
}

#--- SUMMARY -------------------------------------------------------------------

#' @export
summary.svycdiff <- function(object, digits = 4L, ...) {

  cat("\nCDIFF:  \n")

  cdiff <- object$cdiff

  names(cdiff) <- c("CDIFF", "SE", "LCL", "UCL", "P-Value")

  print(round(cdiff, digits))

  cat("\nOutcome Model:  \n")

  summary(object$fit_y)

  cat("\nTreatment Model:  \n")

  summary(object$fit_a)

  cat("\nSelection Model:  \n")

  summary(object$fit_s)
}

#=== END =======================================================================
