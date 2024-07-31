#===============================================================================
#
#  PROGRAM: SIMULATION.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: To quantify the bias and variance associated with the use of
#           various commonly-used weighting schemes in propensity score
#           construction within complex survey design settings against
#           our proposed methods.
#
#  INUPT:   main.R - R script containing main function for the proposed methods
#
#           utils.R - R script containing helper functions for methods
#
#           sim.R - R script containing functions to simulate data
#
#  OUTPUT:  SIMULATION_RESULTS.RData
#
#           An .RData file containing the object `results`, which is a list of
#           eight elements corresponding to the eight simulation settings
#           described in the manuscript. Each list element contains a numeric
#           matrix of size `nsims` x `37`, where each row corresponds to a
#           simulation replicate and each column is a saved metric.
#
#  DETAILS: Methods Being Compared:
#
#             00. Oracle Estimator
#             01. Simple Regression
#             02. Multiple Regression
#             03. IPTW Estimator
#             04. Survey-Weighted Multiple Regression
#             05. IPTW Multiple Regression
#             06. IPTW + Survey Weighted Multiple Regression
#             07. Weighted IPTW + Survey Weighted Multiple Regression
#             08. Naive G-Computation
#             09. Formula 1 - Outcome Modeling and Direct Standardization
#             10. Formula 2a - Inverse Probability Weighting 1
#             11. Formula 2b - Inverse Probability Weighting 2
#
#  UPDATED: 2024-05-25
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- SOURCE NECESSARY PACKAGES -------------------------------------------------

library(pacman)

p_load(progressr, doFuture, numDeriv, betareg, survey, gt, here, update = FALSE)

#--- SOURCE FUNCTIONS ----------------------------------------------------------

source(here("R", "main.R"))

source(here("R", "utils.R"))

source(here("R", "sim.R"))

#--- HELPER FUNCTIONS ----------------------------------------------------------

het_est <- function(dat, fit, weight) {

  smy <- coef(fit)

  est <- with(dat, smy[2] + mean(X*weight)*smy[4])

  return(est)
}

het_err <- function(dat, fit, weight) {

  vcv <- vcov(summary(fit))

  err <- with(dat,

    sqrt(vcv[2,2] + (mean(X*weight)^2)*vcv[4,4] + 2*mean(X*weight)*vcv[2,4]))

  return(err)
}

U_Oracle <- function(theta, U_data) {

  #--- DATA

  yy <- U_data$yy
  aa <- U_data$aa
  xx <- U_data$xx
  ss <- U_data$ss

  S <- U_data$S

  dat <- data.frame(aa, yy, xx)

  dat1 <- data.frame(aa = 1, yy, xx)
  dat0 <- data.frame(aa = 0, yy, xx)

  #--- PARAMETERS

  gamma <- theta[1:4]
  ATE   <- theta[5]

  #--- SCORE EQUATIONS

  #- Outcome Model

  OR_f <- as.formula("~aa*xx")

  Y_AX <- model.matrix(dat, object = OR_f)

  m1 <- model.matrix(dat1, object = OR_f)%*%gamma
  m0 <- model.matrix(dat0, object = OR_f)%*%gamma

  U_Y <- Y_AX*c(yy - Y_AX%*%gamma)*S

  #- ATE

  U_ATE <- ATE - (m1 - m0)

  return(cbind(U_Y, U_ATE))
}

U_IPTW <- function(theta, U_data) {

  #--- DATA

  yy <- U_data$yy
  aa <- U_data$aa
  xx <- U_data$xx
  ss <- U_data$ss

  dat <- data.frame(aa, yy, xx)

  #--- PARAMETERS

  tau <- theta[1:2]
  ATE <- theta[3]

  #--- SCORE EQUATIONS

  #- Propensity Model

  A_f <- as.formula("~xx")

  A_X <- model.matrix(dat, object = A_f)

  P_A1_cond_X <- plogis(A_X%*%tau)

  U_A <- A_X*c(aa - P_A1_cond_X)

  #- ATE

  U_ATE <- ATE - ((aa*yy/P_A1_cond_X) - ((1 - aa)*yy/(1 - P_A1_cond_X)))

  return(cbind(U_A, U_ATE))
}

U_naiveG <- function(theta, U_data) {

  #--- DATA

  yy <- U_data$yy
  aa <- U_data$aa
  xx <- U_data$xx
  ss <- U_data$ss

  dat <- data.frame(aa, yy, xx)

  dat1 <- data.frame(aa = 1, yy, xx)
  dat0 <- data.frame(aa = 0, yy, xx)

  #--- PARAMETERS

  gamma <- theta[1:4]
  ATE   <- theta[5]

  #--- SCORE EQUATIONS

  OR_f <- as.formula("~aa*xx")

  Y_AX <- model.matrix(dat, object = OR_f)

  m1 <- model.matrix(dat1, object = OR_f)%*%gamma
  m0 <- model.matrix(dat0, object = OR_f)%*%gamma

  U_Y <- Y_AX*c(yy - Y_AX%*%gamma)

  U_ATE <- ATE - (m1 - m0)*(1/ss)/sum(1/ss)

  return(cbind(U_Y, U_ATE))
}

G <- function(theta, U_data, U_func) {

  U <- get(U_func)

  return(apply(U(theta, U_data), 2, sum))
}

#=== GLOBAL VARIABLES ==========================================================

tau_X_vec  <- 0:1

beta_A_vec <- 0:1

beta_X_vec <- 0:1

settings <- expand.grid(tau_X_vec, beta_A_vec, beta_X_vec)

nsims <- 200

N <- 100000

#=== SIMULATION ================================================================

plan(multisession)

results <- vector("list", nrow(settings))

tic <- proc.time()

for(i in 1:length(results)) {

  cat("\nSetting", i, "---\n\n")

  #--- GENERATE SUPER POPULATION -----------------------------------------------

  set.seed(i)

  tau_X  <- settings[i, 1]

  beta_A <- settings[i, 2]

  beta_X <- settings[i, 3]

  dat <- simdat(N, X_dist = "continuous", S_known = FALSE,

    tau_0 = -1, tau_X = tau_X,

    beta_0 = -4.5, beta_A = beta_A, beta_X = beta_X,

    hetero = TRUE, alpha_0 = 1, alpha_X = 1, alpha_A = 1, alpha_AX = 0.1)

  ate <- dat$ATE[1]

  #--- RUN OVER NSIMS RESAMPLINGS ----------------------------------------------

  idx <- 1:nsims

  with_progress({

    p <- progressor(along = idx)

    sim_results <- foreach(sim = idx,

      .combine = "rbind", .errorhandling = "pass",

      .options.future = list(seed = TRUE)) %dofuture% {

      p()

      #--- SAMPLE FROM SUPER POPULATION ----------------------------------------

      S <- rbinom(N, 1, dat$P_S_cond_AX)

      samp <- dat[S == 1,]

      U_dat <- with(samp,

        data.frame(aa = A, yy = Y, xx = X, ss = P_S_cond_AX))

      U_dat_Oracle <- with(dat,

        data.frame(aa = A, yy = Y, xx = X, ss = P_S_cond_AX, S = S))

      #--- FIT COMMON MODELS ---------------------------------------------------

      samp$s_wt <- 1 / samp$P_S_cond_AX

      #-- Pr(A = 1 | X, S = 1): Within-Sample Propensity Score

      fit_a <- glm(A ~ X, family = "quasibinomial", data = samp)

      P_A1_cond_X <- predict(fit_a, newdata = samp, type = "response")

      samp$ps_wt <- ifelse(samp$A == 1,

        1/P_A1_cond_X, 1/(1 - P_A1_cond_X))

      samp$comb_wt <- with(samp, s_wt*ps_wt)

      #-- Pr(A = 1 | X): Weighted Propensity Score

      des_wtd_a <- svydesign(

        ids = ~ 1, weights = 1 / samp$P_S_cond_AX, data = samp)

      wtd_fit_a <- svyglm(A ~ X, des_wtd_a, family = "quasibinomial",)

      wtd_P_A1_cond_X <- predict(wtd_fit_a, newdata = samp, type = "response")

      samp$wtd_ps_wt <- ifelse(samp$A == 1,

        1/wtd_P_A1_cond_X , 1/(1 - wtd_P_A1_cond_X ))

      samp$wtd_comb_wt <- with(samp, s_wt*wtd_ps_wt)

      #--- FIT COMPARISON MODELS -----------------------------------------------

      #-- 00. Oracle Estimator

      fit_00 <- glm(Y ~ A*X, data = samp)

      est_00 <- mean(

        predict(fit_00, newdata = transform(dat, A = 1)) -

        predict(fit_00, newdata = transform(dat, A = 0)))

      ests_00 <- c(coef(fit_00), ate = est_00)

      halfmeat_00 <- U_Oracle(theta = ests_00, U_dat_Oracle)

      bread_00 <- jacobian(G,

        ests_00, U_data = U_dat_Oracle, U_func = "U_Oracle")

      IF_00 <- halfmeat_00%*%t(solve(-bread_00))

      err_00 <- sqrt(sum(IF_00[, ncol(IF_00)]^2))

      cov_00 <- ifelse(est_00 - qnorm(0.975)*err_00 < ate &

        est_00 + qnorm(0.975)*err_00 > ate, 1, 0)

      #-- 01. Simple Regression

      fit_01 <- glm(Y ~ A, data = samp)

      est_01 <- coef(summary(fit_01))["A", "Estimate"]

      err_01 <- coef(summary(fit_01))["A", "Std. Error"]

      cov_01 <- ifelse(est_01 - qnorm(0.975)*err_01 < ate &

        est_01 + qnorm(0.975)*err_01 > ate, 1, 0)

      #-- 02. Multiple Regression

      est_02 <- het_est(samp, fit_00, 1)

      err_02 <- het_err(samp, fit_00, 1)

      cov_02 <- ifelse(est_02 - qnorm(0.975)*err_02 < ate &

        est_02 + qnorm(0.975)*err_02 > ate, 1, 0)

      #-- 03. IPTW Estimator

      est_03 <- with(samp, mean((A*Y/P_A1_cond_X) -

        ((1 - A)*Y/(1 - P_A1_cond_X))))

      ests_03 <- c(coef(fit_a), ate = est_03)

      halfmeat_03 <- U_IPTW(theta = ests_03, U_dat)

      bread_03 <- jacobian(G, ests_03, U_data = U_dat, U_func = "U_IPTW")

      IF_03 <- halfmeat_03%*%t(solve(-bread_03))

      err_03 <- sqrt(sum(IF_03[, ncol(IF_03)]^2))

      cov_03 <- ifelse(est_03 - qnorm(0.975)*err_03 < ate &

        est_03 + qnorm(0.975)*err_03 > ate, 1, 0)

      #-- 04. Survey-Weighted Multiple Regression

      fit_04 <- glm(Y ~ A*X, data = samp, weights = samp$s_wt)

      est_04 <- het_est(samp, fit_04, 1)

      err_04 <- het_err(samp, fit_04, 1)

      cov_04 <- ifelse(est_04 - qnorm(0.975)*err_04 < ate &

        est_04 + qnorm(0.975)*err_04 > ate, 1, 0)

      #-- 05. IPTW Multiple Regression

      fit_05 <- glm(Y ~ A*X, data = samp, weights = samp$ps_wt)

      est_05 <- het_est(samp, fit_05, 1)

      err_05 <- het_err(samp, fit_05, 1)

      cov_05 <- ifelse(est_05 - qnorm(0.975)*err_05 < ate &

        est_05 + qnorm(0.975)*err_05 > ate, 1, 0)

      #-- 06. IPTW + Survey Weighted Multiple Regression

      fit_06 <- glm(Y ~ A*X, data = samp, weights = samp$comb_wt)

      est_06 <- het_est(samp, fit_06, 1)

      err_06 <- het_err(samp, fit_06, 1)

      cov_06 <- ifelse(est_06 - qnorm(0.975)*err_06 < ate &

        est_06 + qnorm(0.975)*err_06 > ate, 1, 0)

      #-- 07. Weighted IPTW + Survey Weighted Multiple Regression

      fit_07 <- glm(Y ~ A*X, data = samp, weights = samp$wtd_comb_wt)

      est_07 <- het_est(samp, fit_07, 1)

      err_07 <- het_err(samp, fit_07, 1)

      cov_07 <- ifelse(est_07 - qnorm(0.975)*err_07 < ate &

        est_07 + qnorm(0.975)*err_07 > ate, 1, 0)

      #-- 08. Naive G-Computation

      est_08 <- het_est(samp, fit_00, with(samp, s_wt/sum(s_wt)))

      ests_08 <- c(coef(fit_00), ate = est_08)

      halfmeat_08 <- U_naiveG(theta = ests_08, U_dat)

      bread_08 <- jacobian(G, ests_08, U_data = U_dat, U_func = "U_naiveG")

      IF_08 <- halfmeat_08%*%t(solve(-bread_08))

      err_08 <- sqrt(sum(IF_08[, ncol(IF_08)]^2))

      cov_08 <- ifelse(est_08 - qnorm(0.975)*err_08 < ate &

        est_08 + qnorm(0.975)*err_08 > ate, 1, 0)

      #-- 09. Identification Formula 1

      fit_09 <- svycdiff(samp, "OM",

        A ~ X, P_S_cond_AX ~ A + X, Y ~ A*X, "gaussian")

      est_09 <- fit_09$cdiff[1]

      err_09 <- fit_09$cdiff[2]

      cov_09 <- ifelse(fit_09$cdiff[3] < ate & fit_09$cdiff[4] > ate, 1, 0)

      #-- 10. Identification Formula 2A

      fit_10 <- svycdiff(samp, "IPW1",

        A ~ X, P_S_cond_AX ~ A + X, Y ~ 1, NULL)

      est_10 <- fit_10$cdiff[1]

      err_10 <- fit_10$cdiff[2]

      cov_10 <- ifelse(fit_10$cdiff[3] < ate & fit_10$cdiff[4] > ate, 1, 0)

      #-- 11. Identification Formula 2B

      fit_11 <- svycdiff(samp, "IPW2",

        A ~ X, P_S_cond_AX ~ A + X, Y ~ 1, NULL)

      est_11 <- fit_11$cdiff[1]

      err_11 <- fit_11$cdiff[2]

      cov_11 <- ifelse(fit_11$cdiff[3] < ate & fit_11$cdiff[4] > ate, 1, 0)

      #--- STORE SIMULATION RESULTS --------------------------------------------

      sim_result <- c(

        ate = ate,

        est_00 = est_00, err_00 = err_00, cov_00 = cov_00,
        est_01 = est_01, err_01 = err_01, cov_01 = cov_01,
        est_02 = est_02, err_02 = err_02, cov_02 = cov_02,
        est_03 = est_03, err_03 = err_03, cov_03 = cov_03,
        est_04 = est_04, err_04 = err_04, cov_04 = cov_04,
        est_05 = est_05, err_05 = err_05, cov_05 = cov_05,
        est_06 = est_06, err_06 = err_06, cov_06 = cov_06,
        est_07 = est_07, err_07 = err_07, cov_07 = cov_07,
        est_08 = est_08, err_08 = err_08, cov_08 = cov_08,
        est_09 = est_09, err_09 = err_09, cov_09 = cov_09,
        est_10 = est_10, err_10 = err_10, cov_10 = cov_10,
        est_11 = est_11, err_11 = err_11, cov_11 = cov_11)

      return(sim_result)
    }
  })

  results[[i]] <- sim_results
}

toc <- proc.time(); toc - tic

#=== END =======================================================================
