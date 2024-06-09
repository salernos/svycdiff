#===============================================================================
#
#  PROGRAM: SIM_SENSITIVITY.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: To perform a sensitivity analysis in order to understand the
#           robustness of our outcome-model based approach to potential
#           model misspecification.
#
#  INUPT:   svycdiff.R - R script containing functions for the proposed methods
#
#           utils.R - R script containing helper functions for methods
#
#           sim.R - R script containing functions to simulate data
#
#  OUTPUT:  SENSITIVITY_RESULTS.RData
#
#           An .RData file containing the object `results`, which is a list of
#           five elements corresponding to the five sensitivity analysis
#           settings described in the manuscript. Each list element contains a
#           numeric matrix of size `nsims` x `16`, where each row corresponds
#           to a simulation replicate and each column is a saved metric.
#
#  DETAILS: Methods Being Compared:
#
#             00.  Oracle Estimator
#             09.  Identification Formula 1
#             09a. Identification Formula 1, Misspecified
#             10.  Identification Formula 2A
#             11.  Identification Formula 2B
#
#  UPDATED: 2024-05-25
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- SOURCE FUNCTIONS ----------------------------------------------------------

source(here("R", "svycdiff.R"))

source(here("R", "utils.R"))

source(here("R", "sim.R"))

#--- SOURCE NECESSARY PACKAGES -------------------------------------------------

library(pacman)

p_load(progressr, doFuture, numDeriv, betareg, survey, gt, update = F)

#--- HELPER FUNCTIONS ----------------------------------------------------------

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

G <- function(theta, U_data, U_func) {

  U <- get(U_func)

  return(apply(U(theta, U_data), 2, sum))
}

#=== GLOBAL VARIABLES ==========================================================

gamma_AX_vec <- c(0, 0.01, 0.05, 0.1, 0.5)

nsims <- 200

N <- 100000

#=== SIMULATION ================================================================

plan(multisession)

results <- vector("list", length(gamma_AX_vec))

tic <- proc.time()

for(i in 1:length(results)) {

  cat("\nSetting", i, "---\n\n")

  #--- GENERATE SUPER POPULATION -----------------------------------------------

  set.seed(i)

  gamma_AX  <- gamma_AX_vec[i]

  dat <- simdat(N, X_dist = "continuous", S_known = F,

    tau_0 = -1, tau_X = 1,

    beta_0 = -4.5, beta_A = 1, beta_X = 1,

    hetero = T, alpha_0 = 1, alpha_X = 1, alpha_A = 1, alpha_AX = gamma_AX)

  ate <- dat$ATE[1]

  #--- RUN OVER NSIMS RESAMPLINGS ----------------------------------------------

  idx <- 1:nsims

  with_progress({

    p <- progressor(along = idx)

    sim_results <- foreach(sim = idx,

      .combine = "rbind", .errorhandling = "pass",

      .options.future = list(seed = T)) %dofuture% {

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

      #-- 09. Identification Formula 1

      fit_09 <- svycdiff(samp, "OM",

        A ~ X, P_S_cond_AX ~ A + X, Y ~ A * X, "gaussian")

      est_09 <- fit_09$cdiff[1]

      err_09 <- fit_09$cdiff[2]

      cov_09 <- ifelse(fit_09$cdiff[3] < ate & fit_09$cdiff[4] > ate, 1, 0)

      #-- 09. Identification Formula 1, Misspecified

      fit_09a <- svycdiff(samp, "OM",

        A ~ X, P_S_cond_AX ~ A + X, Y ~ A + X, "gaussian")

      est_09a <- fit_09a$cdiff[1]

      err_09a <- fit_09a$cdiff[2]

      cov_09a <- ifelse(fit_09a$cdiff[3] < ate & fit_09a$cdiff[4] > ate, 1, 0)

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

        est_00  = est_00,  err_00  = err_00,  cov_00  = cov_00,
        est_09  = est_09,  err_09  = err_09,  cov_09  = cov_09,
        est_09a = est_09a, err_09a = err_09a, cov_09a = cov_09a,
        est_10  = est_10,  err_10  = err_10,  cov_10  = cov_10,
        est_11  = est_11,  err_11  = err_11,  cov_11  = cov_11)

      return(sim_result)
    }
  })

  results[[i]] <- sim_results
}

toc <- proc.time(); toc - tic

save(results, file = here("inst", "SENSITIVITY_RESULTS.RData"))

#=== END =======================================================================
