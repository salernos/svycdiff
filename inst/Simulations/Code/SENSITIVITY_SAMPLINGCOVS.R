#===============================================================================
#
#  PROGRAM: SENSITIVITY_SAMPLINGCOVS.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Sensitivity analysis to evaluate the performance of the proposed
#           methods against standard analytic methods and the oracle estimator
#           when we omit unobserved covariates that influence sample selection
#           but not the treatment or outcome models from our modeling.
#
#  INPUT:   main.R    - R script implementing the proposed methods
#           helpers.R - R script containing helper functions for the methods
#           sim.R     - R script containing functions to simulate data
#
#  OUTPUT:  SENSITIVITY_SAMPLINGCOVS.RData
#
#           An .RData file containing the object `results`, which is a tibble
#           of size (nsims x nsettings) x 43, where each row corresponds to a
#           simulation replicate for a particular setting and each column is
#           metric for a particular method. Metrics are the true ACD, and the
#           estimated ACD (est_), the estimated standard error (err_) and the
#           coverage (cov_) for each method. Methods being compared are:
#
#             00. Oracle Estimator
#             01. Simple Regression
#             02. Multiple Regression
#             03. IPTW Estimator
#             04. Survey-Weighted Multiple Regression
#             05. IPTW Multiple Regression
#             06. IPTW + Survey Weighted Multiple Regression
#             07. Weighted IPTW + Survey Weighted Multiple Regression
#             08. Proposed Outcome Modeling Estimator
#             09. Proposed Inverse Probability Weighting Estimator 1
#             10. Proposed Inverse Probability Weighting Estimator 2
#             11. Proposed Doubly-Robust Estimator
#
#           The remaining columns contain simulation-specific parameters such
#           as the simulation setting id (id), the population size (N), the
#           number of X covariates (p), and the sizes of the main effects in
#           the true propensity (tau_X) and selection (beta_A, beta_X) models.
#
#  UPDATED: 2025-07-22
#
#===============================================================================

#=== SETUP =====================================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(pacman)

p_load(progressr, furrr, MASS, numDeriv, betareg, survey, gt, here, tidyverse,

  update = FALSE)

#--- SOURCE PACKAGE FUNCTIONS --------------------------------------------------

source(here("R", "main.R"))

source(here("R", "helpers.R"))

source(here("R", "sim.R"))

#--- HELPER FUNCTIONS ----------------------------------------------------------

#-- BUILD FORMULAS AND ORACLE MATRICES

build_forms <- function(pop_dat) {

  x_cols <- grep("^X", names(pop_dat), value = TRUE)

  u_cols <- grep("^U", names(pop_dat), value = TRUE)

  #- Formulas

  form_a <- as.formula(paste0("A ~ ", paste(x_cols, collapse = " + ")))

  form_s <- as.formula(paste0("pS ~ A + ",

    paste(c(x_cols, u_cols), collapse = " + ")))

  form_y <- as.formula(paste0("Y ~ A + ",

    paste(x_cols, collapse = " + "), " + ",

    paste0("A:", x_cols, collapse = " + ")))

  #- Oracle Design Matrices

  OR_f <- as.formula(paste0("~ A + ",

    paste(x_cols, collapse = " + "), " + ",

    paste0("A:", x_cols, collapse = " + ")))

  X_or  <- model.matrix(OR_f, data = pop_dat)

  X1_or <- model.matrix(OR_f, data = transform(pop_dat, A = 1))

  X0_or <- model.matrix(OR_f, data = transform(pop_dat, A = 0))

  ate <- pop_dat$CDIFF[1]

  list(

    form_a = form_a,
    form_s = form_s,
    form_y = form_y,
    X_or   = X_or,
    X1_or  = X1_or,
    X0_or  = X0_or,
    ate    = ate
  )
}

#-- RUN ONE MONTE CARLO REPLICATE

run_one_sim <- function(sim_i, pop_dat, forms) {

  set.seed(sim_i)

  #- Sample

  S_ind <- rbinom(nrow(pop_dat), 1, pop_dat$pS)

  samp <- pop_dat[S_ind == 1, ]

  samp$s_wt <- 1 / samp$pS

  #- Common Pre-Computations

  form_a <- forms$form_a
  form_s <- forms$form_s
  form_y <- forms$form_y

  X_or  <- forms$X_or
  X1_or <- forms$X1_or
  X0_or <- forms$X0_or

  ate <- forms$ate

  #- Within‐Sample Propensity Score

  A_X <- model.matrix(form_a, samp)

  fit_a <- glm(form_a, data = samp, family = quasibinomial())

  pA <- predict(fit_a, samp, type = "response")

  samp$ps_wt <- ifelse(samp$A == 1, 1 / pA, 1 / (1 - pA))

  samp$comb_wt <- samp$s_wt * samp$ps_wt

  #- Weighted Propensity Score

  des_wtd_a <- svydesign(ids = ~ Cluster, strata = ~ Strata,

    weights = samp$s_wt, data = samp)

  wtd_fit_a <- svyglm(form_a, design = des_wtd_a, family = quasibinomial())

  pA_wtd  <- predict(wtd_fit_a, samp, type = "response")

  samp$wtd_ps_wt <- ifelse(samp$A == 1, 1 / pA_wtd, 1 / (1 - pA_wtd))

  samp$wtd_comb_wt <- samp$s_wt * samp$wtd_ps_wt

  #- Estimators

  # 00. Oracle Estimator

  U_Oracle <- function(theta) {

    p_dims <- ncol(X_or)

    gamma  <- theta[1:p_dims]

    CDIFF  <- theta[p_dims + 1]

    U_Y <- X_or * c(pop_dat$Y - X_or %*% gamma) * S_ind

    U_C <- CDIFF - ((X1_or - X0_or) %*% gamma)

    cbind(U_Y, U_C)
  }

  fit_00 <- glm(form_y, data = samp, family = gaussian())

  est_00 <- mean(

    predict(fit_00, newdata = transform(pop_dat, A = 1)) -

    predict(fit_00, newdata = transform(pop_dat, A = 0)))

  ests_00 <- c(coef(fit_00), CDIFF = est_00)

  halfmeat_00 <- U_Oracle(ests_00)

  bread_00 <- numDeriv::jacobian(function(t) colSums(U_Oracle(t)), ests_00)

  IF_00 <- halfmeat_00 %*% t(solve(-bread_00))

  err_00 <- sqrt(sum(IF_00[, ncol(IF_00)]^2))

  cov_00 <- as.integer(est_00 - 1.96 * err_00 < ate &

    est_00 + qnorm(0.975) * err_00 > ate)

  # 01. Simple Regression

  fit_01 <- glm(Y ~ A, data = samp)

  est_01 <- coef(fit_01)["A"]

  err_01 <- summary(fit_01)$coef["A", "Std. Error"]

  cov_01 <- as.integer(est_01 - qnorm(0.975) * err_01 < ate &

    est_01 + qnorm(0.975) * err_01 > ate)

  # 02. Multiple Regression

  het_est <- function(fit, design, data, treat = "A") {

    coefs <- coef(fit)

    cn <- names(coefs)

    ints <- grep(paste0("^", treat, ":"), cn, value = TRUE)

    ints <- unique(c(ints, grep(paste0(":", treat, "$"), cn, value = TRUE)))

    if (length(ints) == 0) {

      stop("No interactions of the form ", treat, ":X found in model.")
    }

    x_vars <- gsub(paste0("(^", treat, ":|:", treat, "$)"), "", ints)

    wts <- weights(design)

    wts <- if (length(wts)==1) rep(wts, nrow(data)) else wts

    wts_norm <- wts / sum(wts)

    Xmeans <- vapply(x_vars,

      function(xj) { sum(data[[xj]] * wts_norm) }, numeric(1))

    cont <- setNames(rep(0, length(coefs)), cn)

    cont[treat] <- 1

    cont[ints]  <- Xmeans

    svycontrast(fit, cont)
  }

  des_02 <- svydesign(ids = ~ 1, data = samp)

  fit_02 <- svyglm(form_y, design = des_02)

  het_02 <- het_est(fit_02, des_02, samp, treat = "A")

  est_02 <- as.numeric(het_02)

  err_02 <- sqrt(vcov(het_02)[1, 1])

  cov_02 <- as.integer(est_02 - qnorm(0.975) * err_02 < ate &

    est_02 + qnorm(0.975) * err_02 > ate)

  # 03. IPTW Estimator

  U_IPTW <- function(theta) {

    p_dims <- ncol(A_X)

    tau <- theta[1:p_dims]

    CDIFF <- theta[p_dims + 1]

    pA <- plogis(A_X %*% tau)

    U_A <- A_X * c(samp$A - pA)

    U_C <- CDIFF - (samp$A * samp$Y / pA - (1 - samp$A) * samp$Y / (1 - pA))

    cbind(U_A, U_C)
  }

  est_03 <- with(samp, mean(A * Y / pA - (1 - A) * Y / (1 - pA)))

  ests_03 <- c(coef(fit_a), CDIFF = est_03)

  halfmeat_03 <- U_IPTW(ests_03)

  bread_03 <- numDeriv::jacobian(function(t) colSums(U_IPTW(t)), ests_03)

  IF_03 <- halfmeat_03 %*% t(solve(-bread_03))

  err_03 <- sqrt(sum(IF_03[, ncol(IF_03)]^2))

  cov_03 <- as.integer(est_03 - qnorm(0.975) * err_03 < ate &

    est_03 + qnorm(0.975) * err_03 > ate)

  # 04. Survey-Weighted Multiple Regression

  des_04 <- svydesign(ids = ~ Cluster, strata = ~ Strata,

    weights = 1 / samp$pS, data = samp, nest = TRUE)

  fit_04 <- svyglm(form_y, design = des_04)

  het_04 <- het_est(fit_04, des_04, samp, treat = "A")

  est_04 <- as.numeric(het_04)

  err_04 <- sqrt(vcov(het_04)[1, 1])

  cov_04 <- as.integer(est_04 - qnorm(0.975) * err_04 < ate &

    est_04 + qnorm(0.975) * err_04 > ate)

  # 05. IPTW-Weighted Regression

  des_05 <- svydesign(ids = ~ 1, weights = samp$ps_wt, data = samp)

  fit_05 <- svyglm(form_y, design = des_05)

  het_05 <- het_est(fit_05, des_05, samp, treat = "A")

  est_05 <- as.numeric(het_05)

  err_05 <- sqrt(vcov(het_05)[1, 1])

  cov_05 <- as.integer(est_05 - qnorm(0.975) * err_05 < ate &

    est_05 + qnorm(0.975) * err_05 > ate)

  # 06. IPTW + Survey-Weighted Regression

  des_06 <- svydesign(ids = ~ Cluster, strata = ~ Strata,

    weights = samp$comb_wt, data = samp, nest = TRUE)

  fit_06 <- svyglm(form_y, design = des_06)

  het_06 <- het_est(fit_06, des_06, samp, treat = "A")

  est_06 <- as.numeric(het_06)

  err_06 <- sqrt(vcov(het_06)[1, 1])

  cov_06 <- as.integer(est_06 - qnorm(0.975) * err_06 < ate &

    est_06 + qnorm(0.975) * err_06 > ate)

  # 07. Weighted IPTW + Survey-Weighted Regression

  des_07 <- svydesign(ids = ~ Cluster, strata = ~ Strata,

    weights = samp$wtd_comb_wt, data = samp, nest = TRUE)

  fit_07 <- svyglm(form_y, design = des_07)

  het_07 <- het_est(fit_07, des_07, samp, treat="A")

  est_07 <- as.numeric(het_07)

  err_07 <- sqrt(vcov(het_07)[1, 1])

  cov_07 <- as.integer(est_07 - qnorm(0.975) * err_07 < ate &

    est_07 + qnorm(0.975) * err_07 > ate)

  # 08. Proposed Outcome Modeling Estimator

  fit_08 <- svycdiff(samp, "OM", form_a, form_s, form_y,

    y_fam = "gaussian", strata = "Strata", cluster = "Cluster")

  est_08 <- fit_08$cdiff[1]

  err_08 <- fit_08$cdiff[2]

  cov_08 <- as.integer(fit_08$cdiff[3] < ate & fit_08$cdiff[4] > ate)

  # 09. Proposed Inverse Probability Weighting Estimator 1

  fit_09 <- svycdiff(samp, "IPW1", form_a, form_s, Y ~ 1,

    y_fam = NULL, strata = "Strata", cluster = "Cluster")

  est_09 <- fit_09$cdiff[1]

  err_09 <- fit_09$cdiff[2]

  cov_09 <- as.integer(fit_09$cdiff[3] < ate & fit_09$cdiff[4] > ate)

  # 10. Proposed Inverse Probability Weighting Estimator 2

  fit_10 <- svycdiff(samp, "IPW2", form_a, form_s, Y ~ 1,

    y_fam = NULL, strata = "Strata", cluster = "Cluster")

  est_10 <- fit_10$cdiff[1]

  err_10 <- fit_10$cdiff[2]

  cov_10 <- as.integer(fit_10$cdiff[3] < ate & fit_10$cdiff[4] > ate)

  # 11. Proposed Doubly-Robust Estimator

  fit_11 <- svycdiff(samp, "DR", form_a, form_s, form_y,

    y_fam = "gaussian", strata = "Strata", cluster = "Cluster")

  est_11 <- fit_11$cdiff[1]

  err_11 <- fit_11$cdiff[2]

  cov_11 <- as.integer(fit_11$cdiff[3] < ate & fit_11$cdiff[4] > ate)

  #- Combine Into a Named Vector

  tibble::tibble(ate = ate,

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
    est_11 = est_11, err_11 = err_11, cov_11 = cov_11
  )
}

#=== SIMULATION STUDY ==========================================================

#--- SIMULATION SETTINGS -------------------------------------------------------

nsims <- 500

settings <- expand.grid(

  N           = 1e6,
  p           = 2,
  q           = 2,     #- Adding Sampling Covariates

  #- Design Variables

  n_strat     = 10,
  n_clust     = 20,
  sigma_strat = 0.35,
  sigma_clust = 0.15,
  X_fam       = "gaussian",

  #- Propensity Model

  tau_0       = 0,
  tau_A       = 0.5,
  tau_X_val   = c(0, 1),
  tau_X12     = 0,

  #- Selection Model

  beta_0      = -8,
  beta_A      = c(0, 1),
  beta_X_val  = c(0, 1),
  beta_U_val  = 1,

  #- Outcome Model

  Y_fam       = "gaussian",
  alpha_0     = 1,
  alpha_A     = 1,
  alpha_X_val = 1,
  alpha_AX    = 0.1,

  stringsAsFactors = FALSE) |>

  mutate(id = row_number())

#--- SIMULATION ----------------------------------------------------------------

plan(multisession)

handlers("txtprogressbar")

total_steps <- nrow(settings) * nsims

with_progress({

  pr <- progressor(steps = total_steps)

  results <- future_pmap_dfr(

    settings,

    function(N, p, q, n_strat, n_clust, sigma_strat, sigma_clust, X_fam,

      tau_0, tau_A, tau_X_val, tau_X12, beta_0, beta_A, beta_X_val, beta_U_val,

      Y_fam, alpha_0, alpha_A, alpha_X_val, alpha_AX, id) {

      #-- BUILD COEFFICIENT VECTORS

      tau_X   <- rep(tau_X_val,   p)
      beta_X  <- rep(beta_X_val,  p)
      beta_U  <- rep(beta_U_val,  q)
      alpha_X <- rep(alpha_X_val, p)

      #-- GENERATE POPULATION

      pop_dat <- simdat(

        N           = N,
        p           = p,
        q           = q,
        n_strat     = n_strat,
        n_clust     = n_clust,
        sigma_strat = sigma_strat,
        sigma_clust = sigma_clust,
        X_fam       = X_fam,
        tau_0       = tau_0,
        tau_A       = tau_A,
        tau_X       = tau_X,
        tau_X12     = tau_X12,
        beta_0      = beta_0,
        beta_A      = beta_A,
        beta_X      = beta_X,
        beta_U      = beta_U,
        Y_fam       = Y_fam,
        alpha_0     = alpha_0,
        alpha_A     = alpha_A,
        alpha_X     = alpha_X,
        alpha_AX    = alpha_AX
      )

      #-- BUILD FORMULAS AND ORACLE MATRICES

      forms <- build_forms(pop_dat)

      #-- RUN REPLICATES IN PARALLEL

      sims <- future_map_dfr(

        seq_len(nsims),

        function(sim_i) {

          pr()

          run_one_sim(sim_i, pop_dat, forms)
        },

        .options = furrr_options(seed = TRUE)
      )

      sims |>

        mutate(

          id = id,

          N = N,
          p = p,

          tau_X  = tau_X_val,
          beta_A = beta_A,
          beta_X = beta_X_val
        )
    },

    .options = furrr_options(seed = TRUE)
  )
})

#--- UNCOMMENT TO SAVE RESULTS -------------------------------------------------

# save(results, file = here("inst", "Simulations", "Results",
#
#   "SENSITIVITY_SAMPLINGCOVS.RData"))

#=== SUMMARIZE SIMULATION RESULTS =============================================

# NOTE: Assumes `results` is a tibble with one row per setting id × replicate

#-- DEFINE METHOD LABELS

method_codes <- sprintf("%02d", 0:11)

method_labels <- c(

  "Oracle Estimator",
  "Simple Regression",
  "Multiple Regression",
  "IPTW Estimator",
  "Survey-Weighted Multiple Regression",
  "IPTW Multiple Regression",
  "IPTW + Survey-Weighted Multiple Regression",
  "Weighted IPTW + Survey-Weighted Multiple Regression",
  "Outcome Modeling and Direct Standardization",
  "Inverse Probability Weighting 1",
  "Inverse Probability Weighting 2",
  "Augmented Inverse Probability Weighting"
)

#-- PIVOT TO LONG FORMAT

results_long <- results |>

  group_by(id) |>

  mutate(replicate = row_number()) |>

  ungroup() |>

  pivot_longer(

    cols = matches("^(est|err|cov)_\\d\\d$"),
    names_to = c("metric", "method_code"),
    names_sep = "_") |>

  pivot_wider(names_from = metric, values_from = value) |>

  mutate(method = factor(method_code, method_codes, method_labels))

#-- COMPUTE SUMMARY METRICS BY `id` AND `method`

results_summary <- results_long |>

  group_by(id, method) |>

  summarise(

    `True ATE`       = first(ate),
    `Estimated ATE`  = mean(est),
    `Relative Bias`  = mean((est - `True ATE`) / `True ATE`),
    `Analytic SE`    = mean(err),
    `Monte Carlo SE` = sd(est),
    `MSE`            = mean((est - `True ATE`)^2),
    `Coverage`       = mean(cov),

    .groups = "drop")

#-- CREATE GT TABLE

results_gt <- results_summary |>

  arrange(id, method) |>

  gt(groupname_col = "id") |>

  tab_header(

    title = "Full Simulation Results",

    subtitle = paste0("Number of Simulations = ", nsims)) |>

  fmt_number(columns = c(`True ATE`:`Coverage`), decimals = 3)

results_gt

#-- CREATE LATEX TABLE

results_gt |>

  as_latex() |>

  cat()

#=== END =======================================================================
