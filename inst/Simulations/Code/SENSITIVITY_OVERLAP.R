#===============================================================================
#
#  PROGRAM: SENSITIVITY_OVERLAP.R
#
#  AUTHORS: Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Sensitivity analysis to evaluate the performance of the proposed
#           methods against the oracle estimator for varying degrees of
#           overlap (good versus poor) in the propensity score distribution.
#
#  INPUT:   main.R    - R script implementing the proposed methods
#           helpers.R - R script containing helper functions for the methods
#           sim.R     - R script containing functions to simulate data
#
#  OUTPUT:  SENSITIVITY_OVERLAP_SKNOWN.RData
#
#           An .RData file containing the object `results`, which is a tibble
#           of size (nsims x nsettings) x 20, where each row corresponds to a
#           simulation replicate for a particular setting and each column is
#           metric for a particular method. Metrics are the true ACD, and the
#           estimated ACD (est_), the estimated standard error (err_) and the
#           coverage (cov_) for each method. Methods being compared are:
#
#             00. Oracle Estimator
#             01. Proposed Outcome Modeling Estimator
#             02. Proposed Inverse Probability Weighting Estimator 1
#             03. Proposed Inverse Probability Weighting Estimator 2
#             04. Proposed Doubly-Robust Estimator
#
#           The remaining columns contain simulation-specific parameters such
#           as the simulation setting id (id), the population size (N), the
#           number of X covariates (p), and the degree of overlap in the true
#           propensity score model (tau_A).
#
#  UPDATED: 2025-07-22
#
#===============================================================================

#=== SETUP =====================================================================

#--- LOAD NECESSARY PACKAGES ---------------------------------------------------

library(pacman)

p_load(progressr, furrr, MASS, numDeriv, betareg, survey, gt, here, patchwork,

  tidyverse, update = FALSE)

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

  form_s <- as.formula("pS ~ 1")

  form_y <- as.formula(paste0("Y ~ A + ",

    paste(x_cols, collapse = " + "), " + ",

    paste0("A:", x_cols, collapse = " + ")))

  #- Oracle Design Matrices

  OR_f <- as.formula(paste0("Y ~ A + ",

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
    OR_f   = OR_f,
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

  fit_00 <- glm(forms$OR_f, data = samp, family = gaussian())

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

  # 01. Proposed Outcome Modeling Estimator

  fit_01 <- svycdiff(samp, "OM", form_a, form_s, form_y,

    y_fam = "gaussian", strata = "Strata", cluster = "Cluster")

  est_01 <- fit_01$cdiff[1]

  err_01 <- fit_01$cdiff[2]

  cov_01 <- as.integer(fit_01$cdiff[3] < ate & fit_01$cdiff[4] > ate)

  # 02. Proposed Inverse Probability Weighting Estimator 1

  fit_02 <- svycdiff(samp, "IPW1", form_a, form_s, Y ~ 1,

    y_fam = NULL, strata = "Strata", cluster = "Cluster")

  est_02 <- fit_02$cdiff[1]

  err_02 <- fit_02$cdiff[2]

  cov_02 <- as.integer(fit_02$cdiff[3] < ate & fit_02$cdiff[4] > ate)

  # 03. Proposed Inverse Probability Weighting Estimator 2

  fit_03 <- svycdiff(samp, "IPW2", form_a, form_s, Y ~ 1,

    y_fam = NULL, strata = "Strata", cluster = "Cluster")

  est_03 <- fit_03$cdiff[1]

  err_03 <- fit_03$cdiff[2]

  cov_03 <- as.integer(fit_03$cdiff[3] < ate & fit_03$cdiff[4] > ate)

  # 04. Proposed Doubly-Robust Estimator

  fit_04 <- svycdiff(samp, "DR", form_a, form_s, form_y,

    y_fam = "gaussian", strata = "Strata", cluster = "Cluster")

  est_04 <- fit_04$cdiff[1]

  err_04 <- fit_04$cdiff[2]

  cov_04 <- as.integer(fit_04$cdiff[3] < ate & fit_04$cdiff[4] > ate)

  #- Combine Into a Named Vector

  tibble::tibble(ate = ate,

    est_00 = est_00, err_00 = err_00, cov_00 = cov_00,
    est_01 = est_01, err_01 = err_01, cov_01 = cov_01,
    est_02 = est_02, err_02 = err_02, cov_02 = cov_02,
    est_03 = est_03, err_03 = err_03, cov_03 = cov_03,
    est_04 = est_04, err_04 = err_04, cov_04 = cov_04
  )
}

#=== SIMULATION STUDY ==========================================================

#--- SIMULATION SETTINGS -------------------------------------------------------

nsims <- 500

settings <- expand.grid(

  N           = 1e6,
  p           = 2,
  q           = 0,

  #- Design Variables

  n_strat     = 10,
  n_clust     = 20,
  sigma_strat = 0.35,
  sigma_clust = 0.15,
  X_fam       = "gaussian",

  #- Propensity Model

  tau_0       = 0,
  tau_A       = c(0.5, 2),
  tau_X_val   = 1,
  tau_X12     = 0,

  #- Selection Model

  beta_0      = -8,
  beta_A      = 1,
  beta_X_val  = 1,
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

          tau_A  = tau_A
        )
    },

    .options = furrr_options(seed = TRUE)
  )
})

#--- UNCOMMENT TO SAVE RESULTS -------------------------------------------------

# save(results, file = here("inst", "Simulations", "Results",
#
#   "SENSITIVITY_OVERLAP_SKNOWN.RData"))

#=== SUMMARIZE SIMULATION RESULTS ==============================================

#--- EXAMPLE PROPENSITY SCORE DISTRIBUTIONS ------------------------------------

dat_good <- simdat(N = 1e6, p = 2, q = 0, n_strat = 10, n_clust = 20,

  sigma_strat = 0.35, sigma_clust = 0.15, X_fam = "gaussian",

  tau_0 = 0, tau_A = 0.5, tau_X = rep(1, p), tau_X12 = 0,

  beta_0 = -8, beta_A = 1, beta_X = rep(1, p), beta_U = rep(0, q),

  Y_fam = "gaussian", alpha_0 = 1, alpha_A = 1, alpha_X = rep(1, p),

  alpha_AX = 0.1) |>

  mutate(A = factor(A, c(0,1),

    c("Comparison Group (A = 0)", "Group of Interest (A = 1)")))

dat_poor <- simdat(N = 1e6, p = 2, q = 0, n_strat = 10, n_clust = 20,

  sigma_strat = 0.35, sigma_clust = 0.15, X_fam = "gaussian",

  tau_0 = 0, tau_A = 2, tau_X = rep(1, p), tau_X12 = 0,

  beta_0 = -8, beta_A = 1, beta_X = rep(1, p), beta_U = rep(0, q),

  Y_fam = "gaussian", alpha_0 = 1, alpha_A = 1, alpha_X = rep(1, p),

  alpha_AX = 0.1) |>

  mutate(A = factor(A, c(0,1),

    c("Comparison Group (A = 0)", "Group of Interest (A = 1)")))

plt_overlap1 <- ggplot(dat_good, aes(x = pA, group = A, fill = A)) +

  theme_bw() +

  geom_histogram(position = "identity", alpha = 0.5, color = "black") +

  scale_fill_manual(values = c("#AA4AC4", "#00C1D5")) +

  labs(

    title = "Good Overlap Scenario", fill = "Group Membership:",

    x = "\nPropensity Score", y = "Frequency\n")

plt_overlap2 <- ggplot(dat_poor, aes(x = pA, group = A, color = A, fill = A)) +

  theme_bw() +

  geom_histogram(position = "identity", alpha = 0.5, color = "black") +

  scale_fill_manual(values = c("#AA4AC4", "#00C1D5")) +

  labs(

    title = "Poor Overlap Scenario", fill = "Group Membership:",

    x = "\nPropensity Score", y = "Frequency\n")

plt_overlap <- plt_overlap1 + plt_overlap2 +

  plot_annotation(tag_levels = "A") +

  plot_layout(guides = "collect", axes = "collect") &

  theme(legend.position = "bottom",

    axis.title = element_text(size = 14),

    axis.text  = element_text(size = 12))

plt_overlap

#--- UNCOMMENT TO SAVE PLOTS ---------------------------------------------------

# ggsave(here("inst", "Simulations", "Figures", "plt_overlap.png"),
#
#   plt_overlap, height=6, width=12, units="in")

#--- SIMULATION RESULTS --------------------------------------------------------

# NOTE: Assumes `results` is a tibble with one row per setting id × replicate

#-- DEFINE METHOD LABELS

method_codes <- sprintf("%02d", 0:4)

method_labels <- c(

  "Oracle Estimator",
  "Outcome Modeling and Direct Standardization",
  "Inverse Probability Weighting 1",
  "Inverse Probability Weighting 2",
  "Augmented Inverse Probability Weighting"
)

#-- PIVOT TO LONG FORMAT

results_long <- results |>

  mutate(

    Setting = if_else(tau_A == 0.5, "Good Overlap", "Poor Overlap")) |>

  group_by(id) |>

  mutate(replicate = row_number()) |>

  ungroup() |>

  pivot_longer(

    cols = matches("^(est|err|cov)_\\d\\d$"),
    names_to = c("metric", "method_code"),
    names_sep = "_") |>

  pivot_wider(names_from = metric, values_from = value) |>

  mutate(Method = factor(method_code, method_codes, method_labels))

#-- COMPUTE SUMMARY METRICS BY `id` AND `method`

results_summary <- results_long |>

  group_by(Setting, Method) |>

  summarise(

    `True ATE`       = first(ate),
    `Estimated ATE`  = mean(est),
    `Relative Bias`  = mean((est - `True ATE`) / `True ATE`),
    `Analytic SE`    = mean(err),
    `Monte Carle SE` = sd(est),
    `MSE`            = mean((est - `True ATE`)^2),
    `Coverage`       = mean(cov),

    .groups = "drop")

#-- CREATE GT TABLE

results_gt <- results_summary |>

  arrange(Setting, Method) |>

  gt(groupname_col = "Setting") |>

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
