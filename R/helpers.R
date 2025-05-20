
#=== HELPERS ===================================================================

#--- 1. PREP DATA --------------------------------------------------------------

# NOTES: Validate Inputs, Create df0/1, Build Model Matrices, Detect `S_known`

prep_data <- function(df, a_form, s_form, y_form) {

  #–– ASSERTIONS ON INPUTS

  stopifnot(is.data.frame(df))

  stopifnot(inherits(a_form, "formula"),
            inherits(s_form, "formula"),
            inherits(y_form, "formula"))

  A_var <- all.vars(a_form)[1]
  S_var <- all.vars(s_form)[1]
  Y_var <- all.vars(y_form)[1]

  stopifnot(A_var %in% names(df), Y_var %in% names(df), S_var %in% names(df))

  aa <- df[[A_var]]
  ss <- df[[S_var]]
  yy <- df[[Y_var]]

  stopifnot(is.numeric(aa), all(aa %in% c(0, 1)))

  S_known <- ifelse(s_form[[3]] == 1, TRUE, FALSE)

  #-- COUTERFACTUAL DATASETS

  df0 <- df1 <- df

  df0[[A_var]] <- 0L
  df1[[A_var]] <- 1L

  #-- PRECOMPUTE DESIGN MATRICES

  X_a  <- model.matrix(a_form, df)
  X_s  <- model.matrix(s_form, df)
  X_s1 <- model.matrix(s_form, df1)
  X_s0 <- model.matrix(s_form, df0)
  X_y  <- model.matrix(y_form, df)
  X_y1 <- model.matrix(y_form, df1)
  X_y0 <- model.matrix(y_form, df0)

  #-- RETURN

  list(

    df = df, df0 = df0, df1 = df1,

    A = aa, Y = yy, S = ss, S_known = S_known,

    X_a = X_a,
    X_s = X_s, X_s1 = X_s1, X_s0 = X_s0,
    X_y = X_y, X_y1 = X_y1, X_y0 = X_y0
  )
}

#--- 2. FIT MODELS -------------------------------------------------------------

# NOTES: Fit Every Model Depending on `id_form` and `S_known`

#' @importFrom stats quasibinomial

fit_models <- function(prep, id_form, a_form, s_form, y_form, y_fam = NULL) {

  df <- prep$df

  #-- PROPENSITY MODELS

  #- Within-Sample Propensity Model

  fit_a <- glm(a_form, data = df, family = quasibinomial())

  #- Weighted Propensity Model

  des <- survey::svydesign(ids = ~1, weights = 1/prep$S, data = df)

  wtd_fit_a <- survey::svyglm(a_form, design = des, family = quasibinomial())

  #-- SELECTION MODEL

  fit_s <- if (prep$S_known) NULL else betareg::betareg(s_form, data = df)

  #-- OUTCOME MODEL (ONLY OM & DR)

  fit_y <- NULL

  if (id_form %in% c("OM", "DR")) {

    stopifnot(!is.null(y_fam), y_fam %in% c("gaussian", "binomial", "poisson"))

    fit_y <- glm(y_form, data = df, family = y_fam)
  }

  #-- RETURN

  list(fit_a = fit_a, wtd_fit_a = wtd_fit_a, fit_s = fit_s, fit_y = fit_y)
}

#--- 3. COMPUTE TARGETS --------------------------------------------------------

# NOTES: Pull Out Predictions and Compute Point Estimates

compute_targets <- function(models, prep, id_form) {

  df <- prep$df

  aa <- prep$A
  yy <- prep$Y

  S_known <- prep$S_known

  #-- PROPENSITY SCORES

  P_A1  <- predict(models$fit_a,     newdata = df, type = "response")
  wP_A1 <- predict(models$wtd_fit_a, newdata = df, type = "response")

  #-- SELECTION PROBABILITIES

  if (S_known) {

    P_S1_A1 <- df$P_S_cond_A1X
    P_S1_A0 <- df$P_S_cond_A0X

  } else {

    P_S1_A1 <- predict(models$fit_s, newdata = prep$df1, type = "response")
    P_S1_A0 <- predict(models$fit_s, newdata = prep$df0, type = "response")
  }

  P_S_X <- P_S1_A1 * wP_A1 + P_S1_A0 * (1 - wP_A1)

  P_S <- length(prep$S) / sum(1 / prep$S)

  #-- TARGET

  cdiff <- switch(id_form,

    OM = {

      m1 <- predict(models$fit_y, newdata = prep$df1, type = "response")
      m0 <- predict(models$fit_y, newdata = prep$df0, type = "response")

      mean((m1 - m0) / P_S_X) * P_S
    },

    IPW1 = {

      mean((aa * yy / wP_A1 / P_S1_A1) -

       ((1 - aa) * yy / (1 - wP_A1) / P_S1_A0)) * P_S
    },

    IPW2 = {

      mean(((aa * yy / P_A1) - ((1 - aa) * yy / (1 - P_A1))) / P_S_X) * P_S
    },

    DR = {

      m1 <- predict(models$fit_y, newdata = prep$df1, type = "response")
      m0 <- predict(models$fit_y, newdata = prep$df0, type = "response")

      mu1 <- m1 + aa * (yy - m1) / P_A1
      mu0 <- m0 + (1 - aa) * (yy - m0) / (1 - P_A1)

      mean((mu1 - mu0) / P_S_X) * P_S
    },

    stop("invalid id_form")
  )

  #-- RETURN ALL ESTIMATES

  coefs <- switch(id_form,

    OM = c(coef(models$wtd_fit_a),

      if (!S_known) coef(models$fit_s) else NULL,

      coef(models$fit_y)),

    IPW1 = c(coef(models$wtd_fit_a),

      if (!S_known) coef(models$fit_s) else NULL),

    IPW2 = c(coef(models$fit_a), coef(models$wtd_fit_a),

      if (!S_known) coef(models$fit_s) else NULL),

    DR = c(coef(models$fit_a), coef(models$wtd_fit_a),

      if (!S_known) coef(models$fit_s) else NULL,

      coef(models$fit_y)),

    stop("invalid id_form")
  )

  c(coefs, cdiff = cdiff)
}

#--- 4. MAKE U CLOSURE ---------------------------------------------------------

# NOTES: Construct U() Closure for All `id_form`

make_uclosure <- function(prep, models, id_form, y_fam = NULL) {

  with(prep, {

    #-- COEFFICIENT BLOCK LENGTHS

    coeff_lengths <- switch(id_form,

      OM = c(

        nw = length(coef(models$wtd_fit_a)),

        ns = if(!S_known) length(coef(models$fit_s)) else 0,

        ny = length(coef(models$fit_y))),

      IPW1 = c(

        nw = length(coef(models$wtd_fit_a)),

        ns = if(!S_known) length(coef(models$fit_s)) else 0),

      IPW2 = c(

        na = length(coef(models$fit_a)),

        nw = length(coef(models$wtd_fit_a)),

        ns = if(!S_known) length(coef(models$fit_s)) else 0),

      DR = c(

        na = length(coef(models$fit_a)),

        nw = length(coef(models$wtd_fit_a)),

        ns = if(!S_known) length(coef(models$fit_s)) else 0,

        ny = length(coef(models$fit_y)))
    )

    #-- INDICES

    idx <- cumsum(c(0, coeff_lengths, 1))

    #-- SCORE EQUATIONS

    function(theta) {

      #- Split `theta` into `blocks` and `CDIFF`

      params <- theta[1:sum(coeff_lengths)]

      CDIFF <- theta[length(theta)]

      #- Extract Parameter Blocks

      pos <- split(params, rep(names(coeff_lengths), coeff_lengths))

      #- Helpers

      aa <- A; yy <- Y; ss <- S

      split_sel <- function(v) list(beta = v[-length(v)], phi = v[length(v)])

      #- Compute Selection Scores

      get_sel_scores <- function(beta, phi) {

        P_S_AX <- plogis(X_s %*% beta)

        ss_star <- qlogis(ss)

        mu_star <- digamma(P_S_AX * phi) - digamma((1 - P_S_AX) * phi)

        W <- P_S_AX * (1 - P_S_AX)

        U_S <- phi * sweep(X_s, 1, W, `*`) * c(ss_star - mu_star)

        U_Phi <- (P_S_AX * (ss_star - mu_star) + log(1 - ss) -

          digamma((1 - P_S_AX) * phi) + digamma(phi))

        list(U_S = U_S, U_Phi = U_Phi)
      }

      #- Score Equations by `id_form`

      if (id_form == "OM") {

        tau_w <- pos$nw;

        gamma <- pos$ny

        wP <- plogis(X_a %*% tau_w)

        U_wtdA <- X_a * c((aa - wP) / ss)

        if (!S_known) {

          sel <- split_sel(pos$ns)

          sel_scores <- get_sel_scores(sel$beta, sel$phi)
        }

        mu_Y <- switch(y_fam,

          gaussian = X_y %*% gamma,

          binomial = plogis(X_y %*% gamma),

          poisson = exp(X_y %*% gamma))

        U_Y  <- X_y * c(yy - mu_Y)

        m1 <- switch(y_fam,

          gaussian = X_y1 %*% gamma,

          binomial = plogis(X_y1 %*% gamma),

          poisson = exp(X_y1 %*% gamma))

        m0 <- switch(y_fam,

          gaussian = X_y0 %*% gamma,

          binomial = plogis(X_y0 %*% gamma),

          poisson = exp(X_y0 %*% gamma))

        PS1A1 <- if(S_known) P_S_cond_A1X else plogis(X_s1 %*% sel$beta)

        PS1A0 <- if(S_known) P_S_cond_A0X else plogis(X_s0 %*% sel$beta)

        P_S_X  <- PS1A1 * wP + PS1A0 * (1 - wP)

        P_S <- length(S) / sum(1 / S)

        U_C <- CDIFF - ((m1 - m0) / P_S_X) * P_S

        if (S_known) {

          cbind(U_wtdA, U_Y, U_C)

        } else {

          cbind(U_wtdA, sel_scores$U_S, sel_scores$U_Phi, U_Y, U_C)
        }

      } else if (id_form == "IPW1") {

        tau_w <- pos$nw

        wP <- plogis(X_a %*% tau_w)

        U_wtdA <- X_a * c((aa - wP) / ss)

        if (!S_known) {

          sel <- split_sel(pos$ns)

          sel_scores <- get_sel_scores(sel$beta, sel$phi)
        }

        PS1A1 <- if(S_known) P_S_cond_A1X else plogis(X_s1 %*% sel$beta)

        PS1A0 <- if(S_known) P_S_cond_A0X else plogis(X_s0 %*% sel$beta)

        P_S <- length(S) / sum(1 / S)

        U_C <- CDIFF - ((aa * yy / wP / PS1A1) -

          ((1 - aa) * yy / (1 - wP) / PS1A0)) * P_S

        if (S_known) {

          cbind(U_wtdA, U_C)

        } else {

          cbind(U_wtdA, sel_scores$U_S, sel_scores$U_Phi, U_C)
        }

      } else if (id_form == "IPW2") {

        tau <- pos$na

        tau_w <- pos$nw

        P <- plogis(X_a %*% tau)

        U_A <- X_a * c(aa - P)

        wP <- plogis(X_a %*% tau_w)

        U_wtdA <- X_a * c((aa - wP) / ss)

        if (!S_known) {

          sel <- split_sel(pos$ns)

          sel_scores <- get_sel_scores(sel$beta, sel$phi)
        }

        PS1A1 <- if(S_known) P_S_cond_A1X else plogis(X_s1 %*% sel$beta)

        PS1A0 <- if(S_known) P_S_cond_A0X else plogis(X_s0 %*% sel$beta)

        P_S_X <- PS1A1 * wP + PS1A0 * (1 - wP)

        P_S <- length(S) / sum(1 / S)

        U_C <- CDIFF - ((aa * yy / P) -

          ((1 - aa) * yy / (1 - P))) / P_S_X * P_S

        if (S_known) {

          cbind(U_A, U_wtdA, U_C)

        } else {

          cbind(U_A, U_wtdA, sel_scores$U_S, sel_scores$U_Phi, U_C)
        }

      } else if (id_form == "DR") {

        tau <- pos$na

        tau_w <- pos$nw

        gamma <- pos$ny

        P <- plogis(X_a %*% tau)

        U_A <- X_a * c(aa - P)

        wP <- plogis(X_a %*% tau_w)

        U_wtdA <- X_a * c((aa - wP) / ss)

        if (!S_known) {

          sel <- split_sel(pos$ns)

          sel_scores <- get_sel_scores(sel$beta, sel$phi)
        }

        mu_Y <- switch(y_fam,

          gaussian = X_y %*% gamma,

          binomial = plogis(X_y %*% gamma),

          poisson = exp(X_y %*% gamma))

        U_Y <- X_y * c(yy - mu_Y)

        m1 <- switch(y_fam,

          gaussian = X_y1 %*% gamma,

          binomial = plogis(X_y1 %*% gamma),

          poisson = exp(X_y1 %*% gamma))

        m0 <- switch(y_fam,

          gaussian = X_y0 %*% gamma,

          binomial = plogis(X_y0 %*% gamma),

          poisson = exp(X_y0 %*% gamma))

        mu1 <- m1 + aa * (yy - m1) / P

        mu0 <- m0 + (1 - aa) * (yy - m0) / (1 - P)

        PS1A1 <- if(S_known) P_S_cond_A1X else plogis(X_s1 %*% sel$beta)

        PS1A0 <- if(S_known) P_S_cond_A0X else plogis(X_s0 %*% sel$beta)

        P_S_X <- PS1A1 * wP + PS1A0 * (1 - wP)

        P_S <- length(S) / sum(1 / S)

        U_C <- CDIFF - (mu1 - mu0) / P_S_X * P_S

        if (S_known) {

          cbind(U_A, U_wtdA, U_Y, U_C)

        } else {

          cbind(U_A, U_wtdA, sel_scores$U_S, sel_scores$U_Phi, U_Y, U_C)
        }

      } else {

        stop("invalid id_form")
      }
    }
  })
}

#--- 5. COMPUTE VARIANCE-COVARIANCE MATRIX -------------------------------------

# NOTES: Numerical Derivative + Cluster-Robust

compute_vcov <- function(u_func, est, prep, strata = NULL, cluster = NULL) {

  #-- DERIVATIVES

  ests <- as.numeric(est)

  halfmeat <- u_func(ests)

  bread <- numDeriv::jacobian(func = function(t) colSums(u_func(t)), x = ests)

  #-- INFLUENCE FUNCTION

  IF <- tcrossprod(halfmeat, solve(-bread))

  #-- VARIANCE-COVARIANCE MATRIX

  n <- nrow(prep$df)
  r <- ncol(IF)

  st_vec <- if (!is.null(strata))  prep$df[[strata]]  else factor(rep(1, n))
  cl_vec <- if (!is.null(cluster)) prep$df[[cluster]] else factor(seq_len(n))

  cl_id <- interaction(st_vec, cl_vec, drop = TRUE)

  U_clust <- rowsum(IF, group = cl_id)

  strata_cl <- sapply(strsplit(rownames(U_clust), "\\."), `[`, 1)

  V_mat <- matrix(0, r, r)

  for (k in unique(strata_cl)) {

    idx <- which(strata_cl == k)
    U_k <- U_clust[idx, , drop=FALSE]
    M_k <- nrow(U_k)

    if (M_k > 1) {

      Ubar   <- colMeans(U_k)
      D      <- sweep(U_k, 2, Ubar, "-")
      V_k    <- crossprod(D)
      V_mat  <- V_mat + (M_k/(M_k - 1)) * V_k
    }
  }

  V_mat
}
