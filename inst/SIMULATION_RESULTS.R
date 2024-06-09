#===============================================================================
#
#  PROGRAM: SIMULATION_RESULTS.R
#
#  AUTHOR:  Stephen Salerno (ssalerno@fredhutch.org)
#
#  PURPOSE: Summarize simulation study results for manuscript.
#
#  INPUT:   SIMULATION_RESULTS.RData - Results from main simulation study
#
#           An .RData file containing the object `results`, which is a list of
#           eight elements corresponding to the eight simulation settings
#           described in the manuscript. Each list element contains a numeric
#           matrix of size `nsims` x `37`, where each row corresponds to a
#           simulation replicate and each column is a saved metric.
#
#           SENSITIVITY_RESULTS.RData - Results from sensitivity analysis
#
#           An .RData file containing the object `results`, which is a list of
#           five elements corresponding to the five sensitivity analysis
#           settings described in the manuscript. Each list element contains a
#           numeric matrix of size `nsims` x `16`, where each row corresponds
#           to a simulation replicate and each column is a saved metric.
#
#  OUTPUT:  diagram.png - Conceptual diagrams for methods/data examples
#
#           simulation_results.html - Table of summarized simulation results
#
#           simulation results.png - Figure of summarized simulation results
#
#           sensitivity_results.html - Table of sensitivity analysis results
#
#  UPDATED: 2024-05-25
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- NECESSARY PACKAGES --------------------------------------------------------

library(pacman)

p_load(here, gt, ggsci, patchwork, ggraph, tidygraph, tidyverse, update = F)

#--- SIMULATION PARAMETERS -----------------------------------------------------

nsims <- 200; N <- 100000

tau_X_vec <- beta_A_vec <- beta_X_vec <- 0:1

settings <- expand.grid(tau_X_vec, beta_A_vec, beta_X_vec)

names(settings) <- c("tau_X", "beta_A", "beta_X")

mthds <- c("Oracle Estimator", "Simple Regression", "Multiple Regression",

  "IPTW Estimator", "Survey-Weighted Multiple Regression",

  "IPTW Multiple Regression", "IPTW + Survey-Weighted Multiple Regression",

  "Weighted IPTW + Survey-Weighted Multiple Regression",

  "Naively Weighted G-Computation",

  "Outcome Modeling and Direct Standardization",

  "Inverse Probability Weighting 1", "Inverse Probability Weighting 2")

metrics <- c("Setting", "Method", "True ATE", "Estimated ATE", "Bias",

  "Analytic SE", "Monte Carlo SE", "MSE", "Coverage")

#=== RESULTS ===================================================================

#--- SUMMARIZE RESULTS ---------------------------------------------------------

load(here("inst", "SIMULATION_RESULTS.RData"))

bias_cols <- seq(2, 35, 3)

err_cols  <- seq(3, 36, 3)

cov_cols  <- seq(4, 37, 3)

results_summary <- lapply(1:8, function(i) {

  ate   <- mean(results[[i]][, 1])                        #-- TRUE ATE

  est   <- colMeans(results[[i]][, bias_cols])            #-- ESTIMATED ATE

  bias  <- colMeans((results[[i]][, bias_cols] -          #-- BIAS

    results[[i]][, 1]) / results[[i]][, 1])

  anerr <- colMeans(results[[i]][, err_cols])             #-- ANALYTIC SE

  mcerr <- apply(results[[i]][, bias_cols], 2, sd)        #-- MONTE CARLO SE

  mse   <- colMeans((results[[i]][, bias_cols] -          #-- MSE

    results[[i]][, 1])^2)

  cov   <- colMeans(results[[i]][, cov_cols])             #-- COVERAGE

  return(data.frame(i = rep(i, 12),                       #-- COMBINE

    ate = rep(ate, 12), est = est,

    bias = bias, anerr = anerr,

    mcerr = mcerr, mse = mse, cov = cov))
})

results_summary_all <- do.call(rbind, results_summary) |>

  as_tibble() |>

  mutate(method = rep(mthds, 8) |> factor(levels = mthds)) |>

  select(i, method, everything())

names(results_summary_all) <- metrics

results_summary_all

#---TABLE OF ALL RESULTS -------------------------------------------------------

sim_all_gt <- results_summary_all |>

  mutate(Setting = paste("Setting", Setting)) |>

  gt(id = "one", groupname_col = "Setting") |>

  tab_header(

    title = "Full Simulation Results for Each Method (nsims = 200)",

    subtitle = "across different Data Generating Mechanism Settings") |>

  cols_align(align = "left", columns = `Method`) |>

  fmt_number(`True ATE`:`Coverage`, decimals = 3) |>

  tab_style(

    style = cell_text(weight = "bold"),

    locations = cells_row_groups()) |>

  tab_style(

    style = cell_text(weight = "bold"),

    locations = cells_column_labels())

sim_all_gt

gtsave(sim_all_gt, here("inst", "simulation_results.html"))

gt_labs <- apply(settings, 1, function(i) {

  paste(c("\U1D70F\U2093 =", "\U03B2\U2090 =", "\U03B2\U2093 ="), i) }) |>

  t() |> data.frame() |> unite(col = "labs", sep = ", ") |> pull(labs)

gt_labs

facet_labs <- apply(settings, 1, function(i) {

  paste(c("tau[X] ==", "beta[A] ==", "beta[X] =="), i) }) |> t() |>

  data.frame() |> unite(col = "labs", sep = "~") |> pull(labs)

facet_labs

#=== FIGURES ===================================================================

#--- SIMULATION SETTING DAGS ---------------------------------------------------

dags <- tibble(

  from = c(rep(c("X", "A"), each=8), rep(c("X", "A"), each=4), rep("X", 4)),

  to = c(rep("Y", 16), rep("S", 8), rep("A", 4)),

  setting = c(rep(1:8, 2), 5:8, 3:4, 7:8, seq(2, 8, 2))) |>

  as_tbl_graph()

dag_plt <- ggraph(dags, layout = "manual", x = c(0,1,2,1), y = c(1,0,1,2)) +

  theme_void() +

  geom_edge_link(aes(color = factor(to)), width = 1.1,

    start_cap = circle(5, 'mm'), end_cap = circle(5, 'mm')) +

  geom_node_point(aes(color = name), size = 10, alpha = 0.25, stroke = 2) +

  geom_node_text(aes(label = name, color = name), fontface = "bold") +

  facet_edges(~ setting, ncol = 8) +

  scale_edge_color_manual(values = c("#2F65A7", "#75988d", "#9A3324")) +

  scale_color_manual(values = c("#2F65A7", "#9A3324", "#989C97", "#75988d")) +

  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +

  theme(

    legend.position = "none",

    strip.background = element_blank(),

    strip.text = element_blank())

dag_plt

#--- CONCEPTUAL DIAGRAMS -------------------------------------------------------

dag_8 <- tibble(

  from = c("X", "A", "X", "A", "X"),

  to   = c("Y", "Y", "S", "S", "A"),

  setting = 8) |>

  as_tbl_graph()

dag_acd <- tibble(

  from = c("SES", "Race", "SES", "Race", "SES"),

  to   = c("Telomere\nLength", "Telomere\nLength", "Sample\nSelection",

    "Sample\nSelection", "Race"),

  setting = 9) |>

  as_tbl_graph()

dag_ate <- tibble(

  from = c("SES", "Environmental\nExposure", "SES", "Environmental\nExposure",

    "SES"),

  to   = c("Telomere\nLength", "Telomere\nLength", "Sample\nSelection",

    "Sample\nSelection", "Environmental\nExposure"),

  setting = 10) |>

  as_tbl_graph()

dag_plt_8 <- ggraph(dag_8, layout = "manual", x = c(0,1,2,1), y = c(1,0,1,2)) +

  theme_void() +

  geom_edge_link(aes(color = factor(to)), width = 1.1,

    start_cap = circle(17.5, 'mm'), end_cap = circle(17.5, 'mm')) +

  geom_node_point(aes(color = name), size = 40, alpha = 0.25, stroke = 2) +

  geom_node_text(aes(label = name, color = name), fontface = "bold") +

  facet_edges(~ setting, ncol = 8) +

  scale_edge_color_manual(values = c("#2F65A7", "#75988d", "#9A3324")) +

  scale_color_manual(values = c("#2F65A7", "#9A3324", "#989C97", "#75988d")) +

  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +

  theme(

    legend.position = "none",

    strip.background = element_blank(),

    strip.text = element_blank())

dag_plt_8

dag_plt_acd <- ggraph(dag_acd, layout = "manual",

  x = c(0,1,2,1), y = c(1,0,1,2)) +

  theme_void() +

  geom_edge_link(aes(color = factor(to)), width = 1.1,

    start_cap = circle(17.5, 'mm'), end_cap = circle(17.5, 'mm')) +

  geom_node_point(aes(color = name), size = 40, alpha = 0.25, stroke = 2) +

  geom_node_text(aes(label = name, color = name), fontface = "bold") +

  facet_edges(~ setting, ncol = 8) +

  scale_edge_color_manual(values = c("#2F65A7", "#75988d", "#9A3324")) +

  scale_color_manual(values = c("#2F65A7", "#9A3324", "#989C97", "#75988d")) +

  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +

  theme(

    legend.position = "none",

    strip.background = element_blank(),

    strip.text = element_blank())

dag_plt_acd

dag_plt_ate <- ggraph(dag_ate, layout = "manual",

    x = c(0,1,2,1), y = c(1,0,1,2)) +

  theme_void() +

  geom_edge_link(aes(color = factor(to)), width = 1.1,

    arrow = arrow(type = "closed", length = unit(5, 'mm')),

    start_cap = circle(17.5, 'mm'), end_cap = circle(17.5, 'mm')) +

  geom_node_point(aes(color = name), size = 40, alpha = 0.25, stroke = 2) +

  geom_node_text(aes(label = name, color = name), fontface = "bold") +

  facet_edges(~ setting, ncol = 8) +

  scale_edge_color_manual(values = c("#2F65A7", "#75988d", "#9A3324")) +

  scale_color_manual(values = c("#2F65A7", "#9A3324", "#989C97", "#75988d")) +

  expand_limits(x = c(-0.5, 2.5), y = c(-0.5, 2.5)) +

  theme(

    legend.position = "none",

    strip.background = element_blank(),

    strip.text = element_blank())

dag_plt_ate

diagram <- dag_plt_8 + dag_plt_acd + dag_plt_ate +

  plot_annotation(tag_levels = 'A')

ggsave(here("inst", "diagram.png"), diagram, height = 5, width = 15)

#--- SIMULATION SETTING TABLE --------------------------------------------------

rect_dat <- data.frame(

  xstr = seq(0.5, 3.5, 1), xend = seq(1.5, 4.5, 1), col = paste(1:4))

simset_tab <- settings |> as_tibble() |>

  mutate(

    `Propensity Model` = paste0("tau[0] == -1~tau[X] == ", tau_X),

    `Selection Model` = paste0("beta[0] == -4.5~beta[A] == ",

      beta_A, "~beta[X] == ", beta_X),

    `Outcome Model` =

      "gamma[0] == 1~gamma[A] == 1~gamma[X] == 1~gamma[AX] == 0.1",

    `Bias Expected` = c("None", "Confounding", "None", "Confounding",

        "Selection", "Confounding + Selection", "Selection",

        "Confounding + Selection")) |>

  select(-(tau_X:beta_X)) |>

  pivot_longer(everything()) |>

  mutate(

    name = factor(name) |>

      fct_relevel("Propensity Model", "Selection Model",

        "Outcome Model", "Bias Expected"),

    y = rep(1, 32),

    grp = rep(1:8, each = 4))

simset_plt <- ggplot(data = simset_tab) + theme_bw() + coord_flip() +

  facet_grid(~ grp, labeller = labeller(grp = ~ paste("Setting", .x))) +

  scale_x_discrete(1:4, expand = c(0,0), limits = rev) +

  scale_fill_manual(values = c("gray", "#75988d", "#9A3324", "#2F65A7")) +

  geom_rect(data = rect_dat, alpha = 0.25,

    aes(ymin = -Inf, ymax = Inf, xmin = xstr, xmax = xend, fill = col)) +

  geom_text(data = simset_tab, aes(x = name, y = y, label = value),

    parse = T, size = 3) +

  geom_vline(xintercept = c(1.5, 2.5, 3.5)) +

  theme(

    legend.position = "none",

    strip.background = element_blank(),

    strip.text = element_text(face = "bold", size = 12),

    axis.text = element_text(face = "bold", size = 10),

    panel.grid = element_blank(),

    axis.title = element_blank(),

    axis.text.x = element_blank(),

    axis.ticks.x = element_blank())

simset_plt

#--- PLOTS OF METRICS ----------------------------------------------------------

simbias_plt <- results_summary_all |>

  select(Setting, Method, Bias, `Monte Carlo SE`) |>

  rename(se = `Monte Carlo SE`) |>

  mutate(grp = "Bias") |>

  ggplot(aes(x = Bias, y = Method, xmin = Bias - se, xmax = Bias + se)) +

    theme_bw() +

    facet_grid(grp ~ Setting) +

    geom_vline(xintercept = 0, lty = 2) +

    geom_errorbar(color = '#00274C', width = 0.1, linewidth = 1.1) +

    geom_point(color = '#00274C', size = 3) +

    scale_x_continuous(labels = scales::percent_format()) +

    scale_y_discrete(limits = rev) +

    labs(x = "\nPercent Bias", y = NULL) +

    theme(

      strip.background = element_blank(),

      strip.text.x = element_blank(),

      strip.text.y = element_text(face = "bold", size = 12),

      axis.title = element_text(face = "bold", size = 12),

      axis.text = element_text(face = "bold", size = 10))

simbias_plt

simcov_plt <- results_summary_all |>

  select(Setting, Method, Coverage) |>

  mutate(grp = "Coverage") |>

  ggplot(aes(x = Coverage, y = Method)) +

    theme_bw() +

    facet_grid(grp ~ Setting) +

    geom_vline(xintercept = 0.95, lty = 2) +

    geom_point(color = '#00274C', size = 3) +

    scale_x_continuous(labels = scales::percent_format()) +

    scale_y_discrete(limits = rev) +

    labs(x = "\nCoverage", y = NULL) +

    theme(

      strip.background = element_blank(),

      strip.text.x = element_blank(),

      strip.text.y = element_text(face = "bold", size = 12),

      axis.title = element_text(face = "bold", size = 12),

      axis.text = element_text(face = "bold", size = 10))

simcov_plt

#--- COMBINE ALL PLOTS ---------------------------------------------------------

sim_plot <- dag_plt / simset_plt / simbias_plt / simcov_plt +

  plot_layout(heights = c(1, 1, 2, 2)) +

  plot_annotation(tag_levels = 'A')

sim_plot

ggsave(here("inst", "simulation_results.png"), height = 10, width = 20)

#=== SENSITIVITY ANALYSIS ======================================================

load(here("inst", "SENSITIVITY_RESULTS.RData"))

nsims <- 200; N <- 100000

mthds_sens <- c("Oracle Estimator",

  "Outcome Modeling and Direct Standardization",

  "Outcome Modeling and Direct Standardization, Misspecified",

  "Inverse Probability Weighting 1",

  "Inverse Probability Weighting 2")

metrics_sens <- c("Setting", "Method", "True ATE", "Estimated ATE", "Bias",

  "Analytic SE", "Monte Carlo SE", "MSE", "Coverage")

bias_cols <- seq(2, 14, 3)

err_cols  <- seq(3, 15, 3)

cov_cols  <- seq(4, 16, 3)

results_sens <- lapply(1:5, function(i) {

  ate   <- mean(results[[i]][, 1])                        #-- TRUE ATE

  est   <- colMeans(results[[i]][, bias_cols])            #-- ESTIMATED ATE

  bias  <- colMeans((results[[i]][, bias_cols] -          #-- BIAS

    results[[i]][, 1]) / results[[i]][, 1])

  anerr <- colMeans(results[[i]][, err_cols])             #-- ANALYTIC SE

  mcerr <- apply(results[[i]][, bias_cols], 2, sd)        #-- MONTE CARLO SE

  mse   <- colMeans((results[[i]][, bias_cols] -          #-- MSE

    results[[i]][, 1])^2)

  cov   <- colMeans(results[[i]][, cov_cols])             #-- COVERAGE

  return(data.frame(i = i,                                #-- COMBINE

    ate = ate, est = est,

    bias = bias, anerr = anerr,

    mcerr = mcerr, mse = mse, cov = cov))
})

results_sens_all <- do.call(rbind, results_sens) |>

  as_tibble() |>

  mutate(method = rep(mthds_sens, 5) |> factor(levels = mthds_sens)) |>

  select(i, method, everything())

names(results_sens_all) <- metrics_sens

results_sens_all

sens_all_gt <- results_sens_all |>

  mutate(Setting = paste("Setting", Setting)) |>

  gt(id = "one", groupname_col = "Setting") |>

  tab_header(

    title = "Full Simulation Results for Each Method (nsims = 200)",

    subtitle = "across different treatment heterogeneity levels") |>

  cols_align(align = "left", columns = `Method`) |>

  fmt_number(`True ATE`:`Coverage`, decimals = 3) |>

  tab_style(

    style = cell_text(weight = "bold"),

    locations = cells_row_groups()) |>

  tab_style(

    style = cell_text(weight = "bold"),

    locations = cells_column_labels())

sens_all_gt

gtsave(sens_all_gt, here("inst", "sensitivity_results.html"))

#=== END =======================================================================
