#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 5a_model_validation.R
#
# Comprehensive model validation and comparison
# Performs LOO-CV comparisons, posterior predictive checks, and parameter
# comparison across fitted models.
#
# Input: results/c2_run_20260626/model_output/<model>/{posterior_draws.rds, loo.rds,
#        diagnostics.rds} (widened rds bundle, post Phase 5 W1/W2)
# Output: results/loo_* rds, results/ppc_*.pdf, results/model_fit_metrics.csv
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(posterior)
library(loo)
library(bayesplot)

source("scripts/posterior_helpers.R")

# Null-coalesce helper
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("MODEL VALIDATION AND COMPARISON (Standalone)\n")
cat("===========================================\n\n")

# Find all fitted models by looking at the local April mirror.
model_dirs <- list.dirs(APRIL_RUN, recursive = FALSE, full.names = FALSE)
model_dirs <- sort(model_dirs[!grepl("^_", model_dirs)])
model_dirs <- model_dirs[sapply(model_dirs, function(m)
  file.exists(file.path(APRIL_RUN, m, "posterior_draws.rds")))]

cat("Found", length(model_dirs), "models with widened draws:\n")
for (m in model_dirs) cat("  -", m, "\n")

# 1. LOAD LOO RESULTS ---------------------------------------------------------
loo_results <- list()
cat("\nLoading loo.rds for each model...\n")
for (model_name in model_dirs) {
  lo <- tryCatch(load_loo(model_name),
                 error = function(e) { cat("  ", model_name, ":", conditionMessage(e), "\n"); NULL })
  if (!is.null(lo)) {
    loo_results[[model_name]] <- lo
    n_bad <- sum(lo$diagnostics$pareto_k > 0.7)
    cat(sprintf("  %-25s  elpd=%7.1f  p_loo=%5.1f  n_k>0.7=%d\n",
                model_name,
                lo$estimates["elpd_loo", "Estimate"],
                lo$estimates["p_loo", "Estimate"],
                n_bad))
  }
}

dir.create("results", recursive = TRUE, showWarnings = FALSE)

# LOO comparison across models
if (length(loo_results) > 1) {
  cat("\n\nLOO-CV Model Comparison (best to worst):\n")
  loo_comp <- loo_compare(loo_results)
  print(loo_comp)
  saveRDS(loo_results, "results/loo_results.rds")
  saveRDS(loo_comp,    "results/loo_comparison.rds")
}

# 2. POSTERIOR PREDICTIVE CHECKS ---------------------------------------------
cat("\n\nPerforming posterior predictive checks...\n")
for (model_name in names(loo_results)) {
  cat("  PPC for", model_name, "\n")
  stan_data <- load_stan_data(model_name)
  draws     <- load_draws(model_name)
  y_obs <- stan_data$d2H_wax
  if (is.null(y_obs)) { cat("    No d2H_wax in stan_data\n"); next }

  y_rep <- as_draws_matrix(subset_draws(draws, variable = "d2H_rep"))
  if (nrow(y_rep) == 0) next

  pdf(paste0("results/ppc_", model_name, ".pdf"), width = 10, height = 6)
  print(ppc_dens_overlay(y_obs, y_rep[1:min(100, nrow(y_rep)), , drop = FALSE]))
  print(ppc_scatter_avg(y_obs, y_rep))
  dev.off()
}

# 3. KEY PARAMETER COMPARISON ------------------------------------------------
cat("\n\nComparing key parameters across models (back-transformed to ‰)...\n")

param_comparison <- map_df(names(loo_results), function(model_name) {
  draws     <- load_draws(model_name)
  stan_data <- load_stan_data(model_name)

  d2h_mean <- stan_data$d2H_wax_mean_original %||% stan_data$scaling_params$d2H_mean
  d2h_sd   <- stan_data$d2H_wax_sd_original   %||% stan_data$scaling_params$d2H_sd

  get_summary <- function(var) {
    d <- as.numeric(as_draws_matrix(subset_draws(draws, variable = var)))
    list(mean = mean(d), sd = sd(d))
  }
  b0    <- get_summary("beta_0")
  boipc <- get_summary("beta_oipc")
  sg    <- get_summary("sigma")

  data.frame(
    model        = model_name,
    intercept    = b0$mean * d2h_sd + d2h_mean,
    intercept_sd = b0$sd   * d2h_sd,
    slope        = boipc$mean,
    slope_sd     = boipc$sd,
    sigma        = sg$mean * d2h_sd,
    sigma_sd     = sg$sd   * d2h_sd,
    stringsAsFactors = FALSE
  )
})

write.csv(param_comparison, "results/parameter_comparison.csv", row.names = FALSE)

# 4. MODEL FIT METRICS --------------------------------------------------------
cat("\n\nComputing model fit metrics...\n")

fit_metrics_list <- list()
for (model_name in names(loo_results)) {
  res <- tryCatch({
    draws     <- load_draws(model_name)
    stan_data <- load_stan_data(model_name)
    y_obs <- stan_data$d2H_wax
    if (is.null(y_obs) || length(y_obs) == 0) stop("no d2H_wax")

    mu_mat <- as_draws_matrix(subset_draws(draws, variable = "mu"))
    y_pred <- colMeans(mu_mat)
    if (length(y_pred) != length(y_obs)) stop("length mismatch")

    ss_res <- sum((y_obs - y_pred)^2)
    ss_tot <- sum((y_obs - mean(y_obs))^2)
    r_squared <- 1 - ss_res / ss_tot
    rmse <- sqrt(mean((y_obs - y_pred)^2))

    loo_obj <- loo_results[[model_name]]
    data.frame(
      model = model_name,
      r_squared = r_squared,
      rmse = rmse,
      looic = loo_obj$estimates["looic", "Estimate"],
      elpd  = loo_obj$estimates["elpd_loo", "Estimate"],
      p_loo = loo_obj$estimates["p_loo", "Estimate"],
      n_obs = length(y_obs),
      n_bad_k = sum(loo_obj$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat("  Error with", model_name, ":", e$message, "\n"); NULL
  })
  if (!is.null(res)) fit_metrics_list[[model_name]] <- res
}

if (length(fit_metrics_list) > 0) {
  fit_metrics <- bind_rows(fit_metrics_list) %>% arrange(looic)
  write.csv(fit_metrics, "results/model_fit_metrics.csv", row.names = FALSE)
  cat("\nModel fit metrics saved for", nrow(fit_metrics), "models\n")
} else {
  cat("\nERROR: No fit metrics could be calculated\n")
}

cat("\n\nModel comparison complete! Results saved in results/\n")
