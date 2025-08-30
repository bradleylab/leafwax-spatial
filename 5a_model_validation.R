#!/usr/bin/env Rscript
# 5a_model_comparison.R - Comprehensive model comparison

library(tidyverse)
library(cmdstanr)
library(posterior)
library(loo)
library(bayesplot)

# Set data directory
DATA_DIR <- "prepared_data_consolidated"

cat("MODEL VALIDATION AND COMPARISON (Standalone)\n")
cat("===========================================\n\n")

# Find all fitted models
model_dirs <- list.dirs("model_output", recursive = FALSE, full.names = FALSE)
model_dirs <- model_dirs[model_dirs != ""]
model_dirs <- sort(model_dirs)

cat("Found", length(model_dirs), "completed models:\n")
for (m in model_dirs) {
  cat("  -", m, "\n")
}

# Initialize storage
loo_results <- list()
model_summaries <- list()

cat("\nComputing LOO-CV for all models...\n")

# Process each model
for (model_name in model_dirs) {
  cat("  Processing", model_name, "... ")
  
  fit_file <- file.path("model_output", model_name, "fit.rds")
  loo_file <- file.path("model_output", model_name, "loo.rds")
  stan_data_file <- file.path(DATA_DIR, paste0("stan_data_", model_name, ".rds"))
  
  if (!file.exists(fit_file) || !file.exists(stan_data_file)) {
    cat("SKIPPED (missing files)\n")
    next
  }
  
  # Check for existing LOO
  if (file.exists(loo_file)) {
    loo_results[[model_name]] <- readRDS(loo_file)
    cat("(loaded existing)\n")
  } else {
    # Load fit and compute LOO
    fit <- readRDS(fit_file)
    
    tryCatch({
      log_lik <- fit$draws("log_lik", format = "array")
      loo_results[[model_name]] <- loo(log_lik, cores = 2)
      
      # Save for future use
      saveRDS(loo_results[[model_name]], loo_file)
      
      # Check for problematic observations
      pareto_k <- loo_results[[model_name]]$diagnostics$pareto_k
      n_bad <- sum(pareto_k > 0.7)
      if (n_bad > 0) {
        cat("WARNING:", n_bad, "observations with k > 0.7\n")
      } else {
        cat("OK\n")
      }
      
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
    })
  }
}

# LOO comparison
if (length(loo_results) > 1) {
  cat("\n\nLOO-CV Model Comparison:\n")
  loo_comp <- loo_compare(loo_results)
  print(loo_comp)
  
  # Save
  saveRDS(loo_results, "results/loo_results.rds")
  saveRDS(loo_comp, "results/loo_comparison.rds")
}

# Posterior predictive checks
cat("\n\nPerforming posterior predictive checks...\n")

for (model_name in names(loo_results)) {
  cat("  PPC for", model_name, "\n")
  
  fit <- readRDS(file.path("model_output", model_name, "fit.rds"))
  stan_data <- readRDS(file.path(DATA_DIR, paste0("stan_data_", model_name, ".rds")))
  
  # Get observed values
  y_obs <- stan_data$d2H_wax
  
  if (is.null(y_obs)) {
    cat("    Warning: Cannot find observed values for PPC\n")
    next
  }
  
  # Extract posterior predictions
  y_rep <- fit$draws("d2H_rep", format = "matrix")
  
  if (!is.null(y_rep) && nrow(y_rep) > 0) {
    # Create PPC plot
    pdf(paste0("results/ppc_", model_name, ".pdf"), width = 10, height = 6)
    
    # Density overlay
    print(ppc_dens_overlay(y_obs, y_rep[1:min(100, nrow(y_rep)), ]))
    
    # Scatter plot
    print(ppc_scatter_avg(y_obs, y_rep))
    
    dev.off()
  }
}

# Compare key parameters
cat("\n\nComparing key parameters across models...\n")

param_comparison <- map_df(names(loo_results), function(model_name) {
  fit <- readRDS(file.path("model_output", model_name, "fit.rds"))
  stan_data <- readRDS(file.path(DATA_DIR, paste0("stan_data_", model_name, ".rds")))
  
  # Get scaling parameters
  d2h_mean <- stan_data$d2H_wax_mean_original
  d2h_sd <- stan_data$d2H_wax_sd_original
  
  # Extract parameters
  params <- fit$summary(variables = c("beta_0", "beta_oipc", "sigma"))
  
  data.frame(
    model = model_name,
    intercept = params[params$variable == "beta_0", "mean"] * d2h_sd + d2h_mean,
    intercept_sd = params[params$variable == "beta_0", "sd"] * d2h_sd,
    slope = params[params$variable == "beta_oipc", "mean"],
    slope_sd = params[params$variable == "beta_oipc", "sd"],
    sigma = params[params$variable == "sigma", "mean"] * d2h_sd,
    sigma_sd = params[params$variable == "sigma", "sd"] * d2h_sd,
    stringsAsFactors = FALSE
  )
})

write.csv(param_comparison, "results/parameter_comparison.csv", row.names = FALSE)

# Model fit metrics
cat("\n\nComputing model fit metrics...\n")

fit_metrics_list <- list()

for (model_name in names(loo_results)) {
  tryCatch({
    fit <- readRDS(file.path("model_output", model_name, "fit.rds"))
    stan_data <- readRDS(file.path(DATA_DIR, paste0("stan_data_", model_name, ".rds")))
    
    # Get observed values
    y_obs <- stan_data$d2H_wax
    
    if (is.null(y_obs) || length(y_obs) == 0) {
      cat("  Warning: Cannot find observed values for", model_name, "\n")
      next
    }
    
    # Get fitted values (mu)
    mu_summary <- fit$summary(variables = "mu")
    y_pred <- mu_summary$mean
    
    if (length(y_pred) != length(y_obs)) {
      cat("  Warning: Prediction length mismatch for", model_name, "\n")
      next
    }
    
    # Calculate R-squared
    ss_res <- sum((y_obs - y_pred)^2)
    ss_tot <- sum((y_obs - mean(y_obs))^2)
    r_squared <- 1 - ss_res/ss_tot
    
    # Get LOO info
    loo_obj <- loo_results[[model_name]]
    
    fit_metrics_list[[model_name]] <- data.frame(
      model = model_name,
      r_squared = r_squared,
      rmse = sqrt(mean((y_obs - y_pred)^2)),
      looic = loo_obj$estimates["looic", "Estimate"],
      elpd = loo_obj$estimates["elpd_loo", "Estimate"],
      p_loo = loo_obj$estimates["p_loo", "Estimate"],
      n_obs = length(y_obs),
      n_bad_k = sum(loo_obj$diagnostics$pareto_k > 0.7),
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    cat("  Error with", model_name, ":", e$message, "\n")
  })
}

# Combine and save results
if (length(fit_metrics_list) > 0) {
  fit_metrics <- bind_rows(fit_metrics_list)
  fit_metrics <- fit_metrics %>% arrange(looic)
  write.csv(fit_metrics, "results/model_fit_metrics.csv", row.names = FALSE)
  
  cat("\nModel fit metrics saved for", nrow(fit_metrics), "models\n")
} else {
  cat("\nERROR: No fit metrics could be calculated\n")
}

cat("\n\nModel comparison complete! Results saved in results/\n")