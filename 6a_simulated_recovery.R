# run_simulated_recovery.R
# Phase 3: Simulated Data Recovery Tests
#
# Generates synthetic d2H_wax data with KNOWN parameters, fits the model,
# and checks whether it recovers the truth. Three scenarios:
#   3a. Uniform slope (0.7), no spatial slope variation
#   3b. Spatially varying slope (tropics=0.5, high-lat=0.9, global mean=0.7)
#   3c. Confounding stress test (uniform slope=0.7, intercept correlated with OIPC)
#
# Usage: Rscript run_simulated_recovery.R [scenario]
#   scenario: "3a", "3b", "3c", or "all" (default)
#
# Expected runtime per scenario: ~2-3 hrs on r6i.8xlarge

library(cmdstanr)
library(posterior)
library(loo)
library(MASS)  # for mvrnorm

# In the container CmdStan lives at a fixed path; set it explicitly because
# apptainer --containall does not reliably propagate CMDSTAN to cmdstanr's
# discovery. Guarded so local runs keep their own auto-detected CmdStan.
container_cmdstan <- "/root/.cmdstan/cmdstan-2.36.0"
if (dir.exists(container_cmdstan)) set_cmdstan_path(container_cmdstan)

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
scenarios <- if (length(args) > 0 && args[1] != "all") args[1] else c("3a", "3b", "3c")

cat("=== SIMULATED DATA RECOVERY ===\n")
cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# Run from repo root

# Load prepared data and config
stan_data <- readRDS("prepared_data/stan_data_baseline_veg_sp.rds")
config <- readRDS("prepared_data/config_baseline_veg_sp.rds")

# Load existing posteriors for realistic parameter values
orig_draws <- readRDS("model_output/baseline_veg_sp/posterior_draws.rds")

# Compile the model (with fixed range thresholds)
cat("Compiling Stan model...\n")
model <- cmdstan_model("4d_leaf_wax_spatial_model.stan")
cat("Model compiled.\n\n")

# ─── Helper: compute weighted OIPC using estimated lambda ───
compute_oipc_weighted <- function(stan_data, lambda_km) {
  scales <- stan_data$distance_scales
  weights <- exp(-scales / lambda_km)
  weights <- weights / sum(weights)
  stan_data$oipc_values %*% weights
}

# ─── Helper: generate GP surface at observation locations ───
generate_gp_surface <- function(stan_data, ls_std, sigma_gp, seed = NULL) {
  # Compute Matern 3/2 kernel at knots
  n_knots <- stan_data$n_pp_knots
  N <- stan_data$N
  sqrt3 <- sqrt(3)
  
  knot_coords <- stan_data$knot_coords
  obs_coords <- stan_data$coords
  
  # Knot-knot covariance
  K_knots <- matrix(0, n_knots, n_knots)
  for (i in 1:n_knots) {
    K_knots[i, i] <- 1.0  # eta^2 = 1
    if (i < n_knots) for (j in (i+1):n_knots) {
      d <- sqrt(sum((knot_coords[i,] - knot_coords[j,])^2))
      scaled <- sqrt3 * d / ls_std
      K_knots[i, j] <- (1 + scaled) * exp(-scaled)
      K_knots[j, i] <- K_knots[i, j]
    }
  }
  K_knots <- K_knots + diag(1e-4, n_knots)
  
  # Draw knot values
  if (!is.null(seed)) set.seed(seed)
  knot_vals <- mvrnorm(1, mu = rep(0, n_knots), Sigma = sigma_gp^2 * K_knots)
  
  # Cross-covariance (obs x knots)
  K_cross <- matrix(0, N, n_knots)
  for (i in 1:N) {
    for (j in 1:n_knots) {
      d <- sqrt(sum((obs_coords[i,] - knot_coords[j,])^2))
      scaled <- sqrt3 * d / ls_std
      K_cross[i, j] <- (1 + scaled) * exp(-scaled)
    }
  }
  
  # Project to observations
  gp_at_obs <- K_cross %*% solve(K_knots, knot_vals)
  return(as.numeric(gp_at_obs))
}

# ─── Helper: fit model and save results ───
fit_and_save <- function(stan_data_sim, scenario_name, true_params, model, config) {
  outdir <- file.path("model_output", paste0("simrecovery_", scenario_name))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Save true parameters
  saveRDS(true_params, file.path(outdir, "true_params.rds"))
  
  # Save the synthetic stan data
  saveRDS(stan_data_sim, file.path(outdir, "stan_data_synthetic.rds"))
  
  cat("Fitting scenario", scenario_name, "...\n")
  fit_start <- Sys.time()
  
  fit <- model$sample(
    data = stan_data_sim,
    seed = 314,
    chains = config$chains,
    parallel_chains = config$chains,
    iter_warmup = as.integer(config$iter * 0.5),
    iter_sampling = as.integer(config$iter * 0.5),
    adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth,
    refresh = 200
  )
  
  fit_time <- as.numeric(difftime(Sys.time(), fit_start, units = "mins"))
  cat("Fitting completed in", round(fit_time, 1), "minutes\n")
  
  # Extract draws
  draws <- fit$draws(format = "df")
  key_params <- c("beta_0", "beta_oipc", "sigma", "lambda_decay",
                   "sigma_intercept_spatial", "sigma_slope_spatial",
                   "ls_intercept_km", "var_spatial_intercept", "var_spatial_slope",
                   "min_oipc_slope", "max_oipc_slope", "sd_oipc_slope",
                   "rmse", "r_squared",
                   "beta_c4", "beta_tree", "beta_shrub", "beta_grass",
                   "beta_oipc_x_c4", "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
  available <- intersect(key_params, names(draws))
  saveRDS(draws[, c(available, ".chain", ".iteration", ".draw")],
          file.path(outdir, "posterior_draws.rds"))
  
  # Save diagnostics
  saveRDS(fit$diagnostic_summary(), file.path(outdir, "diagnostics.rds"))
  
  # Save LOO
  saveRDS(fit$loo(), file.path(outdir, "loo.rds"))
  
  # Save runtime
  saveRDS(list(
    scenario = scenario_name,
    elapsed_mins = fit_time,
    start_time = fit_start,
    end_time = Sys.time(),
    git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                          error = function(e) NA_character_,
                          warning = function(w) NA_character_)
  ), file.path(outdir, "runtime_info.rds"))
  
  # Report recovery
  cat("\n  --- Recovery check for", scenario_name, "---\n")
  for (p in names(true_params)) {
    if (p %in% names(draws)) {
      post_mean <- mean(draws[[p]])
      post_q <- quantile(draws[[p]], c(0.1, 0.9))
      recovered <- true_params[[p]] >= post_q[1] && true_params[[p]] <= post_q[2]
      cat(sprintf("  %-25s true=%.3f  posterior=%.3f [%.3f, %.3f]  %s\n",
                  p, true_params[[p]], post_mean, post_q[1], post_q[2],
                  ifelse(recovered, "RECOVERED", "MISSED")))
    }
  }
  cat("\n")
  
  return(draws)
}

# ─── Use realistic parameter values from existing posteriors ───
# Take posterior means as "realistic" values for simulation.
# orig_draws is a posterior::draws_array — extract by name.
lambda_km <- mean(posterior::extract_variable(orig_draws, "lambda_decay"))
ls_km <- mean(posterior::extract_variable(orig_draws, "ls_intercept_km"))
coord_scale_km <- mean(stan_data$coord_scaling) * 111.0
ls_std <- ls_km / coord_scale_km  # Convert to standardized units

sigma_int_orig <- mean(posterior::extract_variable(orig_draws,
                                                   "sigma_intercept_spatial"))
sigma_int_std <- sigma_int_orig / stan_data$d2H_wax_sd_original  # Standardized

# Compute weighted OIPC for all observations
oipc_weighted <- compute_oipc_weighted(stan_data, lambda_km)

# Residual sigma (standardized)
sigma_resid <- 0.3  # Moderate residual noise


# ═══════════════════════════════════════════════════════════════
# SCENARIO 3a: Uniform slope, no spatial slope variation
# ═══════════════════════════════════════════════════════════════
if ("3a" %in% scenarios) {
  cat("━━━ SCENARIO 3a: Uniform slope = 0.7, no spatial slope variation ━━━\n")
  
  true_beta_0 <- 0.0
  true_beta_oipc <- 0.7
  
  # Generate intercept GP surface (realistic spatial structure)
  gp_intercept <- generate_gp_surface(stan_data, ls_std, sigma_int_std, seed = 100)
  
  # Generate d2H_wax: intercept + slope*OIPC + noise
  mu_sim <- (true_beta_0 + gp_intercept) + true_beta_oipc * oipc_weighted
  noise <- rnorm(stan_data$N, 0, sqrt(sigma_resid^2 + stan_data$d2H_wax_err^2))
  d2h_sim <- mu_sim + noise
  
  # Replace d2H_wax in stan_data
  sd_3a <- stan_data
  sd_3a$d2H_wax <- as.numeric(d2h_sim)
  
  true_params_3a <- list(beta_oipc = 0.7, beta_0 = 0.0, sigma = sigma_resid)
  
  fit_and_save(sd_3a, "3a", true_params_3a, model, config)
}


# ═══════════════════════════════════════════════════════════════
# SCENARIO 3b: Spatially varying slope
# ═══════════════════════════════════════════════════════════════
if ("3b" %in% scenarios) {
  cat("━━━ SCENARIO 3b: Spatially varying slope (tropics=0.5, high-lat=0.9) ━━━\n")
  
  true_beta_0 <- 0.0
  true_beta_oipc_global <- 0.7
  
  # Create latitude-dependent slope
  # Use raw latitude (unstandardized) for intuitive interpretation
  lat_raw <- stan_data$latitude
  # Linear gradient: slope = 0.5 at equator, 0.9 at poles
  # slope = 0.5 + 0.4 * (|lat| / 90)
  true_slope_spatial <- 0.5 + 0.4 * (abs(lat_raw) / 90)
  cat("  Slope range: [", round(min(true_slope_spatial), 3), ",", 
      round(max(true_slope_spatial), 3), "]\n")
  cat("  Slope mean:", round(mean(true_slope_spatial), 3), "\n")
  
  # Generate intercept GP (same as 3a for comparability)
  gp_intercept <- generate_gp_surface(stan_data, ls_std, sigma_int_std, seed = 100)
  
  # Generate d2H_wax with spatially varying slope
  mu_sim <- (true_beta_0 + gp_intercept) + true_slope_spatial * oipc_weighted
  noise <- rnorm(stan_data$N, 0, sqrt(sigma_resid^2 + stan_data$d2H_wax_err^2))
  d2h_sim <- mu_sim + noise
  
  sd_3b <- stan_data
  sd_3b$d2H_wax <- as.numeric(d2h_sim)
  
  # The global mean slope should be ~0.7
  true_params_3b <- list(beta_oipc = mean(true_slope_spatial), beta_0 = 0.0, sigma = sigma_resid)
  
  fit_and_save(sd_3b, "3b", true_params_3b, model, config)
}


# ═══════════════════════════════════════════════════════════════
# SCENARIO 3c: Intercept confounding stress test
# ═══════════════════════════════════════════════════════════════
if ("3c" %in% scenarios) {
  cat("━━━ SCENARIO 3c: Confounding stress test (intercept correlated with OIPC) ━━━\n")
  
  true_beta_0 <- 0.0
  true_beta_oipc <- 0.7
  
  # Create an intercept surface that is strongly correlated with OIPC
  # This mimics the worst case: the intercept GP could "steal" from the slope
  # Use mean OIPC across scales as the spatial pattern
  oipc_spatial_pattern <- rowMeans(stan_data$oipc_values)
  
  # Scale to have realistic intercept variance (~sigma_int_std)
  # intercept = alpha * OIPC_pattern + small GP noise
  alpha <- 0.3  # ~25% of posterior sigma_intercept_std, gives effective OLS slope ~1.0  # Strong correlation
  gp_noise <- generate_gp_surface(stan_data, ls_std, sigma_int_std * 0.3, seed = 200)
  confounding_intercept <- alpha * oipc_spatial_pattern + gp_noise
  
  cat("  Correlation between intercept and OIPC:", 
      round(cor(confounding_intercept, oipc_weighted), 3), "\n")
  
  # Generate d2H_wax
  mu_sim <- (true_beta_0 + confounding_intercept) + true_beta_oipc * oipc_weighted
  noise <- rnorm(stan_data$N, 0, sqrt(sigma_resid^2 + stan_data$d2H_wax_err^2))
  d2h_sim <- mu_sim + noise
  
  sd_3c <- stan_data
  sd_3c$d2H_wax <- as.numeric(d2h_sim)
  
  true_params_3c <- list(beta_oipc = 0.7, beta_0 = 0.0, sigma = sigma_resid)
  
  fit_and_save(sd_3c, "3c", true_params_3c, model, config)
}


cat("\n=== ALL SCENARIOS COMPLETE ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
