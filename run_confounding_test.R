# run_confounding_test.R
# Phase 3c: Spatial Confounding Simulation Study
#
# Tests whether the Bayesian spatial model can separate a spatially-varying
# intercept from the OIPC slope when both share spatial structure.
#
# Design follows Paciorek (2010, Biostatistics) and Dupont et al. (2022):
#   - Generate confounding intercept as a GP with specified correlation (rho)
#     to the OIPC spatial pattern
#   - True slope is fixed at beta_oipc = 0.7
#   - Vary rho across scenarios to test different confounding strengths
#
# Scenarios:
#   3c_rho00: rho = 0.0 (no confounding; intercept independent of OIPC)
#   3c_rho03: rho = 0.3 (weak confounding)
#   3c_rho05: rho = 0.5 (moderate confounding; standard in literature)
#   3c_empirical: rho = 0.45 (matched to empirical residual-OIPC correlation)
#
# References:
#   Paciorek (2010) Biostatistics 11(4):601-15. doi:10.1093/biostatistics/kxq024
#   Dupont, Wood & Augustin (2022) Biometrics 78(4):1279-90. doi:10.1111/biom.13656
#   Hodges & Reich (2010) Am Stat 64(4):325-34. doi:10.1198/tast.2010.10052
#
# Usage: Rscript run_confounding_test.R [scenario]
#   scenario: "rho00", "rho03", "rho05", "empirical", or "all" (default)
#
# Expected runtime: ~6 hrs per scenario on r6i.8xlarge (solo)

library(cmdstanr)
library(posterior)
library(loo)
library(MASS)  # for mvrnorm

set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] != "all") {
  scenarios <- args[1]
} else {
  scenarios <- c("rho00", "rho03", "rho05", "empirical")
}

cat("=== SPATIAL CONFOUNDING SIMULATION (Paciorek 2010 framework) ===\n")
cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

setwd("/home/shared/leafwax_spatial/spatial_leafwax_model")

# Load prepared data and compile model
stan_data <- readRDS("prepared_data_consolidated/stan_data_baseline_veg_sp.rds")
config <- readRDS("prepared_data_consolidated/config_baseline_veg_sp.rds")
orig_draws <- readRDS("model_output/baseline_veg_sp_fixed_range/posterior_draws.rds")

cat("Compiling Stan model...\n")
model <- cmdstan_model("4d_leaf_wax_spatial_model.stan")
cat("Model compiled.\n\n")

# ─── Realistic parameter values from fitted model ───
lambda_km <- mean(orig_draws$lambda_decay)
ls_km <- mean(orig_draws$ls_intercept_km)
coord_scale_km <- mean(stan_data$coord_scaling) * 111.0
ls_std <- ls_km / coord_scale_km

sigma_int_std <- mean(orig_draws$sigma_intercept_spatial) / stan_data$d2H_wax_sd_original
sigma_resid <- 0.3

# Compute weighted OIPC
scales <- stan_data$distance_scales
weights <- exp(-scales / lambda_km)
weights <- weights / sum(weights)
oipc_weighted <- as.numeric(stan_data$oipc_values %*% weights)

# ─── True parameters (same for all scenarios) ───
TRUE_BETA_OIPC <- 0.7
TRUE_BETA_0 <- 0.0

# ─── Generate Matern 3/2 covariance matrix at observation locations ───
# For generating the confounding intercept as a proper GP
cat("Computing Matern 3/2 covariance at observation locations...\n")
N <- stan_data$N
obs_coords <- stan_data$coords
sqrt3 <- sqrt(3)

# This is expensive (818 x 818) but only done once
K_obs <- matrix(0, N, N)
for (i in 1:N) {
  K_obs[i, i] <- 1.0
  if (i < N) {
    for (j in (i + 1):N) {
      d <- sqrt(sum((obs_coords[i, ] - obs_coords[j, ])^2))
      scaled <- sqrt3 * d / ls_std
      K_obs[i, j] <- (1 + scaled) * exp(-scaled)
      K_obs[j, i] <- K_obs[i, j]
    }
  }
}
K_obs <- K_obs + diag(1e-4, N)  # jitter for numerical stability
cat("Covariance matrix computed.\n\n")

# ─── Generate correlated GP surface (Paciorek 2010 method) ───
# Given a target correlation rho between the GP surface and OIPC,
# construct: Z = rho * (sigma_z/sigma_x) * X_proj + sqrt(1-rho^2) * Z_indep
# where X_proj is OIPC projected onto the GP basis and Z_indep is an
# independent GP draw. This ensures cor(Z, X) = rho while Z retains
# proper GP spatial structure.
generate_confounding_intercept <- function(oipc, K_obs, rho, sigma_z, seed) {
  set.seed(seed)
  N <- length(oipc)

  # Standardize OIPC to unit variance for clean correlation control
  oipc_std <- (oipc - mean(oipc)) / sd(oipc)

  # Draw an independent GP surface
  z_indep <- mvrnorm(1, mu = rep(0, N), Sigma = sigma_z^2 * K_obs)

  # Construct Z with target correlation rho to OIPC
  # Z = rho * sigma_z * oipc_std + sqrt(1 - rho^2) * z_indep
  # This gives cor(Z, oipc) ≈ rho (exact for large N)
  z_confound <- rho * sigma_z * oipc_std + sqrt(1 - rho^2) * z_indep

  # Verify actual correlation
  actual_rho <- cor(z_confound, oipc)
  cat(sprintf("  Target rho=%.2f, achieved rho=%.3f\n", rho, actual_rho))

  return(z_confound)
}

# ─── Fit and save helper ───
fit_and_save <- function(stan_data_sim, scenario_name, true_params, metadata,
                         model, config) {
  outdir <- file.path("model_output", paste0("confounding_", scenario_name))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Save everything needed for reproducibility
  saveRDS(true_params, file.path(outdir, "true_params.rds"))
  saveRDS(metadata, file.path(outdir, "experiment_metadata.rds"))
  saveRDS(stan_data_sim, file.path(outdir, "stan_data_synthetic.rds"))

  cat("  Fitting model...\n")
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
  cat("  Completed in", round(fit_time, 1), "minutes\n")

  # Extract draws
  draws <- fit$draws(format = "df")
  key_params <- c(
    "beta_0", "beta_oipc", "sigma", "lambda_decay",
    "sigma_intercept_spatial", "sigma_slope_spatial",
    "ls_intercept_km", "var_spatial_intercept", "var_spatial_slope",
    "min_oipc_slope", "max_oipc_slope", "sd_oipc_slope",
    "rmse", "r_squared"
  )
  available <- intersect(key_params, names(draws))
  saveRDS(draws[, c(available, ".chain", ".iteration", ".draw")],
    file.path(outdir, "posterior_draws.rds"))
  saveRDS(fit$diagnostic_summary(), file.path(outdir, "diagnostics.rds"))
  saveRDS(fit$loo(), file.path(outdir, "loo.rds"))
  saveRDS(list(
    scenario = scenario_name,
    elapsed_mins = fit_time,
    start_time = fit_start,
    end_time = Sys.time(),
    git_commit = system("git rev-parse HEAD", intern = TRUE)
  ), file.path(outdir, "runtime_info.rds"))

  # Recovery check
  cat("\n  --- Recovery check ---\n")
  beta_post <- draws$beta_oipc
  bias <- mean(beta_post) - true_params$beta_oipc
  coverage_80 <- true_params$beta_oipc >= quantile(beta_post, 0.1) &
    true_params$beta_oipc <= quantile(beta_post, 0.9)
  coverage_95 <- true_params$beta_oipc >= quantile(beta_post, 0.025) &
    true_params$beta_oipc <= quantile(beta_post, 0.975)

  cat(sprintf("  beta_oipc: true=%.3f  posterior=%.3f [%.3f, %.3f]\n",
    true_params$beta_oipc, mean(beta_post),
    quantile(beta_post, 0.025), quantile(beta_post, 0.975)))
  cat(sprintf("  Bias: %.4f\n", bias))
  cat(sprintf("  80%% CI coverage: %s\n", ifelse(coverage_80, "YES", "NO")))
  cat(sprintf("  95%% CI coverage: %s\n", ifelse(coverage_95, "YES", "NO")))

  # What would OLS give on this data?
  ols <- lm(stan_data_sim$d2H_wax ~ oipc_weighted)
  cat(sprintf("  OLS slope (for comparison): %.3f\n", coef(ols)[2]))
  cat("\n")

  return(list(
    bias = bias,
    coverage_80 = coverage_80,
    coverage_95 = coverage_95,
    post_mean = mean(beta_post),
    post_sd = sd(beta_post),
    ols_slope = coef(ols)[2]
  ))
}


# ─── Run scenarios ───
rho_values <- c(rho00 = 0.0, rho03 = 0.3, rho05 = 0.5, empirical = 0.45)
results <- list()

for (scenario in scenarios) {
  rho <- rho_values[scenario]
  if (is.na(rho)) {
    cat("Unknown scenario:", scenario, "- skipping\n")
    next
  }

  cat(sprintf("━━━ Scenario 3c_%s: rho = %.2f ━━━\n", scenario, rho))

  # Generate confounding intercept
  z_confound <- generate_confounding_intercept(
    oipc = oipc_weighted,
    K_obs = K_obs,
    rho = rho,
    sigma_z = sigma_int_std,
    seed = 300 + round(rho * 100)  # Deterministic seed per rho
  )

  # Generate synthetic d2H_wax
  mu_sim <- (TRUE_BETA_0 + z_confound) + TRUE_BETA_OIPC * oipc_weighted
  noise <- rnorm(N, 0, sqrt(sigma_resid^2 + stan_data$d2H_wax_err^2))
  d2h_sim <- mu_sim + noise

  # Replace d2H_wax
  sd_sim <- stan_data
  sd_sim$d2H_wax <- as.numeric(d2h_sim)

  true_params <- list(
    beta_oipc = TRUE_BETA_OIPC,
    beta_0 = TRUE_BETA_0,
    sigma = sigma_resid,
    rho = rho,
    sigma_z = sigma_int_std
  )

  metadata <- list(
    scenario = scenario,
    rho = rho,
    sigma_z = sigma_int_std,
    sigma_resid = sigma_resid,
    true_beta_oipc = TRUE_BETA_OIPC,
    ls_std = ls_std,
    lambda_km = lambda_km,
    design = "Paciorek (2010) correlated GP framework",
    description = paste0(
      "Confounding intercept generated as a GP with Matern 3/2 kernel ",
      "(same kernel and length scale as the analysis model), with ",
      "correlation rho=", rho, " to the OIPC spatial pattern. ",
      "True beta_oipc=", TRUE_BETA_OIPC, ". ",
      "Following Paciorek (2010, Biostatistics 11(4):601-15)."
    ),
    references = c(
      "Paciorek (2010) Biostatistics 11(4):601-15. doi:10.1093/biostatistics/kxq024",
      "Dupont, Wood & Augustin (2022) Biometrics 78(4):1279-90. doi:10.1111/biom.13656",
      "Hodges & Reich (2010) Am Stat 64(4):325-34. doi:10.1198/tast.2010.10052"
    )
  )

  results[[scenario]] <- fit_and_save(
    sd_sim, paste0("3c_", scenario), true_params, metadata, model, config
  )
}

# ─── Summary table ───
cat("\n=== CONFOUNDING SIMULATION SUMMARY ===\n")
cat("True beta_oipc = 0.700 in all scenarios\n\n")
cat(sprintf("%-12s %6s %8s %8s %6s %6s %8s\n",
  "Scenario", "rho", "Post.Mean", "Post.SD", "Bias", "Cov80", "OLS.Slope"))
for (sc in names(results)) {
  r <- results[[sc]]
  cat(sprintf("%-12s %6.2f %8.3f %8.3f %6.4f %6s %8.3f\n",
    sc, rho_values[sc], r$post_mean, r$post_sd, r$bias,
    ifelse(r$coverage_80, "YES", "NO"), r$ols_slope))
}

# Save summary
saveRDS(results, "model_output/confounding_summary.rds")

cat("\n=== ALL CONFOUNDING SCENARIOS COMPLETE ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
