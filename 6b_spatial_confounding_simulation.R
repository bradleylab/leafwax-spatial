# run_confounding_test_v2.R
# Phase 3c: Spatial Confounding Simulation Study (CORRECTED)
#
# Fixes from v1:
#   1. Re-standardizes simulated d2H_wax (mean=0, sd=1) so model priors are
#      correctly calibrated. v1 fed unstandardized data (sd up to 2.0) into a
#      model with priors tuned for sd=1, producing inflated slope estimates.
#   2. Uses the SAME z_indep draw across all scenarios so only rho varies.
#      v1 used different seeds per scenario, confounding the comparison.
#   3. Rescales d2H_wax_err to match new standardization.
#   4. Computes true beta_oipc in re-standardized space for correct evaluation.
#
# Design follows Paciorek (2010, Biostatistics) and Dupont et al. (2022).
# See validation_log/NOTE_spatial_confounding_theory.md for full background.
#
# Usage: Rscript run_confounding_test_v2.R [scenario]
#   scenario: "rho00", "rho03", "rho05", "empirical", or "all" (default)

library(cmdstanr)
library(posterior)
library(loo)
library(MASS)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[1] != "all") {
  scenarios <- args[1]
} else {
  scenarios <- c("rho00", "rho03", "rho05", "empirical")
}

cat("=== SPATIAL CONFOUNDING SIMULATION v2 (Paciorek 2010 framework) ===\n")
cat("Scenarios:", paste(scenarios, collapse = ", "), "\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")

# Run from repo root

SRC_MODEL <- "baseline_sp"
draws_path     <- file.path("model_output", SRC_MODEL, "posterior_draws.rds")
stan_data_path <- file.path("prepared_data",
                            paste0("stan_data_", SRC_MODEL, ".rds"))
config_path    <- file.path("prepared_data",
                            paste0("config_", SRC_MODEL, ".rds"))

for (p in c(draws_path, stan_data_path, config_path)) {
  if (!file.exists(p)) {
    stop("Required input missing: ", p,
         "\n  Source model = '", SRC_MODEL, "'. Fit it first.")
  }
}

cat("Source model for confounding simulation:", SRC_MODEL, "\n")
stan_data <- readRDS(stan_data_path)
config <- readRDS(config_path)
orig_draws <- readRDS(draws_path)

cat("Compiling Stan model...\n")
model <- cmdstan_model("4d_leaf_wax_spatial_model.stan")
cat("Model compiled.\n\n")

# ─── Parameters from fitted model ───
# orig_draws is a posterior::draws_array — extract by name.
lambda_km <- mean(posterior::extract_variable(orig_draws, "lambda_decay"))
ls_km <- mean(posterior::extract_variable(orig_draws, "ls_intercept_km"))
coord_scale_km <- mean(stan_data$coord_scaling) * 111.0
ls_std <- ls_km / coord_scale_km

# sigma_intercept_spatial is back-transformed to original units in generated
# quantities, so dividing by d2H_wax_sd_original gives standardized units
sigma_int_std <- mean(posterior::extract_variable(orig_draws,
                                                  "sigma_intercept_spatial")) /
  stan_data$d2H_wax_sd_original
sigma_resid <- 0.3

cat("GP intercept amplitude (standardized):", round(sigma_int_std, 4), "\n")
cat("GP length scale (km):", round(ls_km, 1), "\n")
cat("Lambda decay (km):", round(lambda_km, 1), "\n\n")

# Compute weighted OIPC using posterior mean lambda
scales <- stan_data$distance_scales
weights <- exp(-scales / lambda_km)
weights <- weights / sum(weights)
oipc_weighted <- as.numeric(stan_data$oipc_values %*% weights)

TRUE_BETA_OIPC_UNSTD <- 0.7
TRUE_BETA_0 <- 0.0
N <- stan_data$N

# ─── Build Matern 3/2 covariance at observation locations ───
cat("Computing Matern 3/2 covariance at observation locations...\n")
obs_coords <- stan_data$coords
sqrt3 <- sqrt(3)

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
K_obs <- K_obs + diag(1e-4, N)
cat("Covariance matrix computed.\n\n")

# ─── Generate z_indep ONCE, reuse across scenarios ───
cat("Drawing independent GP (shared across all scenarios)...\n")
set.seed(42)
z_indep <- mvrnorm(1, mu = rep(0, N), Sigma = sigma_int_std^2 * K_obs)
cat("  z_indep: mean=", round(mean(z_indep), 4),
    " sd=", round(sd(z_indep), 4),
    " cor(z_indep, oipc)=", round(cor(z_indep, oipc_weighted), 4), "\n\n")

# ─── Confounding intercept generator ───
generate_confounding_intercept <- function(oipc, z_indep, rho, sigma_z) {
  oipc_std <- (oipc - mean(oipc)) / sd(oipc)
  z_confound <- rho * sigma_z * oipc_std + sqrt(1 - rho^2) * z_indep
  actual_rho <- cor(z_confound, oipc)
  cat(sprintf("  Target rho=%.2f, achieved rho=%.3f\n", rho, actual_rho))
  return(z_confound)
}

# ─── Fit helper: re-standardize synthetic d2H_wax ───
fit_and_save <- function(stan_data_orig, d2h_sim, scenario_name,
                         true_beta_unstd, metadata, model, config) {
  outdir <- file.path("model_output", paste0("confounding_v2_", scenario_name))
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Re-standardize simulated d2H_wax to mean=0, sd=1
  sim_mean <- mean(d2h_sim)
  sim_sd <- sd(d2h_sim)
  d2h_std <- (d2h_sim - sim_mean) / sim_sd
  cat(sprintf("  Raw d2h_sim: mean=%.4f sd=%.4f\n", sim_mean, sim_sd))
  cat(sprintf("  After re-standardization: mean=%.4f sd=%.4f\n",
              mean(d2h_std), sd(d2h_std)))

  # True causal slope in re-standardized space
  true_beta_std <- true_beta_unstd / sim_sd
  cat(sprintf("  True beta_oipc (unstd space): %.4f\n", true_beta_unstd))
  cat(sprintf("  True beta_oipc (re-std space): %.4f\n", true_beta_std))

  # Build modified stan_data
  sd_sim <- stan_data_orig
  sd_sim$d2H_wax <- as.numeric(d2h_std)

  # Rescale measurement errors to new standardization
  # Original d2H_wax_err is in units of (original permil / d2H_wax_sd_original)
  # New standardization divides by sim_sd additionally
  sd_sim$d2H_wax_err <- stan_data_orig$d2H_wax_err / sim_sd

  # Update back-transformation parameters for generated quantities
  sd_sim$d2H_wax_sd_original <- sim_sd * stan_data_orig$d2H_wax_sd_original
  sd_sim$d2H_wax_mean_original <- sim_mean * stan_data_orig$d2H_wax_sd_original +
    stan_data_orig$d2H_wax_mean_original

  # Save true params and metadata
  true_params <- list(
    beta_oipc_unstd = true_beta_unstd,
    beta_oipc_std = true_beta_std,
    beta_0 = TRUE_BETA_0,
    sigma_unstd = sigma_resid,
    sigma_std = sigma_resid / sim_sd,
    rho = metadata$rho,
    sigma_z = sigma_int_std,
    sim_mean = sim_mean,
    sim_sd = sim_sd
  )

  saveRDS(true_params, file.path(outdir, "true_params.rds"))
  saveRDS(metadata, file.path(outdir, "experiment_metadata.rds"))
  saveRDS(sd_sim, file.path(outdir, "stan_data_synthetic.rds"))

  cat("  Fitting model...\n")
  fit_start <- Sys.time()

  # Cap max_treedepth at 12 for simulation runs. Higher values (14) cause
  # infinite-length trajectories when confounding creates funnel geometry.
  # 12 = 4096 leapfrog steps max, vs 14 = 16384. Treedepth warnings are
  # expected and documented; they don't invalidate the recovery test.
  sim_treedepth <- min(config$max_treedepth, 12)
  cat(sprintf("  max_treedepth: %d (config: %d)\n",
              sim_treedepth, config$max_treedepth))

  fit <- model$sample(
    data = sd_sim,
    seed = 314,
    chains = config$chains,
    parallel_chains = config$chains,
    iter_warmup = as.integer(config$iter * 0.5),
    iter_sampling = as.integer(config$iter * 0.5),
    adapt_delta = config$adapt_delta,
    max_treedepth = sim_treedepth,
    refresh = 100
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
    git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE),
                          error = function(e) NA_character_,
                          warning = function(w) NA_character_),
    version = "v2_restandardized"
  ), file.path(outdir, "runtime_info.rds"))

  # Recovery check (in re-standardized space)
  cat("\n  --- Recovery check (re-standardized space) ---\n")
  beta_post <- draws$beta_oipc
  bias <- mean(beta_post) - true_beta_std
  ci80 <- quantile(beta_post, c(0.1, 0.9))
  ci95 <- quantile(beta_post, c(0.025, 0.975))
  coverage_80 <- true_beta_std >= ci80[1] & true_beta_std <= ci80[2]
  coverage_95 <- true_beta_std >= ci95[1] & true_beta_std <= ci95[2]

  cat(sprintf("  beta_oipc: true=%.4f  posterior=%.4f [%.4f, %.4f]\n",
    true_beta_std, mean(beta_post), ci95[1], ci95[2]))
  cat(sprintf("  Bias: %.4f\n", bias))
  cat(sprintf("  80%% CI coverage: %s\n", ifelse(coverage_80, "YES", "NO")))
  cat(sprintf("  95%% CI coverage: %s\n", ifelse(coverage_95, "YES", "NO")))

  # OLS on re-standardized data
  ols <- lm(d2h_std ~ oipc_weighted)
  ols_slope <- coef(ols)[2]
  cat(sprintf("  OLS slope (re-std): %.4f  (true causal: %.4f)\n",
    ols_slope, true_beta_std))
  cat(sprintf("  OLS bias from confounding: %.4f\n",
    ols_slope - true_beta_std))
  cat(sprintf("  Spatial model bias: %.4f\n", bias))
  cat(sprintf("  Confounding absorbed: %.1f%%\n",
    100 * (1 - bias / (ols_slope - true_beta_std))))
  cat("\n")

  return(list(
    true_beta_std = true_beta_std,
    bias = bias,
    coverage_80 = coverage_80,
    coverage_95 = coverage_95,
    post_mean = mean(beta_post),
    post_sd = sd(beta_post),
    ols_slope = ols_slope,
    sim_sd = sim_sd
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

  cat(sprintf("\n━━━ Scenario 3c_%s: rho = %.2f ━━━\n", scenario, rho))

  # Generate confounding intercept (same z_indep for all scenarios)
  z_confound <- generate_confounding_intercept(
    oipc = oipc_weighted,
    z_indep = z_indep,
    rho = rho,
    sigma_z = sigma_int_std
  )

  # Generate synthetic d2H_wax (in unstandardized simulation space)
  mu_sim <- (TRUE_BETA_0 + z_confound) + TRUE_BETA_OIPC_UNSTD * oipc_weighted
  noise <- rnorm(N, 0, sqrt(sigma_resid^2 + stan_data$d2H_wax_err^2))
  d2h_sim <- mu_sim + noise

  metadata <- list(
    scenario = scenario,
    rho = rho,
    sigma_z = sigma_int_std,
    sigma_resid = sigma_resid,
    true_beta_oipc_unstd = TRUE_BETA_OIPC_UNSTD,
    ls_std = ls_std,
    lambda_km = lambda_km,
    version = "v2",
    fixes = c(
      "Re-standardized simulated d2H_wax to mean=0 sd=1",
      "Same z_indep across all scenarios (seed=42)",
      "Rescaled d2H_wax_err for new standardization"
    ),
    design = "Paciorek (2010) correlated GP framework"
  )

  results[[scenario]] <- fit_and_save(
    stan_data, d2h_sim, paste0("3c_", scenario),
    TRUE_BETA_OIPC_UNSTD, metadata, model, config
  )
}

# ─── Summary table ───
cat("\n=== CONFOUNDING SIMULATION v2 SUMMARY ===\n")
cat("True beta_oipc = 0.700 (unstandardized); varies per scenario after",
    "re-standardization\n\n")
cat(sprintf("%-12s %5s %9s %9s %9s %7s %6s %6s %9s\n",
  "Scenario", "rho", "True.Std", "Post.Mean", "Post.SD",
  "Bias", "Cov80", "Cov95", "OLS.Slope"))
for (sc in names(results)) {
  r <- results[[sc]]
  cat(sprintf("%-12s %5.2f %9.4f %9.4f %9.4f %7.4f %6s %6s %9.4f\n",
    sc, rho_values[sc], r$true_beta_std, r$post_mean, r$post_sd,
    r$bias, ifelse(r$coverage_80, "YES", "NO"),
    ifelse(r$coverage_95, "YES", "NO"), r$ols_slope))
}

saveRDS(results, "model_output/confounding_v2_summary.rds")

cat("\n=== ALL CONFOUNDING v2 SCENARIOS COMPLETE ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
