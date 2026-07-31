#!/usr/bin/env Rscript
# validate_chordal_prep.R
#
# Local (no-fit) end-to-end validation of the chordal prep path on the REAL frozen
# sediment: builds stan_data in memory for representative models and asserts the
# chordal coords / km regularization / length-scale-prior structure. This catches
# integration bugs in prepare_stan_data() that the unit-level acceptance test
# (which only exercises the helpers) cannot. Prep is data preparation, not a fit,
# so it runs locally; the MCMC fits still go to C2.
#
# Runs in memory only — does NOT write to prepared_data/ (that dir has symlinks to
# a retired location; a real run preps into a fresh dir per README_chordal_run.md).
#
# Run from the leafwax_working root:  Rscript scripts/validate_chordal_prep.R

suppressWarnings(suppressMessages({
  source("0_load_config.R")
  source(CONFIG$scripts$spatial_functions)   # 4a: prepare_stan_data, validate_stan_data, ...
  library(tidyverse)
}))

R_EARTH <- 6371
fails <- 0
check <- function(name, cond, detail = "") {
  ok <- isTRUE(cond)
  if (!ok) fails <<- fails + 1
  cat(sprintf("  [%s] %s%s\n", if (ok) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" — ", detail) else ""))
}

# --- Load real sediment + replicate 4b's scaling + jitter setup ----------------
sediment <- readRDS(CONFIG$input_data)
cat("Loaded", nrow(sediment), "sediment records from", CONFIG$input_data, "\n")

set.seed(12345)
jit <- CONFIG$coordinate_jitter
sediment <- sediment %>%
  mutate(longitude = longitude + rnorm(n(), 0, jit),
         latitude  = latitude  + rnorm(n(), 0, jit))

SCALING_PARAMS <- list(
  d2H_mean = mean(sediment$d2H_wax, na.rm = TRUE), d2H_sd = sd(sediment$d2H_wax, na.rm = TRUE),
  oipc_mean = mean(sediment$oipc_d2h20, na.rm = TRUE), oipc_sd = sd(sediment$oipc_d2h20, na.rm = TRUE),
  elev_mean = mean(sediment$elevation_gmted, na.rm = TRUE), elev_sd = sd(sediment$elevation_gmted, na.rm = TRUE),
  c4_mean = CONFIG$c4_standardization$mean, c4_sd = CONFIG$c4_standardization$sd
)
if (isTRUE(CONFIG$climate_standardization$compute_from_data) && "annual_precip" %in% names(sediment)) {
  SCALING_PARAMS$precip_mean <- mean(sediment$annual_precip, na.rm = TRUE)
  SCALING_PARAMS$precip_sd   <- sd(sediment$annual_precip, na.rm = TRUE)
}
has_pft <- all(c("pft_tree", "pft_shrub", "pft_grass") %in% names(sediment))

prep_one <- function(model_name) {
  cfg <- CONFIG$model_configs[[model_name]]
  suppressWarnings(suppressMessages(
    prepare_stan_data(
      data = sediment,
      include_c4 = cfg$include_c4,
      include_pft = ifelse(is.null(cfg$include_pft), FALSE, cfg$include_pft) && has_pft,
      include_gp = cfg$include_gp,
      include_elevation = cfg$include_elevation,
      include_precip = ifelse(is.null(cfg$include_precip), FALSE, cfg$include_precip),
      include_temp = FALSE, include_vpd = FALSE, include_soil = FALSE,
      n_pp_knots = cfg$n_pp_knots,
      SCALING_PARAMS = SCALING_PARAMS, has_pft_columns = has_pft,
      apply_range_factor = cfg$apply_range_factor
    )
  ))
}

assert_spatial <- function(sd, model, expect_rf) {
  cat("== ", model, " ==\n")
  # validate_stan_data must pass (it now checks chordal dims, norms, tau > 0, ls bounds)
  vd <- tryCatch({ validate_stan_data(sd); TRUE }, error = function(e) { cat("   validate ERROR:", conditionMessage(e), "\n"); FALSE })
  check("validate_stan_data passes", vd)
  check("coords are N x 3", all(dim(sd$coords) == c(sd$N, 3)))
  check("knot_coords are n_pp_knots x 3", all(dim(sd$knot_coords) == c(sd$n_pp_knots, 3)))
  cn <- sqrt(rowSums(sd$coords^2)); kn <- sqrt(rowSums(sd$knot_coords^2))
  check("site chordal norms ~ 6371 km", max(abs(cn - R_EARTH)) < 1,
        sprintf("max|norm-R| = %.2e", max(abs(cn - R_EARTH))))
  check("knot chordal norms ~ 6371 km", max(abs(kn - R_EARTH)) < 1)
  check("tau_spatial_slope strictly > 0", all(sd$tau_spatial_slope > 0),
        sprintf("min = %.3f", min(sd$tau_spatial_slope)))
  check("tau_spatial_intercept strictly > 0", all(sd$tau_spatial_intercept > 0))
  check("oipc_range_at_knots >= 0 and finite", all(sd$oipc_range_at_knots >= 0) && all(is.finite(sd$oipc_range_at_knots)))
  check("ls prior fields present + ordered",
        is.finite(sd$ls_log_lower) && is.finite(sd$ls_log_upper) && sd$ls_log_lower < sd$ls_log_upper &&
        is.finite(sd$ls_prior_mean_log) && sd$ls_prior_sd_log > 0)
  check("ls bounds match config (871 / 6434 km)",
        abs(exp(sd$ls_log_lower) - 870.8) < 1 && abs(exp(sd$ls_log_upper) - 6434.1) < 1)
  # apply_range_factor behavior: when OFF, slope tau == intercept tau (base) everywhere
  if (!expect_rf) {
    check("range_factor OFF -> tau_slope == tau_intercept",
          max(abs(sd$tau_spatial_slope - sd$tau_spatial_intercept)) < 1e-12)
    check("apply_range_factor flag == 0", sd$apply_range_factor == 0L)
  } else {
    check("range_factor ON -> some slope tau < intercept tau (shrinkage active)",
          any(sd$tau_spatial_slope < sd$tau_spatial_intercept - 1e-9))
    check("apply_range_factor flag == 1", sd$apply_range_factor == 1L)
  }
}

# baseline_sp (rf ON), full_interact_sp (rf ON, all covariates), baseline_sp_rfoff (rf OFF)
assert_spatial(prep_one("baseline_sp"),        "baseline_sp",        expect_rf = TRUE)
assert_spatial(prep_one("full_interact_sp"),   "full_interact_sp",   expect_rf = TRUE)
assert_spatial(prep_one("baseline_sp_rfoff"),  "baseline_sp_rfoff",  expect_rf = FALSE)

# non-spatial sanity: builds without a GP, coords still N x 3, no crash
cat("== baseline (non-spatial) ==\n")
sdn <- prep_one("baseline")
check("non-spatial: builds + coords N x 3", all(dim(sdn$coords) == c(sdn$N, 3)))
check("non-spatial: n_pp_knots == 1 placeholder", sdn$n_pp_knots == 1)

cat("\n")
if (fails == 0) cat("ALL CHORDAL PREP VALIDATION CHECKS PASSED\n") else {
  cat(sprintf("%d CHECK(S) FAILED\n", fails)); quit(status = 1)
}
