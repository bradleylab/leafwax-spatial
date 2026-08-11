#!/usr/bin/env Rscript
# check_stan_r_parity.R
#
# Numerical Stan <-> R prediction parity for the chordal GP. The
# acceptance test only checks R-internal self-consistency; this checks that the
# leafwax package's predict_spatial_dual_gp() (chordal branch) reproduces the
# spatial field the Stan model actually produced, on the SAME posterior draws.
#
# Stan saves, per calibration site, alpha_spatial = beta_0 + intercept_GP and
# beta_oipc_spatial = beta_oipc + slope_GP (transformed parameters). The package
# predicts the GP deviations at those same sites from the knot effects; we compare
#   package$intercept  vs  (alpha_spatial   - beta_0)
#   package$slope      vs  (beta_oipc_spatial - beta_oipc).
#
# Reads saved RDS + the package source only (no fit, no cmdstanr). Run after the
# pilot fit on the pulled-back output:
#   Rscript scripts/check_stan_r_parity.R \
#     model_output/baseline_sp prepared_data/stan_data_baseline_sp.rds [n_draws]

suppressWarnings(suppressMessages({
  library(posterior)
  source("../leafwax-pkg/R/spatial_interpolation.R")  # predict_spatial_dual_gp (+ helpers)
}))

args <- commandArgs(trailingOnly = TRUE)
pilot_dir <- if (length(args) >= 1) args[[1]] else "model_output/baseline_sp"
stan_data_f <- if (length(args) >= 2) args[[2]] else "prepared_data/stan_data_baseline_sp.rds"
n_draws <- if (length(args) >= 3) as.integer(args[[3]]) else 50L
TOL <- 1e-3   # standardized-response units; jitter/solve differences are ~1e-6

fails <- 0
gate <- function(name, ok, detail = "") {
  if (!isTRUE(ok)) fails <<- fails + 1
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" — ", detail) else ""))
}

sd <- readRDS(stan_data_f)
if (isTRUE(sd$include_gp != 1)) stop("stan_data is not a spatial (GP) model.")
draws <- posterior::as_draws_df(readRDS(file.path(pilot_dir, "posterior_draws.rds")))
if (nrow(draws) > n_draws) draws <- draws[round(seq.int(1, nrow(draws), length.out = n_draws)), ]
nd <- nrow(draws)

# Site + knot coordinates in DEGREES (what the package consumes)
coords_deg <- cbind(sd$longitude, sd$latitude)
knot_deg   <- sd$knot_coords_deg
d2H_sd <- if (!is.null(sd$scaling_params$d2H_sd)) sd$scaling_params$d2H_sd else sd$d2H_wax_sd_original
scaling <- list(d2H_sd = d2H_sd)

# The chordal deposit must be tagged; the raw draws here are chordal-fit, so call
# the chordal branch explicitly (mirrors load_posteriors stamping spatial_metric).
pkg <- predict_spatial_dual_gp(coords_deg, knot_deg, draws, scaling, metric = "chordal")

# Stan-side per-site field (n_draws x N), minus the global terms.
alpha_cols <- grep("^alpha_spatial\\[", names(draws), value = TRUE)
slope_cols <- grep("^beta_oipc_spatial\\[", names(draws), value = TRUE)
gate("Stan alpha_spatial / beta_oipc_spatial present",
     length(alpha_cols) == nrow(coords_deg) && length(slope_cols) == nrow(coords_deg),
     sprintf("alpha=%d slope=%d sites=%d", length(alpha_cols), length(slope_cols), nrow(coords_deg)))

stan_alpha <- as.matrix(draws[, alpha_cols])       # nd x N
stan_slope <- as.matrix(draws[, slope_cols])
beta0 <- draws[["beta_0"]]
betaoipc <- draws[["beta_oipc"]]
stan_int_dev   <- stan_alpha - matrix(beta0,    nd, ncol(stan_alpha))
stan_slope_dev <- stan_slope - matrix(betaoipc, nd, ncol(stan_slope))

int_diff   <- max(abs(pkg$intercept - stan_int_dev))
slope_diff <- max(abs(pkg$slope     - stan_slope_dev))
# relative to the field's own scale
int_scale   <- max(abs(stan_int_dev)); slope_scale <- max(abs(stan_slope_dev))
cat(sprintf("  intercept field: max|pkg - stan| = %.2e (field max |dev| = %.3f)\n", int_diff, int_scale))
cat(sprintf("  slope field:     max|pkg - stan| = %.2e (field max |dev| = %.3f)\n", slope_diff, slope_scale))

gate("intercept GP parity within tol", int_diff < TOL, sprintf("%.2e", int_diff))
gate("slope GP parity within tol", slope_diff < TOL, sprintf("%.2e", slope_diff))

cat("\n")
if (fails == 0) {
  cat("STAN<->R PARITY PASSED — package chordal prediction reproduces the Stan field.\n")
} else {
  cat(sprintf("STAN<->R PARITY: %d check(s) failed — investigate before trusting package predictions.\n", fails))
  quit(status = 1)
}
