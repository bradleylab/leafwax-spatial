#!/usr/bin/env Rscript
# test_chordal_acceptance.R
#
# Fast, dependency-light acceptance tests for the chordal-metric refit
# (implementation spec v2, section 7 — the subset that does not require the full
# prep pipeline or a Stan fit). Validates:
#   1. chordal coordinate geometry (norm = R; known-pair distances)
#   2. compute_spatial_tau replicates the former in-Stan regularization logic
#   3. mPP kriging self-consistency (prediction at knots recovers knot effects)
#   4. length-scale prior/bounds match the config faithful-translation targets
#
# The heavier checks (build real stan_data from the frozen sediment + a small
# Stan smoke fit, Stan-vs-R prediction parity) run in the pre-launch gate with a
# fresh prepared-data directory, immediately before the HPC batch.
#
# Run:  Rscript scripts/test_chordal_acceptance.R   (from the leafwax_working root)

suppressWarnings(suppressMessages({
  source("4a_spatial_functions.R")                 # lonlat_to_chordal(lon, lat, R), compute_spatial_tau, ...
  # 4a's lonlat_to_chordal takes two vectors; the package's takes one lon/lat
  # matrix. They live in separate namespaces in production; in this test both
  # are sourced into the global env, so capture the 4a (two-vector) version
  # before the package definition masks it.
  chordal_from_lonlat <- lonlat_to_chordal
  source("../leafwax-pkg/R/spatial_interpolation.R")  # predict_one_gp_mpp, matern32, pair_distances, lonlat_to_chordal(matrix)
  library(yaml)
}))

R_EARTH <- 6371
fails <- 0
check <- function(name, cond, detail = "") {
  status <- if (isTRUE(cond)) "PASS" else "FAIL"
  if (!isTRUE(cond)) fails <<- fails + 1
  cat(sprintf("  [%s] %s%s\n", status, name, if (nzchar(detail)) paste0(" — ", detail) else ""))
}
gc_km <- function(lon1, lat1, lon2, lat2) {  # great-circle km, for comparison
  geosphere::distHaversine(c(lon1, lat1), c(lon2, lat2), r = R_EARTH * 1000) / 1000
}

cat("== 1. Chordal geometry ==\n")
xyz <- chordal_from_lonlat(c(0, 0, 90, -30, 137), c(0, 1, 0, 45, -25))
norms <- sqrt(rowSums(xyz^2))
check("3-D norm == R for all points", max(abs(norms - R_EARTH)) < 1e-6,
      sprintf("max|norm-R| = %.2e", max(abs(norms - R_EARTH))))

d_1deg <- as.numeric(dist(chordal_from_lonlat(c(0, 0), c(0, 1))))
check("1 deg latitude ~ 111.19 km chordal", abs(d_1deg - 111.19) < 0.05,
      sprintf("chordal = %.3f km, gc = %.3f km", d_1deg, gc_km(0, 0, 0, 1)))

d_90 <- as.numeric(dist(chordal_from_lonlat(c(0, 90), c(0, 0))))
check("90 deg equator chordal = 2R sin(45) = 9010.5 km",
      abs(d_90 - 2 * R_EARTH * sin(pi / 4)) < 1e-3,
      sprintf("chordal = %.1f km, gc = %.1f km (chordal < gc as expected)", d_90, gc_km(0, 0, 90, 0)))

# continental scale (~30 deg): chordal should track great circle to ~1.5%
d_30 <- as.numeric(dist(chordal_from_lonlat(c(0, 30), c(0, 0))))
gc_30 <- gc_km(0, 0, 30, 0)
check("30 deg deviation from great circle < 1.5%", abs(d_30 - gc_30) / gc_30 < 0.015,
      sprintf("chordal = %.1f km, gc = %.1f km, dev = %.2f%%", d_30, gc_30, 100 * (gc_30 - d_30) / gc_30))

cat("== 2. compute_spatial_tau replicates former in-Stan logic ==\n")
dens <- c(0, 5, 15, 20)
rng  <- c(0.10, 0.50, 1.00, 0.00)   # max = 1.0
# base tau: d==0 -> .50 ; d<10 -> .50+.30*d/10 ; else .80
exp_base   <- c(0.50, 0.65, 0.80, 0.80)
# range_factor vs max=1.0: <0.25 -> .2 ; <0.60 -> .5 ; else 1.0
exp_rf     <- c(0.2, 0.5, 1.0, 0.2)
tau_on  <- compute_spatial_tau(dens, rng, apply_range_factor = TRUE)
tau_off <- compute_spatial_tau(dens, rng, apply_range_factor = FALSE)
check("base_tau matches piecewise density rule", max(abs(tau_on$base_tau - exp_base)) < 1e-12)
check("range_factor matches 0.25/0.60 thresholds", max(abs(tau_on$range_factor - exp_rf)) < 1e-12)
check("tau_intercept == base_tau (range_factor slope-only)",
      max(abs(tau_on$tau_intercept - exp_base)) < 1e-12)
check("tau_slope == base_tau * range_factor when enabled",
      max(abs(tau_on$tau_slope - exp_base * exp_rf)) < 1e-12)
check("tau_slope == base_tau when range_factor disabled",
      max(abs(tau_off$tau_slope - exp_base)) < 1e-12)
# max_range == 0 edge: all range_factor default to 1.0 (matches Stan)
tau_zero <- compute_spatial_tau(c(1, 2), c(0, 0), apply_range_factor = TRUE)
check("all-zero OIPC range -> range_factor = 1.0 (Stan edge case)",
      all(tau_zero$range_factor == 1.0))

cat("== 2b. Exact-boundary behaviour (thresholds and radius operators) ==\n")
# base_tau at the density boundary: d < 10 is strict, so d == 10 -> 0.80
bt <- compute_spatial_tau(c(9, 10, 11), c(1, 1, 1), apply_range_factor = FALSE)$base_tau
check("base_tau at d = 9/10/11 = 0.77/0.80/0.80",
      max(abs(bt - c(0.50 + 0.30 * 0.9, 0.80, 0.80))) < 1e-12,
      sprintf("got %.3f/%.3f/%.3f", bt[1], bt[2], bt[3]))
# range_factor thresholds are strict <: r == 0.25*max -> 0.5 (not 0.2);
# r == 0.60*max -> 1.0 (not 0.5). Append a knot at 1.0 so max = 1.0 and the
# boundaries land exactly at 0.25 and 0.60; test the first four.
rf <- compute_spatial_tau(c(1, 1, 1, 1, 1), c(0.2499, 0.25, 0.5999, 0.60, 1.0),
                          apply_range_factor = TRUE)$range_factor
check("range_factor at r = 0.2499/0.25/0.5999/0.60 (max=1) = 0.2/0.5/0.5/1.0",
      max(abs(rf[1:4] - c(0.2, 0.5, 0.5, 1.0))) < 1e-12,
      sprintf("got %.1f/%.1f/%.1f/%.1f", rf[1], rf[2], rf[3], rf[4]))
# density count uses <= (inclusive): an obs exactly at the radius is counted.
knot0 <- chordal_from_lonlat(0, 0)
obs1  <- chordal_from_lonlat(0, 1)                 # 111.194 km from knot
d1    <- as.numeric(dist(rbind(knot0, obs1)))
check("density count <= radius: obs exactly at radius is counted",
      calculate_knot_data_density_km(knot0, obs1, radius_km = d1, verbose = FALSE) == 1 &&
      calculate_knot_data_density_km(knot0, obs1, radius_km = d1 * 0.999, verbose = FALSE) == 0)
# OIPC-range count uses strict < (and requires > 5 inside): 7 obs exactly at the
# radius are all excluded (range 0); nudging the radius outward includes them.
obs7 <- chordal_from_lonlat(rep(0, 7), rep(1, 7))
r_at  <- compute_oipc_range_at_knots(knot0, obs7, 1:7, radius_km = d1)
r_out <- compute_oipc_range_at_knots(knot0, obs7, 1:7, radius_km = d1 * 1.001)
check("OIPC-range uses strict <: obs exactly at radius excluded (range 0)",
      r_at == 0 && abs(r_out - 6) < 1e-9,
      sprintf("at-radius = %.3f, just-outside = %.3f", r_at, r_out))

cat("== 3. mPP kriging self-consistency (predict at knots recovers knot effects) ==\n")
cat("   NOTE: this checks R-internal consistency only; it would pass for any\n")
cat("   internally consistent metric. The definitive Stan-vs-package prediction\n")
cat("   parity check runs at the pre-launch gate on a real fit.\n")
set.seed(314)
n_knots <- 20
knot_lonlat <- cbind(runif(n_knots, -160, 160), runif(n_knots, -50, 60))
z <- matrix(rnorm(n_knots), nrow = 1)          # 1 draw
sigma <- 1.0
ls_km <- 2000
pred_at_knots <- predict_one_gp_mpp(knot_lonlat, knot_lonlat, z, sigma, ls_km,
                                    metric = "chordal")
target <- sigma * z[1, ]
check("prediction at knot locations ~ sigma*z (jitter-limited)",
      max(abs(pred_at_knots[1, ] - target)) < 1e-2,
      sprintf("max abs diff = %.2e", max(abs(pred_at_knots[1, ] - target))))
# length scale used directly in km: doubling ls must change off-knot predictions
off <- cbind(0, 0)
p1 <- predict_one_gp_mpp(off, knot_lonlat, z, sigma, 1000, metric = "chordal")
p2 <- predict_one_gp_mpp(off, knot_lonlat, z, sigma, 4000, metric = "chordal")
check("length scale (km) affects off-knot prediction", abs(p1[1, 1] - p2[1, 1]) > 1e-6)

cat("== 4. Length-scale prior / bounds vs config (nominal km re-expression) ==\n")
cfg <- yaml::read_yaml("config.yaml")$gp_length_scale
check("exp(prior_mean_log_km) ~ 2367 km", abs(exp(cfg$prior_mean_log_km) - 2367) < 2,
      sprintf("= %.1f km", exp(cfg$prior_mean_log_km)))
check("lower_km ~ 870.8", abs(cfg$lower_km - 870.8) < 0.5)
check("upper_km ~ 6434.1", abs(cfg$upper_km - 6434.1) < 0.5)
check("prior_sd_log == 0.4", abs(cfg$prior_sd_log - 0.4) < 1e-9)

cat("== 5. Metric dispatch (predict a posterior with the metric it was fit under) ==\n")
set.seed(7)
K <- 6; nd <- 4
knots5 <- cbind(runif(K, -120, 120), runif(K, -40, 50))
zi <- as.data.frame(matrix(rnorm(nd * K), nd, K)); names(zi) <- paste0("z_intercept_spatial[", 1:K, "]")
zs <- as.data.frame(matrix(rnorm(nd * K), nd, K)); names(zs) <- paste0("z_slope_spatial[", 1:K, "]")
draws5 <- cbind(zi, zs,
                sigma_intercept_spatial = rep(2, nd),  # permil-scale; /d2H_sd inside
                sigma_slope_spatial = rep(0.1, nd),
                ls_intercept_km = rep(2000, nd))
scaling5 <- list(lon_mean = 0, lon_sd = 80, lat_mean = 10, lat_sd = 28, d2H_sd = 20)
newpt <- matrix(c(-90, 38), nrow = 1)
# (a) omitting metric must error (fail-safe, not a silent default)
err <- tryCatch({ predict_spatial_dual_gp(newpt, knots5, draws5, scaling5); "no_error" },
                error = function(e) "errored")
check("dual_gp omitting metric errors (fail-safe)", identical(err, "errored"))
err_one <- tryCatch({ predict_one_gp_mpp(newpt, knots5, matrix(rnorm(K), 1, K), 1, 2000); "no_error" },
                    error = function(e) "errored")
check("one_gp_mpp omitting metric errors (fail-safe)", identical(err_one, "errored"))
# (b) the two metrics give different predictions (the silent mismatch codex measured)
ch <- predict_spatial_dual_gp(newpt, knots5, draws5, scaling5, metric = "chordal")
st <- predict_spatial_dual_gp(newpt, knots5, draws5, scaling5, metric = "standardized")
check("chordal vs standardized predictions differ",
      abs(mean(ch$slope) - mean(st$slope)) > 1e-6,
      sprintf("|Delta slope mean| = %.4f", abs(mean(ch$slope) - mean(st$slope))))
# (c) standardized depends on lon/lat scaling; chordal does not — confirms the
# dispatch really routes to different geometries.
scaling5b <- modifyList(scaling5, list(lon_sd = 40))
st2 <- predict_spatial_dual_gp(newpt, knots5, draws5, scaling5b, metric = "standardized")
ch2 <- predict_spatial_dual_gp(newpt, knots5, draws5, scaling5b, metric = "chordal")
check("standardized depends on lon/lat scaling; chordal does not",
      abs(mean(st$slope) - mean(st2$slope)) > 1e-6 &&
      abs(mean(ch$slope) - mean(ch2$slope)) < 1e-12)

cat("\n")
if (fails == 0) {
  cat("ALL CHORDAL ACCEPTANCE CHECKS PASSED\n")
} else {
  cat(sprintf("%d CHECK(S) FAILED\n", fails))
  quit(status = 1)
}
