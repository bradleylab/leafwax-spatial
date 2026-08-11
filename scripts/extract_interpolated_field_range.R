# extract_interpolated_field_range.R
#
# Computes the global interpolated range of the spatial intercept field
# alpha_spatial(s) and the spatial slope field beta_oipc_spatial(s) for
# each spatial model using the reported map interpolation method (Matérn
# 3/2 kernel with the full conditional GP posterior mean), vectorized and
# applied to all 9 spatial models, with results saved to CSV. These are the
# values used for the reported intercept and slope field ranges
# ("interpolated intercept field ranges from X to Y ‰ globally" and
# "interpolated slope field ranges from X to Y").
#
# Method (chordal metric, matches Figure_04_spatial_maps.R):
#   1. Build a 2°×2° global grid (lon=-179..179, lat=-89..89).
#   2. Convert the grid to 3-D chordal coordinates (km) with
#      lonlat_to_chordal(); knots are stan_data$knot_coords, already 3-D
#      chordal km. Euclidean distance over these coords IS chordal km, so the
#      length scale (ls_intercept_km / ls_slope_km) enters directly in km with
#      no standardized-coordinate conversion.
#   3. Build K_knot_knot via Matérn 3/2 on chordal distance, add the model
#      jitter (1e-4), and solve the GP conditional mean.
#   4. Per draw d, predict the field at every grid cell:
#        pred_z(s)   = k_gk(s)' (K_kk + 1e-4 I)^{-1} z_knots
#        alpha_d(s)  = pred_z_int(s) × sigma_intercept_spatial[d]  (per mil)
#        slope_d(s)  = beta_oipc[d] + pred_z_slp(s) × sigma_slope_spatial[d]
#      The posterior-mean field is E[alpha(s)] = mean_d alpha_d(s) (and
#      likewise for the slope), accumulated as a running sum over a
#      deterministic subsample of draws. This replaces the former
#      plug-in-of-posterior-means approximation (interpolating the mean
#      z/sigma/rho, a plug-in surface rather than the true posterior-mean
#      field).
#   5. Output: per-model min/max/SD of the posterior-mean field, plus the
#      across-model envelope for prose.
#
# Reproducibility: deterministic given SEED (fixed draw subsample).

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(fields)
})

source("scripts/posterior_helpers.R")
.load_saved_model_config <- load_config
source("4a_spatial_functions.R")  # lonlat_to_chordal()

# Posterior-mean field via a deterministic draw subsample (see Method above).
N_DRAWS <- 1000
SEED <- 42

GRID_LON_BY <- 2
GRID_LAT_BY <- 2

# Set LEAFWAX_OUTPUT_DIR to override the generated-output directory.
.output_dir <- Sys.getenv(
  "LEAFWAX_OUTPUT_DIR",
  unset = "model_analysis/spatial_pattern_diagnostics"
)
dir.create(.output_dir, recursive = TRUE, showWarnings = FALSE)
OUT_CSV <- file.path(.output_dir, "interpolated_field_ranges.csv")

spatial_models <- c(
  "baseline_sp", "baseline_env_sp", "baseline_veg_sp",
  "full_sp", "full_interact_sp", "elevation_only_sp",
  "elevation_c4_sp", "c4_only_sp", "elevation_c4_interact_sp"
)

matern32 <- function(d, alpha = 1, rho = 1) {
  s <- sqrt(3) * d / rho
  alpha^2 * (1 + s) * exp(-s)
}

interpolate_one <- function(m) {
  draws  <- load_draws(m)
  sd_obj <- load_stan_data(m)
  cfg <- .load_saved_model_config(m)

  # Knots are 3-D chordal km (stan_data$knot_coords); use directly.
  knot_coords <- sd_obj$knot_coords
  n_knots <- nrow(knot_coords)

  # 2°×2° global grid → 3-D chordal km (same metric as the fitted model)
  lon_seq <- seq(-179, 179, by = GRID_LON_BY)
  lat_seq <- seq(-89,   89, by = GRID_LAT_BY)
  pred_grid <- expand.grid(lon = lon_seq, lat = lat_seq)
  pred_coords <- lonlat_to_chordal(pred_grid$lon, pred_grid$lat)  # K x 3
  n_grid <- nrow(pred_grid)

  # Chordal distances (km) are constant across draws; only rho varies.
  D_kk <- fields::rdist(knot_coords)               # 125 x 125
  D_gk <- fields::rdist(pred_coords, knot_coords)  # K x 125

  # Per-draw parameters (draw rows in a common order across variables)
  z_int_mat   <- as_draws_matrix(subset_draws(draws, variable = "z_intercept_spatial"))
  z_slp_mat   <- as_draws_matrix(subset_draws(draws, variable = "z_slope_spatial"))
  sigma_int_v <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "sigma_intercept_spatial")))
  sigma_slp_v <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "sigma_slope_spatial")))
  beta_oipc_v <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "beta_oipc")))
  ls_int_v    <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "ls_intercept_km")))
  ls_slp_v    <- as.numeric(as_draws_matrix(subset_draws(draws, variable = "ls_slope_km")))
  stopifnot(ncol(z_int_mat) == n_knots, ncol(z_slp_mat) == n_knots)

  # Deterministic subsample of draws (all draws if fewer than N_DRAWS)
  n_total <- nrow(z_int_mat)
  set.seed(SEED)
  draw_idx <- if (n_total > N_DRAWS) sort(sample.int(n_total, N_DRAWS)) else seq_len(n_total)

  # Running-sum accumulation of the per-draw fields → posterior-mean field
  alpha_sum <- numeric(n_grid)
  slope_sum <- numeric(n_grid)
  for (d in draw_idx) {
    rho_int <- ls_int_v[d]
    rho_slp <- ls_slp_v[d]

    K_kk_int <- matern32(D_kk, alpha = 1, rho = rho_int) + diag(1e-4, n_knots)
    K_kk_slp <- matern32(D_kk, alpha = 1, rho = rho_slp) + diag(1e-4, n_knots)
    K_gk_int <- matern32(D_gk, alpha = 1, rho = rho_int)
    K_gk_slp <- matern32(D_gk, alpha = 1, rho = rho_slp)

    # as.numeric(): a draws_matrix row stays a 1xK draws_matrix, which solve()
    # rejects (b must be K x n); coerce to a plain length-K numeric vector.
    pred_z_int <- as.vector(K_gk_int %*% solve(K_kk_int, as.numeric(z_int_mat[d, ])))
    pred_z_slp <- as.vector(K_gk_slp %*% solve(K_kk_slp, as.numeric(z_slp_mat[d, ])))

    alpha_sum <- alpha_sum + pred_z_int * sigma_int_v[d]
    slope_sum <- slope_sum + (beta_oipc_v[d] + pred_z_slp * sigma_slp_v[d])
  }
  n_used <- length(draw_idx)
  alpha_bar <- alpha_sum / n_used   # E[alpha(s)] per grid cell (per mil)
  slope_bar <- slope_model_to_physical(
    slope_sum / n_used,
    cfg$scaling_params
  ) # E[slope(s)] per grid cell, per mil wax per per mil precipitation

  tibble(
    model = m,
    n_grid = n_grid,
    n_knots = n_knots,
    n_draws = n_used,
    alpha_min_permil = min(alpha_bar),
    alpha_max_permil = max(alpha_bar),
    alpha_sd_permil  = sd(alpha_bar),
    slope_min = min(slope_bar),
    slope_max = max(slope_bar),
    slope_sd  = sd(slope_bar),
    ls_intercept_km = mean(ls_int_v[draw_idx]),
    sigma_intercept_permil = mean(sigma_int_v[draw_idx])
  )
}

cat("Interpolating", length(spatial_models),
    sprintf("spatial models on a %.0f°×%.0f° global grid…\n",
            GRID_LON_BY, GRID_LAT_BY))
out <- lapply(spatial_models, function(m) {
  cat("  ", m, "…\n")
  tryCatch(interpolate_one(m),
           error = function(e) {
             message("  failed: ", conditionMessage(e))
             NULL
           })
}) |> bind_rows()

write.csv(out, OUT_CSV, row.names = FALSE)

cat("\nWrote", OUT_CSV, "\n\n")
print(out)

cat(sprintf(
  "\n**Per-model interpolated alpha_spatial(s) ranges from %.0f to +%.0f ‰** (across %d models).\n",
  min(out$alpha_min_permil), max(out$alpha_max_permil), nrow(out)))
cat(sprintf(
  "**Per-model interpolated slope field ranges from %.2f to %.2f.**\n",
  min(out$slope_min), max(out$slope_max)))

# Also produce per-model alpha range (for the figure-caption-style sentence)
cat("\nPer-model alpha_spatial range:\n")
for (i in seq_len(nrow(out))) {
  cat(sprintf("  %-26s %.0f to +%.0f ‰  (sigma_intercept=%.1f)\n",
              out$model[i], out$alpha_min_permil[i], out$alpha_max_permil[i],
              out$sigma_intercept_permil[i]))
}
