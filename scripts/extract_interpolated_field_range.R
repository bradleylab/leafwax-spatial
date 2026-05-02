# extract_interpolated_field_range.R
#
# Computes the global interpolated range of the spatial intercept field
# alpha_spatial(s) and the spatial slope field beta_oipc_spatial(s) for
# each spatial model, using the SAME interpolation method as
# manuscript/figure_code/Figure_04_spatial_maps.R (Matérn 3/2 kernel +
# full conditional GP posterior mean), but vectorized and applied to all
# 9 spatial models, with results saved to CSV.
#
# This factors out a calculation that previously existed only inside
# Figure_04's render path. The numbers below are the same numbers
# Figure_04 prints to stdout when it builds Panel A (intercept) and
# Panel B (slope), and are the values quoted in the manuscript prose
# ("interpolated intercept field ranges from X to Y ‰ globally";
# "interpolated slope field ranges from X to Y").
#
# Method (matches Figure_04 lines 100–210):
#   1. Build a 2°×2° global grid (lon=-179..179, lat=-89..89).
#   2. Standardize grid coords with the same lon/lat mean and sd that
#      stan_prep used (recovered from sd_obj$longitude / sd_obj$latitude).
#   3. Convert ls_intercept_km / ls_slope_km from km to standardized
#      coordinate units using the formula in Figure_04 lines 117–120:
#        ls_lon_std = ls_km / (km_per_deg_lon × lon_sd)
#        ls_lat_std = ls_km / (km_per_deg_lat × lat_sd)
#        ls_std    = sqrt(ls_lon_std² + ls_lat_std²) / sqrt(2)
#      where km_per_deg_lon = 111 × cos(lat_mean); km_per_deg_lat = 111.
#   4. Build K_knot_knot via Matérn 3/2 in standardized space, invert.
#   5. For each grid point, build k_vec and predict
#        pred_z = k_vec' × K_inv × z_intercept_knots
#        alpha(s) = pred_z × sigma_intercept_spatial   (per mil)
#        beta_oipc(s) = beta_oipc + (k_vec' × K_inv × z_slope_knots)
#                       × sigma_slope_spatial          (unitless slope)
#   6. Output: per-model min/max/SD of the interpolated grid values, plus
#      the across-model envelope for prose.
#
# Reproducibility: deterministic — uses posterior mean summaries.

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(fields)
})

source("scripts/posterior_helpers.R")

OUT_CSV <- "model_analysis/spatial_pattern_diagnostics/interpolated_field_ranges.csv"
GRID_LON_BY <- 2
GRID_LAT_BY <- 2

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
  summ <- load_summaries(m)
  sd_obj <- load_stan_data(m)

  knot_coords_std <- sd_obj$knot_coords  # already standardized in stan_prep
  n_knots <- nrow(knot_coords_std)

  # Recover the standardization params that stan_prep used: mean+sd of
  # the observation lon/lat (these were what stan_prep used to
  # standardize knot_coords).
  lon_orig <- sd_obj$longitude
  lat_orig <- sd_obj$latitude
  lon_mean <- mean(lon_orig); lon_sd <- sd(lon_orig)
  lat_mean <- mean(lat_orig); lat_sd <- sd(lat_orig)

  z_intercept_knots <- summ$mean[startsWith(summ$variable, "z_intercept_spatial[")]
  z_slope_knots     <- summ$mean[startsWith(summ$variable, "z_slope_spatial[")]
  sigma_intercept   <- summ$mean[summ$variable == "sigma_intercept_spatial"]
  sigma_slope       <- summ$mean[summ$variable == "sigma_slope_spatial"]
  beta_oipc         <- summ$mean[summ$variable == "beta_oipc"]
  ls_intercept_km   <- summ$mean[summ$variable == "ls_intercept_km"]
  ls_slope_km       <- summ$mean[summ$variable == "ls_slope_km"]

  # Build 2°×2° grid (matches Figure_04)
  lon_seq <- seq(-179, 179, by = GRID_LON_BY)
  lat_seq <- seq(-89,   89, by = GRID_LAT_BY)
  pred_grid <- expand.grid(lon = lon_seq, lat = lat_seq)
  pred_lon_std <- (pred_grid$lon - lon_mean) / lon_sd
  pred_lat_std <- (pred_grid$lat - lat_mean) / lat_sd
  pred_coords_std <- cbind(pred_lon_std, pred_lat_std)

  # Length-scale conversion (matches Figure_04 lines 117-120 + 168-170)
  km_per_deg_lon <- 111 * cos(lat_mean * pi / 180)
  km_per_deg_lat <- 111
  ls_to_std <- function(ls_km) {
    sqrt((ls_km / (km_per_deg_lon * lon_sd))^2 +
         (ls_km / (km_per_deg_lat * lat_sd))^2) / sqrt(2)
  }
  rho_int <- ls_to_std(ls_intercept_km)
  rho_slp <- ls_to_std(ls_slope_km)

  # Knot×knot kernels (vectorized)
  D_kk <- fields::rdist(knot_coords_std)
  K_kk_int <- matern32(D_kk, alpha = 1, rho = rho_int) + diag(1e-9, n_knots)
  K_kk_slp <- matern32(D_kk, alpha = 1, rho = rho_slp) + diag(1e-9, n_knots)
  K_kk_int_inv <- solve(K_kk_int)
  K_kk_slp_inv <- solve(K_kk_slp)

  # Grid×knot kernels
  D_gk <- fields::rdist(pred_coords_std, knot_coords_std)
  K_gk_int <- matern32(D_gk, alpha = 1, rho = rho_int)
  K_gk_slp <- matern32(D_gk, alpha = 1, rho = rho_slp)

  # Predicted standardized z's at each grid cell
  z_int_grid <- as.vector(K_gk_int %*% K_kk_int_inv %*% z_intercept_knots)
  z_slp_grid <- as.vector(K_gk_slp %*% K_kk_slp_inv %*% z_slope_knots)

  # Back-transform: alpha (per mil), slope (unitless)
  alpha_grid <- z_int_grid * sigma_intercept
  slope_grid <- beta_oipc + z_slp_grid * sigma_slope

  tibble(
    model = m,
    n_grid = nrow(pred_grid),
    n_knots = n_knots,
    alpha_min_permil = min(alpha_grid),
    alpha_max_permil = max(alpha_grid),
    alpha_sd_permil  = sd(alpha_grid),
    slope_min = min(slope_grid),
    slope_max = max(slope_grid),
    slope_sd  = sd(slope_grid),
    ls_intercept_km = ls_intercept_km,
    sigma_intercept_permil = sigma_intercept
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

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
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
