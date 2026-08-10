# regen_manuscript_numbers.R
#
# Recompute numeric quantities that appear in the manuscript or supplement
# from the run selected by LEAFWAX_RUN_DIR. For the reported chordal results,
# select results/c2_run_20260728_chordal/model_output.
#
# Run from repo root:
#   Rscript scripts/regen_manuscript_numbers.R
#
# Reproducibility: deterministic — no Bayesian inference, just summaries
# of saved posterior draws + classical OLS/Moran/variogram fits with
# fixed seeds where stochastic.

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(tibble)
  library(spdep)
  library(gstat)
  library(sp)
})

source("scripts/posterior_helpers.R")

# Set LEAFWAX_OUTPUT_DIR to override the generated-output directory. Only the
# output path changes; every calculation is unchanged.
.output_dir <- Sys.getenv("LEAFWAX_OUTPUT_DIR", unset = "model_analysis/reported_outputs")
dir.create(.output_dir, recursive = TRUE, showWarnings = FALSE)
OUT_PATH <- file.path(.output_dir, "REGENERATED_NUMBERS.md")
SEED <- 42
set.seed(SEED)

# Helpers ---------------------------------------------------------------------

fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "g")
fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")
pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f%%"), 100 * x)

bullet <- function(label, val) cat(sprintf("- %s: **%s**\n", label, val))
header <- function(level, txt) cat(sprintf("%s %s\n\n", strrep("#", level), txt))

# Sediment data ---------------------------------------------------------------

sed <- load_sediment()
n_total <- nrow(sed)
elev_range <- range(sed$elevation, na.rm = TRUE)
lat_range  <- range(sed$latitude,  na.rm = TRUE)
n_doi <- length(unique(na.omit(sed$DOI)))
n_compilations <- sed |>
  count(compilation, sort = TRUE) |>
  rename(source = compilation, n = n) |>
  mutate(source = if_else(
    is.na(source) | trimws(source) == "",
    "Not recorded",
    source
  ))
n_measured <- sum(sed$has_measured_elevation, na.rm = TRUE)
oipc_se_range <- range(sed$oipc_se20, na.rm = TRUE)
wax_err_q <- quantile(sed$d2H_wax_err, c(0.05, 0.5, 0.95), na.rm = TRUE)
wax_err_range <- range(sed$d2H_wax_err, na.rm = TRUE)

# OLS Figure 1 fit ------------------------------------------------------------

ols_dat <- sed |> filter(!is.na(d2H_wax), !is.na(oipc_d2h20))
ols <- lm(d2H_wax ~ oipc_d2h20, data = ols_dat)
ols_summ <- summary(ols)
beta0 <- coef(ols)[1]
beta_oipc <- coef(ols)[2]
beta0_se <- ols_summ$coefficients[1, 2]
beta_oipc_se <- ols_summ$coefficients[2, 2]
ols_r2 <- ols_summ$r.squared
ols_rmse <- sqrt(mean(residuals(ols)^2))
ols_sigma <- ols_summ$sigma
ols_n <- nrow(ols_dat)
xbar <- mean(ols_dat$oipc_d2h20)

# 95% PI half-width at xbar (single new prediction, not mean response)
n_ols <- nrow(ols_dat)
sxx <- sum((ols_dat$oipc_d2h20 - xbar)^2)
pi_se <- ols_sigma * sqrt(1 + 1/n_ols + (xbar - xbar)^2 / sxx)
pi_halfwidth <- qt(0.975, df = ols$df.residual) * pi_se

# Classical OLS inverse diagnostic at d2H_wax = -180 --------------------------
# This remains a labelled Fieller/delta-method comparison and is not the
# likelihood-based Bayesian inversion implemented by the leafwax package.

y_target <- -180
xhat <- (y_target - beta0) / beta_oipc

# McClelland approximation (delta method)
var_b0 <- vcov(ols)[1, 1]
var_b1 <- vcov(ols)[2, 2]
cov_b0_b1 <- vcov(ols)[1, 2]
se_xhat_mc <- sqrt(
  var_b0 / beta_oipc^2
  + (xhat^2 * var_b1) / beta_oipc^2
  + 2 * xhat * cov_b0_b1 / beta_oipc^2
  + (ols_sigma^2) / beta_oipc^2
)
mc_ci <- xhat + c(-1, 1) * qt(0.975, ols$df.residual) * se_xhat_mc

# Fieller exact CI
t_crit <- qt(0.975, ols$df.residual)
g <- (t_crit^2 * var_b1) / beta_oipc^2
A <- xhat
B_disc <- (1 - g) * (t_crit^2 / beta_oipc^2) *
  (ols_sigma^2 + var_b0 + 2 * xhat * cov_b0_b1 + xhat^2 * var_b1)
fieller_lo <- (A - sqrt(B_disc)) / (1 - g)
fieller_hi <- (A + sqrt(B_disc)) / (1 - g)

# Spatial coords for spdep ----------------------------------------------------

coords <- ols_dat |>
  select(longitude, latitude) |>
  as.matrix()

# Surface-geodesic nearest-neighbor distances via spdep::knearneigh
knn1 <- knearneigh(coords, k = 1, longlat = TRUE)
nn_dists_km <- nbdists(knn2nb(knn1), coords, longlat = TRUE) |>
  unlist()
nn_mean <- mean(nn_dists_km)
nn_median <- median(nn_dists_km)
pct_within_10  <- mean(nn_dists_km <= 10)
pct_within_100 <- mean(nn_dists_km <= 100)

# Moran's I on OLS residuals --------------------------------------------------

knn8 <- knearneigh(coords, k = 8, longlat = TRUE)
nb8 <- knn2nb(knn8)
lw8 <- nb2listw(nb8, style = "W")
mor_resid <- moran.test(residuals(ols), lw8, alternative = "greater")
mor_raw   <- moran.test(ols_dat$d2H_wax,  lw8, alternative = "greater")

# Multi-distance Moran's I on residuals
multi_d_km <- c(100, 300, 500, 1000, 2000, 3000, 5000)
moran_multi <- lapply(multi_d_km, function(d) {
  nb_d <- dnearneigh(coords, 0, d, longlat = TRUE)
  if (any(card(nb_d) == 0)) {
    # use zero.policy to avoid errors at small d
    lw_d <- nb2listw(nb_d, style = "W", zero.policy = TRUE)
    out <- moran.test(residuals(ols), lw_d, alternative = "greater",
                       zero.policy = TRUE)
  } else {
    lw_d <- nb2listw(nb_d, style = "W")
    out <- moran.test(residuals(ols), lw_d, alternative = "greater")
  }
  list(d = d, I = unname(out$estimate[1]), p = out$p.value)
}) |> bind_rows()

# Variogram on residuals (Matérn 3/2) ----------------------------------------

sp_pts <- ols_dat |>
  mutate(resid = residuals(ols)) |>
  select(longitude, latitude, resid) |>
  as.data.frame()
sp::coordinates(sp_pts) <- ~ longitude + latitude
sp::proj4string(sp_pts) <- sp::CRS("+proj=longlat +datum=WGS84")

vg_emp <- variogram(resid ~ 1, data = sp_pts, cutoff = 5000, width = 100)
vg_init <- vgm(psill = var(residuals(ols)) * 0.5, model = "Mat",
               range = 500, nugget = var(residuals(ols)) * 0.5,
               kappa = 1.5)
vg_fit <- tryCatch(
  fit.variogram(vg_emp, model = vg_init, fit.kappa = FALSE),
  error = function(e) {
    message("variogram fit failed: ", conditionMessage(e))
    NULL
  }
)
if (!is.null(vg_fit)) {
  nugget <- vg_fit$psill[1]
  psill  <- vg_fit$psill[2]
  vg_range <- vg_fit$range[2]
} else {
  nugget <- NA_real_; psill <- NA_real_; vg_range <- NA_real_
}

# Per-model diagnostics -------------------------------------------------------

main_models <- c(
  "baseline", "baseline_env", "baseline_env_sp", "baseline_sp",
  "baseline_veg", "baseline_veg_sp", "c4_only_sp",
  "elevation_c4_interact_sp", "elevation_c4_sp", "elevation_only_sp",
  "full", "full_interact", "full_interact_sp", "full_sp"
)

diag_tbl <- lapply(main_models, function(m) {
  d <- readRDS(file.path(MODEL_RUN_DIR, m, "diagnostics.rds"))
  tibble(
    model = m,
    max_rhat = d$max_rhat,
    min_ess_bulk = d$min_ess_bulk,
    n_divergent = sum(d$summary$num_divergent),
    n_max_treedepth = sum(d$summary$num_max_treedepth),
    min_ebfmi = min(d$summary$ebfmi)
  )
}) |> bind_rows()

# Spatial intercept observation-level stats ----------------------------------

# alpha_spatial[i] is on standardized scale; back-transform with stan_data
spatial_models <- main_models[grepl("_sp$", main_models)]

intercept_stats <- lapply(spatial_models, function(m) {
  pd <- load_draws(m)
  cfg <- load_config(m)
  vars <- variables(pd)
  alpha_vars <- vars[startsWith(vars, "alpha_spatial[")]
  if (length(alpha_vars) == 0) return(NULL)
  alpha_means <- summarise_draws(subset_draws(pd, variable = alpha_vars), mean)
  alpha_pm <- alpha_means$mean * cfg$scaling_params$d2H_sd  # back to per mil
  # site order matches sediment row order
  oipc_orig <- sed$oipc_d2h20
  r <- cor(alpha_pm, oipc_orig)
  tibble(
    model = m,
    r = r,
    r2 = r^2,
    sd_permil = sd(alpha_pm),
    min_permil = min(alpha_pm),
    max_permil = max(alpha_pm)
  )
}) |> bind_rows()

# Spatial slope field stats --------------------------------------------------

slope_stats <- lapply(spatial_models, function(m) {
  pd <- load_draws(m)
  cfg <- load_config(m)
  vars <- variables(pd)
  slope_vars <- vars[startsWith(vars, "beta_oipc_spatial[")]
  if (length(slope_vars) == 0) return(NULL)
  slope_means <- slope_model_to_physical(
    summarise_draws(subset_draws(pd, variable = slope_vars), mean)$mean,
    cfg$scaling_params
  )
  iqr_val <- IQR(slope_means)
  tibble(
    model = m,
    mean_slope = mean(slope_means),
    sd_slope = sd(slope_means),
    iqr_slope = iqr_val,
    min_slope = min(slope_means),
    max_slope = max(slope_means)
  )
}) |> bind_rows()

# 95% prediction interval widths (back-transformed) -------------------------

pi_models <- c("baseline", "baseline_env", "full", "full_interact",
               "baseline_sp", "baseline_env_sp", "full_sp", "full_interact_sp")

pi_stats <- lapply(pi_models, function(m) {
  pd <- load_draws(m)
  cfg <- load_config(m)
  vars <- variables(pd)
  rep_vars <- vars[startsWith(vars, "d2H_rep[")]
  if (length(rep_vars) == 0) return(NULL)
  # convert each draw to per mil, then per-site quantile
  draws_mat <- unclass(as_draws_matrix(subset_draws(pd, variable = rep_vars)))
  draws_permil <- draws_mat * cfg$scaling_params$d2H_sd + cfg$scaling_params$d2H_mean
  q025 <- apply(draws_permil, 2, quantile, 0.025)
  q975 <- apply(draws_permil, 2, quantile, 0.975)
  widths <- q975 - q025
  tibble(
    model = m,
    mean_width = mean(widths),
    sd_width = sd(widths),
    min_width = min(widths),
    max_width = max(widths)
  )
}) |> bind_rows()

# Data-density strata used in Supplementary Section S2.5.2 -------------------
# Assign each observation the density count of its nearest fitted chordal knot,
# then calculate in-sample RMSE from the posterior mean fitted value. This uses
# the exact Stan inputs and posterior draws from baseline_env_sp.
density_stan <- load_stan_data("baseline_env_sp")
density_draws <- load_draws("baseline_env_sp")
nearest_knot <- vapply(seq_len(density_stan$N), function(i) {
  delta <- sweep(density_stan$knot_coords, 2, density_stan$coords[i, ], "-")
  which.min(rowSums(delta^2))
}, integer(1))
site_density <- density_stan$knot_data_density[nearest_knot]
density_group <- cut(
  site_density,
  breaks = c(-Inf, 9, 50, Inf),
  labels = c("Sparse", "Medium", "Dense"),
  right = TRUE
)
mu_vars <- paste0("mu[", seq_len(density_stan$N), "]")
mu_mean <- colMeans(unclass(as_draws_matrix(
  subset_draws(density_draws, variable = mu_vars)
)))
residual_permil <- (density_stan$d2H_wax - mu_mean) *
  density_stan$d2H_wax_sd_original
density_stats <- tibble(
  category = density_group,
  residual_permil = residual_permil
) |>
  group_by(category) |>
  summarise(
    n = n(),
    rmse_permil = sqrt(mean(residual_permil^2)),
    .groups = "drop"
  )

# Residual sigma per spatial model -------------------------------------------

sigma_stats <- lapply(spatial_models, function(m) {
  pd <- load_draws(m)
  cfg <- load_config(m)
  sigma_draws <- as_draws_matrix(subset_draws(pd, variable = "sigma"))
  sigma_permil <- as.numeric(sigma_draws) * cfg$scaling_params$d2H_sd
  tibble(
    model = m,
    sigma_mean = mean(sigma_permil),
    sigma_lo = quantile(sigma_permil, 0.025),
    sigma_hi = quantile(sigma_permil, 0.975)
  )
}) |> bind_rows()

# Models used for the residual-variance and detection-threshold summaries
sigma_stats_focus <- sigma_stats |> filter(model %in% c(
  "baseline_sp", "baseline_env_sp", "full_sp", "full_interact_sp"))

# Detection thresholds -------------------------------------------------------

sigma_resid_med <- median(sigma_stats_focus$sigma_mean)
sigma_anal <- 3
sigma_total <- sqrt(sigma_resid_med^2 + sigma_anal^2)
det_thresh <- function(rho) {
  1.96 * sqrt(2 * sigma_resid_med^2 * (1 - rho) + 2 * sigma_anal^2)
}
det_table <- tibble(
  rho = c(0, 0.5, 0.8, 0.9),
  threshold_permil = sapply(rho, det_thresh)
)

# Detection-threshold endpoints: spatial at rho=0.9 vs non-spatial at rho=0
sigma_baseline_mean <- {
  pd <- load_draws("baseline")
  cfg <- load_config("baseline")
  sigma_draws <- as_draws_matrix(subset_draws(pd, variable = "sigma"))
  mean(as.numeric(sigma_draws) * cfg$scaling_params$d2H_sd)
}
spatial_threshold_rho09 <- 1.96 * sqrt(
  2 * sigma_resid_med^2 * 0.1 + 2 * sigma_anal^2
)
nonspatial_threshold_rho0 <- 1.96 * sqrt(
  2 * sigma_baseline_mean^2 + 2 * sigma_anal^2
)

# WRITE REPORT ----------------------------------------------------------------

sink(OUT_PATH)
header(1, "Regenerated numeric quantities")
cat(sprintf("Model run: `%s`.\n\n", RUN_ID))

header(2, "Dataset")
cat(sprintf("- **n = %s** sites\n", fmt_int(n_total)))
cat(sprintf("- Elevation range: **%s to %s m**\n",
            fmt_int(round(elev_range[1])), fmt_int(round(elev_range[2]))))
cat(sprintf("- Latitude range: %.1f° to %.1f°\n",
            lat_range[1], lat_range[2]))
cat(sprintf("- Unique DOIs: **%d**\n", n_doi))
cat(sprintf("- Sites with measured (not raster) elevation: **%d / %d** (%s)\n",
            n_measured, n_total, pct(n_measured / n_total)))
cat("\n### Compilations\n\n")
cat("| Source | n |\n|---|---:|\n")
for (i in seq_len(nrow(n_compilations))) {
  cat(sprintf("| %s | %d |\n", n_compilations$source[i], n_compilations$n[i]))
}
cat("\n")

header(2, "OLS Figure 1 (point-value fit)")
cat(sprintf("- β₀ = **%.2f ± %.2f**\n", beta0, beta0_se))
cat(sprintf("- β_OIPC = **%.3f ± %.3f**\n",
            beta_oipc, beta_oipc_se))
cat(sprintf("- R² = **%.3f**\n", ols_r2))
cat(sprintf("- RMSE = **%.1f ‰**\n", ols_rmse))
cat(sprintf("- σ̂ (residual SD) = **%.2f ‰**\n", ols_sigma))
cat(sprintf("- n = %d\n", ols_n))
cat(sprintf("- 95%% PI half-width at x̄: **%.1f ‰**\n", pi_halfwidth))
cat("\n### Inverse prediction at δ²H_wax = −180 ‰\n\n")
cat(sprintf("- Point estimate: **x̂ = %.1f ‰**\n", xhat))
cat(sprintf("- McClelland 95%% CI: [%.1f, %.1f]\n",
            mc_ci[1], mc_ci[2]))
if (B_disc > 0 && (1 - g) > 0) {
  cat(sprintf("- Fieller 95%% CI: [%.1f, %.1f]\n",
              fieller_lo, fieller_hi))
} else {
  cat("- Fieller 95% CI: discriminant negative or g >= 1 — undefined\n")
}
cat("\n")

header(2, "Nearest-neighbor distances")
cat(sprintf("- Mean NN distance: **%.1f km**\n", nn_mean))
cat(sprintf("- Median NN distance: **%.1f km**\n", nn_median))
cat(sprintf("- %% within 10 km: **%s**\n", pct(pct_within_10)))
cat(sprintf("- %% within 100 km: **%s**\n", pct(pct_within_100)))

header(2, "Moran's I (residuals of OLS)")
cat(sprintf("- k=8 NN: I = **%.3f**, p ≈ %s\n",
            unname(mor_resid$estimate[1]),
            ifelse(mor_resid$p.value < 1e-3, "0", fmt(mor_resid$p.value))))
cat(sprintf("- k=8 NN raw δ²H_wax: I = **%.3f**\n",
            unname(mor_raw$estimate[1])))
cat("\n| d (km) | I | p |\n|---:|---:|---:|\n")
for (i in seq_len(nrow(moran_multi))) {
  cat(sprintf("| %d | %.3f | %s |\n",
              moran_multi$d[i], moran_multi$I[i],
              ifelse(moran_multi$p[i] < 1e-3, "~0", fmt(moran_multi$p[i]))))
}
cat("\n")

header(2, "Variogram (residuals, Matérn 3/2)")
cat(sprintf("- Nugget: **%.0f ‰²**\n", nugget))
cat(sprintf("- Partial sill: **%.0f ‰²**\n", psill))
cat(sprintf("- Range: **%.0f km**\n", vg_range))
cat("\n*The empirical variogram does not reach a clean sill, so the parametric range is sensitive to functional form. The fitted GP length scale is the model-based spatial-range summary.*\n\n")

header(2, "Per-model diagnostics (S2.4.6)")
cat("| Model | max R̂ | min ESS bulk | divergences | max-treedepth | min EBFMI |\n")
cat("|---|---:|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(diag_tbl))) {
  cat(sprintf("| %s | %.4f | %d | %d | %d | %.3f |\n",
              diag_tbl$model[i],
              diag_tbl$max_rhat[i],
              as.integer(diag_tbl$min_ess_bulk[i]),
              as.integer(diag_tbl$n_divergent[i]),
              as.integer(diag_tbl$n_max_treedepth[i]),
              diag_tbl$min_ebfmi[i]))
}
cat(sprintf("\n**Min ESS across all models = %d** (%s).\n",
            as.integer(min(diag_tbl$min_ess_bulk)),
            diag_tbl$model[which.min(diag_tbl$min_ess_bulk)]))
cat(sprintf("**Min ESS for spatial models = %d** (%s).\n\n",
            as.integer(min(diag_tbl$min_ess_bulk[grepl("_sp$", diag_tbl$model)])),
            diag_tbl$model[grepl("_sp$", diag_tbl$model)][
              which.min(diag_tbl$min_ess_bulk[grepl("_sp$", diag_tbl$model)])]))

header(2, "Spatial intercept — observation-level posterior mean")
cat("| Model | r (vs OIPC) | r² | SD (‰) | Range (‰) |\n|---|---:|---:|---:|---|\n")
for (i in seq_len(nrow(intercept_stats))) {
  cat(sprintf("| %s | %.3f | %.3f | %.1f | %.1f to %.1f |\n",
              intercept_stats$model[i],
              intercept_stats$r[i],
              intercept_stats$r2[i],
              intercept_stats$sd_permil[i],
              intercept_stats$min_permil[i],
              intercept_stats$max_permil[i]))
}
cat(sprintf(
  "\n**Summary:** r = %.2f to %.2f; r² = %.2f to %.2f; SD = %.1f to %.1f ‰; range %.0f to +%.0f ‰.\n\n",
  min(intercept_stats$r), max(intercept_stats$r),
  min(intercept_stats$r2), max(intercept_stats$r2),
  min(intercept_stats$sd_permil), max(intercept_stats$sd_permil),
  min(intercept_stats$min_permil), max(intercept_stats$max_permil)))

header(2, "Spatial slope field — observation-level posterior mean")
cat("| Model | mean | SD | IQR | Range |\n|---|---:|---:|---:|---|\n")
for (i in seq_len(nrow(slope_stats))) {
  cat(sprintf("| %s | %.3f | %.3f | %.3f | %.2f–%.2f |\n",
              slope_stats$model[i],
              slope_stats$mean_slope[i],
              slope_stats$sd_slope[i],
              slope_stats$iqr_slope[i],
              slope_stats$min_slope[i],
              slope_stats$max_slope[i]))
}
cat(sprintf(
  "\n**Across all spatial models:** mean slope %.2f to %.2f, SD ≈ %.2f, IQR ≈ %.2f, range %.2f to %.2f.\n\n",
  min(slope_stats$mean_slope), max(slope_stats$mean_slope),
  median(slope_stats$sd_slope), median(slope_stats$iqr_slope),
  min(slope_stats$min_slope), max(slope_stats$max_slope)))

header(2, "95% prediction interval widths (S2.6.3)")
cat("| Model | Mean width (‰) | SD | Min | Max |\n|---|---:|---:|---:|---:|\n")
for (i in seq_len(nrow(pi_stats))) {
  cat(sprintf("| %s | %.1f | %.1f | %.1f | %.1f |\n",
              pi_stats$model[i],
              pi_stats$mean_width[i],
              pi_stats$sd_width[i],
              pi_stats$min_width[i],
              pi_stats$max_width[i]))
}
sp_widths <- pi_stats |> filter(grepl("_sp$", model))
nsp_widths <- pi_stats |> filter(!grepl("_sp$", model))
cat(sprintf(
  "\n- Spatial mean width: **%.1f ‰**\n- Non-spatial mean width: **%.1f ‰**\n- Spatial reduction: **%s**\n\n",
  mean(sp_widths$mean_width), mean(nsp_widths$mean_width),
  pct(1 - mean(sp_widths$mean_width) / mean(nsp_widths$mean_width))))

header(2, "Nearest-knot density strata (S2.5.2)")
cat("Density is the number of observations within 1,000 chordal km of each observation's nearest fitted knot. RMSE uses the posterior mean fitted value from baseline_env_sp.\n\n")
cat("| Category | Knot-density count | n | In-sample RMSE (per mil) |\n")
cat("|---|---:|---:|---:|\n")
for (i in seq_len(nrow(density_stats))) {
  bounds <- c("<10", "10--50", ">50")[[i]]
  cat(sprintf("| %s | %s | %d | %.1f |\n",
              density_stats$category[i], bounds, density_stats$n[i],
              density_stats$rmse_permil[i]))
}
cat("\n")

header(2, "Residual σ per spatial model")
cat("| Model | σ_residual (‰) | 95% CI |\n|---|---:|---|\n")
for (i in seq_len(nrow(sigma_stats_focus))) {
  cat(sprintf("| %s | %.2f | [%.2f, %.2f] |\n",
              sigma_stats_focus$model[i],
              sigma_stats_focus$sigma_mean[i],
              sigma_stats_focus$sigma_lo[i],
              sigma_stats_focus$sigma_hi[i]))
}
cat(sprintf("\n**Median σ_residual = %.1f ‰**.\n\n", sigma_resid_med))

header(2, "Detection thresholds (§17)")
cat(sprintf("Combined σ_total = √(%.2f² + 3²) = **%.2f ‰**.\n\n",
            sigma_resid_med, sigma_total))
cat("| ρ | Threshold (‰) |\n|---:|---:|\n")
for (i in seq_len(nrow(det_table))) {
  cat(sprintf("| %.1f | %.1f |\n", det_table$rho[i], det_table$threshold_permil[i]))
}
cat(sprintf(
  "\n### Detection-threshold endpoints\n\n- Spatial, ρ=0.9: **%.1f ‰**\n- Non-spatial, ρ=0: **%.1f ‰**; uses σ_baseline = %.2f ‰\n\n",
  spatial_threshold_rho09, nonspatial_threshold_rho0, sigma_baseline_mean))

header(2, "OIPC and measurement uncertainty")
cat(sprintf("- OIPC SE per site (`oipc_se20`): **%.1f to %.1f ‰**\n",
            oipc_se_range[1], oipc_se_range[2]))
cat(sprintf("- δ²H_wax measurement error range: **%.2f to %.1f ‰**\n",
            wax_err_range[1], wax_err_range[2]))
cat(sprintf("- δ²H_wax measurement error quantiles (5/50/95): **%.2f / %.2f / %.2f ‰**\n",
            wax_err_q[1], wax_err_q[2], wax_err_q[3]))

cat("\n---\n\nAll quantities above are reproducible by re-running this script.\n")
sink()

cat("Wrote", OUT_PATH, "\n")
