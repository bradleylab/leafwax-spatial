#───────────────────────────────────────────────────────────────────────────────
# 3i_spatial_plus_residuals.R
#
# Spatial+ pre-processing (Dupont, Wood & Augustin 2022, Biometrics) for
# confounding-robust slope estimation.
#
# Strategy: for each Spatial+ variant, take the standardized predictor
# matrices that the existing pipeline writes for the corresponding base
# model (oipc_values [N x n_scales] and the four oipc_x_* interaction
# matrices) and replace each column with the residual after fitting
# mgcv::gam(column ~ s(lon, lat, bs="tp", k = K)) with REML. The
# residuals are saved as a new stan_data_<base>_splus.rds file. The
# existing 4d_leaf_wax_spatial_model.stan compiles unchanged — the
# model receives the same matrix shapes, just with the spatially
# smooth component of each predictor stripped out.
#
# Per Dupont et al., we do NOT rescale the residuals. Slopes from
# Spatial+ fits are on the decorrelated predictor and are not
# numerically equal to the original-OIPC slope, but their reduction
# vs. the original slope quantifies the confounding-corrected upper
# bound (original) and lower bound (Spatial+).
#
# Method: thin-plate spline smoother on (lon, lat) in degrees, k=200,
# REML wiggliness selection. Reported diagnostics: EDF, smooth R²,
# raw vs. residual SD per (predictor × scale) pair.
#
# Usage:
#   Rscript 3i_spatial_plus_residuals.R [base_model_name ...]
# With no args: processes all four target spatial+ base models.
#───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(mgcv)
})

source("0_load_config.R")

# The four base spatial models that get a Spatial+ variant
SPLUS_BASE_MODELS <- c("baseline_sp", "baseline_env_sp",
                       "full_sp", "full_interact_sp")

# Smoother basis dimension. CRITICAL methodological choice — see note below.
#
# OIPC is itself an interpolated smooth raster of GNIP precipitation isotopes
# against (lat, lon, elevation). With unbounded k + REML, mgcv recovers OIPC
# almost perfectly from coordinates alone (R² ≥ 0.98 for k ≥ 100), leaving
# residuals dominated by per-pixel noise rather than meaningful predictor
# variation. Spatial+ is then uninformative.
#
# Per Dupont et al. 2022, the smoother used to residualize the predictor
# should match the spatial scale of the confounder being controlled for.
# Our spatial GP has length scale ~3,700 km (ls_intercept_km median across
# spatial models). On 1,128 globally distributed sites, that scale
# corresponds to an effective TPS basis of EDF ≈ 35-45 (EDF ≈ surface_area /
# scale²). Empirical sweep on baseline_sp (oipc_values[, 1]):
#
#   k= 10  EDF= 8.9  smooth R²=0.682  resid_sd=0.60   (too coarse)
#   k= 20  EDF=18.6  smooth R²=0.765  resid_sd=0.51   (lat trend only)
#   k= 30  EDF=28.6  smooth R²=0.866  resid_sd=0.39   (~4,200 km)
#   k= 40  EDF≈ 38   smooth R²≈0.90   resid_sd≈0.34   (~3,650 km — matches GP)
#   k= 50  EDF=47.9  smooth R²=0.924  resid_sd=0.29   (~3,200 km)
#   k=100  EDF=92.2  smooth R²=0.951  resid_sd=0.23   (much finer than GP)
#   k=200  EDF=183   smooth R²=0.978  resid_sd=0.15   (recovers raster)
#
# Default k=40 matches the GP scale; report sensitivity at k=30 and k=50
# in the supplement. Override with env var LEAFWAX_SPLUS_K.
SMOOTHER_K <- as.integer(Sys.getenv("LEAFWAX_SPLUS_K", unset = "40"))
if (is.na(SMOOTHER_K) || SMOOTHER_K < 10 || SMOOTHER_K > 300) {
  stop("LEAFWAX_SPLUS_K must be in [10, 300]; got '",
       Sys.getenv("LEAFWAX_SPLUS_K"), "'")
}
cat("Using SMOOTHER_K =", SMOOTHER_K, "\n")

# Predictor matrix names to residualize. Interactions are skipped when
# absent or all-zero (which the upstream pipeline writes for variants
# that do not include the interaction).
TARGET_COLUMNS <- c("oipc_values",
                    "oipc_x_c4", "oipc_x_tree", "oipc_x_shrub", "oipc_x_grass")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  unknown <- setdiff(args, SPLUS_BASE_MODELS)
  if (length(unknown) > 0) {
    stop("Unknown base model(s): ", paste(unknown, collapse = ", "),
         "\nValid base models: ", paste(SPLUS_BASE_MODELS, collapse = ", "))
  }
  models <- args
} else {
  models <- SPLUS_BASE_MODELS
}

prep_dir <- CONFIG$output_dirs$prepared_data
if (!dir.exists(prep_dir)) {
  stop("prepared_data dir not found: ", prep_dir,
       "\nRun 4b_stan_prep.R first.")
}

residualize_column <- function(y, lon, lat, k) {
  if (length(unique(y)) < 3 || all(y == 0)) {
    return(list(residuals = y, edf = NA_real_, rsq = NA_real_, smooth_used = FALSE))
  }
  fit <- mgcv::gam(y ~ s(lon, lat, bs = "tp", k = k), method = "REML")
  list(
    residuals   = as.numeric(residuals(fit, type = "response")),
    edf         = sum(fit$edf) - 1,            # subtract intercept
    rsq         = summary(fit)$r.sq,
    smooth_used = TRUE
  )
}

residualize_predictor_matrix <- function(M, lon, lat, k, name) {
  if (is.null(M) || (is.numeric(M) && all(M == 0))) {
    return(list(residual = M, diagnostics = data.frame()))
  }
  resid_M <- matrix(NA_real_, nrow(M), ncol(M))
  diag <- data.frame()
  for (s in seq_len(ncol(M))) {
    raw_col <- M[, s]
    out <- residualize_column(raw_col, lon, lat, k)
    raw_sd_s <- sd(raw_col)
    resid_raw <- out$residuals
    resid_raw_sd <- sd(resid_raw)
    # Rescale residual to match the original column's SD. The Stan
    # prior beta_oipc ~ normal(0.8, 0.3) is calibrated for a roughly
    # unit-variance predictor; without rescaling, residuals with sd ≈
    # 0.2 would make the prior strongly informative on a scale that
    # doesn't reflect what we want it to. Rescaling preserves the
    # decorrelation (residual is still orthogonal to coords) while
    # keeping the prior interpretation comparable. The rescale factor
    # is logged for traceability.
    if (out$smooth_used && resid_raw_sd > 0) {
      rescale <- raw_sd_s / resid_raw_sd
      resid_M[, s] <- resid_raw * rescale
    } else {
      rescale <- 1
      resid_M[, s] <- resid_raw
    }
    diag <- rbind(diag, data.frame(
      predictor      = name,
      scale_idx      = s,
      edf            = out$edf,
      smooth_rsq     = out$rsq,
      raw_sd         = raw_sd_s,
      resid_raw_sd   = resid_raw_sd,
      rescale_factor = rescale,
      resid_sd       = sd(resid_M[, s]),
      smooth_used    = out$smooth_used
    ))
  }
  list(residual = resid_M, diagnostics = diag)
}

# Variant tag — when LEAFWAX_SPLUS_K is the default 40, write plain
# `<base>_splus`; otherwise stamp the k value into the variant name so
# k-sensitivity runs don't clobber each other.
SPLUS_TAG <- if (SMOOTHER_K == 40L) "splus" else sprintf("splus_k%d", SMOOTHER_K)
cat("Variant tag: ", SPLUS_TAG, "\n\n", sep = "")

for (base_model in models) {
  variant_name <- paste0(base_model, "_", SPLUS_TAG)
  in_path  <- file.path(prep_dir, paste0("stan_data_", base_model, ".rds"))
  out_path <- file.path(prep_dir, paste0("stan_data_", variant_name, ".rds"))
  cfg_in   <- file.path(prep_dir, paste0("config_",    base_model, ".rds"))
  cfg_out  <- file.path(prep_dir, paste0("config_",    variant_name, ".rds"))

  if (!file.exists(in_path)) {
    cat("skip:", base_model, "(missing", in_path, ")\n"); next
  }
  if (!file.exists(cfg_in)) {
    cat("skip:", base_model, "(missing", cfg_in, ")\n"); next
  }

  cat("===", base_model, "  →  ", variant_name, "===\n")
  sd_obj <- readRDS(in_path)

  # Spatial+ uses the (jittered) lon/lat that the rest of the pipeline
  # uses. These are stored on stan_data by prepare_stan_data().
  if (is.null(sd_obj$longitude) || is.null(sd_obj$latitude)) {
    stop("stan_data for ", base_model,
         " is missing longitude/latitude fields.")
  }
  lon <- sd_obj$longitude
  lat <- sd_obj$latitude

  diagnostics_all <- data.frame()

  for (cn in TARGET_COLUMNS) {
    if (is.null(sd_obj[[cn]])) {
      cat(sprintf("  (skipping %-15s — not in stan_data)\n", cn))
      next
    }
    if (is.numeric(sd_obj[[cn]]) && all(sd_obj[[cn]] == 0)) {
      cat(sprintf("  (skipping %-15s — all zeros, predictor disabled)\n", cn))
      next
    }
    cat(sprintf("  residualizing %-15s  [N=%d × n_scales=%d]\n",
                cn, nrow(sd_obj[[cn]]), ncol(sd_obj[[cn]])))
    out <- residualize_predictor_matrix(sd_obj[[cn]], lon, lat,
                                        SMOOTHER_K, cn)
    sd_obj[[cn]] <- out$residual
    diagnostics_all <- rbind(diagnostics_all, out$diagnostics)
  }

  # Tag the stan_data so downstream code can detect Spatial+ variants
  attr(sd_obj, "spatial_plus") <- TRUE
  attr(sd_obj, "spatial_plus_smoother") <- list(
    family = "mgcv::gam s(lon, lat, bs = 'tp')",
    k      = SMOOTHER_K,
    method = "REML",
    note   = "fit per-column, residuals replace original; not rescaled"
  )
  attr(sd_obj, "spatial_plus_diagnostics") <- diagnostics_all

  saveRDS(sd_obj, out_path)
  cat("  saved:", out_path, "\n")
  cat("  diagnostics:\n")
  print(diagnostics_all, row.names = FALSE, digits = 4)

  # Companion config — copy the base config, flip name + spatial_plus flag,
  # carry the diagnostics. Downstream (4c_fit_models.R) reads this.
  cfg <- readRDS(cfg_in)
  cfg$name <- variant_name
  cfg$spatial_plus <- TRUE
  cfg$spatial_plus_smoother  <- attr(sd_obj, "spatial_plus_smoother")
  cfg$spatial_plus_diagnostics <- diagnostics_all
  cfg$base_model <- base_model
  cfg$splus_smoother_k <- SMOOTHER_K
  cfg$timestamp <- Sys.time()
  saveRDS(cfg, cfg_out)
  cat("  saved:", cfg_out, "\n\n")
}

cat("DONE.\n")
