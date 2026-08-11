# regen_precip_space.R
#
# Compute precip-space reconstruction uncertainty and change-detection
# thresholds, propagating slope uncertainty per posterior draw.
#
# Wax-space derivation chain:
#   threshold_wax^(i)(rho_t) = 1.96 *
#     sqrt(2 * sigma_resid^(i)^2 * (1 - rho_t) + 2 * sigma_meas^2)
# (the autocorrelation factor applies only to the residual term;
#  analytical measurement error is independent between samples)
#
# Precip-space (this script):
#   reconstruction_sd_precip^(i)  = sqrt(sigma_resid^(i)^2 + sigma_meas^2)
#                                   / |beta_OIPC^(i)|
#   threshold_precip^(i)(rho_t)   = threshold_wax^(i)(rho_t) / |beta_OIPC^(i)|
#
# For each focus model we report posterior median + 95% credible interval
# across draws of (i) the underlying inputs (beta_OIPC, sigma_resid),
# (ii) the precip-space reconstruction SD, and (iii) the precip-space
# change-detection thresholds at rho_t in {0, 0.5, 0.8, 0.9}.
#
# Run from repo root:
#   Rscript scripts/regen_precip_space.R
# Output: model_analysis/reported_outputs/REGENERATED_NUMBERS_precip_space.md

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(tibble)
})

source("scripts/posterior_helpers.R")

# Set LEAFWAX_OUTPUT_DIR to override the generated-output directory. Only the
# output path changes; every calculation is unchanged.
.output_dir <- Sys.getenv("LEAFWAX_OUTPUT_DIR", unset = "model_analysis/reported_outputs")
dir.create(.output_dir, recursive = TRUE, showWarnings = FALSE)
OUT_PATH   <- file.path(.output_dir, "REGENERATED_NUMBERS_precip_space.md")
SIGMA_MEAS <- 3        # representative analytical SD on wax (per mil)
RHO_GRID   <- c(0, 0.5, 0.8, 0.9)
FOCUS_MODELS <- c("baseline", "baseline_sp", "baseline_env_sp", "full_sp", "full_interact_sp")

q025 <- function(x) unname(quantile(x, 0.025))
q500 <- function(x) unname(quantile(x, 0.500))
q975 <- function(x) unname(quantile(x, 0.975))

per_draw_summaries <- function(model) {
  pd  <- load_draws(model)
  cfg <- load_config(model)
  d2H_sd <- cfg$scaling_params$d2H_sd

  beta_oipc_draws <- slope_model_to_physical(
    as.numeric(as_draws_matrix(subset_draws(pd, variable = "beta_oipc"))),
    cfg$scaling_params
  )
  sigma_draws_std <- as.numeric(as_draws_matrix(subset_draws(pd, variable = "sigma")))
  sigma_draws_pm  <- sigma_draws_std * d2H_sd  # per mil

  # Per draw: total wax-space SD (used only for single-point
  # reconstruction SD), and the autocorrelation-corrected change
  # detection threshold.
  sigma_total_wax <- sqrt(sigma_draws_pm^2 + SIGMA_MEAS^2)
  recon_sd_precip <- sigma_total_wax / abs(beta_oipc_draws)
  threshold_precip <- lapply(RHO_GRID, function(rho) {
    var_diff_wax <- 2 * sigma_draws_pm^2 * (1 - rho) + 2 * SIGMA_MEAS^2
    1.96 * sqrt(var_diff_wax) / abs(beta_oipc_draws)
  })
  names(threshold_precip) <- paste0("rho_", RHO_GRID)

  list(
    model = model,
    beta_oipc = beta_oipc_draws,
    sigma_resid_pm = sigma_draws_pm,
    sigma_total_wax = sigma_total_wax,
    recon_sd_precip = recon_sd_precip,
    threshold_precip = threshold_precip
  )
}

summarise_one <- function(x, label) {
  tibble(
    quantity = label,
    median  = q500(x),
    lo95    = q025(x),
    hi95    = q975(x)
  )
}

build_model_block <- function(s) {
  rows <- bind_rows(
    summarise_one(s$beta_oipc,        "beta_OIPC (global slope)"),
    summarise_one(s$sigma_resid_pm,   "sigma_residual (per mil, wax)"),
    summarise_one(s$sigma_total_wax,  "sigma_total wax (sqrt(sigma_r^2 + 3^2))"),
    summarise_one(s$recon_sd_precip,  "reconstruction SD (precip, single point)"),
    summarise_one(s$threshold_precip$rho_0,   "threshold precip @ rho_t = 0"),
    summarise_one(s$threshold_precip$rho_0.5, "threshold precip @ rho_t = 0.5"),
    summarise_one(s$threshold_precip$rho_0.8, "threshold precip @ rho_t = 0.8"),
    summarise_one(s$threshold_precip$rho_0.9, "threshold precip @ rho_t = 0.9")
  )
  rows$model <- s$model
  rows[, c("model", "quantity", "median", "lo95", "hi95")]
}

cat("Loading posteriors and computing per-draw precip-space quantities ...\n")
all_summaries <- lapply(FOCUS_MODELS, function(m) {
  cat("  ", m, "\n")
  s <- per_draw_summaries(m)
  build_model_block(s)
}) |> bind_rows()

cat("Writing", OUT_PATH, "\n")

sink(OUT_PATH)
cat("# Precipitation-space reconstruction uncertainty and detection thresholds\n\n")
cat(sprintf("Model run: `%s`. Analytical SD assumed = %.1f permil (median observed).\n\n",
            RUN_ID, SIGMA_MEAS))
cat("**Method.** For each posterior draw i we compute\n")
cat("`threshold_precip^(i)(rho_t) = 1.96 * sqrt(2 * sigma_resid^(i)^2 * (1 - rho_t) + 2 * sigma_meas^2)`\n")
cat("`                              / |beta_OIPC^(i)|`.\n")
cat("The autocorrelation factor applies only to the residual term; analytical\n")
cat("measurement error is independent between samples. Reported values are\n")
cat("posterior median and 95 percent credible interval across draws,\n")
cat("propagating slope uncertainty into the precip-space scale.\n\n")

for (m in FOCUS_MODELS) {
  cat(sprintf("## %s\n\n", m))
  cat("| Quantity | median | 95 percent CI |\n|---|---:|---|\n")
  rows <- all_summaries |> dplyr::filter(model == m)
  for (i in seq_len(nrow(rows))) {
    cat(sprintf("| %s | %.2f | [%.2f, %.2f] |\n",
                rows$quantity[i], rows$median[i], rows$lo95[i], rows$hi95[i]))
  }
  cat("\n")
}

cat("---\n\n")
cat("## baseline_env_sp summary\n\n")
focus <- all_summaries |> dplyr::filter(model == "baseline_env_sp")
get_med <- function(q) focus$median[focus$quantity == q]
get_ci  <- function(q) sprintf("[%.0f, %.0f]", focus$lo95[focus$quantity == q],
                               focus$hi95[focus$quantity == q])
cat(sprintf("- Reconstruction uncertainty (single point), precipitation scale: median = **%.0f ‰**, 95 percent CI %s.\n",
            get_med("reconstruction SD (precip, single point)"),
            get_ci("reconstruction SD (precip, single point)")))
for (rho in RHO_GRID) {
  q <- sprintf("threshold precip @ rho_t = %s", rho)
  cat(sprintf("- 95 percent detection threshold for Δδ²H_precip, rho_t = %s: median = **%.0f ‰**, 95 percent CI %s.\n",
              rho, get_med(q), get_ci(q)))
}

sink()
cat("Wrote", OUT_PATH, "\n")
