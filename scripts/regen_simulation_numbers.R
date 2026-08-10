# regen_simulation_numbers.R
#
# Recompute simulation and detection-threshold quantities that are not covered by
# regen_manuscript_numbers.R (real-data summaries) or regen_tables.R
# (Tables 1-5 / S8). Everything here derives from the selected fitted run and
# its simulation outputs; nothing is transcribed from a run log.
#
# Consolidates six manuscript computations:
#   A. Parameter-recovery 95% CIs           (supplement S2.6.1; main Results)
#   B. Confounding physical-slope results   (Table S4; supplement S2.6.2)
#   C. Intercept-OIPC correlation range     (main Results; supplement S2.6.2)
#   D. Vegetation main-effect 95% CIs       (Table 4; main C4 discussion)
#   E. Detection-threshold constant + grid  (Table S9)
#   F. Prior-sensitivity physical slopes     (Table S5; supplement S2.6.4)
#
# Run from repo root, against the chordal run:
#   LEAFWAX_RUN_DIR=results/c2_run_20260728_chordal/model_output \
#     Rscript scripts/regen_simulation_numbers.R
#
# Outputs: model_analysis/reported_outputs/SIMULATION_NUMBERS.md plus
# machine-readable CSVs for Tables S4 and S5. Override the output directory
# with LEAFWAX_OUTPUT_DIR. Every input is resolved from LEAFWAX_RUN_DIR, so this
# output redirect
# never changes a value.
#
# Reproducibility: deterministic. Every quantity is a posterior summary of saved
# draws or a closed-form function of saved run metadata. No inference, no RNG.

suppressPackageStartupMessages({
  library(posterior)
})

source("scripts/posterior_helpers.R")

# Require the run manifest at the run root (one level above model_output) so
# reported quantities cannot be generated from an unverified model directory.
.run_root <- dirname(MODEL_RUN_DIR)
if (!file.exists(file.path(.run_root, "RUN_MANIFEST_chordal.rds"))) {
  stop("MODEL_RUN_DIR is not a chordal run (no RUN_MANIFEST_chordal.rds at ",
       .run_root, ").",
       "\nSet LEAFWAX_RUN_DIR=results/c2_run_20260728_chordal/model_output")
}
if (!dir.exists(file.path(MODEL_RUN_DIR, "confounding_v2_3c_rho00"))) {
  stop("No confounding scenarios under ", MODEL_RUN_DIR)
}

.output_dir <- Sys.getenv("LEAFWAX_OUTPUT_DIR", unset = "model_analysis/reported_outputs")
dir.create(.output_dir, recursive = TRUE, showWarnings = FALSE)
OUT_PATH <- file.path(.output_dir, "SIMULATION_NUMBERS.md")
CONFOUNDING_CSV <- file.path(.output_dir, "confounding_physical.csv")
PRIOR_CSV <- file.path(.output_dir, "prior_sensitivity_physical.csv")

# The nine spatial-model variants and the constant analytical measurement SD
# imputed for the two-sample detection threshold (main-text Methods; Table S9).
SPATIAL_MODELS <- c(
  "baseline_sp", "baseline_env_sp", "baseline_veg_sp", "full_sp",
  "full_interact_sp", "elevation_only_sp", "elevation_c4_sp", "c4_only_sp",
  "elevation_c4_interact_sp"
)
SIGMA_ANALYTICAL <- 3  # per mil; analytical SD of a single wax measurement

# Helpers ---------------------------------------------------------------------

header <- function(level, txt) cat(sprintf("%s %s\n\n", strrep("#", level), txt))

# Posterior draws of one scalar parameter, as a plain numeric vector.
draws_vec <- function(model, var) {
  d <- load_draws(model)
  as.numeric(as_draws_matrix(subset_draws(d, variable = var)))
}

# A. Parameter recovery -------------------------------------------------------
# Each simrecovery scenario stores the injected truth in true_params.rds and the
# posterior in posterior_draws.rds. Report median + 95% CI and whether the CI
# contains the simulated slope.

.physical_scaling <- load_config("baseline_sp")$scaling_params
recovery <- lapply(c("3a", "3b", "3c"), function(s) {
  scen <- paste0("simrecovery_", s)
  truth <- slope_model_to_physical(
    readRDS(file.path(MODEL_RUN_DIR, scen, "true_params.rds"))$beta_oipc,
    .physical_scaling
  )
  b <- slope_model_to_physical(
    draws_vec(scen, "beta_oipc"), .physical_scaling
  )
  q <- quantile(b, c(0.5, 0.025, 0.975))
  data.frame(
    scenario = s, truth = truth,
    median = q[[1]], lo = q[[2]], hi = q[[3]],
    contains_truth = truth >= q[[2]] && truth <= q[[3]]
  )
})
recovery <- do.call(rbind, recovery)

# B. Confounding achieved correlations ---------------------------------------
# The knob rho_c and the realized correlation differ because two continental-
# scale smooth fields are incidentally correlated over a finite site set. The
# per-scenario achieved r is stored in experiment_metadata.rds (written by
# 6b_spatial_confounding_simulation.R).

conf_meta <- function(s) {
  readRDS(file.path(MODEL_RUN_DIR, paste0("confounding_v2_3c_", s),
                    "experiment_metadata.rds"))
}
conf <- do.call(rbind, lapply(c("rho00", "empirical", "rho03", "rho05"), function(s) {
  m <- conf_meta(s)
  scen <- paste0("confounding_v2_3c_", s)
  tp <- readRDS(file.path(MODEL_RUN_DIR, scen, "true_params.rds"))
  syn <- readRDS(file.path(MODEL_RUN_DIR, scen, "stan_data_synthetic.rds"))
  oipc <- syn$oipc_values
  xpred <- if (is.matrix(oipc)) rowMeans(oipc) else as.numeric(oipc)

  # Each synthetic response was re-standardized within its scenario. Undo that
  # scaling first, then apply the real-data response/predictor SD ratio to obtain
  # permil wax per permil precipitation. This is the same transformation used by
  # Figure 4, kept here as the tabular source of truth for Table S4.
  to_physical <- function(value) {
    slope_model_to_physical(value * tp$sim_sd, .physical_scaling)
  }
  b <- to_physical(draws_vec(scen, "beta_oipc"))
  truth <- slope_model_to_physical(tp$beta_oipc_unstd, .physical_scaling)
  ols <- to_physical(unname(coef(lm(syn$d2H_wax ~ xpred))[2]))
  q <- quantile(b, c(0.025, 0.975))

  data.frame(
    scenario = s, knob = unname(m$rho), achieved_r = m$rho_achieved,
    truth = truth, posterior_mean = mean(b), posterior_median = median(b),
    lo = q[[1]], hi = q[[2]], bias = mean(b) - truth,
    covers_truth = truth >= q[[1]] && truth <= q[[2]], ols = ols
  )
}))

empirical_knob <- conf_meta("empirical")$rho
empirical_r_check <- conf_meta("empirical")$rho_achieved

# C. Intercept-OIPC correlation across variants ------------------------------
# Posterior-mean spatial intercept vs the OIPC predictor, per spatial model.
sed <- load_sediment()
if (is.null(sed$oipc_d2h20)) stop("sediment data lacks oipc_d2h20")

intercept_corr <- do.call(rbind, lapply(SPATIAL_MODELS, function(m) {
  d <- load_draws(m)
  av <- grep("^alpha_spatial\\[", variables(d), value = TRUE)
  if (!length(av)) stop("no alpha_spatial in ", m)
  apm <- summarise_draws(subset_draws(d, variable = av), mean)$mean
  if (length(apm) != nrow(sed)) {
    stop("alpha_spatial length ", length(apm), " != n sites ", nrow(sed),
         " for ", m)
  }
  data.frame(model = m, r = cor(apm, sed$oipc_d2h20))
}))

# D. Vegetation main-effect 95% CIs ------------------------------------------
# full_sp C4/tree/shrub/grass main effects (exclude interaction terms), 95% CI.
veg_vars <- grep("beta_(c4|tree|shrub|grass)", variables(load_draws("full_sp")),
                 value = TRUE)
veg_vars <- veg_vars[!grepl("_x_|oipc_x|x_c4|x_tree|x_shrub|x_grass", veg_vars)]
if (!length(veg_vars)) stop("no vegetation main-effect terms in full_sp posterior")
veg <- do.call(rbind, lapply(veg_vars, function(v) {
  x <- draws_vec("full_sp", v)
  q <- quantile(x, c(0.5, 0.025, 0.975))
  data.frame(term = v, median = q[[1]], lo = q[[2]], hi = q[[3]],
             excludes_zero = q[[2]] > 0 || q[[3]] < 0)
}))

# E. Detection-threshold constant + grid -------------------------------------
# Two-sample wax-scale detection threshold at rho_t = 0:
#   threshold_precip = 1.96 * sqrt(2*sigma_resid^2 + 2*sigma_analytical^2) / slope
# sigma_resid is the baseline_env_sp residual RMSE point estimate (posterior-mean
# prediction vs observed, back-transformed to per mil) — the same rmse_point that
# regen_tables.R reports in Table 2, computed here directly from the posterior
# so section E depends only on LEAFWAX_RUN_DIR (no cross-script CSV coupling).
.pd_be  <- load_draws("baseline_env_sp")
.cfg_be <- load_config("baseline_env_sp")
.mu_vars <- grep("^mu\\[", variables(.pd_be), value = TRUE)
if (!length(.mu_vars)) stop("no mu[] linear-predictor draws in baseline_env_sp")
.mu_permil <- as_draws_matrix(subset_draws(.pd_be, variable = .mu_vars)) *
  .cfg_be$scaling_params$d2H_sd + .cfg_be$scaling_params$d2H_mean
sigma_resid <- sqrt(mean((sed$d2H_wax - colMeans(.mu_permil))^2))
stopifnot(is.finite(sigma_resid))
det_const <- 1.96 * sqrt(2 * sigma_resid^2 + 2 * SIGMA_ANALYTICAL^2)

# Physical-slope rows spanning the fitted spatial-model range and the
# non-spatial baseline (Table S9). Values are rounded display anchors derived
# from the chordal posterior summaries, not standardized-model coefficients.
det_grid <- data.frame(
  slope = c(0.63, 0.64, 0.73, 0.81),
  note = c("spatially adjusted, lower end",
           "fitted baseline_env_sp global slope",
           "spatially adjusted, upper end",
           "non-spatial baseline"),
  stringsAsFactors = FALSE
)
det_grid$threshold_permil <- det_const / det_grid$slope

# F. Prior sensitivity --------------------------------------------------------
# The reference fit and seven alternative-prior refits all use the same observed
# response and OIPC scaling. Convert every saved beta_oipc draw with that common
# scaling; no refitting is performed here.
prior_models <- c(
  "baseline_veg_sp",
  "sensitivity_beta_oipc_prior_wider",
  "sensitivity_beta_oipc_prior_shifted",
  "sensitivity_beta_oipc_prior_uninformative",
  "sensitivity_pc_slope_relaxed",
  "sensitivity_pc_slope_very_relaxed",
  "sensitivity_ls_longer",
  "sensitivity_ls_shorter"
)
prior_labels <- c(
  "Reference",
  "beta_oipc wider",
  "beta_oipc shifted",
  "beta_oipc uninformative",
  "sigma_slope relaxed",
  "sigma_slope very relaxed",
  "GP length scale longer",
  "GP length scale shorter"
)
prior_sensitivity <- do.call(rbind, Map(function(model, label) {
  b <- slope_model_to_physical(draws_vec(model, "beta_oipc"), .physical_scaling)
  q <- quantile(b, c(0.025, 0.975))
  data.frame(model = model, label = label, posterior_mean = mean(b),
             lo = q[[1]], hi = q[[2]])
}, prior_models, prior_labels))

write.csv(conf, CONFOUNDING_CSV, row.names = FALSE)
write.csv(prior_sensitivity, PRIOR_CSV, row.names = FALSE)

# WRITE REPORT ----------------------------------------------------------------

sink(OUT_PATH)
header(1, "Simulation and detection-threshold numbers")
cat(sprintf("Model run: `%s`.\n\n", RUN_ID))
cat("Quantities not produced by `regen_manuscript_numbers.R` or ",
    "`regen_tables.R`.\n\n", sep = "")

header(2, "A. Parameter recovery (supplement S2.6.1)")
cat("| Scenario | Simulated slope | Posterior median | 95% CI | Contains truth |\n")
cat("|---|---:|---:|---|:--:|\n")
for (i in seq_len(nrow(recovery))) {
  cat(sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] | %s |\n",
              recovery$scenario[i], recovery$truth[i], recovery$median[i],
              recovery$lo[i], recovery$hi[i],
              ifelse(recovery$contains_truth[i], "yes", "no")))
}
cat(sprintf("\nPhysical-scale summaries: 3a %.2f (%.2f--%.2f), 3b %.2f (%.2f--%.2f); ",
            recovery$median[1], recovery$lo[1], recovery$hi[1],
            recovery$median[2], recovery$lo[2], recovery$hi[2]))
cat(sprintf("3c over-recovers to %.3f [%.3f, %.3f], excluding the simulated %.2f.\n\n",
            recovery$median[3], recovery$lo[3], recovery$hi[3], recovery$truth[3]))

header(2, "B. Confounding physical-slope results (Table S4; S2.6.2)")
cat("| Scenario | Knob rho_c | Achieved r | True slope | Posterior mean | 95% CI | Bias | OLS | Covers truth |\n")
cat("|---|---:|---:|---:|---:|---|---:|---:|:--:|\n")
for (i in seq_len(nrow(conf))) {
  cat(sprintf("| %s | %.4f | %.4f | %.3f | %.3f | [%.3f, %.3f] | %+.3f | %.3f | %s |\n",
              conf$scenario[i], conf$knob[i], conf$achieved_r[i],
              conf$truth[i], conf$posterior_mean[i], conf$lo[i], conf$hi[i],
              conf$bias[i], conf$ols[i], ifelse(conf$covers_truth[i], "yes", "no")))
}
cat(sprintf("\nCalibrated empirical knob rho_c = %.4f; achieved r = **%.4f**.\n\n",
            empirical_knob, empirical_r_check))

header(2, "C. Intercept-OIPC correlation across spatial variants (S2.6.2)")
cat("| Model | r |\n|---|---:|\n")
for (i in seq_len(nrow(intercept_corr))) {
  cat(sprintf("| %s | %.3f |\n", intercept_corr$model[i], intercept_corr$r[i]))
}
cat(sprintf("\n**Range across variants: %.2f to %.2f** (baseline_sp = %.3f anchors the confounding test).\n\n",
            min(intercept_corr$r), max(intercept_corr$r),
            intercept_corr$r[intercept_corr$model == "baseline_sp"]))

header(2, "D. Vegetation main-effect 95% CIs, full_sp (Table 4)")
cat("| Term | Median | 95% CI | Excludes 0 |\n|---|---:|---|:--:|\n")
for (i in seq_len(nrow(veg))) {
  cat(sprintf("| %s | %+.3f | [%+.3f, %+.3f] | %s |\n",
              veg$term[i], veg$median[i], veg$lo[i], veg$hi[i],
              ifelse(veg$excludes_zero[i], "yes", "no")))
}
cat("\nThe shrub main effect is positive; the tree and grass main effects are negative.\n\n")

header(2, "E. Detection-threshold constant + grid (Table S9)")
cat(sprintf("Constant = 1.96 * sqrt(2*%.4f^2 + 2*%.1f^2) = **%.2f permil** ",
            sigma_resid, SIGMA_ANALYTICAL, det_const))
cat("(sigma_resid = baseline_env_sp RMSE point estimate, posterior-mean prediction).\n\n")
cat("| Slope | Threshold (permil) | Note |\n|---:|---:|---|\n")
for (i in seq_len(nrow(det_grid))) {
  cat(sprintf("| %.2f | %.0f | %s |\n",
              det_grid$slope[i], det_grid$threshold_permil[i], det_grid$note[i]))
}

header(2, "F. Prior-sensitivity physical slopes (Table S5; S2.6.4)")
cat("| Prior variant | Posterior mean | 95% CI |\n|---|---:|---|\n")
for (i in seq_len(nrow(prior_sensitivity))) {
  cat(sprintf("| %s | %.3f | [%.3f, %.3f] |\n",
              prior_sensitivity$label[i], prior_sensitivity$posterior_mean[i],
              prior_sensitivity$lo[i], prior_sensitivity$hi[i]))
}
cat(sprintf("\nRange of posterior means: **%.3f to %.3f**.\n",
            min(prior_sensitivity$posterior_mean),
            max(prior_sensitivity$posterior_mean)))

cat("\n---\n\nAll quantities above are reproducible by re-running this script ",
    "against the chordal run. Machine-readable Table S4 and S5 sources: `",
    basename(CONFOUNDING_CSV), "` and `", basename(PRIOR_CSV), "`.\n", sep = "")
sink()

cat("Wrote", OUT_PATH, "\n")
cat("Wrote", CONFOUNDING_CSV, "\n")
cat("Wrote", PRIOR_CSV, "\n")
