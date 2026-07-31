# phase_c_numbers.R
#
# Recompute the manuscript/supplement numbers that were introduced or revised
# during the chordal-metric re-trace and were NOT already covered by
# regen_manuscript_numbers.R (real-data summaries) or regen_tables_v10.R
# (Tables 1-5 / S8). Everything here derives from the committed chordal run and
# its Tier-C simulation outputs; nothing is transcribed from a run log.
#
# Consolidates five ad-hoc computations used while editing the draft:
#   A. Parameter-recovery 95% CIs           (supplement S2.6.1; main Results)
#   B. Confounding achieved correlations    (Table S4; supplement S2.6.2)
#   C. Intercept-OIPC correlation range     (main Results; supplement S2.6.2)
#   D. Vegetation main-effect 95% CIs       (Table 4; main C4 discussion)
#   E. Detection-threshold constant + grid  (Table S9)
#
# Run from repo root, against the chordal run:
#   LEAFWAX_RUN_DIR=results/c2_run_20260728_chordal/model_output \
#     Rscript scripts/phase_c_numbers.R
#
# Output (only): retrace/chordal/PHASE_C_NUMBERS.md (override the output directory
# with LEAFWAX_RETRACE_OUT_DIR). Every INPUT is resolved from LEAFWAX_RUN_DIR
# alone, so the output redirect never changes a value.
#
# Reproducibility: deterministic. Every quantity is a posterior summary of saved
# draws or a closed-form function of saved run metadata. No inference, no RNG.

suppressPackageStartupMessages({
  library(posterior)
})

source("scripts/posterior_helpers.R")  # resolves LEAFWAX_RUN_DIR -> APRIL_RUN

# Both the frozen great-circle run and the chordal run carry Tier-C scenario
# dirs, so a bare-directory check would let a default (frozen) run silently
# produce great-circle numbers. Require the chordal run manifest — written only
# by the chordal build — at the run root (one level above model_output).
.run_root <- dirname(APRIL_RUN)
if (!file.exists(file.path(.run_root, "RUN_MANIFEST_chordal.rds"))) {
  stop("APRIL_RUN is not a chordal run (no RUN_MANIFEST_chordal.rds at ",
       .run_root, ").",
       "\nSet LEAFWAX_RUN_DIR=results/c2_run_20260728_chordal/model_output")
}
if (!dir.exists(file.path(APRIL_RUN, "confounding_v2_3c_rho00"))) {
  stop("No Tier-C confounding scenarios under ", APRIL_RUN)
}

.retrace_out <- Sys.getenv("LEAFWAX_RETRACE_OUT_DIR", unset = "retrace/chordal")
dir.create(.retrace_out, recursive = TRUE, showWarnings = FALSE)
OUT_PATH <- file.path(.retrace_out, "PHASE_C_NUMBERS.md")

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

recovery <- lapply(c("3a", "3b", "3c"), function(s) {
  scen <- paste0("simrecovery_", s)
  truth <- readRDS(file.path(APRIL_RUN, scen, "true_params.rds"))$beta_oipc
  b <- draws_vec(scen, "beta_oipc")
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
# 6b_spatial_confounding_simulation.R). The old great-circle anchor (rho_c=0.45)
# and the calibrated empirical knob map to achieved correlations through the same
# closed form 6b uses; its three inputs (r0, sigma_int_std, sd(z_indep)) are all
# recoverable from the saved scenario metadata (sd(z_indep) is identified by the
# rho03 point and cross-checked at rho05), so no log constant is transcribed.

conf_meta <- function(s) {
  readRDS(file.path(APRIL_RUN, paste0("confounding_v2_3c_", s),
                    "experiment_metadata.rds"))
}
conf <- do.call(rbind, lapply(c("rho00", "empirical", "rho03", "rho05"), function(s) {
  m <- conf_meta(s)
  data.frame(scenario = s, knob = m$rho, achieved_r = m$rho_achieved)
}))

r0 <- conf_meta("rho00")$rho_achieved   # achieved r at knob 0 = cor(z_indep, oipc)
sigma_int <- conf_meta("rho00")$sigma_z # SD of the standardized spatial intercept

# achieved r as a function of knob t and sd(z_indep); mirrors 6b's
# achieved_rho_fn() with cor(oipc_std, oipc) = 1 substituted analytically.
achieved_r_fn <- function(t, sz) {
  a <- t * sigma_int
  b <- sqrt(1 - t^2)
  (a + b * r0 * sz) / sqrt(a^2 + b^2 * sz^2 + 2 * a * b * r0 * sz)
}
# Identify sd(z_indep) from the rho03 scenario, then verify against rho05. Read
# the knobs from metadata so the identification tracks the actual fixed scenarios
# even if they are ever re-run at different knob values.
.k03 <- conf_meta("rho03")$rho
.k05 <- conf_meta("rho05")$rho
sd_z <- uniroot(function(sz) achieved_r_fn(.k03, sz) - conf_meta("rho03")$rho_achieved,
                c(0.1, 5))$root
stopifnot(abs(achieved_r_fn(.k05, sd_z) - conf_meta("rho05")$rho_achieved) < 1e-4)

OLD_ANCHOR_KNOB <- 0.45  # superseded great-circle anchor (supplement S2.6.2)
old_anchor_r <- achieved_r_fn(OLD_ANCHOR_KNOB, sd_z)
empirical_knob <- conf_meta("empirical")$rho
empirical_r_check <- achieved_r_fn(empirical_knob, sd_z)

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
# regen_tables_v10.R reports in Table 2, computed here directly from the posterior
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

# Slope rows spanning the confounding-graded range (Table S9). The lowest row is
# the empirical-confounding simulation truth (re-standardized), the middle rows
# bracket the spatially adjusted real-data range, the top row is the non-spatial
# baseline slope.
det_grid <- data.frame(
  slope = c(0.50, 0.60, 0.62, 0.70, 0.78),
  note = c("empirical-confounding simulation truth (re-std 0.504)",
           "spatially adjusted, lower end",
           "fitted baseline_env_sp global slope",
           "spatially adjusted, upper end",
           "non-spatial baseline"),
  stringsAsFactors = FALSE
)
det_grid$threshold_permil <- det_const / det_grid$slope

# WRITE REPORT ----------------------------------------------------------------

sink(OUT_PATH)
header(1, "Phase-c numbers — chordal re-trace")
cat(sprintf("Generated %s. Run dir: `%s`.\n\n",
            format(Sys.time(), "%Y-%m-%d %H:%M %Z"), APRIL_RUN))
cat("Numbers introduced or revised during the chordal re-trace that are not ",
    "produced by `regen_manuscript_numbers.R` or `regen_tables_v10.R`.\n\n", sep = "")

header(2, "A. Parameter recovery (supplement S2.6.1)")
cat("| Scenario | Simulated slope | Posterior median | 95% CI | Contains truth |\n")
cat("|---|---:|---:|---|:--:|\n")
for (i in seq_len(nrow(recovery))) {
  cat(sprintf("| %s | %.3f | %.3f | [%.3f, %.3f] | %s |\n",
              recovery$scenario[i], recovery$truth[i], recovery$median[i],
              recovery$lo[i], recovery$hi[i],
              ifelse(recovery$contains_truth[i], "yes", "no")))
}
cat("\nMain text quotes 3a as 0.65 (0.55--0.75) and 3b as 0.64 (0.54--0.74); ",
    "3c over-recovers to 0.964 [0.880, 1.047], excluding the simulated 0.70.\n\n", sep = "")

header(2, "B. Confounding achieved correlations (Table S4; S2.6.2)")
cat("| Scenario | Knob rho_c | Achieved r |\n|---|---:|---:|\n")
for (i in seq_len(nrow(conf))) {
  cat(sprintf("| %s | %.4f | %.4f |\n",
              conf$scenario[i], conf$knob[i], conf$achieved_r[i]))
}
cat(sprintf("\n- Identified sd(z_indep) = **%.4f** (from rho03; cross-checked at rho05).\n", sd_z))
cat(sprintf("- Old great-circle anchor rho_c = %.2f -> achieved r = **%.4f** (reported 0.71; superseded).\n",
            OLD_ANCHOR_KNOB, old_anchor_r))
cat(sprintf("- Calibrated empirical knob rho_c = %.4f -> achieved r = **%.4f** (reported 0.41).\n\n",
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
cat("\nUnder chordal the shrub main effect resolves positive (excludes 0) — a flip ",
    "from the frozen run, where it did not. (The other chordal C4 flip, the tree ",
    "interaction, is a Table 5 term, not shown in this main-effects block.)\n\n",
    sep = "")

header(2, "E. Detection-threshold constant + grid (Table S9)")
cat(sprintf("Constant = 1.96 * sqrt(2*%.4f^2 + 2*%.1f^2) = **%.2f permil** ",
            sigma_resid, SIGMA_ANALYTICAL, det_const))
cat("(sigma_resid = baseline_env_sp RMSE point estimate, posterior-mean prediction).\n\n")
cat("| Slope | Threshold (permil) | Note |\n|---:|---:|---|\n")
for (i in seq_len(nrow(det_grid))) {
  cat(sprintf("| %.2f | %.0f | %s |\n",
              det_grid$slope[i], det_grid$threshold_permil[i], det_grid$note[i]))
}

cat("\n---\n\nAll quantities above are reproducible by re-running this script ",
    "against the chordal run.\n", sep = "")
sink()

cat("Wrote", OUT_PATH, "\n")
