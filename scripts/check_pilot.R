#!/usr/bin/env Rscript
# check_pilot.R
#
# Pre-launch gate on the converged baseline_sp chordal pilot:
#   (1) convergence — max R-hat, min ESS, divergences, max-treedepth hits, E-BFMI
#       (all from diagnostics.rds);
#   (2) length-scale sanity — the column must be present and finite, the two
#       reported columns (ls_intercept_km / ls_slope_km, both exp of the single
#       shared log_ls_spatial_km) must agree, and the posterior must NOT pile
#       against either configured bound, else the "nominal km re-expression"
#       prior is too tight.
#
# Bounds are read from the prepared Stan data (ls_log_lower / ls_log_upper), NOT
# hardcoded, so this gate tracks config.yaml automatically. Reads saved RDS only
# (no fit, no cmdstanr).
#
# Usage:  Rscript scripts/check_pilot.R [model_output/baseline_sp] [prepared_data/stan_data_baseline_sp.rds]
# If the stan_data path is omitted it is inferred from the model dir name; if it
# cannot be found the script FAILS (bounds must be authoritative, not guessed).

suppressWarnings(suppressMessages(library(posterior)))
args <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(args) >= 1) args[[1]] else "model_output/baseline_sp"

NEAR_FRAC <- 0.02   # "near a bound" = within 2% of the log-scale span
fails <- 0
gate <- function(name, ok, detail = "") {
  if (!isTRUE(ok)) fails <<- fails + 1
  cat(sprintf("  [%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name,
              if (nzchar(detail)) paste0(" — ", detail) else ""))
}

diag_f  <- file.path(dir, "diagnostics.rds")
draws_f <- file.path(dir, "posterior_draws.rds")
if (!file.exists(diag_f) || !file.exists(draws_f)) {
  stop("Missing diagnostics.rds / posterior_draws.rds in ", dir,
       " — run the pilot fit first (sbatch --array=0 slurm/job_fit_chordal.sh).")
}

# ── Resolve the length-scale bounds from the prepared Stan data ────────────────
model <- basename(dir)
sd_path <- if (length(args) >= 2) args[[2]] else
  file.path("prepared_data", paste0("stan_data_", model, ".rds"))
if (!file.exists(sd_path)) {
  stop("Prepared Stan data not found at '", sd_path, "'. Length-scale bounds must ",
       "come from the fitted stan_data (ls_log_lower/ls_log_upper), not a guess. ",
       "Pass the path as the 2nd argument.")
}
sd <- readRDS(sd_path)
if (is.null(sd$ls_log_lower) || is.null(sd$ls_log_upper)) {
  stop("stan_data at '", sd_path, "' has no ls_log_lower/ls_log_upper — not a ",
       "chordal-era prepared dataset.")
}
LOWER_KM <- exp(sd$ls_log_lower)
UPPER_KM <- exp(sd$ls_log_upper)

# Freshness: the diagnostics/draws being gated must post-date the stan_data they
# were (supposedly) fit against, else we would gate a pilot that predates the
# current prepared data.
cat("== Artifact freshness vs stan_data ==\n")
fresh <- file.mtime(diag_f) > file.mtime(sd_path) && file.mtime(draws_f) > file.mtime(sd_path)
gate("diagnostics.rds and posterior_draws.rds newer than stan_data", fresh,
     sprintf("diag %s, draws %s, stan_data %s",
             format(file.mtime(diag_f), "%H:%M:%S"),
             format(file.mtime(draws_f), "%H:%M:%S"),
             format(file.mtime(sd_path), "%H:%M:%S")))

cat("== Convergence (", dir, ") ==\n", sep = "")
d <- readRDS(diag_f)
max_rhat <- d$max_rhat
min_ess  <- d$min_ess_bulk
ndiv     <- sum(d$summary$num_divergent)
ntree    <- if (!is.null(d$summary$num_max_treedepth)) sum(d$summary$num_max_treedepth) else NA_integer_
ebfmi    <- d$summary$ebfmi
# E-BFMI must be PRESENT and finite for every chain. na.rm would let an all-missing
# ebfmi (e.g. a diagnostics object that never recorded it) pass with count 0.
ebfmi_ok <- !is.null(ebfmi) && length(ebfmi) > 0 && all(is.finite(ebfmi))
n_low_ebfmi <- if (ebfmi_ok) sum(ebfmi < 0.2) else NA_integer_
cat(sprintf("  max R-hat = %.3f | min ESS bulk = %.0f | divergences = %d | max-treedepth hits = %s | E-BFMI<0.2 chains = %s\n",
            max_rhat, min_ess, ndiv,
            ifelse(is.na(ntree), "NA", as.character(ntree)),
            ifelse(is.na(n_low_ebfmi), "NA (missing/nonfinite)", as.character(n_low_ebfmi))))
gate("max R-hat < 1.01", is.finite(max_rhat) && max_rhat < 1.01, sprintf("%.3f", max_rhat))
gate("min ESS bulk > 400", is.finite(min_ess) && min_ess > 400, sprintf("%.0f", min_ess))
gate("no divergent transitions", ndiv == 0, sprintf("%d", ndiv))
gate("no max-treedepth hits", isTRUE(ntree == 0), sprintf("%s", ntree))
gate("E-BFMI present & finite for all chains", ebfmi_ok,
     if (ebfmi_ok) "" else "missing or nonfinite E-BFMI")
gate("no low E-BFMI chains (all >= 0.2)", isTRUE(n_low_ebfmi == 0), sprintf("%s", n_low_ebfmi))

cat("== Length-scale sanity ==\n")
draws <- posterior::as_draws_df(readRDS(draws_f))
have_int <- "ls_intercept_km" %in% names(draws)
have_slp <- "ls_slope_km" %in% names(draws)
gate("length-scale column present", have_int || have_slp,
     if (have_int || have_slp) "" else "neither ls_intercept_km nor ls_slope_km found")
if (!have_int && !have_slp) {
  cat("\nCannot evaluate length-scale gate without a length-scale column.\n")
  quit(status = 1)
}

# ls_intercept_km and ls_slope_km are both exp(log_ls_spatial_km[1]) — the single
# shared length scale reported twice (4d_leaf_wax_spatial_model.stan:441-442).
# If both are present they must be identical draw-for-draw.
if (have_int && have_slp) {
  max_abs_diff <- max(abs(draws[["ls_intercept_km"]] - draws[["ls_slope_km"]]))
  gate("ls_intercept_km == ls_slope_km (shared length scale)",
       is.finite(max_abs_diff) && max_abs_diff < 1e-6, sprintf("max|Δ|=%.2e", max_abs_diff))
}

ls_col <- if (have_int) "ls_intercept_km" else "ls_slope_km"
ls <- draws[[ls_col]]
gate("length-scale draws all finite", all(is.finite(ls)),
     sprintf("%d non-finite", sum(!is.finite(ls))))
ls <- ls[is.finite(ls)]

qs <- quantile(ls, c(0.01, 0.05, 0.5, 0.95, 0.99))
cat(sprintf("  %s quantiles (km): 1%%=%.0f 5%%=%.0f 50%%=%.0f 95%%=%.0f 99%%=%.0f\n",
            ls_col, qs[1], qs[2], qs[3], qs[4], qs[5]))
cat(sprintf("  bounds (km, from stan_data): [%.1f, %.1f]\n", LOWER_KM, UPPER_KM))
loglo <- log(LOWER_KM); loghi <- log(UPPER_KM); span <- loghi - loglo
near_lo <- mean((log(ls) - loglo) / span < NEAR_FRAC)
near_hi <- mean((loghi - log(ls)) / span < NEAR_FRAC)
cat(sprintf("  fraction within %.0f%% of lower bound: %.3f | upper bound: %.3f\n",
            100 * NEAR_FRAC, near_lo, near_hi))
gate("posterior not piling at lower bound (<5% of draws near)", near_lo < 0.05, sprintf("%.3f", near_lo))
gate("posterior not piling at upper bound (<5% of draws near)", near_hi < 0.05, sprintf("%.3f", near_hi))

cat("\n")
if (fails == 0) {
  cat("PILOT GATE PASSED — safe to launch the full 17-fit chordal batch.\n")
} else {
  cat(sprintf("PILOT GATE: %d check(s) failed — review before launching the batch.\n", fails))
  cat("If the length-scale piles at a bound, evaluate and document any change to\n",
      "gp_length_scale bounds in config.yaml, re-prepare, and re-run the pilot.\n")
  quit(status = 1)
}
