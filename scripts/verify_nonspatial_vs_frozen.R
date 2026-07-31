#!/usr/bin/env Rscript
# verify_nonspatial_vs_frozen.R
#
# The non-spatial models are refit fresh (spec v2 §3, revised 2026-07-18). They
# have no GP, so the chordal geometry is inert; with the same data/priors/seed/
# executable they SHOULD target the same posterior as the frozen run. That is a
# claim to TEST, not assume — and whole-file md5 equality cannot test it, because
# the refit deposits intentionally carry extra beta_elev columns (and draw/metadata
# ordering may differ). This script does a shared-PARAMETER equivalence check:
# for the parameters present in BOTH deposits it compares per-parameter posterior
# mean and sd within a tolerance.
#
#   Rscript scripts/verify_nonspatial_vs_frozen.R \
#     --refit  results/c2_run_<date>_chordal/model_output \
#     --frozen results/c2_run_20260626/model_output \
#     [--mean-tol-sd 0.1] [--sd-tol 0.1] [--models baseline,baseline_veg,...]
#
# Exit 0 if every model's shared parameters agree within tolerance; nonzero if any
# exceed it (or a deposit is missing/unreadable). This is an explicitly approved
# heuristic (means within a fraction of a posterior SD; SDs within a factor), NOT
# a formal MCSE-based equivalence test. Informational, run post-fit; not a launch gate.

suppressWarnings(suppressMessages(library(posterior)))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  i <- match(flag, args); if (!is.na(i) && i < length(args)) args[[i + 1]] else default
}
refit_dir  <- get_arg("--refit",  "model_output")
frozen_dir <- get_arg("--frozen", "results/c2_run_20260626/model_output")
# Unit-free equivalence criterion (defined up front). Same seed/data/priors/
# executable should reproduce the shared posterior closely, so we require every
# shared STRUCTURAL parameter to agree in mean to within MEAN_TOL_SD posterior
# SDs, and in SD to within a factor SD_TOL. These are compared per-parameter
# against the frozen parameter's own SD, so they are dimensionless (a raw "1%"
# threshold across parameters with different units, as in an earlier version, is
# not meaningful and is removed).
MEAN_TOL_SD <- as.numeric(get_arg("--mean-tol-sd", "0.1"))   # |Δmean| <= 0.1 * sd_frozen
SD_TOL      <- as.numeric(get_arg("--sd-tol", "0.1"))        # sd_refit/sd_frozen in [0.9, 1.1]
models_arg <- get_arg("--models", "baseline,baseline_veg,baseline_env,full,full_interact")
MODELS <- strsplit(models_arg, ",")[[1]]

# Excluded from the structural comparison: transformed/derived quantities (mu,
# scale_weights, log_lik, lp__) and posterior-predictive draws (d2H_rep, y_rep).
# d2H_rep is stochastic (not deterministic), but it is a predictive replicate, not
# a model parameter, so it is still correctly excluded from a parameter-equivalence
# test. What remains are the structural parameters that define the posterior.
GEN_Q <- c("mu", "d2H_rep", "y_rep", "log_lik", "scale_weights", "lp__")
is_gen_q <- function(v) grepl(paste0("^(", paste(GEN_Q, collapse = "|"), ")(\\[|\\.|$)"), v)

summ <- function(pd) {
  x <- posterior::as_draws_df(readRDS(pd))
  vars <- setdiff(names(x), c(".chain", ".iteration", ".draw"))
  data.frame(var = vars,
             mean = vapply(vars, function(v) mean(x[[v]]), numeric(1)),
             sd   = vapply(vars, function(v) sd(x[[v]]),   numeric(1)),
             row.names = NULL)
}

fails <- 0L
for (m in MODELS) {
  rf <- file.path(refit_dir,  m, "posterior_draws.rds")
  fr <- file.path(frozen_dir, m, "posterior_draws.rds")
  cat("== ", m, " ==\n")
  if (!file.exists(rf) || !file.exists(fr)) {
    cat("  FAIL: missing deposit (refit:", file.exists(rf), " frozen:", file.exists(fr), ")\n")
    fails <- fails + 1L; next
  }
  a <- tryCatch(summ(rf), error = function(e) NULL)
  b <- tryCatch(summ(fr), error = function(e) NULL)
  if (is.null(a) || is.null(b)) { cat("  FAIL: unreadable deposit\n"); fails <- fails + 1L; next }
  # Count refit-only variables on the FULL sets (before any subsetting) — this is
  # where beta_elev etc. legitimately appear only in the refit.
  refit_only <- setdiff(a$var, b$var)
  shared <- intersect(a$var, b$var)
  structural <- shared[!is_gen_q(shared)]          # drop generated quantities
  if (!length(structural)) { cat("  FAIL: no shared structural parameters\n"); fails <- fails + 1L; next }
  ai <- a[match(structural, a$var), ]; bi <- b[match(structural, b$var), ]
  # Parameters deterministically fixed in BOTH runs (e.g. full's zeroed PFT
  # interaction columns beta_oipc_x_{tree,shrub,grass}) have sd == 0 in both, so a
  # ratio-based sd test is 0/0 and meaningless. Treat them as inactive: require
  # their means to agree to a small ABSOLUTE tolerance, and exclude them from the
  # sd-ratio gate. The ratio gate applies only to genuinely varying parameters.
  ABS_ZERO <- 1e-8; ABS_MEAN_TOL <- 1e-6
  inactive <- (bi$sd < ABS_ZERO) & (ai$sd < ABS_ZERO)
  inactive_ok <- !any(inactive) || all(abs(ai$mean[inactive] - bi$mean[inactive]) <= ABS_MEAN_TOL)
  act <- !inactive
  if (any(act)) {
    dmean_sd <- abs(ai$mean[act] - bi$mean[act]) / bi$sd[act]   # bi$sd[act] > ABS_ZERO by construction
    sd_ratio <- ai$sd[act] / bi$sd[act]
    mean_ok <- max(dmean_sd) <= MEAN_TOL_SD
    sd_ok   <- all(sd_ratio >= (1 - SD_TOL) & sd_ratio <= (1 + SD_TOL))
    worst <- structural[act][which.max(dmean_sd)]
    cat(sprintf("  shared=%d structural=%d (active=%d, inactive=%d) refit-only(elev etc.)=%d | max |Δmean|/sd=%.3f at '%s' | sd-ratio [%.3f, %.3f]\n",
                length(shared), length(structural), sum(act), sum(inactive), length(refit_only),
                max(dmean_sd), worst, min(sd_ratio), max(sd_ratio)))
  } else {
    mean_ok <- TRUE; sd_ok <- TRUE
    cat(sprintf("  shared=%d structural=%d (all inactive) refit-only=%d\n",
                length(shared), length(structural), length(refit_only)))
  }
  if (mean_ok && sd_ok && inactive_ok) {
    cat(sprintf("  PASS: active means within %.2f sd, active sd within ±%.0f%%, inactive means agree\n",
                MEAN_TOL_SD, 100 * SD_TOL))
  } else {
    cat(sprintf("  FAIL: %s%s%s\n", if (!mean_ok) "mean drift " else "",
                if (!sd_ok) "sd drift " else "", if (!inactive_ok) "inactive-mean drift" else ""))
    fails <- fails + 1L
  }
}

cat("\n")
if (fails == 0) {
  cat("ALL non-spatial refits agree with the frozen run on shared parameters.\n")
} else {
  cat(fails, "model(s) FAILED shared-parameter equivalence — investigate before treating the",
      "non-spatial refits as reproducing the frozen numbers.\n")
  quit(status = 1)
}
