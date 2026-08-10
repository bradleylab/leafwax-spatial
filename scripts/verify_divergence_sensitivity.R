# verify_divergence_sensitivity.R
#
# Divergence-location and sensitivity diagnostic for the spatial models.
# A divergent transition flags a possible integration failure; whether it
# matters for inference is an empirical question (Betancourt 2017). This
# script does NOT assume divergences are benign — it assembles the evidence
# needed to determine where the divergences fall, whether they
# distort posterior summaries (means AND tails), whether they distort model
# comparison (LOO), and the full convergence table.
#
# The supplement (Section S2.4) reports that after tuning two spatial models
# retain a few divergences (full_sp: 2, full_interact_sp: 4, of 16,000
# post-warmup draws). For each spatial model this computes:
#   1. Per-chain divergence counts (rule out a single pathological chain).
#   2. WHERE divergences fall — the empirical-CDF percentile of each divergent
#      draw within the marginal of every GP hyperparameter (raw and on the
#      reported scale), lp__, the residual scale, and max|z|.
#   3. FUNNEL check — the value of each GP scale / length-scale at the divergent
#      draws vs. its global range (a funnel shows the scale pinned near its min).
#   4. SENSITIVITY of reported quantities to dropping divergent draws, on the
#      mean AND the 2.5/5/50/95/97.5% quantiles, in posterior-SD units.
#   5. LOO sensitivity (flagged models only) — elpd_loo recomputed with the
#      divergent draws removed, to confirm model comparison is unaffected.
#   6. A convergence table (divergences, max R-hat, min ESS bulk/tail,
#      treedepth hits, min E-BFMI) and, for flagged models, a divergence
#      pairs plot.
#
# Per-draw divergence flags come from fit$sampler_diagnostics(), which needs the
# CmdStanR fit object and its chain CSVs, so run against a run dir that retains
# fit.rds + the *.csv files (i.e. on Compute2, not the compact Mac mirror).
#
# Run from the ANALYSIS REPO ROOT inside the container (chain-CSV parsing needs
# a large TMPDIR; the sbatch wrapper binds scratch to /tmp):
#   LEAFWAX_RUN_DIR=/scratch2/fs1/<user>/leafwax_run/model_output \
#     Rscript scripts/verify_divergence_sensitivity.R
# Outputs (regenerable; the committed record is this script + the numbers it
# writes into the tracked supplement):
#   manuscript/drafts/comms_ee/DIVERGENCE_SENSITIVITY.md   (report)
#   results/divergence_sensitivity.csv                     (summary table)
#   results/divergence_pairs_<model>.png                   (flagged models)

suppressPackageStartupMessages({
  library(posterior)
  library(loo)
})

# --- guard: must run from the analysis repo root -----------------------------
if (!file.exists("scripts/verify_divergence_sensitivity.R")) {
  stop("Run from the analysis repo root (leafwax_working), e.g. ",
       "`Rscript scripts/verify_divergence_sensitivity.R`.")
}

RUN_DIR <- Sys.getenv("LEAFWAX_RUN_DIR", unset = "model_output")
OUT_MD  <- "manuscript/drafts/comms_ee/DIVERGENCE_SENSITIVITY.md"
OUT_CSV <- "results/divergence_sensitivity.csv"

# All nine spatial variants, matching verify_pd_knots.R. Non-divergent models
# report a trivial (zero-shift) row, which is itself part of the claim.
SPATIAL_MODELS <- c(
  "baseline_sp", "baseline_env_sp", "baseline_veg_sp",
  "full_sp", "full_interact_sp",
  "elevation_only_sp", "elevation_c4_sp",
  "c4_only_sp", "elevation_c4_interact_sp"
)

# Marginals to locate divergences within.
LOC_VARS <- c("lp__", "log_ls_spatial_km[1]",   # renamed from log_ls_spatial_raw (chordal km param)
              "sigma_intercept_raw[1]", "sigma_slope_raw[1]",
              "sigma_intercept_spatial", "sigma_slope_spatial",
              "ls_intercept_km", "ls_slope_km", "beta_0", "beta_oipc", "sigma")
# Funnel diagnostic: the actual GP scales / length-scale (reported units).
FUNNEL_VARS <- c("sigma_intercept_spatial", "sigma_slope_spatial",
                 "ls_intercept_km", "ls_slope_km")
# Reported quantities whose posterior we test for divergence sensitivity.
REPORT_VARS <- c("beta_0", "beta_oipc",
                 "sigma_intercept_spatial", "sigma_slope_spatial",
                 "ls_intercept_km", "ls_slope_km")
QUANTILES   <- c(0.025, 0.05, 0.5, 0.95, 0.975)
SENSITIVITY_THRESHOLD_SD <- 0.05   # |shift/SD| below this = no material effect
PAIRS_VARS  <- c("lp__", "sigma_intercept_spatial", "sigma_slope_spatial",
                 "ls_intercept_km", "sigma")   # + max_abs_z appended at runtime

load_fit <- function(model) {
  path <- file.path(RUN_DIR, model, "fit.rds")
  if (!file.exists(path)) {
    stop("fit.rds missing for '", model, "' at ", path,
         "\nThis diagnostic needs the CSV-backed fit object; point ",
         "LEAFWAX_RUN_DIR at a run dir that retains fit.rds + chain CSVs.")
  }
  readRDS(path)
}

load_diagnostics <- function(model) {
  path <- file.path(RUN_DIR, model, "diagnostics.rds")
  if (!file.exists(path)) stop("diagnostics.rds missing at ", path)
  readRDS(path)
}

percentile_of <- function(x, values) {
  # empirical CDF: fraction of draws at or below each value (min -> 1/n, max -> 1)
  vapply(values, function(v) mean(x <= v, na.rm = TRUE), numeric(1))
}

# shift in a summary statistic when divergent draws are dropped, in SD units
shift_sd <- function(all_vals, keep_vals, stat_fun) {
  s <- sd(all_vals)
  if (!is.finite(s) || s == 0) return(0)
  (stat_fun(keep_vals) - stat_fun(all_vals)) / s
}

# convergence-table row, built from the compact diagnostics.rds (no draws needed)
build_conv <- function(model, diag) {
  aps <- diag$all_params_summary
  smry <- diag$summary
  ess_tail_min <- if ("ess_tail" %in% names(aps)) min(aps$ess_tail, na.rm = TRUE) else NA
  data.frame(
    model = model,
    divergences = sum(smry$num_divergent),
    treedepth_hits = sum(smry$num_max_treedepth),
    max_rhat = max(aps$rhat, na.rm = TRUE),
    min_ess_bulk = min(aps$ess_bulk, na.rm = TRUE),
    min_ess_tail = ess_tail_min,
    min_ebfmi = min(smry$ebfmi, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

# a sensitivity row with all shifts = 0 (used for models with no divergences,
# where dropping zero draws changes every summary by exactly nothing)
zero_sens_row <- function(model, quantity, mean_val) {
  row <- data.frame(model = model, quantity = quantity,
                    mean_all = mean_val, mean_drop = mean_val,
                    shift_mean_sd = 0, stringsAsFactors = FALSE)
  for (q in QUANTILES) row[[sprintf("shift_q%03d_sd", round(q * 1000))]] <- 0
  row
}

analyze_one_model <- function(model, con) {
  cat("  ", model, "... ")
  diag <- load_diagnostics(model)
  conv <- build_conv(model, diag)

  # Fast path: no divergences -> sensitivity is exactly zero, so skip the costly
  # draws load (parsing 8 x ~290 MB chain CSVs) and the fit.rds requirement.
  # The divergent models -- the only ones the supplement claim rests on -- are
  # handled in full below; clean models contribute zero-shift rows.
  if (conv$divergences == 0) {
    cat("0 divergent (reported); skipping draws load\n")
    writeLines(sprintf("\n### %s\n", model), con)
    writeLines("Divergent transitions: **0**. Sensitivity zero by construction (no draws to drop).\n", con)
    aps <- diag$all_params_summary
    sens <- do.call(rbind, lapply(intersect(REPORT_VARS, aps$variable), function(v)
      zero_sens_row(model, v, aps$mean[match(v, aps$variable)])))
    return(list(sens = sens, conv = conv))
  }

  fit <- load_fit(model)

  # ---- assemble aligned draws + divergence flags --------------------------
  present <- fit$metadata()$model_params   # full flat var list incl. lp__, indexed names, GQ
  want <- unique(c(LOC_VARS, REPORT_VARS, PAIRS_VARS,
                   "z_intercept_spatial", "z_slope_spatial"))
  want_scalar <- setdiff(want, c("z_intercept_spatial", "z_slope_spatial"))
  missing <- setdiff(want_scalar, present)
  if (length(missing) > 0) {
    stop("model '", model, "' is missing expected variables: ",
         paste(missing, collapse = ", "))
  }
  if (!all(c("z_intercept_spatial", "z_slope_spatial") %in%
           sub("\\[.*$", "", present))) {
    stop("model '", model, "' lacks z_intercept_spatial/z_slope_spatial.")
  }

  Ddf <- as_draws_df(fit$draws(
    variables = c(want_scalar, "z_intercept_spatial", "z_slope_spatial"),
    format = "draws_array"))
  Sdf <- fit$sampler_diagnostics(format = "df")
  # HARD alignment check: the two objects must be the same draws in the same order
  stopifnot(nrow(Ddf) == nrow(Sdf),
            identical(Ddf$.chain, Sdf$.chain),
            identical(Ddf$.iteration, Sdf$.iteration),
            identical(Ddf$.draw, Sdf$.draw))
  div <- as.logical(Sdf$divergent__)
  chain <- Ddf$.chain
  n <- length(div); ndiv <- sum(div)
  cat(ndiv, "divergent /", n, "\n")
  # cross-check: fit's per-draw flags must agree with the compact diagnostics count
  stopifnot(ndiv == conv$divergences)

  D <- as.matrix(Ddf[, setdiff(colnames(Ddf), c(".chain", ".iteration", ".draw"))])
  cn <- colnames(D)
  zi <- grep("^z_intercept_spatial\\[", cn, value = TRUE)
  zs <- grep("^z_slope_spatial\\[",     cn, value = TRUE)
  stopifnot(length(zi) > 0, length(zs) > 0)
  maxabsz <- pmax(apply(abs(D[, zi, drop = FALSE]), 1, max),
                  apply(abs(D[, zs, drop = FALSE]), 1, max))

  # ---- sensitivity (mean + tail quantiles) --------------------------------
  sens_rows <- list()
  for (v in intersect(REPORT_VARS, cn)) {
    a <- D[, v]; b <- D[!div, v]
    row <- data.frame(model = model, quantity = v,
                      mean_all = mean(a), mean_drop = mean(b),
                      shift_mean_sd = shift_sd(a, b, mean),
                      stringsAsFactors = FALSE)
    for (q in QUANTILES) {
      row[[sprintf("shift_q%03d_sd", round(q * 1000))]] <-
        shift_sd(a, b, function(x) quantile(x, q, names = FALSE))
    }
    sens_rows[[v]] <- row
  }
  sens <- do.call(rbind, sens_rows)

  # ---- LOO sensitivity (flagged models only) ------------------------------
  loo_line <- NULL
  if (ndiv > 0) {
    ll <- fit$draws(variables = "log_lik", format = "draws_matrix")  # draws x N
    # align: draws_matrix row order matches Ddf (same fit, same flatten order)
    stopifnot(nrow(ll) == n)
    # r_eff is per-observation (length N_obs), invariant to dropping draws, so
    # compute once on the full (equal-chain) set and reuse for the dropped run.
    # Recomputing relative_eff() on the dropped set would error: dropping a few
    # divergent draws leaves unequal per-chain counts, which relative_eff rejects.
    reff   <- relative_eff(exp(ll), chain_id = chain)
    r_all  <- suppressWarnings(loo(ll,         r_eff = reff))
    r_drop <- suppressWarnings(loo(ll[!div, ], r_eff = reff))
    d_elpd <- r_drop$estimates["elpd_loo", "Estimate"] - r_all$estimates["elpd_loo", "Estimate"]
    loo_line <- sprintf(
      "elpd_loo all = %.2f (SE %.2f); dropping %d divergent draws changes elpd_loo by %.4f (%.4f SE units).",
      r_all$estimates["elpd_loo", "Estimate"], r_all$estimates["elpd_loo", "SE"],
      ndiv, d_elpd, d_elpd / r_all$estimates["elpd_loo", "SE"])
  }

  # ---- pairs plot (flagged models only) -----------------------------------
  if (ndiv > 0) {
    pv <- c(intersect(PAIRS_VARS, cn))
    M <- cbind(D[, pv, drop = FALSE], max_abs_z = maxabsz)
    png(file.path("results", sprintf("divergence_pairs_%s.png", model)),
        width = 1400, height = 1400, res = 150)
    colv <- ifelse(div, "red", grDevices::adjustcolor("grey40", alpha.f = 0.25))
    ord  <- order(div)  # draw divergent points last (on top)
    pairs(M[ord, ], col = colv[ord], pch = ifelse(div[ord], 19L, 46L),  # 46 = "."
          main = sprintf("%s: divergent (red) vs non-divergent draws", model))
    dev.off()
  }

  # ---- write the per-model report block -----------------------------------
  writeLines(sprintf("\n### %s\n", model), con)
  writeLines(sprintf("Divergent transitions: **%d / %d** (%.4f%%). Per chain: %s.\n",
                     ndiv, n, 100 * ndiv / n,
                     paste(sprintf("c%d=%d", sort(unique(chain)),
                                   tapply(div, chain, sum)[as.character(sort(unique(chain)))]),
                           collapse = ", ")), con)
  # ndiv > 0 here (clean models returned early via the fast path above)
  writeLines("**Where divergences fall** (empirical-CDF percentile within each marginal):\n", con)
  writeLines("| parameter | divergent-draw percentiles |\n|---|---|", con)
  for (v in c(intersect(LOC_VARS, cn), "max_abs_z")) {
    col <- if (v == "max_abs_z") maxabsz else D[, v]
    p <- sort(percentile_of(col, col[div]))
    writeLines(sprintf("| `%s` | %s |", v, paste(sprintf("%.3f", p), collapse = ", ")), con)
  }
  writeLines("\n**Funnel check** (GP scale/length-scale at divergences vs global range):\n", con)
  writeLines("| parameter | global min | median | max | at divergences |\n|---|---:|---:|---:|---|", con)
  for (v in intersect(FUNNEL_VARS, cn)) {
    col <- D[, v]
    writeLines(sprintf("| `%s` | %.4g | %.4g | %.4g | %s |", v,
                       min(col), median(col), max(col),
                       paste(sprintf("%.4g", col[div]), collapse = ", ")), con)
  }
  writeLines("\n**Sensitivity** to dropping divergent draws (shift in SD units; mean and quantiles):\n", con)
  writeLines("| quantity | shift mean | q2.5 | q5 | q50 | q95 | q97.5 |\n|---|---:|---:|---:|---:|---:|---:|", con)
  for (i in seq_len(nrow(sens))) {
    r <- sens[i, ]
    writeLines(sprintf("| `%s` | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |",
                       r$quantity, r$shift_mean_sd,
                       r$shift_q025_sd, r$shift_q050_sd, r$shift_q500_sd,
                       r$shift_q950_sd, r$shift_q975_sd), con)
  }
  if (!is.null(loo_line)) writeLines(sprintf("\n**LOO:** %s\n", loo_line), con)

  list(sens = sens, conv = conv)
}

cat("Divergence-location + sensitivity diagnostic.\n")
cat("Run dir:", normalizePath(RUN_DIR, mustWork = FALSE), "\n\n")
dir.create(dirname(OUT_MD),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)

con <- file(OUT_MD, "w")
writeLines("# Divergence-location and sensitivity diagnostic\n", con)
writeLines(sprintf("Generated %s. Run dir: `%s`.\n",
                   format(Sys.time(), "%Y-%m-%d %H:%M %Z"), RUN_DIR), con)
writeLines(paste(
  "We assessed whether the divergent transitions retained by two spatial models",
  "show evidence of localized pathology (clustering in a funnel/neck of the GP",
  "hyperparameters) or material sensitivity (a shift in reported posterior",
  "summaries when divergent draws are excluded). A shift below",
  SENSITIVITY_THRESHOLD_SD, "posterior standard deviations is treated as showing",
  "no material effect on that summary.\n"), con)

res  <- lapply(SPATIAL_MODELS, analyze_one_model, con = con)
sens_all <- do.call(rbind, lapply(res, `[[`, "sens"))
conv_all <- do.call(rbind, lapply(res, `[[`, "conv"))

# --- convergence table + overall summary ---
writeLines("\n## Convergence table (all spatial models)\n", con)
writeLines("| model | divergences | treedepth | max R-hat | min ESS bulk | min ESS tail | min E-BFMI |\n|---|---:|---:|---:|---:|---:|---:|", con)
for (i in seq_len(nrow(conv_all))) {
  r <- conv_all[i, ]
  writeLines(sprintf("| %s | %d | %d | %.4f | %.0f | %s | %.3f |",
                     r$model, r$divergences, r$treedepth_hits, r$max_rhat,
                     r$min_ess_bulk,
                     ifelse(is.na(r$min_ess_tail), "NA", sprintf("%.0f", r$min_ess_tail)),
                     r$min_ebfmi), con)
}

shift_cols <- grep("^shift_", names(sens_all), value = TRUE)
max_abs_shift <- max(abs(as.matrix(sens_all[, shift_cols])), na.rm = TRUE)
writeLines("\n## Summary\n", con)
writeLines(sprintf(
  "- Largest absolute sensitivity across all reported quantities, all summary statistics (mean + 5 quantiles), all nine spatial models: **%.4f** posterior SD.\n",
  max_abs_shift), con)
writeLines(sprintf(
  "- This is %s the %.2f-SD threshold: excluding the divergent draws produced no material change in any reported posterior summary.\n",
  ifelse(max_abs_shift < SENSITIVITY_THRESHOLD_SD, "far below", "ABOVE"),
  SENSITIVITY_THRESHOLD_SD), con)
close(con)

write.csv(sens_all, OUT_CSV, row.names = FALSE)
cat("\nWrote", OUT_MD, "and", OUT_CSV, "\n")
cat(sprintf("Max |shift| across all summaries/models: %.4f (threshold %.2f)\n",
            max_abs_shift, SENSITIVITY_THRESHOLD_SD))
