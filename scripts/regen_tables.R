# regen_tables.R
#
# Regenerate manuscript Tables 1--5 from the model run selected by
# LEAFWAX_RUN_DIR. The reported analysis uses
# results/c2_run_20260728_chordal/model_output.
#
# Run from repo root:
#   Rscript scripts/regen_tables.R
#
# Reads the selected run's posterior draws, diagnostics, LOO objects, prepared
# model data, and model-ready calibration data.
#
# Writes generated .tex and companion .csv files to
# model_analysis/reported_outputs/ by default.
#
# Reproducibility: deterministic — only summarises saved posterior draws.

suppressPackageStartupMessages({
  library(posterior)
  library(dplyr)
  library(tibble)
  library(loo)
})

source("scripts/posterior_helpers.R")
source("scripts/spatial_variance_helpers.R")
source("scripts/table_helpers.R")

# Set LEAFWAX_OUTPUT_DIR to override the generated-output directory. Only the
# output path changes; every calculation is unchanged.
.output_dir <- Sys.getenv("LEAFWAX_OUTPUT_DIR", unset = "model_analysis/reported_outputs")
dir.create(.output_dir, recursive = TRUE, showWarnings = FALSE)
.outpath <- function(default) file.path(.output_dir, basename(default))

ALL_MODELS <- c(
  "baseline", "baseline_sp",
  "baseline_env", "baseline_env_sp",
  "baseline_veg", "baseline_veg_sp",
  "full", "full_sp",
  "full_interact", "full_interact_sp",
  "elevation_only_sp", "elevation_c4_sp",
  "c4_only_sp", "elevation_c4_interact_sp"
)

SPATIAL_MODELS <- ALL_MODELS[grepl("_sp$", ALL_MODELS)]

# Predictor strings used in the second column of Table 1.
PREDICTORS <- c(
  baseline                 = "$\\delta^2$H$_p$",
  baseline_sp              = "$\\delta^2$H$_p$ + GP",
  baseline_env             = "$\\delta^2$H$_p$, elev, precip",
  baseline_env_sp          = "$\\delta^2$H$_p$, elev, precip + GP",
  baseline_veg             = "$\\delta^2$H$_p$, PFT, C4, PFT $\\times$ $\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$",
  baseline_veg_sp          = "$\\delta^2$H$_p$, PFT, C4, PFT $\\times$ $\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$ + GP",
  full                     = "$\\delta^2$H$_p$, PFT, C4, elev, precip",
  full_sp                  = "$\\delta^2$H$_p$, PFT, C4, elev, precip + GP",
  full_interact            = "$\\delta^2$H$_p$, PFT $\\times$ $\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$, elev, precip",
  full_interact_sp         = "$\\delta^2$H$_p$, PFT $\\times$ $\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$, elev, precip + GP",
  elevation_only_sp        = "$\\delta^2$H$_p$, elev + GP",
  elevation_c4_sp          = "$\\delta^2$H$_p$, elev, C4 + GP",
  c4_only_sp               = "$\\delta^2$H$_p$, C4 + GP",
  elevation_c4_interact_sp = "$\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$, elev + GP"
)

# Pretty model names with LaTeX-safe underscores
MODEL_LABELS <- c(
  baseline                 = "baseline",
  baseline_sp              = "baseline$\\_{\\text{sp}}$",
  baseline_env             = "baseline\\_env",
  baseline_env_sp          = "baseline\\_env$\\_{\\text{sp}}$",
  baseline_veg             = "baseline\\_veg",
  baseline_veg_sp          = "baseline\\_veg$\\_{\\text{sp}}$",
  full                     = "full",
  full_sp                  = "full$\\_{\\text{sp}}$",
  full_interact            = "full\\_interact",
  full_interact_sp         = "full\\_interact$\\_{\\text{sp}}$",
  elevation_only_sp        = "elevation\\_only$\\_{\\text{sp}}$",
  elevation_c4_sp          = "elevation\\_c4$\\_{\\text{sp}}$",
  c4_only_sp               = "c4\\_only$\\_{\\text{sp}}$",
  elevation_c4_interact_sp = "elevation\\_c4\\_interact$\\_{\\text{sp}}$"
)

# Helper: format mean with [q025, q975] interval ----------------------------

fmt_meanci <- function(mean_val, lo, hi, digits = 1, blank_dash = TRUE) {
  if (blank_dash && (is.na(mean_val) || is.null(mean_val))) return("-")
  fmt <- paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]")
  sprintf(fmt, mean_val, lo, hi)
}

# Per-draw RMSE, R^2, plus posterior summaries of fixed effects -------------

extract_global_params <- function(model) {
  pd <- load_draws(model)
  cfg <- load_config(model)
  sed <- load_sediment()

  d2H_sd   <- cfg$scaling_params$d2H_sd
  d2H_mean <- cfg$scaling_params$d2H_mean

  vars <- variables(pd)
  mu_vars <- vars[startsWith(vars, "mu[")]
  if (length(mu_vars) == 0) {
    stop("No mu[] variables in posterior_draws for ", model)
  }
  mu_mat <- as_draws_matrix(subset_draws(pd, variable = mu_vars))   # [draws, n]
  mu_permil <- mu_mat * d2H_sd + d2H_mean
  y <- sed$d2H_wax

  # Per-draw RMSE and R^2 of linear-predictor against observed y
  rmse_per_draw <- apply(mu_permil, 1, function(m) sqrt(mean((y - m)^2)))
  r2_per_draw   <- apply(mu_permil, 1, function(m) {
    1 - sum((y - m)^2) / sum((y - mean(y))^2)
  })
  # Point-estimate RMSE / R^2 from posterior-mean predictions
  mu_pmean <- colMeans(mu_permil)
  rmse_point <- sqrt(mean((y - mu_pmean)^2))
  r2_point   <- 1 - sum((y - mu_pmean)^2) / sum((y - mean(y))^2)

  # beta_0 back-transformed: beta_0_permil = beta_0 * d2H_sd + d2H_mean
  beta_0_draws    <- as.numeric(as_draws_matrix(subset_draws(pd, variable = "beta_0")))
  beta_0_permil   <- beta_0_draws * d2H_sd + d2H_mean
  beta_oipc_draws <- as.numeric(as_draws_matrix(subset_draws(pd, variable = "beta_oipc")))
  beta_oipc_physical <- slope_model_to_physical(
    beta_oipc_draws, cfg$scaling_params
  )

  # Lambda integration (always present; lambda_decay)
  lambda_draws <- as.numeric(as_draws_matrix(subset_draws(pd, variable = "lambda_decay")))

  # GP scale only for spatial models
  is_spatial <- grepl("_sp$", model)
  if (is_spatial) {
    ls_draws  <- as.numeric(as_draws_matrix(subset_draws(pd, variable = "ls_intercept_km")))
    sig_slope <- slope_model_to_physical(
      as.numeric(as_draws_matrix(subset_draws(pd, variable = "sigma_slope_spatial"))),
      cfg$scaling_params
    )
    sig_int   <- as.numeric(as_draws_matrix(subset_draws(pd, variable = "sigma_intercept_spatial")))
  } else {
    ls_draws  <- NA_real_
    sig_slope <- NA_real_
    sig_int   <- NA_real_
  }

  q <- function(x, p) unname(quantile(x, p))

  tibble(
    model = model,
    rmse_point  = rmse_point,
    rmse_mean   = mean(rmse_per_draw),
    rmse_lo     = q(rmse_per_draw, 0.025),
    rmse_hi     = q(rmse_per_draw, 0.975),
    r2_point    = r2_point,
    r2_mean     = mean(r2_per_draw),
    r2_lo       = q(r2_per_draw, 0.025),
    r2_hi       = q(r2_per_draw, 0.975),
    beta_0_mean = mean(beta_0_permil),
    beta_0_lo   = q(beta_0_permil, 0.025),
    beta_0_hi   = q(beta_0_permil, 0.975),
    beta_oipc_mean = mean(beta_oipc_physical),
    beta_oipc_lo   = q(beta_oipc_physical, 0.025),
    beta_oipc_hi   = q(beta_oipc_physical, 0.975),
    lambda_mean = mean(lambda_draws),
    lambda_lo   = q(lambda_draws, 0.025),
    lambda_hi   = q(lambda_draws, 0.975),
    gp_scale_mean = if (is_spatial) mean(ls_draws) else NA_real_,
    gp_scale_lo   = if (is_spatial) q(ls_draws, 0.025) else NA_real_,
    gp_scale_hi   = if (is_spatial) q(ls_draws, 0.975) else NA_real_,
    knot_slope_sd_mean = if (is_spatial) mean(sig_slope) else NA_real_,
    # sigma_intercept_spatial is already in per-mil units in this Stan
    # parameterization.
    knot_int_sd_mean_permil = if (is_spatial) mean(sig_int) else NA_real_
  )
}

# Marginal spatial-component comparison (Table 3) --------------------------

extract_variance_decomp <- function(model) {
  if (!grepl("_sp$", model)) return(NULL)
  pd <- load_draws(model)
  stan_data <- load_stan_data(model)
  fitted_oipc <- fitted_oipc_draw_matrix(pd, stan_data, model)
  max_weight_error <- attr(fitted_oipc, "scale_weight_max_sum_error")
  summary <- summarize_spatial_components(pd, fitted_oipc, model)

  # These are marginal response-scale component magnitudes, not an additive
  # partition of total response variance. The intercept and slope contributions
  # covary, so their covariance is retained in the CSV for audit but excluded
  # from the two percentages displayed in Figure 2b.
  tibble(
    model = model,
    n_draws = summary$n_draws,
    n_obs = summary$n_obs,
    marginal_intercept_variance_std2 =
      summary$marginal_intercept_variance_std2,
    marginal_slope_variance_std2 =
      summary$marginal_slope_variance_std2,
    twice_covariance_std2 = summary$twice_covariance_std2,
    realized_spatial_variance_std2 =
      summary$realized_spatial_variance_std2,
    residual_variance_std2 = summary$residual_variance_std2,
    marginal_intercept_share_pct =
      summary$marginal_intercept_share_pct,
    marginal_slope_share_pct = summary$marginal_slope_share_pct,
    realized_spatial_share_of_spatial_plus_residual_pct =
      summary$realized_spatial_share_of_spatial_plus_residual_pct,
    scale_weight_max_sum_error_before_normalization = max_weight_error,
    identity_max_abs_error = summary$identity_max_abs_error
  )
}

# Vegetation coefficients (Table 4) — 95% CI per existing footnote ---------

VEG_VARS <- c(
  C4    = "beta_c4",
  tree  = "beta_tree",
  shrub = "beta_shrub",
  grass = "beta_grass",
  precip = "beta_precip"
)

extract_veg_coefs <- function(model) {
  pd <- load_draws(model)
  vars <- variables(pd)

  q <- function(x, p) unname(quantile(x, p))
  out <- list(model = model)
  for (nm in names(VEG_VARS)) {
    v <- VEG_VARS[[nm]]
    if (v %in% vars) {
      d <- as.numeric(as_draws_matrix(subset_draws(pd, variable = v)))
      out[[paste0(nm, "_mean")]] <- mean(d)
      out[[paste0(nm, "_lo")]]   <- q(d, 0.025)
      out[[paste0(nm, "_hi")]]   <- q(d, 0.975)
    } else {
      out[[paste0(nm, "_mean")]] <- NA_real_
      out[[paste0(nm, "_lo")]]   <- NA_real_
      out[[paste0(nm, "_hi")]]   <- NA_real_
    }
  }
  as_tibble(out)
}

# Interaction coefficients (Table 5) — 95% CI ------------------------------
# delta2H_precip x vegetation interactions; present only in interaction
# model variants. Absent terms return NA and render as dashes.

INTERACT_VARS <- c(
  C4    = "beta_oipc_x_c4",
  tree  = "beta_oipc_x_tree",
  shrub = "beta_oipc_x_shrub",
  grass = "beta_oipc_x_grass"
)

extract_interaction_coefs <- function(model) {
  pd <- load_draws(model)
  vars <- variables(pd)

  q <- function(x, p) unname(quantile(x, p))
  out <- list(model = model)
  for (nm in names(INTERACT_VARS)) {
    v <- INTERACT_VARS[[nm]]
    if (v %in% vars) {
      d <- as.numeric(as_draws_matrix(subset_draws(pd, variable = v)))
      out[[paste0(nm, "_mean")]] <- mean(d)
      out[[paste0(nm, "_lo")]]   <- q(d, 0.025)
      out[[paste0(nm, "_hi")]]   <- q(d, 0.975)
    } else {
      out[[paste0(nm, "_mean")]] <- NA_real_
      out[[paste0(nm, "_lo")]]   <- NA_real_
      out[[paste0(nm, "_hi")]]   <- NA_real_
    }
  }
  out$any_interaction <- any(!is.na(c(out$C4_mean, out$tree_mean,
                                      out$shrub_mean, out$grass_mean)))
  as_tibble(out)
}

# ---------- BUILD ALL DATAFRAMES -------------------------------------------

cat("Loading diagnostics...\n")
diag_tbl <- lapply(ALL_MODELS, function(m) {
  d <- readRDS(file.path(MODEL_RUN_DIR, m, "diagnostics.rds"))
  tibble(
    model = m,
    max_rhat = d$max_rhat,
    min_ess  = d$min_ess_bulk
  )
}) |> bind_rows()

cat("Loading LOO objects...\n")
loo_list <- lapply(ALL_MODELS, function(m) load_loo(m))
names(loo_list) <- ALL_MODELS

loo_tbl <- lapply(ALL_MODELS, function(m) {
  l <- loo_list[[m]]
  tibble(
    model = m,
    looic    = l$estimates["looic", "Estimate"],
    se_looic = l$estimates["looic", "SE"],
    p_eff    = l$estimates["p_loo", "Estimate"],
    n_hi_k   = sum(l$diagnostics$pareto_k > 0.7, na.rm = TRUE)
  )
}) |> bind_rows()

# loo_compare ranks and gives elpd_diff with SE; convert to ΔLOOIC
cat("Running loo_compare...\n")
cmp <- loo_compare(loo_list)
# Rows are sorted by elpd_loo descending. Column elpd_diff is relative to row 1.
# elpd_diff column gives Δelpd; ΔLOOIC = -2 × Δelpd; SE column is se_diff (in elpd).
cmp_df <- as.data.frame(cmp)
cmp_df$model <- rownames(cmp_df)
delta_tbl <- tibble(
  model = cmp_df$model,
  delta_looic = -2 * cmp_df$elpd_diff,
  se_delta    = 2 * cmp_df$se_diff
)
# delta_looic is negative-of-Δelpd-(reference - this); since reference is the
# best model and its elpd_diff is 0 with se_diff 0, this gives ΔLOOIC ≥ 0.

table1_df <- ALL_MODELS |>
  tibble(model = _) |>
  left_join(diag_tbl, by = "model") |>
  left_join(loo_tbl,  by = "model") |>
  left_join(delta_tbl, by = "model") |>
  mutate(
    label       = MODEL_LABELS[model],
    predictors  = PREDICTORS[model],
    max_rhat    = sprintf("%.3f", max_rhat),
    min_ess     = sprintf("%d", as.integer(min_ess)),
    looic_str   = sprintf("%.1f", looic),
    se_str      = sprintf("%.1f", se_looic),
    delta_str   = sprintf("%.1f", delta_looic),
    se_delta_str= sprintf("%.1f", se_delta),
    p_eff_str   = sprintf("%.1f", p_eff),
    n_hi_k_str  = sprintf("%d", n_hi_k)
  )

cat("Extracting global params (Table 2)...\n")
table2_rows <- lapply(ALL_MODELS, extract_global_params) |> bind_rows()

cat("Extracting variance decomposition (Table 3)...\n")
table3_rows <- lapply(SPATIAL_MODELS, extract_variance_decomp) |> bind_rows()

cat("Extracting environmental coefficients (Table 4)...\n")
table4_rows <- lapply(ALL_MODELS, extract_veg_coefs) |> bind_rows()

cat("Extracting interaction coefficients (Table 5)...\n")
table5_rows <- lapply(ALL_MODELS, extract_interaction_coefs) |> bind_rows()

# --------- WRITE TABLE 1 (standalone .tex) ---------------------------------

cat("Writing Table 1...\n")

# Build the data frame with the exact column order Table 1 expects
t1_out <- table1_df |>
  transmute(
    Model         = label,
    Predictors    = predictors,
    `Max $\\hat{R}$` = max_rhat,
    `Min ESS`     = min_ess,
    LOOIC         = looic_str,
    SE            = se_str,
    `$\\Delta$LOOIC` = delta_str,
    `SE `         = se_delta_str,
    `$p_{\\text{eff}}$` = p_eff_str,
    `$n_{\\text{hi-k}}$` = n_hi_k_str
  )

# We can't use emit_standalone_tex directly because Table 1 needs landscape +
# specific p{4.5cm} column for predictors. Build the LaTeX by hand to match
# the required table structure.

table1_path <- .outpath("table1_model_performance.tex")
banner <- c(
  "% Auto-generated by scripts/regen_tables.R.",
  "% Do not hand-edit — regenerate via `Rscript scripts/regen_tables.R`.",
  ""
)

t1_lines <- c(
  banner,
  "\\begin{landscape}",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Model performance metrics for all candidate models}",
  "\\label{tab:model-performance}",
  "\\footnotesize",
  "",
  "\\begin{tabular}{lp{4.5cm}cccccccc}",
  "\\toprule",
  "Model & Predictors & Max $\\hat{R}$ & Min ESS & LOOIC & SE & $\\Delta$LOOIC & SE & $p_{\\text{eff}}$ & $n_{\\text{hi-k}}$\\\\",
  "\\midrule"
)

for (i in seq_len(nrow(t1_out))) {
  row <- t1_out[i, ]
  t1_lines <- c(t1_lines, paste0(
    paste(c(row$Model, row$Predictors, row$`Max $\\hat{R}$`, row$`Min ESS`,
            row$LOOIC, row$SE, row$`$\\Delta$LOOIC`, row$`SE `,
            row$`$p_{\\text{eff}}$`, row$`$n_{\\text{hi-k}}$`),
          collapse = " & "),
    "\\\\"
  ))
  # Separate groups of five models for readability.
  if (i %in% c(5, 10) && i < nrow(t1_out)) {
    t1_lines <- c(t1_lines, "\\addlinespace")
  }
}

t1_lines <- c(t1_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{2mm}",
  "\\footnotesize",
  "\\textit{Note:} $\\hat{R}$ = Gelman-Rubin convergence diagnostic; ESS = effective sample size; LOOIC = leave-one-out information criterion; SE = standard error; $\\Delta$LOOIC reported vs the lowest-LOOIC model; $p_{\\text{eff}}$ = effective number of parameters; $n_{\\text{hi-k}}$ = count of observations with Pareto-$k > 0.7$; GP = Gaussian process; $\\delta^2$H$_p$ = $\\delta^2$H$_{\\text{precip}}$; elev = elevation.",
  "\\end{minipage}",
  "\\end{table}",
  "\\end{landscape}"
)

writeLines(t1_lines, table1_path)
cat("  ->", table1_path, "\n")

# --------- WRITE TABLE 2 BODY FRAGMENT -------------------------------------

cat("Writing Table 2 body fragment...\n")

table2_path <- .outpath("table2_global_params_body.tex")
t2_lines <- c(
  "% Auto-generated by scripts/regen_tables.R.",
  "% Do not hand-edit — regenerate via `Rscript scripts/regen_tables.R`.",
  ""
)

for (m in ALL_MODELS) {
  row <- table2_rows |> filter(model == m)
  is_sp <- grepl("_sp$", m)
  rmse  <- sprintf("%.1f", row$rmse_point)
  r2    <- sprintf("%.3f", row$r2_point)
  b0    <- fmt_meanci(row$beta_0_mean,row$beta_0_lo,row$beta_0_hi,digits = 1)
  bo    <- fmt_meanci(row$beta_oipc_mean, row$beta_oipc_lo, row$beta_oipc_hi, digits = 3)
  lam   <- fmt_meanci(row$lambda_mean,row$lambda_lo,row$lambda_hi,digits = 1)
  if (is_sp) {
    gp <- sprintf("%.0f [%.0f, %.0f]",
                  row$gp_scale_mean, row$gp_scale_lo, row$gp_scale_hi)
    ks <- sprintf("%.3f", row$knot_slope_sd_mean)
    ki <- sprintf("%.1f", row$knot_int_sd_mean_permil)
  } else {
    gp <- "-"; ks <- "-"; ki <- "-"
  }
  t2_lines <- c(t2_lines, paste0(
    paste(c(MODEL_LABELS[m], rmse, r2, b0, bo, lam, gp, ks, ki),
          collapse = " & "),
    "\\\\"
  ))
}

writeLines(t2_lines, table2_path)
cat("  ->", table2_path, "\n")

# --------- WRITE TABLE 3 ---------------------------------------------------

cat("Writing Table 3...\n")

table3_path <- .outpath("table3_variance_decomp.tex")
t3_lines <- c(
  "% Auto-generated by scripts/regen_tables.R.",
  "% Do not hand-edit — regenerate via `Rscript scripts/regen_tables.R`.",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Marginal spatial-component variance for spatial models}",
  "\\label{tab:variance_decomp}",
  "\\small",
  "",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Model & Spatial intercept (\\%) & Spatial slope (\\%)\\\\",
  "\\midrule"
)

t3_models_ordered <- c(
  "baseline_sp", "baseline_env_sp", "baseline_veg_sp",
  "full_sp", "full_interact_sp",
  "elevation_only_sp", "elevation_c4_sp",
  "c4_only_sp", "elevation_c4_interact_sp"
)

for (i in seq_along(t3_models_ordered)) {
  m <- t3_models_ordered[i]
  row <- table3_rows |> filter(model == m)
  if (nrow(row) == 0) next
  t3_lines <- c(t3_lines, sprintf(
    "%s & %.1f & %.1f\\\\",
    MODEL_LABELS[m], row$marginal_intercept_share_pct,
    row$marginal_slope_share_pct
  ))
  if (i == 5 && i < length(t3_models_ordered)) {
    t3_lines <- c(t3_lines, "\\addlinespace")
  }
}

t3_lines <- c(t3_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{2mm}",
  "\\footnotesize",
  paste0(
    "\\textit{Note:} For each posterior draw, the spatial-intercept contribution ",
    "and the local slope deviation multiplied by that draw's fitted, scale-weighted ",
    "OIPC predictor were evaluated across the 1,128 calibration sites. The table ",
    "reports each marginal variance as a percentage of their sum. Covariance between ",
    "the two contributions is excluded, so these percentages compare component ",
    "magnitudes and are not an additive partition of total $\\delta^2$H$_{wax}$ variance."
  ),
  "\\end{minipage}",
  "\\end{table}"
)

writeLines(t3_lines, table3_path)
cat("  ->", table3_path, "\n")

# --------- WRITE TABLE 4 ---------------------------------------------------

cat("Writing Table 4...\n")

fmt_veg <- function(mean_val, lo, hi) {
  if (is.na(mean_val)) return("-")
  sprintf("%.3f [%.3f, %.3f]", mean_val, lo, hi)
}

table4_path <- .outpath("table4_environmental.tex")
t4_lines <- c(
  "% Auto-generated by scripts/regen_tables.R.",
  "% Do not hand-edit — regenerate via `Rscript scripts/regen_tables.R`.",
  "",
  "\\begin{landscape}",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Environmental covariate coefficients}",
  "\\label{tab:environmental-coefficients}",
  "\\scriptsize",
  "",
  "\\begin{tabular}{lccccc}",
  "\\toprule",
  "Model & $\\beta_{\\text{C4}}$ & $\\beta_{\\text{tree}}$ & $\\beta_{\\text{shrub}}$ & $\\beta_{\\text{grass}}$ & $\\beta_{\\text{precip}}$\\\\",
  "\\midrule"
)

for (i in seq_along(ALL_MODELS)) {
  m <- ALL_MODELS[i]
  row <- table4_rows |> filter(model == m)
  c4    <- fmt_veg(row$C4_mean,    row$C4_lo,    row$C4_hi)
  tr    <- fmt_veg(row$tree_mean,  row$tree_lo,  row$tree_hi)
  sh    <- fmt_veg(row$shrub_mean, row$shrub_lo, row$shrub_hi)
  gr    <- fmt_veg(row$grass_mean, row$grass_lo, row$grass_hi)
  pr    <- fmt_veg(row$precip_mean,row$precip_lo,row$precip_hi)
  t4_lines <- c(t4_lines, sprintf("%s & %s & %s & %s & %s & %s\\\\",
                                  MODEL_LABELS[m], c4, tr, sh, gr, pr))
  if (i %in% c(5, 10) && i < length(ALL_MODELS)) {
    t4_lines <- c(t4_lines, "\\addlinespace")
  }
}

t4_lines <- c(t4_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{2mm}",
  "\\footnotesize",
  "\\textit{Note:} Coefficient estimates shown as posterior mean [95\\% credible interval]. $\\beta_{\\text{C4}}$ = C4 grass fraction effect; $\\beta_{\\text{tree}}$, $\\beta_{\\text{shrub}}$, $\\beta_{\\text{grass}}$ = plant functional type effects; $\\beta_{\\text{precip}}$ = precipitation effect. Dashes indicate parameters not included in the model. All coefficients on the Stan model's standardized $\\delta^2$H scale.",
  "\\end{minipage}",
  "\\end{table}",
  "\\end{landscape}"
)

writeLines(t4_lines, table4_path)
cat("  ->", table4_path, "\n")

# Table 5: delta2H_precip x vegetation interaction coefficients ------------
# Only interaction-bearing variants are listed (rows where any interaction
# term is present). Portrait table; 4 interaction columns.

table5_path <- .outpath("table5_interactions.tex")
# The documented interaction-model variants (Section S2.3.6) that actually
# estimate the delta2H_precip x vegetation terms. `full`/`full_sp` carry these
# parameters fixed at zero and the elevation_c4 variants do not store them, so
# only the two `*_interact*` variants with estimated interactions are tabulated.
interact_models <- c("full_interact", "full_interact_sp")

t5_lines <- c(
  "% Auto-generated by scripts/regen_tables.R.",
  "% Do not hand-edit — regenerate via `Rscript scripts/regen_tables.R`.",
  "",
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{$\\delta^{2}$H$_{\\text{precip}}\\times$vegetation interaction coefficients}",
  "\\label{tab:interactions}",
  "\\small",
  "",
  "\\begin{tabular}{lcccc}",
  "\\toprule",
  "Model & $\\beta_{\\delta^{2}\\mathrm{H}_{p}\\times\\text{C4}}$ & $\\beta_{\\delta^{2}\\mathrm{H}_{p}\\times\\text{tree}}$ & $\\beta_{\\delta^{2}\\mathrm{H}_{p}\\times\\text{shrub}}$ & $\\beta_{\\delta^{2}\\mathrm{H}_{p}\\times\\text{grass}}$\\\\",
  "\\midrule"
)

for (m in interact_models) {
  row <- table5_rows |> filter(model == m)
  c4 <- fmt_veg(row$C4_mean,    row$C4_lo,    row$C4_hi)
  tr <- fmt_veg(row$tree_mean,  row$tree_lo,  row$tree_hi)
  sh <- fmt_veg(row$shrub_mean, row$shrub_lo, row$shrub_hi)
  gr <- fmt_veg(row$grass_mean, row$grass_lo, row$grass_hi)
  t5_lines <- c(t5_lines, sprintf("%s & %s & %s & %s & %s\\\\",
                                  MODEL_LABELS[m], c4, tr, sh, gr))
}

t5_lines <- c(t5_lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\begin{minipage}{\\linewidth}",
  "\\vspace{2mm}",
  "\\footnotesize",
  "\\textit{Note:} Each coefficient is the change in the $\\delta^{2}$H\\textsubscript{wax}--$\\delta^{2}$H\\textsubscript{precip} slope per unit increase in that vegetation fraction; a positive value indicates a steeper slope in that vegetation type. Estimates shown as posterior mean [95\\% credible interval]. Dashes indicate interaction terms not included in the model. All coefficients on the Stan model's standardized $\\delta^2$H scale.",
  "\\end{minipage}",
  "\\end{table}"
)

writeLines(t5_lines, table5_path)
cat("  ->", table5_path, "\n")

# --------- WRITE COMPANION CSVs ---------------------------------------------
# Each CSV contains the same values as its generated .tex table, providing a
# machine-readable record for independent numeric checks.

cat("Writing companion CSVs to", .output_dir, "...\n")

# Table 1 CSV
t1_csv <- table1_df |>
  transmute(
    model = model,
    predictors = PREDICTORS[model],
    max_rhat = round(as.numeric(max_rhat), 4),
    min_ess  = as.integer(min_ess),
    looic    = round(looic, 1),
    se_looic = round(se_looic, 1),
    delta_looic = round(delta_looic, 1),
    se_delta = round(se_delta, 1),
    p_eff    = round(p_eff, 1),
    n_hi_k   = as.integer(n_hi_k)
  )
write.csv(t1_csv, .outpath("table1_model_performance.csv"),
          row.names = FALSE)

# Table 2 CSV
t2_csv <- table2_rows |>
  transmute(
    model = model,
    rmse_permil_point = rmse_point,
    rmse_permil_mean = rmse_mean,
    rmse_permil_q025 = rmse_lo,
    rmse_permil_q975 = rmse_hi,
    r_squared_point = r2_point,
    r_squared_mean = r2_mean,
    r_squared_q025 = r2_lo,
    r_squared_q975 = r2_hi,
    beta_0_permil_mean = beta_0_mean,
    beta_0_permil_q025 = beta_0_lo,
    beta_0_permil_q975 = beta_0_hi,
    beta_oipc_mean = beta_oipc_mean,
    beta_oipc_q025 = beta_oipc_lo,
    beta_oipc_q975 = beta_oipc_hi,
    lambda_int_km_mean = lambda_mean,
    lambda_int_km_q025 = lambda_lo,
    lambda_int_km_q975 = lambda_hi,
    gp_scale_km_mean = gp_scale_mean,
    gp_scale_km_q025 = gp_scale_lo,
    gp_scale_km_q975 = gp_scale_hi,
    knot_slope_sd_mean = knot_slope_sd_mean,
    knot_intercept_sd_permil_mean = knot_int_sd_mean_permil
  )
write.csv(t2_csv, .outpath("table2_global_params.csv"),
          row.names = FALSE)

# Table 3 CSV. Retain the covariance and realized spatial variance for audit;
# Figure 2 uses only the explicitly named marginal component shares.
t3_csv <- table3_rows
write.csv(t3_csv, .outpath("table3_variance_decomp.csv"),
          row.names = FALSE)

# Table 4 CSV
t4_csv <- table4_rows |>
  transmute(
    model = model,
    beta_c4_mean = C4_mean,     beta_c4_q05 = C4_lo,     beta_c4_q95 = C4_hi,
    beta_tree_mean = tree_mean, beta_tree_q05 = tree_lo, beta_tree_q95 = tree_hi,
    beta_shrub_mean = shrub_mean, beta_shrub_q05 = shrub_lo, beta_shrub_q95 = shrub_hi,
    beta_grass_mean = grass_mean, beta_grass_q05 = grass_lo, beta_grass_q95 = grass_hi,
    beta_precip_mean = precip_mean, beta_precip_q05 = precip_lo, beta_precip_q95 = precip_hi
  )
write.csv(t4_csv, .outpath("table4_environmental.csv"),
          row.names = FALSE)

# Table 5 CSV
t5_csv <- table5_rows |>
  transmute(
    model = model,
    beta_c4_mean = C4_mean,       beta_c4_q025 = C4_lo,       beta_c4_q975 = C4_hi,
    beta_tree_mean = tree_mean,   beta_tree_q025 = tree_lo,   beta_tree_q975 = tree_hi,
    beta_shrub_mean = shrub_mean, beta_shrub_q025 = shrub_lo, beta_shrub_q975 = shrub_hi,
    beta_grass_mean = grass_mean, beta_grass_q025 = grass_lo, beta_grass_q975 = grass_hi
  )
write.csv(t5_csv, .outpath("table5_interactions.csv"),
          row.names = FALSE)

cat("\nAll five tables regenerated from", MODEL_RUN_DIR, "\n")
cat("CSV companions written for independent numeric checks.\n")
