# Comprehensive analysis of full models for Section 3.4
library(tidyverse)
library(posterior)
library(loo)

source("scripts/posterior_helpers.R")
source("manuscript/table_code/table_helpers.R")

# Create output directory
output_dir <- "model_analysis/tables"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("FULL MODEL ANALYSIS FOR SECTION 3.4\n")
cat("====================================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# 1. MODEL PERFORMANCE COMPARISON
#───────────────────────────────────────────────────────────────────────────────

cat("1. MODEL PERFORMANCE METRICS\n")
cat(strrep("-", 60), "\n\n")

# Models to compare — full 14-model roster per manuscript/TABLES.md (Table 1).
# Order matches config.yaml for stable table row ordering.
all_models <- c(
  "baseline", "baseline_sp",
  "baseline_env", "baseline_env_sp",
  "baseline_veg", "baseline_veg_sp",
  "full", "full_sp",
  "full_interact", "full_interact_sp",
  "elevation_only_sp", "elevation_c4_sp",
  "c4_only_sp", "elevation_c4_interact_sp"
)

# Full-model subset used for per-coefficient extraction (section 2).
# Section 1 runs over all 14 models; section 2 only over the four "full" ones.
full_models <- c("full", "full_sp", "full_interact", "full_interact_sp")

performance_metrics <- list()
# Carry loo + summary bundles so Table 1 producer (section 8) can compute
# ΔLOOIC via loo_compare() and pull max_rhat / min_ess_bulk from diagnostics
# without re-reading the rds bundles.
loo_cache         <- list()
diagnostics_cache <- list()

for (model_name in all_models) {
  # Helpers raise a clear error if any piece is missing; wrap in tryCatch
  # so a missing model doesn't abort the whole loop.
  bundle <- tryCatch({
    list(draws     = load_draws(model_name),
         loo       = load_loo(model_name),
         stan_data = load_stan_data(model_name))
  }, error = function(e) {
    cat(sprintf("%-20s: %s\n", model_name, conditionMessage(e)))
    NULL
  })
  if (is.null(bundle)) next

  loo_result <- bundle$loo
  stan_data  <- bundle$stan_data
  loo_cache[[model_name]]         <- loo_result
  diagnostics_cache[[model_name]] <- tryCatch(
    readRDS(file.path(APRIL_RUN, model_name, "diagnostics.rds")),
    error = function(e) NULL
  )

  # Posterior predictive means: y_rep is a draws matrix with one column per
  # observation; colMeans gives the per-obs posterior mean.
  y_rep       <- as_draws_matrix(subset_draws(bundle$draws, variable = "d2H_rep"))
  y_pred_mean <- colMeans(y_rep)
  y_obs       <- stan_data$d2H_wax

  # Calculate R² and RMSE in standardized units
  ss_res <- sum((y_obs - y_pred_mean)^2)
  ss_tot <- sum((y_obs - mean(y_obs))^2)
  r_squared <- 1 - (ss_res / ss_tot)
  rmse <- sqrt(mean((y_obs - y_pred_mean)^2))

  # Convert RMSE to original scale
  rmse_permil <- rmse * stan_data$scaling_params$d2H_sd

  performance_metrics[[model_name]] <- data.frame(
    model = model_name,
    looic = loo_result$estimates["looic", "Estimate"],
    elpd = loo_result$estimates["elpd_loo", "Estimate"],
    p_loo = loo_result$estimates["p_loo", "Estimate"],
    r_squared = r_squared,
    rmse_std = rmse,
    rmse_permil = rmse_permil,
    n = length(y_obs),
    stringsAsFactors = FALSE
  )

  cat(sprintf("%-20s: LOOIC = %8.1f, R² = %.3f, RMSE = %.1f‰\n",
              model_name, performance_metrics[[model_name]]$looic,
              r_squared, rmse_permil))
}

performance_df <- bind_rows(performance_metrics)

# Sort by LOOIC
performance_df <- performance_df %>% arrange(looic)

cat("\nRanked by LOOIC (lower is better):\n")
print(performance_df %>% select(model, looic, r_squared, rmse_permil))

#───────────────────────────────────────────────────────────────────────────────
# 2. EXTRACT COEFFICIENTS FOR FULL MODELS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n2. COEFFICIENT EXTRACTION FOR FULL MODELS\n")
cat(strrep("-", 60), "\n")

coefficient_list <- list()

for (model_name in full_models) {
  cat("\n", model_name, ":\n")

  draws <- tryCatch(load_draws(model_name),
                    error = function(e) { cat("  skip:", conditionMessage(e), "\n"); NULL })
  if (is.null(draws)) next
  config <- load_config(model_name)
  vars_present <- variables(draws)

  # Parameters to extract. Names must match 4d_leaf_wax_spatial_model.stan exactly.
  # Stan exports `beta_precip` (not `beta_precip_amount`) and `beta_oipc_x_*`
  # (not `beta_oipc_*`). There are no `beta_precip_*` interaction terms in the
  # model; those extraction requests were stale and have been removed.
  params_to_extract <- c("beta_oipc", "beta_c4", "beta_tree", "beta_shrub",
                        "beta_grass", "beta_precip")

  # Add OIPC interaction terms if applicable (no precip interactions in Stan)
  if (grepl("interact", model_name)) {
    params_to_extract <- c(params_to_extract,
                          "beta_oipc_x_c4",
                          "beta_oipc_x_tree", "beta_oipc_x_shrub", "beta_oipc_x_grass")
  }

  for (param in params_to_extract) {
    if (!(param %in% vars_present)) next
    d <- as.numeric(as_draws_matrix(subset_draws(draws, variable = param)))

    coef_stats <- data.frame(
      model = model_name,
      parameter = param,
      mean = mean(d),
      median = median(d),
      sd = sd(d),
      lower_95 = quantile(d, 0.025),
      upper_95 = quantile(d, 0.975),
      lower_50 = quantile(d, 0.25),
      upper_50 = quantile(d, 0.75),
      significant = (quantile(d, 0.025) > 0) | (quantile(d, 0.975) < 0),
      stringsAsFactors = FALSE
    )

    coefficient_list[[paste(model_name, param, sep = "_")]] <- coef_stats

    sig_marker <- if(coef_stats$significant) "*" else " "
    cat(sprintf("  %-20s: %7.3f [%7.3f, %7.3f]%s\n",
                param, coef_stats$mean, coef_stats$lower_95,
                coef_stats$upper_95, sig_marker))
  }

  # Elevation enters via B-spline coefficients `beta_elev_bspline` (a vector),
  # not a scalar `beta_elevation`. Summarize mean magnitude across coefficients.
  # beta_elev_bspline is NOT in draws_to_save (Phase 5 W1 widened list); pull
  # its summary from diagnostics.rds where all parameters are summarized.
  summ <- load_summaries(model_name)
  elev_rows <- summ[grepl("^beta_elev_bspline", summ$variable), ]
  if (nrow(elev_rows) > 0) {
    # Use posterior-mean coefficient values; across-coefficient stats only.
    elev_mean_abs <- mean(abs(elev_rows$mean))
    elev_sd       <- sd(elev_rows$mean)
    cat(sprintf("  Elevation (|β| avg) : %7.3f (SD = %.3f)\n", elev_mean_abs, elev_sd))
  }
}

coefficients_df <- bind_rows(coefficient_list)

#───────────────────────────────────────────────────────────────────────────────
# 3. VARIANCE DECOMPOSITION FOR SPATIAL MODELS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n3. VARIANCE DECOMPOSITION ANALYSIS\n")
cat(strrep("-", 60), "\n")

# All 9 spatial models per manuscript/TABLES.md (Table 3 variance decomposition).
spatial_models <- c(
  "baseline_sp",
  "baseline_env_sp",
  "baseline_veg_sp",
  "full_sp",
  "full_interact_sp",
  "elevation_only_sp",
  "elevation_c4_sp",
  "c4_only_sp",
  "elevation_c4_interact_sp"
)

variance_components <- list()

for (model_name in spatial_models) {
  summ <- tryCatch(load_summaries(model_name),
                   error = function(e) { cat(sprintf("%-20s skip: %s\n", model_name, conditionMessage(e))); NULL })
  if (is.null(summ)) next

  # Use Stan's own variance decomposition from generated quantities.
  # 4d_leaf_wax_spatial_model.stan:504-511 defines var_spatial_intercept,
  # var_spatial_slope, var_residual on a consistent scale, and
  # prop_variance_spatial / prop_variance_residual as the final proportions.
  # The prior ad-hoc calculation mixed sigma_intercept_spatial (original ‰)
  # with sigma (standardized) and inflated the spatial share to 100%.
  get_mean <- function(var) {
    row <- summ[summ$variable == var, , drop = FALSE]
    if (nrow(row) != 1) return(NA_real_)
    row$mean
  }

  var_int      <- get_mean("var_spatial_intercept")
  var_slope    <- get_mean("var_spatial_slope")
  var_spatial  <- get_mean("var_total_spatial")
  var_resid    <- get_mean("var_residual")
  var_total    <- get_mean("var_total")
  prop_spatial <- get_mean("prop_variance_spatial")
  prop_resid   <- get_mean("prop_variance_residual")

  # Additional reporting: posterior-mean residual SD in original ‰.
  sigma_resid_orig <- get_mean("sigma_residual_original")

  # Intercept / slope shares expressed two ways:
  #   * as % of total variance (matches what the prior script reported)
  #   * as % of the spatial component only (what Table 3 displays, so the
  #     intercept and slope columns sum to 100 within each spatial model)
  .pct <- function(num, denom) if (is.finite(denom) && denom > 0) 100 * num / denom else NA_real_
  variance_components[[model_name]] <- data.frame(
    model = model_name,
    var_spatial_intercept = var_int,
    var_spatial_slope     = var_slope,
    var_total_spatial     = var_spatial,
    var_residual          = var_resid,
    var_total             = var_total,
    var_intercept_pct     = .pct(var_int,   var_total),
    var_slope_pct         = .pct(var_slope, var_total),
    var_intercept_pct_of_spatial = .pct(var_int,   var_spatial),
    var_slope_pct_of_spatial     = .pct(var_slope, var_spatial),
    var_spatial_total_pct = 100 * prop_spatial,
    var_obs_pct           = 100 * prop_resid,
    sigma_residual_permil = sigma_resid_orig,
    stringsAsFactors = FALSE
  )

  cat(sprintf("\n%s:\n", model_name))
  cat(sprintf("  Spatial variance (total): %.1f%%\n",
              variance_components[[model_name]]$var_spatial_total_pct))
  cat(sprintf("    - Intercept: %.1f%%\n",
              variance_components[[model_name]]$var_intercept_pct))
  cat(sprintf("    - Slope:     %.1f%%\n",
              variance_components[[model_name]]$var_slope_pct))
  cat(sprintf("  Residual variance:        %.1f%% (σ = %.2f ‰)\n",
              variance_components[[model_name]]$var_obs_pct, sigma_resid_orig))
}

variance_df <- bind_rows(variance_components)

# Calculate variance reduction
if ("baseline_sp" %in% names(variance_components)) {
  baseline_spatial_var <- variance_components[["baseline_sp"]]$var_spatial_total_pct

  cat("\nVariance reduction from baseline_sp:\n")
  for (model in c("full_sp", "full_interact_sp")) {
    if (model %in% names(variance_components)) {
      model_spatial_var <- variance_components[[model]]$var_spatial_total_pct
      reduction <- baseline_spatial_var - model_spatial_var
      cat(sprintf("  %s: %.1f%% reduction\n", model, reduction))
    }
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 4. PARAMETER CORRELATION ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n4. PREDICTOR CORRELATION ANALYSIS\n")
cat(strrep("-", 60), "\n")

# Build a per-obs predictor matrix for multicollinearity analysis.
# Environmental covariates in stan_data are 1124×9 matrices (obs × spatial
# scale); use the narrowest-scale column (index 1) as the representative
# point value. Scalar per-obs variables come from the sediment frame.
sd_full <- load_stan_data("full")
sed     <- load_sediment()

.col1 <- function(m) if (is.matrix(m)) m[, 1] else as.numeric(m)

predictors <- data.frame(
  oipc      = .col1(sd_full$oipc_values),
  elevation = .col1(sd_full$elevation_values)
)
# c4 in stan_data: try weighted first, fall back to sediment scalar
if (!is.null(sd_full$c4_values) && is.matrix(sd_full$c4_values)) {
  predictors$c4 <- .col1(sd_full$c4_values)
} else {
  predictors$c4 <- sed$c4_mean_filled
}
if (!is.null(sd_full$precip_values) && is.matrix(sd_full$precip_values)) {
  predictors$precip <- .col1(sd_full$precip_values)
} else if (!is.null(sed$annual_precip)) {
  predictors$precip <- sed$annual_precip
}
if (!is.null(sd_full$pft_tree) && is.matrix(sd_full$pft_tree)) {
  predictors$tree  <- .col1(sd_full$pft_tree)
  predictors$shrub <- .col1(sd_full$pft_shrub)
  predictors$grass <- .col1(sd_full$pft_grass)
}

# Calculate correlation matrix
cor_matrix <- cor(predictors, use = "pairwise.complete.obs")

cat("\nCorrelation matrix:\n")
print(round(cor_matrix, 3))

# Find high correlations
high_cor <- which(abs(cor_matrix) > 0.6 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor) > 0) {
  cat("\nHigh correlations (|r| > 0.6):\n")
  for (i in 1:nrow(high_cor)) {
    if (high_cor[i, 1] < high_cor[i, 2]) {  # Avoid duplicates
      cat(sprintf("  %s vs %s: r = %.3f\n",
                  rownames(cor_matrix)[high_cor[i, 1]],
                  colnames(cor_matrix)[high_cor[i, 2]],
                  cor_matrix[high_cor[i, 1], high_cor[i, 2]]))
    }
  }
} else {
  cat("\nNo high correlations (|r| > 0.6) found between predictors\n")
}

#───────────────────────────────────────────────────────────────────────────────
# 5. COEFFICIENT COMPARISON ACROSS MODEL COMPLEXITY
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n5. COEFFICIENT CHANGES WITH MODEL COMPLEXITY\n")
cat(strrep("-", 60), "\n")

# Compare beta_oipc across models
cat("\nbeta_oipc (OIPC slope) changes:\n")
oipc_coefs <- coefficients_df %>%
  filter(parameter == "beta_oipc") %>%
  select(model, mean, lower_95, upper_95)

if (nrow(oipc_coefs) > 0) {
  print(oipc_coefs)
}

# Compare vegetation effects
cat("\nVegetation effect changes:\n")
veg_params <- c("beta_c4", "beta_tree", "beta_shrub", "beta_grass")

for (param in veg_params) {
  param_data <- coefficients_df %>% filter(parameter == param)
  if (nrow(param_data) > 0) {
    cat(sprintf("\n%s:\n", param))
    for (i in 1:nrow(param_data)) {
      sig_marker <- if(param_data$significant[i]) "*" else " "
      cat(sprintf("  %-20s: %7.3f [%7.3f, %7.3f]%s\n",
                  param_data$model[i], param_data$mean[i],
                  param_data$lower_95[i], param_data$upper_95[i], sig_marker))
    }
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 6. RELATIVE IMPORTANCE
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n6. RELATIVE IMPORTANCE OF PREDICTORS\n")
cat(strrep("-", 60), "\n")

for (model_name in c("full", "full_sp")) {
  model_coefs <- coefficients_df %>%
    filter(model == model_name, !grepl("interact|x", parameter))

  if (nrow(model_coefs) > 0) {
    # Sort by absolute effect size
    model_coefs <- model_coefs %>%
      mutate(abs_mean = abs(mean)) %>%
      arrange(desc(abs_mean))

    cat(sprintf("\n%s - Ranked by absolute effect size:\n", model_name))
    for (i in 1:min(nrow(model_coefs), 6)) {
      sig_marker <- if(model_coefs$significant[i]) "*" else " "
      cat(sprintf("  %2d. %-20s: |β| = %.3f%s\n",
                  i, model_coefs$parameter[i], model_coefs$abs_mean[i], sig_marker))
    }
  }
}

# Check which predictors lose significance
cat("\n\nPredictors that lose significance in full models:\n")
baseline_sig <- coefficients_df %>%
  filter(model %in% c("baseline", "baseline_sp"), significant) %>%
  pull(parameter) %>% unique()

full_nonsig <- coefficients_df %>%
  filter(model %in% c("full", "full_sp"), !significant) %>%
  pull(parameter) %>% unique()

lost_sig <- intersect(baseline_sig, full_nonsig)
if (length(lost_sig) > 0) {
  cat("  ", paste(lost_sig, collapse = ", "), "\n")
} else {
  cat("  None\n")
}

#───────────────────────────────────────────────────────────────────────────────
# 7. SAVE OUTPUTS
#───────────────────────────────────────────────────────────────────────────────

# Save coefficient comparison table
write.csv(coefficients_df,
          file.path(output_dir, "full_model_coefficients.csv"),
          row.names = FALSE)

# Save variance decomposition
write.csv(variance_df,
          file.path(output_dir, "full_model_variance.csv"),
          row.names = FALSE)

# Save performance metrics
write.csv(performance_df,
          file.path(output_dir, "full_model_performance.csv"),
          row.names = FALSE)

# Save correlation matrix
write.csv(cor_matrix,
          file.path(output_dir, "predictor_correlations.csv"))

cat("\n\nAll outputs saved to:", output_dir, "\n")

#───────────────────────────────────────────────────────────────────────────────
# 8. TABLE 1 PRODUCER — manuscript/tables/table1_model_performance.tex
#───────────────────────────────────────────────────────────────────────────────
#
# Columns (matches the hand-authored .tex pre-W6): Model, Predictors, Max Rhat,
# Min ESS, LOOIC, SE, ΔLOOIC, SE(Δ), p_eff, n_hi_k. Numbers come from:
#   loo.rds           — LOOIC, SE, p_loo (= p_eff), pareto_k (n_hi_k = sum k>0.7)
#   diagnostics.rds   — max_rhat, min_ess_bulk
#   loo_compare()     — ΔLOOIC + SE(Δ), measured against the lowest LOOIC in
#                       the set (elpd_diff / se_diff scaled by -2).
#
# `Predictors` is a hand-maintained latex string per model (the model names
# themselves encode structure but the manuscript needs a reader-facing
# predictor list). Model-name latex escaping lives in .latex_model_name().

cat("\n\n8. TABLE 1 PRODUCER\n")
cat(strrep("-", 60), "\n")

predictors_latex <- list(
  baseline                 = "$\\delta^2$H$_p$",
  baseline_sp              = "$\\delta^2$H$_p$ + GP",
  baseline_env             = "$\\delta^2$H$_p$, elev, precip",
  baseline_env_sp          = "$\\delta^2$H$_p$, elev, precip + GP",
  baseline_veg             = "$\\delta^2$H$_p$, PFT, C4",
  baseline_veg_sp          = "$\\delta^2$H$_p$, PFT, C4 + GP",
  full                     = "$\\delta^2$H$_p$, PFT, C4, elev, precip",
  full_sp                  = "$\\delta^2$H$_p$, PFT, C4, elev, precip + GP",
  full_interact            = "$\\delta^2$H$_p$, PFT $\\times$ $\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$, elev, precip",
  full_interact_sp         = "$\\delta^2$H$_p$, PFT $\\times$ $\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$, elev, precip + GP",
  elevation_only_sp        = "$\\delta^2$H$_p$, elev + GP",
  elevation_c4_sp          = "$\\delta^2$H$_p$, elev, C4 + GP",
  c4_only_sp               = "$\\delta^2$H$_p$, C4 + GP",
  elevation_c4_interact_sp = "$\\delta^2$H$_p$, C4 $\\times$ $\\delta^2$H$_p$, elev + GP"
)

.latex_model_name <- function(name) {
  # Mirror the hand-authored table style: trailing `_sp` renders as a math
  # subscript so "baseline_sp" appears as baseline$_{\text{sp}}$; any other
  # underscore is escaped as \_. fixed = TRUE keeps gsub from treating the
  # replacement as a regex backreference.
  if (grepl("_sp$", name)) {
    base <- sub("_sp$", "", name)
    base <- gsub("_", "\\_", base, fixed = TRUE)
    paste0(base, "$\\_{\\text{sp}}$")
  } else {
    gsub("_", "\\_", name, fixed = TRUE)
  }
}

# loo_compare on all cached psis_loo objects; rows come back sorted by
# best → worst with elpd_diff = 0 for the best. Names preserved.
if (length(loo_cache) >= 2) {
  comp <- loo_compare(loo_cache)
  # Rownames are "model1", "model2", ... in loo>=2.4; older loo uses the
  # names directly. Defensive: if rownames look like "model<i>", remap.
  if (all(grepl("^model[0-9]+$", rownames(comp)))) {
    rownames(comp) <- names(loo_cache)[as.integer(sub("model", "", rownames(comp)))]
  }
  comp_df <- as.data.frame(comp)
  comp_df$model <- rownames(comp)
} else {
  comp_df <- data.frame(model = names(loo_cache),
                        elpd_diff = 0, se_diff = 0,
                        stringsAsFactors = FALSE)
}

.hi_k_count <- function(loo_result) {
  # psis_loo puts pareto_k under $diagnostics$pareto_k (loo >= 2.0).
  pk <- loo_result$diagnostics$pareto_k
  if (is.null(pk)) return(NA_integer_)
  sum(pk > 0.7)
}

table1_rows <- lapply(all_models, function(m) {
  if (!(m %in% names(loo_cache))) return(NULL)
  lo   <- loo_cache[[m]]
  diag <- diagnostics_cache[[m]]
  delta_row <- comp_df[comp_df$model == m, , drop = FALSE]
  # ΔLOOIC = -2 * elpd_diff; SE(Δ) = 2 * se_diff. Best model has 0/0.
  delta_looic <- if (nrow(delta_row) == 1) -2 * delta_row$elpd_diff else NA_real_
  se_delta    <- if (nrow(delta_row) == 1)  2 * delta_row$se_diff   else NA_real_

  data.frame(
    model       = .latex_model_name(m),
    predictors  = if (!is.null(predictors_latex[[m]])) predictors_latex[[m]] else "—",
    max_rhat    = if (!is.null(diag$max_rhat))     round(diag$max_rhat, 3)      else NA_real_,
    min_ess     = if (!is.null(diag$min_ess_bulk)) round(diag$min_ess_bulk, 0)  else NA_real_,
    looic       = round(lo$estimates["looic", "Estimate"], 1),
    se_looic    = round(lo$estimates["looic", "SE"], 1),
    delta_looic = round(delta_looic, 1),
    se_delta    = round(se_delta, 1),
    p_eff       = round(lo$estimates["p_loo", "Estimate"], 1),
    n_hi_k      = .hi_k_count(lo),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
table1_df <- bind_rows(table1_rows)

# CSV sibling for the numeric audit (W6e). Human-readable row per model.
write.csv(table1_df,
          file.path(output_dir, "table1_model_performance.csv"),
          row.names = FALSE)

col_names <- c("Model", "Predictors", "Max $\\hat{R}$", "Min ESS",
               "LOOIC", "SE", "$\\Delta$LOOIC", "SE",
               "$p_{\\text{eff}}$", "$n_{\\text{hi-k}}$")

note_text <- paste(
  "$\\hat{R}$ = Gelman-Rubin convergence diagnostic; ESS = effective sample size;",
  "LOOIC = leave-one-out information criterion; SE = standard error;",
  "$\\Delta$LOOIC reported vs the lowest-LOOIC model;",
  "$p_{\\text{eff}}$ = effective number of parameters;",
  "$n_{\\text{hi-k}}$ = count of observations with Pareto-$k > 0.7$;",
  "GP = Gaussian process; $\\delta^2$H$_p$ = $\\delta^2$H$_{\\text{precip}}$; elev = elevation."
)

emit_standalone_tex(
  table1_df,
  path      = "manuscript/tables/table1_model_performance.tex",
  caption   = "Model performance metrics for all candidate models",
  label     = "tab:model_performance",
  col.names = col_names,
  align     = c("l", "p{4.5cm}", rep("c", 8)),
  note      = note_text,
  landscape = TRUE,
  size_macro = "footnotesize",
  source_script = "extract_full_model_analysis.R"
)
cat("wrote manuscript/tables/table1_model_performance.tex\n")

#───────────────────────────────────────────────────────────────────────────────
# 9. TABLE 3 PRODUCER — manuscript/tables/table3_variance_decomp.tex
#───────────────────────────────────────────────────────────────────────────────
#
# Columns: Model, Spatial (%), Residual (%), Intercept (%), Slope (%).
# All 9 spatial models, rows ordered as defined in `spatial_models` above
# (section 3). The intercept / slope columns are percentages of the spatial
# component (so they sum to 100 within a row); the Spatial / Residual
# columns are percentages of total variance (so those two sum to 100).
#
# Source of every number: generated-quantity posterior means pulled from
# diagnostics.rds by section 3. The CLAUDE.md integrity rule applies — no
# re-sampling, no ad-hoc back-transforms. The Phase 5 W5 bug-fix replaced
# the sigma_int² + sigma² mixing with Stan's prop_variance_* outputs; this
# table now carries those values directly.

cat("\n\n9. TABLE 3 PRODUCER\n")
cat(strrep("-", 60), "\n")

# Preserve the `spatial_models` order defined in section 3.
table3_df <- variance_df %>%
  mutate(model_order = match(model, spatial_models)) %>%
  arrange(model_order) %>%
  transmute(
    Model      = sapply(model, .latex_model_name, USE.NAMES = FALSE),
    `Spatial (\\%)`   = sprintf("%.1f", var_spatial_total_pct),
    `Residual (\\%)`  = sprintf("%.1f", var_obs_pct),
    `Intercept (\\%)` = sprintf("%.1f", var_intercept_pct_of_spatial),
    `Slope (\\%)`     = sprintf("%.1f", var_slope_pct_of_spatial)
  )

write.csv(variance_df %>%
            mutate(model_order = match(model, spatial_models)) %>%
            arrange(model_order) %>%
            select(-model_order),
          file.path(output_dir, "table3_variance_decomp.csv"),
          row.names = FALSE)

note3 <- paste(
  "Variance components shown as percentage of total variance.",
  "Spatial = variance explained by spatial Gaussian process;",
  "Residual = unexplained variance;",
  "Intercept/Slope = proportion of spatial variance attributed to",
  "intercept vs.\\ slope components."
)

emit_standalone_tex(
  table3_df,
  path      = "manuscript/tables/table3_variance_decomp.tex",
  caption   = "Variance decomposition for spatial models",
  label     = "tab:variance_decomp",
  align     = c("l", rep("c", 4)),
  note      = note3,
  size_macro = "small",
  source_script = "extract_full_model_analysis.R"
)
cat("wrote manuscript/tables/table3_variance_decomp.tex\n")

#───────────────────────────────────────────────────────────────────────────────
# 10. TABLE 4 PRODUCER — manuscript/tables/table4_environmental.tex
#───────────────────────────────────────────────────────────────────────────────
#
# Columns: Model, beta_C4, beta_tree, beta_shrub, beta_grass, beta_precip.
# Each cell is "mean [lower, upper]" with a 90% credible interval, or "-"
# when the parameter is not included in the model (gated by include_c4,
# include_pft, include_precip in 4d_leaf_wax_spatial_model.stan). The
# predictors are on the Stan model's standardized scale, matching the
# reporting in the prior hand-authored table.
#
# Coefficient draws come from the widened posterior_draws.rds (every β_*
# scalar is in the saved variable list per 4c_fit_models.R:280–289).

cat("\n\n10. TABLE 4 PRODUCER\n")
cat(strrep("-", 60), "\n")

table4_params <- c("beta_c4", "beta_tree", "beta_shrub", "beta_grass", "beta_precip")

.fmt_coef_90 <- function(draws_vec) {
  # Render "mean [q5, q95]" with 3 decimal places matching the prior table.
  q <- quantile(draws_vec, c(0.05, 0.95), names = FALSE)
  sprintf("%.3f [%.3f, %.3f]", mean(draws_vec), q[1], q[2])
}

table4_rows <- lapply(all_models, function(m) {
  draws <- tryCatch(load_draws(m), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  vars_present <- variables(draws)

  row_vals <- sapply(table4_params, function(p) {
    if (!(p %in% vars_present)) return("-")
    d <- as.numeric(as_draws_matrix(subset_draws(draws, variable = p)))
    .fmt_coef_90(d)
  }, USE.NAMES = FALSE)

  data.frame(
    Model                         = .latex_model_name(m),
    `$\\beta_{\\text{C4}}$`       = row_vals[1],
    `$\\beta_{\\text{tree}}$`     = row_vals[2],
    `$\\beta_{\\text{shrub}}$`    = row_vals[3],
    `$\\beta_{\\text{grass}}$`    = row_vals[4],
    `$\\beta_{\\text{precip}}$`   = row_vals[5],
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
})
table4_df <- bind_rows(table4_rows)

# Write a plain-name CSV sibling for the numeric audit (W6e).
table4_csv <- bind_rows(lapply(all_models, function(m) {
  draws <- tryCatch(load_draws(m), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  vars_present <- variables(draws)
  row_list <- list(model = m)
  for (p in table4_params) {
    if (p %in% vars_present) {
      d <- as.numeric(as_draws_matrix(subset_draws(draws, variable = p)))
      q <- quantile(d, c(0.05, 0.95), names = FALSE)
      row_list[[paste0(p, "_mean")]]  <- mean(d)
      row_list[[paste0(p, "_q05")]]   <- q[1]
      row_list[[paste0(p, "_q95")]]   <- q[2]
    } else {
      row_list[[paste0(p, "_mean")]] <- NA_real_
      row_list[[paste0(p, "_q05")]]  <- NA_real_
      row_list[[paste0(p, "_q95")]]  <- NA_real_
    }
  }
  as.data.frame(row_list, stringsAsFactors = FALSE)
}))
write.csv(table4_csv,
          file.path(output_dir, "table4_environmental.csv"),
          row.names = FALSE)

note4 <- paste(
  "Coefficient estimates shown as posterior mean [90\\% credible interval].",
  "$\\beta_{\\text{C4}}$ = C4 grass fraction effect;",
  "$\\beta_{\\text{tree}}$, $\\beta_{\\text{shrub}}$, $\\beta_{\\text{grass}}$",
  "= plant functional type effects;",
  "$\\beta_{\\text{precip}}$ = precipitation effect.",
  "Dashes indicate parameters not included in the model.",
  "All coefficients on the Stan model's standardized $\\delta^2$H scale."
)

emit_standalone_tex(
  table4_df,
  path      = "manuscript/tables/table4_environmental.tex",
  caption   = "Environmental covariate coefficients",
  label     = "tab:environmental",
  align     = c("l", rep("c", 5)),
  note      = note4,
  landscape = TRUE,
  size_macro = "scriptsize",
  source_script = "extract_full_model_analysis.R"
)
cat("wrote manuscript/tables/table4_environmental.tex\n")

#───────────────────────────────────────────────────────────────────────────────
# 11. TABLE 2 PRODUCER — manuscript/tables/table2_global_params_body.tex
#───────────────────────────────────────────────────────────────────────────────
#
# Nine columns: Model, RMSE (‰), R², β_0 (‰), β_OIPC, λ_int (km),
# GP scale (km), Knot Slope SD, Knot Int SD (‰). All rendered as
# "mean [95% CI]"; GP-only columns show "-" for non-spatial models.
#
# Lineage per cell (every one traceable to a draws column):
#   RMSE      sqrt(mean((d2H_wax - mu[draw])^2)) × d2H_sd                  per draw
#   R²        1 - SSres / SStot using per-draw mu                          per draw
#   β_0 (‰)   beta_0[draw] × d2H_sd + d2H_mean                             per draw
#   β_OIPC    beta_oipc[draw]                                              per draw
#   λ_int     lambda_decay[draw] (already km per Stan l.496)               per draw
#   GP scale  ls_intercept_km[draw] (include_gp == 1 only)                 per draw
#   Knot Slope SD  sigma_slope_spatial[draw]                               per draw
#   Knot Int SD (‰) sigma_intercept_spatial[draw] (already ‰ per Stan)    per draw
#
# The producer writes only the body fragment; table2_global_params.tex is
# a hand-maintained wrapper that carries the landscape longtable preamble
# and `\input{}`s this file.

cat("\n\n11. TABLE 2 PRODUCER\n")
cat(strrep("-", 60), "\n")

.fmt_ci95_number <- function(v, digits = 3) {
  q <- quantile(v, c(0.025, 0.975), names = FALSE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"),
          mean(v), q[1], q[2])
}

.get_draws_vec <- function(draws, name) {
  as.numeric(as_draws_matrix(subset_draws(draws, variable = name)))
}

table2_rows <- lapply(all_models, function(m) {
  draws <- tryCatch(load_draws(m), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  sd_m  <- load_stan_data(m)
  vars_present <- variables(draws)

  d2H_sd   <- sd_m$scaling_params$d2H_sd
  d2H_mean <- sd_m$scaling_params$d2H_mean
  y_obs    <- sd_m$d2H_wax

  # Per-draw fit statistics use `mu` (linear-predictor mean, without the
  # observation noise that d2H_rep adds). This matches section 1's point
  # estimate (RMSE / R² evaluated at the posterior mean prediction) when
  # averaged over draws, and gives tight CIs reflecting parameter
  # uncertainty alone — consistent with what the prior hand-authored
  # table 2 reported. `mu` is included in the widened posterior_draws.rds
  # via 4c_fit_models.R:280.
  mu <- as_draws_matrix(subset_draws(draws, variable = "mu"))
  y_mat <- matrix(y_obs, nrow = nrow(mu), ncol = length(y_obs), byrow = TRUE)
  ss_res_draws <- rowSums((y_mat - mu)^2)
  ss_tot <- sum((y_obs - mean(y_obs))^2)
  rmse_draws <- sqrt(rowMeans((y_mat - mu)^2)) * d2H_sd
  r2_draws   <- 1 - ss_res_draws / ss_tot

  # Global intercept on original scale — Stan line 478:
  # intercept_original = beta_0 * d2H_wax_sd_original + d2H_wax_mean_original.
  beta0_draws    <- .get_draws_vec(draws, "beta_0") * d2H_sd + d2H_mean
  beta_oipc_draws <- .get_draws_vec(draws, "beta_oipc")
  lambda_draws   <- .get_draws_vec(draws, "lambda_decay")

  has_gp <- all(c("ls_intercept_km", "sigma_slope_spatial",
                  "sigma_intercept_spatial") %in% vars_present)

  if (has_gp) {
    gp_scale_draws <- .get_draws_vec(draws, "ls_intercept_km")
    knot_slope_sd  <- .get_draws_vec(draws, "sigma_slope_spatial")
    knot_int_sd    <- .get_draws_vec(draws, "sigma_intercept_spatial")
  }

  data.frame(
    Model          = .latex_model_name(m),
    RMSE           = .fmt_ci95_number(rmse_draws,    digits = 1),
    R2             = .fmt_ci95_number(r2_draws,      digits = 3),
    beta_0         = .fmt_ci95_number(beta0_draws,   digits = 1),
    beta_oipc      = .fmt_ci95_number(beta_oipc_draws, digits = 3),
    lambda_int     = .fmt_ci95_number(lambda_draws,  digits = 1),
    gp_scale       = if (has_gp) .fmt_ci95_number(gp_scale_draws, digits = 0) else "-",
    knot_slope_sd  = if (has_gp) sprintf("%.3f", mean(knot_slope_sd))         else "-",
    knot_int_sd    = if (has_gp) sprintf("%.1f", mean(knot_int_sd))            else "-",
    stringsAsFactors = FALSE
  )
})
table2_df <- bind_rows(table2_rows)

# Plain-number CSV for audit. Keeps every posterior summary Table 2 reports.
table2_csv <- bind_rows(lapply(all_models, function(m) {
  draws <- tryCatch(load_draws(m), error = function(e) NULL)
  if (is.null(draws)) return(NULL)
  sd_m  <- load_stan_data(m)
  vars_present <- variables(draws)

  d2H_sd   <- sd_m$scaling_params$d2H_sd
  d2H_mean <- sd_m$scaling_params$d2H_mean
  y_obs    <- sd_m$d2H_wax
  mu       <- as_draws_matrix(subset_draws(draws, variable = "mu"))
  ss_tot   <- sum((y_obs - mean(y_obs))^2)
  y_mat    <- matrix(y_obs, nrow = nrow(mu), ncol = length(y_obs), byrow = TRUE)
  rmse_d   <- sqrt(rowMeans((y_mat - mu)^2)) * d2H_sd
  r2_d     <- 1 - rowSums((y_mat - mu)^2) / ss_tot
  beta0_d  <- .get_draws_vec(draws, "beta_0") * d2H_sd + d2H_mean
  b_oipc_d <- .get_draws_vec(draws, "beta_oipc")
  lam_d    <- .get_draws_vec(draws, "lambda_decay")
  has_gp   <- all(c("ls_intercept_km", "sigma_slope_spatial",
                    "sigma_intercept_spatial") %in% vars_present)
  if (has_gp) {
    gp_d   <- .get_draws_vec(draws, "ls_intercept_km")
    ks_d   <- .get_draws_vec(draws, "sigma_slope_spatial")
    ki_d   <- .get_draws_vec(draws, "sigma_intercept_spatial")
  } else {
    gp_d <- ks_d <- ki_d <- NA_real_
  }
  .summ <- function(v, label) {
    if (all(is.na(v))) return(setNames(list(NA_real_, NA_real_, NA_real_),
                                      paste0(label, c("_mean", "_q025", "_q975"))))
    q <- quantile(v, c(0.025, 0.975), names = FALSE)
    setNames(list(mean(v), q[1], q[2]),
             paste0(label, c("_mean", "_q025", "_q975")))
  }
  as.data.frame(c(list(model = m),
                  .summ(rmse_d,   "rmse_permil"),
                  .summ(r2_d,     "r_squared"),
                  .summ(beta0_d,  "beta_0_permil"),
                  .summ(b_oipc_d, "beta_oipc"),
                  .summ(lam_d,    "lambda_int_km"),
                  .summ(gp_d,     "gp_scale_km"),
                  .summ(ks_d,     "knot_slope_sd"),
                  .summ(ki_d,     "knot_intercept_sd_permil")),
                stringsAsFactors = FALSE)
}))
write.csv(table2_csv,
          file.path(output_dir, "table2_global_params.csv"),
          row.names = FALSE)

# Emit body-only fragment. Wrapper .tex supplies \begin{longtable} + header.
emit_tabular_fragment(
  table2_df,
  path          = "manuscript/tables/table2_global_params_body.tex",
  source_script = "extract_full_model_analysis.R"
)
cat("wrote manuscript/tables/table2_global_params_body.tex\n")

#───────────────────────────────────────────────────────────────────────────────
# SUMMARY REPORT
#───────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("SUMMARY FOR SECTION 3.4\n")
cat(strrep("=", 60), "\n")

# Best performing model
best_model <- performance_df %>% slice(1)
cat("\nBest performing model:", best_model$model, "\n")
cat(sprintf("  LOOIC: %.1f, R²: %.3f, RMSE: %.1f‰\n",
            best_model$looic, best_model$r_squared, best_model$rmse_permil))

# Key findings
cat("\nKey findings:\n")
cat("1. Significant predictors in full_sp:\n")
full_sp_sig <- coefficients_df %>%
  filter(model == "full_sp", significant, !grepl("interact", parameter))
if (nrow(full_sp_sig) > 0) {
  for (i in 1:nrow(full_sp_sig)) {
    cat(sprintf("   - %s: %.3f [%.3f, %.3f]\n",
                full_sp_sig$parameter[i], full_sp_sig$mean[i],
                full_sp_sig$lower_95[i], full_sp_sig$upper_95[i]))
  }
}

cat("\n2. Variance explained by spatial component in full_sp:\n")
if ("full_sp" %in% names(variance_components)) {
  cat(sprintf("   %.1f%% (mostly intercept variation)\n",
              variance_components[["full_sp"]]$var_spatial_total_pct))
}

cat("\n3. Interaction effects:\n")
interact_params <- coefficients_df %>%
  filter(grepl("interact", model), grepl("x|_c4|_tree|_shrub|_grass", parameter),
         significant)
if (nrow(interact_params) > 0) {
  cat("   Significant interactions found:\n")
  for (i in 1:nrow(interact_params)) {
    cat(sprintf("   - %s in %s\n", interact_params$parameter[i], interact_params$model[i]))
  }
} else {
  cat("   No significant interaction effects\n")
}

cat("\nAnalysis complete!\n")