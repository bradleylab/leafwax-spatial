#!/usr/bin/env Rscript
# create_figure1_all_versions.R
# Create multiple versions of Figure 1 with different spatial integration approaches
# Run from figures/ directory

library(tidyverse)
library(posterior)
library(cowplot)

# Source common functions + helper API (rds-only, no fit.rds reads).
source("common_functions.R")
source("posterior_helpers.R")
load_project_config()

cat("\nCREATING MULTIPLE VERSIONS OF FIGURE 1\n")
cat("======================================\n\n")

# Create directories if needed
create_directories()

# Load data from the April run via helpers.
sediment  <- load_sediment()
stan_data <- load_stan_data("baseline")
scaling   <- stan_data$scaling_params

# Widened baseline draws (draws_array) — this replaces the old
# readRDS(..../baseline/fit.rds) + cmdstanr chain-CSV access.
draws_baseline <- load_draws("baseline")
lambda_fitted  <- mean(subset_draws(draws_baseline, variable = "effective_scale_km"))
cat("Fitted effective scale from baseline model:", round(lambda_fitted, 1), "km\n\n")

# Initialize results storage
results_list <- list()

# Function to create consistent plot
create_ols_plot <- function(x, y, title_text, subtitle_text = NULL) {
  # Fit model
  model <- lm(y ~ x)
  
  # Create prediction grid
  pred_df <- tibble(
    x_pred = seq(min(x), max(x), length.out = 200)
  ) %>%
    mutate(
      fit = predict(model, newdata = data.frame(x = x_pred)),
      ci = predict(model, newdata = data.frame(x = x_pred), interval = "confidence"),
      pi = predict(model, newdata = data.frame(x = x_pred), interval = "prediction")
    ) %>%
    mutate(
      ci_lwr = ci[, "lwr"],
      ci_upr = ci[, "upr"],
      pi_lwr = pi[, "lwr"],
      pi_upr = pi[, "upr"]
    ) %>%
    select(-ci, -pi)
  
  # Extract metrics
  msum <- summary(model)
  n_obs <- length(model$residuals)
  r2 <- msum$r.squared
  slope <- coef(model)[2]
  intercept <- coef(model)[1]
  rmse <- sqrt(mean(residuals(model)^2))
  
  # Create subtitle if not provided
  if (is.null(subtitle_text)) {
    subtitle_text <- sprintf("y = %.1f + %.3fx, R² = %.3f, RMSE = %.1f‰, n = %d", 
                            intercept, slope, r2, rmse, n_obs)
  }
  
  # Create plot
  p <- ggplot(data.frame(x = x, y = y), aes(x, y)) +
    theme_minimal(base_size = 14) +
    
    # Prediction interval (95% PI) - lightest blue
    geom_ribbon(
      data = pred_df,
      aes(x = x_pred, ymin = pi_lwr, ymax = pi_upr),
      inherit.aes = FALSE, 
      fill = "#377EB8",
      alpha = 0.15
    ) +
    
    # Confidence interval (95% CI) - darker blue
    geom_ribbon(
      data = pred_df,
      aes(x = x_pred, ymin = ci_lwr, ymax = ci_upr),
      inherit.aes = FALSE, 
      fill = "#377EB8",
      alpha = 0.3
    ) +
    
    # Data points
    geom_point(
      shape = 21,
      fill = "#377EB8",
      color = "white",
      size = 2.5,
      stroke = 0.5,
      alpha = 0.7
    ) +
    
    # Regression line
    geom_line(
      data = pred_df,
      aes(x = x_pred, y = fit),
      inherit.aes = FALSE, 
      color = "#377EB8",
      linewidth = 2
    ) +
    
    # Labels
    labs(
      title = title_text,
      subtitle = subtitle_text,
      x = expression(delta^2 * H[precip] * " (‰)"),
      y = expression(delta^2 * H[wax] * " (‰)")
    ) +
    
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02))
  
  # Return plot and metrics
  return(list(
    plot = p,
    n = n_obs,
    intercept = intercept,
    slope = slope,
    r2 = r2,
    rmse = rmse
  ))
}

# Version i: Point fitting only
cat("Creating version (i): Point fitting only\n")
# Remove NAs for point fitting
valid_idx <- !is.na(sediment$oipc_d2h20)
x_point <- sediment$oipc_d2h20[valid_idx]
y_point <- sediment$d2H_wax[valid_idx]

result_i <- create_ols_plot(x_point, y_point, "OLS regression (point values)")
ggsave("../figures/main_figs/Figure1a_point_fitting.png", result_i$plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("../figures/main_figs/Figure1a_point_fitting.pdf", result_i$plot, width = 10, height = 8, bg = "white")

results_list$point <- data.frame(
  method = "Point values only",
  n = result_i$n,
  intercept = result_i$intercept,
  slope = result_i$slope,
  r_squared = result_i$r2,
  rmse = result_i$rmse,
  spatial_scale = NA
)

# Version ii: Spatial scale of 10 km only
cat("Creating version (ii): 10 km spatial scale\n")
# Find which scale is closest to 10 km
scale_idx <- which.min(abs(stan_data$distance_scales - 10))
cat("  Using scale index", scale_idx, "=", stan_data$distance_scales[scale_idx], "km\n")

oipc_10km <- stan_data$oipc_values[, scale_idx]
oipc_10km_orig <- oipc_10km * scaling$oipc_sd + scaling$oipc_mean
y_all <- sediment$d2H_wax

result_ii <- create_ols_plot(oipc_10km_orig, y_all, "OLS regression (10 km integration)")
ggsave("../figures/main_figs/Figure1b_10km_scale.png", result_ii$plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("../figures/main_figs/Figure1b_10km_scale.pdf", result_ii$plot, width = 10, height = 8, bg = "white")

results_list$km10 <- data.frame(
  method = "10 km integration",
  n = result_ii$n,
  intercept = result_ii$intercept,
  slope = result_ii$slope,
  r_squared = result_ii$r2,
  rmse = result_ii$rmse,
  spatial_scale = 10
)

# Version iii: Equal weights across all scales
cat("Creating version (iii): Equal weights across scales\n")
equal_weights <- rep(1/stan_data$n_scales, stan_data$n_scales)
oipc_equal <- stan_data$oipc_values %*% equal_weights
oipc_equal_orig <- as.numeric(oipc_equal * scaling$oipc_sd + scaling$oipc_mean)

result_iii <- create_ols_plot(oipc_equal_orig, y_all, "OLS regression (equal weights)")
ggsave("../figures/main_figs/Figure1c_equal_weights.png", result_iii$plot, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("../figures/main_figs/Figure1c_equal_weights.pdf", result_iii$plot, width = 10, height = 8, bg = "white")

results_list$equal <- data.frame(
  method = "Equal weights (1-400 km)",
  n = result_iii$n,
  intercept = result_iii$intercept,
  slope = result_iii$slope,
  r_squared = result_iii$r2,
  rmse = result_iii$rmse,
  spatial_scale = mean(stan_data$distance_scales)  # Average scale
)

# Version iv: Bayesian fit with fitted lambda
cat("Creating version (iv): Bayesian baseline model\n")
# Get predictions from Bayesian model (draws_array -> draws_matrix -> colMeans)
y_rep <- as_draws_matrix(subset_draws(draws_baseline, variable = "d2H_rep"))
mu    <- as_draws_matrix(subset_draws(draws_baseline, variable = "mu"))
y_pred <- colMeans(mu)

# Back-transform
y_obs_orig <- stan_data$d2H_wax * scaling$d2H_sd + scaling$d2H_mean
y_pred_orig <- y_pred * scaling$d2H_sd + scaling$d2H_mean

# Get OIPC weighted by fitted scale
scale_weights_mat <- as_draws_matrix(
  subset_draws(draws_baseline, variable = "scale_weights"))
scale_weights_fitted <- colMeans(scale_weights_mat)
oipc_fitted <- stan_data$oipc_values %*% scale_weights_fitted
oipc_fitted_orig <- as.numeric(oipc_fitted * scaling$oipc_sd + scaling$oipc_mean)

# Get Bayesian parameter posterior means
bayes_params <- summarise_draws(
  subset_draws(draws_baseline, variable = c("beta_0", "beta_oipc", "sigma")),
  mean)
bayes_intercept_std <- bayes_params$mean[bayes_params$variable == "beta_0"]
bayes_slope         <- bayes_params$mean[bayes_params$variable == "beta_oipc"]

# Convert intercept to OIPC=0
oipc_mean_std <- -scaling$oipc_mean / scaling$oipc_sd
bayes_intercept <- (bayes_intercept_std + bayes_slope * oipc_mean_std) * scaling$d2H_sd + scaling$d2H_mean

# Calculate R² and RMSE
bayes_r2 <- 1 - var(y_obs_orig - y_pred_orig) / var(y_obs_orig)
bayes_rmse <- sqrt(mean((y_obs_orig - y_pred_orig)^2))

# Per-observation 95% posterior predictive interval on original δ²H scale.
# Build this alongside x / y *before* sorting, then arrange the entire
# frame by x. Previously pred_interval was attached after arrange(), which
# paired each sorted-x row with the unsorted j-th prediction interval —
# a sort-permutation mismatch that made the ribbon zig-zag between
# unrelated observations.
pred_interval <- t(apply(y_rep * scaling$d2H_sd + scaling$d2H_mean, 2, quantile, c(0.025, 0.975)))
plot_df <- data.frame(
  x       = oipc_fitted_orig,
  y       = y_obs_orig,
  y_pred  = y_pred_orig,
  pi_lwr  = pred_interval[, 1],
  pi_upr  = pred_interval[, 2]
) %>% arrange(x)

# Create the Bayesian plot
p_bayes <- ggplot(plot_df, aes(x, y)) +
  theme_minimal(base_size = 14) +

  # 95% posterior predictive interval — wide because it folds in the
  # residual sigma, not a fit-uncertainty envelope.
  geom_ribbon(
    aes(ymin = pi_lwr, ymax = pi_upr),
    fill  = "#377EB8",
    alpha = 0.15
  ) +
  
  # Data points
  geom_point(
    shape = 21,
    fill = "#377EB8",
    color = "white",
    size = 2.5,
    stroke = 0.5,
    alpha = 0.7
  ) +
  
  # Fitted line
  geom_line(
    aes(y = y_pred),
    color = "#377EB8",
    linewidth = 2
  ) +
  
  # Labels
  labs(
    title = "Bayesian baseline model",
    subtitle = sprintf("y = %.1f + %.3fx, R² = %.3f, RMSE = %.1f‰, n = %d, effective scale = %.1f km", 
                      bayes_intercept, bayes_slope, bayes_r2, bayes_rmse, nrow(plot_df), lambda_fitted),
    x = expression(delta^2 * H[precip] * " (‰)"),
    y = expression(delta^2 * H[wax] * " (‰)")
  ) +
  
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12, color = "gray40"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  
  scale_x_continuous(expand = expansion(mult = 0.02)) +
  scale_y_continuous(expand = expansion(mult = 0.02))

ggsave("../figures/main_figs/Figure1d_bayesian_fitted.png", p_bayes, width = 10, height = 8, dpi = 300, bg = "white")
ggsave("../figures/main_figs/Figure1d_bayesian_fitted.pdf", p_bayes, width = 10, height = 8, bg = "white")

# Combined 2x2 panel version with A/B/C/D labels
combined_fig1 <- plot_grid(
  result_i$plot   + labs(title = NULL, subtitle = result_i$plot$labels$subtitle),
  result_ii$plot  + labs(title = NULL, subtitle = result_ii$plot$labels$subtitle),
  result_iii$plot + labs(title = NULL, subtitle = result_iii$plot$labels$subtitle),
  p_bayes         + labs(title = NULL, subtitle = p_bayes$labels$subtitle),
  labels    = c("A", "B", "C", "D"),
  label_size = 16,
  ncol      = 2,
  nrow      = 2,
  align     = "hv"
)
ggsave("../figures/main_figs/Figure_01_combined.png", combined_fig1,
       width = 16, height = 12, dpi = 300, bg = "white")
ggsave("../figures/main_figs/Figure_01_combined.pdf", combined_fig1,
       width = 16, height = 12, bg = "white")

results_list$bayesian <- data.frame(
  method = "Bayesian (fitted weights)",
  n = nrow(plot_df),
  intercept = bayes_intercept,
  slope = bayes_slope,
  r_squared = bayes_r2,
  rmse = bayes_rmse,
  spatial_scale = lambda_fitted
)

# Combine results and save
results_df <- do.call(rbind, results_list)
# Scratch CSV — lands next to the script, not in the tracked figures dir.
write.csv(results_df, "Figure1_comparison_table.csv", row.names = FALSE)

# Print results
cat("\n")
print(results_df)

cat("\nFigures saved to ../figures/main_figs/:\n")
cat("  - Figure1a_point_fitting.png/pdf\n")
cat("  - Figure1b_10km_scale.png/pdf\n")
cat("  - Figure1c_equal_weights.png/pdf\n")
cat("  - Figure1d_bayesian_fitted.png/pdf\n")
cat("  - Figure_01_combined.png/pdf\n")
cat("  - Figure1_comparison_table.csv\n")

# Show scale weights from Bayesian model
cat("\nBayesian model scale weights:\n")
scale_weight_df <- data.frame(
  scale_km = stan_data$distance_scales,
  weight = round(scale_weights_fitted, 4)
)
print(scale_weight_df)