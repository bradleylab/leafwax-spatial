#!/usr/bin/env Rscript
# analyze_spatial_scale_r2.R
# Analyze R² vs spatial integration scale for OLS regression
# Run from figures/ directory

library(tidyverse)
library(scales)

source("posterior_helpers.R")

cat("\nANALYZING R² VS SPATIAL INTEGRATION SCALE\n")
cat("==========================================\n\n")

# Load sediment + Stan data directly from the April 2026 run rather than
# the pre-April `results/3_sediment_ready_for_modeling.rds` / `prepared_data/`
# layout that the archived script used. Stan data for baseline carries the
# same oipc matrix + scaling block the analysis needs.
sediment <- load_sediment()
stan_data <- load_stan_data("baseline")
scaling <- stan_data$scaling_params

# Get available scales
available_scales <- stan_data$distance_scales
cat("Available scales:", paste(available_scales, collapse = ", "), "km\n")

# Define scales to test (0 means point values)
test_scales <- c(0, seq(1, 100, by = 1))
cat("Testing scales from 0 (point) to 100 km\n\n")

# Function to compute weighted OIPC for a given scale
compute_weighted_oipc <- function(target_scale, oipc_matrix, distance_scales) {
  if (target_scale == 0) {
    # Return point values (first scale, which should be closest to 0)
    return(oipc_matrix[, 1])
  }
  
  # Use exponential weighting centered on target scale
  # This mimics what the Bayesian model does but with fixed lambda
  lambda_decay <- 100 / target_scale  # So effective scale = target_scale
  scale_weights <- exp(-lambda_decay * distance_scales / 100)
  scale_weights <- scale_weights / sum(scale_weights)
  
  # Compute weighted average
  return(as.numeric(oipc_matrix %*% scale_weights))
}

# Initialize results
results <- data.frame(
  scale_km = numeric(),
  n_obs = integer(),
  intercept = numeric(),
  slope = numeric(),
  r_squared = numeric(),
  rmse = numeric(),
  aic = numeric(),
  bic = numeric()
)

# For scale 0, use point values
cat("Processing scale 0 km (point values)...\n")
valid_idx <- !is.na(sediment$oipc_d2h20)
if (sum(valid_idx) > 0) {
  x_point <- sediment$oipc_d2h20[valid_idx]
  y_point <- sediment$d2H_wax[valid_idx]
  
  lm_point <- lm(y_point ~ x_point)
  summ_point <- summary(lm_point)
  
  results <- rbind(results, data.frame(
    scale_km = 0,
    n_obs = length(lm_point$residuals),
    intercept = coef(lm_point)[1],
    slope = coef(lm_point)[2],
    r_squared = summ_point$r.squared,
    rmse = sqrt(mean(residuals(lm_point)^2)),
    aic = AIC(lm_point),
    bic = BIC(lm_point)
  ))
}

# Process other scales
cat("Processing scales 1-100 km...\n")
pb <- txtProgressBar(min = 1, max = length(test_scales) - 1, style = 3)

for (i in 2:length(test_scales)) {
  scale <- test_scales[i]
  setTxtProgressBar(pb, i - 1)
  
  # Compute weighted OIPC
  oipc_weighted <- compute_weighted_oipc(scale, stan_data$oipc_values, stan_data$distance_scales)
  
  # Back-transform to original scale
  oipc_orig <- oipc_weighted * scaling$oipc_sd + scaling$oipc_mean
  y_orig <- sediment$d2H_wax
  
  # Fit OLS
  lm_fit <- lm(y_orig ~ oipc_orig)
  summ <- summary(lm_fit)
  
  # Store results
  results <- rbind(results, data.frame(
    scale_km = scale,
    n_obs = length(lm_fit$residuals),
    intercept = coef(lm_fit)[1],
    slope = coef(lm_fit)[2],
    r_squared = summ$r.squared,
    rmse = sqrt(mean(residuals(lm_fit)^2)),
    aic = AIC(lm_fit),
    bic = BIC(lm_fit)
  ))
}
close(pb)

# Find optimal scale
optimal_r2 <- results[which.max(results$r_squared), ]
optimal_aic <- results[which.min(results$aic), ]
optimal_rmse <- results[which.min(results$rmse), ]

cat("\n\nOptimal scales:\n")
cat("  By R²:", optimal_r2$scale_km, "km (R² =", round(optimal_r2$r_squared, 4), ")\n")
cat("  By AIC:", optimal_aic$scale_km, "km (AIC =", round(optimal_aic$aic, 1), ")\n")
cat("  By RMSE:", optimal_rmse$scale_km, "km (RMSE =", round(optimal_rmse$rmse, 2), "‰)\n")

# Create main plot: R² vs scale
p_r2 <- ggplot(results, aes(x = scale_km, y = r_squared)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(size = 2, color = "black") +
  
  # Highlight optimal point
  geom_point(data = optimal_r2, 
             aes(x = scale_km, y = r_squared),
             color = "red", size = 4, shape = 16) +
  geom_text(data = optimal_r2,
            aes(x = scale_km, y = r_squared, 
                label = paste0("Optimal: ", scale_km, " km\nR² = ", round(r_squared, 3))),
            vjust = -1, hjust = 0.5, color = "red", fontface = "bold") +
  
  # Add reference lines
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray50", alpha = 0.5) +
  geom_vline(xintercept = optimal_r2$scale_km, linetype = "dashed", color = "red", alpha = 0.5) +
  
  labs(
    title = "OLS Model Performance vs Spatial Integration Scale",
    subtitle = "Exponential weighting with varying effective scales",
    x = "Spatial Integration Scale (km)",
    y = expression(R^2)
  ) +
  
  scale_x_continuous(breaks = seq(0, 100, by = 10)) +
  scale_y_continuous(labels = number_format(accuracy = 0.001)) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(color = "gray40"),
    panel.grid.minor = element_blank()
  )

# Save the canonical supplement plot (R² vs integration scale) straight
# into manuscript/figures/supplement_figs/. Combined and zoom variants
# were dropped in the figure-output consolidation — the script remains
# self-contained and a reviewer can resurrect them by commenting the
# save calls back in.
ggsave("../figures/supplement_figs/Figure_S4_ols_spatial_scale.png",
       p_r2, width = 10, height = 7, dpi = 300, bg = "white")
ggsave("../figures/supplement_figs/Figure_S4_ols_spatial_scale.pdf",
       p_r2, width = 10, height = 7, device = cairo_pdf)

# Print summary statistics
cat("\n\nSummary of results:\n")
cat("R² range:", round(min(results$r_squared), 4), "-", round(max(results$r_squared), 4), "\n")
cat("R² improvement from point to optimal:", 
    round((optimal_r2$r_squared - results$r_squared[1]) / results$r_squared[1] * 100, 1), "%\n")
cat("Slope range:", round(min(results$slope), 3), "-", round(max(results$slope), 3), "\n")
cat("RMSE range:", round(min(results$rmse), 1), "-", round(max(results$rmse), 1), "‰\n")

# Compare to the Bayesian model result. `effective_scale_km` lives only
# in the diagnostics summary (not the widened posterior_draws), so pull
# its posterior mean from diagnostics.rds via load_summaries().
baseline_summ <- load_summaries("baseline")
bayes_scale <- baseline_summ$mean[baseline_summ$variable == "effective_scale_km"]
if (length(bayes_scale) != 1) bayes_scale <- NA_real_
cat("\nBayesian model effective scale:", round(bayes_scale, 1), "km\n")

# Find R² at Bayesian scale
bayes_scale_r2 <- results %>% 
  filter(abs(scale_km - bayes_scale) == min(abs(scale_km - bayes_scale))) %>%
  pull(r_squared)
cat("OLS R² at Bayesian scale:", round(bayes_scale_r2, 4), "\n")

cat("\nAnalysis complete! Results saved to supplement/\n")