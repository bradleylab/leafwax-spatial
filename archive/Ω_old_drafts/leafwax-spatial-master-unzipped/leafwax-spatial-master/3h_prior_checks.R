#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3h_prior_checks.R
#
# Test prior parameterizations for spatial model components
# Evaluates tau (range) parameters with different knot configurations
# Addresses data sparsity and prior predictive checks
#
# Input: results/3_sediment_ready_for_modeling.rds, Stan model files
# Output: results/prior_checks/ (prior predictive plots and diagnostics)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(fields)

cat("TAU PARAMETERIZATION EVALUATION FOR SPARSE DATA\n")
cat("===============================================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# SETUP
#───────────────────────────────────────────────────────────────────────────────

set.seed(42)

# Physical plausibility ranges
PLAUSIBLE_SLOPES <- list(
  extreme_min = 0.5,
  typical_min = 0.6,
  common_range = c(0.7, 0.9),
  typical_max = 1.0,
  extreme_max = 1.2
)

PLAUSIBLE_INTERCEPTS <- list(
  extreme_min = -40,
  typical_range = c(-20, 20),
  extreme_max = 40
)

#───────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────

# Create Fibonacci lattice knots
create_fibonacci_knots <- function(n_knots = 100) {
  golden_angle <- pi * (3.0 - sqrt(5.0))
  knot_coords <- matrix(NA, n_knots, 2)
  
  for (i in 1:n_knots) {
    theta <- golden_angle * (i - 1)
    z <- 1 - 2 * (i - 0.5) / n_knots
    radius <- sqrt(1 - z^2)
    
    lat <- asin(z) * 180 / pi
    lon <- (theta %% (2 * pi)) * 180 / pi - 180
    
    knot_coords[i, ] <- c(lon, lat)
  }
  
  colnames(knot_coords) <- c("longitude", "latitude")
  return(knot_coords)
}

# Calculate data density at knots
calculate_knot_density <- function(knot_coords, obs_coords, radius_km = 1000) {
  n_knots <- nrow(knot_coords)
  n_obs <- nrow(obs_coords)
  knot_density <- numeric(n_knots)
  
  cat("Calculating data density at", n_knots, "knots...\n")
  cat("Using radius of", radius_km, "km for 'nearby' observations\n")
  
  for (k in 1:n_knots) {
    if (k %% 20 == 0) cat("  Progress:", round(100 * k / n_knots), "%\n")
    dists_km <- rdist.earth(matrix(knot_coords[k,], nrow = 1), obs_coords, miles = FALSE)
    knot_density[k] <- sum(dists_km <= radius_km)
  }
  
  return(knot_density)
}

# Define tau parameterizations optimized for sparse data
tau_parameterizations <- list(
  # Original very tight (for comparison)
  very_tight = function(density) {
    list(
      slope = 0.05 + 0.15 * (1 - exp(-density/20)),
      intercept = 0.10 + 0.20 * (1 - exp(-density/20))
    )
  },
  
  # Moderate base with gradual increase
  moderate_base = function(density) {
    list(
      slope = 0.3 + 0.5 * (1 - exp(-density/30)),
      intercept = 0.4 + 0.6 * (1 - exp(-density/30))
    )
  },
  
  # Square root scaling - smoother transition
  sqrt_scaling = function(density) {
    list(
      slope = 0.25 + 0.65 * sqrt(density/(density + 20)),
      intercept = 0.35 + 0.65 * sqrt(density/(density + 20))
    )
  },
  
  # Log scaling - rapid initial change then plateau
  log_scaling = function(density) {
    list(
      slope = pmin(0.25 + 0.15 * log1p(density), 1.0),
      intercept = pmin(0.35 + 0.15 * log1p(density), 1.0)
    )
  },
  
  # Two-stage: different treatment for true zeros
  two_stage = function(density) {
    list(
      slope = ifelse(density == 0, 0.35, 0.5 + 0.5 * (1 - exp(-density/25))),
      intercept = ifelse(density == 0, 0.45, 0.6 + 0.4 * (1 - exp(-density/25)))
    )
  },
  
  # Adaptive for sparse data
  sparse_adaptive = function(density) {
    # More permissive than original
    n <- length(density)
    slope <- numeric(n)
    intercept <- numeric(n)
    
    for (i in 1:n) {
      d <- density[i]
      if (d == 0) {
        slope[i] <- 0.40
        intercept[i] <- 0.50
      } else if (d < 10) {
        slope[i] <- 0.40 + 0.20 * d/10
        intercept[i] <- 0.50 + 0.20 * d/10
      } else {
        slope[i] <- 0.60 + 0.40 * (1 - exp(-(d-10)/30))
        intercept[i] <- 0.70 + 0.30 * (1 - exp(-(d-10)/30))
      }
    }
    
    if (n == 1) {
      list(slope = slope[1], intercept = intercept[1])
    } else {
      list(slope = slope, intercept = intercept)
    }
  }
)

# Evaluate flexibility function
evaluate_flexibility <- function(tau_func, knot_density, n_sims = 5000) {
  
  n_knots <- length(knot_density)
  unique_densities <- sort(unique(knot_density))
  
  # PC prior parameters
  sigma_slope_rate <- -log(0.05) / 0.3      # P(sigma > 0.3) = 0.05
  sigma_intercept_rate <- -log(0.05) / 30   # P(sigma > 30‰) = 0.05
  
  # Storage
  slope_ranges <- matrix(NA, length(unique_densities), 3)
  intercept_ranges <- matrix(NA, length(unique_densities), 3)
  
  # For each unique density level
  for (d_idx in 1:length(unique_densities)) {
    dens <- unique_densities[d_idx]
    tau_d <- tau_func(dens)
    
    # Simulate slopes
    slopes_sim <- numeric(n_sims)
    intercepts_sim <- numeric(n_sims)
    
    for (i in 1:n_sims) {
      sigma_slope <- rexp(1, rate = sigma_slope_rate)
      sigma_intercept <- rexp(1, rate = sigma_intercept_rate)
      
      z_slope <- rnorm(1, 0, tau_d$slope)
      z_intercept <- rnorm(1, 0, tau_d$intercept)
      
      slopes_sim[i] <- 0.8 + sigma_slope * z_slope
      intercepts_sim[i] <- sigma_intercept * z_intercept
    }
    
    # Get quantiles
    slope_ranges[d_idx, ] <- quantile(slopes_sim, c(0.05, 0.50, 0.95))
    intercept_ranges[d_idx, ] <- quantile(intercepts_sim, c(0.05, 0.50, 0.95))
  }
  
  # Calculate coverage of plausible ranges
  coverage_stats <- list(
    slope_typical_coverage = mean(slope_ranges[,1] <= PLAUSIBLE_SLOPES$typical_min & 
                                   slope_ranges[,3] >= PLAUSIBLE_SLOPES$typical_max),
    slope_extreme_coverage = mean(slope_ranges[,1] <= PLAUSIBLE_SLOPES$extreme_min & 
                                  slope_ranges[,3] >= PLAUSIBLE_SLOPES$extreme_max),
    intercept_typical_coverage = mean(intercept_ranges[,1] <= PLAUSIBLE_INTERCEPTS$typical_range[1] & 
                                      intercept_ranges[,3] >= PLAUSIBLE_INTERCEPTS$typical_range[2])
  )
  
  return(list(
    slope_ranges = slope_ranges,
    intercept_ranges = intercept_ranges,
    coverage = coverage_stats,
    densities = unique_densities
  ))
}

#───────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

# Load actual data
cat("Loading sediment data...\n")
sediment_data <- readRDS("results/3_sediment_ready_for_modeling.rds")
obs_coords <- as.matrix(sediment_data[, c("longitude", "latitude")])
cat("  Loaded", nrow(obs_coords), "observations\n\n")

# Test multiple knot configurations based on optimization results
knot_configs <- list(
  "75_knots" = 75,
  "100_knots" = 100,
  "125_knots" = 125,
  "150_knots" = 150
)

# Store results for comparison
all_results <- list()
all_densities <- list()

for (config_name in names(knot_configs)) {
  n_knots <- knot_configs[[config_name]]
  
  cat("\n", strrep("=", 60), "\n")
  cat("TESTING WITH", n_knots, "KNOTS\n")
  cat(strrep("=", 60), "\n\n")
  
  # Create knots
  knot_coords <- create_fibonacci_knots(n_knots)
  
  # Calculate density with 1000km radius
  knot_density <- calculate_knot_density(knot_coords, obs_coords, radius_km = 1000)
  
  # Store density info
  all_densities[[config_name]] <- knot_density
  
  # Show density distribution
  cat("\nData density distribution (1000km radius):\n")
  density_table <- table(cut(knot_density, 
                             breaks = c(-1, 0, 5, 10, 20, 50, 100, 200, Inf),
                             labels = c("0", "1-5", "6-10", "11-20", "21-50", "51-100", "101-200", "200+")))
  print(density_table)
  
  cat("\nDensity summary statistics:\n")
  cat("  Knots with 0 observations:", sum(knot_density == 0), 
      sprintf("(%.1f%%)", 100 * sum(knot_density == 0) / n_knots), "\n")
  cat("  Mean density:", round(mean(knot_density), 1), "\n")
  cat("  Median density:", median(knot_density), "\n")
  cat("  Max density:", max(knot_density), "\n")
  
  # Evaluate parameterizations
  results <- list()
  
  cat("\nEvaluating tau parameterizations...\n")
  cat("=====================================\n\n")
  
  for (param_name in names(tau_parameterizations)) {
    cat("Testing:", param_name, "\n")
    
    # Show tau values at key densities
    test_densities <- c(0, 5, 10, 20, 50, 100)
    cat("  Tau values at key densities:\n")
    cat("  Density | Tau_slope | Tau_intercept\n")
    cat("  --------|-----------|---------------\n")
    
    for (d in test_densities) {
      tau_d <- tau_parameterizations[[param_name]](d)
      cat(sprintf("  %7d | %9.3f | %13.3f\n", d, tau_d$slope, tau_d$intercept))
    }
    
    # Evaluate flexibility
    eval_result <- evaluate_flexibility(tau_parameterizations[[param_name]], knot_density)
    results[[param_name]] <- eval_result
    
    # Print summary
    cat("\n  Prior predictive ranges (5th, 50th, 95th percentiles):\n")
    cat("  Density | Slope Range          | Intercept Range\n")
    cat("  --------|----------------------|-------------------------\n")
    
    # Show results for actual density values in data
    shown_densities <- c(0, min(knot_density[knot_density > 0]), 
                        median(knot_density[knot_density > 0]), 
                        max(knot_density))
    
    for (d in unique(round(shown_densities))) {
      if (d %in% eval_result$densities) {
        idx <- which(eval_result$densities == d)
        cat(sprintf("  %7d | [%4.2f, %4.2f, %4.2f] | [%5.1f, %5.1f, %5.1f]\n",
                    d,
                    eval_result$slope_ranges[idx, 1],
                    eval_result$slope_ranges[idx, 2],
                    eval_result$slope_ranges[idx, 3],
                    eval_result$intercept_ranges[idx, 1],
                    eval_result$intercept_ranges[idx, 2],
                    eval_result$intercept_ranges[idx, 3]))
      }
    }
    
    cat("\n  Coverage of plausible ranges:\n")
    cat("    Slopes [0.6, 1.0]:", sprintf("%.0f%%", 100 * eval_result$coverage$slope_typical_coverage), "\n")
    cat("    Slopes [0.5, 1.2]:", sprintf("%.0f%%", 100 * eval_result$coverage$slope_extreme_coverage), "\n")
    cat("    Intercepts [-20, 20]:", sprintf("%.0f%%", 100 * eval_result$coverage$intercept_typical_coverage), "\n")
    cat("\n")
  }
  
  all_results[[config_name]] <- results
}

#───────────────────────────────────────────────────────────────────────────────
# COMPARISON AND RECOMMENDATION
#───────────────────────────────────────────────────────────────────────────────

cat("\n", strrep("=", 60), "\n")
cat("COMPARISON ACROSS KNOT CONFIGURATIONS\n")
cat(strrep("=", 60), "\n\n")

# Compare density distributions
cat("Density distribution comparison:\n")
cat("Density range | 50 knots | 100 knots\n")
cat("--------------|----------|----------\n")

density_breaks <- c(-1, 0, 5, 10, 20, 50, 100, 200, Inf)
for (i in 1:(length(density_breaks)-1)) {
  range_label <- ifelse(density_breaks[i] == -1, "0", 
                       paste0(density_breaks[i]+1, "-", density_breaks[i+1]))
  if (density_breaks[i+1] == Inf) range_label <- paste0(density_breaks[i]+1, "+")
  
  count_50 <- sum(all_densities[["50_knots"]] > density_breaks[i] & 
                  all_densities[["50_knots"]] <= density_breaks[i+1])
  count_100 <- sum(all_densities[["100_knots"]] > density_breaks[i] & 
                   all_densities[["100_knots"]] <= density_breaks[i+1])
  
  cat(sprintf("%-13s | %8d | %9d\n", range_label, count_50, count_100))
}

# Score parameterizations for each knot configuration
cat("\n\nScoring results:\n")
cat("================\n\n")

best_overall <- NULL
best_score <- -Inf

for (config_name in names(all_results)) {
  cat("\nConfiguration:", config_name, "\n")
  cat("-----------------\n")
  
  results <- all_results[[config_name]]
  knot_density <- all_densities[[config_name]]
  
  scores <- data.frame(
    parameterization = names(results),
    stringsAsFactors = FALSE
  )
  
  for (param_name in names(results)) {
    res <- results[[param_name]]
    
    # Score based on:
    # 1. Not too tight at zero density (allows some variation)
    zero_idx <- which(res$densities == 0)
    if (length(zero_idx) > 0) {
      zero_slope_width <- res$slope_ranges[zero_idx, 3] - res$slope_ranges[zero_idx, 1]
      scores[scores$parameterization == param_name, "zero_flexibility"] <- 
        ifelse(zero_slope_width < 0.1, 0,
               ifelse(zero_slope_width < 0.2, 5,
                      ifelse(zero_slope_width < 0.3, 10, 5)))
    }
    
    # 2. Good coverage of plausible ranges
    scores[scores$parameterization == param_name, "coverage"] <- 
      (res$coverage$slope_typical_coverage * 10 + 
       res$coverage$slope_extreme_coverage * 5)
    
    # 3. Reasonable behavior at median density
    med_density <- median(knot_density[knot_density > 0])
    if (!is.na(med_density) && med_density > 0) {
      med_idx <- which.min(abs(res$densities - med_density))
      med_slope_range <- res$slope_ranges[med_idx, 3] - res$slope_ranges[med_idx, 1]
      scores[scores$parameterization == param_name, "median_behavior"] <- 
        ifelse(med_slope_range > 0.3 && med_slope_range < 0.6, 10, 5)
    }
    
    # Total score
    scores[scores$parameterization == param_name, "total"] <- 
      rowSums(scores[scores$parameterization == param_name, 2:4], na.rm = TRUE)
  }
  
  # Sort by score
  scores <- scores[order(scores$total, decreasing = TRUE), ]
  print(scores, row.names = FALSE)
  
  # Track best overall
  if (scores$total[1] > best_score) {
    best_score <- scores$total[1]
    best_overall <- list(
      config = config_name,
      param = scores$parameterization[1],
      n_knots = knot_configs[[config_name]]
    )
  }
}

#───────────────────────────────────────────────────────────────────────────────
# FINAL RECOMMENDATION
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 60), "\n")
cat("FINAL RECOMMENDATION\n")
cat(strrep("=", 60), "\n\n")

cat("Best configuration:\n")
cat("  Number of knots:", best_overall$n_knots, "\n")
cat("  Tau parameterization:", best_overall$param, "\n\n")

# Show recommended tau function details
tau_func <- tau_parameterizations[[best_overall$param]]
res <- all_results[[best_overall$config]][[best_overall$param]]

cat("Recommended tau values and resulting prior ranges:\n")
cat("Density | Tau_slope | Tau_intercept | Slope 95% CI    | Intercept 95% CI\n")
cat("--------|-----------|---------------|-----------------|------------------\n")

show_densities <- c(0, 5, 10, 20, 50, 100, 150)
for (d in show_densities) {
  tau_d <- tau_func(d)
  
  # Find closest density in results
  if (d %in% res$densities) {
    idx <- which(res$densities == d)
  } else {
    idx <- which.min(abs(res$densities - d))
  }
  
  cat(sprintf("%7d | %9.3f | %13.3f | [%4.2f, %4.2f] | [%5.1f, %5.1f]\n",
              d, tau_d$slope, tau_d$intercept,
              res$slope_ranges[idx, 1], res$slope_ranges[idx, 3],
              res$intercept_ranges[idx, 1], res$intercept_ranges[idx, 3]))
}

cat("\nWhy this configuration works best:\n")
cat("1. With", best_overall$n_knots, "knots, you have better data coverage\n")
cat("2. The", best_overall$param, "parameterization provides:\n")
cat("   - Reasonable constraint at zero density (not too wild)\n")
cat("   - Flexibility where data exists\n")
cat("   - Smooth transitions (no sudden jumps)\n")

# Additional implementation details for Stan
cat("\n\nImplementation for Stan model:\n")
cat("-------------------------------\n")
cat("// In transformed data block:\n")
cat("vector[n_pp_knots] tau_spatial_slope;\n")
cat("vector[n_pp_knots] tau_spatial_intercept;\n")
cat("if (include_gp == 1) {\n")
cat("  for (k in 1:n_pp_knots) {\n")

# Generate Stan code based on best parameterization
if (best_overall$param == "sparse_adaptive") {
  cat("    real d = knot_data_density[k];\n")
  cat("    if (d == 0) {\n")
  cat("      tau_spatial_slope[k] = 0.40;\n")
  cat("      tau_spatial_intercept[k] = 0.50;\n")
  cat("    } else if (d < 10) {\n")
  cat("      tau_spatial_slope[k] = 0.40 + 0.20 * d/10;\n")
  cat("      tau_spatial_intercept[k] = 0.50 + 0.20 * d/10;\n")
  cat("    } else {\n")
  cat("      tau_spatial_slope[k] = 0.60 + 0.40 * (1 - exp(-(d-10)/30.0));\n")
  cat("      tau_spatial_intercept[k] = 0.70 + 0.30 * (1 - exp(-(d-10)/30.0));\n")
  cat("    }\n")
} else if (best_overall$param == "sqrt_scaling") {
  cat("    real d = knot_data_density[k];\n")
  cat("    tau_spatial_slope[k] = 0.25 + 0.65 * sqrt(d/(d + 20.0));\n")
  cat("    tau_spatial_intercept[k] = 0.35 + 0.65 * sqrt(d/(d + 20.0));\n")
} else {
  cat("    // Implement the ", best_overall$param, " parameterization\n")
}
cat("  }\n")
cat("}\n")

# Create visualization
cat("\nGenerating visualization...\n")

# Plot tau functions for best configuration
density_seq <- seq(0, 200, by = 1)
plot_data <- data.frame()

for (param_name in names(tau_parameterizations)) {
  tau_vals <- tau_parameterizations[[param_name]](density_seq)
  
  plot_data <- rbind(plot_data, data.frame(
    density = density_seq,
    tau_slope = tau_vals$slope,
    tau_intercept = tau_vals$intercept,
    param = param_name,
    is_best = param_name == best_overall$param
  ))
}

p1 <- ggplot(plot_data, aes(x = density, y = tau_slope, 
                             color = param, 
                             linetype = is_best,
                             size = is_best)) +
  geom_line() +
  geom_rug(data = data.frame(x = all_densities[[best_overall$config]]), 
           aes(x = x), sides = "b", alpha = 0.3, inherit.aes = FALSE) +
  scale_x_sqrt(breaks = c(0, 10, 50, 100, 200)) +
  scale_linetype_manual(values = c("FALSE" = "dotted", "TRUE" = "solid")) +
  scale_size_manual(values = c("FALSE" = 0.8, "TRUE" = 1.5)) +
  labs(x = "Number of observations within 1000km",
       y = "Tau (slope regularization)",
       title = paste("Tau parameterizations with", best_overall$n_knots, "knots"),
       subtitle = paste("Best:", best_overall$param, "(solid line)")) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("tau_comparison_sparse_data.pdf", p1, width = 10, height = 6)

# Save recommendation
recommendation <- list(
  best_config = best_overall,
  tau_function = tau_func,
  all_results = all_results,
  all_densities = all_densities,
  timestamp = Sys.time()
)

saveRDS(recommendation, "tau_sparse_data_recommendation.rds")

cat("\nFiles saved:\n")
cat("  - tau_comparison_sparse_data.pdf\n")
cat("  - tau_sparse_data_recommendation.rds\n")

cat("\nAnalysis complete!\n")