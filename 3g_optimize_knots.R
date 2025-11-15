#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3g_optimize_knots.R
#
# Determine optimal number of knots for spatial basis functions
# Balances data coverage, computational efficiency, and spatial resolution
# Tests different knot configurations and evaluates coverage
#
# Input: results/3_sediment_ready_for_modeling.rds
# Output: results/knot_optimization/ (plots and optimal knot selection)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(fields)

cat("OPTIMAL KNOT ANALYSIS FOR GLOBAL COVERAGE\n")
cat("=========================================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────

# Create Fibonacci lattice knots
create_fibonacci_knots <- function(n_knots) {
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

# Analyze knot configuration
analyze_knot_configuration <- function(n_knots, obs_coords, 
                                       radius_km = c(500, 750, 1000, 1500)) {
  
  # Create knots
  knot_coords <- create_fibonacci_knots(n_knots)
  
  # Calculate various metrics
  results <- list(n_knots = n_knots)
  
  # 1. Data coverage at different radii
  coverage_stats <- list()
  for (r in radius_km) {
    knot_density <- numeric(n_knots)
    for (k in 1:n_knots) {
      dists_km <- rdist.earth(matrix(knot_coords[k,], nrow = 1), obs_coords, miles = FALSE)
      knot_density[k] <- sum(dists_km <= r)
    }
    
    coverage_stats[[paste0("r", r)]] <- list(
      n_empty = sum(knot_density == 0),
      n_sparse = sum(knot_density > 0 & knot_density <= 5),
      n_moderate = sum(knot_density > 5 & knot_density <= 20),
      n_dense = sum(knot_density > 20),
      prop_empty = sum(knot_density == 0) / n_knots,
      median_nonzero = ifelse(any(knot_density > 0), 
                              median(knot_density[knot_density > 0]), 
                              NA),
      max_density = max(knot_density)
    )
  }
  results$coverage <- coverage_stats
  
  # 2. Spatial resolution (average distance between knots)
  knot_dists <- rdist.earth(knot_coords, knot_coords, miles = FALSE)
  diag(knot_dists) <- NA
  nearest_neighbor <- apply(knot_dists, 1, min, na.rm = TRUE)
  
  results$resolution <- list(
    mean_nn_dist = mean(nearest_neighbor),
    median_nn_dist = median(nearest_neighbor),
    min_nn_dist = min(nearest_neighbor),
    max_nn_dist = max(nearest_neighbor),
    sd_nn_dist = sd(nearest_neighbor)
  )
  
  # 3. Computational considerations
  results$computational <- list(
    n_kernel_elements = n_knots^2,
    approx_memory_mb = (n_knots^2 * 8) / 1e6,  # rough estimate
    pp_condition = n_knots / nrow(obs_coords)  # knots per observation
  )
  
  # 4. Coverage of specific regions (for paleoclimate)
  # Define key regions
  regions <- list(
    tropics = c(-180, 180, -23.5, 23.5),
    north_temperate = c(-180, 180, 23.5, 66.5),
    south_temperate = c(-180, 180, -66.5, -23.5),
    arctic = c(-180, 180, 66.5, 90),
    antarctic = c(-180, 180, -90, -66.5)
  )
  
  region_coverage <- list()
  for (region_name in names(regions)) {
    bounds <- regions[[region_name]]
    in_region <- knot_coords[,1] >= bounds[1] & knot_coords[,1] <= bounds[2] &
                 knot_coords[,2] >= bounds[3] & knot_coords[,2] <= bounds[4]
    region_coverage[[region_name]] <- sum(in_region)
  }
  results$region_coverage <- region_coverage
  
  # 5. Data clustering analysis
  # Find how observations cluster relative to knots
  obs_to_knot <- matrix(NA, nrow(obs_coords), 1)
  for (i in 1:nrow(obs_coords)) {
    dists <- rdist.earth(matrix(obs_coords[i,], nrow = 1), knot_coords, miles = FALSE)
    obs_to_knot[i] <- min(dists)
  }
  
  results$data_clustering <- list(
    mean_obs_to_knot = mean(obs_to_knot),
    median_obs_to_knot = median(obs_to_knot),
    max_obs_to_knot = max(obs_to_knot),
    prop_within_500km = mean(obs_to_knot <= 500),
    prop_within_1000km = mean(obs_to_knot <= 1000)
  )
  
  return(results)
}

# Calculate theoretical considerations
calculate_theoretical_metrics <- function(n_knots) {
  # Earth's surface area
  earth_surface_km2 <- 510072000
  
  # Average area per knot
  area_per_knot <- earth_surface_km2 / n_knots
  
  # Approximate spacing (assuming uniform distribution)
  approx_spacing <- sqrt(area_per_knot)
  
  # Degrees of freedom relative to typical data size
  typical_n_obs <- 800  # approximate
  dof_ratio <- n_knots / typical_n_obs
  
  return(list(
    area_per_knot_km2 = area_per_knot,
    approx_spacing_km = approx_spacing,
    dof_ratio = dof_ratio
  ))
}

#───────────────────────────────────────────────────────────────────────────────
# LOAD DATA
#───────────────────────────────────────────────────────────────────────────────

cat("Loading sediment data...\n")
sediment_data <- readRDS("results/3_sediment_ready_for_modeling.rds")
obs_coords <- as.matrix(sediment_data[, c("longitude", "latitude")])
n_obs <- nrow(obs_coords)
cat("  Loaded", n_obs, "observations\n\n")

# Show data clustering
cat("Data clustering analysis:\n")
cat("-------------------------\n")

# Simple clustering analysis
obs_dists <- rdist.earth(obs_coords, obs_coords, miles = FALSE)
diag(obs_dists) <- NA

# Count observations within different distances
within_distances <- c(100, 250, 500, 1000, 2000)
for (d in within_distances) {
  nearby_count <- rowSums(obs_dists <= d, na.rm = TRUE)
  cat(sprintf("  Mean obs within %4d km: %5.1f (median: %2.0f, max: %3.0f)\n",
              d, mean(nearby_count), median(nearby_count), max(nearby_count)))
}

# Identify major clusters
# Use simple approach: points with many neighbors
nearby_500km <- rowSums(obs_dists <= 500, na.rm = TRUE)
major_clusters <- which(nearby_500km > quantile(nearby_500km, 0.9))
cat("\n  Major data clusters (>90th percentile density):", length(major_clusters), "locations\n")

#───────────────────────────────────────────────────────────────────────────────
# TEST DIFFERENT KNOT NUMBERS
#───────────────────────────────────────────────────────────────────────────────

# Test range of knot numbers
knot_numbers <- c(25, 30, 40, 50, 60, 75, 100, 125, 150, 200, 250)

cat("\n\nAnalyzing different knot configurations...\n")
cat("==========================================\n\n")

all_results <- list()

for (n_knots in knot_numbers) {
  cat("Testing", n_knots, "knots...\n")
  
  # Analyze configuration
  result <- analyze_knot_configuration(n_knots, obs_coords)
  
  # Add theoretical metrics
  result$theoretical <- calculate_theoretical_metrics(n_knots)
  
  all_results[[as.character(n_knots)]] <- result
  
  # Print summary
  cov_1000 <- result$coverage$r1000
  cat(sprintf("  Coverage (1000km): %d empty (%.1f%%), %d sparse, %d moderate, %d dense\n",
              cov_1000$n_empty, 100 * cov_1000$prop_empty,
              cov_1000$n_sparse, cov_1000$n_moderate, cov_1000$n_dense))
  cat(sprintf("  Resolution: %.0f km mean spacing (%.0f-%.0f km range)\n",
              result$resolution$mean_nn_dist,
              result$resolution$min_nn_dist,
              result$resolution$max_nn_dist))
  cat(sprintf("  Computational: %.1f MB kernel, %.2f knots per observation\n",
              result$computational$approx_memory_mb,
              result$computational$pp_condition))
  cat(sprintf("  Data fit: %.1f%% obs within 500km, %.1f%% within 1000km of knot\n",
              100 * result$data_clustering$prop_within_500km,
              100 * result$data_clustering$prop_within_1000km))
  cat("\n")
}

#───────────────────────────────────────────────────────────────────────────────
# CREATE COMPARISON PLOTS
#───────────────────────────────────────────────────────────────────────────────

cat("\nGenerating comparison plots...\n")

# Extract metrics for plotting
plot_data <- data.frame(
  n_knots = knot_numbers,
  prop_empty_500 = sapply(all_results, function(x) x$coverage$r500$prop_empty),
  prop_empty_1000 = sapply(all_results, function(x) x$coverage$r1000$prop_empty),
  prop_empty_1500 = sapply(all_results, function(x) x$coverage$r1500$prop_empty),
  mean_spacing = sapply(all_results, function(x) x$resolution$mean_nn_dist),
  prop_obs_500 = sapply(all_results, function(x) x$data_clustering$prop_within_500km),
  prop_obs_1000 = sapply(all_results, function(x) x$data_clustering$prop_within_1000km),
  memory_mb = sapply(all_results, function(x) x$computational$approx_memory_mb),
  knots_per_obs = sapply(all_results, function(x) x$computational$pp_condition)
)

# Plot 1: Coverage vs Resolution Trade-off
p1 <- ggplot(plot_data) +
  geom_line(aes(x = n_knots, y = prop_empty_1000 * 100, color = "Empty knots (1000km)"), size = 1.2) +
  geom_line(aes(x = n_knots, y = prop_empty_500 * 100, color = "Empty knots (500km)"), size = 1.2) +
  geom_line(aes(x = n_knots, y = mean_spacing/10, color = "Mean spacing/10 (km)"), size = 1.2, linetype = "dashed") +
  geom_vline(xintercept = c(50, 75, 100), linetype = "dotted", alpha = 0.5) +
  scale_x_continuous(breaks = knot_numbers) +
  scale_color_manual(values = c("Empty knots (1000km)" = "red", 
                                "Empty knots (500km)" = "darkred",
                                "Mean spacing/10 (km)" = "blue")) +
  labs(x = "Number of knots",
       y = "Percentage / Scaled distance",
       title = "Trade-off: Data coverage vs Spatial resolution",
       subtitle = "Lower empty % is better, lower spacing gives better resolution",
       color = "Metric") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Plot 2: Data fit
p2 <- ggplot(plot_data) +
  geom_line(aes(x = n_knots, y = prop_obs_1000 * 100), color = "darkgreen", size = 1.2) +
  geom_line(aes(x = n_knots, y = prop_obs_500 * 100), color = "forestgreen", size = 1.2) +
  geom_hline(yintercept = c(90, 95), linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = c(50, 75, 100), linetype = "dotted", alpha = 0.5) +
  scale_x_continuous(breaks = knot_numbers) +
  scale_y_continuous(limits = c(70, 100)) +
  labs(x = "Number of knots",
       y = "Percentage of observations",
       title = "Observations within distance of nearest knot",
       subtitle = "Higher is better - shows how well knots cover the data") +
  annotate("text", x = 150, y = 98, label = "1000 km", color = "darkgreen") +
  annotate("text", x = 150, y = 88, label = "500 km", color = "forestgreen") +
  theme_minimal()

# Plot 3: Computational cost
p3 <- ggplot(plot_data) +
  geom_line(aes(x = n_knots, y = memory_mb), color = "purple", size = 1.2) +
  geom_point(aes(x = n_knots, y = memory_mb), color = "purple") +
  geom_vline(xintercept = c(50, 75, 100), linetype = "dotted", alpha = 0.5) +
  scale_x_continuous(breaks = knot_numbers) +
  scale_y_log10() +
  labs(x = "Number of knots",
       y = "Kernel memory (MB, log scale)",
       title = "Computational cost",
       subtitle = "Memory requirement for kernel matrices") +
  theme_minimal()

# Combine plots
library(patchwork)
combined_plot <- p1 / p2 / p3
ggsave("optimal_knots_analysis.pdf", combined_plot, width = 10, height = 12)

#───────────────────────────────────────────────────────────────────────────────
# SCORING AND RECOMMENDATION
#───────────────────────────────────────────────────────────────────────────────

cat("\n\nSCORING KNOT CONFIGURATIONS\n")
cat("============================\n\n")

# Score each configuration
scores <- data.frame(n_knots = knot_numbers)

for (i in 1:length(knot_numbers)) {
  n <- knot_numbers[i]
  res <- all_results[[as.character(n)]]
  
  # Score components (0-10 scale)
  
  # 1. Data coverage (want few empty knots at 1000km)
  empty_score <- 10 * (1 - res$coverage$r1000$prop_empty)
  
  # 2. Spatial resolution for paleoclimate (want ~500-1500 km spacing)
  spacing <- res$resolution$mean_nn_dist
  if (spacing < 500) {
    resolution_score <- 5  # Too fine
  } else if (spacing >= 500 && spacing <= 1500) {
    resolution_score <- 10  # Ideal
  } else if (spacing > 1500 && spacing <= 2000) {
    resolution_score <- 7  # Acceptable
  } else {
    resolution_score <- 3  # Too coarse
  }
  
  # 3. Computational feasibility (penalty for too many knots)
  if (n <= 50) {
    comp_score <- 10
  } else if (n <= 100) {
    comp_score <- 8
  } else if (n <= 150) {
    comp_score <- 5
  } else {
    comp_score <- 2
  }
  
  # 4. Data fit (want most observations near a knot)
  fit_score <- 10 * res$data_clustering$prop_within_1000km
  
  # 5. Degrees of freedom (knots should be less than observations)
  if (res$computational$pp_condition < 0.05) {
    dof_score <- 8  # Very few knots relative to data
  } else if (res$computational$pp_condition < 0.1) {
    dof_score <- 10  # Good ratio
  } else if (res$computational$pp_condition < 0.2) {
    dof_score <- 7  # Acceptable
  } else {
    dof_score <- 3  # Too many knots
  }
  
  # Store scores
  scores[i, "empty_penalty"] <- empty_score
  scores[i, "resolution"] <- resolution_score
  scores[i, "computational"] <- comp_score
  scores[i, "data_fit"] <- fit_score
  scores[i, "dof_ratio"] <- dof_score
  
  # Weighted total (adjust weights based on priorities)
  scores[i, "total"] <- 
    0.25 * empty_score +      # Coverage important for inversion
    0.25 * resolution_score + # Resolution critical for paleoclimate
    0.15 * comp_score +       # Computational feasibility
    0.25 * fit_score +        # Must fit the data well
    0.10 * dof_score          # Statistical consideration
}

# Sort by total score
scores <- scores[order(scores$total, decreasing = TRUE), ]

cat("Top configurations:\n")
print(head(scores, 10), row.names = FALSE)

# Get best configuration
best_n <- scores$n_knots[1]
best_result <- all_results[[as.character(best_n)]]

#───────────────────────────────────────────────────────────────────────────────
# DETAILED RECOMMENDATION
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 60), "\n")
cat("RECOMMENDATION: OPTIMAL NUMBER OF KNOTS\n")
cat(strrep("=", 60), "\n\n")

cat("Recommended number of knots:", best_n, "\n\n")

cat("Key metrics for", best_n, "knots:\n")
cat("--------------------------------\n")
cat("Spatial resolution:\n")
cat("  Mean knot spacing:", round(best_result$resolution$mean_nn_dist), "km\n")
cat("  Range:", round(best_result$resolution$min_nn_dist), "-", 
    round(best_result$resolution$max_nn_dist), "km\n")
cat("  Approximate grid:", round(best_result$theoretical$approx_spacing_km), "km\n")

cat("\nData coverage:\n")
cov <- best_result$coverage$r1000
cat("  Empty knots (1000km):", cov$n_empty, sprintf("(%.1f%%)", 100 * cov$prop_empty), "\n")
cat("  Sparse knots (1-5 obs):", cov$n_sparse, "\n")
cat("  Well-covered knots (>20 obs):", cov$n_dense, "\n")

cat("\nData fit:\n")
cat("  Observations within 500km of knot:", 
    sprintf("%.1f%%", 100 * best_result$data_clustering$prop_within_500km), "\n")
cat("  Observations within 1000km of knot:", 
    sprintf("%.1f%%", 100 * best_result$data_clustering$prop_within_1000km), "\n")
cat("  Max distance to nearest knot:", 
    round(best_result$data_clustering$max_obs_to_knot), "km\n")

cat("\nRegional coverage:\n")
for (region in names(best_result$region_coverage)) {
  cat("  ", stringr::str_pad(region, 15), ":", 
      best_result$region_coverage[[region]], "knots\n")
}

cat("\nComputational considerations:\n")
cat("  Kernel matrix size:", round(best_result$computational$approx_memory_mb, 1), "MB\n")
cat("  Knots per observation:", round(best_result$computational$pp_condition, 3), "\n")

# Alternative recommendations
cat("\n\nAlternative configurations to consider:\n")
cat("---------------------------------------\n")
for (i in 2:min(4, nrow(scores))) {
  n_alt <- scores$n_knots[i]
  cat(sprintf("%3d knots (score %.1f):", n_alt, scores$total[i]))
  
  if (n_alt < best_n) {
    cat(" Fewer knots = faster, coarser resolution\n")
  } else {
    cat(" More knots = slower, finer resolution\n")
  }
}

# Save recommendation
recommendation <- list(
  optimal_n_knots = best_n,
  all_results = all_results,
  scores = scores,
  plot_data = plot_data,
  timestamp = Sys.time()
)

saveRDS(recommendation, "optimal_knots_recommendation.rds")

cat("\n\nFiles saved:\n")
cat("  - optimal_knots_analysis.pdf\n")
cat("  - optimal_knots_recommendation.rds\n")

cat("\n\nFor paleoclimate inversions with your data distribution:\n")
cat("  Use", best_n, "knots with larger tau values (see tau analysis)\n")
cat("  This balances global coverage, computational cost, and data fit\n")

cat("\nAnalysis complete!\n")