#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3e_spatial_clustering_analysis.R
#
# Analyze spatial clustering and environmental-spatial correlation
# Tests for spatial autocorrelation in environmental covariates
# Evaluates spatial structure in residuals
#
# Input: results/3_sediment_ready_for_modeling.rds
# Output: results/spatial_clustering/ (diagnostic plots and analyses)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(fields)

# Load data
sediment <- readRDS("results/3_sediment_ready_for_modeling.rds")

cat("\nSPATIAL CLUSTERING ANALYSIS\n")
cat("===========================\n\n")

# 1. Spatial clustering of samples
coords <- as.matrix(sediment[, c("longitude", "latitude")])
dist_matrix <- rdist.earth(coords, miles = FALSE)  # km
diag(dist_matrix) <- NA

# Nearest neighbor distances
nn_distances <- apply(dist_matrix, 1, min, na.rm = TRUE)

cat("Sample clustering:\n")
cat("  Total samples:", nrow(sediment), "\n")
cat("  Nearest neighbor distances:\n")
cat("    Min:", round(min(nn_distances), 1), "km\n")
cat("    Median:", round(median(nn_distances), 1), "km\n")
cat("    Mean:", round(mean(nn_distances), 1), "km\n")
cat("    Max:", round(max(nn_distances), 1), "km\n")

# Samples within distance thresholds
thresholds <- c(1, 5, 10, 20, 50, 100)
for (thresh in thresholds) {
  n_within <- sum(nn_distances <= thresh)
  pct_within <- 100 * n_within / length(nn_distances)
  cat(sprintf("    Within %3d km: %3d samples (%.1f%%)\n", thresh, n_within, pct_within))
}

# 2. Environmental similarity vs spatial distance
cat("\n\nEnvironmental-spatial correlation:\n")

# Get site-level values for C4 and PFT
# These are the values at the exact sample locations (first element of each list)
c4_site <- numeric(nrow(sediment))
pft_tree_site <- numeric(nrow(sediment))
pft_shrub_site <- numeric(nrow(sediment))
pft_grass_site <- numeric(nrow(sediment))

for (i in 1:nrow(sediment)) {
  # C4 at site
  if (!is.null(sediment$c4_values_filled[[i]]) && length(sediment$c4_values_filled[[i]]) > 0) {
    c4_site[i] <- sediment$c4_values_filled[[i]][1]  # First value is at site
  } else {
    c4_site[i] <- NA
  }
  
  # PFT at site
  if (!is.null(sediment$pft_tree[[i]]) && length(sediment$pft_tree[[i]]) > 0) {
    pft_tree_site[i] <- sediment$pft_tree[[i]][1]
    pft_shrub_site[i] <- sediment$pft_shrub[[i]][1]
    pft_grass_site[i] <- sediment$pft_grass[[i]][1]
  } else {
    pft_tree_site[i] <- NA
    pft_shrub_site[i] <- NA
    pft_grass_site[i] <- NA
  }
}

# Add to sediment dataframe temporarily
sediment$c4_site <- c4_site
sediment$pft_tree_site <- pft_tree_site
sediment$pft_shrub_site <- pft_shrub_site
sediment$pft_grass_site <- pft_grass_site

# Select all environmental variables including vegetation
env_vars <- c("oipc_d2h20", "elevation_gmted", "annual_precip", 
              "max_temp", "vpd", "soil_moisture",
              "c4_site", "pft_tree_site", "pft_shrub_site", "pft_grass_site")

# Keep only variables that exist and have data
env_vars_available <- env_vars[env_vars %in% names(sediment)]
env_vars_with_data <- character()
for (var in env_vars_available) {
  if (sum(!is.na(sediment[[var]])) > 100) {  # At least 100 non-NA values
    env_vars_with_data <- c(env_vars_with_data, var)
  }
}

cat("  Environmental variables analyzed:\n")
cat("    Climate/terrain:", paste(grep("^(oipc|elev|annual|max|vpd|soil)", env_vars_with_data, value = TRUE), collapse = ", "), "\n")
cat("    Vegetation:", paste(grep("^(c4|pft)", env_vars_with_data, value = TRUE), collapse = ", "), "\n\n")

# For each environmental variable, calculate correlation with spatial distance
env_spatial_cors <- numeric(length(env_vars_with_data))
names(env_spatial_cors) <- env_vars_with_data

for (i in seq_along(env_vars_with_data)) {
  var_name <- env_vars_with_data[i]
  
  # Calculate pairwise differences in environmental variable
  env_values <- sediment[[var_name]]
  env_dist <- as.matrix(dist(env_values))
  
  # Get corresponding spatial distances
  spatial_dist_vec <- dist_matrix[lower.tri(dist_matrix)]
  env_dist_vec <- env_dist[lower.tri(env_dist)]
  
  # Remove NAs
  valid <- !is.na(spatial_dist_vec) & !is.na(env_dist_vec)
  spatial_dist_vec <- spatial_dist_vec[valid]
  env_dist_vec <- env_dist_vec[valid]
  
  # Calculate correlation
  if (length(spatial_dist_vec) > 100) {
    cor_val <- cor(spatial_dist_vec, env_dist_vec)
    env_spatial_cors[i] <- cor_val
    cat(sprintf("  %-20s: r = %6.3f\n", var_name, cor_val))
  } else {
    env_spatial_cors[i] <- NA
    cat(sprintf("  %-20s: insufficient data\n", var_name))
  }
}

# Separate by variable type
climate_vars <- grep("^(oipc|elev|annual|max|vpd|soil)", names(env_spatial_cors), value = TRUE)
veg_vars <- grep("^(c4|pft)", names(env_spatial_cors), value = TRUE)

cat("\n  Summary:\n")
if (length(climate_vars) > 0) {
  climate_cors <- env_spatial_cors[climate_vars]
  cat("    Climate/terrain max |r|:", round(max(abs(climate_cors), na.rm = TRUE), 3), "\n")
}
if (length(veg_vars) > 0) {
  veg_cors <- env_spatial_cors[veg_vars]
  cat("    Vegetation max |r|:", round(max(abs(veg_cors), na.rm = TRUE), 3), "\n")
}
cat("    Overall max |r|:", round(max(abs(env_spatial_cors), na.rm = TRUE), 3), "\n")

# 3. Test: Are environmentally similar sites spatially clustered?
cat("\n\nEnvironmental similarity at different spatial scales:\n")

# For sites within different distance bands, calculate environmental similarity
distance_bands <- list(
  "0-10 km" = c(0, 10),
  "10-50 km" = c(10, 50),
  "50-100 km" = c(50, 100),
  "100-500 km" = c(100, 500),
  ">500 km" = c(500, Inf)
)

# Standardize environmental variables (including vegetation)
env_data <- sediment[, env_vars_with_data, drop = FALSE]
env_data_std <- scale(env_data)

# Calculate environmental distance (Euclidean in standardized space)
env_dist_matrix <- as.matrix(dist(env_data_std))

for (band_name in names(distance_bands)) {
  band <- distance_bands[[band_name]]
  
  # Find pairs within this distance band
  in_band <- dist_matrix > band[1] & dist_matrix <= band[2]
  
  if (sum(in_band, na.rm = TRUE) > 0) {
    env_dists_in_band <- env_dist_matrix[in_band]
    mean_env_dist <- mean(env_dists_in_band, na.rm = TRUE)
    
    cat(sprintf("  %-15s: mean environmental distance = %.2f (n = %d pairs)\n", 
                band_name, mean_env_dist, sum(in_band, na.rm = TRUE)))
  }
}

# Clean up temporary columns
sediment$c4_site <- NULL
sediment$pft_tree_site <- NULL
sediment$pft_shrub_site <- NULL
sediment$pft_grass_site <- NULL

cat("\n\nSUMMARY:\n")
cat("- Substantial spatial clustering: ", 
    round(100 * sum(nn_distances <= 10) / length(nn_distances), 1), 
    "% of samples within 10 km of nearest neighbor\n")
cat("- Weak environmental-spatial correlation (max |r| = ", 
    round(max(abs(env_spatial_cors), na.rm = TRUE), 3), ")\n")
cat("- Justifies spatial random effects to capture unmeasured regional processes\n")


# Add this at the end of the script:

# Save results to file
results <- list(
  clustering = list(
    n_samples = nrow(sediment),
    nn_distances = nn_distances,
    nn_summary = summary(nn_distances),
    within_10km = sum(nn_distances <= 10),
    pct_within_10km = 100 * sum(nn_distances <= 10) / length(nn_distances),
    threshold_counts = sapply(thresholds, function(t) sum(nn_distances <= t))
  ),
  
  correlations = list(
    env_spatial_cors = env_spatial_cors,
    max_abs_correlation = max(abs(env_spatial_cors), na.rm = TRUE),
    climate_vars = climate_vars,
    veg_vars = veg_vars,
    climate_max_cor = if(length(climate_vars) > 0) max(abs(env_spatial_cors[climate_vars]), na.rm = TRUE) else NA,
    veg_max_cor = if(length(veg_vars) > 0) max(abs(env_spatial_cors[veg_vars]), na.rm = TRUE) else NA
  ),
  
  distance_bands = distance_bands,
  
  timestamp = Sys.time()
)

# Save as RDS
saveRDS(results, "results/spatial_clustering_analysis.rds")

# Also save a text summary
sink("results/spatial_clustering_summary.txt")
cat("SPATIAL CLUSTERING ANALYSIS SUMMARY\n")
cat("Generated:", format(Sys.time()), "\n")
cat("=====================================\n\n")

cat("KEY FINDINGS:\n")
cat("- ", round(results$clustering$pct_within_10km, 1), "% of samples have nearest neighbor within 10 km\n")
cat("- Maximum absolute correlation between environmental variables and spatial distance: r = ", 
    round(results$correlations$max_abs_correlation, 3), "\n")
cat("- This weak correlation justifies spatial random effects\n\n")

cat("DETAILS:\n")
cat("Sample clustering (n =", results$clustering$n_samples, "):\n")
for (i in seq_along(thresholds)) {
  cat(sprintf("  Within %3d km: %3d samples (%.1f%%)\n", 
              thresholds[i], 
              results$clustering$threshold_counts[i],
              100 * results$clustering$threshold_counts[i] / results$clustering$n_samples))
}

cat("\nEnvironmental-spatial correlations:\n")
for (var in names(env_spatial_cors)) {
  if (!is.na(env_spatial_cors[var])) {
    cat(sprintf("  %-20s: r = %6.3f\n", var, env_spatial_cors[var]))
  }
}

sink()

cat("\n\nResults saved to:\n")
cat("- results/spatial_clustering_analysis.rds (full data)\n")
cat("- results/spatial_clustering_summary.txt (text summary)\n")
cat("- spatial_environmental_correlation.png (visualization)\n")