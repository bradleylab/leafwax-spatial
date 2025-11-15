#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3_prep_data.R
#
# Prepare sediment and predictor data for spatial modeling
# Extracts environmental covariates at sediment sites with optimized sampling
# Creates model-ready dataset with all necessary predictors
#
# Input: site_locations.csv, environmental rasters (OIPC, C4, elevation, climate)
# Output: results/3_sediment_ready_for_modeling.rds
#         results/3_sediment_ready_for_modeling.csv
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(sf)
library(terra)

cat("Prep data for Bayesian modeling (OPTIMIZED VERSION WITH CLIMATE)\n")
cat("================================================================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Parameters
#───────────────────────────────────────────────────────────────────────────────
max_radius_deg <- 5
lat_threshold <- 50
elev_threshold <- 1500

# Downsampling parameters - DISABLED for better accuracy
# We now use ALL pixels since aggregation happens in R
MAX_OIPC_PIXELS <- NULL  # No limit
MAX_ELEV_PIXELS <- NULL  # No limit
MAX_PFT_PIXELS <- NULL   # No limit
# C4 already uses all pixels

#───────────────────────────────────────────────────────────────────────────────
# Step 1: Load and clean sediment data (unchanged)
#───────────────────────────────────────────────────────────────────────────────
cat("Step 1: Loading sediment data...\n")

sediment <- read_csv("input_data/global_data_c29.csv", col_types = cols(.default = col_guess())) %>%
  mutate(
    d2H_wax = as.numeric(d2H_wax),
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude),
    d2H_wax_err = as.numeric(d2H_wax_err),
    measured_elevation = as.numeric(elevation)
  ) %>%
  mutate(
    d2H_wax_err = case_when(
      is.na(d2H_wax_err) ~ 3,
      d2H_wax_err == 0 ~ 1,
      TRUE ~ d2H_wax_err
    )
  ) %>%
  filter(
    chain == 29,
    !is.na(d2H_wax),
    !is.na(latitude),
    !is.na(longitude)
  ) %>%
  mutate(
    has_measured_elevation = !is.na(measured_elevation)
  )

cat("  Loaded", nrow(sediment), "sediment records\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 2: Load rasters
#───────────────────────────────────────────────────────────────────────────────
cat("Step 2: Loading rasters...\n")

c4_rast <- rast("results/1_C4_total_mean.tif")
d2h_ann <- rast("input_data/GlobalPrecip/d2h_MA.tif")
d2h_ann_se <- rast("input_data/GlobalPrecip/d2h_se_MA.tif")
elevation_rast <- rast("input_data/elevation_5KMmn_GMTEDmn.tif")
modis_pft <- rast("results/2d_MODIS_PFT_3classes_Downsampled.tif")

# Load TerraClimate rasters
tc_ppt <- rast("results/2f_TerraClimate_ppt_mean_2001_2019.tif")
tc_soil <- rast("results/2f_TerraClimate_soil_mean_2001_2019.tif")
tc_tmax <- rast("results/2f_TerraClimate_tmax_mean_2001_2019.tif")
tc_vpd <- rast("results/2f_TerraClimate_vpd_mean_2001_2019.tif")

cat("  All rasters loaded successfully\n")

# Convert to data frames
cat("  Converting rasters to data frames...\n")
c4_df <- as.data.frame(c4_rast, xy = TRUE, na.rm = FALSE)
names(c4_df)[3] <- "C4_value"

oipc_df <- as.data.frame(d2h_ann, xy = TRUE, na.rm = FALSE)
oipc_se_df <- as.data.frame(d2h_ann_se, xy = TRUE, na.rm = FALSE)
names(oipc_df)[3] <- "OIPC_value"
names(oipc_se_df)[3] <- "OIPC_SE_value"

elevation_df <- as.data.frame(elevation_rast, xy = TRUE, na.rm = FALSE)
names(elevation_df)[3] <- "elevation_value"

pft_df <- as.data.frame(modis_pft, xy = TRUE, na.rm = FALSE)

# Convert TerraClimate rasters to data frames
tc_ppt_df <- as.data.frame(tc_ppt, xy = TRUE, na.rm = FALSE)
tc_soil_df <- as.data.frame(tc_soil, xy = TRUE, na.rm = FALSE)
tc_tmax_df <- as.data.frame(tc_tmax, xy = TRUE, na.rm = FALSE)
tc_vpd_df <- as.data.frame(tc_vpd, xy = TRUE, na.rm = FALSE)
names(tc_ppt_df)[3] <- "ppt_value"
names(tc_soil_df)[3] <- "soil_value"
names(tc_tmax_df)[3] <- "tmax_value"
names(tc_vpd_df)[3] <- "vpd_value"

cat("✓ Converted rasters to data frames\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 3: OPTIMIZED extraction function with smart downsampling
#───────────────────────────────────────────────────────────────────────────────

extract_radius_data_optimized <- function(lon, lat, raster_df, value_col, radius_deg,
                                          max_pixels = NULL, preserve_na_ratio = TRUE) {
  # Calculate distances to all pixels
  dists <- sqrt((raster_df$x - lon)^2 + (raster_df$y - lat)^2)
  
  # Find pixels within radius
  mask <- dists <= radius_deg
  
  if (sum(mask) == 0) {
    return(list(
      n_pixels = 0,
      n_valid_pixels = 0,
      n_na_pixels = 0,
      distances = numeric(0),
      values = numeric(0),
      mean_value = NA_real_,
      median_value = NA_real_,
      was_downsampled = FALSE
    ))
  }
  
  distances <- dists[mask]
  values <- raster_df[[value_col]][mask]
  n_original <- length(values)
  
  # Count NA statistics before any downsampling
  n_na_original <- sum(is.na(values))
  n_valid_original <- sum(!is.na(values))
  na_ratio <- n_na_original / n_original
  
  # Apply downsampling if needed and max_pixels is specified
  was_downsampled <- FALSE
  if (!is.null(max_pixels) && n_original > max_pixels) {
    was_downsampled <- TRUE
    
    # Distance-weighted sampling with inverse square weighting
    weights <- 1 / (distances + 0.01)^2
    
    if (preserve_na_ratio && n_na_original > 0 && n_valid_original > 0) {
      # Sample NA and valid pixels separately to preserve ratio
      na_mask <- is.na(values)
      na_indices <- which(na_mask)
      valid_indices <- which(!na_mask)
      
      # Calculate how many of each type to sample
      n_na_sample <- round(max_pixels * na_ratio)
      n_valid_sample <- max_pixels - n_na_sample
      
      sampled_indices <- c()
      
      # Sample valid pixels
      if (length(valid_indices) > 0 && n_valid_sample > 0) {
        n_to_sample <- min(n_valid_sample, length(valid_indices))
        sampled <- sample(valid_indices, n_to_sample,
                          prob = weights[valid_indices], replace = FALSE)
        sampled_indices <- c(sampled_indices, sampled)
      }
      
      # Sample NA pixels
      if (length(na_indices) > 0 && n_na_sample > 0) {
        n_to_sample <- min(n_na_sample, length(na_indices))
        sampled <- sample(na_indices, n_to_sample,
                          prob = weights[na_indices], replace = FALSE)
        sampled_indices <- c(sampled_indices, sampled)
      }
      
    } else {
      # Simple distance-weighted sampling
      sampled_indices <- sample(1:n_original, max_pixels,
                                prob = weights, replace = FALSE)
    }
    
    # Apply sampling
    distances <- distances[sampled_indices]
    values <- values[sampled_indices]
  }
  
  # Final statistics
  n_valid <- sum(!is.na(values))
  n_na <- sum(is.na(values))
  
  return(list(
    n_pixels = length(values),
    n_valid_pixels = n_valid,
    n_na_pixels = n_na,
    n_pixels_original = n_original,
    distances = distances,
    values = values,
    mean_value = mean(values, na.rm = TRUE),
    median_value = median(values, na.rm = TRUE),
    was_downsampled = was_downsampled
  ))
}

#───────────────────────────────────────────────────────────────────────────────
# Step 4: Extract raster data with OPTIMIZED downsampling
#───────────────────────────────────────────────────────────────────────────────
cat("Step 4: Extracting raster data (with smart downsampling)...\n")
cat("  C4: NO downsampling (preserve coverage)\n")
cat("  OIPC: Max", MAX_OIPC_PIXELS, "pixels per observation\n")
cat("  Elevation: Max", MAX_ELEV_PIXELS, "pixels per observation\n")
cat("  TerraClimate: No downsampling\n\n")

n_points <- nrow(sediment)
progress_interval <- max(1, floor(n_points / 20))

# Extract C4 data - all pixels
cat("Extracting C4 data (keeping all pixels)...\n")
c4_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  c4_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    c4_df, "C4_value", max_radius_deg,
    max_pixels = NULL  # No downsampling for C4
  )
}

# Extract OIPC data
cat("Extracting OIPC data (all pixels)...\n")
oipc_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  oipc_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    oipc_df, "OIPC_value", max_radius_deg,
    max_pixels = NULL
  )
}

# Extract OIPC SE data 
cat("Extracting OIPC SE data (all pixels)...\n")
oipc_se_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  oipc_se_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    oipc_se_df, "OIPC_SE_value", max_radius_deg,
    max_pixels = NULL
  )
}

# Extract elevation data 
cat("Extracting GMTED elevation data (all pixels)...\n")
elevation_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  elevation_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    elevation_df, "elevation_value", max_radius_deg,
    max_pixels = NULL
  )
}

# Extract PFT data 
cat("Extracting MODIS PFT data...\n")
pft_extractions <- vector("list", n_points)

for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  
  # Use same approach but with downsampling
  dists <- sqrt((pft_df$x - sediment$longitude[i])^2 + 
                  (pft_df$y - sediment$latitude[i])^2)
  mask <- dists <= max_radius_deg
  
  if (sum(mask) == 0) {
    pft_extractions[[i]] <- list(
      n_pixels = 0, distances = numeric(0),
      pft_tree = numeric(0), pft_shrub = numeric(0), pft_grass = numeric(0)
    )
    next
  }
  
  distances <- dists[mask]
  pft_data <- pft_df[mask, ]
  
  # Remove pixels where all PFT values are NA
  valid_rows <- !apply(is.na(pft_data[, 3:ncol(pft_data)]), 1, all)
  
  if (sum(valid_rows) == 0) {
    pft_extractions[[i]] <- list(
      n_pixels = 0, distances = numeric(0),
      pft_tree = numeric(0), pft_shrub = numeric(0), pft_grass = numeric(0)
    )
    next
  }
  
  distances <- distances[valid_rows]
  pft_data <- pft_data[valid_rows, ]
  
  # Downsample if needed
  if (FALSE) {    # rm if (nrow(pft_data) > MAX_OIPC_PIXELS) {
    weights <- 1 / (distances + 0.01)^2
    sampled_idx <- sample(1:nrow(pft_data), MAX_OIPC_PIXELS,
                          prob = weights, replace = FALSE)
    distances <- distances[sampled_idx]
    pft_data <- pft_data[sampled_idx, ]
  }
  
  pft_extractions[[i]] <- list(
    n_pixels = length(distances),
    distances = distances,
    pft_tree = pft_data$Tree,
    pft_shrub = pft_data$Shrub,
    pft_grass = pft_data$Grass
  )
}

# Extract TerraClimate precipitation
cat("Extracting TerraClimate precipitation data...\n")
tc_ppt_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  tc_ppt_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    tc_ppt_df, "ppt_value", max_radius_deg,
    max_pixels = NULL
  )
}

# Extract TerraClimate soil moisture
cat("Extracting TerraClimate soil moisture data...\n")
tc_soil_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  tc_soil_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    tc_soil_df, "soil_value", max_radius_deg,
    max_pixels = NULL
  )
}

# Extract TerraClimate temperature
cat("Extracting TerraClimate maximum temperature data...\n")
tc_tmax_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  tc_tmax_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    tc_tmax_df, "tmax_value", max_radius_deg,
    max_pixels = NULL
  )
}

# Extract TerraClimate VPD
cat("Extracting TerraClimate VPD data...\n")
tc_vpd_extractions <- vector("list", n_points)
for (i in 1:n_points) {
  if (i %% progress_interval == 0) {
    cat("  Progress:", round(100 * i / n_points), "%\n")
  }
  tc_vpd_extractions[[i]] <- extract_radius_data_optimized(
    sediment$longitude[i], sediment$latitude[i], 
    tc_vpd_df, "vpd_value", max_radius_deg,
    max_pixels = NULL
  )
}

cat("  Completed all raster extractions\n\n")

# Report downsampling statistics
n_downsampled_oipc <- sum(map_lgl(oipc_extractions, "was_downsampled"))
n_downsampled_elev <- sum(map_lgl(elevation_extractions, "was_downsampled"))

cat("Downsampling summary:\n")
cat("  OIPC: downsampled", n_downsampled_oipc, "of", n_points, "observations\n")
cat("  Elevation: downsampled", n_downsampled_elev, "of", n_points, "observations\n")
cat("  C4: no downsampling applied (preserving coverage)\n")
cat("  TerraClimate: no downsampling applied\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 5: Add extracted data to sediment dataframe
#───────────────────────────────────────────────────────────────────────────────
cat("Step 5: Adding extracted data to dataframe...\n")

# Add C4, OIPC, and elevation data
sediment$c4_distances <- map(c4_extractions, "distances")
sediment$c4_values <- map(c4_extractions, "values")
sediment$c4_n_pixels <- map_int(c4_extractions, "n_pixels")
sediment$c4_n_valid_pixels <- map_int(c4_extractions, "n_valid_pixels")
sediment$c4_n_na_pixels <- map_int(c4_extractions, "n_na_pixels")
sediment$c4_mean <- map_dbl(c4_extractions, "mean_value")

sediment$oipc_distances <- map(oipc_extractions, "distances")
sediment$oipc_values <- map(oipc_extractions, "values")
sediment$oipc_n_pixels <- map_int(oipc_extractions, "n_pixels")
sediment$oipc_mean <- map_dbl(oipc_extractions, "mean_value")

sediment$oipc_se_distances <- map(oipc_se_extractions, "distances")
sediment$oipc_se_values <- map(oipc_se_extractions, "values")
sediment$oipc_se_n_pixels <- map_int(oipc_se_extractions, "n_pixels")
sediment$oipc_se_mean <- map_dbl(oipc_se_extractions, "mean_value")

# Add GMTED elevation data
sediment$elevation_distances <- map(elevation_extractions, "distances")
sediment$elevation_values <- map(elevation_extractions, "values")
sediment$elevation_n_pixels <- map_int(elevation_extractions, "n_pixels")
sediment$elevation_n_valid_pixels <- map_int(elevation_extractions, "n_valid_pixels")
sediment$elevation_mean <- map_dbl(elevation_extractions, "mean_value")

# Add MODIS PFT data (3-category system: Tree, Shrub, Grass)
sediment$pft_distances <- map(pft_extractions, "distances")
sediment$pft_tree <- map(pft_extractions, "pft_tree")
sediment$pft_shrub <- map(pft_extractions, "pft_shrub")
sediment$pft_grass <- map(pft_extractions, "pft_grass")
sediment$pft_n_pixels <- map_int(pft_extractions, "n_pixels")

# Add TerraClimate data
sediment$tc_ppt_distances <- map(tc_ppt_extractions, "distances")
sediment$tc_ppt_values <- map(tc_ppt_extractions, "values")
sediment$tc_ppt_n_pixels <- map_int(tc_ppt_extractions, "n_pixels")
sediment$tc_ppt_mean <- map_dbl(tc_ppt_extractions, "mean_value")

sediment$tc_soil_distances <- map(tc_soil_extractions, "distances")
sediment$tc_soil_values <- map(tc_soil_extractions, "values")
sediment$tc_soil_n_pixels <- map_int(tc_soil_extractions, "n_pixels")
sediment$tc_soil_mean <- map_dbl(tc_soil_extractions, "mean_value")

sediment$tc_tmax_distances <- map(tc_tmax_extractions, "distances")
sediment$tc_tmax_values <- map(tc_tmax_extractions, "values")
sediment$tc_tmax_n_pixels <- map_int(tc_tmax_extractions, "n_pixels")
sediment$tc_tmax_mean <- map_dbl(tc_tmax_extractions, "mean_value")

sediment$tc_vpd_distances <- map(tc_vpd_extractions, "distances")
sediment$tc_vpd_values <- map(tc_vpd_extractions, "values")
sediment$tc_vpd_n_pixels <- map_int(tc_vpd_extractions, "n_pixels")
sediment$tc_vpd_mean <- map_dbl(tc_vpd_extractions, "mean_value")

cat("  Added extraction results to dataframe\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 6: Apply C4 imputation using GMTED elevation
#───────────────────────────────────────────────────────────────────────────────
cat("Step 6: Applying C4 imputation rules using GMTED elevation...\n")

# Function to apply ecological rules using GMTED elevation data
apply_c4_ecological_rules_gmted <- function(values, distances, point_lat, elevation_values, elevation_distances) {
  if (length(values) == 0) return(values)
  
  # For each NA pixel, apply ecological rules
  na_mask <- is.na(values)
  
  if (sum(na_mask) > 0) {
    # Get mean elevation for this point from GMTED data
    if (length(elevation_values) > 0 && !all(is.na(elevation_values))) {
      mean_elevation <- mean(elevation_values, na.rm = TRUE)
    } else {
      mean_elevation <- NA
    }
    
    # Apply rules to NA pixels
    for (i in which(na_mask)) {
      # Use point coordinates for latitude rule
      if (abs(point_lat) >= lat_threshold) {
        values[i] <- 0  # High latitude = no C4
      } else if (!is.na(mean_elevation) && mean_elevation >= elev_threshold) {
        values[i] <- 0  # High elevation = no C4
      }
      # If neither rule applies, leave as NA
    }
  }
  
  return(values)
}

# Apply ecological rules to all points
sediment <- sediment %>%
  mutate(
    # Apply ecological rules using GMTED elevation
    c4_values_filled = pmap(list(c4_values, c4_distances, latitude, elevation_values, elevation_distances), 
                            function(c4_vals, c4_dists, lat, elev_vals, elev_dists) {
                              if (length(c4_vals) > 0) {
                                apply_c4_ecological_rules_gmted(c4_vals, c4_dists, lat, elev_vals, elev_dists)
                              } else {
                                c4_vals
                              }
                            }),
    
    # Recalculate statistics after filling
    c4_n_filled_pixels = map_int(c4_values_filled, ~ sum(!is.na(.x))),
    c4_n_still_na = map_int(c4_values_filled, ~ sum(is.na(.x))),
    c4_mean_filled = map_dbl(c4_values_filled, ~ mean(.x, na.rm = TRUE)),
    
    # Determine C4 source for tracking
    c4_source = case_when(
      c4_n_valid_pixels > 0 ~ "raster_extraction",
      abs(latitude) >= lat_threshold ~ "latitude_rule",
      !is.na(elevation_mean) & elevation_mean >= elev_threshold ~ "elevation_rule_gmted",
      TRUE ~ "missing"
    )
  )

# Summary of C4 processing
pixels_before <- sum(sediment$c4_n_valid_pixels)
pixels_after <- sum(sediment$c4_n_filled_pixels)
pixels_filled <- pixels_after - pixels_before

cat("C4 ecological filling summary (using GMTED elevation):\n")
cat("  Valid pixels before filling:", pixels_before, "\n")
cat("  Valid pixels after filling:", pixels_after, "\n")
cat("  Pixels filled with ecological rules:", pixels_filled, "\n")
cat("  Points with complete C4 data:", sum(sediment$c4_n_still_na == 0), "\n\n")

# C4 source summary
c4_source_summary <- sediment %>% count(c4_source)
cat("C4 source summary:\n")
for (i in 1:nrow(c4_source_summary)) {
  cat("  ", c4_source_summary$c4_source[i], ":", c4_source_summary$n[i], "\n")
}
cat("\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 7: Final data preparation for modeling
#───────────────────────────────────────────────────────────────────────────────
cat("Step 7: Preparing final dataset for modeling...\n")

# Create modeling-ready dataset
sediment_ready <- sediment %>%
  filter(
    !is.na(d2H_wax),
    !is.na(latitude),
    !is.na(longitude),
    !is.na(c4_mean_filled),  # Need C4 data
    !is.na(oipc_mean),       # Need OIPC data
    !is.na(elevation_mean),  # Need GMTED elevation data
    !is.na(tc_ppt_mean),     # Need precipitation data
    !is.na(tc_soil_mean),    # Need soil moisture data
    !is.na(tc_tmax_mean),    # Need temperature data
    !is.na(tc_vpd_mean)      # Need VPD data
  ) %>%
  mutate(
    # Create final variables for modeling
    C4_fraction_5deg = c4_mean_filled,
    oipc_d2h20 = oipc_mean,
    oipc_se20 = pmax(oipc_se_mean, 1.0),
    elevation_gmted = elevation_mean,
    # Add TerraClimate variables
    annual_precip = tc_ppt_mean,      # mm/year
    soil_moisture = tc_soil_mean,     # mm
    max_temp = tc_tmax_mean,          # °C
    vpd = tc_vpd_mean / 10            # Convert hPa to kPa
  ) %>%
  # Remove incomplete records
  filter(
    !is.na(C4_fraction_5deg),
    !is.na(oipc_d2h20),
    !is.na(oipc_se20),
    !is.na(elevation_gmted),
    !is.na(annual_precip),
    !is.na(soil_moisture),
    !is.na(max_temp),
    !is.na(vpd)
  )

cat("  Final dataset ready for modeling\n")
cat("  Records ready for modeling:", nrow(sediment_ready), "out of", nrow(sediment), "\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 8: Identify points still needing AWS lookup (OPTIONAL)
#───────────────────────────────────────────────────────────────────────────────
cat("Step 8: Checking if AWS elevation lookup still needed...\n")

# Check if any points still have C4 gaps that might be filled with better elevation data
points_still_missing_c4 <- sediment_ready %>%
  filter(c4_n_still_na > 0, abs(latitude) < lat_threshold)

if (nrow(points_still_missing_c4) > 0) {
  cat("Points that might benefit from higher-resolution elevation (AWS lookup):", nrow(points_still_missing_c4), "\n")
  cat("  These points have C4 gaps that aren't filled by latitude/GMTED elevation rules\n")
  cat("  Consider AWS lookup for these points if higher precision is needed\n")
  
  # Save these points for optional AWS lookup
  write_csv(points_still_missing_c4 %>% 
              select(latitude, longitude, elevation_mean, c4_n_still_na), 
            "results/3_points_for_optional_aws_lookup.csv")
  cat("  Saved to: results/3_points_for_optional_aws_lookup.csv\n")
} else {
  cat("No points need additional elevation data - GMTED raster sufficient!\n")
}
cat("\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 9: Save final datasets
#───────────────────────────────────────────────────────────────────────────────
cat("Step 9: Saving final datasets...\n")

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
  cat("Created results directory\n")
}

# Save full dataset with all extractions
saveRDS(sediment, "results/3_sediment_with_all_extractions.rds")
cat("Saved full dataset: results/3_sediment_with_all_extractions.rds\n")

# Save modeling-ready dataset
saveRDS(sediment_ready, "results/3_sediment_ready_for_modeling.rds")
cat("Saved modeling dataset: results/3_sediment_ready_for_modeling.rds\n")

# Save CSV version for external use (remove list columns)
write_csv(sediment_ready %>% 
            select(-ends_with("_distances"), -ends_with("_values"), -starts_with("pft_")) %>%
            select(-c4_values_filled),  # Remove list columns for CSV
          "results/3_sediment_ready_for_modeling.csv")
cat("✓ Saved CSV version: results/3_sediment_ready_for_modeling.csv\n")

# Save extraction summary
extraction_summary <- sediment %>%
  summarise(
    total_points = n(),
    c4_pixels_found = sum(c4_n_pixels > 0),
    c4_pixels_filled = sum(c4_n_filled_pixels > 0),
    oipc_pixels_found = sum(oipc_n_pixels > 0),
    elevation_pixels_found = sum(elevation_n_pixels > 0),
    pft_pixels_found = sum(pft_n_pixels > 0),
    tc_ppt_pixels_found = sum(tc_ppt_n_pixels > 0),
    tc_soil_pixels_found = sum(tc_soil_n_pixels > 0),
    tc_tmax_pixels_found = sum(tc_tmax_n_pixels > 0),
    tc_vpd_pixels_found = sum(tc_vpd_n_pixels > 0),
    c4_complete = sum(c4_n_still_na == 0),
    ready_for_modeling = nrow(sediment_ready),
    .groups = "drop"
  )

write_csv(extraction_summary, "results/3_extraction_summary.csv")
cat("Saved summary: results/3_extraction_summary.csv\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 10: Final diagnostics and summary
#───────────────────────────────────────────────────────────────────────────────
cat("FINAL PROCESSING SUMMARY\n")
cat("========================\n")
cat("Total sediment points:", nrow(sediment), "\n")
cat("Ready for modeling:", nrow(sediment_ready), "\n")
cat("Data retention:", round(100 * nrow(sediment_ready) / nrow(sediment), 1), "%\n\n")

cat("EXTRACTION RESULTS:\n")
cat("C4 pixels:", sum(sediment$c4_n_pixels), "(no downsampling)\n")
cat("OIPC pixels:", sum(sediment$oipc_n_pixels), "\n")
cat("Elevation pixels:", sum(sediment$elevation_n_pixels), "\n")
cat("PFT pixels:", sum(sediment$pft_n_pixels), "\n")
cat("TerraClimate ppt pixels:", sum(sediment$tc_ppt_n_pixels), "\n")
cat("TerraClimate soil pixels:", sum(sediment$tc_soil_n_pixels), "\n")
cat("TerraClimate tmax pixels:", sum(sediment$tc_tmax_n_pixels), "\n")
cat("TerraClimate vpd pixels:", sum(sediment$tc_vpd_n_pixels), "\n\n")

cat("C4 IMPUTATION (with GMTED elevation):\n")
for (i in 1:nrow(c4_source_summary)) {
  cat("  ", c4_source_summary$c4_source[i], ":", c4_source_summary$n[i], "\n")
}
cat("\n")

cat("VARIABLE RANGES (final dataset):\n")
ranges <- sediment_ready %>%
  summarise(
    d2H_wax_range = paste(round(min(d2H_wax), 1), "to", round(max(d2H_wax), 1)),
    elevation_range = paste(round(min(elevation_gmted), 0), "to", round(max(elevation_gmted), 0), "m"),
    c4_range = paste(round(min(C4_fraction_5deg), 3), "to", round(max(C4_fraction_5deg), 3)),
    oipc_range = paste(round(min(oipc_d2h20), 1), "to", round(max(oipc_d2h20), 1)),
    precip_range = paste(round(min(annual_precip), 0), "to", round(max(annual_precip), 0), "mm"),
    soil_range = paste(round(min(soil_moisture), 0), "to", round(max(soil_moisture), 0), "mm"),
    temp_range = paste(round(min(max_temp), 1), "to", round(max(max_temp), 1), "°C"),
    vpd_range = paste(round(min(vpd), 2), "to", round(max(vpd), 2), "kPa")
  )
cat("  d2H_wax:", ranges$d2H_wax_range, "‰\n")
cat("  GMTED elevation:", ranges$elevation_range, "\n")
cat("  C4 fraction:", ranges$c4_range, "\n")
cat("  OIPC d2H:", ranges$oipc_range, "‰\n")
cat("  Annual precipitation:", ranges$precip_range, "\n")
cat("  Soil moisture:", ranges$soil_range, "\n")
cat("  Max temperature:", ranges$temp_range, "\n")
cat("  VPD:", ranges$vpd_range, "\n\n")

# Compare measured vs GMTED elevations where available
if (sum(sediment_ready$has_measured_elevation) > 0) {
  comparison_data <- sediment_ready %>%
    filter(has_measured_elevation) %>%
    select(measured_elevation, elevation_gmted) %>%
    filter(!is.na(measured_elevation), !is.na(elevation_gmted))
  
  if (nrow(comparison_data) > 0) {
    correlation <- cor(comparison_data$measured_elevation, comparison_data$elevation_gmted)
    rmse <- sqrt(mean((comparison_data$measured_elevation - comparison_data$elevation_gmted)^2))
    
    cat("MEASURED vs GMTED ELEVATION COMPARISON:\n")
    cat("  Points with both:", nrow(comparison_data), "\n")
    cat("  Correlation:", round(correlation, 3), "\n")
    cat("  RMSE:", round(rmse, 0), "meters\n\n")
  }
}

cat("Ready for Bayesian modeling with climate covariates!\n")
cat("Next step: Run 4b_fit_models.R\n")

# Show sample of final data
cat("\nSample of final dataset:\n")
sample_data <- sediment_ready %>%
  select(compilation, location, latitude, longitude, d2H_wax, 
         elevation_gmted, C4_fraction_5deg, oipc_d2h20,
         annual_precip, soil_moisture, max_temp, vpd) %>%
  head(5)
print(sample_data)

cat("\n✓ Data preparation complete with TerraClimate variables!\n")