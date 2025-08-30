#───────────────────────────────────────────────────────────────────────────────
# 2f_process_terraclimate.R
#
# Process TerraClimate NetCDF files to create 2001-2019 mean rasters
# Creates long-term averages for precipitation, soil moisture, temperature, and VPD
#
# Input: input_data/terraclimate_raw/TerraClimate_*.nc (annual files)
# Output: results/2f_TerraClimate_[variable]_mean_2001_2019.tif
#───────────────────────────────────────────────────────────────────────────────

library(terra)

cat("TERRACLIMATE PROCESSING - COMPUTING 2001-2019 MEANS\n")
cat("==================================================\n\n")

# Set up paths
input_dir <- "input_data/terraclimate_raw"
output_dir <- "results"
dir.create(output_dir, showWarnings = FALSE)

# Variables to process
variables <- c("ppt", "soil", "tmax", "vpd")
start_year <- 2001
end_year <- 2019

# Process each variable
for (var in variables) {
  cat(paste0("Processing ", var, "...\n"))
  
  # Get all files for this variable
  pattern <- paste0("TerraClimate_", var, "_.*\\.nc$")
  files <- list.files(input_dir, pattern = pattern, full.names = TRUE)
  
  # Sort to ensure chronological order
  files <- sort(files)
  
  cat("  Found", length(files), "files\n")
  
  if (length(files) != 19) {
    cat("  WARNING: Expected 19 files, found", length(files), "\n")
  }
  
  # Read first file to check structure
  cat("  Reading files...\n")
  first_rast <- rast(files[1])
  
  # TerraClimate files have 12 layers (one per month)
  # We need to compute annual means first, then overall mean
  
  annual_means <- list()
  
  for (i in seq_along(files)) {
    year <- start_year + i - 1
    cat("    Year", year, "...")
    
    # Read the file (12 monthly layers)
    year_rast <- rast(files[i])
    
    # Compute annual mean from monthly data
    if (var == "ppt") {
      # For precipitation, sum across months to get annual total
      annual_val <- sum(year_rast, na.rm = TRUE)
    } else {
      # For other variables, take mean across months
      annual_val <- mean(year_rast, na.rm = TRUE)
    }
    
    annual_means[[i]] <- annual_val
    cat(" done\n")
  }
  
  # Stack all annual means
  cat("  Computing", start_year, "-", end_year, "mean...\n")
  annual_stack <- rast(annual_means)
  
  # Compute long-term mean
  if (var == "ppt") {
    # For precipitation, we want mean annual precipitation
    longterm_mean <- mean(annual_stack, na.rm = TRUE)
  } else {
    # For other variables, simple mean
    longterm_mean <- mean(annual_stack, na.rm = TRUE)
  }
  
  # Set appropriate names
  names(longterm_mean) <- paste0(var, "_mean_", start_year, "_", end_year)
  
  # Save the result
  output_file <- file.path(output_dir, paste0("2f_TerraClimate_", var, "_mean_2001_2019.tif"))
  writeRaster(longterm_mean, output_file, overwrite = TRUE)
  
  cat("  Saved:", output_file, "\n")
  
  # Report statistics
  stats <- global(longterm_mean, fun = c("min", "max", "mean"), na.rm = TRUE)
  cat("  Statistics:\n")
  cat("    Min:", round(stats$min, 2), "\n")
  cat("    Max:", round(stats$max, 2), "\n")
  cat("    Mean:", round(stats$mean, 2), "\n")
  
  # Add units for clarity
  units <- switch(var,
    ppt = "mm/year",
    soil = "mm",
    tmax = "°C",
    vpd = "hPa (× 0.1 = kPa)"
  )
  cat("    Units:", units, "\n\n")
  
  # Clean up memory
  rm(annual_means, annual_stack, longterm_mean)
  gc()
}

# Create a quick visualization
cat("Creating visualization...\n")
library(viridis)

# Read the processed rasters
ppt <- rast(file.path(output_dir, "2f_TerraClimate_ppt_mean_2001_2019.tif"))
soil <- rast(file.path(output_dir, "2f_TerraClimate_soil_mean_2001_2019.tif"))
tmax <- rast(file.path(output_dir, "2f_TerraClimate_tmax_mean_2001_2019.tif"))
vpd <- rast(file.path(output_dir, "2f_TerraClimate_vpd_mean_2001_2019.tif"))

# Create a 2x2 plot
png(file.path(output_dir, "2f_terraclimate_summary.png"), 
    width = 1200, height = 800, res = 150)

par(mfrow = c(2, 2), mar = c(2, 2, 3, 6))

plot(ppt, main = "Mean Annual Precipitation (mm)", 
     col = viridis(100), axes = FALSE, box = FALSE)
plot(soil, main = "Mean Soil Moisture (mm)", 
     col = viridis(100), axes = FALSE, box = FALSE)
plot(tmax, main = "Mean Max Temperature (°C)", 
     col = viridis(100), axes = FALSE, box = FALSE)
plot(vpd, main = "Mean VPD (hPa)", 
     col = viridis(100), axes = FALSE, box = FALSE)

dev.off()

cat("  Saved visualization: results/2f_terraclimate_summary.png\n\n")

# Check resolution matches other rasters
cat("Checking raster properties...\n")
cat("  Dimensions:", dim(ppt)[1:2], "\n")
cat("  Resolution:", res(ppt), "\n")
cat("  Extent:", as.vector(ext(ppt)), "\n")
cat("  CRS:", crs(ppt, describe = TRUE)$name, "\n\n")

# Compare to C4 raster if it exists
if (file.exists("results/1_C4_total_mean.tif")) {
  c4_rast <- rast("results/1_C4_total_mean.tif")
  cat("Comparison with C4 raster:\n")
  cat("  C4 dimensions:", dim(c4_rast)[1:2], "\n")
  cat("  C4 resolution:", res(c4_rast), "\n")
  
  if (!all(dim(c4_rast)[1:2] == dim(ppt)[1:2])) {
    cat("  ⚠ Warning: Dimensions don't match exactly\n")
    cat("  TerraClimate may need resampling to match other rasters\n")
  } else {
    cat("  ✓ Dimensions match!\n")
  }
}

cat("\nTERRACLIMATE PROCESSING COMPLETE!\n")
cat("================================\n")
cat("Output files:\n")
cat("  - results/2f_TerraClimate_ppt_mean_2001_2019.tif\n")
cat("  - results/2f_TerraClimate_soil_mean_2001_2019.tif\n")
cat("  - results/2f_TerraClimate_tmax_mean_2001_2019.tif\n")
cat("  - results/2f_TerraClimate_vpd_mean_2001_2019.tif\n")
cat("\nThese rasters are ready for use in 3_prep_data.R\n")