#───────────────────────────────────────────────────────────────────────────────
# 1_extract_c4_raster.R
#
# Extract and prepare C4 vegetation fraction data from NetCDF
# Processes C4 distribution data to create a mean 2001-2019 raster
# Applies necessary coordinate transformations and saves for downstream use
#
# Input: C4_distribution_NUS_v2.2.nc (C4 vegetation NetCDF file)
# Output: C4_total_mean.tif (mean C4 fraction raster, 2001-2019)
#         corrected_total_c4.png (visualization)
#───────────────────────────────────────────────────────────────────────────────

library(terra)
has_viridis <- requireNamespace("viridis", quietly = TRUE)

cat("C4 VEGETATION RASTER EXTRACTION\n")
cat("===============================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 1: Load the full C4 NetCDF stack
#───────────────────────────────────────────────────────────────────────────────
cat("Step 1: Loading C4 NetCDF data...\n")

source("0_load_config.R")
c4_all <- rast(CONFIG$data_sources$c4_netcdf)
cat("✓ Loaded C4 NetCDF with", nlyr(c4_all), "layers\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 2: Identify TOTAL C4 layers (annual data 2001-2019)
#───────────────────────────────────────────────────────────────────────────────
cat("Step 2: Identifying TOTAL C4 area layers...\n")

total_layers <- grep("C4_area_years=", names(c4_all))
c4_total_raw <- c4_all[[total_layers]]

cat("Found", length(total_layers), "C4 total area layers\n")

# Report which years are present
layer_names <- names(c4_total_raw)
years_present <- gsub("C4_area_years=", "", layer_names)  # strip prefix
cat("Years present in TOTAL C4 layers:\n")
print(years_present)
cat("\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 3: Data quality check and dimensions
#───────────────────────────────────────────────────────────────────────────────
cat("Step 3: Data quality check...\n")

cat("Raw C4 TOTAL area dimensions:\n")
print(dim(c4_total_raw))
cat("Raw C4 TOTAL area extent:\n")
print(ext(c4_total_raw))
cat("\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 4: Compute 2001–2019 temporal mean
#───────────────────────────────────────────────────────────────────────────────
cat("Step 4: Computing temporal mean (2001-2019)...\n")

c4_total_mean_raw <- mean(c4_total_raw, na.rm = TRUE)

# Check for data coverage issues
cat("✓ Computed temporal mean across", length(years_present), "years\n")
cat("Checking data coverage...\n")
total_cells <- ncell(c4_total_mean_raw)
valid_cells <- global(c4_total_mean_raw, fun = function(x) sum(!is.na(x)))[[1]]
coverage_pct <- round(100 * valid_cells / total_cells, 1)
cat("  Data coverage:", coverage_pct, "% (", valid_cells, "/", total_cells, "cells)\n")

if (coverage_pct < 50) {
  cat("    WARNING: Low data coverage detected!\n")
  cat("  This may indicate processing issues with the NetCDF file\n")
}
cat("\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 5: Apply coordinate transformation (transpose for correct orientation)
#───────────────────────────────────────────────────────────────────────────────
cat("Step 5: Applying coordinate transformations...\n")

# Transpose to correct row/col ordering (required for this specific NetCDF format)
c4_total_mean_fixed <- t(c4_total_mean_raw)

# Set correct global extent
ext(c4_total_mean_fixed) <- c(-180, 180, -90, 90)

# Set correct CRS to WGS84
crs(c4_total_mean_fixed) <- "EPSG:4326"

cat("✓ Applied coordinate transformations\n")
cat("Final extent:", as.vector(ext(c4_total_mean_fixed)), "\n")
cat("Final CRS:", crs(c4_total_mean_fixed), "\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 6: Create visualization for quality check
#───────────────────────────────────────────────────────────────────────────────
cat("Step 6: Creating visualization...\n")

# Create results directory if it doesn't exist
if (!dir.exists("results")) {
  dir.create("results")
  cat("  Created results directory\n")
}

# Save diagnostic plot (skipped if viridis not available, e.g. in container)
if (has_viridis) {
  plot(c4_total_mean_fixed, main = "Mean total C4 fraction (2001-2019)",
       col = viridis::viridis(100))
  png("results/1_corrected_total_c4.png", width = 800, height = 400, res = 150)
  plot(c4_total_mean_fixed, main = "Corrected TOTAL C4 fraction (2001-2019)",
       col = viridis::viridis(100))
  dev.off()
  cat("  Saved visualization: results/1_corrected_total_c4.png\n\n")
} else {
  cat("  Skipping visualization (viridis not installed)\n\n")
}

#───────────────────────────────────────────────────────────────────────────────
# Step 7: Save processed raster for downstream use
#───────────────────────────────────────────────────────────────────────────────
cat("Step 7: Saving processed C4 raster...\n")

writeRaster(c4_total_mean_fixed, "results/1_C4_total_mean.tif", overwrite = TRUE)

cat("  Saved: results/1_C4_total_mean.tif\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 8: Summary statistics
#───────────────────────────────────────────────────────────────────────────────
cat("PROCESSING SUMMARY\n")
cat("==================\n")
cat("Input: C4_distribution_NUS_v2.2.nc\n")
cat("Output: results/1_C4_total_mean.tif\n")
cat("Time period:", min(years_present), "-", max(years_present), "\n")
cat("Spatial extent: Global (-180° to 180°, -90° to 90°)\n")
cat("Coordinate system: WGS84 (EPSG:4326)\n")

# Basic statistics
c4_stats <- global(c4_total_mean_fixed, fun = c("min", "max", "mean"), na.rm = TRUE)
cat("C4 fraction range:", round(c4_stats$min, 4), "to", round(c4_stats$max, 4), "\n")
cat("Global mean C4 fraction:", round(c4_stats$mean, 4), "\n\n")

cat("Ready for prep_data\n")