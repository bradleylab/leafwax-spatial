#───────────────────────────────────────────────────────────────────────────────
# 2d_downsample_modis_to_percentages.R
# 
# Create a downsampled MODIS raster with percentage layers for each land cover type
# Target resolution: similar to OIPC (2083 x 4320)
# Output: Multi-layer raster with % of each IGBP class per pixel
#───────────────────────────────────────────────────────────────────────────────

library(terra)
library(tidyverse)

source("0_load_config.R")


cat("DOWNSAMPLING MODIS TO PERCENTAGE LAYERS\n")
cat("======================================\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 1: Load rasters and check dimensions
#───────────────────────────────────────────────────────────────────────────────

# Load high-resolution MODIS
modis_hires <- rast("results/2_Global_LC_Type1_Modal_2001_2019_WGS84.tif")
cat("High-res MODIS dimensions:", dim(modis_hires)[1:2], "\n")

# Load OIPC for target resolution reference
oipc_ref <- rast(CONFIG$data_sources$oipc_d2h)
cat("OIPC reference dimensions:", dim(oipc_ref)[1:2], "\n")
cat("OIPC resolution:", res(oipc_ref), "\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 2: Calculate aggregation factor
#───────────────────────────────────────────────────────────────────────────────

# Calculate how much to aggregate
target_res <- res(oipc_ref)
current_res <- res(modis_hires)
agg_factor <- round(target_res / current_res)

cat("Aggregation factor:", agg_factor, "\n")
cat("New resolution will be approximately:", current_res * agg_factor, "\n\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 3: Create percentage layers for each IGBP class
#───────────────────────────────────────────────────────────────────────────────

# Define IGBP classes (instead of reading all values which causes memory error)
# Standard IGBP classes 1-17
igbp_classes <- 1:17
cat("Processing IGBP classes:", paste(igbp_classes, collapse = ", "), "\n")
cat("(Using standard IGBP 1-17 range to avoid memory issues)\n\n")

# Function to calculate percentage of each class within aggregated pixels
# This function handles the case where a class might not exist in a pixel
calculate_class_percentage <- function(x, target_class) {
  if (all(is.na(x))) return(0)  # Return 0% instead of NA for missing data
  valid_pixels <- x[!is.na(x)]
  if (length(valid_pixels) == 0) return(0)
  
  percentage <- sum(valid_pixels == target_class) / length(valid_pixels) * 100
  return(round(percentage, 2))  # Round to 2 decimal places to save space
}

# Create percentage rasters for each class
cat("Creating percentage layers (this will take a while)...\n")
percentage_layers <- list()

# Process in chunks to avoid memory issues
for (class_id in igbp_classes) {
  cat("  Processing IGBP class", class_id, "... ")
  start_time <- Sys.time()
  
  # Use terra's aggregate function which handles large rasters efficiently
  pct_layer <- aggregate(modis_hires, 
                         fact = agg_factor,
                         fun = function(x) calculate_class_percentage(x, class_id),
                         filename = paste0("temp_class_", class_id, ".tif"),
                         overwrite = TRUE)
  
  names(pct_layer) <- paste0("IGBP_", sprintf("%02d", class_id), "_pct")
  percentage_layers[[paste0("class_", class_id)]] <- pct_layer
  
  elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "mins")), 1)
  cat("done in", elapsed, "minutes\n")
  
  # Clean up memory
  gc()
}

#───────────────────────────────────────────────────────────────────────────────
# Step 4: Combine into multi-layer raster
#───────────────────────────────────────────────────────────────────────────────

cat("\nCombining into multi-layer raster...\n")

# Load all the temporary files and combine them
temp_files <- paste0("temp_class_", igbp_classes, ".tif")
existing_files <- temp_files[file.exists(temp_files)]

if (length(existing_files) > 0) {
  cat("Loading", length(existing_files), "temporary files...\n")
  temp_rasters <- lapply(existing_files, rast)
  modis_percentages <- rast(temp_rasters)
  
  cat("Final raster dimensions:", dim(modis_percentages)[1:3], "\n")
  cat("Number of layers:", nlyr(modis_percentages), "\n")
} else {
  stop("No temporary files found - aggregation may have failed")
}

#───────────────────────────────────────────────────────────────────────────────
# Step 5: Create PFT aggregated layers 
# We'll just have three - tree, shrub, grass
# C3/C4 will come in via the Luo raster
#───────────────────────────────────────────────────────────────────────────────

cat("\nCreating PFT aggregated layers...\n")

# Define PFT groupings
pft_groups <- list(
  Tree = 1:5,        # All forests
  Shrub = 6:7,       # All shrublands
  Grass = c(8, 9, 10) # Woody savannas + Savannas + Grasslands
  # Note: Removing other categories for now - focus on natural vegetation
)

# Create PFT percentage layers by summing relevant IGBP classes
pft_layers <- list()

for (pft_name in names(pft_groups)) {
  cat("  Creating", pft_name, "layer...\n")
  
  relevant_classes <- pft_groups[[pft_name]]
  class_indices <- paste0("class_", relevant_classes)
  
  # Find which layers exist
  existing_indices <- class_indices[class_indices %in% names(percentage_layers)]
  
  if (length(existing_indices) > 0) {
    if (length(existing_indices) == 1) {
      pft_layer <- percentage_layers[[existing_indices[1]]]
    } else {
      # Sum percentages for multiple classes
      pft_layer <- Reduce("+", percentage_layers[existing_indices])
    }
    names(pft_layer) <- paste0("PFT_", pft_name, "_pct")
    pft_layers[[pft_name]] <- pft_layer
  }
}

# Combine PFT layers
modis_pft <- rast(pft_layers)

#───────────────────────────────────────────────────────────────────────────────
# Step 6: Save outputs
#───────────────────────────────────────────────────────────────────────────────

cat("\nSaving outputs...\n")

# Save full IGBP percentage raster
writeRaster(modis_percentages, 
            "results/2d_MODIS_IGBP_17classes_Downsampled.tif",
            overwrite = TRUE)
cat("✓ Saved: results/2d_MODIS_IGBP_17classes_Downsampled.tif\n")

# Save PFT aggregated raster  
writeRaster(modis_pft,
            "results/2d_MODIS_PFT_3classes_Downsampled.tif", 
            overwrite = TRUE)
cat("✓ Saved: results/2d_MODIS_PFT_3classes_Downsampled.tif\n")

# Clean up temporary files
cat("\nCleaning up temporary files...\n")
temp_files_to_remove <- paste0("temp_class_", igbp_classes, ".tif")
file.remove(temp_files_to_remove[file.exists(temp_files_to_remove)])
cat("✓ Removed temporary files\n")

#───────────────────────────────────────────────────────────────────────────────
# Step 7: Validation and summary
#───────────────────────────────────────────────────────────────────────────────

cat("\nValidation checks...\n")

# Check that percentages sum to ~100% using the saved PFT raster (smaller and more manageable)
tryCatch({
  sample_pixels <- spatSample(modis_pft, size = 1000, na.rm = TRUE, as.df = TRUE)
  if (nrow(sample_pixels) > 0) {
    # Sum all PFT percentages (excluding x,y columns)
    pft_cols <- names(sample_pixels)[!names(sample_pixels) %in% c("x", "y")]
    row_sums <- rowSums(sample_pixels[, pft_cols], na.rm = TRUE)
    cat("Sample pixel PFT percentage sums - Mean:", round(mean(row_sums, na.rm = TRUE), 1), 
        "Range:", round(range(row_sums, na.rm = TRUE), 1), "\n")
  } else {
    cat("Could not sample pixels for validation\n")
  }
}, error = function(e) {
  cat("Validation sampling failed (but files were saved successfully):", e$message, "\n")
})

# File sizes
igbp_size <- tryCatch({
  file.size("results/2d_MODIS_IGBP_17classes_Downsampled.tif") / 1e9
}, error = function(e) NA)

pft_size <- tryCatch({
  file.size("results/2d_MODIS_PFT_3classes_Downsampled.tif") / 1e9
}, error = function(e) NA)

if (!is.na(igbp_size) && !is.na(pft_size)) {
  cat("File sizes - IGBP:", round(igbp_size, 2), "GB, PFT:", round(pft_size, 2), "GB\n")
} else {
  cat("Could not determine file sizes\n")
}
