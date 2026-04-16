#───────────────────────────────────────────────────────────────────────────────
# 2c_reproject_modis.R
#
# Reproject MODIS LC_Type1 data to WGS84 and create temporal modal summary
# Processes annual LC_Type1 rasters to create single modal land cover raster
# for plant functional type analysis
#
# Input: annual_mosaics/Global_LC_Type1_*.tif (annual LC_Type1 in Sinusoidal)
# Output: results/2_Global_LC_Type1_Modal_2001_2019_WGS84.tif (modal LC in WGS84)
#         annual_mosaics_wgs84/Global_LC_Type1_*_latlon.tif (annual WGS84 rasters)
#───────────────────────────────────────────────────────────────────────────────

library(terra)

cat("MODIS LC_TYPE1 REPROJECTION AND MODAL CALCULATION\n")
cat("================================================\n\n")

# Load configuration
source("0_load_config.R")

# Create output directory
dir.create(CONFIG$modis$wgs84_dir, showWarnings = FALSE)

# Get all the sinusoidal files from config path
sin_files <- list.files(CONFIG$modis$mosaic_dir, 
                       pattern = "Global_LC_Type1_.*\\.tif$", 
                       full.names = TRUE)

if (length(sin_files) == 0) {
  stop("No MODIS files found in ", CONFIG$modis$mosaic_dir, 
       "\nPlease run 2b_create_annual_mosaics.sh first.")
}

cat("Found", length(sin_files), "annual LC_Type1 files\n")
print(basename(sin_files))
cat("\n")

# Step 1: Reproject each file to WGS84
cat("Step 1: Reprojecting to WGS84...\n")
wgs84_files <- c()

for (i in 1:length(sin_files)) {
  # Get the input file
  input_file <- sin_files[i]
  
  # Extract year from filename
  year <- sub(".*Global_LC_Type1_(\\d{4})\\.tif", "\\1", basename(input_file))
  
  # Create output filename using config path
  output_file <- file.path(CONFIG$modis$wgs84_dir, 
                          paste0("Global_LC_Type1_", year, "_latlon.tif"))
  
  # Check if output already exists (skip if so)
  if (file.exists(output_file)) {
    cat("  Skipping", year, "(already exists)\n")
    wgs84_files <- c(wgs84_files, output_file)
    next
  }
  
  cat("  Reprojecting", year, "...")
  
  # Load sinusoidal raster
  lc_sin <- rast(input_file)
  
  # Reproject to WGS84 using nearest neighbor (important for categorical data)
  lc_geo <- project(lc_sin, "EPSG:4326", method = "near")
  
  # Save to WGS84 directory
  writeRaster(lc_geo, output_file, overwrite = TRUE, datatype = "INT1U")
  
  wgs84_files <- c(wgs84_files, output_file)
  cat(" Done!\n")
  
  # Clean up memory
  rm(lc_sin, lc_geo)
  gc()
}

cat("✓ All files reprojected to WGS84\n\n")

# Step 2: Create modal (most frequent) land cover raster
cat("Step 2: Computing modal land cover (2001-2019)...\n")

# Create results directory if needed
dir.create("results", showWarnings = FALSE)

modal_output <- "results/2_Global_LC_Type1_Modal_2001_2019_WGS84.tif"

if (file.exists(modal_output)) {
  cat("Modal raster already exists:", modal_output, "\n")
} else {
  cat("Loading all", length(wgs84_files), "annual rasters...\n")
  
  # Load all reprojected files as a raster stack
  lc_stack <- rast(wgs84_files)
  
  cat("Computing modal class for each pixel...\n")
  
  # find modal pixel over time series
  lc_modal <- modal(lc_stack)
  
  # Set proper names and attributes
  names(lc_modal) <- "LC_Type1_Modal"
  
  # Save modal raster
  writeRaster(lc_modal, modal_output, overwrite = TRUE, datatype = "INT1U")
  
  cat("✓ Modal raster saved:", modal_output, "\n")
}

# Step 3: Verification and summary
cat("\nStep 3: Verification...\n")

# Verify the outputs using config path
wgs84_files_final <- list.files(CONFIG$modis$wgs84_dir, 
                               pattern = ".*_latlon\\.tif$", 
                               full.names = TRUE)
cat("Created", length(wgs84_files_final), "WGS84 files in", CONFIG$modis$wgs84_dir, "\n")

# Check modal raster
if (file.exists(modal_output)) {
  modal_raster <- rast(modal_output)
  cat("Modal raster details:\n")
  cat("  Dimensions:", paste(dim(modal_raster), collapse = " x "), "\n")
  cat("  Extent:", paste(round(as.vector(ext(modal_raster)), 2), collapse = ", "), "\n")
  cat("  CRS:", crs(modal_raster), "\n")
  
  # Check land cover classes present
  lc_freq <- freq(modal_raster)
  lc_classes <- sort(lc_freq$value)
  cat("  Land cover classes present:", paste(lc_classes, collapse = ", "), "\n")
  
  rm(modal_raster)
}

cat("\nMODIS PROCESSING COMPLETE!\n")
cat("=========================\n")
cat("Outputs:\n")
cat("  -", CONFIG$modis$wgs84_dir, "(19 annual rasters in WGS84)\n")
cat("  - results/2_Global_LC_Type1_Modal_2001_2019_WGS84.tif (modal summary)\n")
cat("\nReady for next step: 3_prep_data.R\n")