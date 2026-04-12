#───────────────────────────────────────────────────────────────────────────────
# 4b_stan_prep.R
# Prepare all Stan data with spatial weighting and save to files
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)

cat("STAN DATA PREPARATION WITH SPATIAL WEIGHTING\n")
cat("============================================\n\n")

# Source helper functions
source("0_load_config.R")  # Load configuration first
source(config$scripts$spatial_functions)   # Includes both spatial and validation functions

# Check if input file exists
input_file <- CONFIG$input_data
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, "\n",
       "Please run the previous pipeline steps first.")
}

# Load data
cat("Loading sediment data...\n")
sediment <- readRDS(input_file)
cat("✓ Loaded", nrow(sediment), "sediment records\n")

# Check for required data columns from 3_prep_data.R
required_cols <- c("c4_values_filled", "c4_distances",
                   "oipc_values", "oipc_distances",
                   "oipc_se_values",
                   "elevation_values", "elevation_distances",
                   "latitude", "longitude",
                   "d2H_wax", "d2H_wax_err")
missing_cols <- required_cols[!required_cols %in% names(sediment)]

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
       "\n  Re-run 3_prep_data.R.")
}

# Check for climate data columns (optional, used only if include_* flags are TRUE)
climate_site_cols <- c("annual_precip")

missing_site <- climate_site_cols[!climate_site_cols %in% names(sediment)]
if (length(missing_site) > 0) {
  cat("WARNING: Missing climate columns:", paste(missing_site, collapse = ", "), "\n")
}

# Check for PFT columns
pft_columns <- c("pft_tree", "pft_shrub", "pft_grass")
has_pft_columns <- all(pft_columns %in% names(sediment))

# Additional check: verify these are valid list columns with data
if (has_pft_columns) {
  # Check if they're list columns and have data
  for (col in pft_columns) {
    if (!is.list(sediment[[col]]) || length(sediment[[col]]) == 0) {
      has_pft_columns <- FALSE
      break
    }
  }
}

if (!has_pft_columns) {
  cat("WARNING: PFT columns not available or invalid in data\n")
  cat("Models requesting PFT will have it disabled\n\n")
}

# Add coordinate jitter BEFORE any spatial calculations
set.seed(12345)
jitter_amount <- CONFIG$coordinate_jitter
sediment <- sediment %>%
  mutate(
    longitude_original = longitude,
    latitude_original = latitude,
    longitude = longitude + rnorm(n(), 0, jitter_amount),
    latitude = latitude + rnorm(n(), 0, jitter_amount)
  )
cat("✓ Added coordinate jitter (±", jitter_amount, "°) to handle duplicates\n\n")

# Check for required columns
required_arrays <- c("c4_values_filled", "oipc_values", 
                     "oipc_se_values", "elevation_values")
missing_arrays <- required_arrays[!required_arrays %in% names(sediment)]
if (length(missing_arrays) > 0) {
  stop("Missing required pixel arrays: ", paste(missing_arrays, collapse = ", "))
}

# Create scaling parameters from the full dataset
SCALING_PARAMS <- list(
  d2H_mean = mean(sediment$d2H_wax, na.rm = TRUE),
  d2H_sd = sd(sediment$d2H_wax, na.rm = TRUE),
  oipc_mean = mean(sediment$oipc_d2h20, na.rm = TRUE),
  oipc_sd = sd(sediment$oipc_d2h20, na.rm = TRUE),
  elev_mean = mean(sediment$elevation_gmted, na.rm = TRUE),
  elev_sd = sd(sediment$elevation_gmted, na.rm = TRUE),
  c4_mean = CONFIG$c4_standardization$mean,
  c4_sd = CONFIG$c4_standardization$sd
)

# Add precipitation scaling if needed
if (CONFIG$climate_standardization$compute_from_data && "annual_precip" %in% names(sediment)) {
  SCALING_PARAMS$precip_mean = mean(sediment$annual_precip, na.rm = TRUE)
  SCALING_PARAMS$precip_sd = sd(sediment$annual_precip, na.rm = TRUE)
}

# Add other climate scaling if needed
if (CONFIG$climate_standardization$compute_from_data) {
  if ("max_temp" %in% names(sediment)) {
    SCALING_PARAMS$temp_mean = mean(sediment$max_temp, na.rm = TRUE)
    SCALING_PARAMS$temp_sd = sd(sediment$max_temp, na.rm = TRUE)
  }
  if ("vpd" %in% names(sediment)) {
    SCALING_PARAMS$vpd_mean = mean(sediment$vpd, na.rm = TRUE)
    SCALING_PARAMS$vpd_sd = sd(sediment$vpd, na.rm = TRUE)
  }
  if ("soil_moisture" %in% names(sediment)) {
    SCALING_PARAMS$soil_mean = mean(sediment$soil_moisture, na.rm = TRUE)
    SCALING_PARAMS$soil_sd = sd(sediment$soil_moisture, na.rm = TRUE)
  }
}

# Validate scaling parameters
for (param_name in names(SCALING_PARAMS)) {
  if (is.na(SCALING_PARAMS[[param_name]]) || is.infinite(SCALING_PARAMS[[param_name]])) {
    stop("Invalid scaling parameter: ", param_name, " = ", SCALING_PARAMS[[param_name]])
  }
}

cat("Scaling parameters:\n")
cat("  d2H mean:", round(SCALING_PARAMS$d2H_mean, 2), "‰, sd:", round(SCALING_PARAMS$d2H_sd, 2), "‰\n")
cat("  OIPC mean:", round(SCALING_PARAMS$oipc_mean, 2), "‰, sd:", round(SCALING_PARAMS$oipc_sd, 2), "‰\n")
cat("  Elevation mean:", round(SCALING_PARAMS$elev_mean, 0), "m, sd:", round(SCALING_PARAMS$elev_sd, 0), "m\n")
if ("precip_mean" %in% names(SCALING_PARAMS)) {
  cat("  Precipitation mean:", round(SCALING_PARAMS$precip_mean, 0), "mm, sd:", round(SCALING_PARAMS$precip_sd, 0), "mm\n")
}
cat("\n")

# Get model configurations from config
model_configs <- CONFIG$model_configs

# Create directories for prepared data
for (dir_name in CONFIG$output_dirs) {
  dir.create(dir_name, showWarnings = FALSE, recursive = TRUE)
}

# Prepare data for each model
cat("Preparing Stan data for each model configuration...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

for (model_name in names(model_configs)) {
  config <- model_configs[[model_name]]
  
  cat("Preparing:", config$name, "\n")
  
  # Output configuration details
  cat("  Configuration:\n")
  cat("    C4:", config$include_c4, "\n")
  cat("    PFT:", ifelse(is.null(config$include_pft), FALSE, config$include_pft), "\n")
  cat("    Elevation:", config$include_elevation, "\n")
  cat("    Precipitation:", ifelse(is.null(config$include_precip), FALSE, config$include_precip), "\n")
  cat("    Vegetation interactions:", ifelse(is.null(config$include_veg_interactions), FALSE, config$include_veg_interactions), "\n")
  cat("    Spatial GP:", config$include_gp, "\n")
  
  if (config$include_gp) {
    cat("    PP knots:", config$n_pp_knots, "\n")
  }
  
  # Check if already prepared
  output_file <- file.path(CONFIG$output_dirs$prepared_data, paste0("stan_data_", model_name, ".rds"))
  if (file.exists(output_file)) {
    cat("  ✓ Already prepared - skipping\n\n")
    next
  }
  
  # Prepare data
  tryCatch({
    # Time the preparation
    prep_start <- Sys.time()
    
    # Get parameters from config, with defaults
    include_pft <- ifelse(is.null(config$include_pft), FALSE, config$include_pft)
    include_precip <- ifelse(is.null(config$include_precip), FALSE, config$include_precip)
    include_veg_interactions <- ifelse(is.null(config$include_veg_interactions), FALSE, config$include_veg_interactions)
    
    # Check if PFT is requested but not available
    if (include_pft && !has_pft_columns) {
      cat("  WARNING: PFT requested but not available in data - setting to FALSE\n")
      include_pft <- FALSE
    }
    
    # Call prepare_stan_data with parameters from config
    stan_data <- prepare_stan_data(
      data = sediment,
      include_c4 = config$include_c4,
      include_pft = include_pft,
      include_gp = config$include_gp,
      include_elevation = config$include_elevation,
      include_precip = include_precip,
      include_temp = FALSE,   # Always FALSE due to collinearity
      include_vpd = FALSE,    # Always FALSE due to collinearity
      include_soil = FALSE,   # Always FALSE due to collinearity
      n_pp_knots = config$n_pp_knots,
      SCALING_PARAMS = SCALING_PARAMS,
      has_pft_columns = has_pft_columns
    )
    
    # Add vegetation interaction flag to stan_data
    stan_data$include_veg_interactions <- as.integer(include_veg_interactions)
    
    # Validate stan_data before saving
    stan_data <- validate_stan_data(stan_data)
    
    prep_time <- difftime(Sys.time(), prep_start, units = "secs")
    
    # Validate stan_data structure
    required_elements <- c("N", "n_scales", "distance_scales",
                          "d2H_wax", "oipc_values", "elevation_values",
                          "include_veg_interactions")
    missing_elements <- required_elements[!required_elements %in% names(stan_data)]
    if (length(missing_elements) > 0) {
      stop("Missing required stan_data elements: ", paste(missing_elements, collapse = ", "))
    }
    
    # Save prepared data
    saveRDS(stan_data, output_file)
    
    # Also save config info with additional metadata
    config$scaling_params <- SCALING_PARAMS
    config$data_prep_time <- as.numeric(prep_time)
    config$n_observations <- stan_data$N
    config$pft_available <- has_pft_columns
    config$pft_actually_used <- include_pft
    config$timestamp <- Sys.time()
    config$model_version <- "updated_bspline_matern"
    saveRDS(config, file.path(CONFIG$output_dirs$prepared_data, paste0("config_", model_name, ".rds")))
    
    cat("  ✓ Saved to", output_file, "\n")
    cat("  Data dimensions: N =", stan_data$N, ", scales =", stan_data$n_scales, "\n")
    cat("  Preparation time:", round(prep_time, 1), "seconds\n\n")
    
  }, error = function(e) {
    cat("  ✗ ERROR:", e$message, "\n\n")
    # Save error log
    error_log <- list(
      model = model_name,
      error = e$message,
      traceback = traceback(),
      timestamp = Sys.time()
    )
    saveRDS(error_log, file.path(CONFIG$output_dirs$prepared_data, paste0("error_", model_name, ".rds")))
  })
}

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("Data preparation complete!\n")
cat("Prepared data saved in:", CONFIG$output_dirs$prepared_data, "\n\n")

# Summary
cat("Summary of prepared datasets:\n")
successful_models <- c()
failed_models <- c()

for (model_name in names(model_configs)) {
  file_path <- file.path(CONFIG$output_dirs$prepared_data, paste0("stan_data_", model_name, ".rds"))
  if (file.exists(file_path)) {
    size_mb <- file.info(file_path)$size / 1024 / 1024
    config <- readRDS(file.path(CONFIG$output_dirs$prepared_data, paste0("config_", model_name, ".rds")))
    cat(sprintf("  ✓ %-30s %6.1f MB   N=%d   PFT=%s   Interact=%s   Time=%.1fs\n", 
                model_name, size_mb, config$n_observations, 
                ifelse(config$pft_actually_used, "Y", "N"),
                ifelse(config$include_veg_interactions, "Y", "N"),
                config$data_prep_time))
    successful_models <- c(successful_models, model_name)
  } else {
    cat(sprintf("  ✗ %-30s FAILED\n", model_name))
    failed_models <- c(failed_models, model_name)
  }
}

cat("\n")
cat("Successfully prepared:", length(successful_models), "models\n")
if (length(failed_models) > 0) {
  cat("Failed:", length(failed_models), "models\n")
  cat("Failed models:", paste(failed_models, collapse = ", "), "\n")
}

# Create a summary file
summary_info <- list(
  successful_models = successful_models,
  failed_models = failed_models,
  has_pft_data = has_pft_columns,
  n_observations = nrow(sediment),
  scaling_params = SCALING_PARAMS,
  model_version = "updated_bspline_matern",
  timestamp = Sys.time()
)
saveRDS(summary_info, file.path(CONFIG$output_dirs$prepared_data, "preparation_summary.rds"))
cat("\nPreparation summary saved to:", 
    file.path(CONFIG$output_dirs$prepared_data, "preparation_summary.rds"), "\n")