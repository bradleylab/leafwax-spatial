#───────────────────────────────────────────────────────────────────────────────
# 0_load_config.R
#
# Load project configuration from YAML file
# Provides centralized configuration management for all scripts
# Source this at the beginning of each script to access config settings
#
# Input: config.yaml (project configuration file)
# Output: config list object with project settings
#───────────────────────────────────────────────────────────────────────────────

library(yaml)

load_config <- function(config_file = "config.yaml") {
  if (!file.exists(config_file)) {
    stop("Configuration file not found: ", config_file, "\n",
         "Please create config.yaml from the template.")
  }
  
  config <- read_yaml(config_file)
  cat("Loaded configuration from", config_file, "\n")
  
  # Validate essential parameters exist
  required_params <- c("spatial_scales", "model_configs", 
                      "input_data", "output_dirs")
  missing <- setdiff(required_params, names(config))
  if (length(missing) > 0) {
    stop("Missing required config parameters: ", paste(missing, collapse = ", "))
  }
  
  return(config)
}

# Load config globally
CONFIG <- load_config()
config <- CONFIG  # Create lowercase alias for compatibility


# Helper function to get model config
get_model_config <- function(model_name) {
  if (!model_name %in% names(CONFIG$model_configs)) {
    stop("Model '", model_name, "' not found in config.yaml")
  }
  return(CONFIG$model_configs[[model_name]])
}