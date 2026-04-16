#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3b_pre-flight_check.R - Pre-flight checks before Stan modeling
#───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
})

cat("PRE-FLIGHT CHECKS FOR STAN MODELING\n")
cat("===================================\n\n")

# Track any issues
issues <- list()

# 1. Check required files
cat("1. Checking required files...\n")
required_files <- c(
  "results/3_sediment_ready_for_modeling.rds",
  "4d_leaf_wax_spatial_model.stan"
)

for (file in required_files) {
  if (file.exists(file)) {
    cat("  ✓", file, "\n")
  } else {
    cat("  ❌", file, "NOT FOUND\n")
    issues$files <- c(issues$files, file)
  }
}

# 2. Check data integrity
cat("\n2. Checking data integrity...\n")
if (file.exists("results/3_sediment_ready_for_modeling.rds")) {
  sediment <- readRDS("results/3_sediment_ready_for_modeling.rds")
  
  # Basic checks
  cat("  Records:", nrow(sediment), "\n")
  
  # Required columns - UPDATED to include TerraClimate
  required_cols <- c("d2H_wax", "d2H_wax_err", "longitude", "latitude",
                     "c4_values_filled", "c4_distances", 
                     "oipc_values", "oipc_distances", "oipc_se_values",
                     "elevation_values", "elevation_distances",
                     "pft_distances",
                     # TerraClimate columns
                     "tc_ppt_values", "tc_ppt_distances",
                     "tc_soil_values", "tc_soil_distances",
                     "tc_tmax_values", "tc_tmax_distances",
                     "tc_vpd_values", "tc_vpd_distances")
  
  missing_cols <- required_cols[!required_cols %in% names(sediment)]
  if (length(missing_cols) > 0) {
    cat("  ❌ Missing columns:", paste(missing_cols, collapse = ", "), "\n")
    issues$columns <- missing_cols
  } else {
    cat("  ✓ All required columns present\n")
  }
  
  # Check for NA values in key columns - UPDATED
  na_check <- c("d2H_wax", "d2H_wax_err", "longitude", "latitude",
                "annual_precip", "soil_moisture", "max_temp", "vpd")
  for (col in na_check) {
    if (col %in% names(sediment)) {
      n_na <- sum(is.na(sediment[[col]]))
      if (n_na > 0) {
        cat("  ⚠ ", col, "has", n_na, "NA values\n")
      }
    }
  }
  
  # Check pixel arrays - UPDATED to include TerraClimate
  cat("\n  Pixel array summary:\n")
  pixel_cols <- c("c4_values_filled", "oipc_values", "elevation_values",
                  "tc_ppt_values", "tc_soil_values", "tc_tmax_values", "tc_vpd_values")
  
  for (col in pixel_cols) {
    if (col %in% names(sediment)) {
      lengths <- map_int(sediment[[col]], length)
      cat("    ", col, ": total =", sum(lengths), 
          ", mean =", round(mean(lengths), 1),
          ", range =", min(lengths), "-", max(lengths), "\n")
    }
  }
  
  # Check TerraClimate summary variables
  cat("\n  TerraClimate summary statistics:\n")
  if (all(c("annual_precip", "soil_moisture", "max_temp", "vpd") %in% names(sediment))) {
    cat("    Annual precip: ", round(mean(sediment$annual_precip, na.rm = TRUE), 0), 
        " mm (", round(min(sediment$annual_precip, na.rm = TRUE), 0), 
        "-", round(max(sediment$annual_precip, na.rm = TRUE), 0), ")\n", sep = "")
    cat("    Soil moisture: ", round(mean(sediment$soil_moisture, na.rm = TRUE), 0), 
        " mm (", round(min(sediment$soil_moisture, na.rm = TRUE), 0), 
        "-", round(max(sediment$soil_moisture, na.rm = TRUE), 0), ")\n", sep = "")
    cat("    Max temp: ", round(mean(sediment$max_temp, na.rm = TRUE), 1), 
        " °C (", round(min(sediment$max_temp, na.rm = TRUE), 1), 
        "-", round(max(sediment$max_temp, na.rm = TRUE), 1), ")\n", sep = "")
    cat("    VPD: ", round(mean(sediment$vpd, na.rm = TRUE), 2), 
        " kPa (", round(min(sediment$vpd, na.rm = TRUE), 2), 
        "-", round(max(sediment$vpd, na.rm = TRUE), 2), ")\n", sep = "")
  }
}

# 3. Check Stan installation
cat("\n3. Checking Stan installation...\n")
tryCatch({
  cmdstan_path <- cmdstan_path()
  cat("  ✓ CmdStan path:", cmdstan_path, "\n")
  
  version <- cmdstan_version()
  cat("  ✓ CmdStan version:", version, "\n")
}, error = function(e) {
  cat("  ❌ CmdStan not properly installed\n")
  issues$stan <- "CmdStan installation"
})

# 4. Test minimal Stan compilation
cat("\n4. Testing minimal Stan compilation...\n")
if (!any(grepl("stan", names(issues)))) {
  # Create a minimal test
  test_stan <- "
data {
  int<lower=0> N;
  vector[N] y;
}
parameters {
  real mu;
  real<lower=0> sigma;
}
model {
  y ~ normal(mu, sigma);
}
"

  writeLines(test_stan, "test_minimal.stan")
  
  tryCatch({
    test_model <- cmdstan_model("test_minimal.stan", quiet = TRUE)
    cat("  ✓ Basic Stan compilation works\n")
    
    # Test with minimal data
    test_data <- list(N = 10, y = rnorm(10))
    test_fit <- test_model$sample(
      data = test_data,
      chains = 1,
      iter_sampling = 100,
      iter_warmup = 100,
      refresh = 0,
      show_messages = FALSE
    )
    cat("  ✓ Basic Stan sampling works\n")
    
  }, error = function(e) {
    cat("  ❌ Stan compilation/sampling failed:", e$message, "\n")
    issues$stan_compile <- e$message
  })
  
  # Clean up
  if (file.exists("test_minimal.stan")) file.remove("test_minimal.stan")
  if (file.exists("test_minimal")) file.remove("test_minimal")
  # Remove any compiled model files
  temp_files <- list.files(pattern = "test_minimal.*", full.names = TRUE)
  if (length(temp_files) > 0) file.remove(temp_files)
}

# 5. Memory check
cat("\n5. System resources...\n")
if (Sys.info()["sysname"] == "Linux") {
  mem_info <- system("free -h", intern = TRUE)
  cat(paste("  ", mem_info[1:2], collapse = "\n"), "\n")
} else {
  cat("  Memory check only available on Linux\n")
}

# Check R memory limit
cat("  R memory limit:", 
    if (is.infinite(memory.limit())) "unlimited" else paste(memory.limit(), "MB"), 
    "\n")

# 6. Summary
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
if (length(issues) == 0) {
  cat("✅ ALL CHECKS PASSED - Ready for modeling\n")
  cat("\nNext steps:\n")
  cat("  1. Run your Stan model with the updated data\n")
  cat("  2. The data now includes climate covariates:\n")
  cat("     - annual_precip (mm/year)\n")
  cat("     - soil_moisture (mm)\n")
  cat("     - max_temp (°C)\n")
  cat("     - vpd (kPa)\n")
  quit(status = 0)
} else {
  cat("❌ ISSUES FOUND:\n")
  for (issue_type in names(issues)) {
    cat("  -", issue_type, ":", paste(issues[[issue_type]], collapse = ", "), "\n")
  }
  cat("\nPlease fix these issues before proceeding with modeling.\n")
  quit(status = 1)
}