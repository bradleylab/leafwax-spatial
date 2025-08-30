# 5c_spatial_pattern_diagnostics.R
# Analyzes spatial patterns in model intercepts and slopes
# Uses ONLY fitted model files, no dependency on prepared_data/

library(tidyverse)
library(cmdstanr)
library(posterior)
library(viridis)
library(gstat)
library(sp)
library(ape)
library(patchwork)

# Create output directory
output_dir <- "model_analysis/spatial_pattern_diagnostics"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("\nSPATIAL PATTERN DIAGNOSTICS\n")
cat("===========================\n")

#───────────────────────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────

# Extract data from fitted model (no external dependencies)
extract_model_data <- function(fit) {
  # Get metadata
  vars <- fit$metadata()$variables
  
  # Extract essential data that was passed to Stan
  data_list <- list()
  
  # Basic dimensions
  if ("N" %in% vars) {
    N_draws <- fit$draws("beta_0", format = "matrix")
    data_list$N <- ncol(fit$draws("d2H_wax", format = "matrix"))
  }
  
  # Coordinates (if model has spatial component)
  if ("alpha_spatial" %in% vars || "beta_oipc_spatial" %in% vars) {
    # For spatial models, we need to reconstruct coordinates
    # This is tricky without the original data, but we can try to extract from 
    # the spatial parameters if they're saved
    data_list$has_spatial <- TRUE
    data_list$n_locations <- ncol(fit$draws("alpha_spatial", format = "matrix"))
  } else {
    data_list$has_spatial <- FALSE
  }
  
  # Get scaling parameters if available
  if ("d2H_wax_mean_original" %in% vars) {
    data_list$d2H_mean_original <- mean(fit$draws("d2H_wax_mean_original", format = "matrix"))
    data_list$d2H_sd_original <- mean(fit$draws("d2H_wax_sd_original", format = "matrix"))
  }
  
  # Get observed data if saved
  if ("d2H_wax" %in% vars) {
    data_list$d2H_wax <- fit$draws("d2H_wax", format = "matrix")[1,]  # Just first draw
  }
  
  return(data_list)
}

# Spatial autocorrelation test without coordinates
# Uses pairwise differences to infer spatial structure
test_spatial_structure <- function(values) {
  n <- length(values)
  if (n < 20) return(list(has_structure = NA, message = "Too few points"))
  
  # Simple test: is variance of differences between nearby indices
  # smaller than overall variance?
  # This assumes some spatial ordering in the data
  
  overall_var <- var(values)
  
  # Look at differences between adjacent points (assuming some spatial ordering)
  adj_diffs <- diff(values)
  adj_var <- var(adj_diffs)
  
  # Ratio test
  var_ratio <- adj_var / (2 * overall_var)  # Factor of 2 for difference variance
  
  # If ratio << 1, there's spatial structure
  has_structure <- var_ratio < 0.8
  
  return(list(
    has_structure = has_structure,
    var_ratio = var_ratio,
    message = ifelse(has_structure, 
                     "Likely spatial structure detected",
                     "No clear spatial structure")
  ))
}

#───────────────────────────────────────────────────────────────────────────────
# MAIN ANALYSIS FUNCTION
#───────────────────────────────────────────────────────────────────────────────

analyze_spatial_patterns <- function(model_name) {
  cat("\n\nAnalyzing:", model_name, "\n")
  cat(strrep("-", 60), "\n")
  
  # Load fitted model
  fit_file <- file.path("model_output", model_name, "fit.rds")
  if (!file.exists(fit_file)) {
    cat("  Model file not found - skipping\n")
    return(NULL)
  }
  
  fit <- readRDS(fit_file)
  vars <- fit$metadata()$variables
  
  cat("  Model has", length(vars), "variables\n")
  
  # Show spatial-related variables
  spatial_vars <- vars[grepl("spatial|knot|intercept.*km|slope.*km|sigma_intercept|sigma_slope", vars, ignore.case = TRUE)]
  if (length(spatial_vars) > 0) {
    cat("  Spatial-related variables found:", paste(spatial_vars[1:min(10, length(spatial_vars))], collapse = ", "), "\n")
  }
  
  # Check if this is a spatial model - look for multiple possible names
  has_spatial_intercept <- any(c("alpha_spatial", "intercept_spatial", "spatial_intercept") %in% vars)
  has_spatial_slope <- any(c("beta_oipc_spatial", "slope_spatial", "spatial_slope") %in% vars)
  
  # Also check for GP indicators
  has_gp <- any(c("sigma_intercept_spatial", "sigma_slope_spatial", "ls_intercept_km", "ls_slope_km") %in% vars)
  
  # Check if knot-level parameters exist
  has_knot_params <- any(grepl("knot_intercepts|knot_slopes|z_intercept_spatial|z_slope_spatial", vars))
  
  if (!has_spatial_intercept && !has_spatial_slope && !has_gp && !has_knot_params) {
    cat("  No spatial components found - skipping\n")
    cat("  Available variables:", head(vars, 20), "...\n")
    return(NULL)
  }
  
  # If we have GP params but not the spatial fields, try to find them
  if (has_gp && !has_spatial_intercept) {
    # Look for any intercept-related spatial variable
    intercept_vars <- vars[grepl("intercept.*spatial|spatial.*intercept|alpha_spatial", vars, ignore.case = TRUE)]
    if (length(intercept_vars) > 0) {
      cat("  Found intercept variable:", intercept_vars[1], "\n")
      has_spatial_intercept <- TRUE
    }
  }
  
  if (has_gp && !has_spatial_slope) {
    # Look for any slope-related spatial variable
    slope_vars <- vars[grepl("slope.*spatial|spatial.*slope|beta.*spatial", vars, ignore.case = TRUE)]
    if (length(slope_vars) > 0) {
      cat("  Found slope variable:", slope_vars[1], "\n")
      has_spatial_slope <- TRUE
    }
  }
  
  results <- list(model = model_name)
  
  # Get scaling parameters
  d2H_mean_orig <- 1
  d2H_sd_orig <- 1
  if ("d2H_wax_mean_original" %in% vars) {
    d2H_mean_orig <- mean(fit$draws("d2H_wax_mean_original", format = "matrix"))
    d2H_sd_orig <- mean(fit$draws("d2H_wax_sd_original", format = "matrix"))
  }
  
  #─────────────────────────────────────────────────────────────────────────────
  # 1. ANALYZE SPATIAL INTERCEPTS
  #─────────────────────────────────────────────────────────────────────────────
  
  # First check what spatial variables we actually have
  if (has_gp || has_knot_params) {
    cat("\n1. CHECKING FOR SPATIAL PARAMETERS\n")
    
    # Look for intercept-related parameters
    if ("min_oipc_slope" %in% vars && "max_oipc_slope" %in% vars) {
      # We have summary statistics instead of the full spatial fields
      cat("  Found spatial slope summaries\n")
      
      min_slope <- mean(fit$draws("min_oipc_slope", format = "matrix"))
      max_slope <- mean(fit$draws("max_oipc_slope", format = "matrix"))
      mean_slope <- mean(fit$draws("beta_oipc", format = "matrix"))
      sd_slope <- mean(fit$draws("sd_oipc_slope", format = "matrix"))
      
      cat("  Slope range: [", round(min_slope, 3), ",", round(max_slope, 3), "]\n")
      cat("  Slope SD:", round(sd_slope, 3), "\n")
      
      results$slope <- list(
        range = c(min_slope, max_slope),
        mean = mean_slope,
        sd = sd_slope,
        span = max_slope - min_slope,
        spatial_structure = TRUE  # Assumed if GP model
      )
    }
    
    # Check for intercept at zero
    if ("intercept_original" %in% vars) {
      intercept_orig <- mean(fit$draws("intercept_original", format = "matrix"))
      cat("  Global intercept at OIPC=0:", round(intercept_orig, 1), "‰\n")
      
      # For knot-level analysis, check if we have z-scores
      if ("z_intercept_spatial" %in% vars) {
        cat("  Found knot-level intercept parameters\n")
        z_int_draws <- fit$draws("z_intercept_spatial", format = "matrix")
        z_int_means <- colMeans(z_int_draws)
        
        cat("  Knot z-scores range:", round(range(z_int_means), 3), "\n")
        cat("  Knot z-scores SD:", round(sd(z_int_means), 3), "\n")
      }
    }
    
    # Check length scales and variance components
    if ("ls_intercept_km" %in% vars) {
      ls_int <- mean(fit$draws("ls_intercept_km", format = "matrix"))
      ls_slope <- mean(fit$draws("ls_slope_km", format = "matrix"))
      cat("  Length scales: intercept =", round(ls_int), "km, slope =", round(ls_slope), "km\n")
    }
    
    if ("sigma_intercept_spatial" %in% vars) {
      sigma_int <- mean(fit$draws("sigma_intercept_spatial", format = "matrix"))
      sigma_slope <- mean(fit$draws("sigma_slope_spatial", format = "matrix"))
      cat("  Spatial SDs: intercept =", round(sigma_int, 1), "‰, slope =", round(sigma_slope, 3), "\n")
    }
  }
  
  if (has_spatial_intercept) {
    cat("\n1. SPATIAL INTERCEPT ANALYSIS\n")
    
    # Initialize variables
    alpha_means <- NULL
    alpha_sd <- NULL
    
    # Try to extract intercepts - they might be standardized or in knot space
    if ("alpha_spatial" %in% vars) {
      alpha_draws <- fit$draws("alpha_spatial", format = "matrix")
      alpha_means <- colMeans(alpha_draws)
      alpha_sd <- apply(alpha_draws, 2, sd)
    } else if ("z_intercept_spatial" %in% vars) {
      # We have knot-level z-scores, need to reconstruct
      cat("  Note: Using knot-level z-scores for intercept analysis\n")
      
      # Skip if we can't extract properly
      tryCatch({
        z_int_draws <- fit$draws("z_intercept_spatial", format = "matrix")
        
        # Get sigma_intercept to scale properly
        if ("sigma_intercept_spatial" %in% vars) {
          sigma_int_spatial <- mean(fit$draws("sigma_intercept_spatial", format = "matrix"))
        } else if ("sigma_intercept_raw" %in% vars) {
          sigma_int_raw <- mean(fit$draws("sigma_intercept_raw[1]", format = "matrix"))
          sigma_int_spatial <- sigma_int_raw * d2H_sd_orig
        } else {
          sigma_int_spatial <- 1
        }
        
        # Get global intercept
        if ("intercept_original" %in% vars) {
          global_int <- mean(fit$draws("intercept_original", format = "matrix"))
        } else if ("beta_0" %in% vars) {
          beta_0 <- mean(fit$draws("beta_0", format = "matrix"))
          global_int <- beta_0 * d2H_sd_orig + d2H_mean_orig
        } else {
          global_int <- 0
        }
        
        # These are knot-level values, not observation-level
        alpha_means <- colMeans(z_int_draws) * sigma_int_spatial + global_int
        alpha_sd <- apply(z_int_draws, 2, sd) * sigma_int_spatial
        
        cat("  Analyzing", ncol(z_int_draws), "knot intercepts\n")
      }, error = function(e) {
        cat("  Could not extract intercept spatial field\n")
      })
    }
    
    # Check if we successfully extracted intercepts
    if (!is.null(alpha_means)) {
      # Now alpha_means should be in original scale
      alpha_means_orig <- alpha_means
      alpha_sd_orig <- alpha_sd
    } else {
      cat("  Skipping detailed intercept analysis\n")
      has_spatial_intercept <- FALSE
    }
    
    # Basic statistics
    cat("  Range:", round(range(alpha_means_orig), 1), "‰\n")
    cat("  Mean:", round(mean(alpha_means_orig), 1), "‰\n")
    cat("  SD:", round(sd(alpha_means_orig), 1), "‰\n")
    cat("  Total span:", round(diff(range(alpha_means_orig)), 1), "‰\n")
    
    # Test for spatial structure
    spatial_test <- test_spatial_structure(alpha_means_orig)
    cat("  Spatial structure test:", spatial_test$message, "\n")
    if (!is.na(spatial_test$var_ratio)) {
      cat("    Variance ratio:", round(spatial_test$var_ratio, 3), "\n")
    }
    
    # Identify extreme values
    extreme_threshold <- 2  # SDs from mean
    z_scores <- (alpha_means_orig - mean(alpha_means_orig)) / sd(alpha_means_orig)
    n_extreme <- sum(abs(z_scores) > extreme_threshold)
    cat("  Extreme values (>", extreme_threshold, "SD):", n_extreme, 
        "(", round(100*n_extreme/length(z_scores), 1), "%)\n")
    
    # Check for bimodality (might indicate discrete groups)
    density_est <- density(alpha_means_orig)
    peaks <- which(diff(sign(diff(density_est$y))) == -2) + 1
    cat("  Number of modes in distribution:", length(peaks), "\n")
    
    # Store results
    results$intercept <- list(
      range = range(alpha_means_orig),
      mean = mean(alpha_means_orig),
      sd = sd(alpha_means_orig),
      span = diff(range(alpha_means_orig)),
      n_extreme = n_extreme,
      pct_extreme = 100*n_extreme/length(z_scores),
      spatial_structure = spatial_test$has_structure,
      n_modes = length(peaks),
      values = alpha_means_orig  # Store for plotting
    )
    
    # Create diagnostic plot
    p1 <- ggplot(data.frame(intercept = alpha_means_orig), aes(x = intercept)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
      geom_density(color = "red", size = 1) +
      labs(title = paste("Intercept Distribution:", model_name),
           x = "Intercept at OIPC=0 (‰)", y = "Density") +
      theme_minimal()
    
    # Uncertainty vs value plot
    p2 <- ggplot(data.frame(mean = alpha_means_orig, sd = alpha_sd_orig), 
                 aes(x = mean, y = sd)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE) +
      labs(title = "Intercept Uncertainty",
           x = "Intercept (‰)", y = "Posterior SD (‰)") +
      theme_minimal()
    
    # Adjacent differences plot (to check smoothness)
    if (length(alpha_means_orig) > 20) {
      diffs <- diff(sort(alpha_means_orig))
      p3 <- ggplot(data.frame(x = 1:length(diffs), diff = diffs), aes(x = x, y = diff)) +
        geom_point(alpha = 0.6) +
        geom_hline(yintercept = median(diffs), color = "red", linetype = "dashed") +
        labs(title = "Ordered Adjacent Differences",
             x = "Index", y = "Difference (‰)") +
        theme_minimal()
    } else {
      p3 <- NULL
    }
  }
  
  #─────────────────────────────────────────────────────────────────────────────
  # 2. ANALYZE SPATIAL SLOPES
  #─────────────────────────────────────────────────────────────────────────────
  
  if (has_spatial_slope) {
    cat("\n2. SPATIAL SLOPE ANALYSIS\n")
    
    # Try to extract slopes
    if ("beta_oipc_spatial" %in% vars) {
      beta_draws <- fit$draws("beta_oipc_spatial", format = "matrix")
      beta_means <- colMeans(beta_draws)
      beta_sd <- apply(beta_draws, 2, sd)
    } else if ("z_slope_spatial" %in% vars) {
      # We have knot-level z-scores
      cat("  Note: Using knot-level z-scores for slope analysis\n")
      z_slope_draws <- fit$draws("z_slope_spatial", format = "matrix")
      
      # Get sigma_slope to scale properly
      if ("sigma_slope_spatial" %in% vars) {
        sigma_slope_spatial <- mean(fit$draws("sigma_slope_spatial", format = "matrix"))
      } else if ("sigma_slope_raw" %in% vars) {
        sigma_slope_spatial <- mean(fit$draws("sigma_slope_raw[1]", format = "matrix"))
      } else {
        sigma_slope_spatial <- 0.1  # reasonable default
      }
      
      # Get global slope
      if ("beta_oipc" %in% vars) {
        global_slope <- mean(fit$draws("beta_oipc", format = "matrix"))
      } else {
        global_slope <- 0.8  # reasonable default
      }
      
      # These are knot-level values
      beta_means <- colMeans(z_slope_draws) * sigma_slope_spatial + global_slope
      beta_sd <- apply(z_slope_draws, 2, sd) * sigma_slope_spatial
      
      cat("  Analyzing", ncol(z_slope_draws), "knot slopes\n")
    }
    
    # Basic statistics
    cat("  Range:", round(range(beta_means), 3), "\n")
    cat("  Mean:", round(mean(beta_means), 3), "\n")
    cat("  SD:", round(sd(beta_means), 3), "\n")
    cat("  Total span:", round(diff(range(beta_means)), 3), "\n")
    
    # Test for spatial structure
    spatial_test_slope <- test_spatial_structure(beta_means)
    cat("  Spatial structure test:", spatial_test_slope$message, "\n")
    
    # Check for physically implausible values
    n_negative <- sum(beta_means < 0)
    n_too_large <- sum(beta_means > 1.5)
    cat("  Negative slopes:", n_negative, "\n")
    cat("  Slopes > 1.5:", n_too_large, "\n")
    
    # Coefficient of variation
    cv <- sd(beta_means) / mean(beta_means)
    cat("  Coefficient of variation:", round(cv, 3), "\n")
    
    # Store results
    results$slope <- list(
      range = range(beta_means),
      mean = mean(beta_means),
      sd = sd(beta_means),
      span = diff(range(beta_means)),
      n_negative = n_negative,
      n_too_large = n_too_large,
      cv = cv,
      spatial_structure = spatial_test_slope$has_structure,
      values = beta_means  # Store for plotting
    )
    
    # Create diagnostic plots for slopes
    p4 <- ggplot(data.frame(slope = beta_means), aes(x = slope)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgreen", color = "black") +
      geom_density(color = "red", size = 1) +
      geom_vline(xintercept = c(0, 1), linetype = "dashed", alpha = 0.5) +
      labs(title = paste("Slope Distribution:", model_name),
           x = "OIPC Slope", y = "Density") +
      theme_minimal()
    
    p5 <- ggplot(data.frame(mean = beta_means, sd = beta_sd), 
                 aes(x = mean, y = sd)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "loess", se = TRUE) +
      labs(title = "Slope Uncertainty",
           x = "Slope", y = "Posterior SD") +
      theme_minimal()
  }
  
  #─────────────────────────────────────────────────────────────────────────────
  # 3. RELATIONSHIP BETWEEN SLOPES AND INTERCEPTS
  #─────────────────────────────────────────────────────────────────────────────
  
  if (has_spatial_intercept && has_spatial_slope) {
    cat("\n3. SLOPE-INTERCEPT RELATIONSHIP\n")
    
    # Check correlation
    cor_value <- cor(alpha_means_orig, beta_means)
    cat("  Correlation:", round(cor_value, 3), "\n")
    
    # Test if correlation is significant
    cor_test <- cor.test(alpha_means_orig, beta_means)
    cat("  P-value:", format.pval(cor_test$p.value), "\n")
    
    results$slope_intercept_cor <- cor_value
    results$slope_intercept_p <- cor_test$p.value
    
    # Create scatter plot
    p6 <- ggplot(data.frame(intercept = alpha_means_orig, slope = beta_means),
                 aes(x = intercept, y = slope)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = TRUE) +
      labs(title = paste("Slope vs Intercept:", model_name),
           x = "Intercept (‰)", y = "Slope") +
      theme_minimal()
  }
  
  #─────────────────────────────────────────────────────────────────────────────
  # 4. SAVE PLOTS
  #─────────────────────────────────────────────────────────────────────────────
  
  # Combine plots
  if (has_spatial_intercept && has_spatial_slope) {
    if (!is.null(p3)) {
      combined_plot <- (p1 + p2 + p3) / (p4 + p5 + p6)
    } else {
      combined_plot <- (p1 + p2) / (p4 + p5 + p6)
    }
  } else if (has_spatial_intercept) {
    combined_plot <- p1 + p2 + p3
  } else {
    combined_plot <- p4 + p5
  }
  
  ggsave(file.path(output_dir, paste0("spatial_diagnostics_", model_name, ".png")),
         combined_plot, width = 15, height = 10, dpi = 300)
  
  #─────────────────────────────────────────────────────────────────────────────
  # 5. VARIANCE DECOMPOSITION (if available)
  #─────────────────────────────────────────────────────────────────────────────
  
  if ("prop_variance_spatial" %in% vars) {
    cat("\n4. VARIANCE DECOMPOSITION\n")
    
    var_spatial <- mean(fit$draws("prop_variance_spatial", format = "matrix"))
    var_residual <- mean(fit$draws("prop_variance_residual", format = "matrix"))
    
    cat("  Spatial variance:", round(100*var_spatial, 1), "%\n")
    cat("  Residual variance:", round(100*var_residual, 1), "%\n")
    
    results$var_spatial <- var_spatial
    results$var_residual <- var_residual
  }
  
  return(results)
}

#───────────────────────────────────────────────────────────────────────────────
# RUN ANALYSIS FOR ALL MODELS
#───────────────────────────────────────────────────────────────────────────────

# Find all fitted models
model_files <- list.files("model_output", pattern = "fit.rds", 
                         recursive = TRUE, full.names = TRUE)
model_names <- basename(dirname(model_files))

cat("Found", length(model_names), "fitted models\n")

# Analyze each model
all_results <- list()
for (model in model_names) {
  tryCatch({
    results <- analyze_spatial_patterns(model)
    if (!is.null(results)) {
      all_results[[model]] <- results
    }
  }, error = function(e) {
    cat("  Error analyzing", model, ":", e$message, "\n")
  })
}

#───────────────────────────────────────────────────────────────────────────────
# SUMMARY REPORT
#───────────────────────────────────────────────────────────────────────────────

cat("\n\n", strrep("=", 70), "\n")
cat("SPATIAL PATTERN SUMMARY\n")
cat(strrep("=", 70), "\n")

# Create summary dataframe
summary_df <- map_df(all_results, function(res) {
  df <- data.frame(
    model = res$model,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(res$intercept)) {
    df$intercept_span <- res$intercept$span
    df$intercept_sd <- res$intercept$sd
    df$intercept_extreme_pct <- res$intercept$pct_extreme
    df$intercept_spatial <- res$intercept$spatial_structure
    df$intercept_modes <- res$intercept$n_modes
  }
  
  if (!is.null(res$slope)) {
    df$slope_span <- res$slope$span
    df$slope_sd <- res$slope$sd
    df$slope_cv <- res$slope$cv
    df$slope_negative <- res$slope$n_negative
    df$slope_spatial <- res$slope$spatial_structure
  }
  
  if (!is.null(res$slope_intercept_cor)) {
    df$slope_int_cor <- res$slope_intercept_cor
    df$slope_int_p <- res$slope_intercept_p
  }
  
  if (!is.null(res$var_spatial)) {
    df$var_spatial_pct <- round(100 * res$var_spatial, 1)
  }
  
  return(df)
})

# Sort by intercept span (as a measure of spatial variation)
if ("intercept_span" %in% names(summary_df)) {
  summary_df <- summary_df %>% arrange(desc(intercept_span))
}

cat("\nSummary of spatial patterns:\n")
if (nrow(summary_df) > 0) {
  print(as.data.frame(summary_df))
} else {
  cat("No spatial models found\n")
}

# Save summary
write.csv(summary_df, file.path(output_dir, "spatial_pattern_summary.csv"), 
          row.names = FALSE)

# Identify potentially problematic models
cat("\n\nPOTENTIAL ISSUES:\n")
cat(strrep("-", 40), "\n")

# Models with very large intercept spans
if ("intercept_span" %in% names(summary_df)) {
  large_span <- summary_df %>% filter(intercept_span > 150)
  if (nrow(large_span) > 0) {
    cat("\nModels with intercept span > 150‰:\n")
    print(large_span$model)
  }
}

# Models with many extreme values
if ("intercept_extreme_pct" %in% names(summary_df)) {
  high_extreme <- summary_df %>% filter(intercept_extreme_pct > 10)
  if (nrow(high_extreme) > 0) {
    cat("\nModels with >10% extreme intercepts:\n")
    print(high_extreme$model)
  }
}

# Models with negative slopes
if ("slope_negative" %in% names(summary_df)) {
  neg_slopes <- summary_df %>% filter(slope_negative > 0)
  if (nrow(neg_slopes) > 0) {
    cat("\nModels with negative slopes:\n")
    print(neg_slopes$model)
  }
}

# Models with very large slope spans
if ("slope_span" %in% names(summary_df)) {
  large_slope_span <- summary_df %>% filter(slope_span > 0.5)
  if (nrow(large_slope_span) > 0) {
    cat("\nModels with slope span > 0.5:\n")
    print(large_slope_span$model)
  }
}

# Models with high slope CV
if ("slope_cv" %in% names(summary_df)) {
  high_cv <- summary_df %>% filter(slope_cv > 0.2)
  if (nrow(high_cv) > 0) {
    cat("\nModels with slope CV > 0.2:\n")
    print(high_cv$model)
  }
}

# Models with no spatial structure
if ("intercept_spatial" %in% names(summary_df) || "slope_spatial" %in% names(summary_df)) {
  no_spatial <- summary_df
  
  # Check intercept spatial structure if column exists
  if ("intercept_spatial" %in% names(summary_df)) {
    no_spatial <- no_spatial %>% 
      filter(is.na(intercept_spatial) | intercept_spatial == TRUE)
  }
  
  # Check slope spatial structure if column exists  
  if ("slope_spatial" %in% names(summary_df)) {
    no_spatial_slope <- summary_df %>%
      filter(!is.na(slope_spatial) & !slope_spatial)
    
    if (nrow(no_spatial_slope) > 0) {
      cat("\nModels with no detected spatial structure in slopes:\n")
      print(no_spatial_slope$model)
    }
  }
  
  # Check intercept spatial structure
  if ("intercept_spatial" %in% names(summary_df)) {
    no_spatial_int <- summary_df %>%
      filter(!is.na(intercept_spatial) & !intercept_spatial)
    
    if (nrow(no_spatial_int) > 0) {
      cat("\nModels with no detected spatial structure in intercepts:\n")
      print(no_spatial_int$model)
    }
  }
}

cat("\n\nDiagnostics complete! Results saved to:", output_dir, "\n")

# Create final summary plot comparing all models
if (length(all_results) > 1) {
  # Extract all intercept ranges
  intercept_data <- map_df(all_results, function(res) {
    if (!is.null(res$intercept)) {
      data.frame(
        model = res$model,
        value = res$intercept$values,
        type = "Intercept"
      )
    }
  })
  
  slope_data <- map_df(all_results, function(res) {
    if (!is.null(res$slope)) {
      data.frame(
        model = res$model,
        value = res$slope$values,
        type = "Slope"
      )
    }
  })
  
  if (nrow(intercept_data) > 0 || nrow(slope_data) > 0) {
    combined_data <- bind_rows(intercept_data, slope_data)
    
    p_compare <- ggplot(combined_data, aes(x = model, y = value)) +
      geom_boxplot(aes(fill = model)) +
      facet_wrap(~type, scales = "free_y", ncol = 1) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none") +
      labs(title = "Spatial Parameter Distributions Across Models",
           x = "", y = "Value")
    
    ggsave(file.path(output_dir, "model_comparison_spatial.png"),
           p_compare, width = 12, height = 8, dpi = 300)
  }
}

# Save all results for future analysis
saveRDS(all_results, file.path(output_dir, "spatial_diagnostics_full.rds"))