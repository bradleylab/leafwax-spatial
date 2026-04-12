#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# spatial_confounding_check.R
# Quick diagnostic: Check collinearity between spatial effects and δ²H_precip
# This addresses the spatial confounding concern raised in the review
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(cmdstanr)
library(posterior)
library(patchwork)  # For combining plots with / operator

cat("\nSPATIAL CONFOUNDING DIAGNOSTIC (Quick Check)\n")
cat(strrep("=", 50), "\n\n")

# Setup directories
model_dir <- "model_output"
data_dir <- "prepared_data_consolidated"
output_dir <- "model_analysis/spatial_confounding"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get all spatial models (those ending in _sp)
all_models <- list.dirs(model_dir, full.names = FALSE, recursive = FALSE)
spatial_models <- all_models[grepl("_sp$", all_models)]

cat("Found", length(spatial_models), "spatial models:\n")
for (m in spatial_models) cat("  •", m, "\n")
cat("\n")

# Storage for results
results_list <- list()

#═══════════════════════════════════════════════════════════════════════════════
# DIAGNOSTIC 1: Correlation between spatial effect and δ²H_precip
#═══════════════════════════════════════════════════════════════════════════════

for (model_name in spatial_models) {
  
  cat("\n", strrep("═", 70), "\n")
  cat("ANALYZING:", model_name, "\n")
  cat(strrep("═", 70), "\n")
  
  # Load model fit and data
  fit_file <- file.path(model_dir, model_name, "fit.rds")
  data_file <- file.path(data_dir, paste0("stan_data_", model_name, ".rds"))
  
  if (!file.exists(fit_file) | !file.exists(data_file)) {
    cat("  ⚠ Files not found, skipping\n")
    next
  }
  
  fit <- readRDS(fit_file)
  stan_data <- readRDS(data_file)
  
  # Get available variables - use stan_variables not variables!
  all_vars <- fit$metadata()$stan_variables
  
  # ─────────────────────────────────────────────────────────────────────────
  # Extract spatial effect using RESIDUALS approach
  # The spatial effect is what remains after accounting for fixed effects
  # ─────────────────────────────────────────────────────────────────────────
  
  n_sites <- length(stan_data$d2H_wax)
  
  # Strategy: spatial effect = observed - (fixed effects prediction)
  # This captures what the GP is actually explaining
  
  # Get posterior predictions - using mu (fitted values without observation noise)
  if ("mu" %in% all_vars) {
    mu_draws <- fit$draws("mu", format = "matrix")
    y_pred_full <- colMeans(mu_draws)  # Full model prediction (with spatial)
    cat("  ✓ Found mu (fitted values)\n")
  } else if ("d2H_rep" %in% all_vars) {
    # Fallback to d2H_rep if mu not available
    y_rep <- fit$draws("d2H_rep", format = "matrix")
    y_pred_full <- colMeans(y_rep)  # Full model prediction (with spatial)
    cat("  ✓ Found d2H_rep (using as fallback)\n")
  } else {
    cat("  ⚠ Neither mu nor d2H_rep found\n")
    next
  }

  # Get spatial effect directly from alpha_spatial
  spatial_effect <- NULL
  spatial_var_name <- "derived_from_residuals"

  if ("alpha_spatial" %in% all_vars) {
    # Direct extraction of spatial effects
    alpha_spatial_draws <- fit$draws("alpha_spatial", format = "matrix")
    # Take column means to get posterior mean at each site
    spatial_effect <- colMeans(alpha_spatial_draws)
    spatial_var_name <- "alpha_spatial (direct GP contribution)"
    cat("  ✓ Extracted spatial effect from alpha_spatial\n")
    cat(sprintf("    Dimensions: %d draws x %d sites -> %d posterior means\n",
                nrow(alpha_spatial_draws), ncol(alpha_spatial_draws), length(spatial_effect)))
    
  } else {
    # Alternative: compute spatial effect from mu and non-spatial prediction
    cat("  ⚠ No alpha_spatial found, trying alternative methods\n")

    # Check if we can compute it from mu and a baseline model
    if ("mu" %in% all_vars) {
      # The spatial effect could be approximated by comparing to a non-spatial baseline
      # For now, we'll use the residuals as a proxy
      y_obs <- stan_data$d2H_wax

      # Compute model residuals
      residual_full <- y_obs - y_pred_full

      cat("  Computing approximate spatial effect from residuals\n")
      cat("  Note: This is less precise than using alpha_spatial directly\n")

      # The spatial effect is approximately the systematic pattern in residuals
      # This is imperfect but gives us something to work with
      spatial_effect <- residual_full
      spatial_var_name <- "approximate_from_residuals"
    } else {
      cat("  ⚠ Cannot compute spatial effect - missing required variables\n")
      next
    }
  }
  
  if (is.null(spatial_effect)) {
    cat("  ⚠ Could not extract spatial effect\n")
    cat("  Available variables:\n")
    print(head(all_vars, 30))
    next
  }
  
  if (length(spatial_effect) != n_sites) {
    cat("\n  ⚠ Spatial effect has", length(spatial_effect), 
        "values but", n_sites, "sites\n")
    cat("     Skipping detailed analysis\n")
    
    results_list[[model_name]] <- list(
      model = model_name,
      spatial_var = spatial_var_name,
      n_spatial = length(spatial_effect),
      n_sites = n_sites,
      note = "Dimension mismatch"
    )
    next
  }
  
  cat("✓ Extracted spatial effect with", n_sites, "values\n")
  cat("  Method:", spatial_var_name, "\n")
  
  # ─────────────────────────────────────────────────────────────────────────
  # Get δ²H_precip (OIPC) values
  # ─────────────────────────────────────────────────────────────────────────

  # Try different possible names for precipitation isotopes
  d2H_precip <- NULL
  if ("d2H_precip" %in% names(stan_data)) {
    d2H_precip <- stan_data$d2H_precip
  } else if ("oipc_values" %in% names(stan_data)) {
    # oipc_values is a matrix (sites x distance_scales)
    # Use first column (1 km scale) for site-level values
    if (is.matrix(stan_data$oipc_values)) {
      d2H_precip <- stan_data$oipc_values[, 1]  # First column = 1 km scale
      cat("  Using oipc_values[, 1] (1 km scale) as d2H_precip\n")
      cat(sprintf("    Extracted %d site-level values from %d x %d matrix\n",
                  length(d2H_precip), nrow(stan_data$oipc_values), ncol(stan_data$oipc_values)))
    } else {
      d2H_precip <- stan_data$oipc_values
      cat("  Using oipc_values as d2H_precip\n")
    }
  } else if ("OIPC" %in% names(stan_data)) {
    d2H_precip <- stan_data$OIPC
  } else if ("oipc" %in% names(stan_data)) {
    d2H_precip <- stan_data$oipc
  }

  if (is.null(d2H_precip)) {
    cat("  ⚠ Could not find precipitation isotope values\n")
    cat("  Available variables in stan_data:\n")
    cat("    ", paste(names(stan_data)[1:min(20, length(names(stan_data)))], collapse = ", "), "\n")
    next
  }

  # Add dimension check
  cat(sprintf("  Dimensions check: spatial_effect = %d, d2H_precip = %d\n",
              length(spatial_effect), length(d2H_precip)))
  
  # ─────────────────────────────────────────────────────────────────────────
  # Compute correlation
  # ─────────────────────────────────────────────────────────────────────────
  
  cor_result <- cor.test(spatial_effect, d2H_precip, method = "pearson")
  cor_value <- cor_result$estimate
  p_value <- cor_result$p.value
  
  cat("\nCOLLINEARITY RESULTS\n")
  cat(sprintf("  Correlation (r):       %7.3f\n", cor_value))
  cat(sprintf("  P-value:               %.3e\n", p_value))
  cat(sprintf("  R² (shared variance):  %7.3f\n", cor_value^2))
  
  # Interpretation
  abs_cor <- abs(cor_value)
  if (abs_cor > 0.7) {
    cat("  HIGH CONFOUNDING RISK\n")
    cat("    Spatial term strongly correlated with predictor\n")
  } else if (abs_cor > 0.5) {
    cat("  MODERATE CONFOUNDING RISK\n")
    cat("    Consider restricted spatial regression\n")
  } else if (abs_cor > 0.3) {
    cat("  LOW-MODERATE CORRELATION\n")
    cat("    Some shared variance but likely acceptable\n")
  } else {
    cat("  LOW CONFOUNDING RISK\n")
    cat("    Spatial term is largely independent of predictor\n")
  }

  # Determine risk level (outside plot creation for use in results)
  risk_text <- if (abs_cor > 0.7) "HIGH RISK" else if (abs_cor > 0.5) "MODERATE RISK" else "LOW RISK"
  risk_color <- if (abs_cor > 0.7) "red" else if (abs_cor > 0.5) "orange" else "darkgreen"

  # ─────────────────────────────────────────────────────────────────────────
  # Create diagnostic plot (with error handling)
  # ─────────────────────────────────────────────────────────────────────────

  tryCatch({
    df_plot <- data.frame(
      d2H_precip = d2H_precip,
      spatial_effect = spatial_effect,
      longitude = stan_data$longitude,
      latitude = stan_data$latitude
    )

    # Scatter plot with regression line
    p1 <- ggplot(df_plot, aes(x = d2H_precip, y = spatial_effect)) +
      geom_point(alpha = 0.5, size = 2) +
      geom_smooth(method = "lm", color = "red", se = TRUE) +
      labs(
        title = paste(model_name, "- Spatial Confounding Check"),
        subtitle = sprintf("Correlation: r = %.3f, p = %.3e", cor_value, p_value),
        x = "δ²H_precip (‰)",
        y = "Spatial Effect (‰)"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))

    p1 <- p1 + annotate("text", x = Inf, y = Inf,
                       label = risk_text,
                       hjust = 1.1, vjust = 1.5,
                       size = 6, fontface = "bold", color = risk_color)

    ggsave(
      file.path(output_dir, paste0(model_name, "_confounding_scatter.png")),
      p1, width = 8, height = 6, dpi = 300
    )
    cat("  ✓ Saved correlation scatter plot\n")

  }, error = function(e) {
    cat("  ⚠ Warning: Could not create scatter plot:", e$message, "\n")
  })

  # Create data frame for spatial plots (outside tryCatch so it's available for both plots)
  df_plot <- data.frame(
    d2H_precip = d2H_precip,
    spatial_effect = spatial_effect,
    longitude = stan_data$longitude,
    latitude = stan_data$latitude
  )

  # Spatial map showing both variables (with error handling)
  p2 <- NULL
  p3 <- NULL

  tryCatch({
    p2 <- ggplot(df_plot, aes(x = longitude, y = latitude)) +
      geom_point(aes(color = spatial_effect), size = 2.5, alpha = 0.7) +
      scale_color_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, name = "Spatial\nEffect (‰)"
      ) +
      labs(title = paste(model_name, "- Spatial Effect Pattern")) +
      theme_minimal() +
      coord_quickmap()
  }, error = function(e) {
    cat("  ⚠ Warning: Could not create spatial effect plot:", e$message, "\n")
  })

  tryCatch({
    p3 <- ggplot(df_plot, aes(x = longitude, y = latitude)) +
      geom_point(aes(color = d2H_precip), size = 2.5, alpha = 0.7) +
      scale_color_viridis_c(name = "δ²H_precip\n(‰)") +
      labs(title = paste(model_name, "- δ²H_precip Pattern")) +
      theme_minimal() +
      coord_quickmap()
  }, error = function(e) {
    cat("  ⚠ Warning: Could not create d2H_precip plot:", e$message, "\n")
  })

  # Combine plots with error handling
  if (!is.null(p2) && !is.null(p3)) {
    tryCatch({
      # Try to combine plots
      p_combined <- p2 / p3

      ggsave(
        file.path(output_dir, paste0(model_name, "_spatial_patterns.png")),
        p_combined, width = 10, height = 10, dpi = 300
      )
      cat("  ✓ Saved combined spatial pattern plot\n")

    }, error = function(e) {
      cat("  ⚠ Warning: Could not combine spatial plots:", e$message, "\n")

      # Try to save individual plots instead
      if (!is.null(p2)) {
        tryCatch({
          ggsave(
            file.path(output_dir, paste0(model_name, "_spatial_effect_map.png")),
            p2, width = 8, height = 6, dpi = 300
          )
          cat("    Saved individual spatial effect map\n")
        }, error = function(e2) {
          cat("    Could not save spatial effect map:", e2$message, "\n")
        })
      }

      if (!is.null(p3)) {
        tryCatch({
          ggsave(
            file.path(output_dir, paste0(model_name, "_oipc_map.png")),
            p3, width = 8, height = 6, dpi = 300
          )
          cat("    Saved individual OIPC map\n")
        }, error = function(e3) {
          cat("    Could not save OIPC map:", e3$message, "\n")
        })
      }
    })
  }
  
  # ─────────────────────────────────────────────────────────────────────────
  # Store results
  # ─────────────────────────────────────────────────────────────────────────
  
  results_list[[model_name]] <- list(
    model = model_name,
    spatial_var = spatial_var_name,
    n_sites = n_sites,
    correlation = cor_value,
    r_squared = cor_value^2,
    p_value = p_value,
    risk_level = risk_text,
    mean_spatial_effect = mean(spatial_effect),
    sd_spatial_effect = sd(spatial_effect),
    range_spatial_effect = range(spatial_effect)
  )
  
}

#═══════════════════════════════════════════════════════════════════════════════
# CREATE SUMMARY TABLE
#═══════════════════════════════════════════════════════════════════════════════

cat("\n\n", strrep("═", 70), "\n")
cat("SUMMARY ACROSS ALL MODELS\n")
cat(strrep("═", 70), "\n\n")

if (length(results_list) > 0) {
  
  summary_df <- map_dfr(results_list, function(x) {
    data.frame(
      model = x$model,
      spatial_variable = x$spatial_var %||% NA,
      n_sites = x$n_sites,
      n_spatial_values = x$n_spatial %||% x$n_sites,
      correlation = x$correlation %||% NA,
      r_squared = x$r_squared %||% NA,
      p_value = x$p_value %||% NA,
      risk_level = x$risk_level %||% x$note %||% "Unknown"
    )
  })
  
  print(summary_df)
  
  # Save to CSV
  write_csv(summary_df, file.path(output_dir, "spatial_confounding_summary.csv"))
  
  # Create comparison plot
  summary_df_clean <- summary_df %>% filter(!is.na(correlation))
  
  if (nrow(summary_df_clean) > 0) {
    p_summary <- ggplot(summary_df_clean, aes(x = reorder(model, correlation), y = correlation)) +
      geom_col(aes(fill = abs(correlation)), alpha = 0.8) +
      geom_hline(yintercept = c(-0.7, -0.5, -0.3, 0.3, 0.5, 0.7), 
                 linetype = "dashed", alpha = 0.3) +
      scale_fill_gradient(low = "lightblue", high = "darkred", 
                         name = "|r|", limits = c(0, 1)) +
      coord_flip() +
      labs(
        title = "Spatial Confounding Check: All Models",
        subtitle = "Correlation between spatial effect and δ²H_precip",
        x = "Model",
        y = "Pearson Correlation (r)"
      ) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"))
    
    ggsave(
      file.path(output_dir, "confounding_comparison_all_models.png"),
      p_summary, width = 10, height = 6, dpi = 300
    )
  }
  
  # ───────────────────────────────────────────────────────────────────────
  # Print interpretation guide
  # ───────────────────────────────────────────────────────────────────────
  
  cat("\n")
  cat("INTERPRETATION GUIDE\n")
  cat("  |r| > 0.7:        HIGH confounding risk\n")
  cat("  0.5 < |r| <= 0.7: MODERATE risk\n")
  cat("  0.3 < |r| <= 0.5: LOW-MODERATE correlation\n")
  cat("  |r| <= 0.3:       LOW confounding risk\n")
  
} else {
  cat("\n⚠ No results to summarize - check model files\n")
}

cat("\n\n✓ Analysis complete! Results saved to:", output_dir, "\n")
cat("  • Summary table: spatial_confounding_summary.csv\n")
cat("  • Individual plots: *_confounding_scatter.png\n")
cat("  • Spatial patterns: *_spatial_patterns.png\n\n")
