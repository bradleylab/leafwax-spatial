#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 5e_spatial_intercept_correlations_weighted.R
#
# Analyze environmental drivers of spatial intercept patterns using
# SPATIALLY-WEIGHTED environmental variables consistent with model fitting
# Uses exponential decay weighting at the fitted integration scale λ ≈ 8 km
#
# Input: model_output/*/fit.rds (fitted spatial models)
#        results/3_sediment_ready_for_modeling.rds (environmental data)
#        prepared_data_0923/stan_data_*.rds (for weighted values)
# Output: figures/spatial_intercept_correlations_weighted.pdf (multi-panel figure)
#         results/spatial_intercept_correlations_weighted.csv (correlation table)
#───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(posterior)
  library(corrplot)
  library(ggpubr)
  library(cowplot)
  library(viridis)
})

source("scripts/posterior_helpers.R")

cat("SPATIAL INTERCEPT CORRELATION ANALYSIS (SPATIALLY-WEIGHTED)\n")
cat("============================================================\n\n")

# Load configuration
source("0_load_config.R")

# Set working directory if needed
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

#───────────────────────────────────────────────────────────────────────────────
# SPATIAL WEIGHTING FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────

# Function to convert degree distances to km
deg_to_km <- function(distances_deg, center_lat) {
  km_per_deg_lat <- 111.132
  km_per_deg_lon <- 111.132 * cos(center_lat * pi / 180)
  avg_km_per_deg <- (km_per_deg_lat + km_per_deg_lon) / 2
  return(distances_deg * avg_km_per_deg)
}

# Function to compute spatially-weighted mean
compute_weighted_mean <- function(values, distances_deg, scale_km, center_lat) {
  if (length(values) == 0 || all(is.na(values))) {
    return(NA_real_)
  }
  valid <- !is.na(values)
  if (sum(valid) == 0) {
    return(NA_real_)
  }
  vals <- values[valid]
  dists_deg <- distances_deg[valid]

  # Convert degrees to km
  dists_km <- deg_to_km(dists_deg, center_lat)

  # Gaussian weights with exponential decay
  weights <- exp(-(dists_km^2) / (2 * scale_km^2))

  # Weighted mean
  weighted_mean <- sum(vals * weights) / sum(weights)
  return(weighted_mean)
}

#───────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA AND MODEL
#───────────────────────────────────────────────────────────────────────────────

cat("Loading environmental data...\n")
data <- load_sediment()

# Find available spatial models (exclude _prepared_data)
all_models <- list.dirs(APRIL_RUN, full.names = FALSE, recursive = FALSE)
all_models <- all_models[!grepl("^_", all_models)]
spatial_models <- all_models[grepl("_sp$", all_models)]

cat("\nAvailable spatial models:\n")
cat(paste("-", spatial_models), sep = "\n")

if ("full_interact_sp" %in% spatial_models) {
  model_name <- "full_interact_sp"
} else if ("full_sp" %in% spatial_models) {
  model_name <- "full_sp"
} else {
  model_name <- spatial_models[1]
}

cat("\nUsing model:", model_name, "\n")
draws     <- load_draws(model_name)
stan_data <- load_stan_data(model_name)

#───────────────────────────────────────────────────────────────────────────────
# 2. EXTRACT SPATIAL EFFECTS AND LAMBDA
#───────────────────────────────────────────────────────────────────────────────

cat("\nExtracting spatial effects and integration scale...\n")

# Extract lambda (spatial integration scale)
vars_present <- variables(draws)
if ("lambda_decay" %in% vars_present) {
  lambda_km <- mean(as.numeric(
    as_draws_matrix(subset_draws(draws, variable = "lambda_decay"))))
  cat("Fitted λ (integration scale):", round(lambda_km, 2), "km\n")
} else {
  lambda_km <- 8.0
  cat("Lambda not found in model, using default:", lambda_km, "km\n")
}

# Spatial intercept contribution in ‰.
# alpha_spatial[i] is on the standardized scale and already includes beta_0
# (4d_leaf_wax_spatial_model.stan:316). The GP residual per obs, in ‰, is
# (alpha_spatial - beta_0) * d2H_sd. Prior implementation multiplied
# alpha_spatial by sigma_intercept_spatial (already in ‰): mixed-units error.
d2H_sd <- stan_data$scaling_params$d2H_sd

alpha_spatial_draws <- as_draws_matrix(subset(draws, variable = "alpha_spatial"))
beta_0_draws <- as.numeric(as_draws_matrix(subset(draws, variable = "beta_0")))
alpha_spatial_scaled <- sweep(alpha_spatial_draws, 1, beta_0_draws, "-") * d2H_sd

alpha_spatial_means <- colMeans(alpha_spatial_scaled)
n_obs <- length(alpha_spatial_means)

cat("Number of observations:", n_obs, "\n")
cat(sprintf("Range of spatial intercept contribution: %.1f to %.1f ‰\n",
            min(alpha_spatial_means), max(alpha_spatial_means)))
cat(sprintf("SD of spatial intercept contribution:   %.2f ‰\n",
            sd(alpha_spatial_means)))

#───────────────────────────────────────────────────────────────────────────────
# 3. COMPUTE SPATIALLY-WEIGHTED ENVIRONMENTAL VARIABLES
#───────────────────────────────────────────────────────────────────────────────

cat("\nComputing spatially-weighted environmental variables at λ =", round(lambda_km, 2), "km...\n")

# Initialize storage for weighted values
weighted_env <- tibble(
  spatial_effect = alpha_spatial_means,
  longitude = data$longitude[1:n_obs],
  latitude = data$latitude[1:n_obs],
  abs_latitude = abs(data$latitude[1:n_obs])
)

# For variables already in stan_data, use those (they're already weighted)
if ("oipc_values" %in% names(stan_data)) {
  # Unstandardize the OIPC values
  oipc_mean <- mean(data$oipc_d2h20, na.rm = TRUE)
  oipc_sd <- sd(data$oipc_d2h20, na.rm = TRUE)

  # Use first column (corresponds to fitted lambda scale)
  weighted_env$oipc <- stan_data$oipc_values[,1] * oipc_sd + oipc_mean
  cat("  ✓ OIPC (δ²H_precip) - using stan_data values\n")
} else {
  # Compute weighted mean if not available
  weighted_env$oipc <- rep(NA, n_obs)
  for (i in 1:n_obs) {
    if (!is.null(data$oipc_values[[i]])) {
      weighted_env$oipc[i] <- compute_weighted_mean(
        data$oipc_values[[i]],
        data$oipc_distances[[i]],
        lambda_km,
        data$latitude[i]
      )
    }
  }
  cat("  ✓ OIPC (δ²H_precip) - computed weighted means\n")
}

# For elevation (also in stan_data)
if ("elevation_values" %in% names(stan_data)) {
  # Unstandardize
  elev_mean <- mean(data$elevation_gmted, na.rm = TRUE)
  elev_sd <- sd(data$elevation_gmted, na.rm = TRUE)
  weighted_env$elevation <- stan_data$elevation_values[,1] * elev_sd + elev_mean
  cat("  ✓ Elevation - using stan_data values\n")
} else {
  weighted_env$elevation <- rep(NA, n_obs)
  for (i in 1:n_obs) {
    if (!is.null(data$elevation_values[[i]])) {
      weighted_env$elevation[i] <- compute_weighted_mean(
        data$elevation_values[[i]],
        data$elevation_distances[[i]],
        lambda_km,
        data$latitude[i]
      )
    }
  }
  cat("  ✓ Elevation - computed weighted means\n")
}

# For C4 (also in stan_data)
if ("c4_values_filled" %in% names(stan_data)) {
  # Unstandardize using CONFIG values
  c4_mean <- CONFIG$c4_standardization$mean
  c4_sd <- CONFIG$c4_standardization$sd
  weighted_env$c4 <- stan_data$c4_values_filled[,1] * c4_sd + c4_mean
  cat("  ✓ C4 fraction - using stan_data values\n")
} else {
  weighted_env$c4 <- rep(NA, n_obs)
  for (i in 1:n_obs) {
    if (!is.null(data$c4_values_filled[[i]])) {
      weighted_env$c4[i] <- compute_weighted_mean(
        data$c4_values_filled[[i]],
        data$c4_distances[[i]],
        lambda_km,
        data$latitude[i]
      )
    }
  }
  cat("  ✓ C4 fraction - computed weighted means\n")
}

# For PFT (if in stan_data)
if ("pft_grass" %in% names(stan_data) && is.matrix(stan_data$pft_grass)) {
  weighted_env$grass <- stan_data$pft_grass[,1]
  weighted_env$tree <- stan_data$pft_tree[,1]
  weighted_env$shrub <- stan_data$pft_shrub[,1]
  cat("  ✓ PFT percentages - using stan_data values\n")
} else {
  # Compute weighted means for PFTs
  weighted_env$grass <- rep(NA, n_obs)
  weighted_env$tree <- rep(NA, n_obs)
  weighted_env$shrub <- rep(NA, n_obs)

  for (i in 1:n_obs) {
    if (!is.null(data$pft_grass[[i]]) && length(data$pft_grass[[i]]) > 0) {
      weighted_env$grass[i] <- compute_weighted_mean(
        data$pft_grass[[i]],
        data$pft_distances[[i]],
        lambda_km,
        data$latitude[i]
      )
      weighted_env$tree[i] <- compute_weighted_mean(
        data$pft_tree[[i]],
        data$pft_distances[[i]],
        lambda_km,
        data$latitude[i]
      )
      weighted_env$shrub[i] <- compute_weighted_mean(
        data$pft_shrub[[i]],
        data$pft_distances[[i]],
        lambda_km,
        data$latitude[i]
      )
    }
  }
  cat("  ✓ PFT percentages - computed weighted means\n")
}

# For climate variables (NOT in stan_data, need to compute)
# These were excluded from models due to collinearity

# VPD
weighted_env$vpd <- rep(NA, n_obs)
for (i in 1:n_obs) {
  if (!is.null(data$tc_vpd_values[[i]])) {
    weighted_env$vpd[i] <- compute_weighted_mean(
      data$tc_vpd_values[[i]],
      data$tc_vpd_distances[[i]],
      lambda_km,
      data$latitude[i]
    )
  }
}
cat("  ✓ VPD - computed weighted means\n")

# Maximum temperature
weighted_env$max_temp <- rep(NA, n_obs)
for (i in 1:n_obs) {
  if (!is.null(data$tc_tmax_values[[i]])) {
    weighted_env$max_temp[i] <- compute_weighted_mean(
      data$tc_tmax_values[[i]],
      data$tc_tmax_distances[[i]],
      lambda_km,
      data$latitude[i]
    )
  }
}
cat("  ✓ Max temperature - computed weighted means\n")

# Soil moisture
weighted_env$soil <- rep(NA, n_obs)
for (i in 1:n_obs) {
  if (!is.null(data$tc_soil_values[[i]])) {
    weighted_env$soil[i] <- compute_weighted_mean(
      data$tc_soil_values[[i]],
      data$tc_soil_distances[[i]],
      lambda_km,
      data$latitude[i]
    )
  }
}
cat("  ✓ Soil moisture - computed weighted means\n")

# Annual precipitation
weighted_env$precip <- rep(NA, n_obs)
for (i in 1:n_obs) {
  if (!is.null(data$tc_ppt_values[[i]])) {
    weighted_env$precip[i] <- compute_weighted_mean(
      data$tc_ppt_values[[i]],
      data$tc_ppt_distances[[i]],
      lambda_km,
      data$latitude[i]
    )
  }
}
cat("  ✓ Annual precipitation - computed weighted means\n")

# Check for missing values
missing_summary <- weighted_env %>%
  summarise_all(~sum(is.na(.)))

cat("\nMissing values per variable:\n")
print(missing_summary[missing_summary > 0])

# Remove rows with any NA for correlation analysis
analysis_clean <- weighted_env %>%
  drop_na()

cat("\nObservations after removing NAs:", nrow(analysis_clean), "\n")

#───────────────────────────────────────────────────────────────────────────────
# 4. CORRELATION ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

cat("\nPerforming correlation analysis with spatially-weighted variables...\n")
cat("All environmental variables weighted at λ =", round(lambda_km, 2), "km\n\n")

# Select variables for correlation
cor_vars <- c(
  "elevation", "c4", "grass", "tree", "shrub",  # Included in models
  "oipc", "abs_latitude", "max_temp", "vpd", "soil", "precip"  # Excluded
)

# Calculate correlations with spatial effect
correlations <- analysis_clean %>%
  select(all_of(cor_vars)) %>%
  map_dbl(~cor(.x, analysis_clean$spatial_effect, use = "complete.obs"))

# Calculate p-values
p_values <- analysis_clean %>%
  select(all_of(cor_vars)) %>%
  map_dbl(~cor.test(.x, analysis_clean$spatial_effect,
                     method = "pearson", use = "complete.obs")$p.value)

# Create correlation table
cor_table <- tibble(
  Variable = cor_vars,
  Variable_Type = ifelse(Variable %in% c("elevation", "c4", "grass", "tree", "shrub"),
                        "Included in model", "Excluded (collinear)"),
  Correlation = correlations,
  P_value = p_values,
  Significance = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01 ~ "**",
    P_value < 0.05 ~ "*",
    TRUE ~ ""
  ),
  Abs_Correlation = abs(Correlation)
) %>%
  arrange(desc(Abs_Correlation))

# Print summary
cat("\nTop correlations with spatial effects (spatially-weighted at λ =", round(lambda_km, 2), "km):\n")
print(cor_table, n = 15)

# Save correlation table
write_csv(cor_table, "results/spatial_intercept_correlations_weighted.csv")
cat("\nCorrelation table saved to results/spatial_intercept_correlations_weighted.csv\n")

#───────────────────────────────────────────────────────────────────────────────
# 5. VISUALIZATION
#───────────────────────────────────────────────────────────────────────────────

cat("\nCreating visualizations...\n")

# Function to create scatter plot with correlation
create_scatter <- function(data, x_var, x_label, color = "steelblue") {

  # Get correlation and p-value
  cor_test <- cor.test(data[[x_var]], data$spatial_effect,
                       method = "pearson", use = "complete.obs")
  cor_val <- cor_test$estimate
  p_val <- cor_test$p.value

  # Format statistics
  cor_text <- sprintf("r = %.3f", cor_val)
  if (p_val < 0.001) {
    p_text <- "p < 0.001"
  } else {
    p_text <- sprintf("p = %.3f", p_val)
  }

  # Create plot
  p <- ggplot(data, aes(x = .data[[x_var]], y = spatial_effect)) +
    geom_point(alpha = 0.6, size = 2, color = color) +
    geom_smooth(method = "lm", se = TRUE, alpha = 0.2, color = "darkgray") +
    labs(x = x_label,
         y = "Spatial Effect (‰)") +
    annotate("text", x = Inf, y = Inf,
             label = paste(cor_text, "\n", p_text),
             hjust = 1.1, vjust = 1.1, size = 3.5) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )

  return(p)
}

# Create plots for variables EXCLUDED from models (key suspects)
p1 <- create_scatter(analysis_clean, "oipc", "δ²H precip (‰)", "darkred")
p2 <- create_scatter(analysis_clean, "abs_latitude", "Absolute Latitude (°)", "darkred")
p3 <- create_scatter(analysis_clean, "max_temp", "Max Temperature (°C)", "darkred")
p4 <- create_scatter(analysis_clean, "vpd", "VPD (kPa)", "darkred")
p5 <- create_scatter(analysis_clean, "soil", "Soil Moisture", "darkred")
p6 <- create_scatter(analysis_clean, "precip", "Annual Precipitation (mm)", "darkred")

# Create plots for variables INCLUDED in models
p7 <- create_scatter(analysis_clean, "elevation", "Elevation (m)", "darkgreen")
p8 <- create_scatter(analysis_clean, "c4", "C4 Fraction", "darkgreen")
p9 <- create_scatter(analysis_clean, "grass", "Grass PFT (%)", "darkgreen")
p10 <- create_scatter(analysis_clean, "tree", "Tree PFT (%)", "darkgreen")
p11 <- create_scatter(analysis_clean, "shrub", "Shrub PFT (%)", "darkgreen")

# Create correlation heatmap
cor_matrix <- analysis_clean %>%
  select(spatial_effect, all_of(cor_vars)) %>%
  cor(use = "complete.obs")

# Reorder for better visualization
cor_order <- cor_table$Variable
cor_matrix_ordered <- cor_matrix[c("spatial_effect", cor_order),
                                 c("spatial_effect", cor_order)]

# Create heatmap plot
pdf("figures/spatial_intercept_correlations_weighted_heatmap.pdf", width = 10, height = 8)
corrplot(cor_matrix_ordered,
         method = "color",
         type = "upper",
         order = "original",
         addCoef.col = "black",
         number.cex = 0.7,
         tl.cex = 0.8,
         tl.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(100),
         title = paste("Correlations with Spatial Effects (λ =", round(lambda_km, 2), "km)"),
         mar = c(0, 0, 2, 0))
dev.off()

# Combine scatter plots into multi-panel figure
combined_plot <- plot_grid(
  # Excluded variables (suspects for driving spatial patterns)
  plot_grid(p1, p2, p3, p4, p5, p6,
            ncol = 3, nrow = 2,
            labels = c("A", "B", "C", "D", "E", "F"),
            label_size = 12),
  # Included variables
  plot_grid(p7, p8, p9, p10, p11, NULL,
            ncol = 3, nrow = 2,
            labels = c("G", "H", "I", "J", "K", ""),
            label_size = 12),
  ncol = 1,
  rel_heights = c(1, 1)
)

# Add overall title
title <- ggdraw() +
  draw_label(
    paste0("Environmental Correlates of Spatial Effects (", model_name, ")"),
    fontface = 'bold',
    size = 14,
    x = 0.5,
    hjust = 0.5
  ) +
  draw_label(
    paste0("All variables spatially-weighted at λ = ", round(lambda_km, 2), " km (model-estimated integration scale)\n",
           "Red: Variables excluded from model (collinear) | Green: Variables included in model"),
    size = 10,
    x = 0.5,
    y = 0.15,
    hjust = 0.5
  )

final_plot <- plot_grid(
  title,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# Save combined figure
ggsave("figures/spatial_intercept_correlations_weighted.pdf",
       plot = final_plot,
       width = 12, height = 10, dpi = 300)

cat("Figures saved to figures/spatial_intercept_correlations_weighted.pdf\n")

#───────────────────────────────────────────────────────────────────────────────
# 6. SUMMARY STATISTICS
#───────────────────────────────────────────────────────────────────────────────

cat("\n" , strrep("=", 60), "\n")
cat("SUMMARY OF KEY FINDINGS\n")
cat(strrep("=", 60), "\n\n")

# Spatial effect statistics
cat("Spatial Effect Statistics:\n")
cat(sprintf("  Range: %.1f to %.1f ‰\n",
            min(alpha_spatial_means), max(alpha_spatial_means)))
cat(sprintf("  Mean: %.1f ‰\n", mean(alpha_spatial_means)))
cat(sprintf("  SD: %.1f ‰\n", sd(alpha_spatial_means)))

cat("\nIntegration Scale:\n")
cat(sprintf("  λ = %.2f km (fitted from model)\n", lambda_km))
cat("  All environmental variables weighted using exponential decay\n")
cat("  at this spatial scale for consistency with model estimation\n")

# Top correlates among excluded variables
cat("\nTop Correlates (Excluded Variables):\n")
excluded_cors <- cor_table %>%
  filter(Variable_Type == "Excluded (collinear)") %>%
  slice_head(n = 3)

for (i in 1:nrow(excluded_cors)) {
  cat(sprintf("  %s: r = %.3f %s\n",
              excluded_cors$Variable[i],
              excluded_cors$Correlation[i],
              excluded_cors$Significance[i]))
}

# Top correlates among included variables
cat("\nTop Correlates (Included Variables):\n")
included_cors <- cor_table %>%
  filter(Variable_Type == "Included in model") %>%
  slice_head(n = 3)

for (i in 1:nrow(included_cors)) {
  cat(sprintf("  %s: r = %.3f %s\n",
              included_cors$Variable[i],
              included_cors$Correlation[i],
              included_cors$Significance[i]))
}

# Key insight
cat("\nKEY INSIGHT:\n")
if (abs(cor_table$Correlation[cor_table$Variable == "oipc"]) > 0.4) {
  cat("Spatially-weighted δ²H_precip shows strong correlation with spatial effects,\n")
  cat("confirming that the spatial patterns capture basin-integrated precipitation\n")
  cat("isotope gradients that were excluded from fixed effects due to collinearity.\n")
  cat("This methodologically coherent analysis (same spatial scale throughout)\n")
  cat("explains why non-spatial models show inflated δ²H_precip slopes.\n")
} else {
  cat("Spatial effects show complex environmental drivers beyond just δ²H_precip\n")
  cat("when assessed at the basin integration scale.\n")
}

cat("\nMETHODOLOGICAL NOTE:\n")
cat("This analysis uses spatially-weighted environmental variables at λ =", round(lambda_km, 2), "km,\n")
cat("matching the integration scale estimated by the model. This ensures an\n")
cat("apples-to-apples comparison between spatial effects and environmental drivers.\n")

cat("\n" , strrep("=", 60), "\n")
cat("Analysis complete!\n")