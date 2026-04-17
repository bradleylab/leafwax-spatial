#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 5e_spatial_intercept_correlations.R
#
# Analyze environmental drivers of spatial intercept patterns
# Correlates spatial random effects with environmental variables to understand
# what drives the ±60‰ regional baseline fractionation differences
#
# Input: model_output/*/fit.rds (fitted spatial models)
#        results/3_sediment_ready_for_modeling.rds (environmental data)
# Output: figures/spatial_intercept_correlations.pdf (multi-panel figure)
#         results/spatial_intercept_correlations.csv (correlation table)
#───────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(tidyverse)
  library(cmdstanr)
  library(posterior)
  library(corrplot)
  library(ggpubr)
  library(cowplot)
  library(viridis)
})

cat("SPATIAL INTERCEPT CORRELATION ANALYSIS\n")
cat("=======================================\n\n")

# Load configuration
source("0_load_config.R")

# Set working directory if needed
if (!dir.exists("figures")) dir.create("figures", recursive = TRUE)
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

#───────────────────────────────────────────────────────────────────────────────
# 1. LOAD DATA AND MODEL
#───────────────────────────────────────────────────────────────────────────────

cat("Loading environmental data...\n")
data <- readRDS("results/3_sediment_ready_for_modeling.rds")

# Find available spatial models
model_dirs <- list.dirs("model_output", recursive = FALSE)
spatial_models <- model_dirs[grep("_sp", basename(model_dirs))]

cat("\nAvailable spatial models:\n")
cat(paste("-", basename(spatial_models)), sep = "\n")

# Select best model (prefer full_sp or full_interact_sp)
if ("model_output/full_interact_sp" %in% spatial_models) {
  model_path <- "model_output/full_interact_sp/fit.rds"
  model_name <- "full_interact_sp"
} else if ("model_output/full_sp" %in% spatial_models) {
  model_path <- "model_output/full_sp/fit.rds"
  model_name <- "full_sp"
} else {
  # Use first available spatial model
  model_path <- file.path(spatial_models[1], "fit.rds")
  model_name <- basename(spatial_models[1])
}

cat("\nUsing model:", model_name, "\n")
cat("Loading fitted model...\n")
fit <- readRDS(model_path)

#───────────────────────────────────────────────────────────────────────────────
# 2. EXTRACT SPATIAL EFFECTS
#───────────────────────────────────────────────────────────────────────────────

cat("\nExtracting spatial effects (back-transformed to ‰)...\n")

# Stan defines alpha_spatial[i] on the standardized response scale, already
# including the global intercept beta_0 (4d_leaf_wax_spatial_model.stan:316,
# 331-335). The GP intercept *residual* contribution at obs i, in standardized
# units, is alpha_spatial[i] - beta_0. Multiply by stan_data$scaling_params$d2H_sd
# to get per-mille. This replaces a prior units-error that multiplied the raw
# alpha_spatial draws by sigma_intercept_spatial (a different quantity, already
# in ‰).
stan_data_path <- file.path("prepared_data", paste0("stan_data_", model_name, ".rds"))
if (!file.exists(stan_data_path)) {
  stop("stan_data not found at: ", stan_data_path)
}
stan_data <- readRDS(stan_data_path)
d2H_sd <- stan_data$scaling_params$d2H_sd

draws <- fit$draws()
alpha_spatial_draws <- as_draws_matrix(subset(draws, variable = "alpha_spatial"))
beta_0_draws <- as.numeric(as_draws_matrix(subset(draws, variable = "beta_0")))

# Broadcast beta_0 across columns (one observation per column) and convert.
alpha_spatial_residual_std <- sweep(alpha_spatial_draws, 1, beta_0_draws, "-")
alpha_spatial_scaled <- alpha_spatial_residual_std * d2H_sd

alpha_spatial_means <- colMeans(alpha_spatial_scaled)
n_obs <- length(alpha_spatial_means)

cat("Number of observations:", n_obs, "\n")
cat(sprintf("Range of spatial intercept contribution: %.1f to %.1f ‰\n",
            min(alpha_spatial_means), max(alpha_spatial_means)))
cat(sprintf("SD of spatial intercept contribution:   %.2f ‰\n",
            sd(alpha_spatial_means)))

#───────────────────────────────────────────────────────────────────────────────
# 3. PREPARE ENVIRONMENTAL VARIABLES
#───────────────────────────────────────────────────────────────────────────────

cat("\nPreparing environmental variables...\n")

# Create analysis dataframe
analysis_df <- data %>%
  slice(1:n_obs) %>%  # Match number of observations in model
  mutate(
    spatial_effect = alpha_spatial_means,
    abs_latitude = abs(latitude)
  ) %>%
  mutate(
    # Extract means from PFT lists
    grass = map_dbl(pft_grass, ~mean(.x, na.rm = TRUE)),
    tree = map_dbl(pft_tree, ~mean(.x, na.rm = TRUE)),
    shrub = map_dbl(pft_shrub, ~mean(.x, na.rm = TRUE))
  ) %>%
  select(
    # Spatial effect
    spatial_effect,

    # Variables INCLUDED in models
    elevation = elevation_gmted,
    c4 = c4_mean_filled,
    grass,
    tree,
    shrub,

    # Variables EXCLUDED for collinearity (key suspects)
    oipc = oipc_d2h20,  # δ²H_precip
    abs_latitude,
    max_temp,
    vpd,
    soil = soil_moisture,
    precip = annual_precip,

    # Location info for reference
    longitude, latitude
  )

# Check for missing values
missing_summary <- analysis_df %>%
  summarise_all(~sum(is.na(.)))

cat("\nMissing values per variable:\n")
print(missing_summary[missing_summary > 0])

# Remove rows with any NA for correlation analysis
analysis_clean <- analysis_df %>%
  drop_na()

cat("\nObservations after removing NAs:", nrow(analysis_clean), "\n")

#───────────────────────────────────────────────────────────────────────────────
# 4. CORRELATION ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

cat("\nPerforming correlation analysis...\n")

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
cat("\nTop correlations with spatial effects:\n")
print(cor_table, n = 15)

# Save correlation table
write_csv(cor_table, "results/spatial_intercept_correlations.csv")
cat("\nCorrelation table saved to results/spatial_intercept_correlations.csv\n")

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
  p <- ggplot(data, aes_string(x = x_var, y = "spatial_effect")) +
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
pdf("figures/spatial_intercept_correlations_heatmap.pdf", width = 10, height = 8)
corrplot(cor_matrix_ordered,
         method = "color",
         type = "upper",
         order = "original",
         addCoef.col = "black",
         number.cex = 0.7,
         tl.cex = 0.8,
         tl.col = "black",
         col = colorRampPalette(c("blue", "white", "red"))(100),
         title = "Correlations with Spatial Effects",
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
    "Red: Variables excluded from model (collinear) | Green: Variables included in model",
    size = 10,
    x = 0.5,
    y = 0.2,
    hjust = 0.5
  )

final_plot <- plot_grid(
  title,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.08, 1)
)

# Save combined figure
ggsave("figures/spatial_intercept_correlations.pdf",
       plot = final_plot,
       width = 12, height = 10, dpi = 300)

cat("Figures saved to figures/spatial_intercept_correlations.pdf\n")

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
  cat("δ²H_precip shows strong correlation with spatial effects,\n")
  cat("suggesting that the spatial patterns capture precipitation isotope\n")
  cat("gradients that were excluded from the fixed effects due to collinearity.\n")
  cat("This explains why non-spatial models show inflated δ²H_precip slopes.\n")
} else {
  cat("Spatial effects show complex environmental drivers beyond just δ²H_precip.\n")
}

cat("\n" , strrep("=", 60), "\n")
cat("Analysis complete!\n")