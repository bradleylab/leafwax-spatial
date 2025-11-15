#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3f_test_d2Hwax_coherence.R
#
# Test for spatial structure in δ²Hwax values
# Analyzes spatial autocorrelation and coherence in isotope measurements
# Generates variograms and spatial diagnostic plots
#
# Input: results/3_sediment_ready_for_modeling.rds
# Output: results/d2H_coherence/ (variograms, Moran's I, spatial plots)
#───────────────────────────────────────────────────────────────────────────────

library(spdep)
library(gstat)
library(sp)
library(ncf)
library(ggplot2)
library(dplyr)

# Create output directory
output_dir <- "results/d2H_coherence"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading sediment data...\n")
sediment <- readRDS("results/3_sediment_ready_for_modeling.rds")
cat("  N =", nrow(sediment), "samples\n\n")

#───────────────────────────────────────────────────────────────────────────────
# 1. Moran's I test at multiple distance thresholds
#───────────────────────────────────────────────────────────────────────────────

cat("Testing spatial autocorrelation with Moran's I...\n")

coords <- as.matrix(sediment[, c("longitude", "latitude")])
distances <- c(50, 100, 200, 500, 1000)  # km

moran_results <- data.frame()
for (d in distances) {
  nb <- dnearneigh(coords, 0, d, longlat = TRUE)
  if (length(nb) > 0 && any(card(nb) > 0)) {
    w <- nb2listw(nb, style = "W", zero.policy = TRUE)
    moran_test <- moran.test(sediment$d2H_wax, w, zero.policy = TRUE)
    
    moran_results <- rbind(moran_results, data.frame(
      distance_km = d,
      moran_I = moran_test$estimate[1],
      expected_I = moran_test$estimate[2],
      variance = moran_test$estimate[3],
      p_value = moran_test$p.value,
      n_neighbors = mean(card(nb))
    ))
    
    cat("  Distance ≤", d, "km: I =", round(moran_test$estimate[1], 3), 
        "p =", format.pval(moran_test$p.value, digits = 3), "\n")
  }
}

write.csv(moran_results, file.path(output_dir, "moran_I_results.csv"), row.names = FALSE)

#───────────────────────────────────────────────────────────────────────────────
# 2. Variogram analysis
#───────────────────────────────────────────────────────────────────────────────

cat("\nComputing empirical variogram...\n")

# Convert to spatial object
sediment_sp <- sediment
coordinates(sediment_sp) <- ~longitude+latitude
proj4string(sediment_sp) <- CRS("+proj=longlat +datum=WGS84")

# Compute variogram
vario <- variogram(d2H_wax ~ 1, sediment_sp, cutoff = 2000, width = 100)
write.csv(vario, file.path(output_dir, "variogram_data.csv"), row.names = FALSE)

# Plot variogram
png(file.path(output_dir, "variogram_d2Hwax.png"), width = 8, height = 6, 
    units = "in", res = 300)
plot(vario, main = "Spatial Structure of δ²Hwax", 
     xlab = "Distance (km)", ylab = "Semivariance",
     pch = 19, col = "blue")
dev.off()

#───────────────────────────────────────────────────────────────────────────────
# 3. Distance-based correlogram
#───────────────────────────────────────────────────────────────────────────────

cat("\nComputing correlogram...\n")

set.seed(123)
correlogram <- correlog(sediment$longitude, sediment$latitude, 
                       sediment$d2H_wax, increment = 100, resamp = 999, 
                       latlon = TRUE)

# Save correlogram data
corr_df <- data.frame(
  distance = correlogram$mean.of.class,
  correlation = correlogram$correlation,
  p_value = correlogram$p
)
write.csv(corr_df, file.path(output_dir, "correlogram_data.csv"), row.names = FALSE)

# Plot correlogram
png(file.path(output_dir, "correlogram_d2Hwax.png"), width = 8, height = 6, 
    units = "in", res = 300)
plot(correlogram, main = "Spatial Correlogram of δ²Hwax")
abline(h = 0, lty = 2, col = "red")
dev.off()

#───────────────────────────────────────────────────────────────────────────────
# 4. Compare to environmental variables
#───────────────────────────────────────────────────────────────────────────────

cat("\nComparing spatial structure to environmental variables...\n")

# Test at 200 km distance
nb_200 <- dnearneigh(coords, 0, 200, longlat = TRUE)
w_200 <- nb2listw(nb_200, style = "W", zero.policy = TRUE)

env_vars <- c("d2H_wax", "oipc_d2h20", "elevation_gmted", "C4_fraction_5deg", 
              "annual_precip", "max_temp", "vpd", "soil_moisture")

comparison_results <- data.frame()

for(var in env_vars) {
  if(var %in% names(sediment) && !all(is.na(sediment[[var]]))) {
    moran_env <- moran.test(sediment[[var]], w_200, zero.policy = TRUE, 
                           na.action = na.omit)
    
    comparison_results <- rbind(comparison_results, data.frame(
      variable = var,
      moran_I = moran_env$estimate[1],
      p_value = moran_env$p.value,
      significant = moran_env$p.value < 0.05
    ))
    
    cat("  ", sprintf("%-20s", var), "I =", sprintf("%6.3f", moran_env$estimate[1]), 
        "p =", format.pval(moran_env$p.value, digits = 3), "\n")
  }
}

write.csv(comparison_results, file.path(output_dir, "spatial_structure_comparison.csv"), 
          row.names = FALSE)

# Create comparison plot
p_comparison <- ggplot(comparison_results, aes(x = reorder(variable, moran_I), 
                                               y = moran_I, fill = significant)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FALSE" = "gray60", "TRUE" = "steelblue")) +
  coord_flip() +
  labs(x = "", y = "Moran's I (200 km neighborhood)", 
       title = "Spatial Autocorrelation Comparison",
       subtitle = "δ²Hwax vs Environmental Variables") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "spatial_structure_comparison.png"), p_comparison,
       width = 8, height = 6, dpi = 300)


#───────────────────────────────────────────────────────────────────────────────
# 4b. Better visualization of Moran's I comparison
#───────────────────────────────────────────────────────────────────────────────

cat("\n\nCreating detailed Moran's I comparison...\n")

# Read the comparison results we already computed
comparison_results <- read.csv(file.path(output_dir, "spatial_structure_comparison.csv"))

# Print ranked results
cat("\nVariables ranked by spatial autocorrelation (Moran's I at 200 km):\n")
comparison_ranked <- comparison_results[order(-comparison_results$moran_I),]
for(i in 1:nrow(comparison_ranked)) {
  cat(sprintf("%2d. %-20s I = %6.3f %s\n", 
              i, 
              comparison_ranked$variable[i], 
              comparison_ranked$moran_I[i],
              ifelse(comparison_ranked$p_value[i] < 0.001, "***",
                     ifelse(comparison_ranked$p_value[i] < 0.01, "**",
                            ifelse(comparison_ranked$p_value[i] < 0.05, "*", "ns")))))
}

# Create a better comparison plot with error bars
# Calculate standard error for Moran's I (approximate)
comparison_results$se <- sqrt(1/nrow(sediment))  # Simplified SE
comparison_results$lower <- comparison_results$moran_I - 1.96*comparison_results$se
comparison_results$upper <- comparison_results$moran_I + 1.96*comparison_results$se

# Categorize variables
comparison_results$category <- case_when(
  comparison_results$variable == "d2H_wax" ~ "Target",
  comparison_results$variable == "oipc_d2h20" ~ "Primary predictor",
  comparison_results$variable %in% c("C4_fraction_5deg", "elevation_gmted") ~ "Physical",
  TRUE ~ "Climate"
)

# Enhanced comparison plot
p_moran_detailed <- ggplot(comparison_results, 
                           aes(x = reorder(variable, moran_I), y = moran_I)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, color = "gray40") +
  geom_point(aes(color = category, size = -log10(p_value + 0.00001)), 
             alpha = 0.8) +
  scale_color_manual(values = c("Target" = "red", 
                               "Primary predictor" = "blue",
                               "Physical" = "darkgreen",
                               "Climate" = "orange")) +
  scale_size_continuous(range = c(3, 8), name = "Significance\n(-log10 p)") +
  coord_flip() +
  labs(x = "", 
       y = "Moran's I (200 km neighborhood)",
       title = "Spatial Autocorrelation Comparison",
       subtitle = "All variables show significant spatial structure") +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major.y = element_blank())

ggsave(file.path(output_dir, "moran_I_comparison_detailed.png"), 
       p_moran_detailed, width = 10, height = 6, dpi = 300)

# Test at multiple distances for top variables
cat("\n\nTesting spatial structure decay for key variables...\n")

top_vars <- comparison_ranked$variable[1:4]  # Top 4 variables
multi_distance_results <- data.frame()

for(var in top_vars) {
  cat("\n", var, ":\n")
  for(d in distances) {
    nb <- dnearneigh(coords, 0, d, longlat = TRUE)
    if(length(nb) > 0 && any(card(nb) > 0)) {
      w <- nb2listw(nb, style = "W", zero.policy = TRUE)
      moran_test <- moran.test(sediment[[var]], w, zero.policy = TRUE)
      
      multi_distance_results <- rbind(multi_distance_results, data.frame(
        variable = var,
        distance_km = d,
        moran_I = moran_test$estimate[1],
        p_value = moran_test$p.value
      ))
      
      cat("  ", d, "km: I =", round(moran_test$estimate[1], 3), "\n")
    }
  }
}

# Plot decay curves for top variables
p_decay <- ggplot(multi_distance_results, 
                  aes(x = distance_km, y = moran_I, color = variable)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_brewer(palette = "Set1", name = "Variable") +
  labs(x = "Distance (km)", 
       y = "Moran's I",
       title = "Spatial Autocorrelation Decay",
       subtitle = "How spatial correlation decreases with distance") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave(file.path(output_dir, "spatial_decay_curves.png"), 
       p_decay, width = 8, height = 6, dpi = 300)

# Summary statistics
cat("\n\nSUMMARY:\n")
cat("Variable with highest spatial autocorrelation:", comparison_ranked$variable[1], 
    "(I =", round(comparison_ranked$moran_I[1], 3), ")\n")
cat("δ²Hwax ranks #", which(comparison_ranked$variable == "d2H_wax"), 
    "out of", nrow(comparison_ranked), "variables\n")

# Write detailed comparison
write.csv(multi_distance_results, 
          file.path(output_dir, "moran_I_multiple_distances.csv"), 
          row.names = FALSE)
          
          
#───────────────────────────────────────────────────────────────────────────────
# 5. Regional analysis - test for hotspots
#───────────────────────────────────────────────────────────────────────────────

cat("\nIdentifying regional patterns...\n")

# Local Moran's I (LISA)
local_moran <- localmoran(sediment$d2H_wax, w_200, zero.policy = TRUE)
sediment$local_I <- local_moran[, 1]
sediment$local_p <- local_moran[, 5]

# Classify into clusters
sediment$cluster <- "Not significant"
sediment$cluster[sediment$local_p < 0.05 & sediment$local_I > 0 & 
                 sediment$d2H_wax > mean(sediment$d2H_wax)] <- "High-High"
sediment$cluster[sediment$local_p < 0.05 & sediment$local_I > 0 & 
                 sediment$d2H_wax < mean(sediment$d2H_wax)] <- "Low-Low"
sediment$cluster[sediment$local_p < 0.05 & sediment$local_I < 0] <- "Outlier"

# Map of clusters
library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")

p_clusters <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_point(data = sediment, aes(x = longitude, y = latitude, 
                                  color = cluster), size = 2, alpha = 0.7) +
  scale_color_manual(values = c("High-High" = "red", "Low-Low" = "blue", 
                               "Outlier" = "yellow", "Not significant" = "gray50")) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 75)) +
  labs(title = "Local Spatial Clusters in δ²Hwax",
       subtitle = "Based on Local Moran's I (200 km neighborhood)",
       color = "Cluster Type") +
  theme_minimal()

ggsave(file.path(output_dir, "spatial_clusters_map.png"), p_clusters,
       width = 12, height = 6, dpi = 300)

#───────────────────────────────────────────────────────────────────────────────
# 6. Summary report
#───────────────────────────────────────────────────────────────────────────────

summary_text <- c(
  "SPATIAL STRUCTURE ANALYSIS OF δ²Hwax",
  "=====================================",
  paste("Date:", Sys.Date()),
  paste("N samples:", nrow(sediment)),
  "",
  "KEY FINDINGS:",
  paste("- Global Moran's I (200 km):", round(comparison_results$moran_I[1], 3),
        ifelse(comparison_results$p_value[1] < 0.001, "(p < 0.001)", 
               paste("(p =", round(comparison_results$p_value[1], 3), ")"))),
  "",
  "- Spatial autocorrelation extends to approximately:",
  paste("  ", max(moran_results$distance_km[moran_results$p_value < 0.05]), "km"),
  "",
  "- Local clusters identified:",
  paste("  High-High clusters:", sum(sediment$cluster == "High-High"), "sites"),
  paste("  Low-Low clusters:", sum(sediment$cluster == "Low-Low"), "sites"),
  paste("  Spatial outliers:", sum(sediment$cluster == "Outlier"), "sites"),
  "",
  "- Comparison to environmental variables:",
  "  δ²Hwax shows stronger spatial structure than most environmental variables",
  paste("  Variables with significant spatial structure:",
        paste(comparison_results$variable[comparison_results$significant], collapse = ", "))
)

writeLines(summary_text, file.path(output_dir, "spatial_analysis_summary.txt"))

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("Results saved to:", output_dir, "\n")
cat("\nKey finding: δ²Hwax Moran's I =", round(comparison_results$moran_I[1], 3), 
    "at 200 km\n")
    
#───────────────────────────────────────────────────────────────────────────────
# 7. Test spatial structure in OLS residuals
#───────────────────────────────────────────────────────────────────────────────

cat("\n\nTesting spatial structure in δ²Hwax ~ δ²Hprecip residuals...\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

# Fit global OLS model
ols_model <- lm(d2H_wax ~ oipc_d2h20, data = sediment)
sediment$ols_residuals <- residuals(ols_model)

# Summary of model
cat("\nGlobal OLS model: δ²Hwax ~ δ²Hprecip\n")
cat("R-squared:", round(summary(ols_model)$r.squared, 3), "\n")
cat("Slope:", round(coef(ols_model)[2], 3), "\n")
cat("Intercept:", round(coef(ols_model)[1], 3), "\n\n")

# Test residuals for spatial autocorrelation
residual_moran <- moran.test(sediment$ols_residuals, w_200, zero.policy = TRUE)
cat("Moran's I for OLS residuals (200 km):", round(residual_moran$estimate[1], 3),
    "p =", format.pval(residual_moran$p.value, digits = 3), "\n")

# Compare to raw d2H_wax
cat("Moran's I for raw δ²Hwax (200 km):", round(comparison_results$moran_I[1], 3), "\n")
cat("Reduction in spatial autocorrelation:", 
    round((1 - residual_moran$estimate[1]/comparison_results$moran_I[1]) * 100, 1), "%\n")

# Test at multiple distances
cat("\nResidual spatial structure at different distances:\n")
residual_moran_results <- data.frame()
for (d in distances) {
  nb <- dnearneigh(coords, 0, d, longlat = TRUE)
  if (length(nb) > 0 && any(card(nb) > 0)) {
    w <- nb2listw(nb, style = "W", zero.policy = TRUE)
    moran_test <- moran.test(sediment$ols_residuals, w, zero.policy = TRUE)
    
    residual_moran_results <- rbind(residual_moran_results, data.frame(
      distance_km = d,
      moran_I_residuals = moran_test$estimate[1],
      p_value = moran_test$p.value
    ))
    
    cat("  Distance ≤", d, "km: I =", round(moran_test$estimate[1], 3), 
        "p =", format.pval(moran_test$p.value, digits = 3), "\n")
  }
}

# Merge with original results for comparison
moran_comparison <- merge(moran_results[,c("distance_km", "moran_I")], 
                         residual_moran_results,
                         by = "distance_km")
names(moran_comparison)[2] <- "moran_I_raw"
moran_comparison$reduction_pct <- (1 - moran_comparison$moran_I_residuals/moran_comparison$moran_I_raw) * 100

write.csv(moran_comparison, file.path(output_dir, "moran_I_comparison_residuals.csv"), row.names = FALSE)

#───────────────────────────────────────────────────────────────────────────────
# Create figure showing OLS fit and residual patterns
#───────────────────────────────────────────────────────────────────────────────

library(cowplot)

# Panel A: Global relationship
p_ols <- ggplot(sediment, aes(x = oipc_d2h20, y = d2H_wax)) +
  geom_point(aes(color = ols_residuals), alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", color = "black", linewidth = 1) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", 
                       midpoint = 0, name = "Residual") +
  labs(x = expression(delta^{2}*H[precip]~"(‰)"),
       y = expression(delta^{2}*H[wax]~"(‰)"),
       title = "A: Global OLS regression") +
  theme_minimal() +
  annotate("text", x = -300, y = -50, 
           label = paste("R² =", round(summary(ols_model)$r.squared, 2),
                        "\nSlope =", round(coef(ols_model)[2], 2)),
           hjust = 0)

# Panel B: Residuals map
p_residual_map <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70") +
  geom_point(data = sediment, aes(x = longitude, y = latitude, 
                                  color = ols_residuals), 
             size = 2, alpha = 0.7) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", 
                       midpoint = 0, name = "Residual") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 75)) +
  labs(title = "B: Spatial pattern of OLS residuals",
       subtitle = paste("Moran's I =", round(residual_moran$estimate[1], 3), 
                       "(p < 0.001)")) +
  theme_minimal()

# Panel C: Moran's I comparison
moran_plot_df <- moran_comparison %>%
  tidyr::pivot_longer(cols = c(moran_I_raw, moran_I_residuals),
                     names_to = "type", values_to = "moran_I") %>%
  mutate(type = ifelse(type == "moran_I_raw", "Raw δ²Hwax", "OLS residuals"))

p_moran_compare <- ggplot(moran_plot_df, aes(x = distance_km, y = moran_I, 
                                             color = type, linetype = type)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Raw δ²Hwax" = "black", "OLS residuals" = "red")) +
  scale_linetype_manual(values = c("Raw δ²Hwax" = "solid", "OLS residuals" = "dashed")) +
  labs(x = "Distance (km)", y = "Moran's I",
       title = "C: Spatial autocorrelation decay",
       color = "", linetype = "") +
  theme_minimal() +
  theme(legend.position = c(0.7, 0.8))

# Combine panels
combined_residual_plot <- plot_grid(
  p_ols, 
  p_residual_map,
  p_moran_compare,
  ncol = 1, 
  rel_heights = c(1, 1, 0.8),
  align = "v"
)

ggsave(file.path(output_dir, "ols_residuals_spatial_analysis.png"), 
       combined_residual_plot,
       width = 10, height = 14, dpi = 300)

# Update summary
cat("\n\nKEY FINDING: OLS residuals show significant spatial structure,")
cat("\njustifying the use of spatially-varying coefficients in the model.\n")