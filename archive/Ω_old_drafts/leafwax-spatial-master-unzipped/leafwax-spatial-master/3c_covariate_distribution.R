#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3c_covariate_distribution.R
#
# Performs site clustering analysis based on environmental covariates
# Handles collinearity by removing redundant variables before clustering
# Analyzes OIPC-d2Hwax relationships within environmental clusters
#
# Input: results/3_sediment_ready_for_modeling.rds
# Output: results/site_clustering_revised/sediment_with_[k]clusters_revised.rds
#         results/site_clustering_revised/*.pdf (validation plots)
#         results/site_clustering_revised/clustering_summary_revised.rds
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(cluster)
library(factoextra)
library(corrplot)

# Set output directory
output_dir <- "results/site_clustering_revised"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data
cat("REVISED SITE CLUSTERING ANALYSIS\n")
cat("================================\n\n")
sediment <- readRDS("results/3_sediment_ready_for_modeling.rds")
cat("Loaded", nrow(sediment), "sediment records\n\n")

#───────────────────────────────────────────────────────────────────────────────
# 1. IDENTIFY AND REMOVE COLLINEAR VARIABLES
#───────────────────────────────────────────────────────────────────────────────

cat("STEP 1: HANDLING COLLINEARITY\n")
cat("-----------------------------\n\n")

# First, examine all potential clustering variables
all_vars <- sediment %>%
  select(
    # Geographic
    latitude, longitude, elevation_gmted,
    
    # Climate  
    oipc_d2h20, annual_precip, soil_moisture, max_temp, vpd,
    
    # Vegetation
    C4_fraction_5deg,
    
    # Response (we'll decide whether to include this)
    d2H_wax
  ) %>%
  mutate(
    abs_latitude = abs(latitude)
  )

# Check for PFT data
if ("pft_tree" %in% names(sediment)) {
  all_vars <- all_vars %>%
    mutate(
      pft_tree_mean = map_dbl(sediment$pft_tree, ~ mean(.x, na.rm = TRUE)),
      pft_shrub_mean = map_dbl(sediment$pft_shrub, ~ mean(.x, na.rm = TRUE)),
      pft_grass_mean = map_dbl(sediment$pft_grass, ~ mean(.x, na.rm = TRUE))
    )
}

# Calculate correlation matrix (excluding lat/lon)
cor_matrix <- cor(all_vars %>% select(-latitude, -longitude) %>% na.omit())

# Visualize correlations
pdf(file.path(output_dir, "correlations_before_removal.pdf"), width = 10, height = 10)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         title = "Correlations Before Variable Removal")
dev.off()

# Identify high correlations
high_cor_threshold <- 0.8
high_cor_pairs <- which(abs(cor_matrix) > high_cor_threshold & cor_matrix < 1, arr.ind = TRUE)

cat("High correlations (|r| >", high_cor_threshold, "):\n")
if (nrow(high_cor_pairs) > 0) {
  for (i in 1:nrow(high_cor_pairs)) {
    if (high_cor_pairs[i, 1] < high_cor_pairs[i, 2]) {  # Avoid duplicates
      var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
      var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
      r_val <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
      cat(sprintf("  %s ~ %s: r = %.3f\n", var1, var2, r_val))
    }
  }
}

# Decision rules for which variables to keep
cat("\nVariable selection decisions:\n")

# 1. Temperature: keep max_temp, drop abs_latitude (it's derived anyway)
cat("  Temperature axis: KEEP max_temp, DROP abs_latitude\n")

# 2. Moisture: keep annual_precip, drop soil_moisture
cat("  Moisture axis: KEEP annual_precip, DROP soil_moisture\n")

# 3. OIPC: This is correlated with temp but is our key predictor
# Decision: Keep OIPC but acknowledge it partially represents temperature
cat("  Isotopes: KEEP oipc_d2h20 (despite correlation with temp)\n")

# 4. d2H_wax: This is our response variable
# For clustering environmental conditions, we should EXCLUDE it
cat("  Response: DROP d2H_wax (don't cluster on the outcome)\n")

#───────────────────────────────────────────────────────────────────────────────
# 2. CREATE FINAL CLUSTERING DATASET
#───────────────────────────────────────────────────────────────────────────────

cat("\nCreating clustering dataset with independent variables...\n")

# Select final variables for clustering
cluster_vars_final <- sediment %>%
  select(
    # Geographic (for reference, not clustering)
    latitude,
    longitude,
    
    # Independent environmental variables
    elevation_gmted,      # Elevation
    annual_precip,        # Moisture (not soil_moisture)
    max_temp,             # Temperature (not abs_latitude)
    vpd,                  # Atmospheric dryness
    C4_fraction_5deg,     # Vegetation type
    oipc_d2h20           # Isotopic composition of precipitation
  )

# Add PFT data if available
if ("pft_tree" %in% names(sediment)) {
  cluster_vars_final <- cluster_vars_final %>%
    mutate(
      pft_tree_mean = map_dbl(sediment$pft_tree, ~ mean(.x, na.rm = TRUE)),
      pft_shrub_mean = map_dbl(sediment$pft_shrub, ~ mean(.x, na.rm = TRUE)),
      pft_grass_mean = map_dbl(sediment$pft_grass, ~ mean(.x, na.rm = TRUE))
    )
}

# Remove lat/lon for clustering
cluster_data <- cluster_vars_final %>%
  select(-latitude, -longitude) %>%
  na.omit()

# Save row indices for mapping back
complete_rows <- which(complete.cases(cluster_vars_final %>% select(-latitude, -longitude)))

cat("\nFinal variables for clustering:\n")
cat(paste("  -", names(cluster_data)), sep = "\n")
cat("\nObservations with complete data:", nrow(cluster_data), "/", nrow(sediment), "\n")

# Check final correlations
cor_matrix_final <- cor(cluster_data)
pdf(file.path(output_dir, "correlations_after_removal.pdf"), width = 8, height = 8)
corrplot(cor_matrix_final, 
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         title = "Correlations After Variable Removal")
dev.off()

#───────────────────────────────────────────────────────────────────────────────
# 3. DETERMINE OPTIMAL NUMBER OF CLUSTERS
#───────────────────────────────────────────────────────────────────────────────

cat("\nDetermining optimal number of clusters...\n")

# Standardize variables
cluster_scaled <- scale(cluster_data)

set.seed(123)
k_max <- min(15, nrow(cluster_scaled) - 1)

# Elbow method (WSS)
wss <- sapply(1:k_max, function(k) {
  kmeans(cluster_scaled, k, nstart = 25)$tot.withinss
})

# Silhouette method
sil_width <- sapply(2:k_max, function(k) {
  km <- kmeans(cluster_scaled, k, nstart = 25)
  ss <- silhouette(km$cluster, dist(cluster_scaled))
  mean(ss[, 3])
})

# Gap statistic
gap_stat <- clusGap(cluster_scaled, FUN = kmeans, nstart = 25,
                    K.max = k_max, B = 50, verbose = FALSE)

# Plot validation metrics
pdf(file.path(output_dir, "clustering_validation.pdf"), width = 12, height = 4)
par(mfrow = c(1, 3))

# WSS plot
plot(1:k_max, wss, type = "b", pch = 19,
     xlab = "Number of clusters K",
     ylab = "Total within-clusters sum of squares",
     main = "Elbow Method")

# Silhouette plot
plot(2:k_max, sil_width, type = "b", pch = 19,
     xlab = "Number of clusters K",
     ylab = "Average silhouette width",
     main = "Silhouette Method")
abline(v = which.max(sil_width) + 1, lty = 2, col = "red")

# Gap statistic plot
plot(gap_stat, main = "Gap Statistic")

dev.off()

# Determine optimal k
optimal_k_sil <- which.max(sil_width) + 1
optimal_k_gap <- maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"])

cat("\nCluster validation results:\n")
cat("  Silhouette method suggests:", optimal_k_sil, "clusters\n")
cat("  Gap statistic suggests:", optimal_k_gap, "clusters\n")

# Use silhouette method as primary criterion
optimal_k <- optimal_k_sil

cat("  Selected k =", optimal_k, "clusters\n\n")

#───────────────────────────────────────────────────────────────────────────────
# 4. PERFORM K-MEANS CLUSTERING
#───────────────────────────────────────────────────────────────────────────────

set.seed(123)
km_result <- kmeans(cluster_scaled, centers = optimal_k, nstart = 50)

cat("Cluster sizes:\n")
print(table(km_result$cluster))
cat("\n")

#───────────────────────────────────────────────────────────────────────────────
# 5. ADD CLUSTERS BACK TO ORIGINAL DATA
#───────────────────────────────────────────────────────────────────────────────

sediment$site_cluster <- NA
sediment$site_cluster[complete_rows] <- km_result$cluster
sediment_with_clusters <- sediment

# Save
saveRDS(sediment_with_clusters, 
        file.path(output_dir, paste0("sediment_with_", optimal_k, "clusters_revised.rds")))

#───────────────────────────────────────────────────────────────────────────────
# 6. ANALYZE CLUSTER CHARACTERISTICS AND OIPC RELATIONSHIPS
#───────────────────────────────────────────────────────────────────────────────

cat("Analyzing cluster characteristics...\n")

# Environmental characteristics by cluster
cluster_env <- sediment_with_clusters %>%
  filter(!is.na(site_cluster)) %>%
  group_by(site_cluster) %>%
  summarise(
    n = n(),
    # Geographic
    mean_lat = mean(latitude),
    mean_lon = mean(longitude),
    # Environmental
    mean_elev = mean(elevation_gmted),
    mean_precip = mean(annual_precip),
    mean_temp = mean(max_temp),
    mean_vpd = mean(vpd),
    mean_C4 = mean(C4_fraction_5deg),
    mean_oipc = mean(oipc_d2h20),
    # Response
    mean_d2H = mean(d2H_wax),
    sd_d2H = sd(d2H_wax),
    .groups = "drop"
  )

print(cluster_env)

# OIPC relationships by cluster
cat("\n\nAnalyzing OIPC-d2Hwax relationships by cluster...\n")

cluster_slopes <- sediment_with_clusters %>%
  filter(!is.na(site_cluster)) %>%
  group_by(site_cluster) %>%
  summarise(
    n = n(),
    slope = coef(lm(d2H_wax ~ oipc_d2h20))[2],
    intercept = coef(lm(d2H_wax ~ oipc_d2h20))[1],
    r_squared = summary(lm(d2H_wax ~ oipc_d2h20))$r.squared,
    p_value = summary(lm(d2H_wax ~ oipc_d2h20))$coefficients[2,4],
    oipc_range = max(oipc_d2h20) - min(oipc_d2h20),
    .groups = "drop"
  ) %>%
  arrange(site_cluster)

print(cluster_slopes)

# Test for significant interaction
data_for_test <- sediment_with_clusters %>% filter(!is.na(site_cluster))
interaction_model <- lm(d2H_wax ~ oipc_d2h20 * factor(site_cluster), data = data_for_test)
no_interaction_model <- lm(d2H_wax ~ oipc_d2h20 + factor(site_cluster), data = data_for_test)

anova_result <- anova(no_interaction_model, interaction_model)
p_value <- anova_result$`Pr(>F)`[2]

cat("\n\nTest for cluster × OIPC interaction:\n")
print(anova_result)

if (!is.na(p_value) && p_value < 0.05) {
  cat("\n⚠ SIGNIFICANT INTERACTION (p =", round(p_value, 4), ")\n")
  cat("Different clusters show different OIPC slopes!\n")
} else {
  cat("\n✓ No significant interaction (p =", round(p_value, 4), ")\n")
}

#───────────────────────────────────────────────────────────────────────────────
# 7. VARIABLE IMPORTANCE FOR CLUSTERING
#───────────────────────────────────────────────────────────────────────────────

cat("\n\nVariable importance for cluster separation:\n")

var_importance <- data.frame(
  variable = colnames(cluster_scaled),
  f_statistic = NA,
  p_value = NA
)

for (i in 1:ncol(cluster_scaled)) {
  aov_result <- aov(cluster_scaled[,i] ~ factor(km_result$cluster))
  var_importance$f_statistic[i] <- summary(aov_result)[[1]][["F value"]][1]
  var_importance$p_value[i] <- summary(aov_result)[[1]][["Pr(>F)"]][1]
}

var_importance <- var_importance %>%
  arrange(desc(f_statistic))

print(var_importance)

#───────────────────────────────────────────────────────────────────────────────
# 8. VISUALIZATION
#───────────────────────────────────────────────────────────────────────────────

# PCA visualization
pca_result <- prcomp(cluster_scaled, scale = FALSE)

pdf(file.path(output_dir, "clustering_pca_revised.pdf"), width = 10, height = 8)
fviz_pca_biplot(pca_result,
                habillage = km_result$cluster,
                palette = "Set1",
                addEllipses = TRUE,
                ellipse.type = "confidence",
                title = paste("PCA Biplot - Environmental Clustering (k =", optimal_k, ")"),
                repel = TRUE)
dev.off()

# Geographic distribution
pdf(file.path(output_dir, "cluster_geography_revised.pdf"), width = 12, height = 8)
world_map <- map_data("world")

ggplot() +
  geom_polygon(data = world_map, 
               aes(x = long, y = lat, group = group),
               fill = "lightgray", color = "white") +
  geom_point(data = sediment_with_clusters %>% filter(!is.na(site_cluster)),
             aes(x = longitude, y = latitude, color = factor(site_cluster)),
             size = 3, alpha = 0.8) +
  scale_color_brewer(palette = "Set3", name = "Cluster") +
  theme_minimal() +
  labs(title = "Geographic Distribution of Environmental Clusters",
       subtitle = "Clustering based on independent environmental variables only",
       x = "Longitude", y = "Latitude") +
  coord_quickmap()

dev.off()

#───────────────────────────────────────────────────────────────────────────────
# 9. SUMMARY
#───────────────────────────────────────────────────────────────────────────────

cat("\n\nREVISED CLUSTERING SUMMARY\n")
cat("==========================\n\n")

cat("1. COLLINEARITY HANDLING:\n")
cat("   - Removed: abs_latitude (correlated with max_temp)\n")
cat("   - Removed: soil_moisture (correlated with annual_precip)\n")
cat("   - Removed: d2H_wax (response variable)\n")
cat("   - Kept: oipc_d2h20 despite correlation with temp\n\n")

cat("2. CLUSTERING RESULTS:\n")
cat("   - Optimal k =", optimal_k, "\n")
cat("   - Based on environmental conditions only\n")
cat("   - No double-weighting of correlated information\n\n")

cat("3. KEY FINDINGS:\n")
cat("   - Slope range:", round(min(cluster_slopes$slope), 3), "to", 
    round(max(cluster_slopes$slope), 3), "\n")
cat("   - Clusters with R² > 0.3:", sum(cluster_slopes$r_squared > 0.3), "\n")
cat("   - Clusters with p < 0.05:", sum(cluster_slopes$p_value < 0.05), "\n")

# Save summary
summary_list <- list(
  optimal_k = optimal_k,
  cluster_env = cluster_env,
  cluster_slopes = cluster_slopes,
  interaction_p_value = p_value,
  var_importance = var_importance,
  variables_used = names(cluster_data),
  variables_removed = c("abs_latitude", "soil_moisture", "d2H_wax")
)

saveRDS(summary_list, file.path(output_dir, "clustering_summary_revised.rds"))

cat("\n✓ Revised analysis complete!\n")
cat("✓ Results saved to:", output_dir, "\n")