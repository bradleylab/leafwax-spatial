#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# 3d_analyze_collinearity.R
#
# Comprehensive collinearity analysis for spatial leaf wax model predictors
# Evaluates multicollinearity among environmental covariates using VIF and PCA
# Generates diagnostic plots for predictor selection
#
# Input: results/3_sediment_ready_for_modeling.rds
# Output: results/collinearity_analysis/ (diagnostic plots and summaries)
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(corrplot)
library(car)
library(glmnet)
library(cluster)

# Create output directory
output_dir <- "results/collinearity"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data
sediment <- readRDS("results/3_sediment_ready_for_modeling.rds")

#───────────────────────────────────────────────────────────────────────────────
# 1. CORRELATION ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

# Extract main predictors at site level
# First check if PFT columns are lists or numeric
if (is.list(sediment$pft_tree)) {
  # PFT columns are lists, need to extract means
  predictors <- sediment %>%
    select(
      # Climate variables
      oipc_d2h20,
      annual_precip,
      soil_moisture,
      max_temp,
      vpd,
      # Elevation
      elevation_gmted,
      # Vegetation
      C4_fraction_5deg
    ) %>%
    # Add PFT means
    mutate(
      pft_tree_mean = map_dbl(sediment$pft_tree, ~ mean(.x, na.rm = TRUE)),
      pft_shrub_mean = map_dbl(sediment$pft_shrub, ~ mean(.x, na.rm = TRUE)),
      pft_grass_mean = map_dbl(sediment$pft_grass, ~ mean(.x, na.rm = TRUE))
    )
} else {
  # PFT columns are already numeric
  predictors <- sediment %>%
    select(
      # Climate variables
      oipc_d2h20,
      annual_precip,
      soil_moisture,
      max_temp,
      vpd,
      # Elevation
      elevation_gmted,
      # Vegetation
      C4_fraction_5deg,
      starts_with("pft_")
    )
}

# Remove any remaining non-numeric columns
predictors <- predictors %>%
  select(where(is.numeric))

# Correlation matrix
cor_matrix <- cor(predictors, use = "complete.obs")

# Visualize correlations
pdf(file.path(output_dir, "predictor_correlations.pdf"), width = 10, height = 10)
corrplot(cor_matrix, 
         method = "color",
         type = "upper",
         order = "hclust",
         addCoef.col = "black",
         tl.col = "black",
         tl.srt = 45,
         diag = FALSE,
         title = "Predictor Correlations")
dev.off()

# Identify high correlations (|r| > 0.7)
high_cor_pairs <- which(abs(cor_matrix) > 0.7 & cor_matrix != 1, arr.ind = TRUE)
if (nrow(high_cor_pairs) > 0) {
  cat("\nHigh correlation pairs (|r| > 0.7):\n")
  for (i in 1:nrow(high_cor_pairs)) {
    if (high_cor_pairs[i, 1] < high_cor_pairs[i, 2]) {  # Avoid duplicates
      var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
      var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
      r_val <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
      cat(sprintf("  %s ~ %s: r = %.3f\n", var1, var2, r_val))
    }
  }
}

#───────────────────────────────────────────────────────────────────────────────
# 2. VARIANCE INFLATION FACTORS
#───────────────────────────────────────────────────────────────────────────────

# Fit a basic linear model to calculate VIFs
lm_model <- lm(d2H_wax ~ ., data = cbind(d2H_wax = sediment$d2H_wax, predictors))
vif_values <- vif(lm_model)

cat("\nVariance Inflation Factors:\n")
vif_df <- data.frame(
  Variable = names(vif_values),
  VIF = vif_values
) %>%
  arrange(desc(VIF))
print(vif_df)

# Flag problematic VIFs (> 5 or > 10)
cat("\nVariables with VIF > 5:\n")
print(vif_df %>% filter(VIF > 5))

#───────────────────────────────────────────────────────────────────────────────
# 3. PRINCIPAL COMPONENT ANALYSIS
#───────────────────────────────────────────────────────────────────────────────

# Standardize predictors
predictors_scaled <- scale(predictors)

# PCA
pca_result <- prcomp(predictors_scaled, center = FALSE, scale = FALSE)

# Variance explained
var_explained <- summary(pca_result)$importance[2,]
cumvar_explained <- summary(pca_result)$importance[3,]

cat("\nPCA Variance Explained:\n")
for (i in 1:min(5, length(var_explained))) {
  cat(sprintf("  PC%d: %.1f%% (Cumulative: %.1f%%)\n", 
              i, var_explained[i]*100, cumvar_explained[i]*100))
}

# Plot PCA
pdf(file.path(output_dir, "pca_analysis.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

# Scree plot
plot(var_explained, type = "b", 
     xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained",
     main = "PCA Scree Plot")

# Biplot
biplot(pca_result, scale = 0, main = "PCA Biplot")
dev.off()

# Loadings for first 3 PCs
cat("\nPCA Loadings (first 3 components):\n")
print(round(pca_result$rotation[, 1:3], 3))

#───────────────────────────────────────────────────────────────────────────────
# 4. HIERARCHICAL CLUSTERING OF VARIABLES
#───────────────────────────────────────────────────────────────────────────────

# Distance matrix based on correlations
dist_matrix <- as.dist(1 - abs(cor_matrix))
hclust_result <- hclust(dist_matrix, method = "complete")

# Plot dendrogram
pdf(file.path(output_dir, "variable_clustering.pdf"), width = 10, height = 6)
plot(hclust_result, main = "Hierarchical Clustering of Variables",
     xlab = "", sub = "", ylab = "Distance (1 - |correlation|)")
abline(h = 0.3, col = "red", lty = 2)  # Corresponds to |r| = 0.7
dev.off()

#───────────────────────────────────────────────────────────────────────────────
# 5. ELASTIC NET VARIABLE IMPORTANCE
#───────────────────────────────────────────────────────────────────────────────

# Prepare data for glmnet
X <- as.matrix(predictors_scaled)
y <- scale(sediment$d2H_wax)

# Fit elastic net with cross-validation
set.seed(123)
cv_enet <- cv.glmnet(X, y, alpha = 0.5, nfolds = 10)

# Extract coefficients at optimal lambda
coef_matrix <- coef(cv_enet, s = "lambda.min")
coef_values <- as.numeric(coef_matrix)[-1]  # Exclude intercept
coef_names <- rownames(coef_matrix)[-1]  # Exclude intercept

# Plot variable importance
pdf(file.path(output_dir, "elastic_net_importance.pdf"), width = 8, height = 6)
coef_df <- data.frame(
  Variable = coef_names,
  Coefficient = coef_values
) %>%
  filter(abs(Coefficient) > 0) %>%
  arrange(desc(abs(Coefficient)))

if (nrow(coef_df) > 0) {
  p <- ggplot(coef_df, aes(x = reorder(Variable, abs(Coefficient)), y = Coefficient)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Elastic Net Variable Importance",
         x = "Variable", y = "Standardized Coefficient") +
    theme_minimal()
  print(p)
}
dev.off()

#───────────────────────────────────────────────────────────────────────────────
# 6. RECOMMENDATIONS FOR MODEL FITTING
#───────────────────────────────────────────────────────────────────────────────

cat("\n\nCOLLINEARITY ANALYSIS SUMMARY\n")
cat("==============================\n\n")

# Strategy 1: Remove highly correlated variables
cat("Strategy 1: Remove one variable from each highly correlated pair\n")
vars_to_remove <- c()
if (nrow(high_cor_pairs) > 0) {
  # Simple greedy approach - keep variable with lower mean correlation
  for (i in 1:nrow(high_cor_pairs)) {
    if (high_cor_pairs[i, 1] < high_cor_pairs[i, 2]) {
      var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
      var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
      mean_cor1 <- mean(abs(cor_matrix[var1, ]), na.rm = TRUE)
      mean_cor2 <- mean(abs(cor_matrix[var2, ]), na.rm = TRUE)
      if (mean_cor1 > mean_cor2) {
        vars_to_remove <- c(vars_to_remove, var1)
      } else {
        vars_to_remove <- c(vars_to_remove, var2)
      }
    }
  }
  vars_to_remove <- unique(vars_to_remove)
  cat("  Suggested variables to remove:", paste(vars_to_remove, collapse = ", "), "\n")
} else {
  cat("  No highly correlated pairs found\n")
}

# Strategy 2: Use principal components
cat("\nStrategy 2: Use principal components\n")
n_pcs_80 <- which(cumvar_explained >= 0.8)[1]
cat(sprintf("  First %d PCs explain 80%% of variance\n", n_pcs_80))

# Strategy 3: Use regularization
cat("\nStrategy 3: Use regularization in Bayesian framework\n")
cat("  - Use horseshoe or regularized horseshoe priors\n")
cat("  - Use spike-and-slab priors for variable selection\n")
cat("  - Use hierarchical shrinkage for groups of related variables\n")

# Save summary
summary_list <- list(
  correlation_matrix = cor_matrix,
  high_correlations = high_cor_pairs,
  vif_values = vif_df,
  pca_result = pca_result,
  elastic_net_coefs = coef_df,
  recommendations = list(
    vars_to_remove = vars_to_remove,
    n_pcs_for_80_var = n_pcs_80
  )
)

saveRDS(summary_list, file.path(output_dir, "collinearity_analysis_results.rds"))

#───────────────────────────────────────────────────────────────────────────────
# 7. SPATIAL CORRELATION IN PREDICTORS
#───────────────────────────────────────────────────────────────────────────────

# Check if predictors themselves have spatial structure
if (require(gstat, quietly = TRUE) && require(sp, quietly = TRUE)) {
  
  # Convert to spatial object
  coords <- sediment[, c("longitude", "latitude")]
  sp_data <- SpatialPointsDataFrame(coords, predictors)
  
  # Calculate variograms for each predictor
  pdf(file.path(output_dir, "predictor_variograms.pdf"), width = 12, height = 10)
  par(mfrow = c(3, 4))
  
  for (var in names(predictors)) {
    if (sum(!is.na(predictors[[var]])) > 30) {  # Need enough data
      form <- as.formula(paste(var, "~ 1"))
      vgm <- variogram(form, sp_data, cutoff = 2000)
      plot(vgm, main = var, pch = 19)
    }
  }
  dev.off()
  
} else {
  cat("\nNote: Install 'gstat' and 'sp' packages for spatial variogram analysis\n")
}

cat("\nCollinearity analysis complete. Results saved to", output_dir, "\n")

#───────────────────────────────────────────────────────────────────────────────
# 8. CREATE SUMMARY REPORT
#───────────────────────────────────────────────────────────────────────────────

# Write a text summary
summary_file <- file.path(output_dir, "collinearity_summary.txt")
sink(summary_file)

cat("COLLINEARITY ANALYSIS SUMMARY REPORT\n")
cat("====================================\n")
cat("Date:", as.character(Sys.Date()), "\n")
cat("Number of predictors analyzed:", ncol(predictors), "\n")
cat("Number of observations:", nrow(predictors), "\n\n")

cat("PREDICTORS INCLUDED:\n")
cat(paste("-", names(predictors)), sep = "\n")

cat("\n\nHIGH CORRELATIONS (|r| > 0.7):\n")
if (nrow(high_cor_pairs) > 0) {
  for (i in 1:nrow(high_cor_pairs)) {
    if (high_cor_pairs[i, 1] < high_cor_pairs[i, 2]) {
      var1 <- rownames(cor_matrix)[high_cor_pairs[i, 1]]
      var2 <- colnames(cor_matrix)[high_cor_pairs[i, 2]]
      r_val <- cor_matrix[high_cor_pairs[i, 1], high_cor_pairs[i, 2]]
      cat(sprintf("  %s ~ %s: r = %.3f\n", var1, var2, r_val))
    }
  }
} else {
  cat("  None found\n")
}

cat("\n\nVARIANCE INFLATION FACTORS:\n")
print(vif_df)

cat("\n\nPRINCIPAL COMPONENTS:\n")
for (i in 1:min(5, length(var_explained))) {
  cat(sprintf("  PC%d: %.1f%% variance (Cumulative: %.1f%%)\n", 
              i, var_explained[i]*100, cumvar_explained[i]*100))
}

cat("\n\nRECOMMENDATIONS:\n")
cat("1. Variables to consider removing due to collinearity:\n")
if (length(vars_to_remove) > 0) {
  cat("  ", paste(vars_to_remove, collapse = ", "), "\n")
} else {
  cat("   None\n")
}

cat("\n2. Variables with high VIF (>5):\n")
high_vif <- vif_df %>% filter(VIF > 5)
if (nrow(high_vif) > 0) {
  for (i in 1:nrow(high_vif)) {
    cat(sprintf("   %s: VIF = %.2f\n", high_vif$Variable[i], high_vif$VIF[i]))
  }
} else {
  cat("   None\n")
}

cat("\n3. Elastic Net selected variables:\n")
if (exists("coef_df") && nrow(coef_df) > 0) {
  for (i in 1:nrow(coef_df)) {
    cat(sprintf("   %s: coefficient = %.3f\n", 
                coef_df$Variable[i], coef_df$Coefficient[i]))
  }
} else {
  cat("   None selected\n")
}

sink()

cat("\nSummary report written to:", summary_file, "\n")