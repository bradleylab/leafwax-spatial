library(leafwax)
library(ggplot2)
library(posterior)

source("posterior_helpers.R")

# Calculate detection thresholds for various scenarios
confidence_levels <- seq(0.5, 0.99, by = 0.01)
rho_values <- c(0, 0.3, 0.5, 0.7, 0.9)

# Residual σ (original per-mille units) for the spatial full_sp model and
# the non-spatial baseline model. Stan exports `sigma` on the standardized
# scale (4d_leaf_wax_spatial_model.stan:225-ish) and `sigma_residual_original =
# sigma * d2H_wax_sd_original` at line 479. We prefer the exported original-
# scale quantity where available; otherwise back-transform.
.sigma_original <- function(model) {
  draws <- load_draws(model)
  vars  <- variables(draws)
  if ("sigma_residual_original" %in% vars) {
    return(mean(as.numeric(
      as_draws_matrix(subset_draws(draws, variable = "sigma_residual_original")))))
  }
  sigma_std <- mean(as.numeric(
    as_draws_matrix(subset_draws(draws, variable = "sigma"))))
  d2H_sd <- load_stan_data(model)$scaling_params$d2H_sd
  sigma_std * d2H_sd
}

sigma_spatial    <- .sigma_original("full_sp")
sigma_nonspatial <- .sigma_original("baseline")
cat(sprintf("sigma (spatial, full_sp)    = %.2f ‰\n", sigma_spatial))
cat(sprintf("sigma (non-spatial, base)   = %.2f ‰\n", sigma_nonspatial))

sigma_analytical <- 3  # measurement repeatability, fixed
sigma_combined_spatial    <- sqrt(sigma_spatial^2    + sigma_analytical^2)
sigma_combined_nonspatial <- sqrt(sigma_nonspatial^2 + sigma_analytical^2)

# Calculate thresholds
results <- expand.grid(
  confidence = confidence_levels,
  rho = rho_values,
  model = c("Spatial", "Non-spatial")
)

results$threshold <- apply(results, 1, function(row) {
  conf <- as.numeric(row['confidence'])
  rho <- as.numeric(row['rho'])
  model <- row['model']
  
  # Get appropriate sigma
  sigma <- ifelse(model == "Spatial", sigma_combined_spatial, sigma_combined_nonspatial)
  
  # Calculate SD of difference with autocorrelation
  sd_diff <- sigma * sqrt(2 * (1 - rho))
  
  # Get threshold for this confidence level (use 2-sided test)
  qnorm((1 + conf) / 2) * sd_diff
})

# Create the figure
p <- ggplot(results, aes(x = confidence, y = threshold, 
                         color = factor(rho), 
                         linetype = model)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(0.5, 1, 0.1), 
                     labels = scales::percent) +
  scale_y_continuous(breaks = seq(0, 80, 10)) +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"),
                     name = "Autocorrelation (ρ)") +
  scale_linetype_manual(values = c("solid", "dotted"),
                        name = "Model Type",
                        labels = c("Non-spatial" = "Non-spatial", 
                                   "Spatial" = "Spatial")) +
  labs(x = "Confidence Level",
       y = "Detection Threshold (‰ δ²H)") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.box = "vertical") +
  guides(color = guide_legend(order = 1),
         linetype = guide_legend(order = 2, 
                                override.aes = list(size = 1.5)))

# Save as PNG
ggsave("detection_thresholds.png", plot = p, width = 8, height = 6, dpi = 300)

print(p)