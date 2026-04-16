library(leafwax)
library(ggplot2)

# Calculate detection thresholds for various scenarios
confidence_levels <- seq(0.5, 0.99, by = 0.01)
rho_values <- c(0, 0.3, 0.5, 0.7, 0.9)

# For spatial model (sigma = 15.2‰)
sigma_spatial <- 15.2
sigma_analytical <- 3
sigma_combined_spatial <- sqrt(sigma_spatial^2 + sigma_analytical^2)

# For non-spatial model (sigma = 21‰)  
sigma_nonspatial <- 21
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