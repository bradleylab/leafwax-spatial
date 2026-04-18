# Load required libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(cowplot)
library(scales)
library(maps)

# Set up parameters
focal_lon <- -90.3108  # Example location (central US)
focal_lat <- 38.6488
max_dist <- 1000  # km
lambda_values <- c(1, 3, 5, 10, 20, 40, 70, 100, 150, 220, 300, 400)
display_lambdas <- c(20, 70, 300)  # For visualization

# Function to calculate great circle distance
calc_distance <- function(lon1, lat1, lon2, lat2) {
  # Haversine formula
  R <- 6371  # Earth radius in km
  dLat <- (lat2 - lat1) * pi / 180
  dLon <- (lon2 - lon1) * pi / 180
  a <- sin(dLat/2)^2 + cos(lat1 * pi / 180) * cos(lat2 * pi / 180) * sin(dLon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  R * c
}

# Create grid for spatial visualization
grid_res <- 0.5  # degrees
lon_seq <- seq(focal_lon - 15, focal_lon + 15, by = grid_res)
lat_seq <- seq(focal_lat - 10, focal_lat + 10, by = grid_res)
grid <- expand.grid(lon = lon_seq, lat = lat_seq)

# Calculate distances from focal point
grid$distance <- mapply(calc_distance, 
                        focal_lon, focal_lat, 
                        grid$lon, grid$lat)

# Calculate weights for different lambda values
for(lambda in display_lambdas) {
  grid[[paste0("weight_", lambda)]] <- exp(-grid$distance / lambda)
  # Set weights to zero beyond 1000 km
  grid[[paste0("weight_", lambda)]][grid$distance > max_dist] <- 0
}

# Get world map for context
world <- ne_countries(scale = "medium", returnclass = "sf")

# Panel A (formerly Panel B): Cross-section showing decay curves
distances <- seq(0, max_dist, by = 5)
decay_data <- purrr::map_df(lambda_values, function(lambda) {
  data.frame(
    distance = distances,
    weight = exp(-distances / lambda),
    lambda = factor(lambda),
    lambda_num = lambda
  )
})

# Add markers for where weight = 0.37 (1/e)
decay_points <- purrr::map_df(lambda_values, function(lambda) {
  data.frame(
    distance = lambda,
    weight = exp(-1),
    lambda = factor(lambda),
    lambda_num = lambda
  )
})

panel_a <- ggplot(decay_data, aes(x = distance, y = weight)) +
  geom_line(aes(color = lambda_num, group = lambda), size = 0.8) +
  geom_point(data = decay_points, 
             aes(color = lambda_num),
             size = 2, shape = 16) +
  geom_hline(yintercept = exp(-1), linetype = "dotted", color = "gray50") +
  annotate("text", x = 950, y = exp(-1) + 0.05, 
           label = "w = 0.37", size = 3, color = "gray50") +
  scale_color_viridis_c(name = "λ (km)", 
                        breaks = c(10, 50, 100, 200, 300, 400)) +
  scale_x_continuous(breaks = seq(0, 1000, 200)) +
  labs(x = "Distance from site (km)",
       y = "Kernel weight",
       title = "A. Exponential decay kernels") +
  theme_minimal() +
  theme(plot.title = element_text(size = 11, face = "bold"),
        legend.position = "right")

# Panel B (formerly Panel C): 2D heatmaps for three different scales
create_heatmap <- function(lambda_val) {
  weight_col <- paste0("weight_", lambda_val)
  
  p <- ggplot() +
    # Add world for land/water
    geom_sf(data = world, fill = NA, color = "gray70", size = 0.2) +
    # Add the raster heatmap
    geom_raster(data = grid, 
                aes(x = lon, y = lat, fill = .data[[weight_col]]),
                interpolate = TRUE) +
    # Overlay US state borders in white for contrast - thinner
    geom_polygon(data = map_data("state"), 
                 aes(x = long, y = lat, group = group),
                 fill = NA, color = "white", size = 0.15, alpha = 0.7) +
    # Add contour at weight = 0.37 with color matching the lambda value
    geom_contour(data = grid,
                 aes(x = lon, y = lat, z = .data[[weight_col]]),
                 breaks = c(0.37),
                 color = viridis(10)[which.min(abs(lambda_values - lambda_val))], 
                 size = 0.5, linetype = "dashed") +
    scale_fill_viridis(name = "Weight", 
                       limits = c(0, 1),
                       breaks = c(0, 0.37, 0.5, 0.75, 1)) +
    coord_sf(xlim = c(focal_lon - 12, focal_lon + 12),
             ylim = c(focal_lat - 8, focal_lat + 8)) +
    labs(title = paste0("λ = ", lambda_val, " km")) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10),
          legend.position = "none",
          axis.title = element_blank())
  
  return(p)
}

heatmaps <- purrr::map(display_lambdas, create_heatmap)

# Add shared legend for heatmaps
legend_data <- data.frame(
  x = 1:100,
  y = 1:100,
  fill = seq(0, 1, length.out = 100)
)

legend_plot <- ggplot(legend_data, aes(x = x, y = y, fill = fill)) +
  geom_raster() +
  scale_fill_viridis(name = "Weight", 
                     limits = c(0, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.key.height = unit(0.8, "cm"))

legend <- get_legend(legend_plot)

# Arrange heatmaps
panel_b <- plot_grid(
  plot_grid(plotlist = heatmaps, nrow = 1),
  legend,
  rel_widths = c(1, 0.1),
  nrow = 1
)

panel_b_titled <- ggdraw() +
  draw_plot(panel_b, y = 0, height = 0.92) +
  draw_label("B. Spatial weight fields at different scales", 
             x = 0.02, y = 0.96, hjust = 0, size = 11, fontface = "bold")

# Combine all panels
figure_3 <- plot_grid(
  panel_a,
  panel_b_titled,
  ncol = 1,
  rel_heights = c(1, 0.8),
  labels = NULL
)

# Save figure. Output name tracks the supplement numbering used in the
# manuscript (`Figure_S1_spatial_weighting`). The Makefile copies this
# PNG + PDF into manuscript/figures/supplement_figs/. Purely analytic —
# no rds dependencies.
ggsave("../figures/supplement_figs/Figure_S1_spatial_weighting.png",
       figure_3,
       width = 10,
       height = 9,
       dpi = 300,
       bg = "white")
ggsave("../figures/supplement_figs/Figure_S1_spatial_weighting.pdf",
       figure_3,
       width = 10,
       height = 9,
       device = cairo_pdf)