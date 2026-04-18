#───────────────────────────────────────────────────────────────────────────────
# Generate standalone residuals map (Panel B only)
#───────────────────────────────────────────────────────────────────────────────

library(ggplot2)
library(dplyr)
library(rnaturalearth)
library(sf)

# Load sediment frame from the April 2026 run. The Makefile runs this
# script with cwd = manuscript/figure_code/, so the path walks up two
# levels to the repo root.
sediment <- readRDS("../../results/c2_run_20260414/3_sediment_ready_for_modeling.rds")

# Fit OLS model
ols_model <- lm(d2H_wax ~ oipc_d2h20, data = sediment)
sediment$ols_residuals <- residuals(ols_model)

# Load world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create standalone residuals map
p_residual_standalone <- ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray70", linewidth = 0.3) +
  geom_point(data = sediment, 
             aes(x = longitude, y = latitude, color = ols_residuals), 
             size = 2.5, alpha = 0.8) +
  scale_color_gradient2(low = "#2166ac", mid = "gray90", high = "#b2182b", 
                       midpoint = 0, 
                       name = "Residual (‰)",
                       limits = c(-100, 100),  # Adjust based on your data
                       breaks = c(-100, -50, 0, 50, 100),
                       guide = guide_colorbar(
                         title.position = "top",
                         title.hjust = 0.5,
                         barwidth = unit(10, "cm"),
                         barheight = unit(0.4, "cm"),
                         direction = "horizontal"
                       )) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 80), expand = FALSE) +
  labs(title = expression(paste("Spatial pattern of OLS residuals (", delta^{2}, "H"[wax], " ~ ", delta^{2}, "H"[precip], ")")),
       subtitle = paste("Moran's I = 0.614 (p < 0.001); R² =", round(summary(ols_model)$r.squared, 2))) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = 5, b = 5),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    panel.grid = element_line(color = "gray90", linewidth = 0.3),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save as the manuscript's Figure S3. Write a PDF + PNG; the Makefile
# copies them into manuscript/figures/supplement_figs/.
ggsave("../figures/supplement_figs/Figure_S3_residuals.png",
       p_residual_standalone,
       width = 12, height = 7, dpi = 300)
ggsave("../figures/supplement_figs/Figure_S3_residuals.pdf",
       p_residual_standalone,
       width = 12, height = 7, device = cairo_pdf)

# Print summary statistics
cat("\nResiduals summary:\n")
cat("Min:", round(min(sediment$ols_residuals), 1), "‰\n")
cat("Max:", round(max(sediment$ols_residuals), 1), "‰\n")
cat("SD:", round(sd(sediment$ols_residuals), 1), "‰\n")
cat("\nRegions with systematic bias:\n")
cat("- Red areas: δ²Hwax more enriched than predicted\n")
cat("- Blue areas: δ²Hwax more depleted than predicted\n")