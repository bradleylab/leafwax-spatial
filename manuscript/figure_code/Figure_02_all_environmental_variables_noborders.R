#───────────────────────────────────────────────────────────────────────────────
# Eight-panel plot: Sample locations over environmental rasters
# Panel A: OIPC δ2H precipitation
# Panel B: Annual precipitation
# Panel C: Temperature (max)
# Panel D: VPD
# Panel E: Soil moisture
# Panel F: Elevation
# Panel G: C4 fraction
# Panel H: Plant Functional Type (PFT)
#───────────────────────────────────────────────────────────────────────────────

library(terra)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(dplyr)
library(readr)
library(cowplot)
library(scales)

#───────────────────────────────────────────────────────────────────────────────
# 1) Read sediment data
#───────────────────────────────────────────────────────────────────────────────

# April 2026 run: sediment frame only exists as .rds (no CSV export)
samples_rds <- "results/c2_run_20260414/3_sediment_ready_for_modeling.rds"
samples_df  <- as_tibble(readRDS(samples_rds))

# Convert to sf
samples_sf <- st_as_sf(
  samples_df,
  coords = c("longitude", "latitude"),
  crs    = 4326,
  remove = FALSE
)

# Load continent outlines. GCA style rules disallow country boundaries on
# global maps, so dissolve ne_countries() by `continent` to get one
# multipolygon per landmass. Antarctica dropped — none of the sites land
# there and it wastes map real estate.
world_outline <- ne_countries(scale = "medium", returnclass = "sf") |>
  dplyr::filter(continent != "Antarctica") |>
  dplyr::group_by(continent) |>
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
  sf::st_make_valid()

# Style defaults for continent outlines on raster panels. Thin, neutral
# lines that don't compete with viridis fills.
continents_layer <- ggplot2::geom_sf(
  data = world_outline, fill = NA, color = "gray40", linewidth = 0.25
)

#───────────────────────────────────────────────────────────────────────────────
# 2) Panel A: OIPC δ2H precipitation
#───────────────────────────────────────────────────────────────────────────────

oipc_rast <- rast("input_data/GlobalPrecip/d2h_MA.tif")
oipc_df <- as.data.frame(oipc_rast, xy = TRUE)
names(oipc_df)[3] <- "d2H_precip"

panel_a <- ggplot() +
  geom_tile(
    data   = oipc_df,
    aes(x = x, y = y, fill = d2H_precip),
    width  = res(oipc_rast)[1],
    height = res(oipc_rast)[2]
  ) +
  scale_fill_viridis_c(
    name      = expression(delta^{2}*H[precip]~"(‰)"),
    option    = "viridis",
    direction = 1,
    limits    = c(-450, 50),
    breaks    = seq(-400, 0, by = 100),
    guide     = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal"
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = bquote(bold("A:"~delta^{2}*"H precipitation")), x = NULL, y = "Latitude") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )

#───────────────────────────────────────────────────────────────────────────────
# 3) Panel B: Annual precipitation (FIXED with better colors and scale)
#───────────────────────────────────────────────────────────────────────────────

precip_rast <- rast("results/2f_TerraClimate_ppt_mean_2001_2019.tif")
precip_df <- as.data.frame(precip_rast, xy = TRUE)
names(precip_df)[3] <- "annual_precip"

panel_b <- ggplot() +
  geom_tile(
    data   = precip_df,
    aes(x = x, y = y, fill = annual_precip),
    width  = res(precip_rast)[1],
    height = res(precip_rast)[2]
  ) +
  scale_fill_gradientn(
    name = "Annual precip (mm)",
    colors = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e", "#003c30"),
    values = c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1),  # Adjust color positions
    trans = "log10",
    breaks = c(100, 300, 1000, 3000),
    labels = c("100", "300", "1000", "3000"),
    limits = c(50, 6000),  # Set explicit limits
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal",
      label.position = "bottom",
      ticks.colour = "black",
      ticks.linewidth = 0.5
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = "B: Annual precipitation", x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )
   
#───────────────────────────────────────────────────────────────────────────────
# 4) Panel C: Temperature (max)
#───────────────────────────────────────────────────────────────────────────────

temp_rast <- rast("results/2f_TerraClimate_tmax_mean_2001_2019.tif")
temp_df <- as.data.frame(temp_rast, xy = TRUE)
names(temp_df)[3] <- "max_temp"

panel_c <- ggplot() +
  geom_tile(
    data   = temp_df,
    aes(x = x, y = y, fill = max_temp),
    width  = res(temp_rast)[1],
    height = res(temp_rast)[2]
  ) +
  scale_fill_gradientn(
    name = "Max temp (°C)",
    colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", 
               "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
    limits = c(-20, 45),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal"
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = "C: Maximum temperature", x = NULL, y = "Latitude") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )

#───────────────────────────────────────────────────────────────────────────────
# 5) Panel D: VPD
#───────────────────────────────────────────────────────────────────────────────

vpd_rast <- rast("results/2f_TerraClimate_vpd_mean_2001_2019.tif")
vpd_df <- as.data.frame(vpd_rast, xy = TRUE)
names(vpd_df)[3] <- "vpd"

# Panel D: VPD (FIXED)
panel_d <- ggplot() +
  geom_tile(
    data   = vpd_df,
    aes(x = x, y = y, fill = vpd),
    width  = res(vpd_rast)[1],
    height = res(vpd_rast)[2]
  ) +
  scale_fill_gradientn(
    name = "VPD (hPa)",
    colors = rev(c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5")),  # Use rev() to reverse
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal"
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = "D: Vapor pressure deficit", x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )

#───────────────────────────────────────────────────────────────────────────────
# 6) Panel E: Soil moisture
#───────────────────────────────────────────────────────────────────────────────

soil_rast <- rast("results/2f_TerraClimate_soil_mean_2001_2019.tif")
soil_df <- as.data.frame(soil_rast, xy = TRUE)
names(soil_df)[3] <- "soil_moisture"

# Filter out NA values and add a small offset to handle zeros
soil_df <- soil_df %>% 
  filter(!is.na(soil_moisture) & soil_moisture > 0)  # Only keep positive values

panel_e <- ggplot() +
  geom_tile(
    data   = soil_df,
    aes(x = x, y = y, fill = soil_moisture),
    width  = res(soil_rast)[1],
    height = res(soil_rast)[2]
  ) +
  scale_fill_gradientn(
    name = "Soil moisture (mm)",
    colors = c("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e"),
    values = c(0, 0.1, 0.3, 0.5, 0.7, 1.0),  # Adjust color distribution
    limits = c(1, 600),  # Start from 1 instead of 0
    trans = "log10",
    breaks = c(1, 10, 50, 100, 200, 400),
    labels = c("1", "10", "50", "100", "200", "400"),
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal"
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = "E: Soil moisture", x = NULL, y = "Latitude") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )

#───────────────────────────────────────────────────────────────────────────────
# 7) Panel F: Elevation
#───────────────────────────────────────────────────────────────────────────────

elev_rast <- rast("input_data/elevation_5KMmn_GMTEDmn.tif")
elev_df <- as.data.frame(elev_rast, xy = TRUE)
names(elev_df)[3] <- "elevation"
elev_df <- elev_df %>% filter(elevation > 0)

panel_f <- ggplot() +
  geom_tile(
    data   = elev_df,
    aes(x = x, y = y, fill = elevation),
    width  = res(elev_rast)[1],
    height = res(elev_rast)[2]
  ) +
  scale_fill_gradientn(
    name = "Elevation (m)",
    colors = c("#2d4a2b", "#5a7a4f", "#8fa674", "#c5d1a5", "#f4e8c1", 
               "#d4a76a", "#b5653a", "#8b4513", "#ffffff"),
    values = scales::rescale(c(0, 50, 200, 500, 1000, 1500, 2000, 3000, 5000)),
    limits = c(0, 5000),
    breaks = c(0, 1000, 2000, 3000, 4000, 5000),
    oob = squish,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal"
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = "F: Elevation", x = NULL, y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )

#───────────────────────────────────────────────────────────────────────────────
# 8) Panel G: C4 fraction
#───────────────────────────────────────────────────────────────────────────────

c4_rast <- rast("results/1_C4_total_mean.tif")
c4_rast_frac <- c4_rast / 100
c4_df <- as.data.frame(c4_rast_frac, xy = TRUE)
names(c4_df)[3] <- "C4_fraction"

panel_g <- ggplot() +
  geom_tile(
    data   = c4_df,
    aes(x = x, y = y, fill = C4_fraction),
    width  = res(c4_rast_frac)[1],
    height = res(c4_rast_frac)[2]
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    direction = 1,
    limits = c(0, 1),
    name = expression(C[4]~"fraction"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust    = 0.5,
      barwidth       = unit(5, "cm"),
      barheight      = unit(0.25, "cm"),
      direction      = "horizontal"
    )
  ) +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = bquote(bold("G:"~C[4]~"vegetation fraction")), x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    legend.position  = "bottom",
    legend.margin    = margin(t = 1, b = 1),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  )

#───────────────────────────────────────────────────────────────────────────────
# 9) Panel H: Plant Functional Type
#───────────────────────────────────────────────────────────────────────────────

pft_rast <- rast("results/2d_MODIS_PFT_3classes_Downsampled.tif")
pft_df <- as.data.frame(pft_rast, xy = TRUE)

pft_cols <- c("Tree", "Shrub", "Grass")
pft_df$total <- rowSums(pft_df[, pft_cols], na.rm = TRUE)

pft_df <- pft_df %>%
  filter(total > 0) %>%
  mutate(
    tree_norm  = Tree / total,
    shrub_norm = Shrub / total,
    grass_norm = Grass / total
  )

tree_color  <- c(0, 0.5, 0)
shrub_color <- c(0.8, 0.6, 0.2)
grass_color <- c(0.9, 0.9, 0.3)

pft_df <- pft_df %>%
  mutate(
    r = tree_norm * tree_color[1] + shrub_norm * shrub_color[1] + grass_norm * grass_color[1],
    g = tree_norm * tree_color[2] + shrub_norm * shrub_color[2] + grass_norm * grass_color[2],
    b = tree_norm * tree_color[3] + shrub_norm * shrub_color[3] + grass_norm * grass_color[3],
    color = rgb(r, g, b)
  )

panel_h <- ggplot() +
  geom_tile(
    data   = pft_df,
    aes(x = x, y = y, fill = color),
    width  = res(pft_rast)[1],
    height = res(pft_rast)[2]
  ) +
  scale_fill_identity() +
  continents_layer +
  geom_sf(data = samples_sf, color = "red", size = 0.4, alpha = 0.8) +
  coord_sf(xlim = c(-180, 180), ylim = c(-75, 90), expand = FALSE) +
  labs(title = "H: Plant Functional Type", x = "Longitude", y = NULL) +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid.major = element_line(color = "white", linewidth = 0.2),
    panel.background = element_rect(fill = "white"),
    plot.title       = element_text(face = "bold", hjust = 0.02, size = 9)
  ) +
  # Add PFT legend
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = rgb(tree_color[1], tree_color[2], tree_color[3]), col = NA)),
    xmin = -170, xmax = -162, ymin = -70, ymax = -66
  ) +
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = rgb(shrub_color[1], shrub_color[2], shrub_color[3]), col = NA)),
    xmin = -170, xmax = -162, ymin = -65, ymax = -61
  ) +
  annotation_custom(
    grob = grid::rectGrob(gp = grid::gpar(fill = rgb(grass_color[1], grass_color[2], grass_color[3]), col = NA)),
    xmin = -170, xmax = -162, ymin = -60, ymax = -56
  ) +
  annotate("text", x = -160, y = -68, label = "Tree", size = 2, hjust = 0) +
  annotate("text", x = -160, y = -63, label = "Shrub", size = 2, hjust = 0) +
  annotate("text", x = -160, y = -58, label = "Grass", size = 2, hjust = 0)

#───────────────────────────────────────────────────────────────────────────────
# 10) Combine panels in 4x2 grid
#───────────────────────────────────────────────────────────────────────────────

combined_plot <- plot_grid(
  panel_a, panel_b,
  panel_c, panel_d,
  panel_e, panel_f,
  panel_g, panel_h,
  ncol = 2, nrow = 4,
  align = "hv",
  axis = "tblr"
)

print(combined_plot)

#───────────────────────────────────────────────────────────────────────────────
# 11) Save the combined figure
#───────────────────────────────────────────────────────────────────────────────

# Save both PDF and PNG. Use the Cairo device family (cairo_pdf for PDF,
# type="cairo" for PNG) so Unicode glyphs — ‰, δ, subscripts — render on
# Linux hosts where the default bitmap font lacks them.
ggsave(
  "manuscript/figures/main_figs/Figure_02_all_environmental_variables.pdf",
  plot = combined_plot,
  width = 12, height = 16,
  device = cairo_pdf
)
ggsave(
  "manuscript/figures/main_figs/Figure_02_all_environmental_variables.png",
  plot = combined_plot,
  width = 12, height = 16,
  dpi = 300,
  type = "cairo"
)

# Print summary statistics
cat("\nSummary Statistics:\n")
cat("==================\n")
cat("Total sample points:", nrow(samples_df), "\n")
cat("\nVariable ranges at sample sites:\n")
cat("δ²H precipitation:", round(min(samples_df$oipc_d2h20, na.rm=TRUE)), "-", 
    round(max(samples_df$oipc_d2h20, na.rm=TRUE)), "‰\n")
cat("Annual precipitation:", round(min(samples_df$annual_precip, na.rm=TRUE)), "-", 
    round(max(samples_df$annual_precip, na.rm=TRUE)), "mm\n")
cat("Max temperature:", round(min(samples_df$max_temp, na.rm=TRUE), 1), "-", 
    round(max(samples_df$max_temp, na.rm=TRUE), 1), "°C\n")
cat("VPD:", round(min(samples_df$vpd, na.rm=TRUE), 1), "-", 
    round(max(samples_df$vpd, na.rm=TRUE), 1), "hPa\n")
cat("Soil moisture:", round(min(samples_df$soil_moisture, na.rm=TRUE)), "-", 
    round(max(samples_df$soil_moisture, na.rm=TRUE)), "mm\n")
cat("Elevation:", round(min(samples_df$elevation_gmted, na.rm=TRUE)), "-", 
    round(max(samples_df$elevation_gmted, na.rm=TRUE)), "m\n")
cat("C₄ fraction:", round(min(samples_df$C4_fraction_5deg, na.rm=TRUE), 1), "-", 
    round(max(samples_df$C4_fraction_5deg, na.rm=TRUE), 1), "%\n")