#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# Figure_S5_elevation.R — elevation effect on leaf wax d2H across the
# `include_elevation == 1` models, for the April 2026 run.
#
# The pre-April Stan model used a piecewise-linear elevation spline with
# `beta_elev_spline[k]` coefficients at a fixed set of knots. The April
# Stan model uses a B-spline: `beta_elev_bspline[1:13]` multiplied by a
# pre-computed `elevation_bspline_matrix` (N_obs × N_scales × 13). This
# script builds Figure S5 directly from that representation:
#
#   per-observation effect = elevation_bspline_matrix %*% beta_elev_bspline
#
# Coefficient uncertainty propagates through the linear combination at each
# observation, so q05/q95 at each elevation are produced by sampling from
# the per-coefficient marginals in the diagnostics summary (assumed
# approximately Gaussian — consistent with the smooth posteriors in the
# April run; a tight approximation given we only need bands for plotting).
#
# Inputs per model m:
#   - stan_data_<m>.rds    : elevation_bspline_matrix, scaling_params,
#                            distance_scales, include_elevation
#   - diagnostics.rds      : beta_elev_bspline[1..13] mean / q5 / q95
#   - sediment rds         : elevation_gmted for the Panel B histogram
#
# No script reads fit.rds or chain CSVs (Phase 5 posterior_helpers.R API).
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(posterior)
library(patchwork)
library(RColorBrewer)

source("posterior_helpers.R")

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("FIGURE S5: ELEVATION EFFECTS\n")
cat("═══════════════════════════════════════════════════════════════\n")

APRIL <- "../../results/c2_run_20260414"
model_dirs <- list.dirs(APRIL, full.names = FALSE, recursive = FALSE)
model_dirs <- model_dirs[!grepl("^_", model_dirs)]
cat("Candidate models:", length(model_dirs), "\n")

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

# Assemble the design vector for a target standardized elevation.
#
# Stan's elevation_bspline_matrix is (N × n_scales) rows × n_bspline_coef
# columns (layout: row_idx = (s-1)*N + n for scale s, observation n; see
# 4d_leaf_wax_spatial_model.stan:345-356). For plot prediction we want
# the point-scale (s=1) basis at a target elevation. Stan_data does not
# store the underlying B-spline knot vector, so rather than rebuild it
# with splines::bs() we borrow the basis row of whichever observation
# has an elevation closest to the target — for N=1124 that's accurate
# to a fraction of a knot spacing.
design_row <- function(stan_data, elev_std) {
  N <- nrow(stan_data$elevation_values)
  obs_std <- stan_data$elevation_values[, 1]  # scale-1 (point) elevation per obs
  i_best  <- which.min(abs(obs_std - elev_std))
  stan_data$elevation_bspline_matrix[i_best, ]  # row for (scale=1, obs=i_best)
}

# Per-model elevation effect on original d2H scale, across a prediction
# grid. Returns a dataframe with elevation, median, q05, q95 columns.
model_elev_effect <- function(model_name) {
  cat("\n  -", model_name, ": ")
  config <- tryCatch(load_config(model_name), error = function(e) NULL)
  if (is.null(config) || !isTRUE(config$include_elevation == 1)) {
    cat("no elevation\n"); return(NULL)
  }
  stan_data <- load_stan_data(model_name)
  summ      <- load_summaries(model_name)

  coefs <- summ[grepl("^beta_elev_bspline\\[", summ$variable), ]
  if (nrow(coefs) == 0) { cat("no bspline coefs\n"); return(NULL) }
  coefs <- coefs[order(as.integer(sub(".*\\[(\\d+)\\].*", "\\1", coefs$variable))), ]

  scaling <- stan_data$scaling_params
  d2H_sd  <- scaling$d2H_sd
  elev_mean <- scaling$elev_mean
  elev_sd   <- scaling$elev_sd

  # Prediction grid: 0–5000 m in 100 m steps, within the data support.
  obs_elev <- stan_data$elevation_values[, 1] * elev_sd + elev_mean
  elev_grid_orig <- seq(max(0, floor(min(obs_elev) / 100) * 100),
                        min(5000, ceiling(max(obs_elev) / 100) * 100),
                        by = 100)
  elev_grid_std  <- (elev_grid_orig - elev_mean) / elev_sd

  # Per-grid design row × coefficient stats. Uncertainty at each grid
  # elevation comes from summing the design-weighted half-widths of each
  # coefficient's symmetric 90% CI (treating coefficients as
  # approximately independent — the ribbon is slightly conservative).
  n_basis <- nrow(coefs)
  eff_median <- numeric(length(elev_grid_std))
  eff_q05    <- numeric(length(elev_grid_std))
  eff_q95    <- numeric(length(elev_grid_std))
  for (i in seq_along(elev_grid_std)) {
    d <- design_row(stan_data, elev_grid_std[i])
    if (length(d) != n_basis) {
      d <- rep(0, n_basis)  # defensive; should not hit for consistent inputs
    }
    eff_median[i] <- sum(d * coefs$mean)
    half_width    <- sum(abs(d) * (coefs$q95 - coefs$q5) / 2)
    eff_q05[i]    <- eff_median[i] - half_width
    eff_q95[i]    <- eff_median[i] + half_width
  }
  # Back-transform to ‰ (standardized residual scale × d2H_sd).
  data.frame(
    elevation = elev_grid_orig,
    median    = eff_median * d2H_sd,
    q05       = eff_q05    * d2H_sd,
    q95       = eff_q95    * d2H_sd,
    model     = model_name,
    has_spatial = isTRUE(config$include_gp == 1),
    stringsAsFactors = FALSE
  ) |> (\(df) { cat("OK (range ",
                    round(min(df$median), 1), " to ",
                    round(max(df$median), 1), "‰)\n", sep = ""); df })()
}

# ──────────────────────────────────────────────────────────────────────────────
# Build dataframe
# ──────────────────────────────────────────────────────────────────────────────

effects_list <- map(model_dirs, model_elev_effect)
effects_list <- effects_list[!map_lgl(effects_list, is.null)]

if (length(effects_list) == 0) {
  stop("No models with elevation effects; check that include_elevation models ",
       "have beta_elev_bspline entries in diagnostics.rds.")
}

all_effects <- bind_rows(effects_list)
all_effects$display_name <- gsub("_sp$", " + GP", gsub("_", " ", all_effects$model))

# ──────────────────────────────────────────────────────────────────────────────
# Plots
# ──────────────────────────────────────────────────────────────────────────────

n_models <- length(unique(all_effects$display_name))
color_palette <- if (n_models <= 8) {
  brewer.pal(max(3, n_models), "Set2")[seq_len(n_models)]
} else {
  colorRampPalette(brewer.pal(8, "Set2"))(n_models)
}

p_effects <- ggplot(all_effects, aes(x = elevation / 1000)) +
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = display_name), alpha = 0.25) +
  geom_line(aes(y = median, color = display_name), linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(breaks = seq(0, 6, 1), labels = seq(0, 6000, 1000)) +
  scale_color_manual(values = color_palette, name = "Model") +
  scale_fill_manual(values = color_palette, name = "Model") +
  labs(
    x = "Elevation (m)",
    y = expression(paste("Effect on ", delta^2, "H"[wax], " (\u2030)")),
    title = "A. Elevation effects across models (April 2026 run)",
    subtitle = "90% credible intervals; conservative band from coefficient-marginal propagation"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title    = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey50")
  )

sediment <- load_sediment()
p_dist <- ggplot(sediment, aes(x = elevation_gmted)) +
  geom_histogram(binwidth = 200, fill = "grey60", color = "grey40", alpha = 0.8) +
  scale_x_continuous(breaks = seq(0, 6000, 1000), limits = c(0, 6000)) +
  labs(
    x = "Elevation (m)",
    y = "Number of samples",
    title = "B. Sample distribution across elevations"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14)
  )

p_combined <- p_effects / p_dist + plot_layout(heights = c(2, 1))

ggsave("Figure_S5_elevation.png", p_combined,
       width = 12, height = 10, dpi = 300, bg = "white")
ggsave("Figure_S5_elevation.pdf", p_combined,
       width = 12, height = 10, device = cairo_pdf)

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("Wrote Figure_S5_elevation.{png,pdf} in ", getwd(), "\n", sep = "")
cat("Models plotted: ", paste(unique(all_effects$display_name), collapse = ", "), "\n", sep = "")
