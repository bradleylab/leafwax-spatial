#!/usr/bin/env Rscript
#───────────────────────────────────────────────────────────────────────────────
# Figure_S2_pairwise_correlations.R — pairwise Pearson correlations among
# model predictors, restricted to the April 2026 run's 1131-site sediment
# frame. Equivalent to the correlation-plot block in
# 3d_analyze_collinearity.R (pre-April layout, writes to
# results/collinearity/) — lifted into figure_code/ so the supplement
# figure regenerates via `make figures` without pulling in the full
# collinearity pipeline.
#───────────────────────────────────────────────────────────────────────────────

library(dplyr)
library(purrr)
library(corrplot)

source("posterior_helpers.R")

sediment <- load_sediment()

# Same predictor set that 3d_analyze_collinearity.R uses (handles the
# list-column / numeric PFT representations both ways).
if (is.list(sediment$pft_tree)) {
  predictors <- sediment %>%
    select(oipc_d2h20, annual_precip, soil_moisture, max_temp, vpd,
           elevation_gmted, C4_fraction_5deg) %>%
    mutate(
      pft_tree_mean  = map_dbl(sediment$pft_tree,  ~ mean(.x, na.rm = TRUE)),
      pft_shrub_mean = map_dbl(sediment$pft_shrub, ~ mean(.x, na.rm = TRUE)),
      pft_grass_mean = map_dbl(sediment$pft_grass, ~ mean(.x, na.rm = TRUE))
    )
} else {
  predictors <- sediment %>%
    select(oipc_d2h20, annual_precip, soil_moisture, max_temp, vpd,
           elevation_gmted, C4_fraction_5deg,
           starts_with("pft_"))
}
predictors <- predictors %>% select(where(is.numeric))

cor_matrix <- cor(predictors, use = "complete.obs")

# Emit PDF (cairo_pdf for Unicode) + PNG at the manuscript supplement name.
cairo_pdf("../figures/supplement_figs/Figure_S2_pairwise_correlations.pdf", width = 10, height = 10)
corrplot(cor_matrix,
         method    = "color",
         type      = "upper",
         order     = "hclust",
         addCoef.col = "black",
         tl.col    = "black",
         tl.srt    = 45,
         diag      = FALSE,
         title     = "Pairwise predictor correlations (April 2026 run)",
         mar       = c(0, 0, 2, 0))
invisible(dev.off())

png("../figures/supplement_figs/Figure_S2_pairwise_correlations.png",
    width = 10, height = 10, units = "in", res = 300)
corrplot(cor_matrix,
         method    = "color",
         type      = "upper",
         order     = "hclust",
         addCoef.col = "black",
         tl.col    = "black",
         tl.srt    = 45,
         diag      = FALSE,
         title     = "Pairwise predictor correlations (April 2026 run)",
         mar       = c(0, 0, 2, 0))
invisible(dev.off())

cat("Wrote Figure_S2_pairwise_correlations.{png,pdf}\n")
cat("Predictor count: ", ncol(predictors),
    "; high |r| > 0.7 pairs:\n", sep = "")
hi <- which(abs(cor_matrix) > 0.7 & abs(cor_matrix) < 1, arr.ind = TRUE)
for (i in seq_len(nrow(hi))) {
  if (hi[i, 1] < hi[i, 2]) {
    cat(sprintf("  %-20s ~ %-20s  r = %+.3f\n",
                rownames(cor_matrix)[hi[i, 1]],
                colnames(cor_matrix)[hi[i, 2]],
                cor_matrix[hi[i, 1], hi[i, 2]]))
  }
}
