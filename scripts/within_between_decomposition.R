# Decompose the pooled wax--precipitation isotope relationship into regional
# and within-region components using the model-ready calibration data.
#
# Run from any working directory:
#   Rscript scripts/within_between_decomposition.R
#
# Set LEAFWAX_OUTPUT_DIR to override the default generated-output directory.

suppressPackageStartupMessages({
  library(dplyr)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1L) {
  stop("Could not resolve this script's path from the Rscript command line")
}
script_path <- normalizePath(sub("^--file=", "", file_arg))
repo_root <- normalizePath(file.path(dirname(script_path), ".."))

data_path <- file.path(repo_root, "results", "3_sediment_ready_for_modeling.rds")
if (!file.exists(data_path)) {
  stop("Model-ready calibration data not found: ", data_path)
}

output_dir <- Sys.getenv(
  "LEAFWAX_OUTPUT_DIR",
  unset = file.path(repo_root, "model_analysis", "reported_outputs")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

calibration <- readRDS(data_path)
stopifnot(nrow(calibration) == 1128L)

calibration$region <- dplyr::case_when(
  calibration$longitude < -30 ~ "Americas",
  calibration$longitude >= -30 & calibration$longitude < 60 &
    calibration$latitude > 35 ~ "Europe",
  calibration$longitude >= -30 & calibration$longitude < 60 &
    calibration$latitude <= 35 ~ "Africa",
  calibration$longitude >= 60 & calibration$longitude < 140 &
    calibration$latitude > -10 ~ "Asia",
  calibration$longitude >= 140 | calibration$latitude < -10 ~ "Oceania"
)
stopifnot(!anyNA(calibration$region))

wax <- calibration$d2H_wax

decompose_relationship <- function(predictor, label) {
  keep <- is.finite(predictor) & is.finite(wax)
  x <- predictor[keep]
  y <- wax[keep]
  region <- calibration$region[keep]

  regional_means <- data.frame(x = x, y = y, region = region) |>
    group_by(region) |>
    summarise(
      n = dplyr::n(),
      x_mean = mean(x),
      y_mean = mean(y),
      .groups = "drop"
    )

  x_mean <- mean(x)
  y_mean <- mean(y)
  total_covariance <- mean((x - x_mean) * (y - y_mean))
  total_predictor_variance <- mean((x - x_mean)^2)
  weights <- regional_means$n / sum(regional_means$n)
  between_covariance <- sum(
    weights *
      (regional_means$x_mean - x_mean) *
      (regional_means$y_mean - y_mean)
  )
  between_predictor_variance <- sum(
    weights * (regional_means$x_mean - x_mean)^2
  )

  data.frame(
    predictor = label,
    n = length(y),
    global_slope = unname(coef(lm(y ~ x))[2]),
    between_slope_unweighted = unname(
      coef(lm(y_mean ~ x_mean, data = regional_means))[2]
    ),
    between_slope_nweighted = unname(
      coef(lm(y_mean ~ x_mean, data = regional_means,
              weights = n))[2]
    ),
    within_slope = unname(coef(lm(y ~ x + factor(region)))[2]),
    pct_predictor_variance_between =
      100 * between_predictor_variance / total_predictor_variance,
    pct_covariance_between = 100 * between_covariance / total_covariance,
    pct_covariance_within =
      100 * (total_covariance - between_covariance) / total_covariance
  )
}

decomposition <- rbind(
  decompose_relationship(calibration$d2H_precip, "OIPC point value"),
  decompose_relationship(calibration$oipc_mean, "OIPC multi-scale mean")
)

regional_slopes <- calibration |>
  filter(is.finite(oipc_mean), is.finite(d2H_wax)) |>
  group_by(region) |>
  summarise(
    n = dplyr::n(),
    slope = coef(lm(d2H_wax ~ oipc_mean))[2],
    r2 = summary(lm(d2H_wax ~ oipc_mean))$r.squared,
    d2Hprecip_min = min(oipc_mean),
    d2Hprecip_max = max(oipc_mean),
    .groups = "drop"
  ) |>
  arrange(region)

decomposition_path <- file.path(output_dir, "within_between_decomposition.csv")
regional_slopes_path <- file.path(output_dir, "regional_slopes.csv")
write.csv(decomposition, decomposition_path, row.names = FALSE)
write.csv(regional_slopes, regional_slopes_path, row.names = FALSE)

cat("Within/between decomposition written to:\n", decomposition_path, "\n")
cat("Regional slopes written to:\n", regional_slopes_path, "\n")
