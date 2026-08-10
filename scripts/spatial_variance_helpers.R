# Helpers for response-scale summaries of the fitted spatial intercept and
# spatial slope contributions. These functions operate on saved posterior draws;
# they do not refit the model.

as_plain_draw_matrix <- function(draws, variable) {
  out <- unclass(as_draws_matrix(subset_draws(draws, variable = variable)))
  attr(out, "draws") <- NULL
  out
}

row_variance <- function(x) {
  n <- ncol(x)
  if (n < 2L) stop("Variance calculation requires at least two observations")
  (rowSums(x * x) - rowSums(x)^2 / n) / (n - 1)
}

row_covariance <- function(x, y) {
  if (!identical(dim(x), dim(y))) {
    stop("Covariance inputs must have equal dimensions")
  }
  n <- ncol(x)
  if (n < 2L) stop("Covariance calculation requires at least two observations")
  (rowSums(x * y) - rowSums(x) * rowSums(y) / n) / (n - 1)
}

fitted_oipc_draw_matrix <- function(draws, stan_data, model = "model") {
  scale_weights <- as_plain_draw_matrix(draws, "scale_weights")
  n_scales <- stan_data$n_scales

  if (ncol(scale_weights) != n_scales) {
    stop("Scale-weight draws do not align with stan_data for ", model)
  }
  if (ncol(stan_data$oipc_values) != n_scales) {
    stop("OIPC predictor matrix does not align with stan_data for ", model)
  }

  scale_weight_sums <- rowSums(scale_weights)
  max_sum_error <- max(abs(scale_weight_sums - 1))
  if (!is.finite(max_sum_error) || max_sum_error > 1e-4) {
    stop("Saved scale weights differ materially from a simplex for ", model)
  }

  # Saved posterior values have finite decimal precision. After checking that
  # the discrepancy is numerical only, restore exact unit row sums before
  # reconstructing the predictor used by each draw.
  scale_weights <- scale_weights / scale_weight_sums
  fitted_oipc <- scale_weights %*% t(stan_data$oipc_values)
  attr(fitted_oipc, "scale_weight_max_sum_error") <- max_sum_error
  fitted_oipc
}

summarize_spatial_components <- function(draws, oipc_predictor,
                                         model = "model") {
  alpha <- as_plain_draw_matrix(draws, "alpha_spatial")
  local_slope <- as_plain_draw_matrix(draws, "beta_oipc_spatial")
  global_slope <- as.numeric(
    as_draws_matrix(subset_draws(draws, variable = "beta_oipc"))
  )
  residual_sd <- as.numeric(
    as_draws_matrix(subset_draws(draws, variable = "sigma"))
  )

  if (!identical(dim(alpha), dim(local_slope))) {
    stop("Spatial intercept and slope arrays do not align for ", model)
  }
  if (length(global_slope) != nrow(alpha) ||
      length(residual_sd) != nrow(alpha)) {
    stop("Scalar posterior draws do not align with spatial arrays for ", model)
  }

  slope_deviation <- sweep(local_slope, 1, global_slope, FUN = "-")
  if (is.null(dim(oipc_predictor))) {
    if (length(oipc_predictor) != ncol(alpha)) {
      stop("OIPC predictor vector does not align with observations for ", model)
    }
    slope_contribution <- sweep(
      slope_deviation, 2, oipc_predictor, FUN = "*"
    )
  } else {
    if (!identical(dim(oipc_predictor), dim(alpha))) {
      stop("OIPC predictor draws do not align with spatial arrays for ", model)
    }
    slope_contribution <- slope_deviation * oipc_predictor
  }

  var_intercept <- row_variance(alpha)
  var_slope <- row_variance(slope_contribution)
  covariance <- row_covariance(alpha, slope_contribution)
  var_spatial_direct <- row_variance(alpha + slope_contribution)
  var_spatial_identity <- var_intercept + var_slope + 2 * covariance
  identity_error <- max(abs(var_spatial_direct - var_spatial_identity))
  if (!is.finite(identity_error) || identity_error > 1e-8) {
    stop("Variance/covariance identity check failed for ", model)
  }

  mean_intercept <- mean(var_intercept)
  mean_slope <- mean(var_slope)
  mean_twice_covariance <- 2 * mean(covariance)
  mean_spatial_realized <- mean(var_spatial_direct)
  mean_residual <- mean(residual_sd^2)
  marginal_sum <- mean_intercept + mean_slope

  data.frame(
    n_draws = nrow(alpha),
    n_obs = ncol(alpha),
    marginal_intercept_variance_std2 = mean_intercept,
    marginal_slope_variance_std2 = mean_slope,
    twice_covariance_std2 = mean_twice_covariance,
    realized_spatial_variance_std2 = mean_spatial_realized,
    residual_variance_std2 = mean_residual,
    marginal_intercept_share_pct = 100 * mean_intercept / marginal_sum,
    marginal_slope_share_pct = 100 * mean_slope / marginal_sum,
    realized_spatial_share_of_spatial_plus_residual_pct =
      100 * mean_spatial_realized / (mean_spatial_realized + mean_residual),
    identity_max_abs_error = identity_error,
    check.names = FALSE
  )
}
