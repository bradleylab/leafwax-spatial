#───────────────────────────────────────────────────────────────────────────────
# 4a_spatial_functions.R
# Helper functions for spatial data preparation AND validation
# UPDATED: B-splines, consistent great circle distances, coordinate scaling
#───────────────────────────────────────────────────────────────────────────────

library(fields)
library(geosphere)
library(splines)

# Load configuration
source("0_load_config.R")

#───────────────────────────────────────────────────────────────────────────────
# SPATIAL FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────

# Gaussian weighting function based on distance
gaussian_weights <- function(dists, scale) {
  exp(- (dists^2) / (2 * scale^2))
}

# Function to convert degree distances to km using Haversine formula
deg_to_km <- function(distances_deg, center_lat) {
  # More accurate conversion using latitude-dependent scaling
  km_per_deg_lat <- 111.132
  km_per_deg_lon <- 111.132 * cos(center_lat * pi / 180)
  # For small distances, approximate as average
  avg_km_per_deg <- (km_per_deg_lat + km_per_deg_lon) / 2
  return(distances_deg * avg_km_per_deg)
}

# Function to compute weighted mean with safe defaults
compute_weighted_mean <- function(values, distances_deg, scale_km, center_lat) {
  if (length(values) == 0 || all(is.na(values))) {
    return(NA_real_)
  }
  valid <- !is.na(values)
  if (sum(valid) == 0) {
    return(NA_real_)
  }
  vals <- values[valid]
  dists_deg <- distances_deg[valid]
  
  # Convert degrees to km
  dists_km <- deg_to_km(dists_deg, center_lat)
  
  # Gaussian weights
  weights <- exp(-(dists_km^2) / (2 * scale_km^2))
  weight_sum <- sum(weights)
  
  # Return NA if no neighbors within reasonable distance
  if (weight_sum < 1e-10) {
    return(NA_real_)
  }
  return(sum(weights * vals) / weight_sum)
}

# Function to compute weighted interaction between two variables
compute_weighted_interaction <- function(values1, values2, distances_deg, scale_km, center_lat) {
  # Handle NULL inputs
  if (is.null(values1) || is.null(values2) || is.null(distances_deg)) {
    return(NA_real_)
  }
  
  # Handle empty or all-NA inputs
  if (length(values1) == 0 || length(values2) == 0 || 
      length(distances_deg) == 0 ||
      all(is.na(values1)) || all(is.na(values2))) {
    return(NA_real_)
  }
  
  # Ensure equal lengths
  if (length(values1) != length(values2) || length(values1) != length(distances_deg)) {
    warning("Input vectors have different lengths in compute_weighted_interaction")
    return(NA_real_)
  }
  
  # Get valid pairs where both variables are non-NA
  valid <- !is.na(values1) & !is.na(values2) & !is.na(distances_deg)
  if (sum(valid) == 0) {
    return(NA_real_)
  }
  
  # Compute interaction at each location
  interaction_values <- values1[valid] * values2[valid]
  dists_deg <- distances_deg[valid]
  
  # Convert degrees to km
  dists_km <- deg_to_km(dists_deg, center_lat)
  
  # Gaussian weights
  weights <- exp(-(dists_km^2) / (2 * scale_km^2))
  weight_sum <- sum(weights)
  
  # Return NA if no neighbors within reasonable distance
  if (weight_sum < 1e-10) {
    return(NA_real_)
  }
  
  return(sum(weights * interaction_values) / weight_sum)
}

# Convert distances in degrees to great circle distances in km
compute_great_circle_distances <- function(lon1, lat1, lon2, lat2) {
  # Use geosphere for accurate great circle distances
  # Returns distances in meters, so divide by 1000 for km
  if (length(lon2) == 1) {
    # Single point to single point
    return(distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000)
  } else {
    # Single point to multiple points
    p1 <- matrix(c(lon1, lat1), nrow = 1)
    p2 <- cbind(lon2, lat2)
    return(distHaversine(p1, p2) / 1000)
  }
}

# Compute spatially weighted averages of predictors at multiple distance scales
spatial_average_at_scales <- function(predictors_df, distance_matrix_km, scales_km) {
  stopifnot(nrow(predictors_df) == nrow(distance_matrix_km))
  
  lapply(scales_km, function(scale) {
    weights <- gaussian_weights(distance_matrix_km, scale)
    diag(weights) <- 0  # exclude self-influence
    weights <- weights / rowSums(weights)  # normalize weights
    
    weighted_avg <- weights %*% as.matrix(predictors_df)
    colnames(weighted_avg) <- paste0(colnames(predictors_df), "_", scale, "km")
    as.data.frame(weighted_avg)
  })
}

# Generate B-spline basis matrix
generate_bspline_basis <- function(x, n_knots, degree = 3, boundary_knots = NULL) {
  # x: vector of values to evaluate spline at
  # n_knots: number of interior knots
  # degree: degree of B-spline (3 = cubic)
  # boundary_knots: optional 2-element vector for boundary knots
  
  if (is.null(boundary_knots)) {
    boundary_knots <- range(x, na.rm = TRUE)
  }
  
  # Create interior knots
  interior_knots <- seq(boundary_knots[1], boundary_knots[2], 
                       length.out = n_knots + 2)[2:(n_knots + 1)]
  
  # Generate B-spline basis
  basis <- bs(x, knots = interior_knots, degree = degree, 
              Boundary.knots = boundary_knots, intercept = TRUE)
  
  return(basis)
}

# Calculate minimum distance from each point to nearest observation
calculate_min_dist_to_data <- function(coords_all, coords_obs) {
  # coords_all: matrix of all coordinates (N x 2)
  # coords_obs: matrix of observation coordinates (N_obs x 2)
  
  n_all <- nrow(coords_all)
  min_dists <- numeric(n_all)
  
  # For each point, find minimum distance to any observation
  for (i in 1:n_all) {
    # Calculate great circle distances from point i to all observations
    dists <- compute_great_circle_distances(
      coords_all[i,1], coords_all[i,2],
      coords_obs[,1], coords_obs[,2]
    )
    
    # If this point is an observation, exclude self-distance
    self_match <- which(coords_obs[,1] == coords_all[i,1] & 
                       coords_obs[,2] == coords_all[i,2])
    if (length(self_match) > 0) {
      dists[self_match] <- Inf
    }
    
    min_dists[i] <- min(dists)
  }
  
  # Handle case where a point matches itself (set to small value instead of Inf)
  min_dists[is.infinite(min_dists)] <- 0
  
  return(min_dists)
}

# Calculate data density at each knot with exponential regularization
# Based on Lindgren et al. (2011) JRSS-B
calculate_knot_data_density <- function(knot_coords, obs_coords, radius = 0.2, verbose = TRUE) {
  # knot_coords: matrix of knot coordinates (n_knots x 2) in standardized space
  # obs_coords: matrix of observation coordinates (n_obs x 2) in standardized space
  # radius: radius in standardized units for counting nearby observations
  
  n_knots <- nrow(knot_coords)
  n_obs <- nrow(obs_coords)
  knot_density <- numeric(n_knots)
  
  if (verbose) {
    cat("  Calculating data density at", n_knots, "knots...\n")
    cat("  Using radius of", radius, "standardized units for 'nearby' observations\n")
    cat("  Note: Exponential regularization will be applied (Lindgren et al., 2011)\n")
  }
  
  # For each knot, count observations within radius
  for (k in 1:n_knots) {
    # Calculate Euclidean distances in standardized space
    dists <- sqrt((knot_coords[k,1] - obs_coords[,1])^2 + 
                  (knot_coords[k,2] - obs_coords[,2])^2)
    
    # Count observations within radius
    knot_density[k] <- sum(dists <= radius)
  }
  
  if (verbose) {
    # Summary statistics
    cat("    Density range: [", min(knot_density), 
        ", ", max(knot_density), "] observations per knot\n")
    cat("    Mean density:", round(mean(knot_density), 1), "\n")
    cat("    Median density:", median(knot_density), "\n")
    
    # Show tau values that will result from exponential regularization
    tau_base <- 0.5
    tau_values <- numeric(n_knots)
    for (k in 1:n_knots) {
      density_factor <- exp(-knot_density[k] / CONFIG$gp_regularization$density_scaling)
      tau_values[k] <- tau_base * (0.5 + 0.5 * density_factor)
    }
    cat("\n    Resulting tau values (Lindgren et al., 2011):\n")
    cat("    Tau range: [", round(min(tau_values), 3), 
        ", ", round(max(tau_values), 3), "]\n")
  }
  
  return(knot_density)
}

# Predictive Process helper function
select_pp_knots <- function(coords, n_knots, method = "cover.design", min_dist_quantile = 0.05) {
  coords_unique <- unique(coords)
  
  if (nrow(coords_unique) < nrow(coords)) {
    cat("  Note:", nrow(coords) - nrow(coords_unique), "duplicate locations removed\n")
  }
  
  if (nrow(coords_unique) < n_knots) {
    warning(paste("Only", nrow(coords_unique), "unique locations available.",
                  "Using all unique locations as knots."))
    n_knots <- nrow(coords_unique)
    knot_coords <- coords_unique
  } else {
    if (method == "kmeans") {
      set.seed(123)
      km <- kmeans(coords_unique, centers = n_knots, nstart = 25)
      knot_coords <- km$centers
    } else if (method == "cover.design") {
      # First, ensure minimum spacing between candidate points
      dists <- as.matrix(rdist.earth(coords_unique, miles = FALSE))
      diag(dists) <- NA
      min_dist <- quantile(dists[upper.tri(dists)], min_dist_quantile, na.rm = TRUE)
      
      # Iteratively remove points that are too close
      keep <- rep(TRUE, nrow(coords_unique))
      for (i in 1:(nrow(coords_unique)-1)) {
        if (keep[i]) {
          too_close <- which(dists[i, ] < min_dist & keep)
          if (length(too_close) > 0) {
            keep[too_close] <- FALSE
          }
        }
      }
      
      coords_filtered <- coords_unique[keep, ]
      cat("  Filtered to", nrow(coords_filtered), "well-separated locations\n")
      
      # Now apply cover.design to the filtered coordinates
      if (nrow(coords_filtered) < n_knots) {
        cat("  Using all", nrow(coords_filtered), "filtered locations as knots\n")
        knot_coords <- coords_filtered
      } else {
        knot_coords <- cover.design(coords_filtered, n_knots)$design
      }
    } else if (method == "regular_global") {
      # Create deterministic equal-area global grid using Fibonacci sphere
      golden_angle <- pi * (3.0 - sqrt(5.0))
      
      # Generate exactly n_knots points
      knot_coords <- matrix(NA, n_knots, 2)
      
      for (i in 1:n_knots) {
        # Fibonacci sphere algorithm - deterministic for given n_knots
        theta <- golden_angle * (i - 1)
        z <- 1 - 2 * (i - 0.5) / n_knots
        radius <- sqrt(1 - z^2)
        
        lat <- asin(z) * 180 / pi
        lon <- (theta %% (2 * pi)) * 180 / pi - 180
        
        knot_coords[i, ] <- c(lon, lat)
      }
      
      cat("  Created deterministic equal-area global grid with", n_knots, "knots\n")
      
      # Calculate actual spacing
      if (n_knots > 1) {
        all_dists <- fields::rdist.earth(knot_coords, knot_coords, miles = FALSE)
        diag(all_dists) <- NA
        mean_min_dist <- mean(apply(all_dists, 1, min, na.rm = TRUE))
        cat("  Mean nearest-neighbor distance:", round(mean_min_dist), "km\n")
      }
    } else {
      stop("Unknown knot selection method: ", method)
    }
  }
  
  # Ensure proper column names
  if (!is.matrix(knot_coords)) {
    knot_coords <- as.matrix(knot_coords)
  }
  colnames(knot_coords) <- c("longitude", "latitude")
  
  return(knot_coords)
}

#───────────────────────────────────────────────────────────────────────────────
# CONTINENT FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────
# Function to add continent indicators to data
add_continent_indicators <- function(data) {
  data %>%
    mutate(
      # Define continents based on coordinates
      continent = case_when(
        longitude > 60 & longitude < 150 & latitude > -10 & latitude < 60 ~ "Asia",
        longitude > -20 & longitude < 60 & latitude > 35 ~ "Europe", 
        longitude > -20 & longitude < 60 & latitude < 35 ~ "Africa",
        longitude < -30 & longitude > -170 & latitude > 15 ~ "North America",
        longitude < -30 & longitude > -90 & latitude < 15 ~ "South America",
        longitude > 110 & latitude < -10 ~ "Australia",
        TRUE ~ "Other"
      ),
      # Create binary indicators
      is_asia = as.integer(continent == "Asia"),
      is_africa = as.integer(continent == "Africa"),
      is_europe = as.integer(continent == "Europe"),
      is_north_america = as.integer(continent == "North America"),
      is_south_america = as.integer(continent == "South America"),
      is_australia = as.integer(continent == "Australia")
    )
}

#───────────────────────────────────────────────────────────────────────────────
# VALIDATION FUNCTIONS
#───────────────────────────────────────────────────────────────────────────────

# Check kernel matrix condition and fix if needed
check_and_fix_kernel <- function(K, matrix_name = "kernel", 
                                 jitter = 1e-4, 
                                 max_cond = 1e6) {
  
  cat("Checking", matrix_name, "matrix...\n")
  
  # Ensure symmetry first
  K_sym <- 0.5 * (K + t(K))
  
  # Check initial condition number
  eigen_vals <- eigen(K_sym, only.values = TRUE)$values
  min_eval <- min(eigen_vals)
  max_eval <- max(eigen_vals)
  cond_num <- max_eval / max(abs(min_eval), 1e-10)
  
  cat("  Eigenvalue range: [", min_eval, ",", max_eval, "]\n")
  cat("  Condition number:", cond_num, "\n")
  
  # Check if matrix needs fixing
  needs_fixing <- (min_eval < jitter || cond_num > max_cond)
  
  if (!needs_fixing) {
    cat("  ✓ Matrix is well-conditioned\n")
    # Verify Cholesky works
    tryCatch({
      L <- chol(K_sym)
      cat("  ✓ Cholesky decomposition successful\n")
    }, error = function(e) {
      stop(matrix_name, " matrix is not positive definite")
    })
    return(K_sym)
  }
  
  # Add jitter iteratively if needed
  cat("  Matrix needs conditioning (min eval:", min_eval, "< jitter:", jitter, 
      "or cond:", cond_num, "> max:", max_cond, ")\n")
  
  jitter_added <- as.numeric(jitter)
  n_attempts <- 0
  max_attempts <- 10
  
  while (n_attempts < max_attempts) {
    # Add jitter
    n_dim <- nrow(K_sym)
    K_sym <- K_sym + diag(rep(jitter_added, n_dim))
    K_sym <- 0.5 * (K_sym + t(K_sym))
    
    # Recompute condition
    eigen_vals <- eigen(K_sym, only.values = TRUE)$values
    min_eval <- min(eigen_vals)
    max_eval <- max(eigen_vals)
    cond_num <- max_eval / max(abs(min_eval), 1e-10)
    
    cat("  Added jitter:", jitter_added, "- new condition number:", cond_num, 
        ", min eval:", min_eval, "\n")
    
    # Check if matrix is now acceptable
    is_well_conditioned <- (min_eval > 1e-10 && cond_num < max_cond)
    
    if (is_well_conditioned) {
      cat("  ✓ Matrix successfully conditioned\n")
      break
    }
    
    # Increase jitter for next iteration
    if (cond_num > 1e8) {
      jitter_added <- jitter_added * 10
    } else {
      jitter_added <- jitter_added * 2
    }
    
    n_attempts <- n_attempts + 1
  }
  
  # Final check
  if (n_attempts > 0 && (cond_num > max_cond || min_eval < 1e-10)) {
    warning(matrix_name, " matrix could not be adequately conditioned after ", 
            n_attempts, " attempts. Final condition number: ", round(cond_num), 
            ", min eigenvalue: ", min_eval)
  }
  
  # Verify Cholesky works
  tryCatch({
    L <- chol(K_sym)
    cat("  ✓ Cholesky decomposition successful\n")
  }, error = function(e) {
    stop(matrix_name, " matrix is not positive definite after conditioning")
  })
  
  return(K_sym)
}

# Validate stan_data before fitting
validate_stan_data <- function(stan_data) {
  cat("\nValidating Stan data...\n")
  
  # Check dimensions
  cat("Checking dimensions...\n")
  expected_dims <- list(
    oipc_values = c(stan_data$N, stan_data$n_scales),
    elevation_values = c(stan_data$N, stan_data$n_scales),
    c4_values_filled = c(stan_data$N, stan_data$n_scales)
  )
  
  # Add climate matrices if they should exist
  if (stan_data$include_precip == 1) {
    expected_dims$annual_precip <- c(stan_data$N, stan_data$n_scales)
  }
  if (stan_data$include_temp == 1) {
    expected_dims$max_temp <- c(stan_data$N, stan_data$n_scales)
  }
  if (stan_data$include_vpd == 1) {
    expected_dims$vpd <- c(stan_data$N, stan_data$n_scales)
  }
  if (stan_data$include_soil == 1) {
    expected_dims$soil_moisture <- c(stan_data$N, stan_data$n_scales)
  }
  
  for (var_name in names(expected_dims)) {
    if (var_name %in% names(stan_data)) {
      actual_dim <- dim(stan_data[[var_name]])
      expected_dim <- expected_dims[[var_name]]
      if (!all(actual_dim == expected_dim)) {
        stop(var_name, " has wrong dimensions: ", 
             paste(actual_dim, collapse = "x"), " instead of ",
             paste(expected_dim, collapse = "x"))
      }
    }
  }
  
  # Check for non-finite values
  cat("Checking for non-finite values...\n")
  matrix_vars <- c("oipc_values", "elevation_values", "c4_values_filled",
                   "pft_tree", "pft_shrub", "pft_grass",
                   "oipc_x_c4", "oipc_x_tree", "oipc_x_shrub", "oipc_x_grass",
                   "annual_precip", "soil_moisture", "max_temp", "vpd",
                   "temp_x_c4", "vpd_x_c4", "precip_x_tree", "precip_x_grass")
  
  for (var_name in matrix_vars) {
    if (var_name %in% names(stan_data) && is.matrix(stan_data[[var_name]])) {
      n_na <- sum(is.na(stan_data[[var_name]]))
      n_inf <- sum(is.infinite(stan_data[[var_name]]))
      n_nan <- sum(is.nan(stan_data[[var_name]]))
      
      if (n_na + n_inf + n_nan > 0) {
        stop(var_name, " contains non-finite values: ",
             n_na, " NA, ", n_inf, " Inf, ", n_nan, " NaN")
      }
    }
  }
  
  # Check required elements
  cat("Checking required elements...\n")
  required_elements <- c("N", "n_scales", "distance_scales", 
                        "d2H_wax", "d2H_wax_err",
                        "longitude", "latitude",
                        "oipc_values", "elevation_values",
                        "n_basis_knots", "spline_degree",
                        "include_c4", "include_pft", "include_gp",
                        "include_elevation", "include_precip", "include_temp",
                        "include_vpd", "include_soil",
                        "estimate_lambda", "lambda_fixed", "lambda_prior_mean_log",
                        "lambda_prior_sd_log", "use_pc_prior_intercept", "pc_prior_intercept_u",
                         "pc_prior_intercept_alpha", "use_pc_prior_slope", "pc_prior_slope_u", 
                         "pc_prior_slope_alpha")
  
  # Check for coordinates if GP is included
  if (stan_data$include_gp == 1) {
    required_elements <- c(required_elements, "coords", "knot_coords", "n_pp_knots", 
                          "knot_data_density", "density_scaling", "coord_scaling")
  }
  
  # Check for interaction matrices
  interaction_elements <- c("oipc_x_c4", "oipc_x_tree", "oipc_x_shrub", "oipc_x_grass")
  required_elements <- c(required_elements, interaction_elements)
  
  # Add climate interactions if needed
  if (stan_data$include_temp == 1 && stan_data$include_c4 == 1) {
    required_elements <- c(required_elements, "temp_x_c4")
  }
  if (stan_data$include_vpd == 1 && stan_data$include_c4 == 1) {
    required_elements <- c(required_elements, "vpd_x_c4")
  }
  if (stan_data$include_precip == 1 && stan_data$include_pft == 1) {
    required_elements <- c(required_elements, "precip_x_tree", "precip_x_grass")
  }
  
  missing <- setdiff(required_elements, names(stan_data))
  if (length(missing) > 0) {
    stop("Missing required stan_data elements: ", paste(missing, collapse = ", "))
  }
  
  # Check GP configuration if applicable
  if (stan_data$include_gp == 1) {
    cat("Checking GP configuration...\n")
    
    # Check coordinate dimensions
    if (!all(dim(stan_data$coords) == c(stan_data$N, 2))) {
      stop("coords has wrong dimensions")
    }
    
    if (!all(dim(stan_data$knot_coords) == c(stan_data$n_pp_knots, 2))) {
      stop("knot_coords has wrong dimensions")
    }
    
    # Check coordinate scaling
    if (length(stan_data$coord_scaling) != 2) {
      stop("coord_scaling must have 2 elements (lon_sd, lat_sd)")
    }
    
    # Check knot data density
    cat("Checking knot data density for regularizing prior...\n")
    if (length(stan_data$knot_data_density) != stan_data$n_pp_knots) {
      stop("knot_data_density has wrong length")
    }
    if (any(stan_data$knot_data_density < 0)) {
      stop("knot_data_density contains negative values")
    }
    cat("  Knot density range: [", min(stan_data$knot_data_density),
        ", ", max(stan_data$knot_data_density), "] observations\n")
  }
  
  # Check B-spline configuration if elevation included
  if (stan_data$include_elevation == 1) {
    cat("Checking B-spline configuration...\n")
    if (!("elevation_bspline_matrix" %in% names(stan_data))) {
      stop("elevation_bspline_matrix missing for B-spline model")
    }
    expected_cols <- stan_data$n_basis_knots + stan_data$spline_degree + 1
    if (ncol(stan_data$elevation_bspline_matrix) != expected_cols) {
      stop("elevation_bspline_matrix has wrong number of columns")
    }
  }
  
  cat("✓ All validation checks passed\n\n")
  return(stan_data)
}

#───────────────────────────────────────────────────────────────────────────────
# MAIN DATA PREPARATION FUNCTION
#───────────────────────────────────────────────────────────────────────────────
# Main data preparation function with spatial weighting
prepare_stan_data <- function(data, include_c4 = TRUE, include_pft = TRUE, include_gp = TRUE,
                              include_elevation = TRUE, include_precip = FALSE, include_temp = FALSE,
                              include_vpd = FALSE, include_soil = FALSE,
                              n_pp_knots = 100, SCALING_PARAMS, has_pft_columns) {
  
  N <- nrow(data)
  
  # Use config for distance scales
  distance_scales <- CONFIG$spatial_scales
  n_scales <- length(distance_scales)
  
  cat("Pre-computing spatially weighted values at", n_scales, "distance scales:\n")
  cat("  Scales (km):", paste(distance_scales, collapse = ", "), "\n")
  cat("  Note: Scale weights will be fitted by the model\n\n")
  
  # STANDARDIZE THE DATA FIRST
  cat("Standardizing data...\n")
  
  # Safe mean function that returns default if all NA
  safe_mean <- function(x, default = 0) {
    if (all(is.na(x))) {
      return(default)
    } else {
      return(mean(x, na.rm = TRUE))
    }
  }
  
  # Safe SD function that returns 1 if all NA or only one unique value
  safe_sd <- function(x, default = 1) {
    if (all(is.na(x))) {
      return(default)
    }
    x_clean <- x[!is.na(x)]
    if (length(unique(x_clean)) <= 1) {
      return(default)
    }
    return(sd(x, na.rm = TRUE))
  }
  
  # Create standardized versions with safe operations
  data_std <- data %>%
    mutate(
      d2H_wax_std = (d2H_wax - SCALING_PARAMS$d2H_mean) / SCALING_PARAMS$d2H_sd,
      d2H_wax_err_std = d2H_wax_err / SCALING_PARAMS$d2H_sd,
    )
  
  # Initialize storage matrices
  c4_weighted_matrix <- matrix(NA, N, n_scales)
  oipc_weighted_matrix <- matrix(NA, N, n_scales)
  oipc_se_weighted_matrix <- matrix(NA, N, n_scales)
  elev_weighted_matrix <- matrix(NA, N, n_scales)
  pft_weighted_array <- array(NA, dim = c(N, 3, n_scales))
  
  # Initialize interaction matrices
  oipc_x_c4_matrix <- matrix(NA, N, n_scales)
  oipc_x_tree_matrix <- matrix(NA, N, n_scales)
  oipc_x_shrub_matrix <- matrix(NA, N, n_scales)
  oipc_x_grass_matrix <- matrix(NA, N, n_scales)
  
  # Initialize climate matrices if any climate variable is requested
  include_any_climate <- include_precip || include_temp || include_vpd || include_soil
  
  if (include_any_climate) {
    cat("Pre-computing climate variable matrices...\n")
    
    # Initialize climate matrices
    precip_weighted_matrix <- matrix(NA, N, n_scales)
    soil_weighted_matrix <- matrix(NA, N, n_scales)
    temp_weighted_matrix <- matrix(NA, N, n_scales)
    vpd_weighted_matrix <- matrix(NA, N, n_scales)
    
    # Climate interactions
    temp_x_c4_matrix <- matrix(NA, N, n_scales)
    vpd_x_c4_matrix <- matrix(NA, N, n_scales)
    precip_x_tree_matrix <- matrix(NA, N, n_scales)
    precip_x_grass_matrix <- matrix(NA, N, n_scales)
  }
  
  # Progress tracking
  progress_interval <- max(1, floor(N / CONFIG$progress_interval_divisor))
  
  # Compute for each scale
  for (s in 1:n_scales) {
    scale_km <- distance_scales[s]
    cat("Computing for", scale_km, "km scale...\n")
    
    for (i in 1:N) {
      if (i %% progress_interval == 0) {
        cat("  Progress:", round(100 * i / N), "%\n")
      }
      
      # C4
      if (include_c4) {
        c4_weighted_matrix[i, s] <- compute_weighted_mean(
          data$c4_values_filled[[i]], 
          data$c4_distances[[i]], 
          scale_km,
          data$latitude[i]
        )
      }
      
      # OIPC
      oipc_weighted_matrix[i, s] <- compute_weighted_mean(
        data$oipc_values[[i]], 
        data$oipc_distances[[i]], 
        scale_km,
        data$latitude[i]
      )
      oipc_se_weighted_matrix[i, s] <- compute_weighted_mean(
        data$oipc_se_values[[i]], 
        data$oipc_distances[[i]], 
        scale_km,
        data$latitude[i]
      )
      
      # Elevation
      elev_weighted_matrix[i, s] <- compute_weighted_mean(
        data$elevation_values[[i]], 
        data$elevation_distances[[i]], 
        scale_km,
        data$latitude[i]
      )
      
      # PFT with safe operations
      if (include_pft && has_pft_columns && length(data$pft_distances[[i]]) > 0) {
        dists_km <- deg_to_km(data$pft_distances[[i]], data$latitude[i])
        weights <- exp(-(dists_km^2) / (2 * scale_km^2))
        weight_sum <- sum(weights)
        
        if (weight_sum > 1e-10) {
          if (!is.null(data$pft_tree[[i]]) && length(data$pft_tree[[i]]) > 0) {
            valid <- !is.na(data$pft_tree[[i]])
            if (sum(valid) > 0) {
              weighted_sum <- sum(weights[valid] * data$pft_tree[[i]][valid])
              pft_weighted_array[i, 1, s] <- weighted_sum / sum(weights[valid])
            }
          }
          if (!is.null(data$pft_shrub[[i]]) && length(data$pft_shrub[[i]]) > 0) {
            valid <- !is.na(data$pft_shrub[[i]])
            if (sum(valid) > 0) {
              weighted_sum <- sum(weights[valid] * data$pft_shrub[[i]][valid])
              pft_weighted_array[i, 2, s] <- weighted_sum / sum(weights[valid])
            }
          }
          if (!is.null(data$pft_grass[[i]]) && length(data$pft_grass[[i]]) > 0) {
            valid <- !is.na(data$pft_grass[[i]])
            if (sum(valid) > 0) {
              weighted_sum <- sum(weights[valid] * data$pft_grass[[i]][valid])
              pft_weighted_array[i, 3, s] <- weighted_sum / sum(weights[valid])
            }
          }
        }
      }
      
      # Compute interactions at each scale
      if (include_c4) {
        oipc_x_c4_matrix[i, s] <- compute_weighted_interaction(
          data$oipc_values[[i]], 
          data$c4_values_filled[[i]], 
          data$oipc_distances[[i]],  # Assuming same grid for OIPC and C4
          scale_km,
          data$latitude[i]
        )
      }
      
      if (include_pft && has_pft_columns && length(data$pft_distances[[i]]) > 0) {
        # OIPC × Tree
        if (!is.null(data$pft_tree[[i]]) && length(data$pft_tree[[i]]) > 0) {
          oipc_x_tree_matrix[i, s] <- compute_weighted_interaction(
            data$oipc_values[[i]], 
            data$pft_tree[[i]], 
            data$oipc_distances[[i]],  # Assuming aligned grids
            scale_km,
            data$latitude[i]
          )
        }
        
        # OIPC × Shrub
        if (!is.null(data$pft_shrub[[i]]) && length(data$pft_shrub[[i]]) > 0) {
          oipc_x_shrub_matrix[i, s] <- compute_weighted_interaction(
            data$oipc_values[[i]], 
            data$pft_shrub[[i]], 
            data$oipc_distances[[i]], 
            scale_km,
            data$latitude[i]
          )
        }
        
        # OIPC × Grass  
        if (!is.null(data$pft_grass[[i]]) && length(data$pft_grass[[i]]) > 0) {
          oipc_x_grass_matrix[i, s] <- compute_weighted_interaction(
            data$oipc_values[[i]], 
            data$pft_grass[[i]], 
            data$oipc_distances[[i]], 
            scale_km,
            data$latitude[i]
          )
        }
      }
      
      # Climate variables
      if (include_any_climate) {
        # Main climate effects
        if (include_precip) {
          precip_weighted_matrix[i, s] <- compute_weighted_mean(
            data$tc_ppt_values[[i]], 
            data$tc_ppt_distances[[i]], 
            scale_km,
            data$latitude[i]
          )
        }
        
        if (include_soil) {
          soil_weighted_matrix[i, s] <- compute_weighted_mean(
            data$tc_soil_values[[i]], 
            data$tc_soil_distances[[i]], 
            scale_km,
            data$latitude[i]
          )
        }
        
        if (include_temp) {
          temp_weighted_matrix[i, s] <- compute_weighted_mean(
            data$tc_tmax_values[[i]], 
            data$tc_tmax_distances[[i]], 
            scale_km,
            data$latitude[i]
          )
        }
        
        if (include_vpd) {
          vpd_weighted_matrix[i, s] <- compute_weighted_mean(
            data$tc_vpd_values[[i]], 
            data$tc_vpd_distances[[i]], 
            scale_km,
            data$latitude[i]
          )
        }
        
        # Climate interactions
        if (include_temp && include_c4) {
          temp_x_c4_matrix[i, s] <- compute_weighted_interaction(
            data$tc_tmax_values[[i]], 
            data$c4_values_filled[[i]], 
            data$tc_tmax_distances[[i]],
            scale_km,
            data$latitude[i]
          )
        }
        
        if (include_vpd && include_c4) {
          vpd_x_c4_matrix[i, s] <- compute_weighted_interaction(
            data$tc_vpd_values[[i]], 
            data$c4_values_filled[[i]], 
            data$tc_vpd_distances[[i]],
            scale_km,
            data$latitude[i]
          )
        }
        
        if (include_precip && include_pft && has_pft_columns && length(data$pft_distances[[i]]) > 0) {
          if (!is.null(data$pft_tree[[i]]) && length(data$pft_tree[[i]]) > 0) {
            precip_x_tree_matrix[i, s] <- compute_weighted_interaction(
              data$tc_ppt_values[[i]], 
              data$pft_tree[[i]], 
              data$tc_ppt_distances[[i]],
              scale_km,
              data$latitude[i]
            )
          }
          
          if (!is.null(data$pft_grass[[i]]) && length(data$pft_grass[[i]]) > 0) {
            precip_x_grass_matrix[i, s] <- compute_weighted_interaction(
              data$tc_ppt_values[[i]], 
              data$pft_grass[[i]], 
              data$tc_ppt_distances[[i]],
              scale_km,
              data$latitude[i]
            )
          }
        }
      }
    }
  }
  
  # STANDARDIZE THE AGGREGATED VALUES
  cat("\nStandardizing aggregated values...\n")
  
  # Standardize OIPC values
  oipc_weighted_matrix_std <- (oipc_weighted_matrix - SCALING_PARAMS$oipc_mean) / SCALING_PARAMS$oipc_sd
  oipc_se_weighted_matrix_std <- oipc_se_weighted_matrix / SCALING_PARAMS$oipc_sd
  
  # Standardize C4 values using config parameters
  if (include_c4) {
    c4_mean <- CONFIG$c4_standardization$mean
    c4_sd <- CONFIG$c4_standardization$sd
    c4_weighted_matrix_std <- (c4_weighted_matrix - c4_mean) / c4_sd
    cat("  C4 standardized with configured scaling - mean:", c4_mean, "%, sd:", c4_sd, "%\n")
  } else {
    c4_weighted_matrix_std <- matrix(0, N, n_scales)
  }
  
  # Standardize elevation values (convert to km first, then standardize)
  elev_weighted_matrix_km <- elev_weighted_matrix / 1000
  elev_mean_km <- SCALING_PARAMS$elev_mean / 1000
  elev_sd_km <- SCALING_PARAMS$elev_sd / 1000
  elev_weighted_matrix_std <- (elev_weighted_matrix_km - elev_mean_km) / elev_sd_km
  
  # Zero out elevation if not included in model
  if (!include_elevation) {
    cat("  Elevation excluded from model - zeroing elevation matrix\n")
    elev_weighted_matrix_std <- matrix(0, N, n_scales)
  }
  
  # Standardize climate variables if included
  if (include_any_climate) {
    # Get scaling parameters from data or config
    if (CONFIG$climate_standardization$compute_from_data) {
      # Compute from the site-level means
      if (include_precip) {
        SCALING_PARAMS$precip_mean <- safe_mean(data$annual_precip)
        SCALING_PARAMS$precip_sd <- safe_sd(data$annual_precip)
      }
      if (include_temp) {
        SCALING_PARAMS$temp_mean <- safe_mean(data$max_temp)
        SCALING_PARAMS$temp_sd <- safe_sd(data$max_temp)
      }
      if (include_vpd) {
        SCALING_PARAMS$vpd_mean <- safe_mean(data$vpd)
        SCALING_PARAMS$vpd_sd <- safe_sd(data$vpd)
      }
      if (include_soil) {
        SCALING_PARAMS$soil_mean <- safe_mean(data$soil_moisture)
        SCALING_PARAMS$soil_sd <- safe_sd(data$soil_moisture)
      }
    }
    
    # Standardize matrices
    if (include_precip) {
      precip_weighted_matrix_std <- (precip_weighted_matrix - SCALING_PARAMS$precip_mean) / SCALING_PARAMS$precip_sd
    } else {
      precip_weighted_matrix_std <- matrix(0, N, n_scales)
    }
    
    if (include_soil) {
      soil_weighted_matrix_std <- (soil_weighted_matrix - SCALING_PARAMS$soil_mean) / SCALING_PARAMS$soil_sd
    } else {
      soil_weighted_matrix_std <- matrix(0, N, n_scales)
    }
    
    if (include_temp) {
      temp_weighted_matrix_std <- (temp_weighted_matrix - SCALING_PARAMS$temp_mean) / SCALING_PARAMS$temp_sd
    } else {
      temp_weighted_matrix_std <- matrix(0, N, n_scales)
    }
    
    if (include_vpd) {
      vpd_weighted_matrix_std <- (vpd_weighted_matrix - SCALING_PARAMS$vpd_mean) / SCALING_PARAMS$vpd_sd
    } else {
      vpd_weighted_matrix_std <- matrix(0, N, n_scales)
    }
    
    # Standardize interactions
    if (include_temp && include_c4) {
      temp_x_c4_matrix_std <- temp_x_c4_matrix / (SCALING_PARAMS$temp_sd * CONFIG$c4_standardization$sd)
    } else {
      temp_x_c4_matrix_std <- matrix(0, N, n_scales)
    }
    
    if (include_vpd && include_c4) {
      vpd_x_c4_matrix_std <- vpd_x_c4_matrix / (SCALING_PARAMS$vpd_sd * CONFIG$c4_standardization$sd)
    } else {
      vpd_x_c4_matrix_std <- matrix(0, N, n_scales)
    }
    
    if (include_precip && include_pft) {
      precip_x_tree_matrix_std <- precip_x_tree_matrix / SCALING_PARAMS$precip_sd
      precip_x_grass_matrix_std <- precip_x_grass_matrix / SCALING_PARAMS$precip_sd
    } else {
      precip_x_tree_matrix_std <- matrix(0, N, n_scales)
      precip_x_grass_matrix_std <- matrix(0, N, n_scales)
    }
  } else {
    # Create empty matrices if no climate variables
    precip_weighted_matrix_std <- matrix(0, N, n_scales)
    soil_weighted_matrix_std <- matrix(0, N, n_scales)
    temp_weighted_matrix_std <- matrix(0, N, n_scales)
    vpd_weighted_matrix_std <- matrix(0, N, n_scales)
    temp_x_c4_matrix_std <- matrix(0, N, n_scales)
    vpd_x_c4_matrix_std <- matrix(0, N, n_scales)
    precip_x_tree_matrix_std <- matrix(0, N, n_scales)
    precip_x_grass_matrix_std <- matrix(0, N, n_scales)
  }
  
  # Standardize interaction matrices
  # Note: These are products of variables, so scale appropriately
  if (include_c4) {
    # OIPC×C4: product of two standardized variables
    oipc_x_c4_matrix_std <- oipc_x_c4_matrix / (SCALING_PARAMS$oipc_sd * CONFIG$c4_standardization$sd)
  } else {
    oipc_x_c4_matrix_std <- matrix(0, N, n_scales)
  }
  
  if (include_pft) {
    # OIPC×PFT: OIPC is standardized, PFT is proportion (0-1)
    oipc_x_tree_matrix_std <- oipc_x_tree_matrix / SCALING_PARAMS$oipc_sd
    oipc_x_shrub_matrix_std <- oipc_x_shrub_matrix / SCALING_PARAMS$oipc_sd
    oipc_x_grass_matrix_std <- oipc_x_grass_matrix / SCALING_PARAMS$oipc_sd
  } else {
    oipc_x_tree_matrix_std <- matrix(0, N, n_scales)
    oipc_x_shrub_matrix_std <- matrix(0, N, n_scales)
    oipc_x_grass_matrix_std <- matrix(0, N, n_scales)
  }
  
  # Handle missing values by interpolating across scales
  cat("\nHandling missing values...\n")
  
  for (i in 1:N) {
    # For each variable, fill NAs with the mean across scales for that observation
    # If all scales are NA, use 0 (population mean after standardization)
    
    # C4
    if (include_c4 && any(is.na(c4_weighted_matrix_std[i, ]))) {
      non_na <- c4_weighted_matrix_std[i, !is.na(c4_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      c4_weighted_matrix_std[i, is.na(c4_weighted_matrix_std[i, ])] <- fill_value
    }
    
    # OIPC
    if (any(is.na(oipc_weighted_matrix_std[i, ]))) {
      non_na <- oipc_weighted_matrix_std[i, !is.na(oipc_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      oipc_weighted_matrix_std[i, is.na(oipc_weighted_matrix_std[i, ])] <- fill_value
    }
    
    # OIPC SE
    if (any(is.na(oipc_se_weighted_matrix_std[i, ]))) {
      non_na <- oipc_se_weighted_matrix_std[i, !is.na(oipc_se_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0.1)
      oipc_se_weighted_matrix_std[i, is.na(oipc_se_weighted_matrix_std[i, ])] <- fill_value
    }
    
    # Elevation
    if (any(is.na(elev_weighted_matrix_std[i, ]))) {
      non_na <- elev_weighted_matrix_std[i, !is.na(elev_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      elev_weighted_matrix_std[i, is.na(elev_weighted_matrix_std[i, ])] <- fill_value
    }
    
    # PFT
    if (include_pft) {
      for (p in 1:3) {
        if (any(is.na(pft_weighted_array[i, p, ]))) {
          non_na <- pft_weighted_array[i, p, !is.na(pft_weighted_array[i, p, ])]
          fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0.33)
          pft_weighted_array[i, p, is.na(pft_weighted_array[i, p, ])] <- fill_value
        }
      }
    }
    
    # Climate variables
    if (include_precip && any(is.na(precip_weighted_matrix_std[i, ]))) {
      non_na <- precip_weighted_matrix_std[i, !is.na(precip_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      precip_weighted_matrix_std[i, is.na(precip_weighted_matrix_std[i, ])] <- fill_value
    }
    
    if (include_temp && any(is.na(temp_weighted_matrix_std[i, ]))) {
      non_na <- temp_weighted_matrix_std[i, !is.na(temp_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      temp_weighted_matrix_std[i, is.na(temp_weighted_matrix_std[i, ])] <- fill_value
    }
    
    if (include_vpd && any(is.na(vpd_weighted_matrix_std[i, ]))) {
      non_na <- vpd_weighted_matrix_std[i, !is.na(vpd_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      vpd_weighted_matrix_std[i, is.na(vpd_weighted_matrix_std[i, ])] <- fill_value
    }
    
    if (include_soil && any(is.na(soil_weighted_matrix_std[i, ]))) {
      non_na <- soil_weighted_matrix_std[i, !is.na(soil_weighted_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      soil_weighted_matrix_std[i, is.na(soil_weighted_matrix_std[i, ])] <- fill_value
    }
    
    # Handle missing values for interactions
    if (include_c4 && any(is.na(oipc_x_c4_matrix_std[i, ]))) {
      non_na <- oipc_x_c4_matrix_std[i, !is.na(oipc_x_c4_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      oipc_x_c4_matrix_std[i, is.na(oipc_x_c4_matrix_std[i, ])] <- fill_value
    }
    
    if (include_pft) {
      # Tree interaction
      if (any(is.na(oipc_x_tree_matrix_std[i, ]))) {
        non_na <- oipc_x_tree_matrix_std[i, !is.na(oipc_x_tree_matrix_std[i, ])]
        fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
        oipc_x_tree_matrix_std[i, is.na(oipc_x_tree_matrix_std[i, ])] <- fill_value
      }
      
      # Shrub interaction
      if (any(is.na(oipc_x_shrub_matrix_std[i, ]))) {
        non_na <- oipc_x_shrub_matrix_std[i, !is.na(oipc_x_shrub_matrix_std[i, ])]
        fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
        oipc_x_shrub_matrix_std[i, is.na(oipc_x_shrub_matrix_std[i, ])] <- fill_value
      }
      
      # Grass interaction
      if (any(is.na(oipc_x_grass_matrix_std[i, ]))) {
        non_na <- oipc_x_grass_matrix_std[i, !is.na(oipc_x_grass_matrix_std[i, ])]
        fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
        oipc_x_grass_matrix_std[i, is.na(oipc_x_grass_matrix_std[i, ])] <- fill_value
      }
    }
    
    # Climate interactions
    if (include_temp && include_c4 && any(is.na(temp_x_c4_matrix_std[i, ]))) {
      non_na <- temp_x_c4_matrix_std[i, !is.na(temp_x_c4_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      temp_x_c4_matrix_std[i, is.na(temp_x_c4_matrix_std[i, ])] <- fill_value
    }
    
    if (include_vpd && include_c4 && any(is.na(vpd_x_c4_matrix_std[i, ]))) {
      non_na <- vpd_x_c4_matrix_std[i, !is.na(vpd_x_c4_matrix_std[i, ])]
      fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
      vpd_x_c4_matrix_std[i, is.na(vpd_x_c4_matrix_std[i, ])] <- fill_value
    }
    
    if (include_precip && include_pft) {
      if (any(is.na(precip_x_tree_matrix_std[i, ]))) {
        non_na <- precip_x_tree_matrix_std[i, !is.na(precip_x_tree_matrix_std[i, ])]
        fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
        precip_x_tree_matrix_std[i, is.na(precip_x_tree_matrix_std[i, ])] <- fill_value
      }
      
      if (any(is.na(precip_x_grass_matrix_std[i, ]))) {
        non_na <- precip_x_grass_matrix_std[i, !is.na(precip_x_grass_matrix_std[i, ])]
        fill_value <- ifelse(length(non_na) > 0, mean(non_na), 0)
        precip_x_grass_matrix_std[i, is.na(precip_x_grass_matrix_std[i, ])] <- fill_value
      }
    }
  }
  
  # Ensure minimum values for numerical stability
  cat("\nEnsuring minimum values for numerical stability...\n")
  
  # OIPC SE should never be exactly zero
  oipc_se_weighted_matrix_std[oipc_se_weighted_matrix_std < 1e-4] <- 1e-4
  
  # Replace any remaining NaN values with appropriate defaults
  oipc_weighted_matrix_std[is.nan(oipc_weighted_matrix_std)] <- 0
  oipc_se_weighted_matrix_std[is.nan(oipc_se_weighted_matrix_std)] <- 0.1
  elev_weighted_matrix_std[is.nan(elev_weighted_matrix_std)] <- 0
  if (include_c4) {
    c4_weighted_matrix_std[is.nan(c4_weighted_matrix_std)] <- 0
    oipc_x_c4_matrix_std[is.nan(oipc_x_c4_matrix_std)] <- 0
  }
  if (include_pft) {
    pft_weighted_array[is.nan(pft_weighted_array)] <- 0.33
    oipc_x_tree_matrix_std[is.nan(oipc_x_tree_matrix_std)] <- 0
    oipc_x_shrub_matrix_std[is.nan(oipc_x_shrub_matrix_std)] <- 0
    oipc_x_grass_matrix_std[is.nan(oipc_x_grass_matrix_std)] <- 0
  }
  if (include_any_climate) {
    precip_weighted_matrix_std[is.nan(precip_weighted_matrix_std)] <- 0
    temp_weighted_matrix_std[is.nan(temp_weighted_matrix_std)] <- 0
    vpd_weighted_matrix_std[is.nan(vpd_weighted_matrix_std)] <- 0
    soil_weighted_matrix_std[is.nan(soil_weighted_matrix_std)] <- 0
    temp_x_c4_matrix_std[is.nan(temp_x_c4_matrix_std)] <- 0
    vpd_x_c4_matrix_std[is.nan(vpd_x_c4_matrix_std)] <- 0
    precip_x_tree_matrix_std[is.nan(precip_x_tree_matrix_std)] <- 0
    precip_x_grass_matrix_std[is.nan(precip_x_grass_matrix_std)] <- 0
  }
  
  # Check for any remaining NAs or NaNs
  cat("\nChecking for NaN/NA values in matrices...\n")
  
  check_matrix <- function(mat, name) {
    n_na <- sum(is.na(mat))
    n_nan <- sum(is.nan(mat))
    n_inf <- sum(is.infinite(mat))
    if (n_na > 0 || n_nan > 0 || n_inf > 0) {
      cat("  WARNING:", name, "has", n_na, "NA,", n_nan, "NaN, and", n_inf, "Inf values\n")
      return(FALSE)
    }
    return(TRUE)
  }
  
  all_good <- TRUE
  all_good <- check_matrix(oipc_weighted_matrix_std, "OIPC") && all_good
  all_good <- check_matrix(oipc_se_weighted_matrix_std, "OIPC SE") && all_good
  all_good <- check_matrix(elev_weighted_matrix_std, "Elevation") && all_good
  if (include_c4) {
    all_good <- check_matrix(c4_weighted_matrix_std, "C4") && all_good
    all_good <- check_matrix(oipc_x_c4_matrix_std, "OIPC×C4") && all_good
  }
  if (include_pft) {
    for (p in 1:3) {
      pft_name <- c("Tree", "Shrub", "Grass")[p]
      all_good <- check_matrix(pft_weighted_array[, p, ], paste("PFT", pft_name)) && all_good
    }
    all_good <- check_matrix(oipc_x_tree_matrix_std, "OIPC×Tree") && all_good
    all_good <- check_matrix(oipc_x_shrub_matrix_std, "OIPC×Shrub") && all_good
    all_good <- check_matrix(oipc_x_grass_matrix_std, "OIPC×Grass") && all_good
  }
  if (include_any_climate) {
    if (include_precip) all_good <- check_matrix(precip_weighted_matrix_std, "Precipitation") && all_good
    if (include_temp) all_good <- check_matrix(temp_weighted_matrix_std, "Temperature") && all_good
    if (include_vpd) all_good <- check_matrix(vpd_weighted_matrix_std, "VPD") && all_good
    if (include_soil) all_good <- check_matrix(soil_weighted_matrix_std, "Soil Moisture") && all_good
    if (include_temp && include_c4) all_good <- check_matrix(temp_x_c4_matrix_std, "Temp×C4") && all_good
    if (include_vpd && include_c4) all_good <- check_matrix(vpd_x_c4_matrix_std, "VPD×C4") && all_good
    if (include_precip && include_pft) {
      all_good <- check_matrix(precip_x_tree_matrix_std, "Precip×Tree") && all_good
      all_good <- check_matrix(precip_x_grass_matrix_std, "Precip×Grass") && all_good
    }
  }
  
  if (all_good) {
    cat("  ✓ All matrices are clean (no NaN/NA/Inf values)\n")
  }
  
  # Generate B-spline basis for elevation if included
  elevation_bspline_matrix <- NULL
  n_basis_knots <- CONFIG$n_elevation_knots
  
  if (include_elevation) {
    cat("\nGenerating B-spline basis for elevation...\n")
    
    # Get unique elevation values across all scales
    all_elev_values <- as.vector(elev_weighted_matrix_std)
    elev_range <- range(all_elev_values[!is.infinite(all_elev_values)], na.rm = TRUE)
    
    # Handle edge case where all elevations are the same
    if (abs(elev_range[2] - elev_range[1]) < 1e-10) {
      cat("  WARNING: All elevations are identical - using constant basis\n")
      # Create a simple constant basis
      elevation_bspline_matrix <- matrix(1, N * n_scales, 1)
      n_basis_knots <- 0  # No interior knots
    } else {
      # Generate B-spline basis for each observation at each scale
      basis_list <- list()
      for (s in 1:n_scales) {
        basis_s <- generate_bspline_basis(
          elev_weighted_matrix_std[, s],
          n_knots = n_basis_knots,
          degree = CONFIG$elevation_spline_degree,
          boundary_knots = elev_range
        )
        basis_list[[s]] <- basis_s
      }
      
      # Stack all bases (will be reshaped in Stan)
      elevation_bspline_matrix <- do.call(rbind, basis_list)
      elevation_bspline_matrix <- scale(elevation_bspline_matrix, center = TRUE, scale = FALSE)

      
      cat("  Created B-spline basis with", n_basis_knots, "interior knots\n")
      cat("  Degree:", CONFIG$elevation_spline_degree, "\n")
      cat("  Elevation range (standardized):", round(elev_range[1], 3), "to", round(elev_range[2], 3), "\n")
      cat("  Basis dimensions:", nrow(elevation_bspline_matrix), "x", ncol(elevation_bspline_matrix), "\n")
    }
  } else {
    # IMPORTANT: When elevation is not included, we still need valid dimensions
    cat("\nElevation not included - creating dummy B-spline matrix\n")
    n_basis_knots <- 0  # Must be >= 0 for Stan constraint
    # Create a dummy matrix with the expected number of columns
    # n_basis_knots + spline_degree + 1 = 0 + 3 + 1 = 4
    n_cols <- n_basis_knots + CONFIG$elevation_spline_degree + 1
    elevation_bspline_matrix <- matrix(0, nrow = N * n_scales, ncol = n_cols)
  }
  
  # Add continent indicators
  data_with_continents <- add_continent_indicators(data)
  
  # Prepare coordinates and compute scaling
  coords_std <- matrix(0, N, 2)
  pp_knot_coords <- matrix(0, 1, 2)
  knot_data_density <- numeric(1)
  coord_scaling <- c(1, 1)  # Default
  
  if (include_gp) {
    cat("\nSelecting", n_pp_knots, "knots for Predictive Process approximation...\n")
    
    # Get standardized coordinates and scaling factors
    lon_mean <- safe_mean(data$longitude)
    lon_sd <- safe_sd(data$longitude)
    lat_mean <- safe_mean(data$latitude)
    lat_sd <- safe_sd(data$latitude)
    
    # Store scaling for length scale conversion
    coord_scaling <- c(lon_sd, lat_sd)
    
    coords_std <- cbind(
      (data$longitude - lon_mean) / lon_sd,
      (data$latitude - lat_mean) / lat_sd
    )
    
    # Select knots using configured method
    if (CONFIG$knot_selection_method == "regular_global") {
      # For global grid, we need to work in original coordinate space
      coords_orig <- cbind(data$longitude, data$latitude)
      
      # Get knots in original space
      knot_coords_orig <- select_pp_knots(coords_orig, n_pp_knots, 
                                          method = CONFIG$knot_selection_method,
                                          min_dist_quantile = CONFIG$min_dist_quantile)
      
      # Now standardize the knot coordinates
      pp_knot_coords <- cbind(
        (knot_coords_orig[,1] - lon_mean) / lon_sd,
        (knot_coords_orig[,2] - lat_mean) / lat_sd
      )
      colnames(pp_knot_coords) <- c("longitude", "latitude")
    } else {
      # For other methods, work in standardized space
      pp_knot_coords <- select_pp_knots(coords_std, n_pp_knots, 
                                        method = CONFIG$knot_selection_method,
                                        min_dist_quantile = CONFIG$min_dist_quantile)
    }
    
    cat("✓ Knots selected using space-filling design\n")
    
    # Calculate knot data density for regularizing prior
    cat("\nCalculating data density at knots for GP regularization...\n")
    
    # Determine radius based on config
    density_radius <- CONFIG$gp_regularization$density_radius_std
    density_scaling <- CONFIG$gp_regularization$density_scaling
    
    knot_data_density <- calculate_knot_data_density(pp_knot_coords, coords_std, 
                                                      radius = density_radius, 
                                                      verbose = TRUE)
  }
  
  # Prepare Stan data
  stan_data <- list(
    N = N,
    n_scales = n_scales,
    distance_scales = distance_scales,
    
    # Use standardized response
    d2H_wax = data_std$d2H_wax_std,
    d2H_wax_err = data_std$d2H_wax_err_std,
    
    # Coordinates for kernel computation
    coords = coords_std,
    coord_scaling = coord_scaling,  # NEW: for proper length scale conversion
    
    # Original coordinates (don't standardize these)
    longitude = data$longitude,
    latitude = data$latitude,
    
    # Standardized matrices - RENAMED TO MATCH STAN MODEL
    oipc_values = oipc_weighted_matrix_std,
    oipc_se_values = oipc_se_weighted_matrix_std,
    elevation_values = elev_weighted_matrix_std,
    c4_values_filled = if (include_c4) c4_weighted_matrix_std else matrix(0, N, n_scales),
    
    # PFT matrices - RENAMED TO MATCH STAN MODEL
    pft_tree = if (include_pft) pft_weighted_array[, 1, ] else matrix(0, N, n_scales),
    pft_shrub = if (include_pft) pft_weighted_array[, 2, ] else matrix(0, N, n_scales),
    pft_grass = if (include_pft) pft_weighted_array[, 3, ] else matrix(0, N, n_scales),
    
    # Climate matrices - MUST HAVE 0 COLUMNS WHEN NOT INCLUDED
    annual_precip = if (include_precip) precip_weighted_matrix_std else matrix(0, N, 0),
    soil_moisture = if (include_soil) soil_weighted_matrix_std else matrix(0, N, 0),
    max_temp = if (include_temp) temp_weighted_matrix_std else matrix(0, N, 0),
    vpd = if (include_vpd) vpd_weighted_matrix_std else matrix(0, N, 0),
    
    # Pre-computed interaction matrices
    oipc_x_c4 = if (include_c4) oipc_x_c4_matrix_std else matrix(0, N, n_scales),
    oipc_x_tree = if (include_pft) oipc_x_tree_matrix_std else matrix(0, N, n_scales),
    oipc_x_shrub = if (include_pft) oipc_x_shrub_matrix_std else matrix(0, N, n_scales),
    oipc_x_grass = if (include_pft) oipc_x_grass_matrix_std else matrix(0, N, n_scales),
    
    # Climate interactions
    temp_x_c4 = temp_x_c4_matrix_std,
    vpd_x_c4 = vpd_x_c4_matrix_std,
    precip_x_tree = precip_x_tree_matrix_std,
    precip_x_grass = precip_x_grass_matrix_std,
    
    # B-spline basis for elevation
    elevation_bspline_matrix = elevation_bspline_matrix,
    n_basis_knots = n_basis_knots,
    spline_degree = CONFIG$elevation_spline_degree,
    
    # GP configuration WITHOUT kernel matrices
    n_pp_knots = if (include_gp) n_pp_knots else 1,
    knot_coords = pp_knot_coords,
    knot_data_density = knot_data_density,
    density_scaling = if (include_gp) density_scaling else 1,
    
    # Model configuration
    include_c4 = as.integer(include_c4),
    include_pft = as.integer(include_pft),
    include_gp = as.integer(include_gp),
    include_elevation = as.integer(include_elevation),
    
    # Individual climate flags
    include_precip = as.integer(include_precip),
    include_temp = as.integer(include_temp),
    include_vpd = as.integer(include_vpd),
    include_soil = as.integer(include_soil),
    
    # Lambda control parameters - UPDATED for direct km interpretation
    estimate_lambda = as.integer(CONFIG$lambda_control$estimate),
    lambda_fixed = CONFIG$lambda_control$fixed_value,
    lambda_prior_mean_log = CONFIG$lambda_control$prior_mean_log,
    lambda_prior_sd_log = CONFIG$lambda_control$prior_sd_log,
    
    # PC prior parameters
    use_pc_prior_intercept = as.integer(CONFIG$pc_prior_intercept$use_pc_prior),
	pc_prior_intercept_u = CONFIG$pc_prior_intercept$u_permil / SCALING_PARAMS$d2H_sd,
	pc_prior_intercept_alpha = CONFIG$pc_prior_intercept$alpha,

	use_pc_prior_slope = as.integer(CONFIG$pc_prior_slope$use_pc_prior),
	pc_prior_slope_u = CONFIG$pc_prior_slope$u_value,  # No scaling - already unitless
	pc_prior_slope_alpha = CONFIG$pc_prior_slope$alpha,
    
    # Prior predictive checks
    d2H_wax_sd_original = SCALING_PARAMS$d2H_sd,
    d2H_wax_mean_original = SCALING_PARAMS$d2H_mean,
    
    # Store scaling parameters for back-transformation
    scaling_params = SCALING_PARAMS,
    
    # Continent indicators (for interaction effects)
    is_asia = data_with_continents$is_asia,
    is_africa = data_with_continents$is_africa,
    is_europe = data_with_continents$is_europe,
    is_north_america = data_with_continents$is_north_america,
    is_south_america = data_with_continents$is_south_america,
    is_australia = data_with_continents$is_australia
  )
  
  # Final validation
  cat("\nFinal validation of stan_data...\n")
  final_check_passed <- TRUE
  
  for (name in names(stan_data)) {
    if (is.numeric(stan_data[[name]]) || is.matrix(stan_data[[name]])) {
      n_na <- sum(is.na(stan_data[[name]]))
      n_nan <- sum(is.nan(stan_data[[name]]))
      n_inf <- sum(is.infinite(stan_data[[name]]))
      
      if (n_na > 0 || n_nan > 0 || n_inf > 0) {
        cat("  ERROR:", name, "contains problematic values!\n")
        cat("    NA count:", n_na, "\n")
        cat("    NaN count:", n_nan, "\n")
        cat("    Inf count:", n_inf, "\n")
        final_check_passed <- FALSE
      }
    }
  }
  
  if (final_check_passed) {
    cat("  ✓ All stan_data elements are clean\n")
  } else {
    cat("  ⚠ Some stan_data elements contain problematic values\n")
  }
  
  cat("\nPre-aggregation complete!\n")
  cat("  Observations:", N, "\n")
  cat("  Distance scales:", n_scales, "\n")
  cat("  B-spline knots:", n_basis_knots, "(interior)\n")
  cat("  Spline degree:", CONFIG$elevation_spline_degree, "\n")
  cat("  Data standardized: YES\n")
  cat("  Lambda control: ", ifelse(CONFIG$lambda_control$estimate, "ESTIMATED", 
                                    paste0("FIXED at ", CONFIG$lambda_control$fixed_value, " km")), "\n")
  if (include_gp) {
    cat("  Predictive Process knots:", n_pp_knots, "\n")
    cat("  GP kernel: Matérn 3/2\n")
    cat("  GP length scale: ESTIMATED FROM DATA\n")
    cat("  GP regularization: Density-based prior with scaling factor", density_scaling, "\n")
    cat("  Coordinate scaling (lon, lat SD):", round(coord_scaling, 3), "\n")
    cat("  PC prior: ", ifelse(CONFIG$pc_prior$use_pc_prior, 
                               paste0("YES (P(σ > ", CONFIG$pc_prior$u_permil, "‰) = ", CONFIG$pc_prior$alpha, ")"),
                               "NO (using half-normal)"), "\n")
  }
  cat("  Climate variables included:\n")
  if (include_precip) cat("    - Precipitation\n")
  if (include_temp) cat("    - Temperature\n")
  if (include_vpd) cat("    - VPD\n")
  if (include_soil) cat("    - Soil moisture\n")
  if (!include_any_climate) cat("    - None\n")
  cat("\n")
  
  return(stan_data)
}