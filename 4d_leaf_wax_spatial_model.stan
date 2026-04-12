// 4d_leaf_wax_spatial_model.stan
// Hierarchical model with spatially-varying OIPC relationship

functions {
  // Matérn 3/2 covariance kernel
  matrix cov_matern32(array[] vector x, real alpha, real rho) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    real sqrt3 = sqrt(3.0);
    
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i+1):N) {
        real dist = distance(x[i], x[j]);
        real scaled_dist = sqrt3 * dist / rho;
        K[i, j] = sq_alpha * (1 + scaled_dist) * exp(-scaled_dist);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
  
  // Cross-covariance for Matérn 3/2
  matrix cov_matern32_cross(array[] vector x1, array[] vector x2, real alpha, real rho) {
    int N1 = size(x1);
    int N2 = size(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    real sqrt3 = sqrt(3.0);
    
    for (i in 1:N1) {
      for (j in 1:N2) {
        real dist = distance(x1[i], x2[j]);
        real scaled_dist = sqrt3 * dist / rho;
        K[i, j] = sq_alpha * (1 + scaled_dist) * exp(-scaled_dist);
      }
    }
    return K;
  }
}

data {
  // Sample size
  int<lower=1> N;
  
  // Response variable (standardized)
  vector[N] d2H_wax;
  vector<lower=0>[N] d2H_wax_err;
  
  // Spatial coordinates (standardized)
  array[N] vector[2] coords;  // lon, lat for each observation
  vector[2] coord_scaling;    // SD of lon, lat for length scale conversion
  
  // Multi-scale predictors
  int<lower=1> n_scales;
  vector<lower=0>[n_scales] distance_scales;  // in km
  
  // Model configuration flags (MUST BE BEFORE conditional declarations)
  int<lower=0, upper=1> include_c4;
  int<lower=0, upper=1> include_pft;
  int<lower=0, upper=1> include_gp;
  int<lower=0, upper=1> include_elevation;
  int<lower=0, upper=1> include_precip;
  int<lower=0, upper=1> include_temp;
  int<lower=0, upper=1> include_vpd;
  int<lower=0, upper=1> include_soil;
  int<lower=0, upper=1> include_veg_interactions;
  
  // Core predictor matrices (N x n_scales)
  matrix[N, n_scales] oipc_values;
  matrix<lower=0>[N, n_scales] oipc_se_values;
  
  // Non-collinear predictors based on analysis
  matrix[N, n_scales] elevation_values;
  matrix[N, n_scales] c4_values_filled;
  
  // Climate matrices - conditional on flags
  matrix[N, include_precip ? n_scales : 0] annual_precip;
  matrix[N, include_temp ? n_scales : 0] max_temp;
  matrix[N, include_vpd ? n_scales : 0] vpd;
  matrix[N, include_soil ? n_scales : 0] soil_moisture;
  
  // PFT matrices
  matrix[N, n_scales] pft_tree;
  matrix[N, n_scales] pft_shrub;
  matrix[N, n_scales] pft_grass;
  
  // Pre-computed interaction matrices
  matrix[N, n_scales] oipc_x_c4;
  matrix[N, n_scales] oipc_x_tree;
  matrix[N, n_scales] oipc_x_shrub;
  matrix[N, n_scales] oipc_x_grass;
  
  // Climate interaction matrices (not used in current model)
  matrix[N, n_scales] temp_x_c4;
  matrix[N, n_scales] vpd_x_c4;
  matrix[N, n_scales] precip_x_tree;
  matrix[N, n_scales] precip_x_grass;
  
  // B-spline configuration (MUST BE BEFORE elevation_bspline_matrix)
  int<lower=0> n_basis_knots;  // interior knots
  int<lower=1> spline_degree;
  
  // B-spline basis for elevation (uses n_basis_knots and spline_degree)
  matrix[N * n_scales, n_basis_knots + spline_degree + 1] elevation_bspline_matrix;
  
  // GP configuration (MUST BE BEFORE arrays that use n_pp_knots)
  int<lower=1> n_pp_knots;
  
  // GP knot locations and density
  array[n_pp_knots] vector[2] knot_coords;
  vector<lower=0>[n_pp_knots] knot_data_density;
  real<lower=0> density_scaling;
  
  // Lambda control (spatial scale weighting)
  int<lower=0, upper=1> estimate_lambda;
  real<lower=0> lambda_fixed;
  real lambda_prior_mean_log;
  real<lower=0> lambda_prior_sd_log;
  
  // PC prior parameters for spatial effects
  int<lower=0, upper=1> use_pc_prior_intercept;
  real<lower=0> pc_prior_intercept_u;
  real<lower=0, upper=1> pc_prior_intercept_alpha;
  
  int<lower=0, upper=1> use_pc_prior_slope;
  real<lower=0> pc_prior_slope_u;
  real<lower=0, upper=1> pc_prior_slope_alpha;
  
  // Scaling parameters for back-transformation
  real d2H_wax_mean_original;
  real<lower=0> d2H_wax_sd_original;
}

transformed data {
  // Number of B-spline coefficients
  int n_bspline_coef = n_basis_knots + spline_degree + 1;
  
  // Average coordinate scaling for length scale conversion (MOVED UP)
  real coord_scale_km = mean(coord_scaling) * 111.0;  // Approx km per degree
  
  // Pre-compute density-based regularization factors
  vector[n_pp_knots] tau_spatial_slope;
  vector[n_pp_knots] tau_spatial_intercept;
  
  // Compute local OIPC range at each knot
  vector[n_pp_knots] oipc_range_at_knots;
  
  if (include_gp == 1) {
    // First, compute OIPC ranges WITHOUT using oipc_weighted
    // We'll use the raw oipc_values at scale 0 as a proxy
    real search_radius_km = 1000.0;  
    real search_radius = search_radius_km / coord_scale_km;
    
    for (k in 1:n_pp_knots) {
      real local_oipc_min = positive_infinity();
      real local_oipc_max = negative_infinity();
      int local_count = 0;
      
      for (n in 1:N) {
        real dist = distance(coords[n], knot_coords[k]);
        if (dist < search_radius) {
          // Use first scale as proxy for OIPC values
          local_oipc_min = fmin(local_oipc_min, oipc_values[n, 1]);
          local_oipc_max = fmax(local_oipc_max, oipc_values[n, 1]);
          local_count += 1;
        }
      }
      
      // Set range to 0 if too few points
      oipc_range_at_knots[k] = local_count > 5 ? 
                                local_oipc_max - local_oipc_min : 0.0;
    }
    
    // Compute global OIPC range for data-relative thresholds
    real max_oipc_range = max(oipc_range_at_knots);
    
    // Now compute tau values with OIPC range adjustment
    for (k in 1:n_pp_knots) {
      real d = knot_data_density[k];
      // First: existing density-based tau
      if (d == 0) {
        tau_spatial_slope[k] = 0.50; 
        tau_spatial_intercept[k] = 0.50; 
      } else if (d < 10) {
        tau_spatial_slope[k] = 0.50 + 0.30 * d/10;
        tau_spatial_intercept[k] = 0.50 + 0.30 * d/10;
      } else {
        tau_spatial_slope[k] = 0.80; 
        tau_spatial_intercept[k] = 0.80;   
      }
      
      // Multiply slope tau by OIPC range factor
      // Thresholds are relative to the global OIPC range at knots
      // (fixed thresholds in standardized units were unreachable;
      //  see validation_log/01_correlation_and_regularization_diagnostic.md)
      real range_factor;
      if (oipc_range_at_knots[k] < 0.25 * max_oipc_range) {  // < 25% of global range
        range_factor = 0.2;
      } else if (oipc_range_at_knots[k] < 0.60 * max_oipc_range) {  // 25-60% of global range
        range_factor = 0.5;
      } else {                                                        // > 60% of global range
        range_factor = 1.0;
      }
      
      tau_spatial_slope[k] = tau_spatial_slope[k] * range_factor;  // Slope only
    }
  } else {
    oipc_range_at_knots = rep_vector(0.0, n_pp_knots);
    tau_spatial_slope = rep_vector(1.0, n_pp_knots);
    tau_spatial_intercept = rep_vector(1.0, n_pp_knots);
  }
  
  // PC prior rate parameter
  real pc_prior_lambda_intercept = -log(pc_prior_intercept_alpha) / pc_prior_intercept_u;
  real pc_prior_lambda_slope = -log(pc_prior_slope_alpha) / pc_prior_slope_u;
}

parameters {
  // Global intercept and slope
  real beta_0;
  real beta_oipc;  // Global OIPC effect
  
  // B-spline coefficients for elevation
  vector[include_elevation ? n_bspline_coef : 0] beta_elev_bspline;
  array[include_elevation ? 1 : 0] real<lower=0> tau_elev_bspline_raw;
  
  // Non-collinear fixed effects (based on analysis)
  array[include_c4 ? 1 : 0] real beta_c4_raw;
  array[include_precip ? 1 : 0] real beta_precip_raw;
  
  // PFT main effects
  array[include_pft ? 1 : 0] real beta_tree_raw;
  array[include_pft ? 1 : 0] real beta_shrub_raw;
  array[include_pft ? 1 : 0] real beta_grass_raw;
  
  // Vegetation interaction effects (conditional)
  array[include_veg_interactions && include_c4 ? 1 : 0] real beta_oipc_x_c4_raw;
  array[include_veg_interactions && include_pft ? 1 : 0] real beta_oipc_x_tree_raw;
  array[include_veg_interactions && include_pft ? 1 : 0] real beta_oipc_x_shrub_raw;
  array[include_veg_interactions && include_pft ? 1 : 0] real beta_oipc_x_grass_raw;
  
  // Exponential decay for scale weighting (in km)
  array[estimate_lambda ? 1 : 0] real<lower=1, upper=400> lambda_decay_raw;
  
  // GP parameters for spatially-varying intercept and slope
  array[include_gp ? 1 : 0] real<lower=-2, upper=0> log_ls_spatial_raw;  // Single length scale for slope and intercept
  array[include_gp ? 1 : 0] real<lower=0> sigma_intercept_raw;   // SD of intercept GP
  array[include_gp ? 1 : 0] real<lower=0> sigma_slope_raw;       // SD of slope GP
  
  // Spatial random effects at knots
  vector[include_gp ? n_pp_knots : 0] z_intercept_spatial;
  vector[include_gp ? n_pp_knots : 0] z_slope_spatial;
  
  // Residual variance (no separate nugget)
  real<lower=0> sigma;
}

transformed parameters {
  // Get lambda_decay value (now directly in km)
  real lambda_decay = estimate_lambda == 1 ? lambda_decay_raw[1] : lambda_fixed;
  
  // Compute scale weights - no normalization needed
  vector[n_scales] scale_weights;
  for (i in 1:n_scales) {
    scale_weights[i] = exp(-distance_scales[i] / lambda_decay);
  }
  scale_weights = scale_weights / sum(scale_weights);
  
  // Compute weighted predictors
  vector[N] oipc_weighted = oipc_values * scale_weights;
  vector[N] oipc_se_weighted = oipc_se_values * scale_weights;
  vector[N] elev_weighted = elevation_values * scale_weights;
  vector[N] c4_weighted = include_c4 ? c4_values_filled * scale_weights : rep_vector(0, N);
  vector[N] precip_weighted = include_precip ? annual_precip * scale_weights : rep_vector(0, N);
  
  // PFT weighted predictors
  vector[N] tree_weighted = include_pft ? pft_tree * scale_weights : rep_vector(0, N);
  vector[N] shrub_weighted = include_pft ? pft_shrub * scale_weights : rep_vector(0, N);
  vector[N] grass_weighted = include_pft ? pft_grass * scale_weights : rep_vector(0, N);
  
  // Interaction weighted predictors
  vector[N] oipc_x_c4_weighted = include_veg_interactions && include_c4 ? 
                                 oipc_x_c4 * scale_weights : rep_vector(0, N);
  vector[N] oipc_x_tree_weighted = include_veg_interactions && include_pft ? 
                                   oipc_x_tree * scale_weights : rep_vector(0, N);
  vector[N] oipc_x_shrub_weighted = include_veg_interactions && include_pft ? 
                                    oipc_x_shrub * scale_weights : rep_vector(0, N);
  vector[N] oipc_x_grass_weighted = include_veg_interactions && include_pft ? 
                                    oipc_x_grass * scale_weights : rep_vector(0, N);
  
  // Initialize parameters - direct conditional access
  real beta_c4 = include_c4 == 1 ? beta_c4_raw[1] : 0;
  real beta_precip = include_precip == 1 ? beta_precip_raw[1] : 0;
  real tau_elev_bspline = include_elevation == 1 ? tau_elev_bspline_raw[1] : 1.0;
  
  // PFT parameters
  real beta_tree = include_pft == 1 ? beta_tree_raw[1] : 0;
  real beta_shrub = include_pft == 1 ? beta_shrub_raw[1] : 0;
  real beta_grass = include_pft == 1 ? beta_grass_raw[1] : 0;
  
  // Interaction parameters
  real beta_oipc_x_c4 = (include_veg_interactions == 1 && include_c4 == 1) ? 
                        beta_oipc_x_c4_raw[1] : 0;
  real beta_oipc_x_tree = (include_veg_interactions == 1 && include_pft == 1) ? 
                          beta_oipc_x_tree_raw[1] : 0;
  real beta_oipc_x_shrub = (include_veg_interactions == 1 && include_pft == 1) ? 
                           beta_oipc_x_shrub_raw[1] : 0;
  real beta_oipc_x_grass = (include_veg_interactions == 1 && include_pft == 1) ? 
                           beta_oipc_x_grass_raw[1] : 0;
  
  // Spatially-varying intercept and slope
  vector[N] alpha_spatial = rep_vector(beta_0, N);  // Intercept at each location
  vector[N] beta_oipc_spatial = rep_vector(beta_oipc, N);  // Slope at each location
  
  if (include_gp == 1) {
    real ls_spatial = exp(log_ls_spatial_raw[1]);
    real sigma_intercept = sigma_intercept_raw[1];
    real sigma_slope = sigma_slope_raw[1];
    
    // Compute GP kernel matrices ONCE
    matrix[n_pp_knots, n_pp_knots] K_knots = cov_matern32(knot_coords, 1.0, ls_spatial);
    matrix[n_pp_knots, n_pp_knots] K_knots_jitter = K_knots + 
                                                     diag_matrix(rep_vector(1e-4, n_pp_knots));
    matrix[N, n_pp_knots] K_cross = cov_matern32_cross(coords, knot_coords, 1.0, ls_spatial);
    
    // Spatial effects with density-based regularization
    vector[n_pp_knots] knot_intercepts = sigma_intercept * z_intercept_spatial;
    vector[n_pp_knots] knot_slopes = sigma_slope * z_slope_spatial;
    
    // Project to observation locations
    alpha_spatial += K_cross * mdivide_left_spd(K_knots_jitter, knot_intercepts);
    beta_oipc_spatial += K_cross * mdivide_left_spd(K_knots_jitter, knot_slopes);
}
  
  // Build location-specific linear predictor
  vector[N] mu;
  for (n in 1:N) {
    mu[n] = alpha_spatial[n] + beta_oipc_spatial[n] * oipc_weighted[n];
    
    // Add elevation effect using B-splines
    if (include_elevation == 1 && n_bspline_coef > 0) {
      // Extract the appropriate B-spline basis for this observation
      // The basis is stacked, so we need to average across scales
      real elev_effect = 0;
      for (s in 1:n_scales) {
        int row_idx = (s - 1) * N + n;
        elev_effect += scale_weights[s] * dot_product(
          elevation_bspline_matrix[row_idx, :], 
          beta_elev_bspline
        );
      }
      mu[n] += elev_effect;
    }
    
    // Add non-collinear effects
    if (include_c4 == 1) {
      mu[n] += beta_c4 * c4_weighted[n];
    }
    
    if (include_precip == 1) {
      mu[n] += beta_precip * precip_weighted[n];
    }
    
    // Add PFT main effects
    if (include_pft == 1) {
      mu[n] += beta_tree * tree_weighted[n] + 
               beta_shrub * shrub_weighted[n] + 
               beta_grass * grass_weighted[n];
    }
    
    // Add vegetation interaction effects
    if (include_veg_interactions == 1) {
      if (include_c4 == 1) {
        mu[n] += beta_oipc_x_c4 * oipc_x_c4_weighted[n];
      }
      if (include_pft == 1) {
        mu[n] += beta_oipc_x_tree * oipc_x_tree_weighted[n] +
                 beta_oipc_x_shrub * oipc_x_shrub_weighted[n] +
                 beta_oipc_x_grass * oipc_x_grass_weighted[n];
      }
    }
  }
}

model {
  // Priors
  beta_0 ~ normal(0, 5);
  beta_oipc ~ normal(0.8, 0.3);  // Expect slope around 0.8
  
  // B-spline coefficients with smoothing
  if (include_elevation == 1) {
    tau_elev_bspline_raw[1] ~ normal(0, 1);
    beta_elev_bspline[1] ~ normal(0, 2);  // First coefficient
    for (k in 2:n_bspline_coef) {
      beta_elev_bspline[k] ~ normal(beta_elev_bspline[k-1], tau_elev_bspline);
    }
  }
  
  // Non-collinear effects
  if (include_c4 == 1) {
    beta_c4_raw[1] ~ normal(0, 2);
  }
  
  if (include_precip == 1) {
    // Weakly informative, matching other effect prior scales
    beta_precip_raw[1] ~ normal(0, 0.5);
  }
  
  // PFT main effects
  if (include_pft == 1) {
    beta_tree_raw[1] ~ normal(0, 2);
    beta_shrub_raw[1] ~ normal(0, 2);
    beta_grass_raw[1] ~ normal(0, 2);
  }
  
  // Vegetation interaction effects
  if (include_veg_interactions == 1) {
    if (include_c4 == 1) {
      beta_oipc_x_c4_raw[1] ~ normal(0, 0.5);  // Moderate interaction expected
    }
    if (include_pft == 1) {
      beta_oipc_x_tree_raw[1] ~ normal(0, 0.5);
      beta_oipc_x_shrub_raw[1] ~ normal(0, 0.5);
      beta_oipc_x_grass_raw[1] ~ normal(0, 0.5);
    }
  }
  
  // Scale weighting - updated prior for direct km interpretation
  if (estimate_lambda == 1) {
    lambda_decay_raw[1] ~ lognormal(lambda_prior_mean_log, lambda_prior_sd_log);
  }
  
  // GP hyperparameters
 if (include_gp == 1) {
    // Length scale prior
    log_ls_spatial_raw[1] ~ normal(-1.0, 0.4);
    
    // Spatial SD priors
    sigma_intercept_raw[1] ~ exponential(pc_prior_lambda_intercept);
    sigma_slope_raw[1] ~ exponential(pc_prior_lambda_slope);
    
    // Spatially-varying effects with density regularization
    for (k in 1:n_pp_knots) {
      z_intercept_spatial[k] ~ normal(0, tau_spatial_intercept[k]);
      z_slope_spatial[k] ~ normal(0, tau_spatial_slope[k]);
    }
}  

// Residual variance (combines nugget and residual)
sigma ~ normal(0, 2);

// Likelihood: OIPC measurement error propagated as beta^2 * SE_x^2
for (n in 1:N) {
    real total_var = square(d2H_wax_err[n]) + square(beta_oipc_spatial[n] * oipc_se_weighted[n]) + square(sigma);
    real total_sd = sqrt(total_var);
    d2H_wax[n] ~ normal(mu[n], total_sd);
}
}

generated quantities {
  // Posterior predictive checks
  vector[N] d2H_rep;
  vector[N] log_lik;
  
  for (n in 1:N) {
    // Same error propagation as likelihood block above
    real total_var = square(d2H_wax_err[n]) + square(beta_oipc_spatial[n] * oipc_se_weighted[n]) + square(sigma);
    real total_sd = sqrt(total_var);
    d2H_rep[n] = normal_rng(mu[n], total_sd);
    log_lik[n] = normal_lpdf(d2H_wax[n] | mu[n], total_sd);
  }
  
  // Back-transform parameters
  real intercept_original = beta_0 * d2H_wax_sd_original + d2H_wax_mean_original;
  real sigma_residual_original = sigma * d2H_wax_sd_original;
  
  // Summary of spatial variation
  real sigma_intercept_spatial = 0;
  real sigma_slope_spatial = 0;
  real ls_intercept_km = 0;
  real ls_slope_km = 0;
  
  if (include_gp == 1) {
    sigma_intercept_spatial = sigma_intercept_raw[1] * d2H_wax_sd_original;
    sigma_slope_spatial = sigma_slope_raw[1];
    // Proper conversion using coordinate scaling
    ls_intercept_km = exp(log_ls_spatial_raw[1]) * coord_scale_km;
    ls_slope_km = exp(log_ls_spatial_raw[1]) * coord_scale_km;
  }
  
  // Effective scale (lambda directly in km now)
  real effective_scale_km = lambda_decay;
  
  // Range of OIPC slopes across space
  real min_oipc_slope = min(beta_oipc_spatial);
  real max_oipc_slope = max(beta_oipc_spatial);
  real sd_oipc_slope = sd(beta_oipc_spatial);
  
  // Variance decomposition (simplified without nugget)
  real var_spatial_intercept = include_gp ? variance(alpha_spatial - beta_0) : 0;
  real var_spatial_slope = include_gp ? variance((beta_oipc_spatial - beta_oipc) .* oipc_weighted) : 0;
  real var_total_spatial = var_spatial_intercept + var_spatial_slope;
  real var_residual = square(sigma);
  real var_total = var_total_spatial + var_residual;
  
  real prop_variance_spatial = var_total > 0 ? var_total_spatial / var_total : 0;
  real prop_variance_residual = var_total > 0 ? var_residual / var_total : 0;
  
  // Model fit metrics
  real mse = mean(square(d2H_wax - mu));
  real rmse = sqrt(mse);
  real r_squared = 1 - variance(d2H_wax - mu) / variance(d2H_wax);
  
  // Diagnostic outputs for OIPC range analysis
  vector[n_pp_knots] knot_oipc_ranges = oipc_range_at_knots;
  vector[n_pp_knots] tau_final_slope = tau_spatial_slope;
  vector[n_pp_knots] tau_final_intercept = tau_spatial_intercept;
  
  // Range factor threshold diagnostics
  real range_threshold_low = include_gp ? 0.25 * max(oipc_range_at_knots) : 0;
  real range_threshold_high = include_gp ? 0.60 * max(oipc_range_at_knots) : 0;

}
