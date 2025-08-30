#───────────────────────────────────────────────────────────────────────────────
# Generate comprehensive summary table for revised clusters
# Run after the revised clustering analysis
#───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(openxlsx)
library(broom)

# Load the revised clustered data
sediment_with_clusters <- readRDS("results/site_clustering_revised/sediment_with_14clusters_revised.rds")

# Function to format confidence intervals
format_ci <- function(estimate, lower, upper, digits = 2) {
  sprintf("%.2f [%.2f, %.2f]", 
          round(estimate, digits), 
          round(lower, digits), 
          round(upper, digits))
}

# Function to describe cluster based on characteristics
describe_cluster_detailed <- function(row) {
  lat <- as.numeric(row["mean_latitude"])
  elev <- as.numeric(row["mean_elevation"])
  temp <- as.numeric(row["mean_max_temp"])
  precip <- as.numeric(row["mean_annual_precip"])
  c4 <- as.numeric(row["mean_C4_fraction"])
  
  # Hemisphere
  if (lat < -20) {
    hemisphere <- "Southern"
  } else if (lat > 20) {
    hemisphere <- "Northern"
  } else {
    hemisphere <- "Tropical"
  }
  
  # Latitude zone
  if (abs(lat) < 10) {
    lat_zone <- "equatorial"
  } else if (abs(lat) < 23.5) {
    lat_zone <- "tropical"
  } else if (abs(lat) < 35) {
    lat_zone <- "subtropical"
  } else if (abs(lat) < 50) {
    lat_zone <- "temperate"
  } else if (abs(lat) < 60) {
    lat_zone <- "boreal"
  } else {
    lat_zone <- "arctic"
  }
  
  # Elevation
  if (elev < 200) {
    elev_desc <- "lowland"
  } else if (elev < 500) {
    elev_desc <- "low elevation"
  } else if (elev < 1500) {
    elev_desc <- "mid-elevation"
  } else if (elev < 2500) {
    elev_desc <- "highland"
  } else {
    elev_desc <- "high mountain"
  }
  
  # Temperature
  if (temp < 0) {
    temp_desc <- "very cold"
  } else if (temp < 10) {
    temp_desc <- "cold"
  } else if (temp < 15) {
    temp_desc <- "cool"
  } else if (temp < 20) {
    temp_desc <- "mild"
  } else if (temp < 25) {
    temp_desc <- "warm"
  } else if (temp < 30) {
    temp_desc <- "hot"
  } else {
    temp_desc <- "very hot"
  }
  
  # Precipitation
  if (precip < 200) {
    precip_desc <- "arid"
  } else if (precip < 400) {
    precip_desc <- "semi-arid"
  } else if (precip < 800) {
    precip_desc <- "moderate rainfall"
  } else if (precip < 1500) {
    precip_desc <- "wet"
  } else if (precip < 2000) {
    precip_desc <- "very wet"
  } else {
    precip_desc <- "extremely wet"
  }
  
  # C4 vegetation
  if (c4 < 5) {
    veg_desc <- "C3-dominated"
  } else if (c4 < 20) {
    veg_desc <- "mixed C3/C4"
  } else if (c4 < 40) {
    veg_desc <- "significant C4"
  } else {
    veg_desc <- "C4-dominated"
  }
  
  # Combine descriptions
  full_desc <- paste0(hemisphere, " ", lat_zone, " ", elev_desc, " - ", 
                      temp_desc, ", ", precip_desc, " (", veg_desc, ")")
  
  return(full_desc)
}

#───────────────────────────────────────────────────────────────────────────────
# CALCULATE COMPREHENSIVE CLUSTER SUMMARIES
#───────────────────────────────────────────────────────────────────────────────

cat("Calculating revised cluster summaries...\n")

cluster_summary <- data.frame()

for (k in sort(unique(sediment_with_clusters$site_cluster[!is.na(sediment_with_clusters$site_cluster)]))) {
  
  # Get data for this cluster
  cluster_data <- sediment_with_clusters %>%
    filter(site_cluster == k)
  
  # Fit linear model for OIPC relationship
  if (nrow(cluster_data) > 3) {
    lm_fit <- lm(d2H_wax ~ oipc_d2h20, data = cluster_data)
    lm_summary <- summary(lm_fit)
    lm_confint <- confint(lm_fit)
    
    # Extract model statistics
    slope <- coef(lm_fit)[2]
    intercept <- coef(lm_fit)[1]
    slope_se <- lm_summary$coefficients[2, 2]
    intercept_se <- lm_summary$coefficients[1, 2]
    slope_p <- lm_summary$coefficients[2, 4]
    r_squared <- lm_summary$r.squared
    adj_r_squared <- lm_summary$adj.r.squared
    residual_se <- lm_summary$sigma
    
    # F-statistic and model p-value
    if (!is.null(lm_summary$fstatistic)) {
      f_stat <- lm_summary$fstatistic[1]
      model_p <- pf(lm_summary$fstatistic[1], 
                    lm_summary$fstatistic[2], 
                    lm_summary$fstatistic[3], 
                    lower.tail = FALSE)
    } else {
      f_stat <- NA
      model_p <- NA
    }
  } else {
    # Too few points for regression
    slope <- NA
    intercept <- NA
    slope_p <- NA
    r_squared <- NA
    adj_r_squared <- NA
    residual_se <- NA
    f_stat <- NA
    model_p <- NA
    lm_confint <- matrix(NA, 2, 2)
  }
  
  # Create summary row
  summary_row <- data.frame(
    # Cluster identification
    cluster = k,
    n_sites = nrow(cluster_data),
    
    # Geographic summary
    mean_latitude = mean(cluster_data$latitude),
    sd_latitude = sd(cluster_data$latitude),
    mean_longitude = mean(cluster_data$longitude),
    sd_longitude = sd(cluster_data$longitude),
    
    # Environmental means
    mean_elevation = mean(cluster_data$elevation_gmted),
    mean_annual_precip = mean(cluster_data$annual_precip),
    mean_max_temp = mean(cluster_data$max_temp),
    mean_vpd = mean(cluster_data$vpd),
    mean_C4_fraction = mean(cluster_data$C4_fraction_5deg),
    
    # Environmental SDs
    sd_elevation = sd(cluster_data$elevation_gmted),
    sd_annual_precip = sd(cluster_data$annual_precip),
    sd_max_temp = sd(cluster_data$max_temp),
    sd_vpd = sd(cluster_data$vpd),
    sd_C4_fraction = sd(cluster_data$C4_fraction_5deg),
    
    # Isotope values
    mean_oipc_d2h20 = mean(cluster_data$oipc_d2h20),
    sd_oipc_d2h20 = sd(cluster_data$oipc_d2h20),
    range_oipc = max(cluster_data$oipc_d2h20) - min(cluster_data$oipc_d2h20),
    mean_d2H_wax = mean(cluster_data$d2H_wax),
    sd_d2H_wax = sd(cluster_data$d2H_wax),
    
    # Regression statistics
    slope = slope,
    slope_CI = ifelse(is.na(slope), "NA", 
                      format_ci(slope, lm_confint[2,1], lm_confint[2,2])),
    slope_p_value = slope_p,
    intercept = intercept,
    intercept_CI = ifelse(is.na(intercept), "NA",
                          format_ci(intercept, lm_confint[1,1], lm_confint[1,2])),
    r_squared = r_squared,
    adj_r_squared = adj_r_squared,
    residual_se = residual_se,
    f_statistic = f_stat,
    model_p_value = model_p
  )
  
  # Add PFT data if available
  if ("pft_tree" %in% names(cluster_data)) {
    summary_row$mean_pft_tree <- mean(sapply(cluster_data$pft_tree, mean, na.rm = TRUE))
    summary_row$mean_pft_shrub <- mean(sapply(cluster_data$pft_shrub, mean, na.rm = TRUE))
    summary_row$mean_pft_grass <- mean(sapply(cluster_data$pft_grass, mean, na.rm = TRUE))
  }
  
  # Bind to results
  cluster_summary <- rbind(cluster_summary, summary_row)
}

# Add descriptions
cluster_summary$description <- apply(cluster_summary, 1, describe_cluster_detailed)

# Add significance flag
cluster_summary$significant <- ifelse(cluster_summary$slope_p_value < 0.05, "Yes", "No")

# Reorder for readability
cluster_summary <- cluster_summary %>%
  select(
    # Identity
    cluster, description, n_sites,
    
    # Geography
    mean_latitude, sd_latitude, mean_longitude, sd_longitude,
    
    # Environment
    mean_elevation, mean_annual_precip, mean_max_temp, mean_vpd, 
    mean_C4_fraction,
    
    # Isotopes
    mean_oipc_d2h20, sd_oipc_d2h20, range_oipc,
    mean_d2H_wax, sd_d2H_wax,
    
    # Regression results
    slope_CI, intercept_CI, r_squared, adj_r_squared,
    slope_p_value, model_p_value, significant, residual_se,
    
    # Everything else
    everything()
  )

#───────────────────────────────────────────────────────────────────────────────
# CREATE FORMATTED EXCEL OUTPUT
#───────────────────────────────────────────────────────────────────────────────

cat("Creating Excel file...\n")

# Create workbook
wb <- createWorkbook()

# Add main summary sheet
addWorksheet(wb, "Revised Cluster Summary")
writeData(wb, "Revised Cluster Summary", cluster_summary)

# Format headers
headerStyle <- createStyle(
  fontSize = 12,
  fontColour = "#FFFFFF",
  fgFill = "#4472C4",
  halign = "center",
  valign = "center",
  textDecoration = "bold",
  border = "TopBottomLeftRight"
)

addStyle(wb, "Revised Cluster Summary", headerStyle, 
         rows = 1, cols = 1:ncol(cluster_summary), 
         gridExpand = TRUE)

# Format numbers
numberStyle <- createStyle(numFmt = "0.00")
addStyle(wb, "Revised Cluster Summary", numberStyle,
         rows = 2:(nrow(cluster_summary)+1),
         cols = which(sapply(cluster_summary, is.numeric)),
         gridExpand = TRUE)

# Highlight significant relationships
significantStyle <- createStyle(fgFill = "#E7F5E7")
sig_rows <- which(cluster_summary$significant == "Yes") + 1
if (length(sig_rows) > 0) {
  addStyle(wb, "Revised Cluster Summary", significantStyle,
           rows = sig_rows,
           cols = which(names(cluster_summary) %in% c("slope_CI", "r_squared", "significant")),
           gridExpand = TRUE)
}

# Highlight strong relationships (R² > 0.3)
strongStyle <- createStyle(fgFill = "#FFE6E6")
strong_rows <- which(cluster_summary$r_squared > 0.3) + 1
if (length(strong_rows) > 0) {
  addStyle(wb, "Revised Cluster Summary", strongStyle,
           rows = strong_rows,
           cols = which(names(cluster_summary) == "r_squared"),
           gridExpand = TRUE)
}

# Adjust column widths
setColWidths(wb, "Revised Cluster Summary", 
             cols = 1:ncol(cluster_summary), 
             widths = "auto")

# Add interpretation sheet
interpretation_data <- data.frame(
  Cluster = 1:14,
  Short_Name = c(
    "Tropical lowland rainforest",
    "Subtropical highland",
    "Arctic/Subarctic",
    "Temperate highland steppe",
    "Southern dry subtropics",
    "East Asian monsoon",
    "Tibetan Plateau",
    "Southern tropical dry",
    "Central Asian cold steppe",
    "Equatorial Africa",
    "Arctic maritime",
    "Eastern US temperate",
    "East Asian mountains",
    "African highlands"
  ),
  Key_Features = c(
    "Hot, extremely wet, low elevation tropical sites",
    "Warm subtropical mountains with high rainfall",
    "Cold arctic sites with moderate moisture",
    "Mid-latitude dry highlands",
    "Southern hemisphere arid/semi-arid regions",
    "Monsoon-influenced lowlands",
    "High elevation cold dry plateau",
    "Dry tropical/subtropical S. hemisphere",
    "Cold continental highlands",
    "Hot equatorial sites with high C4",
    "Maritime arctic with strong OIPC relationship",
    "Temperate eastern North America",
    "High mountains of East Asia",
    "African highlands with mixed C3/C4"
  ),
  OIPC_Relationship = c(
    "Weak positive",
    "Very weak",
    "Moderate positive",
    "Strong positive",
    "None",
    "Weak negative (unreliable)",
    "Moderate positive",
    "Weak positive",
    "Very weak",
    "Strong positive",
    "Very strong positive",
    "Very weak",
    "Weak positive",
    "Moderate positive"
  )
)

addWorksheet(wb, "Cluster Interpretation")
writeData(wb, "Cluster Interpretation", interpretation_data)

# Add metadata sheet
metadata <- data.frame(
  Information = c(
    "Analysis Date",
    "Total Sites",
    "Sites with Clusters",
    "Number of Clusters",
    "Variables Used for Clustering",
    "Variables Removed",
    "",
    "Key Improvements:",
    "- No double-weighting of correlated variables",
    "- Clustering based on environment only (not d2H_wax)",
    "- More balanced importance across variable types",
    "",
    "Color Coding:",
    "- Green: Significant slopes (p<0.05)",
    "- Pink: Strong relationships (R²>0.3)"
  ),
  Value = c(
    as.character(Sys.Date()),
    as.character(nrow(sediment_with_clusters)),
    as.character(sum(!is.na(sediment_with_clusters$site_cluster))),
    "14",
    "elevation, precip, temp, VPD, C4, OIPC, PFTs",
    "abs_latitude, soil_moisture, d2H_wax",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    ""
  )
)

addWorksheet(wb, "Metadata")
writeData(wb, "Metadata", metadata)

# Save the workbook
output_file <- "results/site_clustering_revised/revised_cluster_summary_table.xlsx"
saveWorkbook(wb, output_file, overwrite = TRUE)

cat("\n✓ Revised cluster summary table saved to:", output_file, "\n")

#───────────────────────────────────────────────────────────────────────────────
# PRINT SUMMARY
#───────────────────────────────────────────────────────────────────────────────

cat("\nREVISED CLUSTER CHARACTERISTICS:\n")
cat("================================\n\n")

# Print brief description of each cluster
for (i in 1:nrow(cluster_summary)) {
  cat(sprintf("Cluster %2d (n=%3d): %s\n", 
              cluster_summary$cluster[i],
              cluster_summary$n_sites[i],
              cluster_summary$description[i]))
  
  if (cluster_summary$r_squared[i] > 0.3) {
    cat(sprintf("           → Strong OIPC relationship: R²=%.3f, slope=%s\n",
                cluster_summary$r_squared[i],
                cluster_summary$slope_CI[i]))
  } else if (cluster_summary$significant[i] == "Yes") {
    cat(sprintf("           → Significant but weak: R²=%.3f\n",
                cluster_summary$r_squared[i]))
  } else if (cluster_summary$r_squared[i] < 0.01) {
    cat("           → No OIPC relationship\n")
  }
}

cat("\n✓ Analysis complete!\n")