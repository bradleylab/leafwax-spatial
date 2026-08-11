#!/usr/bin/env Rscript
# 2g_data_audit.R
#
# Data audit (detection only) for the leaf-wax n-C29 dD calibration compilation.
# Non-destructive: reads input_data/global_data_c29.csv and writes review tables
# to data/audit/ -- merges/removes nothing. The frozen dataset is built downstream
# by 2i_freeze_calibration.R once the duplicate + archive decisions are settled.
#
# Methods:
#   1a  archive split for "Sediment": coordinate land/ocean heuristic
#       (in ocean -> marine sediment; on land -> lake sediment; near coast -> flag).
#   2   conservative, DOI-aware duplicate detection: within-source multiples are
#       genuine replicates (kept); only CROSS-source SAME-DOI coordinate matches
#       whose d2H_wax agree within tolerance are proposed as re-ingestion
#       duplicates; everything else -> review. Nothing auto-dropped.
#
# Run: Rscript 2g_data_audit.R   (from repo root)
# Out: data/audit/{duplicate_candidates,duplicate_decisions,archive_type_proposal,
#                  exclusion_log,provenance_summary}.csv  (+ console report)

suppressWarnings(suppressMessages({
  library(dplyr); library(sf); library(rnaturalearth)
}))

WAX_TOL <- 5          # permil; cross-source rows within this are "same sample" candidates
COORD_DP <- 2         # ~1.1 km grouping resolution
COAST_BUFFER_KM <- 25 # sediment sites within this of the coastline -> ambiguous

out_dir <- "data/audit"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
raw <- read.csv("input_data/global_data_c29.csv", stringsAsFactors = FALSE)
raw$obs_id <- sprintf("obs_%04d", seq_len(nrow(raw)))
raw$compilation_clean <- ifelse(is.na(raw$compilation) | raw$compilation == "",
                                "direct primary ingest", raw$compilation)

# ---------------------------------------------------------------------------
# 1a. Archive-type heuristic (lake vs marine for "Sediment")
# ---------------------------------------------------------------------------
land  <- ne_countries(scale = "medium", returnclass = "sf") |> st_make_valid() |> st_union()
coast <- ne_coastline(scale = "medium", returnclass = "sf")
pts <- st_as_sf(raw, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

on_land <- lengths(st_intersects(pts, land)) > 0
dist_coast_km <- as.numeric(st_distance(pts, st_union(coast))) / 1000

raw$archive_type <- dplyr::case_when(
  raw$sample_type == "Soil" ~ "soil",
  raw$sample_type == "Sediment" &  on_land ~ "lake sediment",
  raw$sample_type == "Sediment" & !on_land ~ "marine sediment",
  TRUE ~ "other/unknown"
)
raw$archive_ambiguous <- raw$sample_type == "Sediment" & dist_coast_km < COAST_BUFFER_KM
raw$dist_coast_km <- round(dist_coast_km, 1)

archive_prop <- raw |>
  select(obs_id, source, compilation_clean, latitude, longitude,
         sample_type, archive_type, archive_ambiguous, dist_coast_km) |>
  arrange(desc(archive_ambiguous), archive_type)
write.csv(archive_prop, file.path(out_dir, "archive_type_proposal.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Duplicate detection (conservative)
# ---------------------------------------------------------------------------
raw$coord_key <- paste(round(raw$latitude, COORD_DP), round(raw$longitude, COORD_DP), raw$chain)

grp <- raw |>
  group_by(coord_key) |>
  summarise(n = dplyr::n(), n_sources = n_distinct(source), n_dois = n_distinct(DOI),
            wax_spread = ifelse(n > 1, max(d2H_wax) - min(d2H_wax), 0),
            obs_ids = paste(obs_id, collapse = ";"), .groups = "drop")

multi <- grp |> filter(n > 1)
# Conservative + DOI-aware: only SAME-DOI cross-source matches with agreeing
# d2H_wax are treated as re-ingestion duplicates. Different-DOI coincidences at
# the same location are genuine-different-study candidates -> review (keep both).
multi$class <- dplyr::case_when(
  multi$n_sources == 1                                                       ~ "within-source replicates (keep all)",
  multi$n_sources > 1 & multi$n_dois == 1 & multi$wax_spread <= WAX_TOL      ~ "CROSS-SOURCE same-DOI duplicate",
  multi$n_sources > 1 & multi$n_dois == 1 & multi$wax_spread >  WAX_TOL      ~ "cross-source same-DOI, values differ (review)",
  TRUE                                                                       ~ "cross-source different-DOI (review)"
)
write.csv(multi[order(-multi$n_sources, -multi$n), ],
          file.path(out_dir, "duplicate_candidates.csv"), row.names = FALSE)

# Proposed decisions, one row per observation in a multi-row group.
dec <- raw |>
  filter(coord_key %in% multi$coord_key) |>
  left_join(multi |> select(coord_key, class, n, n_sources, wax_spread), by = "coord_key") |>
  group_by(coord_key) |>
  mutate(
    canonical = obs_id == obs_id[which.min(ifelse(is.na(d2H_wax_err), 999, d2H_wax_err))],
    proposed_action = dplyr::case_when(
      grepl("within-source", class)              ~ "keep (genuine replicate)",
      grepl("same-DOI duplicate", class) &  canonical ~ "keep (canonical of duplicate set)",
      grepl("same-DOI duplicate", class) & !canonical ~ "DROP (same-DOI re-ingestion)",
      TRUE                                       ~ "REVIEW (decide manually)"
    )
  ) |>
  ungroup() |>
  select(coord_key, obs_id, source, compilation_clean, DOI, latitude, longitude,
         chain, d2H_wax, d2H_wax_err, archive_type, class, n, n_sources,
         wax_spread, canonical, proposed_action) |>
  arrange(coord_key, desc(canonical))
write.csv(dec, file.path(out_dir, "duplicate_decisions.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Provenance summary + exclusion-criteria log (no exclusions applied yet)
# ---------------------------------------------------------------------------
prov <- raw |>
  group_by(compilation_clean) |>
  summarise(n_obs = dplyr::n(), n_dois = n_distinct(DOI),
            n_soil = sum(archive_type == "soil"),
            n_lake = sum(archive_type == "lake sediment"),
            n_marine = sum(archive_type == "marine sediment"), .groups = "drop") |>
  arrange(desc(n_obs))
write.csv(prov, file.path(out_dir, "provenance_summary.csv"), row.names = FALSE)

excl <- data.frame(
  criterion = c("compound", "coordinates", "d2H_wax", "blank compilation",
                "cross-source duplicate"),
  rule = c("keep n-C29 only (already enforced upstream)",
           "require finite latitude & longitude",
           "require finite d2H_wax",
           "relabel compilation as 'direct primary ingest' (NOT excluded)",
           "drop non-canonical members of confirmed SAME-DOI re-ingestion sets only"),
  rows_affected = c(NA,
                    sum(!is.finite(raw$latitude) | !is.finite(raw$longitude)),
                    sum(!is.finite(raw$d2H_wax)),
                    sum(raw$compilation_clean == "direct primary ingest"),
                    sum(dec$proposed_action == "DROP (same-DOI re-ingestion)")),
  status = "REQUIRES MANUAL REVIEW"
)
write.csv(excl, file.path(out_dir, "exclusion_log.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Console report
# ---------------------------------------------------------------------------
cat("=== CEE data audit (DETECTION ONLY — nothing merged) ===\n")
cat("raw rows:", nrow(raw), "\n\n")
cat("--- archive-type proposal (1a heuristic) ---\n")
print(table(raw$archive_type))
cat("Sediment rows flagged ambiguous (within", COAST_BUFFER_KM, "km of coast):",
    sum(raw$archive_ambiguous), "\n\n")
cat("--- duplicate candidates (tol =", WAX_TOL, "permil, coord", COORD_DP, "dp) ---\n")
print(table(multi$class))
cat("\nProposed actions across grouped observations:\n")
print(table(dec$proposed_action))
cat("\nNET if all proposals accepted: drop",
    sum(dec$proposed_action == "DROP (same-DOI re-ingestion)"),
    "of", nrow(raw), "rows ->",
    nrow(raw) - sum(dec$proposed_action == "DROP (same-DOI re-ingestion)"), "canonical.\n")
cat("\nReview files in data/audit/:\n  duplicate_candidates.csv\n  duplicate_decisions.csv",
    "\n  archive_type_proposal.csv\n  provenance_summary.csv\n  exclusion_log.csv\n")
cat("\nNOTHING removed. Frozen file built only after sign-off.\n")
