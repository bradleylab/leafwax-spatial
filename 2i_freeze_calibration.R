#!/usr/bin/env Rscript
# 2i_freeze_calibration.R
#
# Builds the audited leaf-wax n-C29 delta-2H calibration dataset from the source
# compilation and its versioned curation table. The source compilation
# (input_data/global_data_c29.csv) is read-only and untouched.
#
# Decisions applied:
#   - One value-matched re-ingestion duplicate is excluded: Garcin 2012 Barombi
#     Mbo. Other same-DOI coordinate matches are distinct samples.
#   - Archive classes combine a coordinate-based initial classification with
#     source-specific review documented in input_data/calibration_curation_v1.csv.
#   - Review flags are retained in the curation table but are not exclusions.
#
# Run:  Rscript 2i_freeze_calibration.R   (from repo root)
# Out:  input_data/leafwax_d2h_c29_calibration_v1.{csv,rds}  (+ .parquet if arrow)
#       input_data/leafwax_d2h_c29_calibration_v1_dictionary.md
#       input_data/leafwax_d2h_c29_calibration_v1_exclusions.csv

suppressWarnings(suppressMessages(library(dplyr)))

DATASET <- "leafwax_d2h_c29_calibration_v1"
EXPECTED_SOURCE_MD5 <- "5c451a2ecc4db4f51fbfd15cffe2bc99"
ALLOWED_ARCHIVE_CLASSES <- c(
  "soil", "lake sediment", "marine sediment", "fluvial sediment"
)

repo_root <- "."   # scripts run from the repository root (see README)
frozen_dir <- Sys.getenv(
  "LEAFWAX_FROZEN_DATA_DIR",
  unset = file.path(repo_root, "input_data")
)
dir.create(frozen_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Source compilation + stable row ids ------------------------------------
source_path <- file.path(repo_root, "input_data", "global_data_c29.csv")
source_md5 <- unname(tools::md5sum(source_path))
if (!identical(source_md5, EXPECTED_SOURCE_MD5)) {
  stop(
    "global_data_c29.csv does not match the curated source version. ",
    "Expected MD5 ", EXPECTED_SOURCE_MD5, "; found ", source_md5, "."
  )
}

raw <- read.csv(source_path,
                stringsAsFactors = FALSE, check.names = FALSE)
raw$obs_id <- sprintf("obs_%04d", seq_len(nrow(raw)))
n_raw <- nrow(raw)

# ---- Versioned scientific curation decisions --------------------------------
curation_path <- file.path(repo_root, "input_data", "calibration_curation_v1.csv")
curation <- read.csv(curation_path, stringsAsFactors = FALSE, check.names = FALSE)

required_curation_columns <- c(
  "obs_id", "include_in_calibration", "exclusion_reason", "archive_class",
  "archive_assignment", "review_flag"
)
missing_curation_columns <- setdiff(required_curation_columns, names(curation))
if (length(missing_curation_columns) > 0) {
  stop("Missing curation columns: ", paste(missing_curation_columns, collapse = ", "))
}
if (!identical(curation$obs_id, raw$obs_id)) {
  stop("calibration_curation_v1.csv must contain every source obs_id in source-row order.")
}
if (anyDuplicated(curation$obs_id)) {
  stop("calibration_curation_v1.csv contains duplicate obs_id values.")
}
if (anyNA(curation$include_in_calibration)) {
  stop("include_in_calibration must be TRUE or FALSE for every source row.")
}
unexpected_archive_classes <- setdiff(unique(curation$archive_class), ALLOWED_ARCHIVE_CLASSES)
if (length(unexpected_archive_classes) > 0) {
  stop("Unexpected archive classes: ", paste(unexpected_archive_classes, collapse = ", "))
}
if (any(!curation$include_in_calibration & !nzchar(curation$exclusion_reason))) {
  stop("Every excluded source row must have an exclusion_reason.")
}

# ---- Assemble ----------------------------------------------------------------
d <- raw %>%
  left_join(curation, by = "obs_id")
stopifnot(nrow(d) == n_raw)   # joins must not fan out

# ---- Apply exclusions --------------------------------------------------------
d <- d %>% mutate(
  .finite_ok = is.finite(latitude) & is.finite(longitude) & is.finite(d2H_wax),
  .is_drop   = !include_in_calibration)
n_nonfinite <- sum(!d$.finite_ok)
drop_ids <- d$obs_id[d$.is_drop]

frozen <- d %>% filter(.finite_ok & !.is_drop)
n_flagged <- sum(nzchar(frozen$review_flag))

# Final column set: 13 core fields (identity/provenance, coordinates, predictor,
# response, audited medium). The curation table retains the assignment basis,
# exclusion decisions, and review flags. The unreliable raw sample_type, the
# constant chain (==29), and the near-empty d2H_precip_year are omitted.
col_order <- c("obs_id", "source", "compilation", "location", "DOI",
               "latitude", "longitude", "elevation",
               "d2H_precip", "d2H_precip_err", "d2H_wax", "d2H_wax_err",
               "archive_class")
frozen <- frozen[, col_order]

# ---- Write data --------------------------------------------------------------
csv_path <- file.path(frozen_dir, paste0(DATASET, ".csv"))
rds_path <- file.path(frozen_dir, paste0(DATASET, ".rds"))
write.csv(frozen, csv_path, row.names = FALSE)
saveRDS(frozen, rds_path)
parquet_ok <- FALSE
try({ if (requireNamespace("arrow", quietly = TRUE)) {
  arrow::write_parquet(frozen, file.path(frozen_dir, paste0(DATASET, ".parquet")))
  parquet_ok <- TRUE } }, silent = TRUE)

# Integrity stamp: written here so it can never drift from the csv it describes
# (R 4.4 has no tools::sha256; use digest, the common dep). Format matches
# `shasum -a 256`: hash + two spaces + bare filename.
sha_ok <- FALSE
try({ if (requireNamespace("digest", quietly = TRUE)) {
  csv_sha <- digest::digest(csv_path, algo = "sha256", file = TRUE)
  writeLines(sprintf("%s  %s.csv", csv_sha, DATASET),
             file.path(frozen_dir, paste0(DATASET, ".sha256")))
  sha_ok <- TRUE } }, silent = TRUE)

# ---- Exclusion log -----------------------------------------------------------
excl_rows <- d %>% filter(.is_drop | !.finite_ok) %>%
  mutate(exclusion_reason = case_when(
    .is_drop ~ exclusion_reason,
    TRUE     ~ "non-finite latitude/longitude/d2H_wax")) %>%
  select(obs_id, source, compilation, DOI, latitude, longitude, d2H_wax, exclusion_reason)
write.csv(excl_rows, file.path(frozen_dir, paste0(DATASET, "_exclusions.csv")), row.names = FALSE)

# ---- Data dictionary ---------------------------------------------------------
dict <- c(
  "# Audited leaf-wax calibration data dictionary", "",
  "Audited global leaf-wax *n*-C29 alkane δ²H calibration compilation. Built by",
  "`2i_freeze_calibration.R` from `input_data/global_data_c29.csv`",
  "(read-only) and `input_data/calibration_curation_v1.csv`.",
  "",
  sprintf("- Source rows: %d", n_raw),
  sprintf("- Duplicates dropped (same-DOI re-ingestion): %d", length(drop_ids)),
  sprintf("- Non-finite lat/lon/d2H_wax dropped: %d", n_nonfinite),
  sprintf("- **Audited calibration rows: %d**", nrow(frozen)),
  "",
  "## Columns (13)", "",
  "- `obs_id` -- stable id, `obs_%04d` by source-compilation row order.",
  "- `source`, `compilation`, `location`, `DOI` -- provenance as compiled.",
  "- `latitude`, `longitude`, `elevation` -- site coordinates (WGS84) and metres.",
  "- `d2H_precip`, `d2H_precip_err` -- precipitation δ²H predictor (‰, VSMOW).",
  "- `d2H_wax`, `d2H_wax_err` -- leaf-wax *n*-C29 δ²H response (‰, VSMOW).",
  "- `archive_class` -- final depositional class: soil / lake sediment /",
  "  marine sediment / **fluvial sediment**. Per-row assignments and review",
  "  flags are recorded in `input_data/calibration_curation_v1.csv`.",
  "",
  "Omitted from the audited calibration file:",
  "`chain` (constant 29), `d2H_precip_year` (near-empty), `sample_type` (a coarse",
  "source field not used after archive-class curation), and the curation columns",
  "`archive_assignment`, `review_flag`, and `exclusion_reason`.",
  "",
  "## archive_class composition", "",
  paste0("- ", names(table(frozen$archive_class)), ": ", as.integer(table(frozen$archive_class))),
  "",
  "## Archive-class assignments", "",
  "- Repasch et al. 2020 PANGAEA deposit -> fluvial (whole Rio Bermejo suspended/bed/floodplain system).",
  "- Gaviria-Lugo et al. 2023 -> Table 1 'Sediment type' (Riverine/Soils/Marine);",
  "  GeoB#### = MARUM marine core-tops.",
  "- Gensel et al. 2022 -> per-sample sub-environment from PANGAEA",
  "  doi:10.1594/PANGAEA.935586 (upper reach/floodplain/swamp/delta -> fluvial; lake -> lake).",
  "- Hren & Brandon 2026 -> soil (leaf-wax delta-2H dataset; Fig. 2 'soil sample location map').",
  "- Bai 2014 / Schwab 2015 / Lu 2020 / Feng 2019 / Jaeschke 2018 -> soil (whole-study);",
  "  leaf-wax delta-2H samples are terrestrial soils, confirmed from each paper full",
  "  text. DOI-keyed override (avoids Bai 2011/2014 and",
  "  Gaviria-Lugo/Lu homograph leakage).",
  "- Garcin et al. 2012 -> lake sediment for the coastal Debundscha row the coordinate",
  "  heuristic mis-filed as marine (Garcin = 11 lake surface sediments).",
  "",
  "## Retained curation flags", "",
  "These records remain in the audited calibration dataset; the flags document",
  "classification ambiguities rather than exclusions.", "",
  "- Hren sediment-tagged drainage rows (9): soil vs fluvial unresolved.",
  "- Gaviria zonal-mean rows (3): possible redundant aggregates of the 26 riverbed sites.",
  "",
  "## Parameters (from 2g_data_audit.R)", "",
  "- duplicate coordinate resolution: 2 dp (~1.1 km); wax tolerance: 5 permil.",
  "- coast ambiguity buffer: 25 km.",
  "",
  sprintf("Companion files: `%s_exclusions.csv`, `%s.sha256`.", DATASET, DATASET),
  "The public repository tracks the CSV. Local RDS and Parquet formats are generated when their dependencies are available."
)
writeLines(dict, file.path(frozen_dir, paste0(DATASET, "_dictionary.md")))

# ---- Console summary ---------------------------------------------------------
cat(sprintf("AUDITED: %s\n  %d source -> %d retained (-%d excluded, -%d nonfinite)\n\n",
            DATASET, n_raw, nrow(frozen), length(drop_ids), n_nonfinite))
cat("columns:", ncol(frozen), "-", paste(names(frozen), collapse = ", "), "\n")
cat("\narchive_class:\n"); print(table(frozen$archive_class))
cat("\nrows with a review flag (retained in the curation table):", n_flagged, "\n")
cat("parquet written:", parquet_ok, " | sha256 written:", sha_ok, "\n")
