#!/usr/bin/env Rscript
# 2i_freeze_calibration.R
#
# Builds the frozen leaf-wax n-C29 delta-2H calibration dataset from the audited
# compilation. This is the FROZEN step: it applies the operator-approved audit
# decisions (duplicate drops + sample-level archive classification) to a new,
# versioned file. The source compilation (input_data/global_data_c29.csv) is
# read-only and untouched.
#
# Decisions applied:
#   - Drop the value-matched re-ingestion duplicates (operator_decision in
#     data/audit/duplicate_review.csv). After grounding every same-DOI candidate
#     against the primary per-sample tables (NotebookLM), only one survives:
#     Garcin 2012 Barombi Mbo. The 2g heuristic's other re-ingestion flags were
#     reversed (distinct transect/cluster samples sharing a 2-dp coordinate).
#   - Archive class: sample-level overrides (data/audit/archive_overrides_proposal.csv,
#     each grounded in the source paper / its data archive) for Repasch /
#     Gaviria-Lugo / Gensel / Hren & Brandon, plus NLM-confirmed whole-study soil
#     reclassification of Bai 2014 / Schwab 2015 / Lu 2020 / Feng 2019 /
#     Jaeschke 2018 and a Garcin 2012 coastal-lake de-marine; coordinate heuristic
#     for all others (data/audit/archive_type_proposal.csv). 4th class: "fluvial sediment".
#   - Flags (not exclusions): near-coast lake/marine ambiguity; Hren
#     sediment-tagged drainage rows; Gaviria zonal-mean rows.
#
# Run:  Rscript 2i_freeze_calibration.R   (from repo root, after 2g + 2h)
# Out:  data/frozen/leafwax_d2h_c29_calibration_v1.{csv,rds}  (+ .parquet if arrow)
#       data/frozen/leafwax_d2h_c29_calibration_v1_dictionary.md
#       data/frozen/leafwax_d2h_c29_calibration_v1_exclusions.csv

suppressWarnings(suppressMessages(library(dplyr)))

DATASET <- "leafwax_d2h_c29_calibration_v1"

repo_root <- "."   # scripts run from the repository root (see README)
audit <- file.path(repo_root, "data", "audit")
frozen_dir <- file.path(repo_root, "data", "frozen")
dir.create(frozen_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Base compilation + obs_id (same convention as 2g_data_audit.R) ---------
raw <- read.csv(file.path(repo_root, "input_data", "global_data_c29.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
raw$obs_id <- sprintf("obs_%04d", seq_len(nrow(raw)))
n_raw <- nrow(raw)

# ---- Decisions ---------------------------------------------------------------
# Duplicate drops come from the human-reviewed, value-matched decisions
# (operator_decision in duplicate_review.csv), NOT the 2g heuristic
# proposed_action. Value-matching each same-DOI candidate against the primary
# per-sample tables (via the NotebookLM "Leafwax" notebook) reversed all but one
# of the heuristic's re-ingestion flags: the colliding rows are distinct
# transect/cluster/land-use samples that merely round to the same 2-dp
# coordinate. Only Garcin 2012 Barombi Mbo (a single crater lake at that
# coordinate) is a genuine re-ingestion. See HANDOFF.md and
# Lessons/coordinate-collision-not-proof-of-duplicate.md.
rev <- read.csv(file.path(audit, "duplicate_review.csv"), stringsAsFactors = FALSE)
drop_ids <- rev$obs_id[grepl("^DROP", rev$operator_decision)]

ov <- read.csv(file.path(audit, "archive_overrides_proposal.csv"), stringsAsFactors = FALSE) %>%
  select(obs_id, archive_class_ov = archive_class, archive_evidence = evidence)

heur <- read.csv(file.path(audit, "archive_type_proposal.csv"), stringsAsFactors = FALSE) %>%
  select(obs_id, archive_type_heur = archive_type, archive_ambiguous, dist_coast_km)

stopifnot(!anyDuplicated(ov$obs_id), !anyDuplicated(heur$obs_id))

# ---- Assemble ----------------------------------------------------------------
d <- raw %>%
  left_join(ov, by = "obs_id") %>%
  left_join(heur, by = "obs_id")
stopifnot(nrow(d) == n_raw)   # joins must not fan out

# Final archive class: override wins; "FLAG"/"review" override labels are mapped
# to their underlying clean class and surfaced via archive_flag instead.
d <- d %>% mutate(
  archive_source = ifelse(!is.na(archive_class_ov), "source-override", "coord-heuristic"),
  archive_class = case_when(
    is.na(archive_class_ov)                              ~ archive_type_heur,
    grepl("^soil", archive_class_ov)                     ~ "soil",
    grepl("zonal mean", archive_class_ov)                ~ "fluvial sediment",
    TRUE                                                 ~ archive_class_ov
  )
)

# Flags (semicolon-joined; not exclusions)
flag <- function(cond, txt) ifelse(cond, txt, "")
join_flags <- function(...) {
  m <- cbind(...)
  apply(m, 1, function(r) paste(r[r != ""], collapse = "; "))
}
d <- d %>% mutate(
  archive_flag = join_flags(
    flag(archive_source == "coord-heuristic" & archive_ambiguous %in% TRUE,
         "near-coast: lake/marine ambiguous"),
    flag(!is.na(archive_class_ov) & grepl("sediment-tagged", archive_class_ov),
         "Hren sediment-tagged drainage row: soil vs fluvial unresolved"),
    flag(!is.na(archive_class_ov) & grepl("zonal mean", archive_class_ov),
         "Gaviria zonal-mean aggregate: possible redundant with riverbed sites")
  )
)

# ---- Apply exclusions --------------------------------------------------------
d <- d %>% mutate(
  .finite_ok = is.finite(latitude) & is.finite(longitude) & is.finite(d2H_wax),
  .is_drop   = obs_id %in% drop_ids)
n_nonfinite <- sum(!d$.finite_ok)

frozen <- d %>% filter(.finite_ok & !.is_drop)
n_flagged <- sum(frozen$archive_flag != "")   # capture before trimming columns

# Final column set: 13 core fields (identity/provenance, coordinates, predictor,
# response, audited medium). Audit-process metadata used to DERIVE archive_class
# (archive_source/evidence/flag/ambiguous, dist_coast_km), the unreliable raw
# sample_type, the constant chain (==29), and the near-empty d2H_precip_year are
# dropped from the released file; they remain in data/audit/ for provenance.
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
    .is_drop ~ "duplicate: non-canonical member of same-DOI re-ingestion set",
    TRUE     ~ "non-finite latitude/longitude/d2H_wax")) %>%
  select(obs_id, source, compilation, DOI, latitude, longitude, d2H_wax, exclusion_reason)
write.csv(excl_rows, file.path(frozen_dir, paste0(DATASET, "_exclusions.csv")), row.names = FALSE)

# ---- Data dictionary ---------------------------------------------------------
dict <- c(
  sprintf("# %s -- data dictionary", DATASET), "",
  "Frozen global leaf-wax n-C29 alkane delta-2H calibration compilation. Built by",
  "`2i_freeze_calibration.R` from `input_data/global_data_c29.csv`",
  "(read-only) with the audited duplicate + archive decisions applied. Journal-neutral name.",
  "",
  sprintf("- Source rows: %d", n_raw),
  sprintf("- Duplicates dropped (same-DOI re-ingestion): %d", length(drop_ids)),
  sprintf("- Non-finite lat/lon/d2H_wax dropped: %d", n_nonfinite),
  sprintf("- **Frozen rows: %d**", nrow(frozen)),
  "",
  "## Columns (13)",
  "- `obs_id` -- stable id, `obs_%04d` by source-compilation row order.",
  "- `source`, `compilation`, `location`, `DOI` -- provenance as compiled.",
  "- `latitude`, `longitude`, `elevation` -- site coordinates (WGS84) and metres.",
  "- `d2H_precip`, `d2H_precip_err` -- precipitation delta-2H predictor (permil, VSMOW).",
  "- `d2H_wax`, `d2H_wax_err` -- leaf-wax n-C29 delta-2H response (permil, VSMOW).",
  "- `archive_class` -- final depositional class: soil / lake sediment /",
  "  marine sediment / **fluvial sediment** (the audit's 4th class). Derivation",
  "  recorded below; see data/audit/ for the per-row evidence/flags.",
  "",
  "Dropped from the released file (retained in data/audit/ for provenance):",
  "`chain` (constant 29), `d2H_precip_year` (near-empty), `sample_type` (raw binary,",
  "unreliable -- superseded by `archive_class`), and the audit-process metadata",
  "`archive_source` / `archive_evidence` / `archive_flag` / `archive_ambiguous` /",
  "`dist_coast_km` used to derive `archive_class`.",
  "",
  "## archive_class composition",
  paste0("- ", names(table(frozen$archive_class)), ": ", as.integer(table(frozen$archive_class))),
  "",
  "## Source-override provenance",
  "- Repasch et al. 2021 -> fluvial (whole Rio Bermejo suspended/bed/floodplain system).",
  "- Gaviria-Lugo et al. 2023 -> Table 1 'Sediment type' (Riverine/Soils/Marine);",
  "  GeoB#### = MARUM marine core-tops.",
  "- Gensel et al. 2022 -> per-sample sub-environment from PANGAEA",
  "  doi:10.1594/PANGAEA.935586 (upper reach/floodplain/swamp/delta -> fluvial; lake -> lake).",
  "- Hren & Brandon 2026 -> soil (leaf-wax delta-2H dataset; Fig. 2 'soil sample location map').",
  "- Bai 2014 / Schwab 2015 / Lu 2020 / Feng 2019 / Jaeschke 2018 -> soil (whole-study);",
  "  leaf-wax delta-2H samples are terrestrial soils, confirmed from each paper full",
  "  text via NotebookLM (2026-06-23). DOI-keyed override (avoids Bai 2011/2014 and",
  "  Gaviria-Lugo/Lu homograph leakage).",
  "- Garcin et al. 2012 -> lake sediment for the coastal Debundscha row the coordinate",
  "  heuristic mis-filed as marine (Garcin = 11 lake surface sediments, NLM).",
  "",
  "## Flags (review items, NOT exclusions)",
  "- Hren sediment-tagged drainage rows (9): soil vs fluvial unresolved.",
  "- Gaviria zonal-mean rows (3): possible redundant aggregates of the 26 riverbed sites.",
  "",
  "## Parameters (from 2g_data_audit.R)",
  "- duplicate coordinate resolution: 2 dp (~1.1 km); wax tolerance: 5 permil.",
  "- coast ambiguity buffer: 25 km.",
  "",
  sprintf("Companion files: `%s_exclusions.csv`, `%s.sha256`.", DATASET, DATASET),
  if (!parquet_ok) "(parquet not written -- `arrow` package unavailable; csv + rds only.)" else
    "Formats: csv, rds, parquet."
)
writeLines(dict, file.path(frozen_dir, paste0(DATASET, "_dictionary.md")))

# ---- Console summary ---------------------------------------------------------
cat(sprintf("FROZEN: %s\n  %d source -> %d frozen (-%d dup, -%d nonfinite)\n\n",
            DATASET, n_raw, nrow(frozen), length(drop_ids), n_nonfinite))
cat("columns:", ncol(frozen), "-", paste(names(frozen), collapse = ", "), "\n")
cat("\narchive_class:\n"); print(table(frozen$archive_class))
cat("\nrows with an audit flag (dropped from file; counted from data/audit/):", n_flagged, "\n")
cat("parquet written:", parquet_ok, " | sha256 written:", sha_ok, "\n")
