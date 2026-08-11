#!/usr/bin/env Rscript

# Audit the source and DOI counts reported in Supplementary Table 1.
#
# This script deliberately produces a data artifact rather than manuscript
# LaTeX. It checks that every retained model record maps to exactly one table
# row and records the distinct primary-publication DOIs underlying each row.
# Run from the analysis-repository root:
#   Rscript scripts/audit_compilation_sources.R

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) args[[1L]] else "results/3_sediment_ready_for_modeling.csv"
output_path <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "scripts/reference_outputs/compilation_source_counts.csv"
}

required_columns <- c("source", "compilation", "DOI")
records <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
missing_columns <- setdiff(required_columns, names(records))
if (length(missing_columns) > 0L) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

normalize_doi <- function(x) {
  x <- trimws(tolower(x))
  x <- sub("^https?://(dx\\.)?doi\\.org/", "", x)
  x[nchar(x) == 0L] <- NA_character_
  x
}

records$DOI <- normalize_doi(records$DOI)
if (anyNA(records$DOI)) {
  stop(sum(is.na(records$DOI)), " retained records have no primary-publication DOI")
}

# The fit-time model-ready file retained an obsolete provenance label for the
# 27 Río Bermejo records. The hydrogen-isotope measurements come from the 2020
# PANGAEA deposit (10.1594/PANGAEA.925616), and the associated peer-reviewed
# publication is Dosch et al. (2024; 10.5194/esurf-12-907-2024), not Repasch et
# al. (2021). This correction changes source metadata only; coordinates and all
# modeled isotope/environmental values are unchanged.
is_rio_bermejo <- records$source == "Repasch et al. 2021" &
  records$DOI == "10.1038/s41561-021-00845-7"
if (sum(is_rio_bermejo) != 27L) {
  stop("Expected 27 fit-time Río Bermejo records; found ", sum(is_rio_bermejo))
}
records$source[is_rio_bermejo] <- "Repasch et al. 2020"
records$DOI[is_rio_bermejo] <- "10.5194/esurf-12-907-2024"

compilation_labels <- c(
  "ETH Zurich: doi:10.3929/ethz-b-000412154" = "ETH Zurich",
  "Hren & Brandon 2026 Zenodo: 17209981" = "Hren and Brandon (2026)",
  "Liu" = "Liu and An (2019)",
  "McFarlin" = "McFarlin et al. (2019)",
  "Manual extraction" = "Manual extraction from primary literature",
  "Ladd" = "Ladd et al. (2021)",
  "Supplement Table S1 (soil)" = "Roy and Sanyal (2022) Table S1 (soil)",
  "Wyoming Data Repository doi:10.15786/20126483" = "Wyoming Data Repository",
  "Struck et al. 2020 DataSheet1" = "Struck et al. (2020) DataSheet1",
  "NOAA WDS Paleoclimatology: study 27790" = "NOAA WDS Paleoclimatology study 27790",
  "PANGAEA: doi:10.1594/PANGAEA.859572" = "PANGAEA 10.1594/PANGAEA.859572",
  "Corcoran et al. 2022 NSF PAR preprint + Elsevier supplementary (mmc1.xlsx)" = "Corcoran et al. (2022)",
  "Roy & Sanyal 2022 Table 2" = "Roy and Sanyal (2022) Table 2"
)

direct_labels <- c(
  "Gaviria-Lugo et al. 2023" = "Gaviria-Lugo et al. (2023)",
  "Gaviria-Lugo et al., 2023" = "Gaviria-Lugo et al. (2023)",
  "Repasch et al. 2020" = "Repasch et al. (2020) PANGAEA deposit",
  "Gensel et al. 2022" = "Gensel et al. (2022)",
  "Nieto-Moreno et al. 2016" = "Nieto-Moreno et al. (2016)",
  "Wang et al., 2023" = "Wang et al. (2023)"
)

records$table_source <- unname(compilation_labels[records$compilation])
is_direct <- is.na(records$compilation) | trimws(records$compilation) == ""
records$table_source[is_direct] <- unname(direct_labels[records$source[is_direct]])

if (anyNA(records$table_source)) {
  unmapped <- unique(records[is.na(records$table_source), c("source", "compilation")])
  print(unmapped, row.names = FALSE)
  stop(nrow(unmapped), " source/compilation combinations do not map to a table row")
}

table_order <- unique(c(unname(compilation_labels), unname(direct_labels)))
if (!setequal(unique(records$table_source), table_order)) {
  stop("Mapped table rows differ from the 18 expected Supplementary Table 1 rows")
}

summarize_source <- function(label) {
  subset <- records[records$table_source == label, , drop = FALSE]
  dois <- sort(unique(subset$DOI))
  data.frame(
    table_source = label,
    n_samples = nrow(subset),
    n_primary_dois = length(dois),
    primary_dois = paste(dois, collapse = "; "),
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, lapply(table_order, summarize_source))
audit <- rbind(
  audit,
  data.frame(
    table_source = "TOTAL UNIQUE",
    n_samples = nrow(records),
    n_primary_dois = length(unique(records$DOI)),
    primary_dois = paste(sort(unique(records$DOI)), collapse = "; "),
    stringsAsFactors = FALSE
  )
)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(audit, output_path, row.names = FALSE, na = "")

message("Wrote ", output_path)
message(
  "Audited ", nrow(records), " retained records, ",
  length(table_order), " table rows, and ",
  length(unique(records$DOI)), " unique primary-publication DOIs."
)
