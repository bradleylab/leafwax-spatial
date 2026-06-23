#!/usr/bin/env Rscript
# 2h_archive_overrides.R
#
# Proposes a sample-level archive-type override for the four sources whose
# leaf-wax δ²H samples are NOT cleanly captured by the coordinate heuristic +
# the compilation's binary `sample_type` (Sediment/Soil). Each assignment is
# grounded in an authoritative source, recorded in `evidence`:
#   - Repasch et al. 2021  : whole source is the Río Bermejo fluvial system.
#   - Gaviria-Lugo et al. 2023 : Table 1 "Sediment type" column (Riverine /
#       Soils / Marine), recovered verbatim via NotebookLM (2026-06-22).
#   - Gensel et al. 2022   : per-sample sub-environment from the PANGAEA data
#       set doi:10.1594/PANGAEA.935586 (downloaded to data/audit/external/).
#   - Hren & Brandon 2026  : leaf-wax δ²H samples are soils (paper Fig. 2,
#       "Soil sample location map for n-alkane δ²H data"); 9 rows tagged by
#       river drainage are flagged for the operator, not silently reclassified.
#
# PROPOSAL ONLY. Non-destructive: reads input_data/global_data_c29.csv (read
# only) + the PANGAEA file, writes data/audit/archive_overrides_proposal.csv.
# Nothing is applied to the compilation. Re-runnable.
#
# Run:  Rscript 2h_archive_overrides.R   (from repo root, after 2g_data_audit.R)

suppressWarnings(suppressMessages(library(dplyr)))

repo_root <- "."   # scripts run from the repository root (see README)

d <- read.csv(file.path(repo_root, "input_data", "global_data_c29.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
# obs_id by input row order -- identical convention to 2g_data_audit.R, so the
# freeze step can join overrides on obs_id.
d$obs_id <- sprintf("obs_%04d", seq_len(nrow(d)))

norm <- function(x) tolower(gsub("[^a-z0-9]", "", x, ignore.case = TRUE))

# ---- Gensel: MK id -> sub-environment, from the PANGAEA data set -------------
# Self-contained: download the data set (doi:10.1594/PANGAEA.935586) if absent.
pan <- file.path(repo_root, "data", "audit", "external", "gensel_2021_pangaea.tsv")
if (!file.exists(pan)) {
  dir.create(dirname(pan), showWarnings = FALSE, recursive = TRUE)
  download.file("https://doi.pangaea.de/10.1594/PANGAEA.935586?format=textfile",
                pan, quiet = TRUE, mode = "wb")
}
L <- readLines(pan, warn = FALSE)
endc <- which(trimws(L) == "*/")[1]
gtab <- read.delim(text = paste(L[(endc + 1):length(L)], collapse = "\n"),
                   sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
ev_col <- grep("Event", names(gtab), value = TRUE)[1]
sc_col <- grep("Sample comment", names(gtab), value = TRUE)[1]
gmap <- setNames(gtab[[sc_col]], gtab[[ev_col]])   # "MK01-10" -> "upper reach" etc.

# ---- Gaviria-Lugo Table 1 soil sites (closed list of 12) ---------------------
# Riverine sites use a "Rio/Quebrada/Estero" prefix; marine sites are "GeoB...";
# these 12 bare names are the topsoils. (Gaviria-Lugo et al. 2023, Table 1.)
gav_soil <- norm(c("Choros", "Talca A", "Talca B", "Cajon delMaipo",
                   "SanAntonio - Maipo", "Rapel", "BioBio A", "BioBio B",
                   "CalleCalle", "Bueno", "Maullin A", "Maullin B"))

classify <- function(source, location, sample_type) {
  loc <- location
  if (grepl("Repasch", source)) {
    return(c("fluvial sediment", "Repasch 2021: entire source is Rio Bermejo (suspended load, bedload, floodplain)."))
  }
  if (grepl("Gaviria", source)) {
    if (grepl("^GeoB", loc))                 return(c("marine sediment", "Gaviria 2023 Table 1: GeoB MARUM core-top = Marine."))
    if (grepl("zone .*mean", loc, ignore.case = TRUE))
      return(c("review (zonal mean)", "Gaviria 2023: zonal-mean row, underlying medium not a single site."))
    if (norm(loc) %in% gav_soil)             return(c("soil", "Gaviria 2023 Table 1: bare catchment name = Soils (topsoil)."))
    return(c("fluvial sediment", "Gaviria 2023 Table 1: Rio/Quebrada/Estero site = Riverine (riverbed)."))
  }
  if (grepl("Gensel", source)) {
    mk <- sub(".*\\(([^)]+)\\).*", "\\1", loc)        # extract "MKxx-y"
    sub_env <- gmap[[mk]]
    if (is.null(sub_env) || is.na(sub_env)) return(c("review (MK unmatched)", paste0("Gensel: MK id '", mk, "' not found in PANGAEA 935586.")))
    if (sub_env == "lake")                  return(c("lake sediment", paste0("Gensel PANGAEA 935586: ", mk, " = lake (Lake St Lucia).")))
    return(c("fluvial sediment", paste0("Gensel PANGAEA 935586: ", mk, " = ", sub_env, " (river/wetland).")))
  }
  if (grepl("Hren", source)) {
    # All Hren rows are drainage-labeled; the 9 sample_type=Sediment rows are
    # the catchment/river-tagged ones worth a second look. The other 182 soils.
    if (identical(sample_type, "Sediment"))
      return(c("soil (FLAG: sediment-tagged)", "Hren & Brandon 2026: leaf-wax δ²H = soils; this row is sample_type=Sediment (river-drainage) -- operator to confirm soil vs fluvial."))
    return(c("soil", "Hren & Brandon 2026: leaf-wax δ²H dataset is soils (Fig. 2)."))
  }
  c(NA_character_, NA_character_)
}

targets <- d %>% filter(grepl("Repasch|Gaviria|Gensel|Hren", source))
res <- t(mapply(classify, targets$source, targets$location, targets$sample_type))
targets$archive_class <- res[, 1]
targets$evidence <- res[, 2]

out <- targets %>%
  mutate(heuristic_sample_type = sample_type) %>%
  select(obs_id, source, compilation, location, latitude, longitude, chain,
         d2H_wax, heuristic_sample_type, archive_class, evidence) %>%
  arrange(source, location)

write.csv(out, file.path(repo_root, "data", "audit", "archive_overrides_proposal.csv"),
          row.names = FALSE)

# ---- Verification ------------------------------------------------------------
cat("Wrote data/audit/archive_overrides_proposal.csv  (", nrow(out), "rows )\n\n")
cat("=== archive_class by source ===\n")
print(addmargins(table(out$source, out$archive_class)))
cat("\n=== Gaviria-Lugo expected 26 fluvial / 12 soil / 29 marine / 3 review ===\n")
print(table(out$archive_class[grepl("Gaviria", out$source)]))
cat("\n=== Gensel split (lake vs fluvial) ===\n")
print(table(out$archive_class[grepl("Gensel", out$source)]))
