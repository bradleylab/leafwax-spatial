# Audited leaf-wax calibration data dictionary

Audited global leaf-wax *n*-C29 alkane δ²H calibration compilation. Built by
`2i_freeze_calibration.R` from `input_data/global_data_c29.csv`
(read-only) and `input_data/calibration_curation_v1.csv`.

- Source rows: 1136
- Duplicates dropped (same-DOI re-ingestion): 1
- Non-finite lat/lon/d2H_wax dropped: 0
- **Audited calibration rows: 1135**

## Columns (13)

- `obs_id` -- stable id, `obs_%04d` by source-compilation row order.
- `source`, `compilation`, `location`, `DOI` -- provenance as compiled.
- `latitude`, `longitude`, `elevation` -- site coordinates (WGS84) and metres.
- `d2H_precip`, `d2H_precip_err` -- precipitation δ²H predictor (‰, VSMOW).
- `d2H_wax`, `d2H_wax_err` -- leaf-wax *n*-C29 δ²H response (‰, VSMOW).
- `archive_class` -- final depositional class: soil / lake sediment /
  marine sediment / **fluvial sediment**. Per-row assignments and review
  flags are recorded in `input_data/calibration_curation_v1.csv`.

Omitted from the audited calibration file:
`chain` (constant 29), `d2H_precip_year` (near-empty), `sample_type` (a coarse
source field not used after archive-class curation), and the curation columns
`archive_assignment`, `review_flag`, and `exclusion_reason`.

## archive_class composition

- fluvial sediment: 75
- lake sediment: 490
- marine sediment: 52
- soil: 518

## Archive-class assignments

- Repasch et al. 2020 PANGAEA deposit -> fluvial (whole Rio Bermejo suspended/bed/floodplain system).
- Gaviria-Lugo et al. 2023 -> Table 1 'Sediment type' (Riverine/Soils/Marine);
  GeoB#### = MARUM marine core-tops.
- Gensel et al. 2022 -> per-sample sub-environment from PANGAEA
  doi:10.1594/PANGAEA.935586 (upper reach/floodplain/swamp/delta -> fluvial; lake -> lake).
- Hren & Brandon 2026 -> soil (leaf-wax delta-2H dataset; Fig. 2 'soil sample location map').
- Bai 2014 / Schwab 2015 / Lu 2020 / Feng 2019 / Jaeschke 2018 -> soil (whole-study);
  leaf-wax delta-2H samples are terrestrial soils, confirmed from each paper full
  text. DOI-keyed override (avoids Bai 2011/2014 and
  Gaviria-Lugo/Lu homograph leakage).
- Garcin et al. 2012 -> lake sediment for the coastal Debundscha row the coordinate
  heuristic mis-filed as marine (Garcin = 11 lake surface sediments).

## Retained curation flags

These records remain in the audited calibration dataset; the flags document
classification ambiguities rather than exclusions.

- Hren sediment-tagged drainage rows (9): soil vs fluvial unresolved.
- Gaviria zonal-mean rows (3): possible redundant aggregates of the 26 riverbed sites.

## Parameters (from 2g_data_audit.R)

- duplicate coordinate resolution: 2 dp (~1.1 km); wax tolerance: 5 permil.
- coast ambiguity buffer: 25 km.

Companion files: `leafwax_d2h_c29_calibration_v1_exclusions.csv`, `leafwax_d2h_c29_calibration_v1.sha256`.
The public repository tracks the CSV. Local RDS and Parquet formats are generated when their dependencies are available.
