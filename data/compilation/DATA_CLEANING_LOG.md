# Data Cleaning Log: Global C29 n-Alkane d2H Dataset

**Input file:** `global_data_c29_merged.csv` (1092 rows)
**Output file:** `global_data_c29_cleaned.csv` (1002 rows)
**Date:** 2026-04-06
**Script:** `clean_dataset.py`

---

## 1. Typo Correction: Wang et al. (2017)

One transcription error was identified and corrected.

| Field | Original | Corrected |
|---|---|---|
| Source | Wang et al., 2017 | Wang et al., 2017 |
| Location | Tibetan plateau | Tibetan plateau |
| Latitude | 22.7 | 22.7 |
| Longitude | 101.3 | 101.3 |
| Elevation (m) | 956 | 956 |
| d2H_wax (permil) | **-288.0** | **-188.0** |

**Evidence chain:**

1. The Hren & Brandon (2026) Zenodo dataset (DOI: 10.5281/zenodo.17209981), Table S3, identifies this as sample XIS-13-03 and reports d2H C29 = -188.2 permil.
2. The site is at 956 m elevation in subtropical Yunnan province, China.
3. Flanking sites from the same transect report d2H_wax values of -213 permil (931 m, lat 22.3) and -210 permil (1400 m, lat 22.2).
4. A value of -288 permil would be more depleted than sites above 4000 m on the Tibetan Plateau in the same dataset (e.g., -241 permil at 4226 m). This is physically implausible for a low-elevation subtropical site.
5. The error is consistent with a transcription typo where the leading digit "1" was replaced with "2" (i.e., -188 became -288).

---

## 2. Deduplication

### 2.1 Definitions

- **Cross-study duplicate:** The same physical sample appearing in two different publications. One study reused data from another. Identified by coordinates matching within 0.01 degrees AND d2H_wax matching within 2 permil AND different source/reference fields.
- **Within-study duplicate:** The same measurement appearing multiple times within one study's dataset, likely from data entry errors or overlapping compilation entries. Identified by same source, coordinates within 0.01 degrees, and d2H_wax within 2 permil.
- **Same-site replicate:** Multiple independent samples collected at the same location, each with a distinct d2H_wax value. These capture real biological and analytical variability and are retained in the dataset.

### 2.2 Matching Criteria

| Parameter | Threshold |
|---|---|
| Coordinate tolerance | 0.01 degrees (~1.1 km) |
| d2H_wax tolerance | 2 permil |
| Replicate threshold | >5 permil range within cluster |

### 2.3 Cross-Study Duplicates Removed (4 rows)

These are cases where the same sample was published in two different studies. The row from the original measurement source was retained.

| Removed Source | Removed Location | Lat | Lon | d2H_wax | Kept Source | Kept Row |
|---|---|---|---|---|---|---|
| Garcin et al., 2012 | MANE | 5.03 | 9.83 | -184.83 | Garcin et al 2012 | 189 |
| Sachse et al 2004 | LPM | 40.93 | 15.61 | -169.0 | Leider et al 2013 | 365 |
| Sachse et al 2004 | HZM | 50.12 | 6.88 | -198.0 | Mugler et al 2008 | 488 |
| Xia et al 2008 | NC | 30.64 | 90.6 | -246.0 | Mugler et al 2008 | 489 |

Notes:
- Garcin et al 2012 / Garcin et al., 2012 represent the same study with slightly different citation formatting in the two compilation sources. The row with d2H_wax_err reported (4.0 permil) was retained.
- Sachse et al. (2004) data for Lago di Monticchio and Holzmaar were later re-reported in Leider et al. (2013) and Mugler et al. (2008), respectively.
- Xia et al. (2008) Tibet data overlaps with Mugler et al. (2008). The Mugler row (which reports d2H_wax_err = 36.4) was retained.

### 2.4 Within-Study Duplicates Removed (73 rows)

These are rows from the same study with coordinates within 0.01 degrees and d2H_wax within 2 permil. For each cluster, one representative row was retained.

| Cluster | Source | Location | Rows Removed | d2H Range (permil) | Kept Row |
|---|---|---|---|---|---|
| 1 | Bai et al., 2011 | Mt Kunlun | 2 | 1.0 | 23 |
| 2 | Bai et al., 2014 | Mt Qilian | 1 | 0.0 | 41 |
| 3 | Bai et al., 2014 | Mt Qilian | 1 | 2.0 | 46 |
| 4 | Bai et al., 2014 | Mt Qilian | 1 | 0.0 | 48 |
| 5 | Bai et al., 2014 | Mt Qilian | 1 | 1.0 | 56 |
| 6 | Bai et al., 2014 | Mt Qilian | 1 | 1.0 | 57 |
| 7 | Bai et al., 2014 | Mt Qilian | 2 | 2.0 | 63 |
| 8 | Bai et al., 2014 | Mt Qilian | 1 | 0.0 | 68 |
| 9 | Chikaraishi and Naraoka, 2006 | Gunma Prefecture | 2 | 2.0 | 75 |
| 10 | Daniels et al 2017 | S6/S7 | 1 | 0.0 | 80 |
| 11 | Daniels et al 2017 | I7/I5 | 1 | 1.0 | 83 |
| 12 | Daniels et al 2017 | I8/I1 | 1 | 0.0 | 84 |
| 13 | Daniels et al 2017 | FOG1/FOG2 | 1 | 2.0 | 86 |
| 14 | Daniels et al 2017 | I3/I4/I6 | 2 | 1.0 | 89 |
| 15 | Feng et al., 2019 | JS7/JS5 | 1 | 0.0 | 136 |
| 16 | Feng et al., 2019 | JH3/JH4 | 1 | 0.0 | 147 |
| 17 | Feng et al., 2019 | LJ3/LJ4 | 1 | 0.0 | 149 |
| 18 | Feng et al., 2019 | LJ10/LJ14/LJ15 | 2 | 0.0 | 158 |
| 19 | Feng et al., 2019 | LJ11/LJ12 | 1 | 2.0 | 159 |
| 21 | Herrmann et al., 2017 | WRZ | 1 | 2.0 | 211 |
| 22 | Herrmann et al., 2017 | SRZ | 1 | 2.0 | 254 |
| 23 | Jaeschke et al., 2018 | VI sites | 2 | 3.0 | 278 |
| 27 | Jia et al., 2008 | Mt Gongga | 1 | 1.0 | 310 |
| 28 | Jia et al., 2008 | Mt Gongga | 6 | 1.0 | 313 |
| 29 | Jia et al., 2008 | Mt Gongga | 1 | 1.0 | 336 |
| 30 | Jia et al., 2008 | Mt Gongga | 1 | 1.0 | 340 |
| 32 | Lu et al., 2020 | JY-04/JY-05 | 1 | 1.0 | 446 |
| 33 | Lu et al., 2020 | FS-01/FS-03 | 1 | 1.0 | 451 |
| 34 | Lu et al., 2020 | CB-50/CB-02 | 1 | 1.0 | 454 |
| 35 | Lu et al., 2020 | ET-02/ET-03 | 1 | 1.0 | 456 |
| 36 | Lu et al., 2020 | ABQ-03/ABQ-05 | 1 | 2.0 | 458 |
| 37 | Lu et al., 2020 | SNQ-02/SNQ-01 | 1 | 1.0 | 463 |
| 38 | Lu et al., 2020 | XM-02/XM-01 | 1 | 1.0 | 465 |
| 41 | Pautler et al., 2014 | Christie Mine | 2 | 1.0 | 497 |
| 43 | Pautler et al., 2014 | Goldbottom Creek | 1 | 1.0 | 517 |
| 44 | Pautler et al., 2014 | Kluane Lake | 1 | 1.0 | 520 |
| 45 | Peterse et al., 2009 | Mt. Kilimanjaro | 4 | 5.0 | 524 |
| 47 | Peterse et al., 2009 | Mt. Kilimanjaro | 3 | 1.0 | 529 |
| 48 | Polissar and Freeman 2010 | Laguna Misteque/Royal C | 1 | 1.8 | 557 |
| 50 | Seki et al., 2010 | Hokkaido, Japan | 1 | 1.0 | 591 |
| 51 | Seki et al., 2010 | Hokkaido, Japan | 1 | 0.0 | 592 |
| 53 | Struck et al., 2020 | UGC3/UGC12 | 1 | 0.0 | 624 |
| 54 | Struck et al., 2020 | TLC2/TLC3 | 1 | 1.9 | 635 |
| 55 | Struck et al., 2020 | TLC12/TLC5 | 1 | 0.0 | 639 |
| 56 | van der Veen et al., 2020 | S127 sp1/sp2 | 1 | 0.87 | 667 |
| 57 | Wang et al., 2023 | Techeng T1/T2 | 1 | 0.0 | 827 |
| 58 | Wang et al., 2023 | Techeng T4/T6 | 1 | 0.0 | 830 |
| 59 | Corcoran et al., 2022 | Tarn 2/4 | 3 | 4.0 | 833 |
| 60 | Hren & Brandon 2026 | Himalaya, Brahmaputra | 2 | 1.0 | 845 |
| 61 | Hren & Brandon 2026 | Himalaya, Brahmaputra | 1 | 0.0 | 860 |
| 62 | Hren & Brandon 2026 | Himalaya, Brahmaputra | 1 | 0.0 | 895 |
| 63 | Hren & Brandon 2026 | Himalaya, Brahmaputra | 1 | 1.0 | 903 |
| 64 | Sisley & Wolhowe (2023) | Cascades, Leader Lake | 1 | 1.7 | 1033 |

### 2.5 Deduplication Summary

| Category | Rows Removed |
|---|---|
| Cross-study duplicates | 4 |
| Within-study duplicates | 73 |
| **Total removed** | **77** |

---

## 3. Same-Site Replicates Retained

Seven duplicate clusters were identified where rows from the same source had d2H_wax values spanning >5 permil. These represent independent samples collected at the same (or nearby) location and capture real biological and analytical variability. All rows in these clusters were retained.

Coordinates in compilation datasets are often rounded or jittered during preprocessing, which can cause genuinely distinct sampling points to share identical coordinates. The d2H spread confirms these are independent measurements, not data entry errors.

| Cluster | Source | Location | Rows Retained | d2H Range (permil) |
|---|---|---|---|---|
| 24 | Jaeschke et al., 2018 | Jimma transect (sites I-VI) | 10 | 6.0 |
| 25 | Jia et al., 2008 | Mt Gongga | 8 | 7.0 |
| 26 | Jia et al., 2008 | Mt Gongga | 22 | 11.0 |
| 42 | Pautler et al., 2014 | Tatlow Camp, Quartz Creek | 11 | 13.0 |
| 46 | Peterse et al., 2009 | Mt. Kilimanjaro | 6 | 9.0 |
| 49 | Seki et al., 2010 | Hokkaido, Japan | 14 | 10.0 |
| 52 | Seki et al., 2010 | Hokkaido, Japan | 12 | 13.0 |

**Total replicate rows retained: 83**

These within-site replicates are valuable for quantifying the natural variability of leaf wax d2H at a given location, which is directly relevant to the calibration uncertainty addressed in this study.

---

## 4. Final Dataset Summary

| Metric | Value |
|---|---|
| Rows before cleaning | 1092 |
| Typo corrections | 1 (value changed, row retained) |
| Rows removed (deduplication) | 77 |
| **Final row count** | **1015** |
| Unique publications | 54 |
| Latitude range | -51.70 to 76.85 |
| Longitude range | -179.94 to 168.59 |
| d2H_wax range | -285.0 to -75.0 permil |

### Source Breakdown by Data Origin

| Data Source | Rows |
|---|---|
| original (pre-existing compilations) | 767 |
| v7_new (newly added in v7) | 248 |

### Top 10 Contributing Publications

| Publication | Rows |
|---|---|
| Hren & Brandon 2026 | 143 |
| Gaviria-Lugo et al., 2023 | 67 |
| Herrmann et al., 2017 | 60 |
| Li et al., 2019 | 54 |
| Bai et al., 2014 | 41 |
| Jia et al., 2008 | 38 |
| Struck et al., 2020 | 36 |
| Douglas et al 2012 | 35 |
| Wang et al., 2017 | 35 |
| Sisley & Wolhowe (2023) | 33 |

---

## 5. Removal of Rows with Missing Critical Fields

13 rows were removed because they lacked latitude, longitude, or d2H_wax — fields required by the spatial model.

| Source | Location | Missing field | Rows |
|---|---|---|---|
| Schwab et al. 2015 | Cameroon | longitude | 10 |
| Krull et al. 2006 | Queensland | d2H_wax | 2 |
| Wang et al. 2017 | Tibetan plateau | latitude, longitude | 1 |

All 13 rows were present in the original dataset and would have been filtered during preprocessing (`3_prep_data.R`, line 40: filter for non-NA d2H_wax, latitude, longitude). They are removed here for clarity.

---

## 6. Final Dataset Summary

| Metric | Value |
|---|---|
| Total rows | 1002 |
| From original dataset | 754 |
| New from v7 compilation | 248 |
| Unique publications | ~54 |
| Latitude range | -51.70 to 76.85 |
| Longitude range | -159.37 to 176.60 |
| d2H_wax range | -285.0 to -75.0 permil |
| Rows with d2H_wax_err reported | ~520 |

### Cleaning steps applied (in order)

1. Typo correction: 1 value changed (-288 → -188)
2. Cross-study duplicates removed: 4 rows
3. Within-study exact duplicates removed: 73 rows
4. Rows with missing lat/lon/d2H_wax removed: 13 rows
5. Same-site replicates (different d2H values at shared coordinates): **retained**

**Net change from original 837-row dataset:** +165 rows (248 new sites added, 83 rows removed through deduplication and missing-data filtering).

---

## 7. V8 Integration

**Input:** `modern_leaf_wax_d2H_C29_merged_v8_fixed.csv` (1192 rows)
**Output:** `global_data_c29_final.csv` (1131 rows)
**Script:** `build_final_dataset.py`
**Date:** 2026-04-03

The v8_fixed compilation is a superset of the v7-derived `global_data_c29_cleaned.csv` (1002 rows) with verified values from a re-audit against the primary supplements. It was merged into the model-ready format directly, replacing the v7 pipeline output.

### 7.1 Changes from v7 audit (v8_fixed corrections)

v8 incorporates the following verified corrections from the v7-era `global_data_c29_cleaned.csv`. Values reflect what is now present in the final dataset; full row-level before/after diffs are preserved in the v8_fixed file and its audit flags column.

| Study | Nature of correction |
|---|---|
| Nieto-Moreno et al. | δ²H n-C29, latitude/longitude, and precipitation values re-extracted from the supplement; audit flags cleared. |
| Li et al. (2019) | δ²H n-C29 values verified against Table S1; small number of typo/transposition fixes. |
| Repasch et al. | Precipitation δ²H reconciled to the supplement's reported values; site coordinates verified. |
| Corcoran et al. | Same-site replicate structure preserved; per-sample δ²H values re-read from the source table. |
| Struck et al. | δ²H n-C29 corrected for several UGC/TLC transect rows. |

The Wang et al. (2017) Tibetan plateau row (lat 22.7, lon 101.3, 956 m) still carries δ²H n-C29 = -288.0‰ in `v8_fixed`, which was flagged in Section 1 as a transcription typo (correct value -188.2‰ per Hren & Brandon 2026 Zenodo Table S3). The value was NOT edited during this merge because the scope of this step is column mapping and deduplication only; the typo correction should be re-applied in the v8 source file before final model fitting.

### 7.2 Column mapping

The v8_fixed columns were renamed/mapped to the schema consumed by `3_prep_data.R`:

| v8 column | Final column |
|---|---|
| Reference | source |
| Compilation | compilation |
| Location | location |
| Latitude (°N) | latitude |
| Longitude (°E) | longitude |
| Elevation (m) | elevation |
| δ²H Precipitation (‰ VSMOW) | d2H_precip |
| δ²H Precipitation SD (±‰) | d2H_precip_err |
| (none) | d2H_precip_year (left blank) |
| (constant) | chain = 29 |
| δ²H n-C29 (‰ VSMOW) | d2H_wax |
| δ²H n-C29 SD (±‰) | d2H_wax_err |
| Origin | data_source |
| DOI | DOI |
| Audit Flags | audit_flags |
| Sample Type | sample_type |

The first twelve columns are the R pipeline contract; the last four are preserved for audit and will be ignored by `3_prep_data.R`.

### 7.3 Filtering and deduplication on v8_fixed

Applied in order:

1. **Missing required fields:** 8 rows dropped (2 missing latitude, 8 missing longitude, 0 missing d2H_wax). 1192 → 1184.
2. **Chain check:** all 1184 rows verified as C29.
3. **Clustering:** 47 multi-row clusters identified using coordinate tolerance 0.01° and d2H tolerance 2‰.
4. **Cross-study duplicates removed:** 32 rows (different source, within tolerances — keeper chosen from rows with `data_source == "User compilation"` when available).
5. **Within-study duplicates removed:** 21 rows (same source, within tolerances — keeper chosen from rows with reported `d2H_wax_err` when available).
6. **Same-site replicates retained:** 45 rows across 5 clusters where the d2H range exceeded 5‰, reflecting genuine biological/analytical variability.

### 7.4 Final summary (v8 integration)

| Metric | v7 cleaned (previous) | v8 final (this step) |
|---|---|---|
| Starting rows | 1092 | 1192 |
| Dropped for missing fields | 13 | 8 |
| Cross-study duplicates removed | 4 | 32 |
| Within-study duplicates removed | 73 | 21 |
| Same-site replicates retained | 83 rows / 7 clusters | 45 rows / 5 clusters |
| Final row count | 1002 | **1131** |
| Latitude range | -51.70 to 76.85 | -51.70 to 76.85 |
| Longitude range | -159.37 to 176.60 | -179.94 to 168.59 |
| d2H_wax range (‰) | -285.0 to -75.0 | -288.0 to -75.0 |
| Unique publications | ~54 | 74 |

The larger count of cross-study duplicates in v8 (32 vs 4) reflects the broader compilation: v8 merges the "User compilation", "Computer v6", and "Gap-fill v8" sources, so rows that were previously in only one source now appear in multiple and get flagged.

### 7.5 data_source composition of final dataset

| data_source (Origin) | Rows |
|---|---|
| Both | 402 |
| Computer v6 | 377 |
| User compilation | 288 |
| Gap-fill v8 | 64 |
| **Total** | **1131** |

### 7.6 Geographic composition (coarse bounding boxes)

| Region | Rows |
|---|---|
| Asia | 581 |
| N America | 193 |
| S America | 153 |
| Africa | 143 |
| Europe | 30 |
| Other (boundary/ocean-basin) | 15 |
| Oceania | 10 |
| Arctic (lat ≥ 72°N) | 6 |

Region labels are coarse lat/lon bounding boxes and are for audit only; they do not affect modelling.
