"""
Build the final model-ready dataset from the verified v8_fixed CSV.

Pipeline:
    v8_fixed.csv  -> column mapping  -> drop missing -> deduplicate -> final.csv

Dedup rules (match clean_dataset.py semantics):
    - Within-study duplicates: same source, coords within 0.01 deg, d2H within
      2 permil -> keep one (the row with the smallest index, preferring
      non-missing d2H_wax_err).
    - Cross-study duplicates: different sources, coords within 0.01 deg, d2H
      within 2 permil -> keep a row whose data_source is "User compilation"
      (fall back to first row if none).
    - Same-site replicates: same coords, d2H range > 5 permil -> retain all.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/_for_GCA")
V8_PATH = BASE / "modern_leaf_wax_d2H_C29_merged_v8_fixed.csv"
OUTPUT_PATH = BASE / "global_data_c29_final.csv"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
COORD_TOL_DEG = 0.01
D2H_TOL_PERMIL = 2.0
REPLICATE_THRESHOLD_PERMIL = 5.0
CHAIN_LENGTH = 29

MODEL_COLUMNS = [
    "source",
    "compilation",
    "location",
    "latitude",
    "longitude",
    "elevation",
    "d2H_precip",
    "d2H_precip_err",
    "d2H_precip_year",
    "chain",
    "d2H_wax",
    "d2H_wax_err",
]

EXTRA_COLUMNS = [
    "data_source",
    "DOI",
    "audit_flags",
    "sample_type",
]

COLUMN_RENAME = {
    "Reference": "source",
    "Compilation": "compilation",
    "Location": "location",
    "Latitude (\u00b0N)": "latitude",
    "Longitude (\u00b0E)": "longitude",
    "Elevation (m)": "elevation",
    "\u03b4\u00b2H Precipitation (\u2030 VSMOW)": "d2H_precip",
    "\u03b4\u00b2H Precipitation SD (\u00b1\u2030)": "d2H_precip_err",
    "\u03b4\u00b2H n-C29 (\u2030 VSMOW)": "d2H_wax",
    "\u03b4\u00b2H n-C29 SD (\u00b1\u2030)": "d2H_wax_err",
    "Origin": "data_source",
    "Audit Flags": "audit_flags",
    "Sample Type": "sample_type",
}

REGIONS = {
    "N America": dict(lat_min=15.0, lat_max=72.0, lon_min=-170.0, lon_max=-50.0),
    "S America": dict(lat_min=-56.0, lat_max=15.0, lon_min=-90.0, lon_max=-34.0),
    "Europe": dict(lat_min=34.0, lat_max=72.0, lon_min=-25.0, lon_max=45.0),
    "Africa": dict(lat_min=-35.0, lat_max=34.0, lon_min=-20.0, lon_max=55.0),
    "Asia": dict(lat_min=-10.0, lat_max=72.0, lon_min=45.0, lon_max=180.0),
    "Oceania": dict(lat_min=-50.0, lat_max=-10.0, lon_min=110.0, lon_max=180.0),
    "Arctic": dict(lat_min=72.0, lat_max=90.0, lon_min=-180.0, lon_max=180.0),
}


def assign_region(lat: float, lon: float) -> str:
    """Coarse region label based on lat/lon bounding boxes."""
    if pd.isna(lat) or pd.isna(lon):
        return "Unknown"
    if lat >= 72.0:
        return "Arctic"
    for name, box in REGIONS.items():
        if name == "Arctic":
            continue
        if (box["lat_min"] <= lat <= box["lat_max"]
                and box["lon_min"] <= lon <= box["lon_max"]):
            return name
    return "Other"


# ---------------------------------------------------------------------------
# Step 1: load and map columns
# ---------------------------------------------------------------------------
df_raw = pd.read_csv(V8_PATH)
n_start = len(df_raw)
print(f"Starting rows (v8_fixed): {n_start}")

df = df_raw.rename(columns=COLUMN_RENAME).copy()
df["chain"] = CHAIN_LENGTH
df["d2H_precip_year"] = np.nan  # no source for year field in v8

# Coerce numerics
for col in ["latitude", "longitude", "elevation",
            "d2H_precip", "d2H_precip_err", "d2H_wax", "d2H_wax_err"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

n_after_map = len(df)
print(f"After column mapping: {n_after_map}")

# ---------------------------------------------------------------------------
# Step 2: drop rows missing required fields
# ---------------------------------------------------------------------------
missing_mask = (
    df["latitude"].isna()
    | df["longitude"].isna()
    | df["d2H_wax"].isna()
)
n_missing_lat = df["latitude"].isna().sum()
n_missing_lon = df["longitude"].isna().sum()
n_missing_wax = df["d2H_wax"].isna().sum()
n_dropped_missing = int(missing_mask.sum())

df = df.loc[~missing_mask].reset_index(drop=True)
print(
    f"Dropped rows with missing lat/lon/d2H_wax: {n_dropped_missing} "
    f"(missing latitude={n_missing_lat}, longitude={n_missing_lon}, "
    f"d2H_wax={n_missing_wax})"
)

# ---------------------------------------------------------------------------
# Step 3: verify chain == 29
# ---------------------------------------------------------------------------
non_c29 = (df["chain"] != CHAIN_LENGTH).sum()
assert non_c29 == 0, f"Found {non_c29} rows with chain != {CHAIN_LENGTH}"
print(f"All {len(df)} rows are C29 (verified)")

# ---------------------------------------------------------------------------
# Step 4: cluster-based deduplication
# ---------------------------------------------------------------------------
# Build clusters via a union-find over pairs that satisfy:
#   |lat_i - lat_j| <= COORD_TOL_DEG
#   |lon_i - lon_j| <= COORD_TOL_DEG
#   |d2H_wax_i - d2H_wax_j| <= D2H_TOL_PERMIL
# We restrict the scan to rows with identical rounded coordinates (3 dp) to
# keep this O(n * avg_bucket_size) while still catching everything within
# 0.01 deg. To be safe we expand each rounded key to its 9 neighbors.
n = len(df)
parent = list(range(n))


def find(x: int) -> int:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


def union(a: int, b: int) -> None:
    ra, rb = find(a), find(b)
    if ra != rb:
        parent[rb] = ra


lat_key = (df["latitude"] / COORD_TOL_DEG).round().astype(int).to_numpy()
lon_key = (df["longitude"] / COORD_TOL_DEG).round().astype(int).to_numpy()
d2h = df["d2H_wax"].to_numpy()

buckets: dict[tuple[int, int], list[int]] = defaultdict(list)
for i in range(n):
    buckets[(lat_key[i], lon_key[i])].append(i)

seen_pairs = 0
for i in range(n):
    li, oi = lat_key[i], lon_key[i]
    for dla in (-1, 0, 1):
        for dlo in (-1, 0, 1):
            key = (li + dla, oi + dlo)
            if key not in buckets:
                continue
            for j in buckets[key]:
                if j <= i:
                    continue
                if abs(df.at[i, "latitude"] - df.at[j, "latitude"]) > COORD_TOL_DEG:
                    continue
                if abs(df.at[i, "longitude"] - df.at[j, "longitude"]) > COORD_TOL_DEG:
                    continue
                if abs(d2h[i] - d2h[j]) > D2H_TOL_PERMIL:
                    continue
                union(i, j)
                seen_pairs += 1

# Group by cluster root
clusters: dict[int, list[int]] = defaultdict(list)
for i in range(n):
    clusters[find(i)].append(i)

multi_clusters = {r: rows for r, rows in clusters.items() if len(rows) > 1}
print(f"Multi-row clusters (>=2 rows within tolerances): {len(multi_clusters)}")

# Apply the dedup rules
rows_to_remove: set[int] = set()
n_cross_study = 0
n_within_study = 0
replicates_kept_rows = 0
replicate_clusters = 0


def pick_cross_study_keeper(sub: pd.DataFrame) -> int:
    """Prefer a User compilation row; otherwise first row."""
    user_rows = sub[sub["data_source"] == "User compilation"]
    if len(user_rows) > 0:
        return int(user_rows.index[0])
    return int(sub.index[0])


def pick_within_study_keeper(sub: pd.DataFrame) -> int:
    """Prefer rows with a reported d2H_wax_err; otherwise first row."""
    with_err = sub[sub["d2H_wax_err"].notna()]
    if len(with_err) > 0:
        return int(with_err.index[0])
    return int(sub.index[0])


for rows in multi_clusters.values():
    sub = df.loc[rows]
    d2h_range = float(sub["d2H_wax"].max() - sub["d2H_wax"].min())
    sources = sub["source"].unique()

    if len(sources) > 1:
        keeper = pick_cross_study_keeper(sub)
        for idx in rows:
            if idx != keeper:
                rows_to_remove.add(idx)
                n_cross_study += 1
        continue

    if d2h_range > REPLICATE_THRESHOLD_PERMIL:
        replicate_clusters += 1
        replicates_kept_rows += len(rows)
        continue

    keeper = pick_within_study_keeper(sub)
    for idx in rows:
        if idx != keeper:
            rows_to_remove.add(idx)
            n_within_study += 1

df_clean = df.drop(index=sorted(rows_to_remove)).reset_index(drop=True)

print(f"Cross-study duplicates removed: {n_cross_study}")
print(f"Within-study duplicates removed: {n_within_study}")
print(
    f"Same-site replicate clusters retained: {replicate_clusters} "
    f"({replicates_kept_rows} rows)"
)

# ---------------------------------------------------------------------------
# Step 5: drop any remaining missing required fields (belt-and-suspenders)
# ---------------------------------------------------------------------------
final_missing = (
    df_clean["latitude"].isna()
    | df_clean["longitude"].isna()
    | df_clean["d2H_wax"].isna()
)
if final_missing.any():
    print(f"WARNING: dropping {int(final_missing.sum())} residual rows with missing fields")
    df_clean = df_clean.loc[~final_missing].reset_index(drop=True)

# ---------------------------------------------------------------------------
# Step 6: order columns and write out
# ---------------------------------------------------------------------------
ordered_cols = MODEL_COLUMNS + EXTRA_COLUMNS
for col in ordered_cols:
    if col not in df_clean.columns:
        df_clean[col] = np.nan
df_out = df_clean[ordered_cols].copy()
df_out.to_csv(OUTPUT_PATH, index=False)
print(f"\nWrote {OUTPUT_PATH} with {len(df_out)} rows x {len(df_out.columns)} cols")

# ---------------------------------------------------------------------------
# Step 7: summary report
# ---------------------------------------------------------------------------
print("\n" + "=" * 64)
print("SUMMARY")
print("=" * 64)
print(f"Starting rows                    : {n_start}")
print(f"After column mapping             : {n_after_map}")
print(f"Dropped for missing lat/lon/wax  : {n_dropped_missing}")
print(f"Cross-study duplicates removed   : {n_cross_study}")
print(f"Within-study duplicates removed  : {n_within_study}")
print(
    f"Same-site replicates kept        : {replicates_kept_rows} rows "
    f"in {replicate_clusters} clusters"
)
print(f"Final row count                  : {len(df_out)}")

print("\nBreakdown by data_source:")
for src, count in df_out["data_source"].value_counts(dropna=False).items():
    label = src if pd.notna(src) else "(missing)"
    print(f"  {label:<20} {count}")

print("\nGeographic coverage:")
regions = df_out.apply(
    lambda r: assign_region(r["latitude"], r["longitude"]), axis=1
)
for name in ["N America", "S America", "Europe", "Asia", "Africa",
             "Oceania", "Arctic", "Other", "Unknown"]:
    c = int((regions == name).sum())
    if c > 0:
        print(f"  {name:<12} {c}")

print("\nRanges:")
print(
    f"  latitude : {df_out['latitude'].min():.2f} to "
    f"{df_out['latitude'].max():.2f}"
)
print(
    f"  longitude: {df_out['longitude'].min():.2f} to "
    f"{df_out['longitude'].max():.2f}"
)
print(
    f"  d2H_wax  : {df_out['d2H_wax'].min():.1f} to "
    f"{df_out['d2H_wax'].max():.1f} permil"
)
print(f"  unique sources (publications): {df_out['source'].nunique()}")

print("\nDone.")
