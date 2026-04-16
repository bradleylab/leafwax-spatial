"""
Cross-check Wang et al. 2017 rows across datasets and detect duplicates in merged dataset.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from itertools import combinations

BASE = Path("/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/_for_GCA")

# ============================================================
# TASK 1: Wang et al. 2017 cross-check
# ============================================================
print("=" * 80)
print("TASK 1: Wang et al. 2017 cross-check")
print("=" * 80)

# Load datasets
original = pd.read_csv(BASE / "global_data_c29_ORIGINAL.csv")
zenodo = pd.read_csv(BASE / "hren_brandon_2026_table_s3.csv")
# Convert Zenodo lat/lon to numeric (some rows have "-" or other non-numeric)
zenodo["Latitude"] = pd.to_numeric(zenodo["Latitude"], errors="coerce")
zenodo["Longitude"] = pd.to_numeric(zenodo["Longitude"], errors="coerce")
zenodo["d2H C29"] = pd.to_numeric(zenodo["d2H C29"], errors="coerce")
merged = pd.read_csv(BASE / "global_data_c29_merged.csv")

print(f"\nOriginal dataset: {len(original)} rows, columns: {list(original.columns)}")
print(f"Zenodo dataset: {len(zenodo)} rows, columns: {list(zenodo.columns)}")
print(f"Merged dataset: {len(merged)} rows, columns: {list(merged.columns)}")

# Extract Wang rows from each dataset
wang_orig = original[original["source"].str.contains("Wang", case=False, na=False)].copy()
wang_zenodo = zenodo[zenodo["Source"].str.contains("Wang", case=False, na=False)].copy()
wang_merged = merged[merged["source"].str.contains("Wang", case=False, na=False)].copy()

print(f"\n--- Wang row counts ---")
print(f"Original: {len(wang_orig)} Wang rows")
print(f"Zenodo:   {len(wang_zenodo)} Wang rows")
print(f"Merged:   {len(wang_merged)} Wang rows")

# Show Wang sources
print(f"\nWang sources in original: {wang_orig['source'].unique()}")
print(f"Wang sources in zenodo: {wang_zenodo['Source'].unique()}")
print(f"Wang sources in merged: {wang_merged['source'].unique()}")

# Zenodo columns: Longitude (idx 2), Latitude (idx 3), d2H C29 (idx 5)
print("\n--- Matching Zenodo Wang rows to Original Wang rows (within 0.1 deg) ---")

MATCH_TOL = 0.1

matched_zenodo_to_orig = []
unmatched_zenodo = []

for zi, zrow in wang_zenodo.iterrows():
    z_lat = zrow["Latitude"]
    z_lon = zrow["Longitude"]
    z_d2h = zrow["d2H C29"]

    # Find matches in original
    matches = []
    for oi, orow in wang_orig.iterrows():
        o_lat = orow["latitude"]
        o_lon = orow["longitude"]
        o_d2h = orow["d2H_wax"]

        if abs(z_lat - o_lat) < MATCH_TOL and abs(z_lon - o_lon) < MATCH_TOL:
            matches.append((oi, orow, abs(z_d2h - o_d2h) if pd.notna(z_d2h) and pd.notna(o_d2h) else None))

    if matches:
        for oi, orow, d2h_diff in matches:
            matched_zenodo_to_orig.append({
                "zenodo_idx": zi,
                "orig_idx": oi,
                "zenodo_sample": zrow.get("Sample", ""),
                "orig_location": orow.get("location", ""),
                "zenodo_lat": z_lat,
                "zenodo_lon": z_lon,
                "orig_lat": orow["latitude"],
                "orig_lon": orow["longitude"],
                "zenodo_d2H": z_d2h,
                "orig_d2H": orow["d2H_wax"],
                "d2H_diff": d2h_diff,
            })
    else:
        unmatched_zenodo.append({
            "zenodo_idx": zi,
            "sample": zrow.get("Sample", ""),
            "lat": z_lat,
            "lon": z_lon,
            "d2H_C29": z_d2h,
            "source": zrow["Source"],
        })

print(f"\nZenodo Wang rows matching original Wang rows: {len(matched_zenodo_to_orig)}")
if matched_zenodo_to_orig:
    match_df = pd.DataFrame(matched_zenodo_to_orig)
    print(match_df.to_string(index=False))

    # Check d2H agreement
    discrepancies = match_df[match_df["d2H_diff"] > 2]
    if len(discrepancies) > 0:
        print(f"\n  ** {len(discrepancies)} pairs with d2H discrepancy > 2 permil:")
        print(discrepancies.to_string(index=False))
    else:
        print("\n  All matched pairs agree within 2 permil on d2H.")

print(f"\nZenodo Wang rows with NO match in original: {len(unmatched_zenodo)}")
if unmatched_zenodo:
    unmatched_df = pd.DataFrame(unmatched_zenodo)
    print(unmatched_df.to_string(index=False))

# Now check unmatched Zenodo Wang rows against merged dataset
print("\n--- Checking unmatched Zenodo Wang rows against merged dataset ---")

truly_new = []
for row in unmatched_zenodo:
    z_lat = row["lat"]
    z_lon = row["lon"]
    z_d2h = row["d2H_C29"]

    found_in_merged = False
    for mi, mrow in merged.iterrows():
        m_lat = mrow["latitude"]
        m_lon = mrow["longitude"]
        if abs(z_lat - m_lat) < MATCH_TOL and abs(z_lon - m_lon) < MATCH_TOL:
            found_in_merged = True
            m_d2h = mrow["d2H_wax"]
            d2h_diff = abs(z_d2h - m_d2h) if pd.notna(z_d2h) and pd.notna(m_d2h) else None
            print(f"  Zenodo row {row['zenodo_idx']} ({row['sample']}) at ({z_lat}, {z_lon}) "
                  f"matches merged row {mi} ({mrow['source']}) — "
                  f"d2H: zenodo={z_d2h}, merged={m_d2h}, diff={d2h_diff}")
            break

    if not found_in_merged:
        truly_new.append(row)

print(f"\nTruly new Zenodo Wang rows (not in original or merged): {len(truly_new)}")
if truly_new:
    print(pd.DataFrame(truly_new).to_string(index=False))

# ============================================================
# TASK 2: Full duplicate detection on merged dataset
# ============================================================
print("\n" + "=" * 80)
print("TASK 2: Full duplicate detection on merged dataset")
print("=" * 80)

print(f"\nMerged dataset: {len(merged)} rows")

# Columns of interest
REPORT_COLS = ["source", "location", "latitude", "longitude", "d2H_wax", "d2H_wax_err", "data_source"]

# Build duplicate pairs
dup_pairs = []
n = len(merged)

lats = merged["latitude"].values
lons = merged["longitude"].values
d2h = merged["d2H_wax"].values

print("Scanning for duplicate pairs...")

for i in range(n):
    for j in range(i + 1, n):
        lat_diff = abs(lats[i] - lats[j])
        lon_diff = abs(lons[i] - lons[j])

        if pd.isna(lats[i]) or pd.isna(lats[j]) or pd.isna(lons[i]) or pd.isna(lons[j]):
            continue

        d2h_i = d2h[i]
        d2h_j = d2h[j]
        d2h_diff = abs(d2h_i - d2h_j) if pd.notna(d2h_i) and pd.notna(d2h_j) else None

        is_dup = False

        # Criterion 1: lat/lon within 0.01 AND d2H within 2
        if lat_diff < 0.01 and lon_diff < 0.01:
            if d2h_diff is not None and d2h_diff <= 2:
                is_dup = True

        # Criterion 2: lat/lon within 0.05 AND d2H exact match
        if lat_diff < 0.05 and lon_diff < 0.05:
            if d2h_diff is not None and d2h_diff == 0:
                is_dup = True

        if is_dup:
            dup_pairs.append((i, j))

print(f"Found {len(dup_pairs)} duplicate pairs")

# Also find rows with same location but different d2H (>5 permil)
print("\n--- Rows with similar lat/lon but DIFFERENT d2H (>5 permil) ---")
conflict_pairs = []
for i in range(n):
    for j in range(i + 1, n):
        if pd.isna(lats[i]) or pd.isna(lats[j]) or pd.isna(lons[i]) or pd.isna(lons[j]):
            continue
        lat_diff = abs(lats[i] - lats[j])
        lon_diff = abs(lons[i] - lons[j])

        if lat_diff < 0.01 and lon_diff < 0.01:
            d2h_i = d2h[i]
            d2h_j = d2h[j]
            if pd.notna(d2h_i) and pd.notna(d2h_j):
                d2h_diff = abs(d2h_i - d2h_j)
                if d2h_diff > 5:
                    conflict_pairs.append((i, j, d2h_diff))

print(f"Found {len(conflict_pairs)} conflict pairs (same location, d2H diff > 5 permil)")
for i, j, diff in conflict_pairs:
    row_i = merged.iloc[i]
    row_j = merged.iloc[j]
    print(f"  Rows {i} & {j}: d2H_diff={diff:.1f}")
    print(f"    [{i}] {row_i['source']}, {row_i.get('location','')}, "
          f"({row_i['latitude']:.4f}, {row_i['longitude']:.4f}), "
          f"d2H={row_i['d2H_wax']}, err={row_i['d2H_wax_err']}, ds={row_i['data_source']}")
    print(f"    [{j}] {row_j['source']}, {row_j.get('location','')}, "
          f"({row_j['latitude']:.4f}, {row_j['longitude']:.4f}), "
          f"d2H={row_j['d2H_wax']}, err={row_j['d2H_wax_err']}, ds={row_j['data_source']}")

# Cluster duplicates using union-find
parent = list(range(n))

def find(x: int) -> int:
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

def union(x: int, y: int) -> None:
    px, py = find(x), find(y)
    if px != py:
        parent[px] = py

for i, j in dup_pairs:
    union(i, j)

# Group into clusters
from collections import defaultdict
clusters: dict[int, list[int]] = defaultdict(list)
for i, j in dup_pairs:
    # Only include rows that are part of duplicate pairs
    pass

# Rebuild: collect all rows involved in any dup pair
dup_rows = set()
for i, j in dup_pairs:
    dup_rows.add(i)
    dup_rows.add(j)

clusters = defaultdict(list)
for r in sorted(dup_rows):
    root = find(r)
    clusters[root].append(r)

print(f"\n--- Duplicate clusters ---")
print(f"Total duplicate clusters: {len(clusters)}")
total_remove = sum(len(members) - 1 for members in clusters.values())
print(f"Total rows that would be removed: {total_remove}")

# For each cluster, recommend which row to keep
flag_records = []

for cluster_id, (root, members) in enumerate(sorted(clusters.items(), key=lambda x: x[1][0])):
    print(f"\n  Cluster {cluster_id + 1} ({len(members)} rows):")

    # Score each row: prefer "original" > "v7_new", prefer has d2H_wax_err, prefer more metadata
    best_idx = members[0]
    best_score = -999

    for idx in members:
        row = merged.iloc[idx]
        score = 0
        # Prefer original
        if row.get("data_source") == "original":
            score += 10
        elif str(row.get("data_source", "")).startswith("v7"):
            score += 5
        # Prefer has d2H_wax_err
        if pd.notna(row.get("d2H_wax_err")):
            score += 3
        # Prefer more non-null columns
        non_null_count = row.notna().sum()
        score += non_null_count * 0.1

        if score > best_score:
            best_score = score
            best_idx = idx

    for idx in members:
        row = merged.iloc[idx]
        keep = "KEEP" if idx == best_idx else "REMOVE"
        print(f"    [{idx}] {keep} | {row['source']}, {row.get('location','')}, "
              f"({row['latitude']:.4f}, {row['longitude']:.4f}), "
              f"d2H={row['d2H_wax']}, err={row['d2H_wax_err']}, ds={row['data_source']}")

        flag_records.append({
            "row_index": idx,
            "cluster_id": cluster_id + 1,
            "action": keep,
            "source": row["source"],
            "location": row.get("location", ""),
            "latitude": row["latitude"],
            "longitude": row["longitude"],
            "d2H_wax": row["d2H_wax"],
            "d2H_wax_err": row["d2H_wax_err"],
            "data_source": row["data_source"],
        })

# Save flagged duplicates
flag_df = pd.DataFrame(flag_records)
out_path = BASE / "duplicate_flags.csv"
flag_df.to_csv(out_path, index=False)
print(f"\n--- Summary ---")
print(f"Total duplicate clusters: {len(clusters)}")
print(f"Total rows flagged: {len(flag_records)}")
print(f"Total rows to remove: {total_remove}")
print(f"Saved to: {out_path}")
