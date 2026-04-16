"""Merge original C29 leaf wax dataset with v7 Perplexity-generated dataset.

Keeps original rows as authoritative; adds new v7 rows that have
valid coordinates and d2H_wax but no coordinate match in the original.
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent
ORIGINAL_PATH = BASE_DIR / "global_data_c29_ORIGINAL.csv"
V7_PATH = BASE_DIR / "modern_leaf_wax_d2H_C29_merged_v7.csv"
OUTPUT_PATH = BASE_DIR / "global_data_c29_merged.csv"

COORD_TOL = 0.05  # degrees
WAX_DIFF_THRESHOLD = 5.0  # permil

# ---------------------------------------------------------------------------
# Column mapping: v7 -> original format
# ---------------------------------------------------------------------------
V7_COL_MAP: dict[str, str] = {
    "Reference": "source",
    "Compilation": "compilation",
    "Location": "location",
    "Latitude (°N)": "latitude",
    "Longitude (°E)": "longitude",
    "Elevation (m)": "elevation",
    "δ²H n-C29 (‰ VSMOW)": "d2H_wax",
    "δ²H n-C29 SD (±‰)": "d2H_wax_err",
    "δ²H Precipitation (‰ VSMOW)": "d2H_precip",
    "δ²H Precipitation SD (±‰)": "d2H_precip_err",
}

# Extra v7 columns to preserve for new rows
V7_EXTRA_COLS: list[str] = [
    "Sample Type",
    "Functional Type",
    "MAP (mm/yr)",
    "MAT (°C)",
    "DOI",
    "Audit Flags",
]


def load_original(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    print(f"Original dataset: {len(df)} rows")
    return df


def load_v7(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    print(f"V7 dataset: {len(df)} rows")
    return df


def map_v7_columns(v7: pd.DataFrame) -> pd.DataFrame:
    """Rename v7 columns to match original schema and add chain=29."""
    mapped = v7.rename(columns=V7_COL_MAP).copy()
    mapped["chain"] = 29
    # Carry through extra columns under original names
    for col in V7_EXTRA_COLS:
        if col in v7.columns:
            mapped[col] = v7[col]
    return mapped


def find_coordinate_matches(
    orig: pd.DataFrame, v7: pd.DataFrame, tol: float
) -> tuple[np.ndarray, np.ndarray, list[dict]]:
    """For each v7 row, find if any original row matches within tol.

    Returns:
        matched_v7_mask: boolean array, True if v7 row has a match in original
        matched_orig_indices: for each matched v7 row, the index of best orig match
        discrepancies: list of dicts for d2H_wax differences > WAX_DIFF_THRESHOLD
    """
    orig_lat = orig["latitude"].values
    orig_lon = orig["longitude"].values
    v7_lat = v7["latitude"].values
    v7_lon = v7["longitude"].values

    matched_v7_mask = np.zeros(len(v7), dtype=bool)
    matched_orig_indices = np.full(len(v7), -1, dtype=int)
    discrepancies: list[dict] = []

    for i in range(len(v7)):
        if pd.isna(v7_lat[i]) or pd.isna(v7_lon[i]):
            continue
        lat_diff = np.abs(orig_lat - v7_lat[i])
        lon_diff = np.abs(orig_lon - v7_lon[i])
        within_tol = (lat_diff < tol) & (lon_diff < tol)
        if within_tol.any():
            # Pick closest match
            dist = lat_diff + lon_diff
            dist[~within_tol] = np.inf
            best_idx = int(np.argmin(dist))
            matched_v7_mask[i] = True
            matched_orig_indices[i] = best_idx

            # Check d2H_wax discrepancy
            v7_wax = v7["d2H_wax"].iloc[i]
            orig_wax = orig["d2H_wax"].iloc[best_idx]
            if pd.notna(v7_wax) and pd.notna(orig_wax):
                diff = abs(v7_wax - orig_wax)
                if diff > WAX_DIFF_THRESHOLD:
                    discrepancies.append(
                        {
                            "v7_row": i,
                            "orig_row": best_idx,
                            "source": v7["source"].iloc[i],
                            "location": v7["location"].iloc[i],
                            "lat": v7_lat[i],
                            "lon": v7_lon[i],
                            "v7_d2H_wax": v7_wax,
                            "orig_d2H_wax": orig_wax,
                            "diff": diff,
                        }
                    )

    return matched_v7_mask, matched_orig_indices, discrepancies


def check_schwab_cameroon(orig: pd.DataFrame, v7: pd.DataFrame) -> None:
    """Report on Schwab et al. 2015 Cameroon sites."""
    orig_schwab = orig[orig["source"].str.contains("Schwab", case=False, na=False)]
    v7_schwab = v7[v7["source"].str.contains("Schwab", case=False, na=False)]

    print("\n--- Schwab et al. 2015 check ---")
    print(f"  Original: {len(orig_schwab)} rows")
    print(f"  V7: {len(v7_schwab)} rows")

    if len(v7_schwab) > 0:
        missing_lon = v7_schwab[v7_schwab["longitude"].isna()]
        print(f"  V7 rows missing longitude: {len(missing_lon)}")
    print(
        f"  All {len(orig_schwab)} original Schwab rows will be kept as authoritative."
    )


def classify_latitude_band(lat: float) -> str:
    """Classify latitude into geographic region."""
    if pd.isna(lat):
        return "Unknown"
    if lat > 60:
        return "Arctic/Subarctic (>60N)"
    if lat > 30:
        return "Northern mid-latitudes (30-60N)"
    if lat > 0:
        return "Northern tropics (0-30N)"
    if lat > -30:
        return "Southern tropics (0-30S)"
    if lat > -60:
        return "Southern mid-latitudes (30-60S)"
    return "Antarctic (<60S)"


def main() -> None:
    # Load data
    orig = load_original(ORIGINAL_PATH)
    v7_raw = load_v7(V7_PATH)

    # Map v7 columns
    v7 = map_v7_columns(v7_raw)

    # Check Schwab sites
    check_schwab_cameroon(orig, v7)

    # Find coordinate matches
    print("\nMatching v7 rows to original by coordinates...")
    matched_mask, matched_indices, discrepancies = find_coordinate_matches(
        orig, v7, COORD_TOL
    )

    n_matched = matched_mask.sum()
    print(f"  V7 rows matched to original: {n_matched}")
    print(f"  V7 rows NOT matched (candidate new): {(~matched_mask).sum()}")

    # Report discrepancies
    if discrepancies:
        print(f"\n--- d2H_wax discrepancies > {WAX_DIFF_THRESHOLD} permil ---")
        for d in discrepancies:
            print(
                f"  {d['source']}, {d['location']} "
                f"({d['lat']:.3f}, {d['lon']:.3f}): "
                f"v7={d['v7_d2H_wax']:.1f}, orig={d['orig_d2H_wax']:.1f}, "
                f"diff={d['diff']:.1f}"
            )
    else:
        print("\n  No d2H_wax discrepancies > 5 permil found.")

    # Identify new rows from v7
    unmatched_v7 = v7[~matched_mask].copy()

    # Filter: must have valid lat, lon, d2H_wax
    valid_new = unmatched_v7[
        unmatched_v7["latitude"].notna()
        & unmatched_v7["longitude"].notna()
        & unmatched_v7["d2H_wax"].notna()
    ].copy()

    dropped = len(unmatched_v7) - len(valid_new)
    print(f"\n  Unmatched v7 rows dropped (missing lat/lon/d2H_wax): {dropped}")
    print(f"  Valid new rows to add: {len(valid_new)}")

    # Build merged dataset
    # Original rows with provenance
    orig_merged = orig.copy()
    orig_merged["data_source"] = "original"

    # New v7 rows with provenance
    # Select columns that match original schema + extras
    orig_cols = list(orig.columns)
    new_rows = valid_new[[c for c in orig_cols if c in valid_new.columns]].copy()
    new_rows["data_source"] = "v7_new"

    # Add extra v7 columns
    for col in V7_EXTRA_COLS:
        if col in valid_new.columns:
            new_rows[col] = valid_new[col].values

    # Also add extras as NaN to original rows for consistent schema
    for col in V7_EXTRA_COLS:
        if col not in orig_merged.columns:
            orig_merged[col] = np.nan

    # Concatenate
    merged = pd.concat([orig_merged, new_rows], ignore_index=True)

    print(f"\n{'=' * 60}")
    print("MERGE SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Total rows in merged dataset: {len(merged)}")
    print(f"  Rows from original: {(merged['data_source'] == 'original').sum()}")
    print(f"  New rows from v7:   {(merged['data_source'] == 'v7_new').sum()}")

    # Geographic breakdown of new rows
    new_only = merged[merged["data_source"] == "v7_new"]
    if len(new_only) > 0:
        new_only = new_only.copy()
        new_only["region"] = new_only["latitude"].apply(classify_latitude_band)
        print("\n  Geographic breakdown of new rows:")
        for region, count in new_only["region"].value_counts().sort_index().items():
            print(f"    {region}: {count}")

        print("\n  Top sources in new rows:")
        for source, count in new_only["source"].value_counts().head(15).items():
            print(f"    {source}: {count}")

    # Save
    merged.to_csv(OUTPUT_PATH, index=False)
    print(f"\nSaved merged dataset to: {OUTPUT_PATH}")
    print(f"  Columns: {list(merged.columns)}")


if __name__ == "__main__":
    main()
