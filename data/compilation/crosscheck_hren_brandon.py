"""Cross-check Hren & Brandon 2026 rows in merged dataset against Zenodo source.

Identifies:
- Value discrepancies (d2H_wax differs >3 permil)
- Missing rows (in Zenodo but not in merged)
- Potentially hallucinated rows (in merged but not in Zenodo)
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent
MERGED_PATH = BASE_DIR / "global_data_c29_merged.csv"
ZENODO_PATH = BASE_DIR / "hren_brandon_2026_table_s3.csv"

COORD_TOL = 0.05  # degrees
WAX_DIFF_THRESHOLD = 3.0  # permil


def load_merged_hren(path: Path) -> pd.DataFrame:
    """Load merged dataset and filter to Hren/Brandon rows."""
    df = pd.read_csv(path)
    mask = df["source"].str.contains("Hren|Brandon", case=False, na=False)
    hren = df[mask].copy().reset_index(drop=True)
    print(f"Merged dataset: {len(df)} total rows")
    print(f"Hren/Brandon rows in merged: {len(hren)}")
    return hren


def load_zenodo(path: Path) -> pd.DataFrame:
    """Load Zenodo source data for Hren & Brandon 2026."""
    df = pd.read_csv(path)
    # Ensure numeric types for coordinate and isotope columns
    for col in ["Latitude", "Longitude", "d2H C29", "d2H C27", "d2H C31"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    print(f"\nZenodo dataset: {len(df)} total rows")

    # Report flags
    flagged = df[df["Flag"].notna()]
    unflagged = df[df["Flag"].isna()]
    print(f"  Unflagged rows: {len(unflagged)}")
    print(f"  Flagged rows: {len(flagged)}")
    if len(flagged) > 0:
        flag_counts = flagged["Flag"].value_counts()
        for flag, count in flag_counts.items():
            print(f"    '{flag}': {count}")

    # Report rows with valid C29 data
    has_c29 = df["d2H C29"].notna()
    unflagged_c29 = df[df["Flag"].isna() & has_c29]
    print(f"  Unflagged rows with valid d2H C29: {len(unflagged_c29)}")

    return df


def match_datasets(
    merged_hren: pd.DataFrame, zenodo: pd.DataFrame, tol: float
) -> tuple[np.ndarray, np.ndarray, list[int], list[int]]:
    """Match merged Hren rows to Zenodo rows by coordinates.

    Returns:
        merged_matched: bool mask for merged rows that matched
        zenodo_matched: bool mask for zenodo rows that matched
        merged_to_zenodo: mapping from merged index to zenodo index (-1 if unmatched)
        zenodo_to_merged: mapping from zenodo index to merged index (-1 if unmatched)
    """
    m_lat = merged_hren["latitude"].values
    m_lon = merged_hren["longitude"].values
    z_lat = zenodo["Latitude"].values
    z_lon = zenodo["Longitude"].values

    merged_matched = np.zeros(len(merged_hren), dtype=bool)
    zenodo_matched = np.zeros(len(zenodo), dtype=bool)
    merged_to_zenodo = np.full(len(merged_hren), -1, dtype=int)
    zenodo_to_merged = np.full(len(zenodo), -1, dtype=int)

    for i in range(len(merged_hren)):
        if pd.isna(m_lat[i]) or pd.isna(m_lon[i]):
            continue
        lat_diff = np.abs(z_lat - m_lat[i])
        lon_diff = np.abs(z_lon - m_lon[i])
        within_tol = (lat_diff < tol) & (lon_diff < tol)
        if within_tol.any():
            dist = lat_diff + lon_diff
            dist[~within_tol] = np.inf
            best_idx = int(np.argmin(dist))
            merged_matched[i] = True
            zenodo_matched[best_idx] = True
            merged_to_zenodo[i] = best_idx
            zenodo_to_merged[best_idx] = i

    return merged_matched, zenodo_matched, merged_to_zenodo, zenodo_to_merged


def main() -> None:
    merged_hren = load_merged_hren(MERGED_PATH)
    zenodo = load_zenodo(ZENODO_PATH)

    if len(merged_hren) == 0:
        print("\nNo Hren/Brandon rows found in merged dataset. Nothing to cross-check.")
        return

    # Match
    print(f"\n{'=' * 60}")
    print("COORDINATE MATCHING")
    print(f"{'=' * 60}")
    merged_matched, zenodo_matched, m2z, z2m = match_datasets(
        merged_hren, zenodo, COORD_TOL
    )

    n_merged_matched = merged_matched.sum()
    n_zenodo_matched = zenodo_matched.sum()
    print(
        f"  Merged Hren rows matched to Zenodo: {n_merged_matched} / {len(merged_hren)}"
    )
    print(f"  Zenodo rows matched to merged: {n_zenodo_matched} / {len(zenodo)}")

    # --- d2H_wax value comparison ---
    print(f"\n{'=' * 60}")
    print(f"VALUE COMPARISON (threshold: >{WAX_DIFF_THRESHOLD} permil)")
    print(f"{'=' * 60}")

    value_discrepancies: list[dict] = []
    for i in range(len(merged_hren)):
        if not merged_matched[i]:
            continue
        z_idx = m2z[i]
        merged_wax = merged_hren["d2H_wax"].iloc[i]
        zenodo_wax = zenodo["d2H C29"].iloc[z_idx]

        if pd.isna(merged_wax) or pd.isna(zenodo_wax):
            continue

        diff = abs(merged_wax - zenodo_wax)
        if diff > WAX_DIFF_THRESHOLD:
            value_discrepancies.append(
                {
                    "merged_source": merged_hren["source"].iloc[i],
                    "merged_location": merged_hren["location"].iloc[i],
                    "zenodo_sample": zenodo["Sample"].iloc[z_idx],
                    "lat": merged_hren["latitude"].iloc[i],
                    "lon": merged_hren["longitude"].iloc[i],
                    "merged_d2H_wax": merged_wax,
                    "zenodo_d2H_C29": zenodo_wax,
                    "diff": diff,
                    "zenodo_flag": zenodo["Flag"].iloc[z_idx],
                }
            )

    if value_discrepancies:
        print(f"  Found {len(value_discrepancies)} discrepancies:")
        for d in value_discrepancies:
            flag_str = (
                f" [FLAG: {d['zenodo_flag']}]" if pd.notna(d["zenodo_flag"]) else ""
            )
            print(
                f"    {d['zenodo_sample']} ({d['lat']:.4f}, {d['lon']:.4f}): "
                f"merged={d['merged_d2H_wax']:.1f}, zenodo={d['zenodo_d2H_C29']:.1f}, "
                f"diff={d['diff']:.1f}{flag_str}"
            )
    else:
        print(
            "  No value discrepancies found. All matched values agree within 3 permil."
        )

    # --- Zenodo rows NOT in merged (missing) ---
    print(f"\n{'=' * 60}")
    print("ZENODO ROWS NOT IN MERGED DATASET")
    print(f"{'=' * 60}")

    zenodo_missing = zenodo[~zenodo_matched].copy()
    zenodo_missing_unflagged = zenodo_missing[zenodo_missing["Flag"].isna()]
    zenodo_missing_flagged = zenodo_missing[zenodo_missing["Flag"].notna()]
    zenodo_missing_has_c29 = zenodo_missing_unflagged[
        zenodo_missing_unflagged["d2H C29"].notna()
    ]

    print(f"  Total Zenodo rows not matched: {len(zenodo_missing)}")
    print(f"    Unflagged: {len(zenodo_missing_unflagged)}")
    print(f"    Flagged: {len(zenodo_missing_flagged)}")
    print(f"    Unflagged with valid C29: {len(zenodo_missing_has_c29)}")

    if len(zenodo_missing_has_c29) > 0:
        print("\n  Unflagged Zenodo rows with C29 data missing from merged:")
        for _, row in zenodo_missing_has_c29.iterrows():
            print(
                f"    {row['Sample']} ({row['Latitude']:.4f}, {row['Longitude']:.4f}): "
                f"d2H_C29={row['d2H C29']:.1f}, "
                f"Source: {str(row['Source'])[:60]}"
            )

    # --- Merged Hren rows NOT in Zenodo (potentially hallucinated) ---
    print(f"\n{'=' * 60}")
    print("MERGED HREN/BRANDON ROWS NOT IN ZENODO (potentially hallucinated)")
    print(f"{'=' * 60}")

    merged_unmatched = merged_hren[~merged_matched].copy()
    print(f"  Total unmatched: {len(merged_unmatched)}")
    if len(merged_unmatched) > 0:
        for _, row in merged_unmatched.iterrows():
            lat_str = f"{row['latitude']:.4f}" if pd.notna(row["latitude"]) else "NaN"
            lon_str = f"{row['longitude']:.4f}" if pd.notna(row["longitude"]) else "NaN"
            wax_str = f"{row['d2H_wax']:.1f}" if pd.notna(row["d2H_wax"]) else "NaN"
            ds = row.get("data_source", "unknown")
            print(
                f"    {row['source']}, {row['location']} "
                f"({lat_str}, {lon_str}): "
                f"d2H_wax={wax_str}, data_source={ds}"
            )

    # --- Summary ---
    print(f"\n{'=' * 60}")
    print("CROSS-CHECK SUMMARY")
    print(f"{'=' * 60}")
    print(f"  Merged Hren/Brandon rows: {len(merged_hren)}")
    print(f"  Zenodo rows: {len(zenodo)} ({(zenodo['Flag'].isna()).sum()} unflagged)")
    print(f"  Coordinate matches: {n_merged_matched}")
    print(f"  Value discrepancies (>{WAX_DIFF_THRESHOLD}‰): {len(value_discrepancies)}")
    print(
        f"  Zenodo unflagged+C29 rows missing from merged: {len(zenodo_missing_has_c29)}"
    )
    print(f"  Merged rows not in Zenodo (hallucination risk): {len(merged_unmatched)}")

    # Save discrepancy report
    report_path = BASE_DIR / "hren_brandon_crosscheck_report.csv"
    if value_discrepancies:
        pd.DataFrame(value_discrepancies).to_csv(report_path, index=False)
        print(f"\n  Value discrepancy report saved to: {report_path}")


if __name__ == "__main__":
    main()
