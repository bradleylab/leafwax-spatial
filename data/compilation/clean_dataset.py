"""
Data cleaning for global C29 leaf wax dataset.

Tasks:
1. Fix Wang et al. 2017 transcription error (d2H_wax -288 -> -188)
2. Remove cross-study and within-study true duplicates using duplicate_flags.csv
3. Save cleaned dataset and produce summary statistics
"""

from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path("/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/_for_GCA")
MERGED_PATH = BASE / "global_data_c29_merged.csv"
FLAGS_PATH = BASE / "duplicate_flags.csv"
OUTPUT_PATH = BASE / "global_data_c29_cleaned.csv"

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
df = pd.read_csv(MERGED_PATH)
flags = pd.read_csv(FLAGS_PATH)

n_start = len(df)
print(f"Starting row count: {n_start}")

# ---------------------------------------------------------------------------
# Task 1: Fix Wang et al. 2017 typo
# ---------------------------------------------------------------------------
typo_mask = (
    df["source"].str.contains("Wang et al., 2017", na=False)
    & (df["d2H_wax"] == -288.0)
    & (df["latitude"].between(22.6, 22.8))
    & (df["longitude"].between(101.2, 101.4))
)

n_typo = typo_mask.sum()
assert n_typo == 1, f"Expected 1 Wang typo row, found {n_typo}"

old_val = df.loc[typo_mask, "d2H_wax"].values[0]
df.loc[typo_mask, "d2H_wax"] = -188.0
print(f"\nTask 1: Fixed Wang et al. 2017 typo: {old_val} -> -188.0")
print(f"  Row index: {df.index[typo_mask][0]}")

# ---------------------------------------------------------------------------
# Task 2: Remove duplicates using duplicate_flags.csv
# ---------------------------------------------------------------------------
# Classify each cluster
COORD_TOL = 0.01
D2H_TOL = 2.0
D2H_REPLICATE_THRESHOLD = 5.0

removal_log: list[dict] = []
rows_to_remove: list[int] = []
replicates_retained: list[dict] = []

for cluster_id in sorted(flags["cluster_id"].unique()):
    cluster = flags[flags["cluster_id"] == cluster_id].copy()
    sources_in_cluster = cluster["source"].unique()
    d2h_range = cluster["d2H_wax"].max() - cluster["d2H_wax"].min()

    # Rule 1: Cross-study duplicates (different sources) -- remove flagged rows
    if len(sources_in_cluster) > 1:
        keep_rows = cluster[cluster["action"] == "KEEP"]
        remove_rows = cluster[cluster["action"] == "REMOVE"]
        for _, row in remove_rows.iterrows():
            kept_row = keep_rows.iloc[0]
            removal_log.append({
                "cluster_id": cluster_id,
                "removed_row_index": row["row_index"],
                "removed_source": row["source"],
                "removed_location": row["location"],
                "removed_lat": row["latitude"],
                "removed_lon": row["longitude"],
                "removed_d2H_wax": row["d2H_wax"],
                "reason": "cross-study duplicate",
                "kept_row_index": kept_row["row_index"],
                "kept_source": kept_row["source"],
            })
            rows_to_remove.append(row["row_index"])
        continue

    # Same source from here on
    # Rule 2: Same source, d2H range > threshold -- same-site replicates, KEEP ALL
    if d2h_range > D2H_REPLICATE_THRESHOLD:
        for _, row in cluster.iterrows():
            replicates_retained.append({
                "cluster_id": cluster_id,
                "row_index": row["row_index"],
                "source": row["source"],
                "location": row["location"],
                "lat": row["latitude"],
                "lon": row["longitude"],
                "d2H_wax": row["d2H_wax"],
                "d2h_range_in_cluster": d2h_range,
            })
        continue

    # Rule 3: Same source, d2H within threshold -- true duplicates, remove flagged
    keep_rows = cluster[cluster["action"] == "KEEP"]
    remove_rows = cluster[cluster["action"] == "REMOVE"]
    for _, row in remove_rows.iterrows():
        kept_row = keep_rows.iloc[0]
        removal_log.append({
            "cluster_id": cluster_id,
            "removed_row_index": row["row_index"],
            "removed_source": row["source"],
            "removed_location": row["location"],
            "removed_lat": row["latitude"],
            "removed_lon": row["longitude"],
            "removed_d2H_wax": row["d2H_wax"],
            "reason": "within-study duplicate (same source, d2H within 2 permil)",
            "kept_row_index": kept_row["row_index"],
            "kept_source": kept_row["source"],
        })
        rows_to_remove.append(row["row_index"])

# Apply removals
rows_to_remove_set = set(rows_to_remove)
df_cleaned = df.drop(index=list(rows_to_remove_set)).reset_index(drop=True)

n_removed = n_start - len(df_cleaned)
n_cross_study = sum(1 for r in removal_log if r["reason"] == "cross-study duplicate")
n_within_study = sum(
    1 for r in removal_log if r["reason"].startswith("within-study")
)

print(f"\nTask 2: Deduplication summary")
print(f"  Cross-study duplicates removed: {n_cross_study}")
print(f"  Within-study duplicates removed: {n_within_study}")
print(f"  Total rows removed: {n_removed}")
print(f"  Same-site replicates retained (clusters): "
      f"{len(set(r['cluster_id'] for r in replicates_retained))}")
print(f"  Same-site replicate rows retained: {len(replicates_retained)}")

# ---------------------------------------------------------------------------
# Task 3: Save cleaned dataset
# ---------------------------------------------------------------------------
df_cleaned.to_csv(OUTPUT_PATH, index=False)
print(f"\nTask 3: Saved cleaned dataset to {OUTPUT_PATH}")
print(f"  Final row count: {len(df_cleaned)}")

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("SUMMARY STATISTICS")
print("=" * 60)
print(f"  Rows before cleaning: {n_start}")
print(f"  Typo corrections: 1")
print(f"  Rows removed (deduplication): {n_removed}")
print(f"  Rows after cleaning: {len(df_cleaned)}")

print(f"\n  Source breakdown (data_source):")
for src, count in df_cleaned["data_source"].value_counts().items():
    print(f"    {src}: {count}")

print(f"\n  Unique sources (publications): {df_cleaned['source'].nunique()}")
print(f"  Latitude range: {df_cleaned['latitude'].min():.2f} to "
      f"{df_cleaned['latitude'].max():.2f}")
print(f"  Longitude range: {df_cleaned['longitude'].min():.2f} to "
      f"{df_cleaned['longitude'].max():.2f}")
print(f"  d2H_wax range: {df_cleaned['d2H_wax'].min():.1f} to "
      f"{df_cleaned['d2H_wax'].max():.1f}")

# ---------------------------------------------------------------------------
# Build removal log for documentation
# ---------------------------------------------------------------------------
removal_df = pd.DataFrame(removal_log)
print(f"\n  Removal log ({len(removal_df)} entries):")
for _, r in removal_df.iterrows():
    print(f"    Cluster {r['cluster_id']}: removed row {r['removed_row_index']} "
          f"({r['removed_source']}, {r['removed_location']}, "
          f"d2H={r['removed_d2H_wax']}) -- {r['reason']} "
          f"-- kept row {r['kept_row_index']} ({r['kept_source']})")

# Save removal log for documentation script reference
removal_df.to_csv(BASE / "_removal_log.csv", index=False)

replicates_df = pd.DataFrame(replicates_retained)
replicates_df.to_csv(BASE / "_replicates_retained.csv", index=False)

print("\nDone.")
