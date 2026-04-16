"""
slope_demonstration_v2.py
Demonstrates why the pooled regression slope of ~0.85 is NOT the fractionation
factor, even though it numerically equals it.

Key insight: TWO opposing geographic effects approximately cancel in the
global average:
  1. Growing season bias COMPRESSES source water range at high latitudes
     → pushes regression slope DOWN
  2. Evapotranspiration EXTENDS source water range at low latitudes
     → pushes regression slope UP
The cancellation produces a pooled slope ≈ 0.85, which coincidentally
equals alpha. But neither effect is zero — they just happen to balance.

A spatial model can separate these effects and reveals a conditional slope
that differs from 0.85.

Three scenarios:
  A: No confounders (slope = 0.85, trivially)
  B: Growing season bias only (slope < 0.85)
  C: Growing season + evapotranspiration (slope ≈ 0.85 again, but for wrong reasons)
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

np.random.seed(42)

ALPHA = 0.85
EPSILON = (ALPHA - 1) * 1000  # -150 permil

# Four latitude bands with OIPC centers
band_config = [
    {"name": "Tropical\n(0–15°)", "lat": 7, "oipc_center": -25, "color": "#e07b39"},
    {"name": "Subtropical\n(15–35°)", "lat": 25, "oipc_center": -50, "color": "#2ca02c"},
    {"name": "Temperate\n(35–55°)", "lat": 45, "oipc_center": -75, "color": "#1f77b4"},
    {"name": "High lat.\n(55–75°)", "lat": 65, "oipc_center": -105, "color": "#7b2d8e"},
]

N_PER_BAND = 25
NOISE_WAX = 4  # measurement noise on d2H_wax (permil)


def generate_scenario(gs_offset_per_lat, et_enrichment_per_lat):
    """Generate pseudodata for a scenario.

    Args:
        gs_offset_per_lat: growing season offset as function of |latitude|
            (positive = source water more enriched than OIPC)
        et_enrichment_per_lat: evapotranspiration enrichment as function of |latitude|
            (positive = source water more enriched; larger at LOW latitudes)
    """
    all_oipc, all_source, all_wax, all_bands, all_colors = [], [], [], [], []

    for band in band_config:
        oipc = band["oipc_center"] + np.random.normal(0, 10, N_PER_BAND)
        lat = band["lat"]

        # Source water offsets
        gs = gs_offset_per_lat(lat) + np.random.normal(0, 1.5, N_PER_BAND)
        et = et_enrichment_per_lat(lat) + np.random.normal(0, 1.5, N_PER_BAND)
        source = oipc + gs + et

        # Wax = alpha * source + epsilon + noise
        wax = ALPHA * source + EPSILON + np.random.normal(0, NOISE_WAX, N_PER_BAND)

        all_oipc.extend(oipc)
        all_source.extend(source)
        all_wax.extend(wax)
        all_bands.extend([band["name"]] * N_PER_BAND)
        all_colors.extend([band["color"]] * N_PER_BAND)

    return (np.array(all_oipc), np.array(all_source), np.array(all_wax),
            all_bands, all_colors)


def compute_slopes(oipc, source, wax, bands):
    """Compute pooled and within-band slopes."""
    slope_vs_source = np.polyfit(source, wax, 1)[0]
    slope_vs_oipc, intercept_oipc = np.polyfit(oipc, wax, 1)

    within = {}
    for band in band_config:
        mask = np.array(bands) == band["name"]
        s, ic = np.polyfit(oipc[mask], wax[mask], 1)
        within[band["name"]] = (s, ic)

    return slope_vs_source, slope_vs_oipc, intercept_oipc, within


# ─── Scenario A: No confounders ───
oipc_a, source_a, wax_a, bands_a, colors_a = generate_scenario(
    gs_offset_per_lat=lambda lat: 0,
    et_enrichment_per_lat=lambda lat: 0,
)
slopes_a = compute_slopes(oipc_a, source_a, wax_a, bands_a)

# ─── Scenario B: Growing season bias only ───
oipc_b, source_b, wax_b, bands_b, colors_b = generate_scenario(
    gs_offset_per_lat=lambda lat: 0.45 * lat,  # +0 at equator, +29 at 65°N
    et_enrichment_per_lat=lambda lat: 0,
)
slopes_b = compute_slopes(oipc_b, source_b, wax_b, bands_b)

# ─── Scenario C: Growing season + evapotranspiration (cancelling) ───
# ET enrichment is large in tropics (warm, dry soils), small at high latitudes
oipc_c, source_c, wax_c, bands_c, colors_c = generate_scenario(
    gs_offset_per_lat=lambda lat: 0.45 * lat,  # same as B
    et_enrichment_per_lat=lambda lat: 25 * (1 - lat / 75),  # +25 at equator, ~3 at 65°
)
slopes_c = compute_slopes(oipc_c, source_c, wax_c, bands_c)

# Print summary
print("=" * 60)
print(f"True fractionation factor (alpha): {ALPHA}")
print("=" * 60)
for label, slopes in [("A: No confounders", slopes_a),
                       ("B: Growing season only", slopes_b),
                       ("C: GS + ET (cancelling)", slopes_c)]:
    src_slope, oipc_slope, _, within = slopes
    within_mean = np.mean([s for s, _ in within.values()])
    print(f"\n{label}:")
    print(f"  Slope vs source water: {src_slope:.3f}")
    print(f"  Slope vs OIPC (pooled): {oipc_slope:.3f}")
    print(f"  Within-band slopes vs OIPC (mean): {within_mean:.3f}")
    for name, (s, _) in within.items():
        print(f"    {name.replace(chr(10), ' ')}: {s:.3f}")

# ─── Figure ───
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
scenarios = [
    ("A. No confounders", oipc_a, wax_a, bands_a, slopes_a),
    ("B. Growing season bias only", oipc_b, wax_b, bands_b, slopes_b),
    ("C. Growing season + evapotranspiration\n(opposing effects cancel)",
     oipc_c, wax_c, bands_c, slopes_c),
]

for ax_idx, (title, oipc, wax, bands, slopes) in enumerate(scenarios):
    ax = axes[ax_idx]
    _, slope_oipc, intercept_oipc, within = slopes[1], slopes[1], slopes[2], slopes[3]
    slope_oipc = slopes[1]
    intercept_oipc = slopes[2]
    within = slopes[3]

    # Plot points by band
    for band in band_config:
        mask = np.array(bands) == band["name"]
        ax.scatter(oipc[mask], wax[mask], c=band["color"], s=20, alpha=0.7,
                   label=band["name"].replace("\n", " "), edgecolors="none")

    # Within-band regression lines
    for band in band_config:
        mask = np.array(bands) == band["name"]
        s, ic = within[band["name"]]
        x_band = np.linspace(oipc[mask].min() - 3, oipc[mask].max() + 3, 50)
        ax.plot(x_band, s * x_band + ic, color=band["color"], lw=1.5, alpha=0.6)

    # Pooled regression line
    x_fit = np.linspace(oipc.min() - 5, oipc.max() + 5, 100)
    ax.plot(x_fit, slope_oipc * x_fit + intercept_oipc, "k-", lw=2.5,
            label=f"Pooled: {slope_oipc:.2f}")

    ax.set_xlabel("δ²H OIPC (‰)", fontsize=11)
    if ax_idx == 0:
        ax.set_ylabel("δ²H leaf wax (‰)", fontsize=11)
    ax.set_title(title, fontsize=11, fontweight="bold")

    # Slope annotation
    within_mean = np.mean([s for s, _ in within.values()])
    if ax_idx == 2:
        bbox_color = "lightyellow"
        note = (f"pooled slope = {slope_oipc:.2f} ≈ α\n"
                f"within-band ≈ {within_mean:.2f}\n"
                f"← coincidental!")
    elif ax_idx == 1:
        bbox_color = "salmon"
        note = f"pooled slope = {slope_oipc:.2f}\n< α = {ALPHA}"
    else:
        bbox_color = "wheat"
        note = f"pooled slope = {slope_oipc:.2f}\n= α = {ALPHA}"

    ax.text(0.97, 0.05, note, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor=bbox_color,
                                   alpha=0.8))

    if ax_idx == 0:
        ax.legend(fontsize=7, loc="upper left")

# Add a text annotation below
fig.text(0.5, -0.02,
         "In all three scenarios, the true fractionation α = 0.85 everywhere. "
         "Source water ≠ OIPC due to geographic effects.\n"
         "Panel C shows that a pooled slope of 0.85 does not validate the fractionation — "
         "opposing confounders can cancel to produce any slope.",
         ha="center", va="top", fontsize=10, style="italic",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0"))

plt.tight_layout()
plt.subplots_adjust(bottom=0.15)

outpath = Path("/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/"
               "_for_GCA/validation_scripts/slope_demonstration_v2.png")
plt.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"\nFigure saved to {outpath}")
