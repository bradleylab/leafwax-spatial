"""
slope_demonstration.py
Generates synthetic data and figure demonstrating why the regression slope
of d2H_wax vs OIPC differs from the fractionation factor (alpha = 0.85).

Uses growing season bias as the mechanism: plants preferentially use
summer precipitation, which is more enriched than the annual mean,
especially at high latitudes.

Produces a 3-panel figure:
  A: d2H_wax vs source water (slope = alpha = 0.85 at every site)
  B: d2H_wax vs OIPC (slope < 0.85 due to growing season compression)
  C: Same as B but with within-band regression lines showing Simpson's paradox

References:
  Hayes (2001) Rev Mineral Geochem 43:225-77
  Sachse et al. (2012) Annu Rev Earth Planet Sci 40:221-49
  Daniels et al. (2017) GCA 215:105-19
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

np.random.seed(42)

# ─── Parameters ───
ALPHA = 0.85  # Fractionation factor (true, constant everywhere)
EPSILON = (ALPHA - 1) * 1000  # = -150 permil

# Four latitude bands
bands = {
    "Tropical\n(0-15°)": {"lat_center": 7, "oipc_center": -25, "gs_offset": 5},
    "Subtropical\n(15-35°)": {"lat_center": 25, "oipc_center": -50, "gs_offset": 12},
    "Temperate\n(35-55°)": {"lat_center": 45, "oipc_center": -75, "gs_offset": 22},
    "High latitude\n(55-75°)": {"lat_center": 65, "oipc_center": -105, "gs_offset": 32},
}

colors = ["#e07b39", "#2ca02c", "#1f77b4", "#7b2d8e"]

# ─── Generate pseudodata ───
# Within each band, 25 sites with OIPC varying due to elevation/local climate
n_per_band = 25

all_oipc = []
all_source = []
all_wax = []
all_band = []
all_color = []

for i, (name, params) in enumerate(bands.items()):
    # Within-band OIPC variation (elevation, local climate)
    oipc = params["oipc_center"] + np.random.normal(0, 12, n_per_band)

    # Growing season offset (roughly constant within band, with small noise)
    gs_offset = params["gs_offset"] + np.random.normal(0, 2, n_per_band)

    # Source water = OIPC + growing season offset
    source = oipc + gs_offset

    # d2H_wax = alpha * source + epsilon + measurement noise
    wax = ALPHA * source + EPSILON + np.random.normal(0, 4, n_per_band)

    all_oipc.extend(oipc)
    all_source.extend(source)
    all_wax.extend(wax)
    all_band.extend([name] * n_per_band)
    all_color.extend([colors[i]] * n_per_band)

all_oipc = np.array(all_oipc)
all_source = np.array(all_source)
all_wax = np.array(all_wax)

# ─── Compute regression slopes ───
# Overall: d2H_wax vs source water
slope_source, intercept_source = np.polyfit(all_source, all_wax, 1)

# Overall: d2H_wax vs OIPC
slope_oipc, intercept_oipc = np.polyfit(all_oipc, all_wax, 1)

# Within-band: d2H_wax vs OIPC
band_names = list(bands.keys())
within_slopes = {}
within_intercepts = {}
for name in band_names:
    mask = np.array(all_band) == name
    s, ic = np.polyfit(all_oipc[mask], all_wax[mask], 1)
    within_slopes[name] = s
    within_intercepts[name] = ic

print(f"Fractionation factor (alpha): {ALPHA:.3f}")
print(f"Slope of d2H_wax vs source water: {slope_source:.3f}")
print(f"Slope of d2H_wax vs OIPC (pooled): {slope_oipc:.3f}")
print(f"Within-band slopes vs OIPC:")
for name in band_names:
    print(f"  {name.replace(chr(10), ' ')}: {within_slopes[name]:.3f}")

# ─── Figure ───
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Panel A: d2H_wax vs source water
ax = axes[0]
for i, name in enumerate(band_names):
    mask = np.array(all_band) == name
    ax.scatter(all_source[mask], all_wax[mask], c=colors[i], s=25, alpha=0.7,
               label=name.replace("\n", " "), edgecolors="none")

x_fit = np.linspace(min(all_source) - 5, max(all_source) + 5, 100)
ax.plot(x_fit, slope_source * x_fit + intercept_source, "k-", lw=2,
        label=f"Slope = {slope_source:.2f}")

ax.set_xlabel("δ²H source water (‰)", fontsize=11)
ax.set_ylabel("δ²H leaf wax (‰)", fontsize=11)
ax.set_title("A. Wax vs source water\n(the true reactant)", fontsize=12,
             fontweight="bold")
ax.legend(fontsize=8, loc="upper left")
ax.text(0.95, 0.05, f"slope = {slope_source:.2f}\n≈ α = {ALPHA}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.8))

# Panel B: d2H_wax vs OIPC (pooled regression only)
ax = axes[1]
for i, name in enumerate(band_names):
    mask = np.array(all_band) == name
    ax.scatter(all_oipc[mask], all_wax[mask], c=colors[i], s=25, alpha=0.7,
               edgecolors="none")

x_fit = np.linspace(min(all_oipc) - 5, max(all_oipc) + 5, 100)
ax.plot(x_fit, slope_oipc * x_fit + intercept_oipc, "k-", lw=2)

# Draw arrows showing the "shift" for one example point per band
for i, name in enumerate(band_names):
    mask = np.array(all_band) == name
    # Pick median point
    idx = np.where(mask)[0][n_per_band // 2]
    ax.annotate("", xy=(all_oipc[idx], all_wax[idx]),
                xytext=(all_source[idx], all_wax[idx]),
                arrowprops=dict(arrowstyle="->", color=colors[i], lw=1.5,
                                alpha=0.5))

ax.set_xlabel("δ²H OIPC (‰)", fontsize=11)
ax.set_ylabel("δ²H leaf wax (‰)", fontsize=11)
ax.set_title("B. Wax vs OIPC\n(the proxy for the reactant)", fontsize=12,
             fontweight="bold")
ax.text(0.95, 0.05, f"pooled slope = {slope_oipc:.2f}\n< α = {ALPHA}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="salmon", alpha=0.8))

# Panel C: Same as B but with within-band regression lines
ax = axes[2]
for i, name in enumerate(band_names):
    mask = np.array(all_band) == name
    ax.scatter(all_oipc[mask], all_wax[mask], c=colors[i], s=25, alpha=0.7,
               label=name.replace("\n", " "), edgecolors="none")

    # Within-band regression line
    x_band = np.linspace(all_oipc[mask].min() - 3, all_oipc[mask].max() + 3, 50)
    ax.plot(x_band, within_slopes[name] * x_band + within_intercepts[name],
            color=colors[i], lw=2, alpha=0.8)

# Pooled regression (dashed)
x_fit = np.linspace(min(all_oipc) - 5, max(all_oipc) + 5, 100)
ax.plot(x_fit, slope_oipc * x_fit + intercept_oipc, "k--", lw=2, alpha=0.5,
        label=f"Pooled: {slope_oipc:.2f}")

ax.set_xlabel("δ²H OIPC (‰)", fontsize=11)
ax.set_ylabel("δ²H leaf wax (‰)", fontsize=11)
ax.set_title("C. Simpson's paradox\n(within-band vs pooled slopes)", fontsize=12,
             fontweight="bold")

# Annotation
within_mean = np.mean(list(within_slopes.values()))
ax.text(0.95, 0.05,
        f"within-band slopes ≈ {within_mean:.2f}\npooled slope = {slope_oipc:.2f}",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

plt.tight_layout()

outpath = Path("/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/_for_GCA/validation_scripts/slope_demonstration.png")
plt.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"\nFigure saved to {outpath}")

# Also save the pseudodata
import csv
data_path = outpath.with_suffix(".csv")
with open(data_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["band", "oipc", "source_water", "d2h_wax",
                     "growing_season_offset"])
    for i in range(len(all_oipc)):
        writer.writerow([all_band[i].replace("\n", " "),
                         f"{all_oipc[i]:.1f}",
                         f"{all_source[i]:.1f}",
                         f"{all_wax[i]:.1f}",
                         f"{all_source[i] - all_oipc[i]:.1f}"])
print(f"Data saved to {data_path}")
