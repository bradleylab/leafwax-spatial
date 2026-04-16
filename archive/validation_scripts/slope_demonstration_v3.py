"""
slope_demonstration_v3.py
The simplest possible demonstration that the regression slope of d2H_wax vs OIPC
is not the fractionation factor.

Two-site arithmetic: same fractionation (alpha = 0.85), same wax values,
different x-axis (source water vs OIPC) → different slopes.
No statistics, no models — just division.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path

# ─── The two sites ───
ALPHA = 0.85
EPSILON = (ALPHA - 1) * 1000  # -150

sites = {
    "Equatorial lake": {
        "oipc": -20,
        "gs_offset": 5,  # small growing season offset
    },
    "Polar lake": {
        "oipc": -120,
        "gs_offset": 40,  # large: summer precip >> annual mean
    },
}

for name, s in sites.items():
    s["source"] = s["oipc"] + s["gs_offset"]
    s["wax"] = ALPHA * s["source"] + EPSILON

eq = sites["Equatorial lake"]
po = sites["Polar lake"]

wax_range = eq["wax"] - po["wax"]
source_range = eq["source"] - po["source"]
oipc_range = eq["oipc"] - po["oipc"]

slope_vs_source = wax_range / source_range
slope_vs_oipc = wax_range / oipc_range

print(f"Equatorial: OIPC={eq['oipc']}  source={eq['source']}  wax={eq['wax']:.1f}")
print(f"Polar:      OIPC={po['oipc']}  source={po['source']}  wax={po['wax']:.1f}")
print(f"d2H_wax range: {wax_range:.2f}")
print(f"Source water range: {source_range}  → slope = {slope_vs_source:.3f}")
print(f"OIPC range: {oipc_range}  → slope = {slope_vs_oipc:.3f}")

# ─── Figure: 2 panels ───
fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)

marker_props = dict(s=120, zorder=5, edgecolors="black", linewidths=1)

# ─── Panel A: d2H_wax vs source water ───
ax = axes[0]

ax.scatter(eq["source"], eq["wax"], c="#e07b39", marker="o",
           label="Equatorial lake", **marker_props)
ax.scatter(po["source"], po["wax"], c="#7b2d8e", marker="s",
           label="Polar lake", **marker_props)

# Regression line
x_range = np.array([po["source"] - 10, eq["source"] + 10])
ax.plot(x_range, ALPHA * x_range + EPSILON, "k-", lw=2, zorder=3)

# Annotations
ax.annotate(f'Source = {eq["source"]}‰\n(OIPC {eq["oipc"]}‰ + {eq["gs_offset"]}‰ offset)',
            xy=(eq["source"], eq["wax"]), xytext=(eq["source"] + 8, eq["wax"] - 12),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", color="gray"),
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#fff3e0"))

ax.annotate(f'Source = {po["source"]}‰\n(OIPC {po["oipc"]}‰ + {po["gs_offset"]}‰ offset)',
            xy=(po["source"], po["wax"]), xytext=(po["source"] + 12, po["wax"] + 10),
            fontsize=8.5, arrowprops=dict(arrowstyle="->", color="gray"),
            bbox=dict(boxstyle="round,pad=0.3", facecolor="#f3e5f5"))

# Range bracket (vertical)
ax.annotate("", xy=(po["source"] - 6, po["wax"]),
            xytext=(po["source"] - 6, eq["wax"]),
            arrowprops=dict(arrowstyle="<->", color="darkred", lw=1.5))
ax.text(po["source"] - 9, (eq["wax"] + po["wax"]) / 2,
        f"Δwax\n{wax_range:.1f}‰", fontsize=8, ha="right", color="darkred",
        fontweight="bold")

# Range bracket (horizontal)
ax.annotate("", xy=(po["source"], po["wax"] + 5),
            xytext=(eq["source"], po["wax"] + 5),
            arrowprops=dict(arrowstyle="<->", color="navy", lw=1.5))
ax.text((eq["source"] + po["source"]) / 2, po["wax"] + 8,
        f"Δsource = {source_range}‰", fontsize=8, ha="center", color="navy",
        fontweight="bold")

ax.set_xlabel("δ²H source water (‰)", fontsize=12)
ax.set_ylabel("δ²H leaf wax (‰)", fontsize=12)
ax.set_title("A. Wax vs source water", fontsize=13, fontweight="bold")
ax.text(0.97, 0.05, f"slope = {wax_range:.1f} / {source_range} = {slope_vs_source:.2f}\n= α (fractionation factor)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.9),
        fontweight="bold")

ax.set_xlim(-135, 5)

# ─── Panel B: d2H_wax vs OIPC ───
ax = axes[1]

ax.scatter(eq["oipc"], eq["wax"], c="#e07b39", marker="o",
           label="Equatorial lake", **marker_props)
ax.scatter(po["oipc"], po["wax"], c="#7b2d8e", marker="s",
           label="Polar lake", **marker_props)

# Regression line through the two points
slope_line = slope_vs_oipc
int_line = eq["wax"] - slope_line * eq["oipc"]
x_range = np.array([po["oipc"] - 10, eq["oipc"] + 10])
ax.plot(x_range, slope_line * x_range + int_line, "k-", lw=2, zorder=3)

# Also show where the points WOULD be if source = OIPC (ghost points)
eq_ghost_wax = ALPHA * eq["oipc"] + EPSILON
po_ghost_wax = ALPHA * po["oipc"] + EPSILON
ax.scatter(eq["oipc"], eq_ghost_wax, c="#e07b39", marker="o", s=80, alpha=0.25,
           edgecolors="gray", linestyles="dashed", zorder=4)
ax.scatter(po["oipc"], po_ghost_wax, c="#7b2d8e", marker="s", s=80, alpha=0.25,
           edgecolors="gray", linestyles="dashed", zorder=4)

# Ghost regression line (slope = 0.85)
ax.plot(x_range, ALPHA * x_range + EPSILON, "k--", lw=1.5, alpha=0.3, zorder=2,
        label="slope = 0.85 (if source = OIPC)")

# Arrows from ghost to actual (showing the offset effect)
ax.annotate("", xy=(eq["oipc"], eq["wax"]),
            xytext=(eq["oipc"], eq_ghost_wax),
            arrowprops=dict(arrowstyle="->", color="#e07b39", lw=2, alpha=0.6))
ax.text(eq["oipc"] + 3, (eq["wax"] + eq_ghost_wax) / 2,
        f"+{eq['gs_offset']}‰\noffset", fontsize=8, color="#e07b39", alpha=0.8)

ax.annotate("", xy=(po["oipc"], po["wax"]),
            xytext=(po["oipc"], po_ghost_wax),
            arrowprops=dict(arrowstyle="->", color="#7b2d8e", lw=2, alpha=0.6))
ax.text(po["oipc"] + 3, (po["wax"] + po_ghost_wax) / 2,
        f"+{po['gs_offset']}‰\noffset", fontsize=8, color="#7b2d8e", alpha=0.8)

# Range bracket (vertical) - same wax range
ax.annotate("", xy=(po["oipc"] - 6, po["wax"]),
            xytext=(po["oipc"] - 6, eq["wax"]),
            arrowprops=dict(arrowstyle="<->", color="darkred", lw=1.5))
ax.text(po["oipc"] - 9, (eq["wax"] + po["wax"]) / 2,
        f"Δwax\n{wax_range:.1f}‰\n(same!)", fontsize=8, ha="right", color="darkred",
        fontweight="bold")

# Range bracket (horizontal) - larger OIPC range
ax.annotate("", xy=(po["oipc"], po["wax"] + 5),
            xytext=(eq["oipc"], po["wax"] + 5),
            arrowprops=dict(arrowstyle="<->", color="navy", lw=1.5))
ax.text((eq["oipc"] + po["oipc"]) / 2, po["wax"] + 8,
        f"ΔOIPC = {oipc_range}‰ (wider!)", fontsize=8, ha="center", color="navy",
        fontweight="bold")

ax.set_xlabel("δ²H OIPC (‰)", fontsize=12)
ax.set_title("B. Same wax values vs OIPC", fontsize=13, fontweight="bold")
ax.text(0.97, 0.05, f"slope = {wax_range:.1f} / {oipc_range} = {slope_vs_oipc:.2f}\n≠ α  (NOT the fractionation)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=11,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="salmon", alpha=0.9),
        fontweight="bold")
ax.legend(fontsize=8, loc="upper left")

ax.set_xlim(-135, 5)

# Bottom annotation
fig.text(0.5, -0.04,
         "The fractionation α = 0.85 is identical at both sites. The d2H_wax values are identical in both panels.\n"
         "The only difference is the x-axis: source water (what the plant used) vs OIPC (what we measure).\n"
         "Because the OIPC-to-source-water offset is larger at the polar site (+40‰ vs +5‰),\n"
         "the OIPC range (100‰) exceeds the source water range (65‰), and the slope drops from 0.85 to 0.55.",
         ha="center", va="top", fontsize=10, style="italic",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f0"))

plt.tight_layout()
plt.subplots_adjust(bottom=0.22)

outpath = Path("/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/"
               "_for_GCA/validation_scripts/slope_demonstration_v3.png")
plt.savefig(outpath, dpi=200, bbox_inches="tight")
print(f"\nFigure saved to {outpath}")
