#!/usr/bin/env python3
"""
1D profile plot of the CylinderAttraction weight function.

Plots alpha(dist) for the exterior and interior bands, matching the
parameters in app_init.f90 :: CylinderAttraction.
"""

import numpy as np
import matplotlib.pyplot as plt

# --- Parameters (must match app_init.f90) ---
dBandOut  = 0.15
dAlphaOut = 0.6
dBandIn   = 0.04
dAlphaIn  = 0.95

# --- Compute alpha(dist) ---
dist = np.linspace(-dBandIn * 1.2, dBandOut * 1.1, 2000)
alpha = np.zeros_like(dist)

for i, d in enumerate(dist):
    if d > 0.0 and d <= dBandOut:
        alpha[i] = dAlphaOut * ((dBandOut - d) / dBandOut) ** 2
    elif d < 0.0:
        absDist = abs(d)
        if absDist <= dBandIn:
            alpha[i] = -dAlphaIn * (1.0 - absDist / dBandIn)

# --- Plot ---
fig, ax = plt.subplots(figsize=(9, 5))

ax.plot(dist, alpha, 'b-', linewidth=2)
ax.axvline(0, color='r', ls='--', linewidth=1.5, label='cylinder surface (dist=0)')
ax.axhline(0, color='k', ls='-', linewidth=0.5)

# Mark band boundaries
ax.axvline(dBandOut, color='gray', ls='-.', label=f'dBandOut = {dBandOut}')
ax.axvline(-dBandIn, color='orange', ls=':', label=f'-dBandIn = {-dBandIn}')

# Annotate key values
ax.annotate(f'dAlphaOut = {dAlphaOut} (peak)',
            xy=(0.001, dAlphaOut), xytext=(0.04, dAlphaOut + 0.15),
            arrowprops=dict(arrowstyle='->', color='green'),
            fontsize=10, color='green')
ax.annotate(f'-dAlphaIn = {-dAlphaIn}',
            xy=(0.0, -dAlphaIn), xytext=(-0.03, -dAlphaIn - 0.15),
            arrowprops=dict(arrowstyle='->', color='purple'),
            fontsize=10, color='purple')

ax.set_xlabel('dist  (signed radial distance from cylinder surface)', fontsize=12)
ax.set_ylabel('alpha  (attraction weight)', fontsize=12)
ax.set_title('CylinderAttraction weight function', fontsize=14)
ax.legend(loc='upper right', fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_xlim(dist[0], dist[-1])
ax.set_ylim(-1.2, 1.0)

plt.tight_layout()
plt.savefig('attraction_weights_1d.png', dpi=150)
print("Saved: attraction_weights_1d.png")
plt.show()
