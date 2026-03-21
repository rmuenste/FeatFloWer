#!/usr/bin/env python3
"""
2D cross-section heatmap of the CylinderAttraction weight function.

Shows alpha(x,y) as a color map on a slice through the cylinder center,
matching the parameters in app_init.f90 :: CylinderAttraction.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# --- Parameters (must match app_init.f90) ---
dBandOut  = 0.15
dAlphaOut = 0.6
dBandIn   = 0.04
dAlphaIn  = 0.95

# Cylinder geometry (from q2p1_fac3d benchmark)
CylCenter = np.array([0.2, 0.2])
CylRadius = 0.05

# --- Grid ---
margin = dBandOut + 0.02
nx, ny = 800, 800
x = np.linspace(CylCenter[0] - CylRadius - margin,
                CylCenter[0] + CylRadius + margin, nx)
y = np.linspace(CylCenter[1] - CylRadius - margin,
                CylCenter[1] + CylRadius + margin, ny)
X, Y = np.meshgrid(x, y)

dx = X - CylCenter[0]
dy = Y - CylCenter[1]
rxy = np.sqrt(dx**2 + dy**2)
dist = rxy - CylRadius

# --- Compute alpha on grid ---
alpha = np.zeros_like(dist)

# Exterior: quadratic decay from dAlphaOut at surface to 0 at dBandOut
mask = (dist > 0) & (dist <= dBandOut)
alpha[mask] = dAlphaOut * ((dBandOut - dist[mask]) / dBandOut) ** 2

# Interior: linear fade from -dAlphaIn at surface to 0 at dBandIn
mask = (dist < 0) & (np.abs(dist) <= dBandIn)
alpha[mask] = -dAlphaIn * (1.0 - np.abs(dist[mask]) / dBandIn)

# --- Plot ---
fig, ax = plt.subplots(figsize=(8, 8))

# Diverging colormap: blue (interior/push) -> white (zero) -> red (exterior/pull)
norm = TwoSlopeNorm(vmin=-dAlphaIn, vcenter=0, vmax=dAlphaOut)
im = ax.pcolormesh(X, Y, alpha, cmap='RdBu_r', norm=norm, shading='auto')

# Draw cylinder surface
theta = np.linspace(0, 2 * np.pi, 200)
ax.plot(CylCenter[0] + CylRadius * np.cos(theta),
        CylCenter[1] + CylRadius * np.sin(theta),
        'k-', linewidth=2, label='cylinder surface')

# Draw band boundaries
for r, ls, lbl in [
    (CylRadius + dBandOut, '-.', f'dBandOut = {dBandOut}'),
    (CylRadius - dBandIn, ':', f'dBandIn = {dBandIn}'),
]:
    ax.plot(CylCenter[0] + r * np.cos(theta),
            CylCenter[1] + r * np.sin(theta),
            color='gray', ls=ls, linewidth=1, label=lbl)

cbar = fig.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('alpha (attraction weight)', fontsize=12)

ax.set_xlabel('x', fontsize=12)
ax.set_ylabel('y', fontsize=12)
ax.set_title('CylinderAttraction weight field (cross-section)', fontsize=14)
ax.set_aspect('equal')
ax.legend(loc='upper right', fontsize=9)

plt.tight_layout()
plt.savefig('attraction_weights_2d.png', dpi=150)
print("Saved: attraction_weights_2d.png")
plt.show()
