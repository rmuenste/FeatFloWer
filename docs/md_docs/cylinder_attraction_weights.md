# Cylinder Attraction Weight Function

This document describes the attraction weight function used by the
`CylinderAttraction` smoother in `applications/q2p1_fac3d/app_init.f90`.
The smoother concentrates mesh nodes near the cylinder surface to improve
resolution of the immersed boundary.

## Geometry

Let $\mathbf{c} = (c_x, c_y)$ be the cylinder center (projected onto the
$xy$-plane), $R$ the cylinder radius, and $\delta$ the radius offset
(`dRadiusOffset`).  For each mesh node $\mathbf{x}_i = (x_i, y_i, z_i)$
define:

$$
\Delta x = x_i - c_x, \qquad \Delta y = y_i - c_y, \qquad
r_{xy} = \sqrt{\Delta x^2 + \Delta y^2}
$$

### Signed distance

The signed radial distance from the offset surface is:

$$
d = r_{xy} - R - \delta
$$

where $d > 0$ means exterior, $d < 0$ means interior.

### Attraction target

The nearest point on the attraction radius $R_a = R - \delta$ is:

$$
\mathbf{p}_i = \mathbf{c} + R_a \, \frac{(\Delta x,\, \Delta y)}{r_{xy}}
$$

Nodes with $r_{xy} < 10^{-14}$ (essentially at the cylinder axis) are
skipped.

## Weight function $\alpha(d)$

The weight $\alpha$ controls how strongly a node is pulled toward the
attraction target.  It is defined piecewise:

### Exterior band ($0 < d \le B_\mathrm{out}$)

Quadratic decay from the peak at the surface to zero at the outer band
limit:

$$
\alpha(d) = \alpha_\mathrm{out} \left(\frac{B_\mathrm{out} - d}{B_\mathrm{out}}\right)^{\!2}
$$

| Value | At $d = 0^+$ | At $d = B_\mathrm{out}$ |
|-------|-------------|------------------------|
| $\alpha$ | $\alpha_\mathrm{out}$ | $0$ |

### Interior band ($-B_\mathrm{in} \le d < 0$)

Linear decay from the peak at the surface to zero at the inner band limit:

$$
\alpha(d) = \alpha_\mathrm{in} \left(1 - \frac{|d|}{B_\mathrm{in}}\right)
$$

| Value | At $d = 0^-$ | At $d = -B_\mathrm{in}$ |
|-------|-------------|------------------------|
| $\alpha$ | $\alpha_\mathrm{in}$ | $0$ |

### Outside both bands

$$
\alpha(d) = 0 \qquad \text{for } d > B_\mathrm{out} \text{ or } d < -B_\mathrm{in}
$$

## Node displacement

The node position is updated in the $xy$-plane only (the $z$-coordinate is
unchanged):

$$
\mathbf{x}_i \;\leftarrow\; \mathbf{x}_i + \omega\,\alpha(d)\,(\mathbf{p}_i - \mathbf{x}_i)
$$

where $\omega$ is the relaxation parameter (`dOmega`).  The update is
applied for `nIter` iterations, each followed by a boundary
re-parametrization step.

The effective per-iteration displacement magnitude is:

$$
|\Delta \mathbf{x}_i| = \omega\,\alpha(d)\,|\mathbf{p}_i - \mathbf{x}_i|
$$

## Parameters

| Symbol | Parameter | Default | Description |
|--------|-----------|---------|-------------|
| $B_\mathrm{out}$ | `dBandOut` | 0.15 | Outer band width |
| $\alpha_\mathrm{out}$ | `dAlphaOut` | 0.8 | Peak exterior weight (at surface) |
| $B_\mathrm{in}$ | `dBandIn` | 0.04 | Inner band width |
| $\alpha_\mathrm{in}$ | `dAlphaIn` | 0.95 | Peak interior weight (at surface) |
| $\delta$ | `dRadiusOffset` | 0.0 | Inward shift of attraction target and distance origin |
| $\omega$ | `dOmega` | 0.2 | Relaxation factor (passed at call site) |
| — | `nIter` | 8 | Number of attraction iterations (loop in call site) |

## Visualization weight

The `ComputeCylAttractionWeight` subroutine computes the same $\alpha(d)$
for VTK output.  Interior weights are stored as **negative** values to
distinguish push (interior) from pull (exterior) in ParaView:

$$
w(d) = \begin{cases}
\alpha_\mathrm{out}\left(\dfrac{B_\mathrm{out}-d}{B_\mathrm{out}}\right)^{\!2} & 0 < d \le B_\mathrm{out} \\[6pt]
-\alpha_\mathrm{in}\left(1-\dfrac{|d|}{B_\mathrm{in}}\right) & -B_\mathrm{in} \le d < 0 \\[6pt]
0 & \text{otherwise}
\end{cases}
$$
