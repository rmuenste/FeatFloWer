#!/usr/bin/env python3
"""Compute Reynolds/Dean numbers of the ATC slug flow from VTK output.

Usage: python3 compute_Re.py [step]   (default: latest main.NNNNN.pvtu in _vtk/)

Reads all per-rank pieces of one time frame, measures velocity scales in the
DATA (not the BC formula), and forms:
  Re_wall = U_wall_max * d / nu      (imposed-BC sanity check ~ Omega*r)
  Re_bulk = <|u|>_nodes * d / nu     (bulk response scale)
  Re_ax   = max|u_axial| * d / nu    (conveyor-cell strength; u.e_phi about z)
  Re_sec  = max|u_secondary| * d/nu  (in-plane Dean-vortex strength)
  De      = Re_ax * sqrt(a / Rc)     (Dean number of the curved segment)
Geometry (from mesh/centerline, CGS): a=0.2375 cm, d=0.4750 cm, Rc=2.521 cm.
Fluid: ultra-pure water, nu = 0.01 cm^2/s.
Laminar evidence: compare two consecutive frames -> steady => laminar.
"""
import glob, math, re, sys
import numpy as np
import meshio

NU = 0.01          # cm^2/s (water, CGS)
A_BORE = 0.2375    # cm
D_BORE = 2 * A_BORE
R_COIL = 2.521     # cm (centerline radius about drum z-axis)

def frame_steps():
    steps = sorted(int(re.search(r"main\.(\d+)\.pvtu", f).group(1))
                   for f in glob.glob("_vtk/main.*.pvtu"))
    return steps

def load_frame(step):
    pieces = sorted(glob.glob(f"_vtk/res_node_*.{step:05d}.vtu"))
    if not pieces:
        sys.exit(f"no pieces for step {step}")
    pts, vel = [], []
    for p in pieces:
        m = meshio.read(p)
        pts.append(m.points)
        vel.append(m.point_data["Velocity"])
    return np.vstack(pts), np.vstack(vel)

def analyze(step):
    x, u = load_frame(step)
    speed = np.linalg.norm(u, axis=1)
    r = np.hypot(x[:, 0], x[:, 1])
    ephi = np.column_stack([-x[:, 1] / r, x[:, 0] / r, np.zeros(len(r))])
    u_ax = np.einsum("ij,ij->i", u, ephi)          # along drum azimuth ~ tube axis
    u_sec = np.linalg.norm(u - u_ax[:, None] * ephi, axis=1)
    print(f"frame step {step}: {len(x)} points")
    print(f"  max |u|            = {speed.max():8.4f} cm/s")
    print(f"  <|u|> (node mean)  = {speed.mean():8.4f} cm/s")
    print(f"  u_axial min/max    = {u_ax.min():8.4f} / {u_ax.max():8.4f} cm/s")
    print(f"  max |u_secondary|  = {u_sec.max():8.4f} cm/s")
    re_wall = speed.max() * D_BORE / NU
    re_bulk = speed.mean() * D_BORE / NU
    re_ax = max(abs(u_ax.min()), abs(u_ax.max())) * D_BORE / NU
    re_sec = u_sec.max() * D_BORE / NU
    dean = re_ax * math.sqrt(A_BORE / R_COIL)
    print(f"  Re_wall = {re_wall:7.1f}   Re_bulk = {re_bulk:7.1f}"
          f"   Re_ax = {re_ax:7.1f}   Re_sec = {re_sec:7.1f}   De = {dean:7.1f}")
    return u

if __name__ == "__main__":
    steps = frame_steps()
    step = int(sys.argv[1]) if len(sys.argv) > 1 else steps[-1]
    u_last = analyze(step)
    # steadiness (laminar evidence): relative change vs previous frame
    prev = [s for s in steps if s < step]
    if prev:
        _, u_prev = load_frame(prev[-1])
        if u_prev.shape == u_last.shape:
            num = np.linalg.norm(u_last - u_prev)
            den = np.linalg.norm(u_last)
            print(f"  steadiness ||u(t)-u(t-1)||/||u(t)|| = {num/den:.3e} "
                  f"(vs frame {prev[-1]}; << 1 and shrinking => steady/laminar)")
