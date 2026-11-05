#!/usr/bin/env python3
"""
d52_composite_target.py - mass-conserving composite closure target for the
D5.2 annular Couette numerical viscometer (rows d52_v21_einstein,
d52_v22_phi10, d52_v23_phi20, d52_v25f_l4_poscont / _verdict).

The "composite target" is the closure prediction evaluated pointwise on the
MEASURED local volume fraction phi(r,z) of the sphere cloud that is actually
in the instrument (RSA clouds keep 0.5 d surface clearances from both
cylinders, so the wall layers are particle-free and the interior is denser
than the nominal phi; the torque ratio must read below the naive closure).

Construction
------------
1. Logged sphere positions (particle_force.log columns px py pz, or a
   pre-extracted "t id x y z" file) over a plateau window are sampled at a
   fixed cadence (default: 21 snapshots at 1 t.u.).
2. MASS-CONSERVING deposit: each sphere is represented by M Monte Carlo
   sample points drawn uniformly inside the sphere (deterministic seed).
   Each point carries V_sphere/M and is assigned to the (r,z) cell that
   contains it - nothing is clipped, nothing is smeared over a half-radius
   footprint. Points that fall outside the annulus (cannot happen for a
   cloud that respects the walls; counted and reported) are moved to the
   nearest cell so that the deposited volume is exactly N*V_sphere.
   Consequently sum_cells phi_cell*V_cell == N_spheres*V_sphere and the
   global phi is N*V_sphere/V_annulus by construction; the tool prints the
   conservation ratio (must be 1.000000).
3. Closure applied pointwise: Einstein 1+2.5phi, Batchelor
   1+2.5phi+6.2phi^2, Krieger-Dougherty (1-phi/0.64)^(-1.6).
4. Radial layers act in series for the torque. In a Couette cell with a
   radially stratified viscosity the torque is constant across r,
   T = 2 pi H r^3 eta(r) d(omega)/dr, so
       eta_eff = int dr/r^3  /  int dr/(eta(r) r^3)
   (harmonic mean with weights w_j = (r_j^-2 - r_{j+1}^-2)/2). This is the
   default ("--series-weight r3"). "--series-weight uniform" gives the plain
   harmonic mean (equal-thickness layers), kept for comparison with older
   evaluations.
5. Horizontal slices act in parallel (equal heights): arithmetic mean over z.

Usage
-----
  tools/d52_composite_target.py POSFILE [POSFILE ...] --closure einstein \
      --tmin 230 --tmax 250 --dt-snap 1.0 --grid 20 --samples 2000 --seed 1

  When several files carry the same snapshot time (segment seams of a
  chained run) the LAST file on the command line wins (the continuing
  segment).

  --convergence   runs the (2000, 20000) x (20, 40) matrix and prints a
                  table (each with both series weights).
  --anchors       prints the naive closure at the nominal phi and the
                  wall-depletion-only estimate (uniform phi over the
                  occupied band r in [r_i+d/2, r_a-d/2], zero in the
                  0.5 d clearance layers).
"""
import argparse
import math
import sys

import numpy as np

# instrument (D5.2 through-hole Couette cell, d=1 units)
R_I = 5.0
R_A = 10.0
Z_0 = 0.0
Z_1 = 10.0
D_SPHERE = 1.0

CLOSURES = {
    "einstein": lambda p: 1.0 + 2.5 * p,
    "batchelor": lambda p: 1.0 + 2.5 * p + 6.2 * p * p,
    "kd": lambda p: np.power(1.0 - np.minimum(p, 0.6399) / 0.64, -1.6),
}


def read_positions(paths, tmin, tmax, dt_snap, tol=2.5e-4):
    """Return dict time -> (N,3) array. Accepts particle_force.log (14 numeric
    columns: time ip fx fy fz tx ty tz px py pz vx vy vz), pre-extracted
    "t id x y z", or those prefixed by a filename tag column. Later files
    override earlier ones for the same snapshot time."""
    snaps = {}
    for path in paths:
        local = {}
        with open(path) as fh:
            for line in fh:
                s = line.split()
                if not s or s[0].startswith("#"):
                    continue
                try:
                    float(s[0])
                except ValueError:
                    s = s[1:]  # leading tag column
                t = float(s[0])
                if t < tmin - tol or t > tmax + tol:
                    continue
                k = (t - tmin) / dt_snap
                kr = round(k)
                if abs(k - kr) * dt_snap > tol:
                    continue
                if len(s) >= 14:
                    xyz = (float(s[8]), float(s[9]), float(s[10]))
                elif len(s) >= 5:
                    xyz = (float(s[2]), float(s[3]), float(s[4]))
                else:
                    raise ValueError(f"{path}: cannot parse line: {line!r}")
                local.setdefault(kr, []).append(xyz)
        for kr, pts in local.items():
            snaps[kr] = np.asarray(pts)
    return {tmin + k * dt_snap: v for k, v in sorted(snaps.items())}


def sample_sphere_points(rng, n):
    """n points uniform in the unit ball."""
    v = rng.normal(size=(n, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    r = np.cbrt(rng.random(n))
    return v * r[:, None]


def deposit(snaps, grid, samples, seed, radius=D_SPHERE / 2):
    """Mass-conserving MC deposit. Returns (phi[nr,nz], stats)."""
    nr = nz = grid
    r_edges = np.linspace(R_I, R_A, nr + 1)
    z_edges = np.linspace(Z_0, Z_1, nz + 1)
    v_cell = (math.pi * (r_edges[1:] ** 2 - r_edges[:-1] ** 2))[:, None] * np.diff(z_edges)[None, :]
    vol = np.zeros((nr, nz))
    rng = np.random.default_rng(seed)
    v_sphere = 4.0 / 3.0 * math.pi * radius ** 3
    n_spheres = 0
    n_out = 0
    rmin, rmax = np.inf, -np.inf
    for t, pos in snaps.items():
        n_spheres += len(pos)
        rc = np.hypot(pos[:, 0], pos[:, 1])
        rmin, rmax = min(rmin, rc.min()), max(rmax, rc.max())
        off = sample_sphere_points(rng, samples * len(pos)).reshape(len(pos), samples, 3) * radius
        p = pos[:, None, :] + off
        r = np.hypot(p[..., 0], p[..., 1]).ravel()
        z = p[..., 2].ravel()
        ir = np.searchsorted(r_edges, r, side="right") - 1
        iz = np.searchsorted(z_edges, z, side="right") - 1
        out = (ir < 0) | (ir >= nr) | (iz < 0) | (iz >= nz)
        n_out += int(out.sum())
        # conserve: fold any outside point into the nearest cell
        ir = np.clip(ir, 0, nr - 1)
        iz = np.clip(iz, 0, nz - 1)
        np.add.at(vol, (ir, iz), v_sphere / samples)
    # vol accumulates over snapshots: phi is the per-snapshot (time-averaged) field
    phi = vol / v_cell / len(snaps)
    v_annulus = math.pi * (R_A ** 2 - R_I ** 2) * (Z_1 - Z_0)
    n_per_snap = n_spheres / max(len(snaps), 1)
    stats = dict(
        n_snapshots=len(snaps),
        n_spheres_per_snapshot=n_per_snap,
        n_sample_points=n_spheres * samples,
        n_outside=n_out,
        deposited_volume=float(vol.sum()),
        expected_volume=n_spheres * v_sphere,
        conservation=float(vol.sum()) / (n_spheres * v_sphere),
        phi_global=float(vol.sum()) / len(snaps) / v_annulus,
        phi_nominal=n_per_snap * v_sphere / v_annulus,
        phi_max=float(phi.max()),
        r_center_min=float(rmin),
        r_center_max=float(rmax),
        r_edges=r_edges,
    )
    return phi, stats


def series_weights(r_edges, kind):
    if kind == "r3":
        return 0.5 * (r_edges[:-1] ** -2 - r_edges[1:] ** -2)
    if kind == "uniform":
        return np.diff(r_edges)
    raise ValueError(kind)


def composite(phi, r_edges, closure, weight_kind):
    eta = CLOSURES[closure](phi)  # [nr, nz]
    w = series_weights(r_edges, weight_kind)[:, None]
    eta_slice = w.sum() / (w / eta).sum(axis=0)  # series in r, per z slice
    return float(eta_slice.mean()), eta_slice  # parallel (mean) over z


def depletion_anchor(phi_nominal, closure, weight_kind):
    """Uniform phi over the occupied band [R_I+d/2, R_A-d/2], zero in the
    clearance layers - the wall-depletion-only estimate."""
    a, b = R_I + D_SPHERE / 2, R_A - D_SPHERE / 2
    v_ann = math.pi * (R_A ** 2 - R_I ** 2)
    v_band = math.pi * (b ** 2 - a ** 2)
    phi_band = phi_nominal * v_ann / v_band
    eta_band = float(CLOSURES[closure](np.array(phi_band)))
    if weight_kind == "r3":
        f = lambda x, y: 0.5 * (x ** -2 - y ** -2)
    else:
        f = lambda x, y: y - x
    tot = f(R_I, R_A)
    walls = f(R_I, a) + f(b, R_A)
    band = f(a, b)
    return phi_band, eta_band, tot / (walls + band / eta_band)


def run_one(snaps, closure, grid, samples, seed, weight_kind, verbose=True):
    phi, st = deposit(snaps, grid, samples, seed)
    eta_comp, eta_slice = composite(phi, st["r_edges"], closure, weight_kind)
    if verbose:
        print(f"snapshots {st['n_snapshots']}  spheres/snapshot {st['n_spheres_per_snapshot']:.0f}  "
              f"grid {grid}x{grid}  samples/sphere {samples}  seed {seed}  series-weight {weight_kind}")
        print(f"  centre radius range        : {st['r_center_min']:.4f} .. {st['r_center_max']:.4f}")
        print(f"  sample points outside annulus: {st['n_outside']} of {st['n_sample_points']}")
        print(f"  CONSERVATION deposited/expected volume = {st['conservation']:.6f}")
        print(f"  global phi (deposit)       : {st['phi_global']:.6f}   nominal N*V_s/V_annulus = {st['phi_nominal']:.6f}")
        print(f"  local phi max              : {st['phi_max']:.4f}")
        print(f"  per-slice eta range        : {eta_slice.min():.4f} .. {eta_slice.max():.4f}")
        print(f"  COMPOSITE ({closure}, series-{weight_kind} in r, mean over z): {eta_comp:.5f}")
    return eta_comp, st, eta_slice


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("positions", nargs="+")
    ap.add_argument("--closure", choices=CLOSURES, default="einstein")
    ap.add_argument("--tmin", type=float, required=True)
    ap.add_argument("--tmax", type=float, required=True)
    ap.add_argument("--dt-snap", type=float, default=1.0)
    ap.add_argument("--grid", type=int, default=20)
    ap.add_argument("--samples", type=int, default=2000, help="MC points per sphere")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--series-weight", choices=("r3", "uniform"), default="r3")
    ap.add_argument("--convergence", action="store_true")
    ap.add_argument("--anchors", action="store_true")
    ap.add_argument("--measured", type=float, default=None, help="measured eta to compare against")
    args = ap.parse_args()

    snaps = read_positions(args.positions, args.tmin, args.tmax, args.dt_snap)
    if not snaps:
        sys.exit("no snapshots found in the requested window")
    times = list(snaps)
    counts = sorted({len(v) for v in snaps.values()})
    print(f"files: {' '.join(args.positions)}")
    print(f"window t in [{args.tmin}, {args.tmax}] cadence {args.dt_snap}: {len(times)} snapshots "
          f"({times[0]:.3f} .. {times[-1]:.3f}), spheres per snapshot {counts}")

    eta_comp, st, _ = run_one(snaps, args.closure, args.grid, args.samples, args.seed, args.series_weight)

    if args.measured is not None:
        dev = args.measured / eta_comp - 1.0
        print(f"  measured eta {args.measured:.4f} vs composite {eta_comp:.5f}: deviation {100*dev:+.2f} percent")

    if args.anchors:
        naive = float(CLOSURES[args.closure](np.array(st["phi_nominal"])))
        print(f"anchors: naive {args.closure} at nominal phi {st['phi_nominal']:.4f}: {naive:.5f}")
        for wk in ("r3", "uniform"):
            pb, eb, ea = depletion_anchor(st["phi_nominal"], args.closure, wk)
            print(f"  wall-depletion-only ({wk}): phi_band {pb:.4f}, eta_band {eb:.5f} -> eta_eff {ea:.5f}")

    if args.convergence:
        print("\nconvergence table (composite; conservation ratio; local phi max)")
        print(f"{'grid':>6} {'samples':>8} {'weight':>8} {'composite':>10} {'conserv':>9} {'phi_max':>8} {'slice min':>10} {'slice max':>10}")
        for grid in (20, 40):
            for samples in (2000, 20000):
                for wk in ("r3", "uniform"):
                    e, s, es = run_one(snaps, args.closure, grid, samples, args.seed, wk, verbose=False)
                    print(f"{grid:>6} {samples:>8} {wk:>8} {e:>10.5f} {s['conservation']:>9.6f} {s['phi_max']:>8.4f} {es.min():>10.4f} {es.max():>10.4f}")


if __name__ == "__main__":
    main()
