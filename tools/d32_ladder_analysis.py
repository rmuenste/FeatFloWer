#!/usr/bin/env python3
"""D3.2 hindered-settling ladder analysis.

Parses particle_force.log (columns: time ip fx fy fz tx ty tz px py pz vx vy vz)
from walled-column swarm rundirs and extracts, per run:

  - plateau settling speed U = -<vz>_particles averaged over the settled window
  - U/u_t against the single-particle reference u_t
  - measured cloud volume fraction phi_cloud = N*Vp / (A_cs * H_cloud), with
    H_cloud = (z_max - z_min + d) from particle extents (full cross-section
    assumed occupied; walls make this a mild overestimate of cloud volume)
  - inter-particle velocity spread std(vz), window-averaged

Emits one table row per rundir plus a Richardson-Zaki fit
  ln(U/u_t) = n * ln(1 - phi)  ->  n = slope through origin
over all runs (and per-N seed scatter).

Usage:
  python3 tools/d32_ladder_analysis.py q2p1_dns_rundir_d32_n20_s1 [more dirs...] \
      [--window 15 25] [--ut 0.4061] [--area 36.0] [--radius 0.5]
"""
import argparse
import math
import os
import sys
from collections import defaultdict


def parse_run(path, t0, t1):
    """Return per-snapshot (t, mean_vz, std_vz, zmin, zmax, N) within [t0, t1]."""
    log = os.path.join(path, "particle_force.log")
    snaps = []  # (t, [vz...], [pz...])
    cur_t = None
    vz, pz = [], []
    with open(log) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 14:
                continue
            t = float(parts[0])
            if cur_t is None or abs(t - cur_t) > 1e-9:
                if cur_t is not None and t0 <= cur_t <= t1 and vz:
                    snaps.append((cur_t, vz, pz))
                cur_t, vz, pz = t, [], []
            vz.append(float(parts[13]))
            pz.append(float(parts[10]))
        if cur_t is not None and t0 <= cur_t <= t1 and vz:
            snaps.append((cur_t, vz, pz))
    return snaps


def analyze(path, t0, t1, ut, area, radius):
    snaps = parse_run(path, t0, t1)
    if not snaps:
        return None
    d = 2.0 * radius
    vp = 4.0 / 3.0 * math.pi * radius ** 3
    means, spreads, phis = [], [], []
    n_particles = len(snaps[0][1])
    for _, vz, pz in snaps:
        m = sum(vz) / len(vz)
        means.append(m)
        spreads.append(math.sqrt(sum((v - m) ** 2 for v in vz) / len(vz)))
        h = max(pz) - min(pz) + d
        phis.append(len(vz) * vp / (area * h))
    mu = sum(means) / len(means)
    sd = math.sqrt(sum((m - mu) ** 2 for m in means) / len(means))
    return {
        "run": os.path.basename(path.rstrip("/")),
        "N": n_particles,
        "n_snaps": len(snaps),
        "U": -mu,
        "U_std": sd,
        "U_ut": -mu / ut,
        "phi": sum(phis) / len(phis),
        "spread": sum(spreads) / len(spreads),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundirs", nargs="+")
    ap.add_argument("--window", nargs=2, type=float, default=[15.0, 25.0])
    ap.add_argument("--ut", type=float, default=0.4061)
    ap.add_argument("--area", type=float, default=36.0)
    ap.add_argument("--radius", type=float, default=0.5)
    args = ap.parse_args()

    rows = []
    for rd in args.rundirs:
        r = analyze(rd, args.window[0], args.window[1], args.ut, args.area, args.radius)
        if r is None:
            print(f"WARN: no snapshots in window for {rd}", file=sys.stderr)
            continue
        rows.append(r)

    hdr = f"{'run':34s} {'N':>4s} {'snaps':>5s} {'U':>8s} {'U_std':>7s} {'U/ut':>6s} {'phi_cloud':>9s} {'spread':>7s}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(f"{r['run']:34s} {r['N']:4d} {r['n_snaps']:5d} {r['U']:8.4f} "
              f"{r['U_std']:7.4f} {r['U_ut']:6.3f} {r['phi']:9.4f} {r['spread']:7.4f}")

    # RZ fit: ln(U/ut) = n ln(1-phi), least squares through origin
    pts = [(math.log(1.0 - r["phi"]), math.log(r["U_ut"])) for r in rows if r["U_ut"] > 0]
    if len(pts) >= 2:
        sxx = sum(x * x for x, _ in pts)
        sxy = sum(x * y for x, y in pts)
        n_fit = sxy / sxx
        resid = [y - n_fit * x for x, y in pts]
        rms = math.sqrt(sum(e * e for e in resid) / len(resid))
        print(f"\nRZ fit over {len(pts)} runs: n = {n_fit:.2f}  (rms ln-residual {rms:.3f})")
        by_n = defaultdict(list)
        for r in rows:
            by_n[r["N"]].append(r["U_ut"])
        for N in sorted(by_n):
            v = by_n[N]
            print(f"  N={N:4d}: U/ut = {', '.join(f'{x:.3f}' for x in v)}"
                  + (f"  (seed span {max(v)-min(v):.3f})" if len(v) > 1 else ""))


if __name__ == "__main__":
    main()
