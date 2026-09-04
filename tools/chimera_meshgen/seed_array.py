#!/usr/bin/env python3
"""Seed a fixed sphere array in a periodic box and emit the Chimera body table.

Design: chimera-integration-design.md v3, section 7 and the Phase-5
roadmap row ("H_k seeding").  Memory note `static-chimera-closures`:
per-particle atmosphere width

    H_k = min(H_max, 0.5 * nearest surface-to-surface gap of particle k)

computed once at seeding guarantees the paper's non-overlap assumption
(atmosphere k never intersects body j; pairwise-disjoint atmospheres for
Chimera-S) by construction.  Gaps are measured with the periodic
minimum-image convention, so a particle's own periodic images count as
neighbours (a single sphere in a unit cell sees its image at distance L).

Modes:
  sc      simple-cubic lattice of n^3 spheres (n = --n, default 1); the
          Hasimoto (1959) array.  --offset shifts the lattice (e.g. 0,0,0
          puts the single sphere at the box corner so that its atmosphere
          straddles all periodic faces -- the periodic-donor check).
  random  random sequential addition (hard-sphere, minimum surface gap
          --mingap diameters), same algorithm as tools/gen_random_array.py
          but WITHOUT explicit image spheres (the Chimera geometry is
          minimum-image periodic itself).

The sphere radius follows from --phi (solid volume fraction) and the
count, or is given directly with --radius.

Output (--out): the SimPar@ChimeraParticleFile format,
    <count>
    sphere  cx cy cz  R  H_k        (one line per body)
with a '#' header recording L, phi, R, seed, H_max and the minimum gap.

Usage:
  seed_array.py --mode sc --n 1 --phi 0.0193925 --box 1 --hmax 1.0 --out particles.dat
  seed_array.py --mode random --count 8 --phi 0.1 --seed 1 --mingap 0.05 --out particles.dat
"""

import argparse
import math
import random
import sys


def min_image_dist(p, q, L):
    d2 = 0.0
    for c in range(3):
        d = abs(p[c] - q[c])
        d = min(d, L - d)
        d2 += d * d
    return math.sqrt(d2)


def nearest_gap(k, pts, r, L):
    """Nearest surface-to-surface gap of particle k (own images included)."""
    gap = L - 2.0 * r  # own periodic image at distance L
    for j, q in enumerate(pts):
        if j == k:
            continue
        gap = min(gap, min_image_dist(pts[k], q, L) - 2.0 * r)
    return gap


def seed_sc(n, L, offset):
    pts = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                pts.append((((i + 0.5) * L / n + offset[0]) % L,
                            ((j + 0.5) * L / n + offset[1]) % L,
                            ((k + 0.5) * L / n + offset[2]) % L))
    return pts


def seed_random(count, r, L, seed, mingap, maxtries):
    dmin = 2.0 * r * (1.0 + mingap)
    rng = random.Random(seed)
    pts = []
    tries = 0
    while len(pts) < count:
        tries += 1
        if tries > maxtries:
            sys.exit("FAILED: placed %d/%d after %d tries (phi too high for RSA?)"
                     % (len(pts), count, tries))
        p = (rng.random() * L, rng.random() * L, rng.random() * L)
        if all(min_image_dist(p, q, L) >= dmin for q in pts):
            pts.append(p)
    return pts


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mode", choices=("sc", "random"), required=True)
    ap.add_argument("--box", type=float, default=1.0, help="periodic box length L")
    ap.add_argument("--n", type=int, default=1, help="sc: spheres per axis")
    ap.add_argument("--count", type=int, default=1, help="random: number of spheres")
    ap.add_argument("--phi", type=float, help="solid volume fraction")
    ap.add_argument("--radius", type=float, help="sphere radius (overrides --phi)")
    ap.add_argument("--offset", default=None,
                    help="sc: lattice shift 'dx,dy,dz' (default: centred lattice)")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--mingap", type=float, default=0.05,
                    help="random: minimum surface gap in diameters")
    ap.add_argument("--hmax", type=float, default=1.0,
                    help="H_max in units of the radius (H_max = hmax * R)")
    ap.add_argument("--maxtries", type=int, default=2000000)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    L = a.box
    count = a.n ** 3 if a.mode == "sc" else a.count
    if a.radius is not None:
        r = a.radius
        phi = count * 4.0 / 3.0 * math.pi * r ** 3 / L ** 3
    elif a.phi is not None:
        phi = a.phi
        r = (phi * 3.0 * L ** 3 / (4.0 * math.pi * count)) ** (1.0 / 3.0)
    else:
        sys.exit("need --phi or --radius")

    if a.mode == "sc":
        offset = (0.0, 0.0, 0.0)
        if a.offset is not None:
            offset = tuple(float(s) for s in a.offset.split(","))
            # shift such that the FIRST lattice site lands at the offset
            offset = tuple(offset[c] - 0.5 * L / a.n for c in range(3))
        pts = seed_sc(a.n, L, offset)
    else:
        pts = seed_random(count, r, L, a.seed, a.mingap, a.maxtries)

    hmax = a.hmax * r
    gaps = [nearest_gap(k, pts, r, L) for k in range(len(pts))]
    if min(gaps) <= 0.0:
        sys.exit("FAILED: overlapping spheres (min gap %g)" % min(gaps))
    hk = [min(hmax, 0.5 * g) for g in gaps]

    with open(a.out, "w") as f:
        f.write("# Chimera body table (SimPar@ChimeraParticleFile), seed_array.py\n")
        f.write("# mode=%s L=%.16g count=%d phi=%.10g R=%.16g seed=%d "
                "Hmax=%.16g mingap_measured=%.10g Hk_min=%.10g Hk_max=%.10g\n"
                % (a.mode, L, len(pts), phi, r, a.seed, hmax, min(gaps),
                   min(hk), max(hk)))
        f.write("# H_k = min(H_max, 0.5 * nearest surface gap), minimum-image periodic\n")
        f.write("%d\n" % len(pts))
        for p, h in zip(pts, hk):
            f.write("sphere  %.16g %.16g %.16g   %.16g  %.16g\n"
                    % (p[0], p[1], p[2], r, h))
    print("wrote %s: %d spheres, R=%.6g phi=%.6g min gap=%.4g (%.3g d), "
          "H_k in [%.4g, %.4g]" % (a.out, len(pts), r, phi, min(gaps),
                                    min(gaps) / (2 * r), min(hk), max(hk)))


if __name__ == "__main__":
    main(sys.argv[1:])
