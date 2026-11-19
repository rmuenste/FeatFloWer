#!/usr/bin/env python3
"""Initial particle bed of case 11 (laterally periodic box, particles above the floor), written
as particles.xyz (x y z in metres, the file read by pe through example.json xyzFilePath_).

Two modes:
  rsa      loose random sequential addition, centres uniform in [0, LX) x [0, LY) x [r, H0 - r]
           with periodic overlap checks in x and y and a minimum surface gap GAP. The algorithm
           and the order of the random calls are those of the peer's make_bed.py (case_specs/
           10_batch_settling/reference, PERIODIC=1), so a seed gives the same bed as the
           peer's CFDEM run of case 11 (run_fluid.sh: FL_PHI0 0.25 over FL_H0 0.06 m, seed 1).
  lattice  simple cubic lattice of spacing s = (1 + GAP) d (solid fraction pi/6 (1 + GAP)^-3),
           layers from z = r + GAP d upwards, each centre jittered by up to JIT (fraction of the
           free gap) in x, y, z; a dense start for fixed-bed probes. N is reached layer by
           layer; the last layer is filled in random order.

N: given directly (--n) or from the bed of the spec: the number that fills the cross-section
to the static height HS at porosity EPS (--hs 0.025 --eps 0.40: 9 789 for bed R, 27 190 for
bed A). The default is the spec's static bed.

usage: make_bed_s11.py OUT.xyz [--mode rsa|lattice] [--lx 0.015] [--ly 0.015] [--d 0.87e-3]
         [--n N | --hs 0.025 --eps 0.40] [--h0 0.06] [--gap 0.02] [--jit 0.3] [--seed 1]"""
import argparse
import math
import random


def rsa(n, lx, ly, h0, d, gap, seed):
    # verbatim algorithm of make_bed.py (PERIODIC=1); random.seed(seed) before the first call
    random.seed(seed)
    r = d / 2
    dmin2 = (d + gap * d) ** 2
    h = d + gap * d
    grid = {}
    pts = []
    nx = int(lx / h)
    ny = int(ly / h)
    hx = lx / nx
    hy = ly / ny
    tries = 0
    while len(pts) < n:
        tries += 1
        p = (random.uniform(0, lx), random.uniform(0, ly), random.uniform(r, h0 - r))
        c = (min(int(p[0] / hx), nx - 1), min(int(p[1] / hy), ny - 1), int(p[2] / h))
        ok = True
        for i in (-1, 0, 1):
            for j in (-1, 0, 1):
                for k in (-1, 0, 1):
                    ci, cj = (c[0] + i) % nx, (c[1] + j) % ny
                    for q in grid.get((ci, cj, c[2] + k), ()):
                        dx = p[0] - q[0]
                        dy = p[1] - q[1]
                        dx -= lx * round(dx / lx)
                        dy -= ly * round(dy / ly)
                        if dx * dx + dy * dy + (p[2] - q[2]) ** 2 < dmin2:
                            ok = False
        if ok:
            grid.setdefault(c, []).append(p)
            pts.append(p)
        if tries > 200 * n:
            raise SystemExit('rsa: too many tries; lower N or raise H0')
    return pts


def lattice(n, lx, ly, d, gap, jit, seed):
    random.seed(seed)
    r = d / 2
    nx = int(lx / ((1 + gap) * d))
    ny = int(ly / ((1 + gap) * d))
    sx, sy = lx / nx, ly / ny            # periodic: the lattice closes over the period
    sz = (1 + gap) * d
    free_x, free_y, free_z = 0.5 * (sx - d), 0.5 * (sy - d), 0.5 * (sz - d)
    pts = []
    layer = 0
    while len(pts) < n:
        z0 = r + gap * d + layer * sz
        cells = [(i, j) for j in range(ny) for i in range(nx)]
        if n - len(pts) < len(cells):
            random.shuffle(cells)
            cells = cells[:n - len(pts)]
        for i, j in cells:
            pts.append(((i + 0.5) * sx + jit * free_x * random.uniform(-1, 1),
                        (j + 0.5) * sy + jit * free_y * random.uniform(-1, 1),
                        z0 + jit * free_z * random.uniform(-1, 1)))
        layer += 1
    return [((x % lx), (y % ly), z) for x, y, z in pts], nx, ny, sz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--mode', choices=('rsa', 'lattice'), default='rsa')
    ap.add_argument('--lx', type=float, default=0.015)
    ap.add_argument('--ly', type=float, default=0.015)
    ap.add_argument('--d', type=float, default=0.87e-3)
    ap.add_argument('--n', type=int, default=None)
    ap.add_argument('--hs', type=float, default=0.025)
    ap.add_argument('--eps', type=float, default=0.40)
    ap.add_argument('--h0', type=float, default=0.06)
    ap.add_argument('--gap', type=float, default=0.02)
    ap.add_argument('--jit', type=float, default=0.3)
    ap.add_argument('--seed', type=int, default=1)
    a = ap.parse_args()
    vp = math.pi * a.d ** 3 / 6
    n = a.n if a.n is not None else round((1 - a.eps) * a.lx * a.ly * a.hs / vp)
    if a.mode == 'rsa':
        pts = rsa(n, a.lx, a.ly, a.h0, a.d, a.gap, a.seed)
        info = f'rsa over H0 {a.h0:g} m: solid fraction {n * vp / (a.lx * a.ly * a.h0):.4f}'
    else:
        pts, nx, ny, sz = lattice(n, a.lx, a.ly, a.d, a.gap, a.jit, a.seed)
        top = max(p[2] for p in pts) + a.d / 2
        info = (f'lattice {nx} x {ny} per layer, spacing {sz:.4e} m, top {top * 1e3:.2f} mm, '
                f'solid fraction below the top {n * vp / (a.lx * a.ly * top):.4f}')
    with open(a.out, 'w') as f:
        for p in pts:
            f.write('%.9e %.9e %.9e\n' % p)
    print(f'{a.out}: N {n}, d {a.d:g} m, box {a.lx:g} x {a.ly:g} m; {info}; '
          f'static height at eps {a.eps:g}: {n * vp / ((1 - a.eps) * a.lx * a.ly) * 1e3:.2f} mm')


if __name__ == '__main__':
    main()
