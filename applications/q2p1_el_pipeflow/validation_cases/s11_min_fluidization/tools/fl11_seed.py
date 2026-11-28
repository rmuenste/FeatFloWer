#!/usr/bin/env python3
"""Particle positions of a q2p1_el_pipeflow run at a logged time, as a pe seed file (case 11).

Reads the el_particles_rank*.log files of RUNDIR (columns: time ip fx fy fz tx ty tz px py pz vx vy
vz; each rank writes the particles it owns), takes the last logged time not later than --time,
collects the positions of that time from all files and writes them, sorted by the particle index
ip, as "x y z" lines (the tokens of the log, unchanged) to OUT, the format of pe's xyzFilePath_.

Refuses (exit 1, nothing written) when the count at that time differs from --count or a particle
index appears twice. Reports the count, the z range of the centres, the smallest centre gap
|x_i - x_j| - 2 r over all pairs (minimum image in the periodic x and y directions of period
--lx / --ly) and the smallest wall gaps z - r - zmin and zmax - z - r. A negative gap is an
overlap: pe's seeding then needs "seedAllowContact_": true in the json.

usage: fl11_seed.py RUNDIR --time T --out particles.xyz [--count 9789] [--radius 4.35e-4]
                    [--lx 0.015 --ly 0.015 --zmin 0 --zmax 0.1]
       fl11_seed.py --self-test
"""
import argparse
import glob
import os
import sys
import tempfile

import numpy as np


def log_files(rundir):
    files = sorted(glob.glob(os.path.join(rundir, 'el_particles_rank*.log')))
    if not files:
        raise SystemExit(f'no el_particles_rank*.log in {rundir}')
    return files


def pick_time(files, tmax):
    """The last logged time (as printed) not later than tmax over all files (streamed: the logs
    of a protocol run are gigabytes)."""
    best = None
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                tok = line.split(None, 1)[0]
                t = float(tok)
                if t <= tmax + 1e-9 and (best is None or t > best[0]):
                    best = (t, tok)
    if best is None:
        raise SystemExit(f'no logged time <= {tmax}')
    return best[1]


def read_time(files, tk):
    """[(ip, x, y, z, (xs, ys, zs), velocity)] of the lines with time token tk in all files."""
    rows = []
    for f in files:
        with open(f) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                t = line.split()
                if t[0] != tk or len(t) < 14:
                    continue
                rows.append((int(t[1]), float(t[8]), float(t[9]), float(t[10]), (t[8], t[9], t[10]),
                             (float(t[11]), float(t[12]), float(t[13]))))
    return rows


def min_image_gaps(pos, i, j, r, lx, ly):
    """|x_i - x_j| - 2 r for the indices j (minimum image in x, y)"""
    d = pos[j] - pos[i]
    d[:, 0] -= lx*np.round(d[:, 0]/lx)
    d[:, 1] -= ly*np.round(d[:, 1]/ly)
    return np.sqrt((d*d).sum(axis=1)) - 2*r


def min_pair_gap_all_pairs(pos, r, lx, ly):
    """the smallest gap over all pairs, one row of the distance matrix at a time"""
    best, pair = float('inf'), None
    for i in range(len(pos) - 1):
        j = np.arange(i + 1, len(pos))
        g = min_image_gaps(pos, i, j, r, lx, ly)
        m = int(np.argmin(g))
        if g[m] < best:
            best, pair = float(g[m]), (i, int(j[m]))
    return best, pair


def min_pair_gap(pos, r, lx, ly):
    """Smallest |x_i - x_j| - 2 r (minimum image in x, y) and its pair. A cell list of width >= 4 r finds
    every pair closer than 4 r between centres (a bed: all that matters); when its smallest gap is not below
    2 r, a closer pair in cells that are not adjacent is not excluded and all pairs are searched."""
    n = len(pos)
    if n < 2:
        return float('inf'), None
    cut = 4.0*r
    nx, ny = max(1, int(lx // cut)), max(1, int(ly // cut))
    zmin = pos[:, 2].min()
    nz = max(1, int((pos[:, 2].max() - zmin) // cut) + 1)
    cx = np.floor(np.mod(pos[:, 0], lx)/lx*nx).astype(int) % nx
    cy = np.floor(np.mod(pos[:, 1], ly)/ly*ny).astype(int) % ny
    cz = np.minimum(((pos[:, 2] - zmin)//cut).astype(int), nz - 1)
    cells = {}
    for i in range(n):
        cells.setdefault((cx[i], cy[i], cz[i]), []).append(i)
    best, pair = float('inf'), None
    for (a, b, c), idx in cells.items():
        nb = []
        for da in (-1, 0, 1):
            for db in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    k = ((a + da) % nx, (b + db) % ny, c + dc)
                    nb.extend(cells.get(k, []))
        nb = np.unique(np.array(nb))
        for i in idx:
            j = nb[nb > i]
            if len(j) == 0:
                continue
            g = min_image_gaps(pos, i, j, r, lx, ly)
            m = int(np.argmin(g))
            if g[m] < best:
                best, pair = float(g[m]), (i, int(j[m]))
    if best < cut - 2*r:
        return best, pair
    return min_pair_gap_all_pairs(pos, r, lx, ly)


def extract(rundir, tmax, out, count, r, lx, ly, zmin, zmax, report=sys.stdout):
    files = log_files(rundir)
    tk = pick_time(files, tmax)
    rows = read_time(files, tk)
    ips = [row[0] for row in rows]
    dup = len(ips) - len(set(ips))
    print(f'{len(files)} rank files; last logged time <= {tmax}: {tk}; {len(rows)} particles', file=report)
    if dup:
        print(f'REFUSED: {dup} particle index(es) appear twice at t = {tk}', file=report)
        return 1
    if len(rows) != count:
        print(f'REFUSED: {len(rows)} particles at t = {tk}, expected {count}', file=report)
        return 1
    rows.sort(key=lambda row: row[0])
    pos = np.array([[row[1], row[2], row[3]] for row in rows])
    vel = np.array([row[5] for row in rows])
    gap, pair = min_pair_gap(pos, r, lx, ly)
    floor_gap = pos[:, 2].min() - r - zmin
    top_gap = zmax - pos[:, 2].max() - r
    speed = np.sqrt((vel*vel).sum(axis=1))
    print(f'z range of the centres {pos[:, 2].min():.6e} .. {pos[:, 2].max():.6e} m', file=report)
    print(f'smallest pair gap {gap:.4e} m ({100*gap/(2*r):.3f} % of d) between ip {rows[pair[0]][0]} '
          f'and {rows[pair[1]][0]}' if pair else 'no pair', file=report)
    print(f'smallest wall gap: floor {floor_gap:.4e} m ({100*floor_gap/(2*r):.3f} % of d), top {top_gap:.4e} m',
          file=report)
    print(f'speeds at that time: max {speed.max():.3e} m/s, rms {np.sqrt((speed**2).mean()):.3e} m/s '
          f'(not seeded: pe seeds one uniform velocity)', file=report)
    if gap < 0 or floor_gap < 0:
        print('overlaps present: the json needs "seedAllowContact_": true', file=report)
    with open(out, 'w') as f:
        for row in rows:
            f.write(' '.join(row[4]) + '\n')
    print(f'wrote {len(rows)} positions to {out}', file=report)
    return 0


def self_test():
    ok = True
    r = 0.4
    with tempfile.TemporaryDirectory() as d:
        # 2 x 2 x 2 lattice of spacing 1 in a periodic 2 x 2 box: pair gap 1 - 2 r = 0.2, the
        # periodic images at distance 1 as well; particle 8 sits low: floor gap z - r = 0.1
        pts = [(i, j, k) for k in range(2) for j in range(2) for i in range(2)]
        lines = {1: [], 2: []}
        for n, (i, j, k) in enumerate(pts, start=1):
            z = 1.0 + k if n != 8 else 1.0 + k
            lines[1 if n % 2 else 2].append((n, 0.5 + i, 0.5 + j, z))
        for rank, rows in lines.items():
            with open(os.path.join(d, f'el_particles_rank{rank:04d}.log'), 'w') as f:
                f.write('# time ip fx fy fz tx ty tz px py pz vx vy vz\n')
                for t in ('0.100000', '0.200000'):
                    for n, x, y, z in rows:
                        zz = z - 0.5 if (n == 1 and t == '0.200000') else z
                        f.write(f'{t:>14s} {n:7d}' + ' 0.0'*6 + f' {x:.8E} {y:.8E} {zz:.8E} 0.0 0.0 0.0\n')
        out = os.path.join(d, 'p.xyz')
        sink = open(os.devnull, 'w')
        ok &= extract(d, 0.15, out, 8, r, 2.0, 2.0, 0.0, 10.0, sink) == 0
        xyz = np.loadtxt(out)
        ok &= xyz.shape == (8, 3) and abs(xyz[0, 2] - 1.0) < 1e-12
        gap, _ = min_pair_gap(xyz, r, 2.0, 2.0)
        ok &= abs(gap - 0.2) < 1e-9
        # time 0.2: particle 1 moved down by 0.5: gap to its neighbour above is 1.5 - 0.8, floor gap 0.1
        ok &= extract(d, 0.25, out, 8, r, 2.0, 2.0, 0.0, 10.0, sink) == 0
        xyz = np.loadtxt(out)
        ok &= abs(xyz[0, 2] - 0.5) < 1e-12 and abs(xyz[:, 2].min() - r - 0.1) < 1e-9
        # refusals: wrong count, duplicate index
        ok &= extract(d, 0.15, out, 9, r, 2.0, 2.0, 0.0, 10.0, sink) == 1
        with open(os.path.join(d, 'el_particles_rank0003.log'), 'w') as f:
            f.write('      0.100000       3' + ' 0.0'*6 + ' 1.0 1.0 1.0 0.0 0.0 0.0\n')
        ok &= extract(d, 0.15, out, 8, r, 2.0, 2.0, 0.0, 10.0, sink) == 1
        # periodic minimum image: two particles across the x seam at distance 0.9
        pos = np.array([[0.05, 1.0, 5.0], [1.95 - 0.8, 1.0, 5.0], [1.95, 1.0, 5.0]])
        gap, pair = min_pair_gap(pos, r, 2.0, 2.0)
        ok &= abs(gap - (0.1 - 0.8)) < 1e-9 and pair == (0, 2)
        # a sparse configuration: the nearest pair lies in cells that are not adjacent
        gap, pair = min_pair_gap(np.array([[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [50.0, 40.0, 3.0]]), 0.1, 100.0, 100.0)
        ok &= abs(gap - 9.8) < 1e-9 and pair == (0, 1)
        # ... and across the periodic boundary
        gap, pair = min_pair_gap(np.array([[1.0, 0.0, 0.0], [95.0, 0.0, 0.0], [50.0, 40.0, 3.0]]), 0.1, 100.0, 100.0)
        ok &= abs(gap - 5.8) < 1e-9 and pair == (0, 1)
    print('fl11_seed self-test', 'PASSED' if ok else 'FAILED')
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('rundir', nargs='?')
    ap.add_argument('--time', type=float, help='take the last logged time not later than this [s]')
    ap.add_argument('--out', help='output seed file')
    ap.add_argument('--count', type=int, default=9789, help='expected particle count (default 9789)')
    ap.add_argument('--radius', type=float, default=4.35e-4, help='particle radius [m]')
    ap.add_argument('--lx', type=float, default=0.015)
    ap.add_argument('--ly', type=float, default=0.015)
    ap.add_argument('--zmin', type=float, default=0.0)
    ap.add_argument('--zmax', type=float, default=0.1)
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not (a.rundir and a.time is not None and a.out):
        ap.error('RUNDIR, --time and --out are needed')
    return extract(a.rundir, a.time, a.out, a.count, a.radius, a.lx, a.ly, a.zmin, a.zmax)


if __name__ == '__main__':
    sys.exit(main())
