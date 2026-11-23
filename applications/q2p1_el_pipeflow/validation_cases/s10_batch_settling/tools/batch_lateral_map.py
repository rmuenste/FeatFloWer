#!/usr/bin/env python3
"""Lateral maps for the batch-settling case from FeatFloWer particle logs.
usage: batch_lateral_map.py RUNDIR T1 [T2 ...]     (env: LX LY NX ZLO ZHI)
For each requested time: NX x NX map (rows = y index from 0, columns = x index) of the column solid
fraction (centre counts * V_p / column volume) and of the column-mean v_z over ZLO..ZHI (default
40-100 mm), their mean and rms over the columns, the seam columns (index 0 and NX-1 in x or y) against
the others, the position range, and the number of particles within r of a periodic plane.
Also prints the Poisson expectation of the rms of the column solid fraction."""
import os, sys, glob, math, collections
def env(k, v): return float(os.environ.get(k, v))
LX, LY, NX = env('LX', 0.02), env('LY', 0.02), int(env('NX', 6))
ZLO, ZHI = env('ZLO', 0.04), env('ZHI', 0.10); d = 1e-3; vp = math.pi * d**3 / 6
run = sys.argv[1]; want = [float(x) for x in sys.argv[2:]]
frames = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log'))):
    with open(f, errors='replace') as fh:
        for l in fh:
            if l.startswith('#'): continue
            w = l.split()
            if len(w) < 14: continue
            t = round(float(w[0]), 6)
            if any(abs(t - x) < 1e-6 for x in want): frames[t].append((int(w[1]),) + tuple(float(x) for x in w[8:14]))
dx = LX / NX; vcol = dx * dx * (ZHI - ZLO)
def stats(vals):
    m = sum(vals) / len(vals); return m, math.sqrt(sum((v - m)**2 for v in vals) / len(vals))
for t in sorted(frames):
    P = frames[t]; ids = [p[0] for p in P]
    xs = [p[1] for p in P]; ys = [p[2] for p in P]; zs = [p[3] for p in P]
    print('=== t %.4f  N %d  unique ids %d  x [%.6e, %.6e]  y [%.6e, %.6e]  z [%.6e, %.6e]' % (
        t, len(P), len(set(ids)), min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)))
    print('    outside [0,L): x %d  y %d ; within r of x planes %d, of y planes %d ; max|vx| %.5f max|vy| %.5f mean vz all %.6f' % (
        sum(1 for x in xs if x < 0 or x >= LX), sum(1 for y in ys if y < 0 or y >= LY),
        sum(1 for x in xs if x < d / 2 or x > LX - d / 2), sum(1 for y in ys if y < d / 2 or y > LY - d / 2),
        max(abs(p[4]) for p in P), max(abs(p[5]) for p in P), sum(p[6] for p in P) / len(P)))
    n = [[0] * NX for _ in range(NX)]; s = [[0.0] * NX for _ in range(NX)]
    for _, x, y, z, vx, vy, vz in P:
        if ZLO <= z < ZHI:
            i = min(max(int(x / dx), 0), NX - 1); j = min(max(int(y / dx), 0), NX - 1); n[j][i] += 1; s[j][i] += vz
    phi = [[n[j][i] * vp / vcol for i in range(NX)] for j in range(NX)]
    vz = [[s[j][i] / n[j][i] if n[j][i] else float('nan') for i in range(NX)] for j in range(NX)]
    print('  column solid fraction (rows y0..y%d, cols x0..x%d):' % (NX - 1, NX - 1))
    for j in range(NX): print('   ' + ' '.join('%.4f' % v for v in phi[j]))
    print('  column-mean v_z (mm/s):')
    for j in range(NX): print('   ' + ' '.join('%+7.3f' % (1e3 * v) for v in vz[j]))
    fp = [v for r in phi for v in r]; fv = [v for r in vz for v in r]
    mp, rp = stats(fp); mv, rv = stats(fv); ntot = sum(sum(r) for r in n)
    print('  phi: mean %.4f rms %.4f (Poisson rms %.4f)   v_z: mean of columns %.3f mm/s  rms %.3f mm/s  min %.3f max %.3f' % (
        mp, rp, mp / math.sqrt(ntot / NX**2), 1e3 * mv, 1e3 * rv, 1e3 * min(fv), 1e3 * max(fv)))
    def grp(sel):
        nn = sum(n[j][i] for j in range(NX) for i in range(NX) if sel(i, j)); ss = sum(s[j][i] for j in range(NX) for i in range(NX) if sel(i, j))
        k = sum(1 for j in range(NX) for i in range(NX) if sel(i, j)); return nn * vp / (k * vcol), 1e3 * ss / max(nn, 1), k
    e = (0, NX - 1)
    for name, sel in (('x-seam columns (i=0,%d)' % (NX - 1), lambda i, j: i in e), ('non-x-seam', lambda i, j: i not in e),
                      ('y-seam columns (j=0,%d)' % (NX - 1), lambda i, j: j in e), ('non-y-seam', lambda i, j: j not in e),
                      ('any seam ring', lambda i, j: i in e or j in e), ('interior', lambda i, j: i not in e and j not in e)):
        a, b, k = grp(sel); print('  %-26s phi %.4f  particle-mean v_z %+.3f mm/s  (%d columns)' % (name, a, b, k))
    # profiles along x and y (column sums)
    px = [sum(n[j][i] for j in range(NX)) * vp / (NX * vcol) for i in range(NX)]; py = [sum(n[j][i] for i in range(NX)) * vp / (NX * vcol) for j in range(NX)]
    vx_ = [1e3 * sum(s[j][i] for j in range(NX)) / max(sum(n[j][i] for j in range(NX)), 1) for i in range(NX)]
    vy_ = [1e3 * sum(s[j][i] for i in range(NX)) / max(sum(n[j][i] for i in range(NX)), 1) for j in range(NX)]
    print('  phi(x): ' + ' '.join('%.4f' % v for v in px) + '   phi(y): ' + ' '.join('%.4f' % v for v in py))
    print('  vz(x):  ' + ' '.join('%+.3f' % v for v in vx_) + '   vz(y): ' + ' '.join('%+.3f' % v for v in vy_))
