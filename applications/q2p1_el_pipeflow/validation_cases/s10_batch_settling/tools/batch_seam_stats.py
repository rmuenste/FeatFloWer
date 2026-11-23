#!/usr/bin/env python3
"""Seam statistics for a laterally periodic batch-settling run.
usage: batch_seam_stats.py RUNDIR T1 [T2 ...]   (env: LX LY ZLO ZHI)
Particles (el_particles_rank*.log) in ZLO..ZHI, binned by the distance of the centre to the nearest
periodic plane of the axis (0-0.5, 0.5-1, 1-2, 2-3.33, >3.33 mm): count, expected count for a uniform
bed, rms of the velocity component NORMAL to that plane, mean v_z.
Fluid (VTK frame nearest in the series res_node_*.NNNNN.vtu, frame index = step-1): rms of the nodal
velocity component normal to the plane on the plane nodes and on the other nodes, and correlation of
column v_z with the column solid fraction (6 x 6 element columns)."""
import os, sys, glob, math, collections, re
def env(k, v): return float(os.environ.get(k, v))
LX, LY, ZLO, ZHI = env('LX', 0.02), env('LY', 0.02), env('ZLO', 0.04), env('ZHI', 0.10)
run = sys.argv[1]; want = [float(x) for x in sys.argv[2:]]
frames = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log'))):
    with open(f, errors='replace') as fh:
        for l in fh:
            if l.startswith('#'): continue
            w = l.split()
            if len(w) < 14: continue
            t = round(float(w[0]), 6)
            if any(abs(t - x) < 1e-6 for x in want): frames[t].append(tuple(float(x) for x in w[8:14]))
edges = [0, 0.5e-3, 1e-3, 2e-3, 3.3333e-3, 0.0101]
vtk = sorted({int(f.split('.')[-2]) for f in glob.glob(os.path.join(run, '_vtk', 'res_node_0001.*.vtu'))})
def arr(fn):
    txt = open(fn).read(); out = {}
    for m in re.finditer(r'<DataArray[^>]*Name="(\w+)"[^>]*>(.*?)</DataArray>', txt, re.S):
        if m.group(1) in ('Points', 'Velocity'): out[m.group(1)] = [float(x) for x in m.group(2).split()]
    return out
for t in sorted(frames):
    P = [p for p in frames[t] if ZLO <= p[2] < ZHI]
    print('=== t %.4f  particles in %.0f-%.0f mm: %d' % (t, ZLO * 1e3, ZHI * 1e3, len(P)))
    for ax, L, name in ((0, LX, 'x'), (1, LY, 'y')):
        for b in range(len(edges) - 1):
            S = [p for p in P if edges[b] <= min(p[ax], L - p[ax]) < edges[b + 1]]
            if not S: continue
            exp = len(P) * 2 * (edges[b + 1] - edges[b]) / L
            print('  %s planes, distance %.2f-%.2f mm: N %5d (uniform %.0f, ratio %.3f)  rms v_%s %.3f mm/s  mean v_z %+.3f mm/s' % (
                name, edges[b] * 1e3, edges[b + 1] * 1e3, len(S), exp, len(S) / exp, name,
                1e3 * math.sqrt(sum(p[3 + ax]**2 for p in S) / len(S)), 1e3 * sum(p[5] for p in S) / len(S)))
    # column correlation
    n = [[0] * 6 for _ in range(6)]; s = [[0.0] * 6 for _ in range(6)]
    for p in P:
        i = min(int(p[0] / (LX / 6)), 5); j = min(int(p[1] / (LY / 6)), 5); n[j][i] += 1; s[j][i] += p[5]
    a = [n[j][i] for j in range(6) for i in range(6)]; b = [s[j][i] / n[j][i] for j in range(6) for i in range(6)]
    ma = sum(a) / 36; mb = sum(b) / 36
    cov = sum((x - ma) * (y - mb) for x, y in zip(a, b)); va = sum((x - ma)**2 for x in a); vb = sum((y - mb)**2 for y in b)
    print('  correlation of column-mean v_z with column count: %.3f' % (cov / math.sqrt(va * vb)))
    fr = int(round(t / 5e-4)) - 1
    if fr in vtk:
        acc = {}
        for fn in glob.glob(os.path.join(run, '_vtk', 'res_node_*.%05d.vtu' % fr)):
            a_ = arr(fn); Q = a_['Points']; V = a_['Velocity']
            for i in range(len(Q) // 3):
                if not (ZLO - 1e-9 <= Q[3 * i + 2] <= ZHI + 1e-9): continue
                for ax, L in ((0, LX), (1, LY)):
                    on = min(Q[3 * i + ax], L - Q[3 * i + ax]) < 1e-7; q = acc.setdefault((ax, on), [0, 0.0, 0.0]); q[0] += 1; q[1] += V[3 * i + ax]**2; q[2] += V[3 * i + 2]**2
        for ax in (0, 1):
            print('  fluid (VTK frame %d) rms u_%s: on the %s planes %.3f mm/s, other nodes %.3f mm/s ; rms u_z on the planes %.3f, other nodes %.3f' % (
                fr, 'xy'[ax], 'xy'[ax], 1e3 * math.sqrt(acc[(ax, True)][1] / acc[(ax, True)][0]), 1e3 * math.sqrt(acc[(ax, False)][1] / acc[(ax, False)][0]),
                1e3 * math.sqrt(acc[(ax, True)][2] / acc[(ax, True)][0]), 1e3 * math.sqrt(acc[(ax, False)][2] / acc[(ax, False)][0])))
