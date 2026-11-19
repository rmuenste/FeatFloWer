#!/usr/bin/env python3
"""Lateral maps of nodal fields from the FeatFloWer VTK output (ascii .vtu, one piece per rank).
usage: batch_vtk_lateral.py RUNDIR FRAME [FRAME ...]     (env: ZLO ZHI, default 40-100 mm)
For each frame index (res_node_RRRR.FFFFF.vtu): maps over the (x, y) node lines of the z-average
(nodes with ZLO <= z <= ZHI) of the solid fraction 1 - ELEpsilonF, ELAlphaP and the fluid u_z;
duplicates of a node on several ranks are averaged (their spread is reported).
Rows = y node index from 0, columns = x node index; first and last row/column are the periodic planes."""
import os, sys, glob, re, collections, math
ZLO = float(os.environ.get('ZLO', 0.04)); ZHI = float(os.environ.get('ZHI', 0.10))
def arrays(fn):
    txt = open(fn).read(); out = {}
    for m in re.finditer(r'<DataArray[^>]*Name="(\w+)"[^>]*>(.*?)</DataArray>', txt, re.S):
        if m.group(1) in ('Points', 'ELEpsilonF', 'ELAlphaP', 'Velocity', 'connectivity'): out[m.group(1)] = [float(x) for x in m.group(2).split()]
    return out
run = sys.argv[1]
for fr in sys.argv[2:]:
    files = sorted(glob.glob(os.path.join(run, '_vtk', 'res_node_*.%05d.vtu' % int(fr))))
    node = collections.defaultdict(list); cell = collections.defaultdict(list)
    for fn in files:
        a = arrays(fn); P = a['Points']; n = len(P) // 3; C = [int(x) for x in a['connectivity']]
        for i in range(n):
            key = (round(P[3 * i] * 1e6), round(P[3 * i + 1] * 1e6), round(P[3 * i + 2] * 1e6))
            node[key].append((a['Velocity'][3 * i + 2],))
        for e in range(len(C) // 8):
            v = C[8 * e:8 * e + 8]; key = tuple(round(sum(P[3 * i + k] for i in v) / 8 * 1e6) for k in range(3))
            cell[key].append((a['ELEpsilonF'][e], a['ELAlphaP'][e]))
    def maps(D, q, tr):
        xs = sorted({k[0] for k in D}); ys = sorted({k[1] for k in D}); col = collections.defaultdict(lambda: [0, 0.0])
        for k, L in D.items():
            if ZLO * 1e6 - 1 <= k[2] <= ZHI * 1e6 + 1:
                c = col[(k[0], k[1])]; c[0] += 1; c[1] += sum(v[q] for v in L) / len(L)
        return [[tr(col[(x, y)][1] / col[(x, y)][0]) for x in xs] for y in ys]
    print('=== frame %s  files %d  nodes %d  cells %d (max rank copies of a cell %d)' % (fr, len(files), len(node), len(cell), max(len(L) for L in cell.values())))
    for name, M, per in (('solid fraction 1 - ELEpsilonF (element columns)', maps(cell, 0, lambda v: 1 - v), False),
                         ('ELAlphaP (element columns)', maps(cell, 1, lambda v: v), False),
                         ('fluid u_z, mm/s (node lines; first and last row/column are the periodic planes)', maps(node, 0, lambda v: 1e3 * v), True)):
        print('  ' + name + ':')
        for row in M: print('   ' + ' '.join('%8.4f' % v for v in row))
        ny = len(M); nx = len(M[0]); flat = [v for r in M for v in r]; m = sum(flat) / len(flat)
        seam = [M[j][i] for j in range(ny) for i in range(nx) if i in (0, nx - 1) or j in (0, ny - 1)]
        inner = [M[j][i] for j in range(ny) for i in range(nx) if not (i in (0, nx - 1) or j in (0, ny - 1))]
        print('   mean %.5f rms %.5f | outer ring (at the periodic planes) mean %.5f | others mean %.5f' % (
            m, math.sqrt(sum((v - m)**2 for v in flat) / len(flat)), sum(seam) / len(seam), sum(inner) / len(inner))
            + ('' if not per else ' | max |f(x=0)-f(x=L)| %.3e  max |f(y=0)-f(y=L)| %.3e' % (
            max(abs(M[j][0] - M[j][-1]) for j in range(ny)), max(abs(M[0][i] - M[-1][i]) for i in range(nx)))))
