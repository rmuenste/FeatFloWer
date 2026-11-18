#!/usr/bin/env python3
"""Per-step seam report of a laterally periodic batch-settling probe run (30 steps, VTK and
particle log every step).
usage: pk_probe_report.py RUNDIR [RUNDIR ...]        (env: L period, default 0.02)
Per step (from run_slurm.log): initial defect of the first pressure solve (DefInitPres of the
first 'XxX' line of the step), number of multigrid cycles of that solve, maximum norm of the
fluid velocity (last value printed in the step).
Per VTK frame (frame index = step - 1): extrema of the cell field ELDepsFDt over all ranks, and
the centre (mm) of the cells that hold them.
Crossings (el_particles_rank*.log): particle and step whenever a coordinate of a particle jumps by
more than L/2 between consecutive logged steps.
Several run directories are printed side by side in the order given."""
import os, sys, re, glob, collections
L = float(os.environ.get('L', 0.02))

def log_steps(run):
    out = collections.OrderedDict(); step = None; vmax = None; first = None
    with open(os.path.join(run, 'run_slurm.log'), errors='replace') as f:
        for l in f:
            if l.startswith(' maximum norm of velocity:'):
                vmax = float(l.split(':')[1])
            elif l.startswith('XxX'):
                w = l.split(); st = int(w[-1])
                if st not in out:
                    out[st] = [float(w[2].replace('D', 'E')), int(w[1]), None]
                step = st
            elif l.startswith('time:') and step is not None:
                out[step][2] = vmax
    if step is not None and out[step][2] is None: out[step][2] = vmax
    return out

def deps_extrema(run):
    res = {}
    frames = sorted({int(f.split('.')[-2]) for f in glob.glob(os.path.join(run, '_vtk', 'res_node_0001.*.vtu'))})
    for fr in frames:
        lo = (1e99, None); hi = (-1e99, None)
        for fn in glob.glob(os.path.join(run, '_vtk', 'res_node_*.%05d.vtu' % fr)):
            txt = open(fn).read(); a = {}
            for m in re.finditer(r'<DataArray[^>]*Name="(\w+)"[^>]*>(.*?)</DataArray>', txt, re.S):
                if m.group(1) in ('Points', 'ELDepsFDt', 'connectivity'): a[m.group(1)] = m.group(2).split()
            if 'ELDepsFDt' not in a: continue
            P = [float(x) for x in a['Points']]; C = [int(x) for x in a['connectivity']]
            D = [float(x) for x in a['ELDepsFDt']]
            for e, v in enumerate(D):
                if v < lo[0] or v > hi[0]:
                    c = tuple(sum(P[3 * i + k] for i in C[8 * e:8 * e + 8]) / 8 * 1e3 for k in range(3))
                    if v < lo[0]: lo = (v, c)
                    if v > hi[0]: hi = (v, c)
        res[fr + 1] = (lo, hi)
    return res

def crossings(run):
    pos = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(run, 'el_particles_rank*.log')):
        with open(f, errors='replace') as fh:
            for l in fh:
                if l.startswith('#'): continue
                w = l.split()
                if len(w) < 14: continue
                pos[int(w[1])][round(float(w[0]), 6)] = (float(w[8]), float(w[9]))
    ev = []
    for ip, d in pos.items():
        ts = sorted(d)
        for a, b in zip(ts, ts[1:]):
            for k, ax in enumerate('xy'):
                if abs(d[b][k] - d[a][k]) > 0.5 * L: ev.append((b, ip, ax))
    return sorted(ev)

for run in sys.argv[1:]:
    S = log_steps(run); E = deps_extrema(run); X = crossings(run)
    print('=== %s' % run)
    print('# step  DefInitPres  nMG  max|u| [m/s] | ELDepsFDt min [1/s] at (x,y,z mm) | max at (x,y,z mm)')
    for st in S:
        e = E.get(st)
        es = '' if e is None else ' | %10.3e (%5.2f %5.2f %6.2f) | %10.3e (%5.2f %5.2f %6.2f)' % ((e[0][0],) + e[0][1] + (e[1][0],) + e[1][1])
        print('%4d  %10.3e  %3d  %10.4e%s' % (st, S[st][0], S[st][1], S[st][2] if S[st][2] is not None else float('nan'), es))
    dt = 5e-4
    print('# crossings (step = time/dt, particle, axis): ' + ', '.join('%d:%d:%s' % (round(t / dt), ip, ax) for t, ip, ax in X))
    print('# number of crossings %d, distinct particles %d' % (len(X), len({ip for _, ip, _ in X})))
