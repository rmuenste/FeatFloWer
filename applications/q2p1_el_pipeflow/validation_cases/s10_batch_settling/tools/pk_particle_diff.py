#!/usr/bin/env python3
"""Largest particle-state difference between two runs of the same bed (el_particles_rank*.log).
usage: pk_particle_diff.py RUN_A RUN_B T [SX SY]      (env: L period, default 0.02)
Particles are matched by id at logged time T. Positions of RUN_B are moved back by (SX, SY) and
compared modulo the period in x and y; prints the largest and rms difference of position and
velocity, and the rms particle speed for scale."""
import os, sys, glob, math
L = float(os.environ.get('L', 0.02))
def load(run, t):
    P = {}
    for f in glob.glob(os.path.join(run, 'el_particles_rank*.log')):
        with open(f, errors='replace') as fh:
            for l in fh:
                if l.startswith('#'): continue
                w = l.split()
                if len(w) >= 14 and abs(float(w[0]) - t) < 1e-9: P[int(w[1])] = [float(x) for x in w[8:14]]
    return P
t = float(sys.argv[3]); sx = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0; sy = float(sys.argv[5]) if len(sys.argv) > 5 else 0.0
A = load(sys.argv[1], t); B = load(sys.argv[2], t); ids = sorted(set(A) & set(B))
dp = []; dv = []; sp = 0.0
for i in ids:
    a, b = A[i], B[i]
    dx = [(b[0] - sx - a[0] + L / 2) % L - L / 2, (b[1] - sy - a[1] + L / 2) % L - L / 2, b[2] - a[2]]
    dp.append(math.sqrt(sum(x * x for x in dx))); dv.append(math.sqrt(sum((b[k] - a[k])**2 for k in (3, 4, 5))))
    sp += sum(a[k]**2 for k in (3, 4, 5))
n = len(ids); rms = lambda v: math.sqrt(sum(x * x for x in v) / n)
print('t %.4f  particles A %d B %d common %d' % (t, len(A), len(B), n))
print('position difference [m]: max %.3e rms %.3e' % (max(dp), rms(dp)))
print('velocity difference [m/s]: max %.3e rms %.3e   (rms particle speed %.3e)' % (max(dv), rms(dv), math.sqrt(sp / n)))
