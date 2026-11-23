#!/usr/bin/env python3
"""Health and structure diagnostics for the batch-settling case (spec 10) from FeatFloWer particle
logs (el_particles_rank*.log: time ip f(3) t(3) p(3) v(3)).
usage: batch_diag_ff.py RUNDIR [DT_PRINT]      (env: LX LY DZ DX ZLO ZHI)
Per printed frame:
  N, max |v|                         particle count and largest speed (burst / loss check)
  phi(z) in the core [ZLO, ZHI]      mean, std, counting-noise estimate, lag-1 autocorrelation
                                     (bins of DZ, default one fluid element = 2.5 mm; centre counts)
  lateral structure in the core      solid fraction and mean v_z in the wall ring of DX-wide columns
                                     against the interior, and rms over the columns of the column-mean
                                     v_z (lateral convection check)
The core is fixed in space (default 60-120 mm) so that both fronts stay outside it until late."""
import os, sys, glob, math, collections
def env(k, v): return float(os.environ.get(k, v))
LX, LY, DZ, DX = env('LX', 0.02), env('LY', 0.02), env('DZ', 2.5e-3), env('DX', 2.5e-3)
ZLO, ZHI = env('ZLO', 0.06), env('ZHI', 0.12); d = 1e-3; vp = math.pi * d**3 / 6
run = sys.argv[1]; dtp = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
frames = collections.defaultdict(list)
for f in sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log'))):
    with open(f, errors='replace') as fh:
        for l in fh:
            if l.startswith('#'): continue
            w = l.split()
            if len(w) < 14: continue
            frames[round(float(w[0]), 6)].append(tuple(float(x) for x in w[8:14]))
nz = int(round((ZHI - ZLO) / DZ)); nx = int(round(LX / DX)); ny = int(round(LY / DX))
print('# t  N  max|v|  | core phi(z): mean std noise lag1 | phi wall/interior | vz wall  vz interior  rms column vz')
last = -1e9
for t in sorted(frames):
    P = frames[t]
    if t - last < dtp - 1e-9 and t != max(frames): continue
    last = t
    vmax = max(math.sqrt(p[3]**2 + p[4]**2 + p[5]**2) for p in P)
    cnt = [0] * nz; col = collections.defaultdict(lambda: [0, 0.0])
    for x, y, z, vx, vy, vz in P:
        if ZLO <= z < ZHI:
            cnt[min(int((z - ZLO) / DZ), nz - 1)] += 1
            c = col[(min(int(x / DX), nx - 1), min(int(y / DX), ny - 1))]; c[0] += 1; c[1] += vz
    phi = [c * vp / (LX * LY * DZ) for c in cnt]; m = sum(phi) / nz
    sd = math.sqrt(sum((p - m)**2 for p in phi) / nz)
    noise = m / math.sqrt(max(sum(cnt) / nz, 1))
    lag = sum((phi[i] - m) * (phi[i + 1] - m) for i in range(nz - 1)) / max(sum((p - m)**2 for p in phi), 1e-300)
    wall = [0, 0.0, 0]; inner = [0, 0.0, 0]; cm = []
    for i in range(nx):
        for j in range(ny):
            n, s = col[(i, j)]; tgt = wall if i in (0, nx - 1) or j in (0, ny - 1) else inner
            tgt[0] += n; tgt[1] += s; tgt[2] += 1
            if n: cm.append(s / n)
    vcol = DX * DX * (ZHI - ZLO); mall = sum(cm) / len(cm)
    print('%.2f %d %.4f | %.4f %.4f %.4f %+.2f | %.4f %.4f | %.5f %.5f %.5f' % (
        t, len(P), vmax, m, sd, noise, lag,
        wall[0] * vp / (wall[2] * vcol), inner[0] * vp / (inner[2] * vcol),
        wall[1] / max(wall[0], 1), inner[1] / max(inner[0], 1),
        math.sqrt(sum((c - mall)**2 for c in cm) / len(cm))))
