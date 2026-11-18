#!/usr/bin/env python3
"""Front positions for the batch-settling case (spec 10) from FeatFloWer particle logs
(el_particles_rank*.log: time ip fx fy fz tx ty tz px py pz vx vy vz, written every
SimPar@ELParticleLogFreq steps). Definitions and arithmetic are those of the peer's
reference/batch_eval.py (profile, fronts, fits copied unchanged), only the reader differs.
usage: batch_eval_ff.py RUNDIR PHI0 [T_FIT0 T_FIT1]   (env: LX LY H0 DZ)
Profile phi(z): exact sphere-slab overlap in bins of DZ (default 2 mm).
Upper front: highest z where phi crosses phi0/2. Lower front: lowest z above the floor where phi
falls through (phi0 + phi_sed)/2, phi_sed = mean of the sediment plateau (bins below the front,
excluding the floor bin). Speeds by least squares over [T_FIT0, T_FIT1] (default 1 s to the end).
Core velocity: mean v_z of the particles more than 10 mm away from both fronts."""
import os, sys, glob, math, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import batch_reference as ref
def env(k, v): return float(os.environ.get(k, v))
d, LX, LY, H0, DZ = ref.d, env('LX', 0.02), env('LY', 0.02), env('H0', 0.16), env('DZ', 2e-3)
run, phi0 = sys.argv[1], float(sys.argv[2]); r = d / 2; A = LX * LY
nb = int(round(0.20 / DZ))
def cap(h):  # volume of sphere below height h above its bottom, 0<=h<=d
    h = min(max(h, 0.0), d); return math.pi * h * h * (3 * r - h) / 3
def profile(zs):
    v = [0.0] * nb
    for z in zs:
        b0 = max(int((z - r) / DZ), 0); b1 = min(int((z + r) / DZ), nb - 1)
        for b in range(b0, b1 + 1):
            v[b] += cap((b + 1) * DZ - (z - r)) - cap(b * DZ - (z - r))
    return [x / (A * DZ) for x in v]
def fronts(phi):
    zc = [(b + 0.5) * DZ for b in range(nb)]
    up = None
    for b in range(nb - 1, 0, -1):
        if phi[b - 1] >= phi0 / 2 > phi[b]:
            up = zc[b - 1] + (phi[b - 1] - phi0 / 2) / (phi[b - 1] - phi[b]) * DZ; break
    k = 1
    while k < nb and phi[k] > 0.45: k += 1
    sed = sum(phi[1:k - 1]) / max(k - 2, 1) if k > 3 else None
    lo = 0.0
    if sed:
        lev = 0.5 * (phi0 + sed)
        for b in range(1, nb - 1):
            if phi[b] >= lev > phi[b + 1]:
                lo = zc[b] + (phi[b] - lev) / (phi[b] - phi[b + 1]) * DZ; break
    return up, lo, sed
frames = collections.defaultdict(lambda: ([], []))
for f in sorted(glob.glob(os.path.join(run, 'el_particles_rank*.log'))):
    with open(f, errors='replace') as fh:
        for l in fh:
            if l.startswith('#'): continue
            w = l.split()
            if len(w) < 14: continue
            fr = frames[round(float(w[0]), 6)]; fr[0].append(float(w[10])); fr[1].append(float(w[13]))
rows = []; prof = {}
for t in sorted(frames):
    zs, vz = frames[t]; phi = profile(zs); up, lo, sed = fronts(phi); prof[t] = phi
    core = [v for z, v in zip(zs, vz) if up and lo + 0.01 < z < up - 0.01]
    rows.append((t, up, lo, sed, sum(core) / len(core) if core else float('nan'), len(zs)))
def fit(k, t0, t1):
    P = [(r_[0], r_[k]) for r_ in rows if t0 <= r_[0] <= t1 and r_[k] is not None]
    n = len(P); mt = sum(p[0] for p in P) / n; my = sum(p[1] for p in P) / n
    s = sum((p[0] - mt) * (p[1] - my) for p in P) / sum((p[0] - mt)**2 for p in P); return s, my - s * mt
t0 = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0; t1 = float(sys.argv[4]) if len(sys.argv) > 4 else rows[-1][0]
phi0a = rows[0][5] * math.pi * d**3 / 6 / (A * H0)
eps = 1 - phi0a; ur = ref.slip(eps); vu = eps * ur
print('# t  z_upper  z_lower  phi_sed  <vz>_core  N')
for r_ in rows[::max(len(rows) // 12, 1)]: print('%.3f %s %.5f %s %.6f %d' % (r_[0], '%.5f' % r_[1] if r_[1] else 'nan', r_[2], '%.4f' % r_[3] if r_[3] else 'nan', r_[4], r_[5]))
print('phi0 actual %.5f  fit window %.2f-%.2f s  frames %d' % (phi0a, t0, t1, len(rows)))
if sum(1 for r_ in rows if t0 <= r_[0] <= t1) < 3: print('not enough frames in the fit window'); sys.exit(0)
su, iu = fit(1, t0, t1); sl, il = fit(2, t0, t1)
seds = [r_[3] for r_ in rows if r_[3] and t0 <= r_[0] <= t1]; sed = sum(seds) / len(seds) if seds else float('nan')
vc = [r_[4] for r_ in rows if t0 <= r_[0] <= t1]; vcm = sum(vc) / len(vc)
print('upper front speed  %.6f m/s   reference %.6f   ratio %.4f   (z at t=0 extrapolated %.5f)' % (-su, vu, -su / vu, iu))
print('core particle vz   %.6f m/s   reference %.6f   ratio %.4f' % (-vcm, vu, -vcm / vu))
vl = phi0a * vu / (sed - phi0a)
print('lower front speed  %.6f m/s   reference %.6f (measured phi_sed %.4f)   ratio %.4f' % (sl, vl, sed, sl / vl))
print('lower front speed reference at PHIPACK %.2f: %.6f' % (ref.phip, phi0a * vu / (ref.phip - phi0a)))
# suspension plateau between the fronts at the last frame (bins 10 mm away from both fronts)
t = rows[-1][0]; up, lo = rows[-1][1], rows[-1][2]; zc = [(b + 0.5) * DZ for b in range(nb)]
pl = [p for z, p in zip(zc, prof[t]) if up and lo + 0.01 < z < up - 0.01]
if pl: print('suspension plateau phi at t %.2f: mean %.4f min %.4f max %.4f (%d bins)' % (t, sum(pl) / len(pl), min(pl), max(pl), len(pl)))
