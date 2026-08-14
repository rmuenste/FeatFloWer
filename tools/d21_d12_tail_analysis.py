#!/usr/bin/env python3
"""D2.1 free-approach + D1.2 grid-noise extraction from existing ten Cate runs.

D2.1 (sphere-wall approach, free variant): the tail of every ten Cate run is
a sphere settling onto the bottom wall with force and position logged every
step. Measured drag amplification lambda = F_z/(6 pi mu R |v|) vs the
Brenner (1961) exact perpendicular-translation solution lambda_B(eps),
eps = gap/R. Deliverable: the gap (in cells h) where the resolved force
departs from the analytic divergence -- the crossover map feeding the PE
lubrication-threshold rule. Quasi-steady Stokes reference is only honest at
E1 (Re 1.5); E2-E4 are recorded as context, not gated.

D1.2 (grid-crossing noise): in the peak-velocity plateau the sphere moves at
~constant speed across the mesh; force fluctuation there is the intrinsic
FBM re-classification noise. Amplitude (rms of detrended F_z, % of mean)
vs D/h and cells-crossed-per-step.

Wall at z=0 is the full-box bottom; sphere R=0.0075 m. ten Cate Table I
fluids. Runs: q2p1_dns_rundir_{e1..e4}_{l2,l3,l4} (dt=1ms rows).
"""
import math, os, re, sys
import numpy as np

R = 0.0075
FLUID = {  # case: (mu [Pa s], rho_f, u_inf [m/s])  ten Cate Table I
    'e1': (0.373, 970.0, 0.038),
    'e2': (0.212, 965.0, 0.060),
    'e3': (0.113, 962.0, 0.091),
    'e4': (0.058, 960.0, 0.128),
}
DH = {'l2': 11.4, 'l3': 23.9, 'l4': 49.1}
H = {L: 2 * R / dh for L, dh in DH.items()}   # cell size [m]
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def parse(run):
    path = os.path.join(ROOT, run, 'run_slurm.log')
    t, fz, z, vz = [], [], [], []
    pat = re.compile(r'SED_BENCH_(FORCE|POS|VEL)\s+time=\s*([0-9.E+-]+)\s+ip=\s*1\s+'
                     r'([0-9.E+-]+)\s+([0-9.E+-]+)\s+([0-9.E+-]+)')
    rows = {}
    for line in open(path, errors='replace'):
        m = pat.search(line)
        if not m:
            continue
        kind, tt = m.group(1), float(m.group(2))
        rows.setdefault(tt, {})[kind] = float(m.group(5))
    for tt in sorted(rows):
        r = rows[tt]
        if all(k in r for k in ('FORCE', 'POS', 'VEL')):
            t.append(tt); fz.append(r['FORCE']); z.append(r['POS']); vz.append(r['VEL'])
    return (np.array(t), np.array(fz), np.array(z), np.array(vz))


def brenner_lambda(eps):
    """Brenner 1961 exact drag correction, sphere translating toward a wall."""
    if eps <= 0:
        return math.inf
    a = math.acosh(1.0 + eps)
    s = 0.0
    for n in range(1, 6000):
        num = 2 * math.sinh((2 * n + 1) * a) + (2 * n + 1) * math.sinh(2 * a)
        den = 4 * math.sinh((n + 0.5) * a) ** 2 - (2 * n + 1) ** 2 * math.sinh(a) ** 2
        term = n * (n + 1) / ((2 * n - 1.0) * (2 * n + 3.0)) * (num / den - 1.0)
        s += term
        if abs(term) < 1e-12 * abs(s) and n > 10:
            break
    return (4.0 / 3.0) * math.sinh(a) * s


def d21(case, lvl, t, fz, z, vz):
    mu, rho, _ = FLUID[case]
    h = H[lvl]
    gap = z - R                       # wall at z=0
    eps = gap / R
    ok = (vz < -1e-4) & (eps > 1e-4)  # falling, not yet in contact
    lam = fz / (6 * math.pi * mu * R * (-vz))
    lamB = np.array([brenner_lambda(e) for e in eps])
    dev = lam / lamB - 1.0
    # crossover: last gap (approaching) where |deviation| exceeds 20%
    idx = np.where(ok & (eps < 4.0))[0]
    cross_h = None
    for i in idx:
        if abs(dev[i]) > 0.20 and gap[i] < gap[idx[0]]:
            cross_h = gap[i] / h
            break
    # quasi-steadiness diagnostic in the tail (eps<1): inertia term vs force
    tail = ok & (eps < 1.0)
    qs = None
    if tail.sum() > 5:
        dvdt = np.gradient(vz, t)
        m_eff = (1120.0 + 0.5 * rho) * (4. / 3.) * math.pi * R ** 3
        qs = float(np.median(np.abs(m_eff * dvdt[tail]) / np.abs(fz[tail])))
    # sample lambda at eps = 1, 0.5, 0.25, 0.1
    samples = {}
    for e0 in (1.0, 0.5, 0.25, 0.1):
        if (ok & (eps < e0)).any():
            i = np.argmin(np.abs(eps - e0) + 1e9 * (~ok))
            samples[e0] = (float(lam[i]), brenner_lambda(e0), float(gap[i] / h))
    return cross_h, samples, qs, float(eps[ok][-1]) if ok.any() else None


def d12(case, lvl, t, fz, z, vz):
    ipk = np.argmin(vz)
    vpk = vz[ipk]
    win = np.abs(vz - vpk) < 0.02 * abs(vpk)
    # contiguous window around the peak
    i0, i1 = ipk, ipk
    while i0 > 0 and win[i0 - 1]:
        i0 -= 1
    while i1 < len(win) - 1 and win[i1 + 1]:
        i1 += 1
    if i1 - i0 < 20:
        return None
    tw, fw = t[i0:i1], fz[i0:i1]
    c = np.polyfit(tw, fw, 2)
    resid = fw - np.polyval(c, tw)
    dt = np.median(np.diff(tw))
    cells_per_step = abs(vpk) * dt / H[lvl]
    return dict(n=i1 - i0, rms_rel=float(np.std(resid) / abs(np.mean(fw))),
                vjit=float(np.std(np.diff(vz[i0:i1]))), cps=float(cells_per_step))


def main():
    runs = [(c, L) for c in ('e1', 'e2', 'e3', 'e4') for L in ('l2', 'l3', 'l4')]
    print('=== D1.2 grid-crossing noise (peak plateau, dt=1ms runs) ===')
    print('%-8s %6s %8s %12s %10s %8s' % ('run', 'D/h', 'n_win', 'F rms/|F|', 'dv jitter', 'cells/step'))
    for c, L in runs:
        run = 'q2p1_dns_rundir_%s_%s' % (c, L)
        if not os.path.exists(os.path.join(ROOT, run, 'run_slurm.log')):
            continue
        t, fz, z, vz = parse(run)
        if len(t) < 100:
            continue
        r = d12(c, L, t, fz, z, vz)
        if r:
            print('%-8s %6.1f %8d %11.2e%% %10.2e %8.3f'
                  % (c + '_' + L, DH[L], r['n'], 100 * r['rms_rel'], r['vjit'], r['cps']))
    print()
    print('=== D2.1 wall approach vs Brenner (free variant) ===')
    print('run       D/h | lambda meas/Brenner at eps=1/0.5/0.25/0.1 | cross(|dev|>20%) gap/h | QS |F_in/F| | eps_end')
    for c, L in runs:
        run = 'q2p1_dns_rundir_%s_%s' % (c, L)
        if not os.path.exists(os.path.join(ROOT, run, 'run_slurm.log')):
            continue
        t, fz, z, vz = parse(run)
        if len(t) < 100:
            continue
        cross, samples, qs, eend = d21(c, L, t, fz, z, vz)
        ss = '  '.join('%.2f/%.2f' % (v[0], v[1]) for e, v in sorted(samples.items(), reverse=True))
        print('%-8s %6.1f | %s | %s | %s | %s'
              % (c + '_' + L, DH[L], ss,
                 ('%.2f' % cross) if cross else 'none<4R',
                 ('%.3f' % qs) if qs is not None else '-',
                 ('%.3f' % eend) if eend is not None else '-'))


if __name__ == '__main__':
    main()
