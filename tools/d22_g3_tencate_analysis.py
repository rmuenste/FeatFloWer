#!/usr/bin/env python3
"""D2.2 gate G3: lubricated ten Cate bottom approach (paper Fig. 13 analog).

Runs: q2p1_dns_rundir_e{1,2}_l3_g3def (jobs 141665/141666) - the certified D1.3
E1/E2 fixtures with the lubrication add-on ON in the production configuration
(kroupaDeficit, mesh clamp c = 2, saturation eps_c = 0.1 and Vinogradova slip
active). Baselines: the certified lubrication-OFF runs (our "No F_lub" curves).
Reference: digitized PIV curves tc-ref/ref_E{1,2}.dat (3.5-4% peak-speed
digitization audit applies; curves authoritative for shape).

MEASURED (both cases): identical to baseline down to the activation gap 2h
(deficit inert above the band by construction), then a more gradual
deceleration through the film band - u(1h) slower by 15% (E1) / 16% (E2),
u(0.5h) by 12% / 17% - and a finite landing: 2h->rest 334 ms vs 301 (E1),
140 vs 124 ms (E2), i.e. +11-13%, NOT the paper's unbounded
time-to-contact tail (their pure-lubrication force never lands the sphere;
the deficit + saturation + hard contact do). Peak lubrication force reaches
only 0.17-0.23 x buoyant weight: at D/h = 23.9 the resolved FBM already
carries most of the film (G2: FBM alone is only -15/-26% below Brenner), so
the correction is modest BY CONSTRUCTION - unlike ten Cate's LBM, which had
no resolved film below one grid spacing and needed the full Eq.-10 force.
The digitized PIV cannot discriminate at these scales (~11 samples in the
approach window, near-zero velocities at the measurement floor); the
quantitative certification of the model is G2/G2b, this gate certifies the
qualitative Fig.-13 behavior and the absence of the landing pathology.

Output: docs/md_docs/dns_figures/d22_g3_tencate.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 0.0075
H = 6.27151451e-4   # L3 h_min; activation gap = 2H

def series(run, tag='SED_BENCH_VEL', col=7):
    t, v = [], []
    for line in open(os.path.join(run, 'run_slurm.log')):
        if line.startswith(tag):
            p = line.split(); t.append(float(p[2])); v.append(float(p[col]))
    return np.array(t), np.array(v)

CASES = (('E1', 'q2p1_dns_rundir_e1_l3', 'q2p1_dns_rundir_e1_l3_g3def', (3.0, 4.4)),
         ('E2', 'q2p1_dns_rundir_e2_l3', 'q2p1_dns_rundir_e2_l3_g3def', (1.9, 2.8)))

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
for ax, (case, base, g3, (t0, t1)) in zip(axes, CASES):
    ref = np.loadtxt('tc-ref/ref_%s.dat' % case)
    tb, vb = series(base); tg, vg = series(g3)
    zb = series(base, 'SED_BENCH_POS')[1]
    t2h = tb[np.where(zb - R <= 2 * H)[0][0]]

    m = (ref[:, 0] >= t0) & (ref[:, 0] <= t1)
    ax.plot(ref[m, 0], ref[m, 1], 's', ms=5, mfc='none', mec='k', mew=1.1,
            label='Measurement (PIV, digitized)')
    mb = (tb >= t0) & (tb <= t1)
    ax.plot(tb[mb], vb[mb], '--', c='#c44e52', lw=1.6, label='FBM, no $F_{lub}$')
    mg = (tg >= t0) & (tg <= t1)
    ax.plot(tg[mg], vg[mg], '-', c='#55a868', lw=1.6,
            label='FBM + kroupaDeficit (production config)')
    ax.axvline(t2h, color='0.5', ls=':', lw=1)
    ax.text(t2h, ax.get_ylim()[0] * 0.05, ' gap = 2h', rotation=90,
            va='bottom', ha='right', fontsize=8, color='0.4')
    ax.set_xlabel('t [s]')
    ax.set_ylabel('u [m/s]')
    ax.set_title('%s (Re = %s) — bottom approach' %
                 (case, '1.5' if case == 'E1' else '4.1'))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, loc='lower right')

fig.suptitle('G3: lubricated ten Cate bottom approach (Fig. 13 analog) — '
             'jobs 141665/141666, D/h = 23.9', fontsize=11)
out = 'docs/md_docs/dns_figures/d22_g3_tencate.png'
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(out, dpi=150)
print('wrote', out)
