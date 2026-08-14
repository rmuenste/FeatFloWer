#!/usr/bin/env python3
"""D2.1 figure: measured wall-approach drag vs Brenner's exact solution.

Top: lambda = F/(6 pi mu R V) against eps = gap/R (log-log), Brenner exact
curve + the constant-V ladder runs (d21_l3, d21_l4). Points follow the
curve and fall off below the per-level resolution limit.
Bottom: deviation from Brenner against gap measured in CELLS (gap/h) - the
crossover map; -10%/-20% lines and the gap=2h rule marked.

Output: docs/md_docs/dns_figures/d21_brenner_crossover.png
"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d21_d12_tail_analysis import parse, brenner_lambda, H

R, MU, V = 0.0075, 0.373, 0.02
RUNS = [('q2p1_dns_rundir_d21_l3', 'l3', 23.9, '#c44e52', 'o'),
        ('q2p1_dns_rundir_d21_l4', 'l4', 49.1, '#4c72b0', 's')]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 8.4),
                               gridspec_kw={'height_ratios': [3, 2]})

eps_c = np.logspace(math.log10(0.04), math.log10(4.5), 300)
lam_c = [brenner_lambda(e) for e in eps_c]
ax1.loglog(eps_c, lam_c, 'k-', lw=1.8, label='Brenner (1961), exact', zorder=3)
ax1.loglog(eps_c, 1 / eps_c + 0.9713, ls=':', c='0.6', lw=1,
           label=r'lubrication asymptote $\lambda\sim1/\varepsilon$')

for run, lvl, dh, col, mk in RUNS:
    t, fz, z, vz = parse(run)
    gap = z - R
    eps = gap / R
    lam = fz / (6 * math.pi * MU * R * V)
    ok = (t > 0.2) & (eps > 5e-3)
    eps, lam, gap = eps[ok], lam[ok], gap[ok]
    # subsample uniformly in log(eps) so the near-wall tail stays populated
    targets = np.logspace(math.log10(eps.min()), math.log10(eps.max()), 70)
    s = np.unique([np.argmin(np.abs(eps - e0)) for e0 in targets])
    ax1.loglog(eps[s], lam[s], mk, ms=4, mfc='none', mec=col, mew=1.1,
               label='FBM, D/h = %.0f' % dh, zorder=4)
    # crossover marker at gap = 2h
    e2h = 2 * H[lvl] / R
    ax1.axvline(e2h, color=col, ls='--', lw=0.9, alpha=0.6)
    ax1.text(e2h, 1.55, ' 2h (D/h=%.0f)' % dh, rotation=90, va='bottom',
             ha='right', fontsize=8, color=col)
    # bottom panel: deviation vs gap in cells
    lamB = np.array([brenner_lambda(e) for e in eps])
    dev = 100 * (lam / lamB - 1)
    gh = gap / H[lvl]
    ax2.semilogx(gh[s], dev[s], mk, ms=4, mfc='none', mec=col, mew=1.1,
                 label='D/h = %.0f' % dh)

ax1.set_xlabel(r'gap / R   ($\varepsilon$)')
ax1.set_ylabel(r'drag amplification  $\lambda = F_z\,/\,6\pi\mu R V$')
ax1.set_title('Sphere approaching a wall at constant V — measured drag vs Brenner\n'
              r'(E1 fluid $\mu$=0.373, V=0.02 m/s, Re=0.78, $\rho_p/\rho_f=10^6$)',
              fontsize=10)
ax1.legend(fontsize=8, loc='upper right')
ax1.grid(True, which='both', alpha=0.25, lw=0.4)

ax2.axhline(0, color='k', lw=0.8)
ax2.axhline(-10, color='0.5', ls='--', lw=0.8)
ax2.axhline(-20, color='0.5', ls=':', lw=0.8)
ax2.text(30, -9.2, '-10%', fontsize=8, color='0.4')
ax2.text(30, -19.2, '-20%', fontsize=8, color='0.4')
ax2.axvspan(0.5, 2.0, color='orange', alpha=0.12)
ax2.text(1.0, 6, 'sub-grid:\ncorrection\ntakes over', fontsize=8,
         ha='center', color='#a06000')
ax2.set_xlabel('gap / h   (cells)')
ax2.set_ylabel('deviation from Brenner  [%]')
ax2.set_xlim(0.5, 60)
ax2.set_ylim(-45, 12)
ax2.legend(fontsize=8, loc='lower right')
ax2.grid(True, which='both', alpha=0.25, lw=0.4)
ax2.set_title('Crossover map: same data in cell units', fontsize=10)

fig.tight_layout()
out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'docs', 'md_docs', 'dns_figures', 'd21_brenner_crossover.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=150)
print('wrote', out)
