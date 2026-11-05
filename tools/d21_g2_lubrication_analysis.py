#!/usr/bin/env python3
"""D2.2 gate G2: lubricated Brenner wall-approach - overlap/double-count measurement.

Run: q2p1_dns_rundir_d21_l3_g2lub (job 141393) - the certified d21_l3 constant-V
fixture (R = 7.5 mm, V = 0.02 m/s, ten Cate E1 fluid, D/h = 23.9) with the pe
lubrication add-on enabled: Kroupa-2016 wall resistances, bare-resistance
configuration (no slip correction, no h_c saturation), mesh clamp c = 2
(activation gap 2h = 1.2543 mm; first live use of the set_lubrication_mesh_dx
wiring - onset measured at t = 1.563, predicted 1.5623).

KEY FACT (measured): the G2 run's FBM force equals the lubrication-OFF baseline
run to >= 4 digits in every band - the near-rigid particle's trajectory is
unchanged, so the resolved flow is too. Superposition holds, and any model
variant F_lub(gap) can therefore be evaluated OFFLINE against this one dataset:
lambda_sum = lambda_FBM(baseline) + lambda_model(gap).

Measured verdict table (deviation of lambda_sum vs Brenner exact):

  variant                          1h-2h band       <1h band
  FBM alone                        -14.8 %          -25.9 %
  full Kroupa wall C_n, act 2h     +74.6 %  (run)   +67.4 %  (run)
  bare 1/h, act 2h                 +70.1 %          +63.8 %
  Eq.-10 deficit (1/h - 1/2h)       +6.7 %          +23.0 %
  Eq.-10 deficit (1/h - 1/1h)      -14.8 % (off)    -17.8 %

DECISION (DESIGN_SPEC d22_lubrication section 7): the full resistance set
double-counts the FBM-resolved film by ~+75 % in the overlap - the deficit
form (Kroupa minus its value at activation, reducing to ten Cate Eq. 10 in
the leading normal term) is REQUIRED for CFD-coupled runs. The ideal
correction lambda_ideal = Brenner - FBM rises 0.7 -> 4.5 across
gap/h in [0.83, 2): that curve, resolution-collapsed in gap/h (D2.1
collapse finding), is the calibrated deficit D2.1's inverse-collapse idea
anticipated - measured here, one calibration per method.

Output: docs/md_docs/dns_figures/d21_g2_lubrication.png
"""
import math, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d21_d12_tail_analysis import parse, brenner_lambda, H

R, MU, V = 0.0075, 0.373, 0.02
h = H['l3']
DENOM = 6 * math.pi * MU * R * V

def parse_dns_lub(run):
    lub = {}
    for line in open(os.path.join(run, 'run_slurm.log')):
        if line.startswith('DNS_LUB'):
            p = line.split()
            lub[round(float(p[2]), 4)] = float(p[4])
    return lub

t, fz, z, vz = parse('q2p1_dns_rundir_d21_l3_g2lub')
lub = parse_dns_lub('q2p1_dns_rundir_d21_l3_g2lub')
flub = np.array([lub.get(round(tt, 4), 0.0) for tt in t])

gap = z - R
gh = gap / h
eps = gap / R
lamF = fz / DENOM
lamL = flub / DENOM
lamS = lamF + lamL
lamB = np.array([brenner_lambda(e) for e in eps])

act = 2 * h                      # mesh-clamp activation gap (armed, measured onset 1.563)
deficit = np.where(gap < act, R / gap - R / act, 0.0)   # ten Cate Eq.-10 class

def banddev(mask, lam):
    d = 100 * (lam[mask] / lamB[mask] - 1)
    return d.mean(), d.min(), d.max(), mask.sum()

m12 = (gh >= 1) & (gh < 2)
m01 = (gh < 1) & (gap > 0)
print('deviation vs Brenner:  mean [min,max] % (n)')
for name, lam in [('FBM alone', lamF), ('FBM + full Kroupa (run)', lamS),
                  ('FBM + Eq.10 deficit c=2', lamF + deficit)]:
    for mk, lbl in [(m12, '1h-2h'), (m01, '<1h ')]:
        mn, lo, hi, n = banddev(mk, lam)
        print(f'  {name:26s} {lbl}: {mn:+6.1f} [{lo:+6.1f},{hi:+6.1f}] ({n})')

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.2, 8.4),
                               gridspec_kw={'height_ratios': [3, 2]})

eps_c = np.logspace(math.log10(0.04), math.log10(4.5), 300)
ax1.loglog(eps_c, [brenner_lambda(e) for e in eps_c], 'k-', lw=1.8,
           label='Brenner (1961), exact', zorder=3)

ok = (t > 0.2) & (eps > 5e-3)
targets = np.logspace(math.log10(eps[ok].min()), math.log10(eps[ok].max()), 80)
s = np.unique([np.argmin(np.abs(eps - e0)) for e0 in targets])
ax1.loglog(eps[s], lamF[s], 'o', ms=4, mfc='none', mec='#c44e52', mew=1.1,
           label='FBM alone (identical to lub-OFF baseline)', zorder=4)
inb = s[gap[s] < act]
ax1.loglog(eps[inb], lamS[inb], 's', ms=4.5, mfc='none', mec='#8172b2', mew=1.2,
           label='FBM + full Kroupa wall term (measured)', zorder=5)
ax1.loglog(eps[inb], (lamF + deficit)[inb], '^', ms=4.5, mfc='none',
           mec='#55a868', mew=1.2,
           label='FBM + Eq.-10 deficit (offline, superposition)', zorder=5)
ax1.axvline(act / R, color='0.4', ls='--', lw=0.9)
ax1.text(act / R, 1.5, ' activation 2h (mesh clamp, armed)', rotation=90,
         va='bottom', ha='right', fontsize=8, color='0.3')
ax1.set_xlabel(r'gap / R   ($\varepsilon$)')
ax1.set_ylabel(r'$\lambda = F_z\,/\,6\pi\mu R V$')
ax1.set_title('G2: lubricated wall approach, D/h = 23.9 (job 141393)\n'
              'full resistance double-counts the resolved film; deficit form closes it')
ax1.legend(fontsize=8, loc='upper right')
ax1.grid(alpha=0.25, which='both')

for lam, col, mk, lbl in [(lamF, '#c44e52', 'o', 'FBM alone'),
                          (lamS, '#8172b2', 's', '+ full Kroupa'),
                          (lamF + deficit, '#55a868', '^', '+ Eq.-10 deficit')]:
    dev = 100 * (lam / lamB - 1)
    ax2.semilogx(gh[s], dev[s], mk, ms=4, mfc='none', mec=col, mew=1.1, label=lbl)
ax2.axhline(0, color='k', lw=0.8)
ax2.axvspan(gh[gap > 0].min(), 2, color='0.9', zorder=0)
ax2.axvline(2, color='0.4', ls='--', lw=0.9)
ax2.set_xlabel('gap / h  (cells)')
ax2.set_ylabel('deviation from Brenner  [%]')
ax2.set_ylim(-40, 100)
ax2.legend(fontsize=8)
ax2.grid(alpha=0.25, which='both')

out = 'docs/md_docs/dns_figures/d21_g2_lubrication.png'
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.tight_layout()
fig.savefig(out, dpi=150)
print('wrote', out)
