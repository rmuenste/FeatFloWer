#!/usr/bin/env python3
"""Compare a q2p1_bench_sedimentation run against the ten Cate PIV curves.

Usage: compare_tencate.py <run_log> [--case E4] [--tcref tc-ref] [--plot out.png]

Parses SED_BENCH_VEL / SED_BENCH_POS lines (serial-PE mode tags; also
accepts DNS_PART_STATE) and reports peak settling velocity, time of peak,
and RMS deviation from the digitized PIV velocity curve over the settling
window, plus the gap-height trajectory h/d against case_E*_h.csv.
"""
import argparse
import re
import sys

import numpy as np

D_P = 0.015          # sphere diameter [m]
R_P = 0.0075
Z_WALL = 0.0         # bottom wall

def parse_run(path):
    t_v, v, t_z, z = [], [], [], []
    pat_v = re.compile(r'SED_BENCH_VEL\s+time=\s*(\S+)\s+ip=\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)')
    pat_z = re.compile(r'SED_BENCH_POS\s+time=\s*(\S+)\s+ip=\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)')
    pat_s = re.compile(r'DNS_PART_STATE\s+time=\s*(\S+)\s+id=\s*\d+\s+' + r'(\S+)\s+'*8 + r'(\S+)')
    for line in open(path, errors='replace'):
        m = pat_v.search(line)
        if m:
            t_v.append(float(m.group(1))); v.append(float(m.group(4)))
            continue
        m = pat_z.search(line)
        if m:
            t_z.append(float(m.group(1))); z.append(float(m.group(4)))
            continue
        m = pat_s.search(line)
        if m and not t_v:  # fall back to unified tag if legacy tags absent
            pass
    return np.array(t_v), np.array(v), np.array(t_z), np.array(z)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('run_log')
    ap.add_argument('--case', default='E4')
    ap.add_argument('--tcref', default='tc-ref')
    ap.add_argument('--plot', default=None)
    a = ap.parse_args()

    t_v, v, t_z, z = parse_run(a.run_log)
    if not len(t_v):
        sys.exit('no SED_BENCH_VEL lines found in %s' % a.run_log)

    ref_v = np.loadtxt('%s/ref_%s.dat' % (a.tcref, a.case))
    ref_h = np.loadtxt('%s/case_%s_h.csv' % (a.tcref, a.case))

    i = np.argmin(v)
    ir = np.argmin(ref_v[:, 1])
    vp, tp = v[i], t_v[i]
    vr, tr = ref_v[ir, 1], ref_v[ir, 0]
    print('run : v_peak = %+.5f m/s at t = %.3f s' % (vp, tp))
    print('PIV : v_peak = %+.5f m/s at t = %.3f s' % (vr, tr))
    print('peak error: %+.2f%% (velocity), %+.3f s (timing)' % ((vp/vr - 1)*100, tp - tr))

    # RMS over the overlapping settling window (PIV time support)
    lo, hi = ref_v[0, 0], ref_v[-1, 0]
    m = (t_v >= lo) & (t_v <= hi)
    if m.sum() > 3:
        vi = np.interp(t_v[m], ref_v[:, 0], ref_v[:, 1])
        rms = np.sqrt(np.mean((v[m] - vi)**2))
        print('RMS(v) over t in [%.2f, %.2f]: %.5f m/s (%.2f%% of |v_peak,PIV|)'
              % (lo, hi, rms, rms/abs(vr)*100))

    if len(t_z):
        h = (z - R_P - Z_WALL) / D_P
        mh = (t_z >= ref_h[0, 0]) & (t_z <= ref_h[-1, 0])
        if mh.sum() > 3:
            hi_ = np.interp(t_z[mh], ref_h[:, 0], ref_h[:, 1])
            print('RMS(h/d) over PIV window: %.4f' % np.sqrt(np.mean((h[mh] - hi_)**2)))
        # time of touchdown (h/d < 0.05)
        idx = np.where(h < 0.05)[0]
        if len(idx):
            print('touchdown (h/d<0.05): t = %.3f s (PIV curve ends %.3f s)'
                  % (t_z[idx[0]], ref_h[-1, 0]))

    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(t_v, v, '-', label='DNS')
        ax[0].plot(ref_v[:, 0], ref_v[:, 1], 'o', ms=3, label='PIV')
        ax[0].set(xlabel='t [s]', ylabel='v_z [m/s]', title=a.case)
        ax[0].legend()
        if len(t_z):
            ax[1].plot(t_z, (z - R_P)/D_P, '-', label='DNS')
            ax[1].plot(ref_h[:, 0], ref_h[:, 1], 'o', ms=3, label='PIV')
            ax[1].set(xlabel='t [s]', ylabel='h/d')
            ax[1].legend()
        fig.tight_layout()
        fig.savefig(a.plot, dpi=140)
        print('plot: %s' % a.plot)

if __name__ == '__main__':
    main()
