#!/usr/bin/env python3
"""Compare two q2p1_bench_sedimentation runs of the same case, typically the
serial-PE and the parallel-PE build (campaign gate D0-E4-XMODE), on the sphere
trajectory they print.

Usage: xmode_compare.py <log_A> <log_B> [--label-a serial] [--label-b parallel]
                        [--case E4] [--tol-peak 0.1] [--plot out.png]

Both logs are read for the unified DNS_PART_STATE lines (both PE modes,
SimPar@PrintParticleState = Yes); a log without them falls back to the
serial-mode SED_BENCH_VEL / SED_BENCH_POS tags. Reported: peak settling
velocity and its time for both runs, their relative difference against
--tol-peak (the XMODE gate), the velocity and gap-height deviation over the
settling window (from 0.05 s to the first touchdown, h/d < 0.05, minus 0.05 s)
and over the impact/rebound window after it, and, with --case, both peaks
against the printed ten Cate Table II ratio.

The two runs need not share a partition or rank count, so the comparison is
tolerance based, never bitwise; the sample times are matched to 1e-7 s.
"""
import argparse
import re
import sys

import numpy as np

D_P = 0.015          # sphere diameter [m]
R_P = 0.0075
# ten Cate Table II printed u_max/u_inf ratios x Table I u_inf (see compare_tencate.py)
V_PEAK_PRINTED = {'E1': -0.947 * 0.038, 'E2': -0.953 * 0.060,
                  'E3': -0.959 * 0.091, 'E4': -0.955 * 0.128}

PAT_STATE = re.compile(r'DNS_PART_STATE\s+time=\s*(\S+)\s+id=\s*(\d+)\s+'
                       + r'(\S+)\s+' * 8 + r'(\S+)')
PAT_VEL = re.compile(r'SED_BENCH_VEL\s+time=\s*(\S+)\s+ip=\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)')
PAT_POS = re.compile(r'SED_BENCH_POS\s+time=\s*(\S+)\s+ip=\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)')


def parse_run(path, particle_id=None):
    """Return (t, z, vz) arrays of one particle, from DNS_PART_STATE or the
    serial SED_BENCH tags. The z position and z velocity are what the
    sedimentation benchmark compares. The particle id differs between the
    modes (serial mode numbers the FeatFloWer particle list from 1, parallel
    mode prints PE's system id), so without --id the log must hold exactly
    one particle."""
    by_id = {}
    t_v, v_v, t_z, z_z = [], [], [], []
    with open(path, errors='replace') as f:
        for line in f:
            m = PAT_STATE.search(line)
            if m:
                rec = by_id.setdefault(int(m.group(2)), ([], [], []))
                rec[0].append(float(m.group(1)))
                rec[1].append(float(m.group(5)))      # position(3): groups 3,4,5
                rec[2].append(float(m.group(8)))      # velocity(3): groups 6,7,8
                continue
            m = PAT_VEL.search(line)
            if m:
                t_v.append(float(m.group(1))); v_v.append(float(m.group(4)))
                continue
            m = PAT_POS.search(line)
            if m:
                t_z.append(float(m.group(1))); z_z.append(float(m.group(4)))
    if by_id:
        if particle_id is None:
            if len(by_id) != 1:
                sys.exit('%s: %d particles (%s), choose one with --id'
                         % (path, len(by_id), ', '.join(map(str, sorted(by_id)))))
            particle_id = next(iter(by_id))
        if particle_id not in by_id:
            sys.exit('%s: no DNS_PART_STATE lines for id %d' % (path, particle_id))
        t, z, vz = by_id[particle_id]
        return np.array(t), np.array(z), np.array(vz), 'DNS_PART_STATE id %d' % particle_id
    if t_v:
        t_v, v_v = np.array(t_v), np.array(v_v)
        zz = np.interp(t_v, np.array(t_z), np.array(z_z)) if t_z else np.full_like(t_v, np.nan)
        return t_v, zz, v_v, 'SED_BENCH_VEL/POS'
    sys.exit('%s: no DNS_PART_STATE or SED_BENCH_VEL lines' % path)


def align(ta, tb, tol=1e-7):
    """Indices of the sample times common to both runs."""
    ka = np.round(ta / tol).astype(np.int64)
    kb = np.round(tb / tol).astype(np.int64)
    common, ia, ib = np.intersect1d(ka, kb, return_indices=True)
    return ia, ib


def touchdown(t, z, hd_limit=0.05):
    """First time the gap height h/d = (z - R)/D drops below hd_limit, or None."""
    h = (z - R_P) / D_P
    idx = np.where(h < hd_limit)[0]
    return t[idx[0]] if len(idx) else None


def window_stats(t, va, vb, lo, hi, floor=1e-4):
    """Max absolute and mean relative deviation of vb from va on t in [lo, hi].
    The relative mean skips samples where |va| is below floor (the start)."""
    m = (t >= lo) & (t <= hi)
    if not m.any():
        return 0, np.nan, np.nan
    d = np.abs(vb[m] - va[m])
    rel = np.abs(va[m]) > max(floor, 0.0)
    mean_rel = np.mean(d[rel] / np.abs(va[m][rel])) * 100 if rel.any() else np.nan
    return int(m.sum()), d.max(), mean_rel


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('log_a'); ap.add_argument('log_b')
    ap.add_argument('--label-a', default='serial'); ap.add_argument('--label-b', default='parallel')
    ap.add_argument('--case', choices=sorted(V_PEAK_PRINTED), help='ten Cate case for the Table II context')
    ap.add_argument('--tol-peak', type=float, default=0.1, help='XMODE gate on the peak velocity, percent (default 0.1)')
    ap.add_argument('--plot', help='write a v(t) / h(t) comparison figure')
    ap.add_argument('--id-a', type=int, help='particle id to read from log_A (default: the only one)')
    ap.add_argument('--id-b', type=int, help='particle id to read from log_B (default: the only one)')
    a = ap.parse_args()

    ta, za, va, src_a = parse_run(a.log_a, a.id_a)
    tb, zb, vb, src_b = parse_run(a.log_b, a.id_b)
    print('%-9s: %6d samples, t in [%.4f, %.4f] s, source %s' % (a.label_a, len(ta), ta[0], ta[-1], src_a))
    print('%-9s: %6d samples, t in [%.4f, %.4f] s, source %s' % (a.label_b, len(tb), tb[0], tb[-1], src_b))

    ia, ib = align(ta, tb)
    if len(ia) == 0:
        sys.exit('no common sample times')
    t, za, va, zb, vb = ta[ia], za[ia], va[ia], zb[ib], vb[ib]
    print('common sample times: %d' % len(t))

    # peak settling velocity (most negative vz) and its time
    ipa, ipb = np.argmin(va), np.argmin(vb)
    vpa, tpa, vpb, tpb = va[ipa], t[ipa], vb[ipb], t[ipb]
    dpeak = (vpb / vpa - 1) * 100
    print('%-9s: v_peak = %+.6f m/s at t = %.4f s' % (a.label_a, vpa, tpa))
    print('%-9s: v_peak = %+.6f m/s at t = %.4f s' % (a.label_b, vpb, tpb))
    print('peak difference %s vs %s: %+.4f%% (velocity), %+.4f s (timing)' % (a.label_b, a.label_a, dpeak, tpb - tpa))
    if a.case:
        vpr = V_PEAK_PRINTED[a.case]
        print('Table II %s printed v_peak = %+.5f m/s: %s %+.3f%%, %s %+.3f%%'
              % (a.case, vpr, a.label_a, (vpa / vpr - 1) * 100, a.label_b, (vpb / vpr - 1) * 100))

    # windows: settling up to the first touchdown of either run, impact after it
    td = [x for x in (touchdown(t, za), touchdown(t, zb)) if x is not None]
    t_end = t[-1]
    if td:
        t_imp = min(td)
        print('touchdown (h/d < 0.05): %s %s, %s %s'
              % (a.label_a, '%.4f s' % touchdown(t, za) if touchdown(t, za) else 'none',
                 a.label_b, '%.4f s' % touchdown(t, zb) if touchdown(t, zb) else 'none'))
        settle_hi = max(0.05, t_imp - 0.05)
    else:
        settle_hi = t_end
    windows = [('settling', 0.05, settle_hi)]
    if td:
        windows.append(('impact/rebound', settle_hi, t_end))
    print('%-16s %8s %18s %14s %18s' % ('window', 'samples', 'max|dv| [m/s]', 'mean rel v', 'max|dh/d|'))
    for name, lo, hi in windows:
        n, dmax, mrel = window_stats(t, va, vb, lo, hi)
        _, dhmax, _ = window_stats(t, (za - R_P) / D_P, (zb - R_P) / D_P, lo, hi, floor=-1)
        print('%-16s %8d %18.3e %13.3f%% %18.3e   t in [%.3f, %.3f]' % (name, n, dmax, mrel, dhmax, lo, hi))

    verdict = 'PASS' if abs(dpeak) <= a.tol_peak else 'FAIL'
    print('XMODE_GATE v_peak %s vs %s: %+.4f%% (tol %.2f%%) -> %s' % (a.label_b, a.label_a, dpeak, a.tol_peak, verdict))

    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
        ax[0].plot(t, va, label=a.label_a); ax[0].plot(t, vb, '--', label=a.label_b)
        ax[0].set_ylabel('v_z [m/s]'); ax[0].legend(); ax[0].grid(alpha=.3)
        ax[1].plot(t, (za - R_P) / D_P, label=a.label_a); ax[1].plot(t, (zb - R_P) / D_P, '--', label=a.label_b)
        ax[1].set_ylabel('h/d'); ax[1].set_xlabel('t [s]'); ax[1].grid(alpha=.3)
        fig.suptitle('%s vs %s%s' % (a.label_a, a.label_b, ' (%s)' % a.case if a.case else ''))
        fig.tight_layout(); fig.savefig(a.plot, dpi=130)
        print('plot: %s' % a.plot)
    return 0 if verdict == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
