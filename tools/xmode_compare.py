#!/usr/bin/env python3
"""Compare two FeatFloWer DNS runs of the same case, typically the serial-PE and
the parallel-PE build (campaign gate D0-E4-XMODE and its extension to DKT and
the D3.2 settling clouds), on the particle trajectories they print.

Usage: xmode_compare.py <log_A> <log_B> [--label-a serial] [--label-b parallel]
                        [--case E4] [--tol-peak 0.1] [--plot out.png]
                        [--id-a N --id-b M] [--diameter D] [--plateau T0 T1]

Both logs are read for the unified DNS_PART_STATE lines (both PE modes,
SimPar@PrintParticleState = Yes); a log without them falls back to the
serial-mode SED_BENCH_VEL / SED_BENCH_POS tags of q2p1_bench_sedimentation.

One particle per log (ten Cate): peak settling velocity and its time for both
runs, their relative difference against --tol-peak (the XMODE gate), the
velocity and gap-height deviation over the settling window (from 0.05 s to the
first touchdown, h/d < 0.05, minus 0.05 s) and over the impact/rebound window
after it, and, with --case, both peaks against the printed ten Cate Table II
ratio.

Several particles per log (DKT pair, settling clouds): the particles of the two
runs are matched by their first printed position, since the ids differ between
the modes (serial mode numbers the FeatFloWer particle list, parallel mode
prints PE's system id). Reported per particle: peak settling velocity in both
runs and the difference; for the cloud: the mean settling velocity <v_z>(t),
its peak and the same window deviations. The gate is on the cloud-mean peak,
or, with --plateau T0 T1, on the time mean of <v_z> over that window (the D3.2
hindered-settling observable <v_settle>); the per-particle spread is reported,
not gated (close encounters make single trajectories partition sensitive).

The two runs need not share a partition or rank count, so the comparison is
tolerance based, never bitwise; the sample times are matched to 1e-7 s.
"""
import argparse
import re
import sys

import numpy as np

# ten Cate Table II printed u_max/u_inf ratios x Table I u_inf (see compare_tencate.py)
V_PEAK_PRINTED = {'E1': -0.947 * 0.038, 'E2': -0.953 * 0.060,
                  'E3': -0.959 * 0.091, 'E4': -0.955 * 0.128}

PAT_STATE = re.compile(r'DNS_PART_STATE\s+time=\s*(\S+)\s+id=\s*(\d+)\s+'
                       + r'(\S+)\s+' * 8 + r'(\S+)')
PAT_VEL = re.compile(r'SED_BENCH_VEL\s+time=\s*(\S+)\s+ip=\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)')
PAT_POS = re.compile(r'SED_BENCH_POS\s+time=\s*(\S+)\s+ip=\s*\d+\s+(\S+)\s+(\S+)\s+(\S+)')


class Track:
    """Trajectory of one particle: sample times, start position, z, v_z."""

    def __init__(self, pid, t, p0, z, vz):
        self.id, self.t, self.p0, self.z, self.vz = pid, np.array(t), np.array(p0), np.array(z), np.array(vz)


def parse_all(path):
    """Return {id: Track} from the DNS_PART_STATE lines, or the single serial
    SED_BENCH track under id 0 when the unified lines are absent."""
    by_id = {}
    t_v, v_v, t_z, z_z = [], [], [], []
    with open(path, errors='replace') as f:
        for line in f:
            m = PAT_STATE.search(line)
            if m:
                pid = int(m.group(2))
                rec = by_id.setdefault(pid, ([], None, [], []))
                rec[0].append(float(m.group(1)))
                if rec[1] is None:
                    by_id[pid] = (rec[0], [float(m.group(3)), float(m.group(4)), float(m.group(5))], rec[2], rec[3])
                    rec = by_id[pid]
                rec[2].append(float(m.group(5)))      # position(3): groups 3,4,5
                rec[3].append(float(m.group(8)))      # velocity(3): groups 6,7,8
                continue
            m = PAT_VEL.search(line)
            if m:
                t_v.append(float(m.group(1))); v_v.append(float(m.group(4)))
                continue
            m = PAT_POS.search(line)
            if m:
                t_z.append(float(m.group(1))); z_z.append(float(m.group(4)))
    if by_id:
        return {pid: Track(pid, r[0], r[1], r[2], r[3]) for pid, r in by_id.items()}, 'DNS_PART_STATE'
    if t_v:
        t_v, v_v = np.array(t_v), np.array(v_v)
        zz = np.interp(t_v, np.array(t_z), np.array(z_z)) if t_z else np.full_like(t_v, np.nan)
        p0 = [0.0, 0.0, zz[0] if len(zz) else 0.0]
        return {0: Track(0, t_v, p0, zz, v_v)}, 'SED_BENCH_VEL/POS'
    sys.exit('%s: no DNS_PART_STATE or SED_BENCH_VEL lines' % path)


def select(tracks, particle_id, path):
    if particle_id is None:
        if len(tracks) != 1:
            sys.exit('%s: %d particles (%s), choose one with --id'
                     % (path, len(tracks), ', '.join(map(str, sorted(tracks)))))
        return next(iter(tracks.values()))
    if particle_id not in tracks:
        sys.exit('%s: no DNS_PART_STATE lines for id %d' % (path, particle_id))
    return tracks[particle_id]


def match_by_start(ta, tb, scale):
    """Pair the particles of two runs by their first printed position. Every
    particle must have exactly one partner closer than 1e-6 * scale."""
    pairs = []
    used = set()
    for a in ta.values():
        d = {b.id: np.linalg.norm(a.p0 - b.p0) for b in tb.values() if b.id not in used}
        bid = min(d, key=d.get)
        if d[bid] > 1e-6 * scale:
            sys.exit('particle %d of A (start %s) has no partner in B (nearest %.3e away)'
                     % (a.id, a.p0, d[bid]))
        used.add(bid)
        pairs.append((a, tb[bid]))
    return pairs


def align(ta, tb, tol=1e-7):
    """Indices of the sample times common to both runs."""
    ka = np.round(ta / tol).astype(np.int64)
    kb = np.round(tb / tol).astype(np.int64)
    common, ia, ib = np.intersect1d(ka, kb, return_indices=True)
    return ia, ib


def touchdown(t, z, radius, diameter, hd_limit=0.05):
    """First time the gap height h/d = (z - R)/D drops below hd_limit, or None."""
    h = (z - radius) / diameter
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


def report_windows(t, za, va, zb, vb, radius, diameter, la, lb):
    """Print the settling / impact window table for one z, v_z pair; returns (settle_hi, had_touchdown)."""
    tda, tdb = touchdown(t, za, radius, diameter), touchdown(t, zb, radius, diameter)
    td = [x for x in (tda, tdb) if x is not None]
    t_end = t[-1]
    if td:
        print('touchdown (h/d < 0.05): %s %s, %s %s'
              % (la, '%.4f s' % tda if tda is not None else 'none', lb, '%.4f s' % tdb if tdb is not None else 'none'))
        settle_hi = max(0.05, min(td) - 0.05)
    else:
        settle_hi = t_end
    windows = [('settling', 0.05, settle_hi)]
    if td:
        windows.append(('impact/rebound', settle_hi, t_end))
    print('%-16s %8s %18s %14s %18s' % ('window', 'samples', 'max|dv| [m/s]', 'mean rel v', 'max|dh/d|'))
    for name, lo, hi in windows:
        n, dmax, mrel = window_stats(t, va, vb, lo, hi)
        _, dhmax, _ = window_stats(t, (za - radius) / diameter, (zb - radius) / diameter, lo, hi, floor=-1)
        print('%-16s %8d %18.3e %13.3f%% %18.3e   t in [%.3f, %.3f]' % (name, n, dmax, mrel, dhmax, lo, hi))


def peak(t, v):
    i = np.argmin(v)
    return v[i], t[i]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('log_a'); ap.add_argument('log_b')
    ap.add_argument('--label-a', default='serial'); ap.add_argument('--label-b', default='parallel')
    ap.add_argument('--case', choices=sorted(V_PEAK_PRINTED), help='ten Cate case for the Table II context')
    ap.add_argument('--tol-peak', type=float, default=0.1, help='XMODE gate on the peak velocity, percent (default 0.1)')
    ap.add_argument('--plot', help='write a v(t) / h(t) comparison figure')
    ap.add_argument('--id-a', type=int, help='particle id to read from log_A (single-particle mode)')
    ap.add_argument('--id-b', type=int, help='particle id to read from log_B (single-particle mode)')
    ap.add_argument('--diameter', type=float, default=0.015,
                    help='sphere diameter for the gap height h/d (default 0.015 m, ten Cate; DKT column: 1.0)')
    ap.add_argument('--plateau', type=float, nargs=2, metavar=('T0', 'T1'),
                    help='gate on the time mean of the (cloud-mean) settling velocity over [T0, T1] instead of the peak')
    a = ap.parse_args()
    radius, diameter = 0.5 * a.diameter, a.diameter

    tracks_a, src_a = parse_all(a.log_a)
    tracks_b, src_b = parse_all(a.log_b)
    print('%-9s: %d particle(s), source %s' % (a.label_a, len(tracks_a), src_a))
    print('%-9s: %d particle(s), source %s' % (a.label_b, len(tracks_b), src_b))

    single = (len(tracks_a) == 1 and len(tracks_b) == 1) or a.id_a is not None or a.id_b is not None
    if single:
        A = select(tracks_a, a.id_a, a.log_a)
        B = select(tracks_b, a.id_b, a.log_b)
        pairs = [(A, B)]
    else:
        if len(tracks_a) != len(tracks_b):
            sys.exit('particle counts differ: %d vs %d' % (len(tracks_a), len(tracks_b)))
        scale = max(1.0, max(np.abs(tr.p0).max() for tr in tracks_a.values()))
        pairs = match_by_start(tracks_a, tracks_b, scale)
        print('matched %d particle pairs by start position' % len(pairs))

    # common sample times, taken from the first pair (all particles of a run share them)
    ia, ib = align(pairs[0][0].t, pairs[0][1].t)
    if len(ia) == 0:
        sys.exit('no common sample times')
    t = pairs[0][0].t[ia]
    print('%-9s: %6d samples, %-9s: %6d samples, common: %d, t in [%.4f, %.4f] s'
          % (a.label_a, len(pairs[0][0].t), a.label_b, len(pairs[0][1].t), len(t), t[0], t[-1]))

    def cut(tr, idx):
        return tr.z[idx], tr.vz[idx]

    if single:
        A, B = pairs[0]
        za, va = cut(A, ia); zb, vb = cut(B, ib)
        vpa, tpa = peak(t, va); vpb, tpb = peak(t, vb)
        dpeak = (vpb / vpa - 1) * 100
        print('%-9s: v_peak = %+.6f m/s at t = %.4f s (id %d)' % (a.label_a, vpa, tpa, A.id))
        print('%-9s: v_peak = %+.6f m/s at t = %.4f s (id %d)' % (a.label_b, vpb, tpb, B.id))
        print('peak difference %s vs %s: %+.4f%% (velocity), %+.4f s (timing)' % (a.label_b, a.label_a, dpeak, tpb - tpa))
        if a.case:
            vpr = V_PEAK_PRINTED[a.case]
            print('Table II %s printed v_peak = %+.5f m/s: %s %+.3f%%, %s %+.3f%%'
                  % (a.case, vpr, a.label_a, (vpa / vpr - 1) * 100, a.label_b, (vpb / vpr - 1) * 100))
        report_windows(t, za, va, zb, vb, radius, diameter, a.label_a, a.label_b)
        gate_label = 'v_peak'
        plot_series = [(a.label_a, za, va), (a.label_b, zb, vb)]
    else:
        print('%-6s %-6s %14s %14s %10s %10s' % ('id_A', 'id_B', 'v_peak_A', 'v_peak_B', 'dpeak', 'dt_peak'))
        dpeaks = []
        VA, VB, ZA, ZB = [], [], [], []
        for A, B in pairs:
            za, va = cut(A, ia); zb, vb = cut(B, ib)
            vpa, tpa = peak(t, va); vpb, tpb = peak(t, vb)
            dp = (vpb / vpa - 1) * 100 if vpa != 0 else np.nan
            dpeaks.append(dp)
            print('%-6d %-6d %+14.6f %+14.6f %+9.4f%% %+9.4f s' % (A.id, B.id, vpa, vpb, dp, tpb - tpa))
            VA.append(va); VB.append(vb); ZA.append(za); ZB.append(zb)
        dpeaks = np.array(dpeaks)
        print('per-particle peak difference: mean %+.4f%%, max |.| %.4f%% (reported, not gated)'
              % (np.nanmean(dpeaks), np.nanmax(np.abs(dpeaks))))
        va, vb = np.mean(VA, axis=0), np.mean(VB, axis=0)
        za, zb = np.mean(ZA, axis=0), np.mean(ZB, axis=0)
        vpa, tpa = peak(t, va); vpb, tpb = peak(t, vb)
        dpeak = (vpb / vpa - 1) * 100
        print('cloud mean <v_z>: %s peak %+.6f m/s at t = %.4f s, %s peak %+.6f m/s at t = %.4f s'
              % (a.label_a, vpa, tpa, a.label_b, vpb, tpb))
        print('cloud-mean peak difference %s vs %s: %+.4f%% (velocity), %+.4f s (timing)'
              % (a.label_b, a.label_a, dpeak, tpb - tpa))
        # windows on the cloud mean; the "touchdown" is the lowest particle's
        zmin_a, zmin_b = np.min(ZA, axis=0), np.min(ZB, axis=0)
        report_windows(t, zmin_a, va, zmin_b, vb, radius, diameter, a.label_a, a.label_b)
        gate_label = 'cloud-mean v_peak'
        plot_series = [(a.label_a, za, va), (a.label_b, zb, vb)]

    dgate = dpeak
    if a.plateau:
        lo, hi = a.plateau
        m = (t >= lo) & (t <= hi)
        if not m.any():
            sys.exit('no samples in the plateau window [%g, %g]' % (lo, hi))
        ma, mb = va[m].mean(), vb[m].mean()
        dgate = (mb / ma - 1) * 100
        print('plateau mean of %s over t in [%.3f, %.3f] (%d samples): %s %+.6f m/s, %s %+.6f m/s, difference %+.4f%%'
              % ('<v_z>' if not single else 'v_z', lo, hi, int(m.sum()), a.label_a, ma, a.label_b, mb, dgate))
        gate_label = 'plateau mean %s' % ('<v_z>' if not single else 'v_z')
    verdict = 'PASS' if abs(dgate) <= a.tol_peak else 'FAIL'
    print('XMODE_GATE %s %s vs %s: %+.4f%% (tol %.2f%%) -> %s'
          % (gate_label, a.label_b, a.label_a, dgate, a.tol_peak, verdict))

    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
        styles = ['-', '--']
        for (lab, z, v), st in zip(plot_series, styles):
            ax[0].plot(t, v, st, label=lab)
            ax[1].plot(t, (z - radius) / diameter, st, label=lab)
        ax[0].set_ylabel('v_z [m/s]' if single else '<v_z> [m/s]'); ax[0].legend(); ax[0].grid(alpha=.3)
        ax[1].set_ylabel('h/d' if single else 'mean h/d'); ax[1].set_xlabel('t [s]'); ax[1].grid(alpha=.3)
        fig.suptitle('%s vs %s%s' % (a.label_a, a.label_b, ' (%s)' % a.case if a.case else ''))
        fig.tight_layout(); fig.savefig(a.plot, dpi=130)
        print('plot: %s' % a.plot)
    return 0 if verdict == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
