#!/usr/bin/env python3
"""D2.3 DKT pair-tilt analysis.

Parses particle_force.log (columns: time ip fx fy fz tx ty tz px py pz vx vy vz)
from a two-particle DKT rundir and reports, per run:

  - t_kiss: first time the center separation reaches contact (sep <= d + tol)
  - tilt(t): angle of the pair axis from vertical, deg, in the UNFOLDED
    0-180 convention (matches the site exporter dkt_export_series.py):
    the reference direction is particle 2 above particle 1, so tilt > 90
    means the pair has swung past horizontal and exchanged roles.
    (Versions before 2026-08-23 folded to 0-90: tilt_folded =
    min(tilt, 180 - tilt); datasheet rows written with the folded
    convention are annotated in place.)
  - tilt at kiss and at end of run
  - exponential growth rate of the tilt over a fit window (default final 3 t.u.)
    with the corresponding doubling time
  - end-state per-particle vz (a leader/trailer velocity split signals an
    active tumble; identical values signal a locked rigid pair)

Used for datasheet row d23_omegafix_rerun (frictional stall retired as an
hcaf_angvel_reset artifact; see practitioner's guide section 6, v2.2).

Usage:
  python3 tools/d23_tilt_analysis.py q2p1_dns_rundir_dkt_offset \
      q2p1_dns_rundir_dkt_offset_nofric q2p1_dns_rundir_dkt_offset_omegafix \
      [--diameter 1.0] [--fit-window 3.0] [--print-every 2.5]
"""
import argparse
import math
import os
import sys


def load(path):
    rows = {}
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#'):
                continue
            p = ln.split()
            if len(p) < 14:
                continue
            t = float(p[0])
            ip = int(p[1])
            rows.setdefault(t, {})[ip] = {
                'P': (float(p[8]), float(p[9]), float(p[10])),
                'V': (float(p[11]), float(p[12]), float(p[13])),
            }
    return {t: r for t, r in rows.items() if len(r) == 2}


def tilt_deg(r):
    (x1, y1, z1), (x2, y2, z2) = r[1]['P'], r[2]['P']
    return math.degrees(math.atan2(math.hypot(x2 - x1, y2 - y1), z2 - z1))


def separation(r):
    (x1, y1, z1), (x2, y2, z2) = r[1]['P'], r[2]['P']
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)


def analyze(path, d, fit_window, print_every):
    rows = load(os.path.join(path, 'particle_force.log'))
    ts = sorted(rows)
    if not ts:
        print(f'{path}: no two-particle records', file=sys.stderr)
        return
    t_kiss = next((t for t in ts if separation(rows[t]) <= d + 5e-4 * d), None)
    t_end = ts[-1]

    name = os.path.basename(path.rstrip('/'))
    print(f'== {name}')
    if print_every > 0:
        print('    t     tilt(deg)   sep-d')
        prev = None
        for t in ts:
            if prev is None or t - prev >= print_every - 1e-6:
                print(f'  {t:6.2f}  {tilt_deg(rows[t]):9.2f}  {separation(rows[t]) - d:8.4f}')
                prev = t

    kiss_txt = f'{t_kiss:.3f} (tilt {tilt_deg(rows[t_kiss]):.2f} deg)' if t_kiss else 'no contact'
    print(f'  t_kiss = {kiss_txt};  tilt(t={t_end:.2f}) = {tilt_deg(rows[t_end]):.2f} deg')

    sel = [t for t in ts if t >= t_end - fit_window and tilt_deg(rows[t]) > 0]
    if len(sel) > 2:
        xs, ys = sel, [math.log(tilt_deg(rows[t])) for t in sel]
        mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
        denom = sum((x - mx) ** 2 for x in xs)
        slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / denom
        dbl = math.log(2) / slope if slope > 1e-12 else float('inf')
        print(f'  ln-tilt growth over final {fit_window:g} t.u.: {slope:.3f}/t.u.'
              f'  (doubling {dbl:.1f} t.u.)')
    vz1, vz2 = rows[t_end][1]['V'][2], rows[t_end][2]['V'][2]
    print(f'  end vz: p1 {vz1:.4f}  p2 {vz2:.4f}  (split {abs(vz1 - vz2):.4f})')


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('rundirs', nargs='+')
    ap.add_argument('--diameter', type=float, default=1.0)
    ap.add_argument('--fit-window', type=float, default=3.0)
    ap.add_argument('--print-every', type=float, default=2.5,
                    help='tilt-table print interval in t.u.; 0 disables the table')
    args = ap.parse_args()
    for rd in args.rundirs:
        analyze(rd, args.diameter, args.fit_window, args.print_every)


if __name__ == '__main__':
    main()
