#!/usr/bin/env python3
"""Inlet velocity table of case 11 for SimPar@ELInletVelocityTable (rows "t U", superficial
velocity in m/s; EL_INLET_TABLE: linear between rows, constant beyond the ends).

Protocol of the spec (reference/inlet_table.py of case_specs/11_min_fluidization, same
parameters and defaults): U0 from t = 0 to TPREP + TSTEP (preparation, then step 0), then steps
k = 1..n, n = round((U0 - UMIN)/DU): a linear transition of TRANS seconds starting at
t_k = TPREP + k TSTEP from U0 - (k-1) DU to U0 - k DU, held until t_k + TSTEP. End
TPREP + (n+1) TSTEP (7.5 s with the defaults).

usage: make_inlet_table.py OUT [--u0 0.30] [--du 0.02] [--umin 0.04] [--tprep 0.5] [--tstep 0.5]
                               [--trans 0.02] [--const U] [--steps STEPS_OUT]
  --const U   one row "0 U": constant inlet velocity (probes)
  --steps F   also write the step list "t_start t_end U" (the evaluation windows of the post tool)
The rows equal the "csv" output of the spec's inlet_table.py to the printed digits (%.6g)."""
import argparse


def protocol(u0, du, umin, tprep, tstep, trans):
    n = int(round((u0 - umin) / du))
    pts = [(0.0, u0)]
    steps = [(tprep, tprep + tstep, u0)]
    for k in range(1, n + 1):
        t = tprep + k * tstep
        u = u0 - k * du
        pts += [(t, u0 - (k - 1) * du), (t + trans, u)]
        steps.append((t, t + tstep, u))
    tend = tprep + (n + 1) * tstep
    pts.append((tend, pts[-1][1]))
    return pts, steps, tend


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--u0', type=float, default=0.30)
    ap.add_argument('--du', type=float, default=0.02)
    ap.add_argument('--umin', type=float, default=0.04)
    ap.add_argument('--tprep', type=float, default=0.5)
    ap.add_argument('--tstep', type=float, default=0.5)
    ap.add_argument('--trans', type=float, default=0.02)
    ap.add_argument('--const', type=float, default=None)
    ap.add_argument('--steps', default=None)
    a = ap.parse_args()
    if a.const is not None:
        pts, steps, tend = [(0.0, a.const)], [], None
        head = f'# case 11 constant inlet U = {a.const:g} m/s (superficial)'
    else:
        pts, steps, tend = protocol(a.u0, a.du, a.umin, a.tprep, a.tstep, a.trans)
        head = (f'# case 11 protocol: U0 {a.u0:g} DU {a.du:g} UMIN {a.umin:g} TPREP {a.tprep:g} '
                f'TSTEP {a.tstep:g} TRANS {a.trans:g}; end {tend:g} s; t [s] U [m/s] (superficial)')
    with open(a.out, 'w') as f:
        f.write(head + '\n')
        for t, u in pts:
            f.write('%.6g %.6g\n' % (t, u))
    if a.steps:
        with open(a.steps, 'w') as f:
            f.write('# t_start t_end U\n')
            for s in steps:
                f.write('%.6g %.6g %.6g\n' % s)
    print(f'{a.out}: {len(pts)} rows' + (f', end {tend:g} s, {len(steps)} steps' if tend else ''))


if __name__ == '__main__':
    main()
