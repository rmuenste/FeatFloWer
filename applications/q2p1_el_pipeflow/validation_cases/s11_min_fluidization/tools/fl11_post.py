#!/usr/bin/env python3
"""Case 11 (minimum fluidisation) evaluation of a q2p1_el_pipeflow run.

Reads run_slurm.log of RUNDIR:
  EL_BED_AUDIT       (SimPar@ELBedAuditFreq > 0; EL_BED_AUDIT module): t, step, plane pressures,
                     dp, its hydrostatic part, buoyant weight, floor contact, bed heights,
                     window porosity
  EL_FORCE_BUDGET    (every step): drag, pressure force, grav_buoy summed over the particles
  EL_INLET_TABLE     (echo of SimPar@ELInletVelocityTable): U(t), the protocol of the run
  start-up echo      'EL kernel width factor' (mapping support = factor x d), 'EL fluid gravity',
                     'Global CFD domain extents' (z of the inlet plane); the deck
                     _data/q2p1_param.dat when the echo lacks the width factor
Per velocity step of the protocol (steps file of make_inlet_table.py --steps, or the plateaus
of the echoed table): mean and std of dp_dyn over the last 0.3 s and over the two sub-windows
0.2-0.35 s and 0.35-0.5 s of the step (flagged when they differ by more than --flag of the
plateau), the force-based dp (drag + pressure force)/A, bed height H = z99 + r, porosity of the
bed window. A step is complete when the audit samples its evaluation window densely (audit
interval at most a tenth of the window, no gap wider than one interval), not merely when the
log reaches its end.

Gates of the spec (section 5, as amended on 2026-11-20):
  F1  plateau dp A / buoyant weight, mean over the steps with dp* > 0.95           (1 +- 0.02)
  F2  fixed branch: dynamic window gradient -<dp/dz> of the audit (bed core; gradz_win_dyn, or
      for older logs gradz_win when the hydrostatic part is zero: dp_hydro = 0 on every audit
      line used) against the closure at the code's own MAPPED gas fraction (eps_kern)  (+-5 %),
      only when every contributing audit window excludes max(2 d, mapping support) at both
      ends (measured: lower end - inlet plane, z99 + r - upper end, and win_margin when the
      line has it; 2 d for binaries before 627caf47). Otherwise the ratios are printed as a
      measurement and the verdict is 'unavailable' with the reason. The ratio at the geometric
      window porosity (rule before the amendment), dp_dyn / H_eq and dp_dyn / H are ungated.
  F3  U_mf = where the fit dp* = a U + b U^2 of the fixed branch (through the origin; steps
      below 0.8 of the plateau) meets dp* = 1, against the closure U_mf at the WHOLE-BED porosity
      1 - N V_p / (A (z99 + r))                                                      (+-5 %);
      the crossing with the measured plateau and the ratios at the core porosity are reported
  F4  U_mf against the measurement 0.182 m/s (reported)
  F5  particle count constant over the run; no centre above mid-column during the plateau of the
      lowest velocity of the ECHOED protocol (the constant part; the 20 ms transition into it is
      reported separately and labelled). Evaluated independently of F1-F3, when the audit
      covers that plateau; a steps file cannot redefine the target.
Every ratio needs a positive finite denominator (area of the inlet plane, buoyant weight,
closure value); otherwise the quantity is printed as nan and the gates that need it are
'unavailable' instead of an exception.
The closures come from fl11_closures.py (shared with core_deficit.py); --closure gidaspow
(default), ergun, wenyu: dp/dz(U, eps) = n F_d(U/eps) / eps with the per-particle forms of
EL_FORCES (model A; identical to the spec's umf_reference.py for ergun and wenyu).

usage: fl11_post.py RUNDIR [--steps FILE] [--d 0.87e-3] [--rhop 860] [--rhof 1.0] [--mu 1.8e-5]
                    [--g 9.81] [--lz 0.10] [--closure gidaspow] [--flag 0.02] [--series OUT.csv]
       fl11_post.py --self-test       synthetic logs: missing inlet plane, empty bed, incomplete
                                      steps file, short window, margin / gradient rules
A run that did not reach the protocol end is evaluated over the steps it covers (partial steps
are listed but not used for the gates)."""
import argparse
import contextlib
import io
import math
import os
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fl11_closures import CLOSURES, closure_grad, closure_umf   # noqa: E402  (re-exported)

NUM = r'(?:[-+]?[0-9.]+(?:[EeDd][-+]?[0-9]+)?|[-+]?NaN)'     # NaN: a value the audit marks unavailable
NAN = float('nan')
TOL_T = 1e-9          # time comparisons (s)
TOL_Z = 2e-8          # height comparisons (m): print precision of the ES14.6 heights of the audit line


def f(x):
    return NAN if 'NaN' in x else float(x.replace('D', 'E').replace('d', 'e'))


def finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def pos(x):
    """x when it is a positive finite number (a usable denominator), else nan"""
    return x if finite(x) and x > 0 else NAN


def div(a, b):
    """a / b with b a positive finite number; nan otherwise (never an exception)"""
    b = pos(b)
    return a / b if finite(b) and finite(a) else NAN


def parse(rundir):
    audit, budget, table, pairs = [], {}, [], {}
    pat_a = re.compile(r'^EL_BED_AUDIT t=\s*(%s) step=\s*(\d+) p_in_out=\s*(%s)\s+(%s) area_in_out=\s*(%s)\s+(%s)'
                       r' dp_hydro_dyn=\s*(%s)\s+(%s)\s+(%s) n=\s*(\d+) weight_z_floor_z=\s*(%s)\s+(%s)'
                       r' z50_95_99_max=\s*(%s)\s+(%s)\s+(%s)\s+(%s) eps_bed=\s*(%s) window=\s*(%s)\s+(%s)'
                       r'(?: gradz_win=\s*(%s))?(?: floor_bal_z=\s*(%s))?(?: eps_kern_slip_win=\s*(%s)\s+(%s))?'
                       r'(?: win_margin=\s*(%s))?(?: gradz_win_dyn=\s*(%s))?'
                       % ((NUM,) + (NUM,) * 7 + (NUM,) * 2 + (NUM,) * 4 + (NUM,) * 3 + (NUM,) * 4 + (NUM,) * 2))
    pat_b = re.compile(r'^EL_FORCE_BUDGET step=\s*(\d+) drag=\s*(%s)\s+(%s)\s+(%s) pressure=\s*(%s)\s+(%s)\s+(%s)'
                       r' lift=.* grav_buoy=\s*(%s)\s+(%s)\s+(%s)' % ((NUM,) * 9))
    for line in open(os.path.join(rundir, 'run_slurm.log'), errors='replace'):
        if line.startswith('EL_BED_AUDIT'):
            m = pat_a.match(line)
            if not m:
                continue
            g = m.groups()
            audit.append(dict(t=f(g[0]), step=int(g[1]), p_in=f(g[2]), p_out=f(g[3]), a_in=f(g[4]),
                              a_out=f(g[5]), dp=f(g[6]), dp_h=f(g[7]), dp_dyn=f(g[8]), n=int(g[9]),
                              wz=f(g[10]), floor=f(g[11]), z50=f(g[12]), z95=f(g[13]), z99=f(g[14]),
                              zmax=f(g[15]), eps=f(g[16]), zlo=f(g[17]), zhi=f(g[18]),
                              gw=f(g[19]) if g[19] else NAN,
                              fbal=f(g[20]) if g[20] else None,
                              ekern=f(g[21]) if g[21] else NAN,
                              skern=f(g[22]) if g[22] else NAN,
                              margin=f(g[23]) if g[23] else None,          # window margin (from 627caf47 on)
                              gw_dyn=f(g[24]) if g[24] else NAN,
                              has_gw_dyn=bool(g[24])))                      # field written (from 627caf47 on)
        elif line.startswith('EL_FORCE_BUDGET'):
            m = pat_b.match(line)
            if m:
                g = m.groups()
                budget[int(g[0])] = (f(g[3]), f(g[6]), f(g[9]))
        elif line.startswith('EL_NEWTON_PAIR'):
            c = line.split()
            pairs[int(c[4])] = (f(c[8]), f(c[12]))       # dp_z, expected_z of the particle momentum
        elif line.lstrip().startswith('EL_INLET_TABLE'):
            c = line.split()
            table.append((f(c[1]), f(c[2])))
    # floor force from the particle momentum balance (EL_BED_AUDIT floor_bal_z; older logs: from the
    # EL_NEWTON_PAIR line of the same step, which carries the same two numbers)
    if len(audit) > 1:
        dt = (audit[1]['t'] - audit[0]['t'])/(audit[1]['step'] - audit[0]['step'])
        for row in audit:
            if row['fbal'] is None:
                pr = pairs.get(row['step'])
                row['fbal'] = (pr[0] - pr[1])/dt if pr else NAN
    return audit, budget, table


def parse_setup(rundir):
    """kernel width factor (and where it came from), fluid gravity, z of the inlet plane"""
    setup = dict(width=None, width_src=None, fluid_gravity=None, zmin=None)
    pat_w = re.compile(r'^EL kernel width factor\s*=\s*(%s)' % NUM)
    pat_g = re.compile(r'^EL fluid gravity\s*=\s*([TF])')
    pat_x = re.compile(r'^Global CFD domain extents:\s*(%s)\s+(%s)\s+(%s)\s+(%s)\s+(%s)\s+(%s)' % ((NUM,) * 6))
    for line in open(os.path.join(rundir, 'run_slurm.log'), errors='replace'):
        if line.startswith('EL_'):
            break                                         # the start-up echo precedes every EL_ stream
        m = pat_w.match(line)
        if m and setup['width'] is None:
            setup['width'], setup['width_src'] = f(m.group(1)), 'log echo'
        m = pat_g.match(line)
        if m and setup['fluid_gravity'] is None:
            setup['fluid_gravity'] = (m.group(1) == 'T')
        m = pat_x.match(line)
        if m and setup['zmin'] is None:
            setup['zmin'] = f(m.group(5))
    if setup['width'] is None:
        deck = os.path.join(rundir, '_data', 'q2p1_param.dat')
        if os.path.exists(deck):
            for line in open(deck, errors='replace'):
                m = re.match(r'^\s*SimPar@ELKernelWidthFactor\s*=\s*(%s)' % NUM, line)
                if m:
                    setup['width'], setup['width_src'] = f(m.group(1)), 'deck'
        if setup['width'] is None:
            setup['width'], setup['width_src'] = 2.5, 'EL_CONFIG default'
    return setup


def u_of_t(table, t):
    if t <= table[0][0]:
        return table[0][1]
    for (t0, u0), (t1, u1) in zip(table, table[1:]):
        if t <= t1:
            return u0 + (u1 - u0) * (t - t0) / (t1 - t0)
    return table[-1][1]


def plateaus(table):
    """steps from the echoed table: intervals of constant U. A step starts where the transition
    into it starts (the convention of the steps file of make_inlet_table.py), so the sub-windows
    0.2-0.35 s and 0.35-0.5 s do not depend on the source of the step list. The preparation is
    merged with the first plateau; main() evaluates the last 0.5 s of it."""
    out = []
    for i, ((t0, u0), (t1, u1)) in enumerate(zip(table, table[1:])):
        if u0 == u1 and t1 - t0 > 0.1:
            ramp = i > 0 and table[i - 1][1] != u0
            out.append((table[i - 1][0] if ramp else t0, t1, u0))
    return out


def lowest_plateau(table):
    """(u, [(t_a, t_b)], [(t_transition_start, t_a)]): the lowest velocity held constant in the
    echoed protocol, its constant intervals and the transitions into them"""
    held = []
    for i, ((t0, u0), (t1, u1)) in enumerate(zip(table, table[1:])):
        if u0 == u1 and t1 - t0 > 0.1:
            held.append((u0, t0, t1, table[i - 1][0] if i > 0 and table[i - 1][1] != u0 else t0))
    if not held:
        return None
    umin = min(h[0] for h in held)
    return (umin, [(h[1], h[2]) for h in held if h[0] == umin], [(h[3], h[1]) for h in held if h[0] == umin])


def stats(rows, key, t0, t1, open_end=False):
    """mean, std, count of rows[key] over t0 <= t <= t1 (t < t1 with open_end: disjoint windows)"""
    eps = 1e-9
    v = [r[key] for r in rows if t0 - eps <= r['t'] and (r['t'] < t1 - eps if open_end else r['t'] <= t1 + eps)]
    if not v:
        return NAN, NAN, 0
    m = sum(v) / len(v)
    s = math.sqrt(sum((x - m) ** 2 for x in v) / len(v))
    return m, s, len(v)


def audit_interval(rows):
    """the audit sampling interval (median spacing of the audit times); nan below two lines"""
    d = sorted(b['t'] - a['t'] for a, b in zip(rows, rows[1:]) if b['t'] > a['t'])
    return d[len(d) // 2] if d else NAN


def covered(rows, t0, t1, dt_a):
    """the audit samples [t0, t1] densely: an audit interval of at most a tenth of the window (ten
    samples or more) and no gap wider than one interval (half an interval of slack), both ends
    included"""
    if not finite(dt_a) or t1 <= t0 or dt_a > (t1 - t0) / 10:
        return False
    ts = [r['t'] for r in rows if t0 - TOL_T <= r['t'] <= t1 + TOL_T]
    if not ts:
        return False
    pts = [t0] + ts + [t1]
    return all(b - a <= 1.5 * dt_a + TOL_T for a, b in zip(pts, pts[1:]))


def window(rows, t0, t1):
    return [r for r in rows if t0 - TOL_T <= r['t'] <= t1 + TOL_T]


def linfit(x, y):
    n = len(x)
    if n < 2:
        return NAN, NAN
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if sxx == 0:
        return NAN, NAN
    s = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx
    return s, my - s * mx


def read_steps(path):
    return [tuple(map(float, l.split()[:3])) for l in open(path) if l.strip() and not l.lstrip().startswith('#')]


def match_steps(listed, protocol):
    """steps of the list that are not plateaus of the echoed protocol, and plateaus of the
    protocol missing from the list (matched on the end time and the velocity)"""
    def same(a, b):
        return abs(a[1] - b[1]) <= 1e-6 and abs(a[2] - b[2]) <= 1e-6 * max(1.0, abs(b[2]))
    foreign = [s for s in listed if not any(same(s, p) for p in protocol)]
    missing = [p for p in protocol if not any(same(s, p) for s in listed)]
    return foreign, missing


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir', nargs='?')
    ap.add_argument('--steps', default=None)
    ap.add_argument('--d', type=float, default=0.87e-3)
    ap.add_argument('--rhop', type=float, default=860.0)
    ap.add_argument('--rhof', type=float, default=1.0)
    ap.add_argument('--mu', type=float, default=1.8e-5)
    ap.add_argument('--g', type=float, default=9.81)
    ap.add_argument('--lz', type=float, default=0.10)
    ap.add_argument('--closure', default='gidaspow', choices=CLOSURES)
    ap.add_argument('--flag', type=float, default=0.02)
    ap.add_argument('--series', default=None)
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.rundir:
        ap.error('RUNDIR is required')

    audit, budget, table = parse(a.rundir)
    if not audit:
        sys.exit('no EL_BED_AUDIT lines (SimPar@ELBedAuditFreq not set?)')
    if not table:
        sys.exit('no EL_INLET_TABLE echo (SimPar@ELInletVelocityTable not set?)')
    setup = parse_setup(a.rundir)
    r = a.d / 2
    vp = math.pi * a.d ** 3 / 6
    area = pos(audit[-1]['a_in'])          # nan without an inlet plane (no element face at zmin)
    for row in audit:
        row['U'] = u_of_t(table, row['t'])
        fb = budget.get(row['step'])
        row['dp_force'] = div(fb[0] + fb[1], area) if fb else NAN
        row['dp_drag'] = div(fb[0], area) if fb else NAN
        row['H'] = row['z99'] + r
    W = pos(-audit[-1]['wz'])              # nan for an empty bed
    wa = div(W, area)
    n0 = audit[0]['n']
    tend = audit[-1]['t']
    dt_a = audit_interval(audit)
    print(f"{os.path.basename(os.path.abspath(a.rundir))}: {len(audit)} audit lines, t {audit[0]['t']:.4f} .. {tend:.4f} s, "
          f"N {n0}, A {area:.6e} m2, buoyant weight {W:.6e} N, W/A {wa:.3f} Pa, "
          f"hydrostatic part in dp {audit[-1]['dp_h']:.4f} Pa (subtracted: dp_dyn)")
    if not finite(wa):
        print('  W/A unavailable (' + ', '.join(x for x, ok in (('no inlet plane: area_in not positive', finite(area)),
                                                                ('empty bed: no buoyant weight', finite(W))) if not ok)
              + '): dp*, F1, F2, F3 unavailable')
    if a.series:
        with open(a.series, 'w') as fh:
            fh.write('t,step,U,dp,dp_hydro,dp_dyn,dp_force,dp_drag,floor_z,z50,z95,z99,zmax,eps_bed,n\n')
            for row in audit:
                fh.write(','.join(str(row[k]) for k in ('t', 'step', 'U', 'dp', 'dp_h', 'dp_dyn', 'dp_force',
                                                        'dp_drag', 'floor', 'z50', 'z95', 'z99', 'zmax', 'eps', 'n')) + '\n')

    # ---- the protocol of the run (echo) and the step list -----------------------------------------
    protocol = plateaus(table)
    if protocol and protocol[0][1] - protocol[0][0] > 0.6:
        # the first plateau of the protocol holds preparation + step 0: evaluate its last 0.5 s
        protocol[0] = (protocol[0][1] - 0.5, protocol[0][1], protocol[0][2])
    if a.steps:
        steps = read_steps(a.steps)
        foreign, missing = match_steps(steps, protocol)
        print(f"\nstep list {a.steps}: {len(steps)} steps; echoed protocol {len(protocol)} plateaus"
              + ('' if not (foreign or missing) else
                 "; NOT the echoed protocol:"
                 + (f" {len(missing)} plateau(s) missing from the list {[round(p[2], 4) for p in missing]}" if missing else '')
                 + (f" {len(foreign)} listed step(s) not in the echo {[round(s[2], 4) for s in foreign]} "
                    "(excluded from the gates)" if foreign else '')))
        steps = [s for s in steps if s not in foreign]
    else:
        steps = protocol
        missing = []
        print(f"\nstep list: the {len(protocol)} plateaus of the echoed protocol")
    proto_end = max((p[1] for p in protocol), default=NAN)
    print(f"protocol end {proto_end:.4f} s, audit to {tend:.4f} s, audit interval {dt_a:.6f} s"
          + ('' if finite(proto_end) and tend >= proto_end - TOL_T else ': the run has not reached the protocol end'))

    print('\n step  U [m/s]  window      dp*      std     sub1    sub2   flag  F/(W/A)  H [mm]  eps_bed  floor/W  dp*+floor  (floor_pe/W)')
    res = []
    for k, (t0, t1, u) in enumerate(steps):
        lo = min(t1 - 0.3, t0 + 0.2)
        complete = covered(audit, lo, t1, dt_a)
        m, s, cnt = stats(audit, 'dp_dyn', t1 - 0.3, t1)
        m1, _, _ = stats(audit, 'dp_dyn', t0 + 0.2, t0 + 0.35, open_end=True)
        m2, _, _ = stats(audit, 'dp_dyn', t0 + 0.35, t0 + 0.5)
        mf, _, _ = stats(audit, 'dp_force', t1 - 0.3, t1)
        h, _, _ = stats(audit, 'H', t1 - 0.3, t1)
        e, _, _ = stats(audit, 'eps', t1 - 0.3, t1)
        fl, _, _ = stats(audit, 'fbal', t1 - 0.3, t1)
        flpe, _, _ = stats(audit, 'floor', t1 - 0.3, t1)
        flag = finite(wa) and abs(m1 - m2) > a.flag * wa
        res.append(dict(k=k, t0=t0, t1=t1, U=u, dp=m, std=s, cnt=cnt, eps=e, H=h, fl=div(fl, W),
                        complete=complete and cnt > 0))
        print(f" {k:3d}  {u:6.3f}  {t1-0.3:5.2f}-{t1:5.2f}  {div(m, wa):7.4f}  {div(s, wa):7.4f}  {div(m1, wa):6.3f}  {div(m2, wa):6.3f}  "
              f"{'DIFF' if flag else '    '}  {div(mf, wa):7.4f}  {h*1e3:6.2f}  {e:7.4f}  {div(fl, W):7.4f}  {div(m, wa) + div(fl, W):8.4f}   ({div(flpe, W):.4f})"
              f"{'' if complete else '  (incomplete)'}")
    print('  dp* = dp_dyn A / W over 0.2-0.5 s of the step; floor from the particle momentum balance; floor_pe = the pe '
          'wall-contact impulse (unrelaxed PGS impulse, 1/0.9 of the applied one)')

    gates_f1_f3(a, res, audit, setup, wa, area, n0, vp, r, dt_a, missing)
    gate_f5(a, audit, table, dt_a)


def gates_f1_f3(a, res, audit, setup, wa, area, n0, vp, r, dt_a, missing):
    """F1, F3, F4 and F2 (each 'unavailable' with the reason when its inputs are not there)"""
    use = [x for x in res if x['complete'] and finite(x['dp'])]
    if missing:
        print(f"\nnote: the step list lacks {len(missing)} plateau(s) of the echoed protocol: F1-F3 use the listed steps only")
    if not finite(wa):
        print('\nF1, F2, F3 unavailable: W/A has no positive finite value')
        return
    if len(use) < 3:
        print(f'\nF1, F2, F3 unavailable: {len(use)} complete step(s) with a measured dp, three needed')
        return
    # ---- F1 (spec rule of 2026-11-20): mean dp* over all steps with dp* > 0.95 ------------------
    flu = [x for x in use if x['dp'] > 0.95 * wa]
    if flu:
        plateau = sum(x['dp'] for x in flu) / len(flu)
        print(f"\nF1 plateau: mean dp* over {len(flu)} steps with dp* > 0.95 (U {min(x['U'] for x in flu):.2f} .. "
              f"{max(x['U'] for x in flu):.2f}): {plateau / wa:.6f}   (gate 1 +- 0.02: {'PASS' if abs(plateau / wa - 1) <= 0.02 else 'FAIL'}); "
              f"dp* + floor/W over these {sum(x['dp']/wa + x['fl'] for x in flu)/len(flu):.4f}")
    else:
        plateau = wa
        print('\nF1 not evaluated: no step with dp* > 0.95; plateau taken as W/A')
    fixed = [x for x in use if x['dp'] < 0.8 * plateau]
    if len(fixed) < 2:
        print(f'F2, F3 unavailable: {len(fixed)} fixed-branch step(s) (dp* < 0.8 of the plateau), two needed')
        return
    efix = sum(x['eps'] for x in fixed) / len(fixed)
    uc = closure_umf(efix, a.d, a.rhop, a.rhof, a.mu, a.g, a.closure)
    # whole-bed porosity 1 - N V_p / (A (z99 + r)): the reference porosity of F3 (spec, 2026-11-20)
    ebed = sum(1 - div(n0 * vp, area * x['H']) for x in fixed) / len(fixed)
    ub = closure_umf(ebed, a.d, a.rhop, a.rhof, a.mu, a.g, a.closure)
    # ---- F3 (agreed with the peer): Ergun form dp* = a U + b U^2 through the origin over the fixed steps
    su2 = sum(x['U']**2 for x in fixed); su3 = sum(x['U']**3 for x in fixed); su4 = sum(x['U']**4 for x in fixed)
    sy1 = sum(x['U']*x['dp']/wa for x in fixed); sy2 = sum(x['U']**2*x['dp']/wa for x in fixed)
    det = su2*su4 - su3*su3
    ca = (sy1*su4 - sy2*su3)/det if det != 0 else NAN
    cb = (su2*sy2 - su3*sy1)/det if det != 0 else NAN

    def cross(level):
        # positive root of b U^2 + a U = level; nan when the fit has no positive crossing
        if not (finite(ca) and finite(cb)): return NAN
        if abs(cb) < 1e-14: return level/ca if ca > 0 else NAN
        disc = ca*ca + 4*cb*level
        if disc < 0: return NAN
        roots = [x for x in ((-ca + math.sqrt(disc))/(2*cb), (-ca - math.sqrt(disc))/(2*cb)) if x > 0]
        return min(roots) if roots else NAN
    u1 = cross(1.0); up = cross(plateau/wa)
    r3 = div(u1, ub)
    print(f"F3 fixed branch {len(fixed)} steps (U {min(x['U'] for x in fixed):.3f} .. {max(x['U'] for x in fixed):.3f}): "
          f"dp* = {ca:.4f} U + {cb:.4f} U^2 (through the origin)")
    print(f"   U_mf at dp* = 1: {u1:.4f} m/s; {a.closure} U_mf at the whole-bed porosity {ebed:.4f}: {ub:.4f} m/s; "
          f"ratio {r3:.4f}   (gate 0.05: {('PASS' if abs(r3-1) <= 0.05 else 'FAIL') if finite(r3) else 'unavailable'})")
    print(f"   reported, not gated: U_mf at the measured plateau {plateau/wa:.4f}: {up:.4f} m/s = {div(up, ub):.4f} of the closure; "
          f"at the core porosity {efix:.4f} the closure gives {uc:.4f} m/s: ratios {div(u1, uc):.4f} / {div(up, uc):.4f}")
    print(f"F4 U_mf / 0.182 (NETL measurement): {u1/0.182:.4f} (dp* = 1) / {up/0.182:.4f} (plateau) (reported, not gated)")
    # ---- secondary (rules before 2026-11-20): F1 over U > 1.2 U_mf of the linear fit, linear-fit U_mf
    pl2 = wa
    umf2 = NAN
    fl2 = []
    for _ in range(3):
        fx2 = [x for x in use if x['dp'] < 0.8 * pl2]
        if len(fx2) >= 2:
            s2, i2 = linfit([x['U'] for x in fx2], [x['dp'] for x in fx2])
            umf2 = (pl2 - i2) / s2 if finite(s2) and s2 else NAN
        fl2 = [x for x in use if (not math.isnan(umf2)) and x['U'] > 1.2 * umf2] or \
              [x for x in use if x['dp'] >= 0.8 * pl2]
        if fl2:
            pl2 = sum(x['dp'] for x in fl2) / len(fl2)
    s, i = linfit([x['U'] for x in fixed], [x['dp'] for x in fixed])
    print(f"   secondary: F1 over the {len(fl2)} steps with U > 1.2 U_mf(linear) {pl2/wa:.4f}; linear fit slope {s:.1f} Pa s/m "
          f"(intercept {i:.2f} Pa), U_mf(linear) {umf2:.4f} m/s = {div(umf2, uc):.4f} of the closure")
    gate_f2(a, fixed, audit, setup, area, n0, vp, r)


def gate_f2(a, fixed, audit, setup, area, n0, vp, r):
    """F2 in three readings: (k) the dynamic window gradient of the audit against the closure at the
    mapped gas fraction (the gate), (w) the same gradient at the geometric window porosity, (e)
    dp_dyn / H_eq with H_eq the height of a uniform bed at eps_bed, (h) dp_dyn / H with H = z99 + r
    against the closure at the mean porosity of that height (end layers included)"""
    reasons = []
    # -- which gradient: the dynamic one, or the total one where the hydrostatic part is zero
    rows_all = [row for x in fixed for row in window(audit, x['t1'] - 0.3, x['t1'])]
    if rows_all and all(row['has_gw_dyn'] for row in rows_all):
        gkey, gtext = 'gw_dyn', 'gradz_win_dyn (dynamic window gradient)'
    elif rows_all and all(row['dp_h'] == 0.0 for row in rows_all):
        gkey = 'gw'
        gtext = ('gradz_win (log before 627caf47 without gradz_win_dyn; hydrostatic part zero: dp_hydro = 0 on all '
                 f"{len(rows_all)} audit lines used" + ('' if setup['fluid_gravity'] is None else
                                                         f", start-up echo 'EL fluid gravity = {'T' if setup['fluid_gravity'] else 'F'}'") + ')')
    else:
        gkey, gtext = 'gw', 'gradz_win (total gradient)'
        reasons.append('no dynamic window gradient (gradz_win_dyn absent) and a nonzero hydrostatic part (dp_hydro)')
    # -- the window margins of every contributing audit line
    required = max(2 * a.d, setup['width'] * a.d)
    tol = max(1e-6 * required, TOL_Z)
    lows, highs, fields = [], [], []
    for row in rows_all:
        if setup['zmin'] is not None:
            lows.append(row['zlo'] - setup['zmin'])
        highs.append(row['z99'] + r - row['zhi'])
        fields.append(row['margin'] if row['margin'] is not None else 2 * a.d)   # binaries before 627caf47: 2 d
    measured = [x for x in lows + highs + fields if x == x]
    mmin = min(measured) if measured and len(measured) == len(lows) + len(highs) + len(fields) else NAN
    margin_ok = finite(mmin) and mmin >= required - tol
    if not margin_ok:
        reasons.append('window inside the mapping support')
    print(f" F2: gradient {gtext}")
    print(f"     audit window margin over {len(rows_all)} lines: min {div(mmin, a.d):.3f} d "
          f"(lower end - inlet plane min {min(lows, default=NAN)/a.d:.3f} d, z99 + r - upper end min {min(highs, default=NAN)/a.d:.3f} d, "
          f"{'win_margin field' if all(row['margin'] is not None for row in rows_all) else 'margin 2 d (binary before 627caf47)'} "
          f"min {min(fields, default=NAN)/a.d:.3f} d); required max(2 d, mapping support {setup['width']:g} d "
          f"[{setup['width_src']}]) = {required/a.d:.3f} d: {'met' if margin_ok else 'NOT met'}")
    print(' F2: step  U      eps_bed  window-grad  closure(eps_bed)  ratio_w | eps_kern ratio_k | dp/H_eq  ratio_e | eps_H   dp/H   closure(eps_H)  ratio_h')
    worst_k, worst_w = 0.0, 0.0
    for x in fixed:
        x['g'], _, _ = stats(audit, gkey, x['t1'] - 0.3, x['t1'])
        x['ek'], _, _ = stats(audit, 'ekern', x['t1'] - 0.3, x['t1'])
    have_k = all(finite(x['g']) and finite(x['ek']) and x['ek'] > 0 for x in fixed)
    have_w = all(finite(x['g']) for x in fixed)
    if not have_k:
        reasons.append('no window gradient and mapped gas fraction on every fixed step (gradz_win, eps_kern_slip_win)')
    for x in fixed:
        heq = div(n0 * vp, area * (1 - x['eps']))
        clo = closure_grad(x['U'], x['eps'], a.d, a.rhof, a.mu, a.closure)
        eh = 1 - div(n0 * vp, area * x['H'])
        cloh = closure_grad(x['U'], eh, a.d, a.rhof, a.mu, a.closure)
        rw = div(x['g'], clo)
        rk = div(x['g'], closure_grad(x['U'], x['ek'], a.d, a.rhof, a.mu, a.closure)) if finite(x['ek']) and x['ek'] > 0 else NAN
        re_ = div(div(x['dp'], heq), clo)
        rh = div(div(x['dp'], x['H']), cloh)
        worst_k = max(worst_k, abs(rk - 1)) if have_k else NAN
        worst_w = max(worst_w, abs(rw - 1)) if have_w else NAN
        print(f"     {x['k']:3d}  {x['U']:.3f}  {x['eps']:.4f}  {x['g']:10.1f}  {clo:10.1f}  {rw:12.4f} | {x['ek']:.4f} {rk:7.4f} | {div(x['dp'], heq):8.1f} {re_:7.4f} | "
              f"{eh:.4f} {div(x['dp'], x['H']):8.1f} {cloh:10.1f}  {rh:9.4f}")
    if not reasons:
        print(f"F2 worst |ratio_k - 1| (window gradient against the closure at the mapped gas fraction) {worst_k:.4f}   "
              f"(gate 0.05: {'PASS' if worst_k <= 0.05 else 'FAIL'})")
    else:
        print(f"F2 measured, not gated: worst |ratio_k - 1| {worst_k:.4f}")
        print(f"F2 verdict unavailable ({'; '.join(reasons)})")
    if have_w:
        print(f"   reported, not gated: worst |ratio_w - 1| at the geometric window porosity (rule before 2026-11-20) {worst_w:.4f}")


def gate_f5(a, audit, table, dt_a):
    """F5 over the lowest plateau of the echoed protocol, independent of the other gates"""
    nset = sorted({row['n'] for row in audit})
    low = lowest_plateau(table)
    if low is None:
        print(f"F5 not evaluated: the echoed table holds no velocity plateau (particle counts seen {nset})")
        return
    u, held, trans = low
    rows = [row for ta, tb in held for row in window(audit, ta, tb)]
    cov = all(covered(audit, ta, tb, dt_a) for ta, tb in held)
    span = ', '.join(f'{ta:.3f}-{tb:.3f} s' for ta, tb in held)
    if not cov:
        have = window(audit, held[0][0], held[-1][1])
        print(f"F5 not evaluated: the audit does not cover the plateau of the lowest protocol velocity {u:.3f} m/s "
              f"({span}; {len(have)} audit lines in it, audit to {audit[-1]['t']:.4f} s; particle counts seen so far {nset})")
        return
    zs = [row['zmax'] for row in rows]
    if not all(finite(z) for z in zs) or nset == [0]:
        print(f"F5 unavailable: no particle centre height on every audit line of the lowest plateau (particle counts {nset})")
        return
    zm = max(zs)
    rows_t = [row for ta, tb in trans for row in window(audit, ta, tb)]
    zt = max((row['zmax'] for row in rows_t if finite(row['zmax'])), default=NAN)
    ok = len(nset) == 1 and nset[0] > 0 and zm < a.lz / 2
    print(f"F5 particle counts seen {nset}; highest centre during the plateau of the lowest velocity "
          f"{u:.3f} m/s ({span}, {len(rows)} audit lines): {zm*1e3:.2f} mm (mid-column {a.lz*500:.1f} mm): "
          f"{'PASS' if ok else 'FAIL'}")
    print(f"   stricter, reported: including the {1e3*sum(tb - ta for ta, tb in trans):.0f} ms transition into it "
          f"({', '.join(f'{ta:.3f}-{tb:.3f} s' for ta, tb in trans)}): highest centre {max(zm, zt)*1e3:.2f} mm")


# ---------------------------------------------------------------------------------------------------
# self-test: synthetic logs
# ---------------------------------------------------------------------------------------------------
def _e(x, w=14, p=6):
    return ('%*s' % (w, 'NaN')) if x != x else '%*.*E' % (w, p, x)


def _write_run(path, protocol, dt_a=0.01, t_end=None, area=2.25e-4, n=9789, width=2.5, margin_d=2.5,
               zlo_d=None, dyn=True, dp_hydro=0.0, d=0.87e-3, steps_file=None):
    """a synthetic run directory: start-up echo, table echo, EL_BED_AUDIT and EL_FORCE_BUDGET lines.
    dp* = min(1, U/0.18) (fixed bed below 0.18 m/s), bed top 24 mm, mapped porosity 0.37"""
    os.makedirs(os.path.join(path, '_data'), exist_ok=True)
    rhop, rhof, g = 860.0, 1.0, 9.81
    vp = math.pi * d ** 3 / 6
    W = n * (rhop - rhof) * vp * g
    tab = [(0.0, protocol[0][1])]
    for (t, u) in protocol[1:]:
        tab += [(t, tab[-1][1]), (t + 0.02, u)]
    t_proto = protocol[-1][0] + 0.5
    tab.append((t_proto, tab[-1][1]))
    t_end = t_proto if t_end is None else t_end
    lines = ['EL kernel                 = deen_poly', 'EL kernel width factor    = ' + _e(width),
             'EL fluid gravity          = ' + ('T' if dp_hydro else 'F')]
    lines += ['  EL_INLET_TABLE %s %s' % (_e(t, 16, 8), _e(u, 16, 8)) for t, u in tab]
    lines += ['Global CFD domain extents: ' + ' '.join(_e(x) for x in (0.0, 0.015, 0.0, 0.015, 0.0, 0.1))]
    k = 1
    while k * dt_a <= t_end + 1e-12:
        t = k * dt_a
        u = u_of_t(tab, t)
        dps = min(1.0, u / 0.18)
        dp = dps * div(W, area) + dp_hydro if (area > 0 and n > 0) else NAN
        z99 = 0.024 if n > 0 else NAN
        m = margin_d * d
        zlo = (zlo_d if zlo_d is not None else margin_d) * d if n > 0 else NAN
        zhi = z99 + d / 2 - m if n > 0 else NAN
        grad = dp / (0.1) if dp == dp else NAN
        gdyn = grad - dp_hydro / 0.1 if dp == dp else NAN
        line = ('EL_BED_AUDIT t= %s step= %d p_in_out= %s %s area_in_out= %s %s dp_hydro_dyn= %s %s %s n= %d '
                'weight_z_floor_z= %s %s z50_95_99_max= %s %s %s %s eps_bed= %s window= %s %s gradz_win= %s '
                'floor_bal_z= %s eps_kern_slip_win= %s %s' % (
                    _e(t), k * 40, _e(dp, 16, 8), _e(0.0, 16, 8), _e(area), _e(area), _e(dp, 16, 8),
                    _e(dp_hydro, 16, 8), _e(dp - dp_hydro if dp == dp else NAN, 16, 8), n,
                    _e(-W if n > 0 else 0.0), _e(0.0), _e(z99 / 2), _e(z99 * 0.96), _e(z99), _e(z99 + 1e-4 if n > 0 else NAN),
                    _e(0.37 if n > 0 else NAN), _e(zlo), _e(zhi), _e(grad, 16, 8), _e((1 - dps) * W, 16, 8),
                    _e(0.37 if n > 0 else NAN, 16, 8), _e(0.1 if n > 0 else NAN, 16, 8)))
        if margin_d is not None and dyn:
            line += ' win_margin= %s gradz_win_dyn= %s' % (_e(m if n > 0 else NAN), _e(gdyn, 16, 8))
        lines.append(line)
        lines.append('EL_FORCE_BUDGET step= %d drag= %s %s %s pressure= %s %s %s lift= 0 0 0 grav_buoy= %s %s %s' % (
            k * 40, _e(0.0), _e(0.0), _e(dps * W * 0.9), _e(0.0), _e(0.0), _e(dps * W * 0.1), _e(0.0), _e(0.0), _e(-W)))
        k += 1
    with open(os.path.join(path, 'run_slurm.log'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    if steps_file is not None:
        with open(os.path.join(path, '_data', 'steps.txt'), 'w') as fh:
            fh.write('# t_start t_end U\n' + ''.join('%.6g %.6g %.6g\n' % s for s in steps_file))
    return path


def self_test():
    proto = [(0.5, 0.30), (1.0, 0.20), (1.5, 0.10), (2.0, 0.05)]          # (start of step, U); end 2.5 s
    full_steps = [(0.5, 1.0, 0.30), (1.0, 1.5, 0.20), (1.5, 2.0, 0.10), (2.0, 2.5, 0.05)]
    cases = [
        ('reference (all gates measurable)', dict(), [], ['F1 plateau', 'F3 fixed branch', 'F2 worst |ratio_k - 1|',
                                                          'F5 particle counts seen [9789]', 'PASS', 'stricter, reported'],
         ['unavailable', 'Traceback']),
        ('missing inlet plane (area_in = 0)', dict(area=0.0), [], ['W/A unavailable (no inlet plane', 'F1, F2, F3 unavailable',
                                                                   'F5 particle counts seen'], ['Traceback']),
        ('empty bed (n = 0)', dict(n=0), [], ['empty bed', 'F1, F2, F3 unavailable', 'F5 unavailable'], ['Traceback']),
        ('incomplete steps file (lowest step missing)', dict(steps_file=full_steps[:3]), ['--steps', '{run}/_data/steps.txt'],
         ['1 plateau(s) missing from the list [0.05]', 'lowest velocity 0.050 m/s', 'lacks 1 plateau'], ['Traceback']),
        ('foreign step in the steps file', dict(steps_file=full_steps + [(2.5, 3.0, 0.01)]), ['--steps', '{run}/_data/steps.txt'],
         ['1 listed step(s) not in the echo [0.01]', 'lowest velocity 0.050 m/s'], ['Traceback', '0.010 m/s']),
        ('short window (run ends inside the lowest plateau)', dict(t_end=2.3), [],
         ['has not reached the protocol end', '(incomplete)', 'F5 not evaluated: the audit does not cover'], ['Traceback']),
        ('sparse audit (gaps wider than one interval)', dict(dt_a=0.4), [], ['F1, F2, F3 unavailable', 'F5 not evaluated'],
         ['Traceback']),
        ('window margin 2 d with a 2.5 d support', dict(margin_d=2.0), [], ['NOT met', 'F2 verdict unavailable (window inside the mapping support)',
                                                                            'F5 particle counts seen'], ['Traceback']),
        ('lower window end at 2 d, field says 2.5 d', dict(zlo_d=2.0), [], ['F2 verdict unavailable (window inside the mapping support)'],
         ['Traceback']),
        ('old log: no win_margin / gradz_win_dyn, zero hydrostatic part', dict(margin_d=2.0, dyn=False), [],
         ['hydrostatic part zero', 'margin 2 d (binary before 627caf47)', 'F2 verdict unavailable (window inside the mapping support)'],
         ['Traceback']),
        ('old log with a 2 d support: margin met, total gradient admitted', dict(margin_d=2.0, dyn=False, width=2.0), [],
         ['hydrostatic part zero', 'F2 worst |ratio_k - 1|'], ['F2 verdict unavailable', 'Traceback']),
        ('old log with a hydrostatic part', dict(margin_d=2.0, dyn=False, width=2.0, dp_hydro=0.981), [],
         ['nonzero hydrostatic part', 'F2 verdict unavailable'], ['Traceback']),
    ]
    nfail = 0
    with tempfile.TemporaryDirectory(prefix='fl11_post_selftest_') as tmp:
        for i, (name, kw, extra, must, mustnot) in enumerate(cases):
            run = _write_run(os.path.join(tmp, 'run%02d' % i), proto, **kw)
            out = io.StringIO()
            try:
                with contextlib.redirect_stdout(out):
                    main([run] + [x.format(run=run) for x in extra])
            except Exception as exc:                     # noqa: BLE001  (a crash is a failed check)
                import traceback
                out.write('Traceback\n' + traceback.format_exc())
            text = out.getvalue()
            bad = [m for m in must if m not in text] + [m for m in mustnot if m in text]
            print(f"self-test {i:2d} {'ok  ' if not bad else 'FAIL'} {name}" + ('' if not bad else f": {bad}"))
            if bad:
                nfail += 1
                print('    ' + text.replace('\n', '\n    '))
    print(f"fl11_post self-test: {len(cases) - nfail}/{len(cases)} passed")
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
