#!/usr/bin/env python3
"""Case 11 (minimum fluidisation) evaluation of a q2p1_el_pipeflow run.

Reads run_slurm.log of RUNDIR:
  EL_BED_AUDIT       (SimPar@ELBedAuditFreq > 0; EL_BED_AUDIT module): t, step, plane pressures,
                     dp, its hydrostatic part, buoyant weight, floor contact, bed heights,
                     window porosity
  EL_FORCE_BUDGET    (every step): drag, pressure force, grav_buoy summed over the particles
  EL_INLET_TABLE     (echo of SimPar@ELInletVelocityTable): U(t), THE protocol of the run
  start-up echo      'EL kernel width factor' (mapping support = factor x d), 'EL fluid gravity',
                     'Global CFD domain extents' (z of the inlet plane and of the top: column
                     height); the deck _data/q2p1_param.dat when the echo lacks the width factor
  time lines         'time: t / t_end | ... | dt: dt' of the first step: start of the run and the
                     evaluation horizon (the deck TimeStep, MaxNumStep, StartSimTime otherwise)

The protocol. The official gates F1-F5 always come from the ECHOED table, each plateau exactly
once, with the interpolation rule of the solver (EL_INLET_TABLE_VALUE: linear between the rows,
the first and the last value held beyond the ends): the inlet history is the table extended by
the hold of its last value up to the evaluation horizon (and of its first value back to the start
of the run). A plateau is an interval of constant U longer than 0.1 s; a step starts where the
transition into it starts (the convention of make_inlet_table.py --steps); the first plateau
(preparation + step 0) is evaluated over its last 0.5 s. Without a horizon a terminal hold has
no known end: F5 is withheld when the lowest velocity is that hold.
A --steps file (make_inlet_table.py --steps) is validated against the echoed protocol in all
three columns and for multiplicity; it never feeds a gate. When it differs it gives a separately
labelled steps-file summary (not gated). Only a log WITHOUT a table echo is evaluated from the
steps file, and then every verdict carries '[protocol from the steps file, not from the echo]'.

Per step: mean and std of dp_dyn over the last 0.3 s and over the two sub-windows 0.2-0.35 s and
0.35-0.5 s of the step (flagged when they differ by more than --flag of the plateau), the
force-based dp (drag + pressure force)/A, bed height H = z99 + r - zmin, porosity of the bed
window. A step is complete when (1) the run reached the end of the step (an audit line at its end
time, within the print precision of the time) and (2) the audit samples the evaluation window at
its own cadence (cadence at most a tenth of the window, no gap wider than one cadence, at least
one distinct sample per cadence interval; tolerance: print precision of the time).

Gates of the spec (section 5, as amended on 2026-11-20):
  F1  plateau dp A / buoyant weight, mean over the steps with dp* > 0.95           (1 +- 0.02)
  F2  fixed branch: dynamic window gradient -<dp/dz> of the audit (bed core; gradz_win_dyn, or
      for older logs gradz_win when the hydrostatic part is zero: dp_hydro = 0 on every audit
      line used) against the closure at the code's own MAPPED gas fraction (eps_kern)  (+-5 %).
      Gated only when every contributing audit window excludes max(2 d, mapping support) at both
      ends (measured on every contributing line: lower end - inlet plane, z99 + r - upper end,
      and win_margin, 2 d for binaries before 627caf47; all three must be finite), every mapped
      sample is valid (finite, 0 < eps < 1: the old -1 sentinel and an unclipped 1 are refused),
      and the closure is positive and the ratio finite on EVERY fixed step. Otherwise the ratios
      are printed as a measurement and the verdict is 'unavailable' with the reason. The ratio at
      the geometric window porosity, dp_dyn / H_eq and dp_dyn / H are ungated.
  F3  U_mf = where the fit dp* = a U + b U^2 of the fixed branch (through the origin; steps
      below 0.8 of the plateau) meets dp* = 1, against the closure U_mf at the WHOLE-BED porosity
      1 - N V_p / (A (z99 + r - zmin))                                               (+-5 %);
      the crossing with the measured plateau and the ratios at the core porosity are reported
  F4  U_mf against the measurement 0.182 m/s (reported)
  F5  (a) particle count constant over the run: the audit covers the run from its start to the
      horizon at its cadence; (b) no centre above the physical mid-column zmin + Lz/2 (both from
      the echoed extents) during the plateau of the lowest velocity of the protocol (the constant
      part, audit coverage required; the 20 ms transition into it is reported separately and
      labelled). Evaluated independently of F1-F3.
Every ratio needs a positive finite denominator (area of the inlet plane, buoyant weight,
closure value); otherwise the quantity is printed as nan and the gates that need it are
'unavailable' instead of an exception.
The closures come from fl11_closures.py (shared with core_deficit.py); --closure gidaspow
(default), ergun, wenyu: dp/dz(U, eps) = n F_d(U/eps) / eps with the per-particle forms of
EL_FORCES (model A; identical to the spec's umf_reference.py for ergun and wenyu).
fl11_movie.py takes its step table from prepare() and step_table() of this module.

usage: fl11_post.py RUNDIR [--steps FILE] [--d 0.87e-3] [--rhop 860] [--rhof 1.0] [--mu 1.8e-5]
                    [--g 9.81] [--lz LZ] [--closure gidaspow] [--flag 0.02] [--series OUT.csv]
       fl11_post.py --self-test       synthetic logs (the counterexamples of reviews r2 and r3)
--lz is only checked against the column height of the echoed extents (which F5 uses)."""
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
TOL_Z = 2e-8          # height comparisons (m): print precision of the ES14.6 heights of the audit line
FROM_STEPS = ' [protocol from the steps file, not from the echo]'


def f(x):
    return NAN if 'NaN' in x else float(x.replace('D', 'E').replace('d', 'e'))


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def pos(x):
    """x when it is a positive finite number (a usable denominator), else nan"""
    return x if finite(x) and x > 0 else NAN


def div(a, b):
    """a / b with b a positive finite number; nan otherwise (never an exception)"""
    b = pos(b)
    return a / b if finite(b) and finite(a) else NAN


def tol_t(t):
    """print precision of a time of the audit line (ES14.6: 7 significant digits)"""
    return 1e-6 * max(1.0, abs(t)) if finite(t) else 1e-6


# ---------------------------------------------------------------------------------------------------
# reading
# ---------------------------------------------------------------------------------------------------
# time line of a step: 'time: t / t_max | itns: istep / nsteps | dt: dt' (D format with 4 digits)
PAT_TIME = re.compile(r'^\s*time:\s*(%s)\s*/\s*(%s)\s*\|\s*itns:\s*(\d+)\s*/\s*(\d+)\s*\|\s*dt:\s*(%s)' % (NUM, NUM, NUM))
TOL_D4 = 5e-4          # relative print precision of the 4-digit D format of the time line


def parse(rundir):
    audit, budget, table, pairs = [], {}, [], {}
    dt_log = NAN
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
                              fbal=f(g[20]) if g[20] else NAN,
                              has_fbal=bool(g[20]),
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
        elif not finite(dt_log) and line.lstrip().startswith('time:'):
            m = PAT_TIME.match(line)
            if m:
                dt_log = f(m.group(5))
    # floor force from the particle momentum balance: EL_BED_AUDIT floor_bal_z; an older log without
    # the field: from the EL_NEWTON_PAIR line of the same step (the same two numbers), only with that
    # line and a positive finite time step (the log's dt, else the audit spacing); nan otherwise
    dt = dt_log
    if not finite(pos(dt)) and len(audit) > 1 and audit[1]['step'] != audit[0]['step']:
        dt = (audit[1]['t'] - audit[0]['t']) / (audit[1]['step'] - audit[0]['step'])
    for row in audit:
        if not row['has_fbal']:
            pr = pairs.get(row['step'])
            row['fbal'] = (pr[0] - pr[1]) / dt if pr and finite(pos(dt)) else NAN
    return audit, budget, table


def parse_setup(rundir):
    """kernel width factor (and where it came from), fluid gravity, extents of the column (z), start
    of the run, evaluation horizon (and where it came from) and the time the run reached.
    Time: the solver stops at min(MaxSimTime, start + MaxNumStep dt). The time lines of the first and
    the last step of the log give start, t_max, step counts and dt with 4 digits; the deck gives
    StartSimTime, MaxSimTime, TimeStep, MaxNumStep exactly. The deck values are used when they agree
    with the log within its print precision (or the log has no time line), the log values otherwise
    (they are what ran). t_run = start + (steps of the last time line) x dt."""
    setup = dict(width=None, width_src=None, fluid_gravity=None, zmin=None, zmax=None,
                 t_start=NAN, horizon=NAN, horizon_src=None, t_run=NAN)
    pat_w = re.compile(r'^EL kernel width factor\s*=\s*(%s)' % NUM)
    pat_g = re.compile(r'^EL fluid gravity\s*=\s*([TF])')
    pat_x = re.compile(r'^Global CFD domain extents:\s*(%s)\s+(%s)\s+(%s)\s+(%s)\s+(%s)\s+(%s)' % ((NUM,) * 6))
    first = last = None
    for line in open(os.path.join(rundir, 'run_slurm.log'), errors='replace'):
        if line.lstrip().startswith('time:'):
            m = PAT_TIME.match(line)
            if m:
                first = (f(m.group(1)), f(m.group(2)), int(m.group(3)), int(m.group(4)), f(m.group(5)))
                break                                     # the first time line follows the start-up echo
        m = pat_w.match(line)
        if m and setup['width'] is None:
            setup['width'], setup['width_src'] = f(m.group(1)), 'log echo'
        m = pat_g.match(line)
        if m and setup['fluid_gravity'] is None:
            setup['fluid_gravity'] = (m.group(1) == 'T')
        m = pat_x.match(line)
        if m and setup['zmin'] is None:
            setup['zmin'], setup['zmax'] = f(m.group(5)), f(m.group(6))
    with open(os.path.join(rundir, 'run_slurm.log'), 'rb') as fh:        # the last time line: the tail
        fh.seek(0, 2)
        fh.seek(max(0, fh.tell() - 4000000))
        tail = fh.read().decode(errors='replace').splitlines()
    for line in reversed(tail):
        if line.lstrip().startswith('time:'):
            m = PAT_TIME.match(line)
            if m:
                last = (f(m.group(1)), f(m.group(2)), int(m.group(3)), int(m.group(4)), f(m.group(5)))
                break
    deck = os.path.join(rundir, '_data', 'q2p1_param.dat')
    keys = {}
    if os.path.exists(deck):
        for line in open(deck, errors='replace'):
            m = re.match(r'^\s*SimPar@(\w+)\s*=\s*(%s)\s*$' % NUM, line.strip())
            if m:
                keys[m.group(1)] = f(m.group(2))
    if setup['width'] is None:
        if 'ELKernelWidthFactor' in keys:
            setup['width'], setup['width_src'] = keys['ELKernelWidthFactor'], 'deck'
        else:
            setup['width'], setup['width_src'] = 2.5, 'EL_CONFIG default'

    def agree(a, b):
        return finite(a) and finite(b) and abs(a - b) <= TOL_D4 * max(abs(a), abs(b), 1.0 if b == 0 else 0.0)
    log_ok = first is not None and all(finite(x) for x in (first[0], first[1], first[4])) and first[4] > 0
    deck_ok = all(k in keys for k in ('TimeStep', 'MaxNumStep')) and finite(pos(keys['TimeStep']))
    dt = keys['TimeStep'] if deck_ok and (not log_ok or agree(keys['TimeStep'], first[4])) else (first[4] if log_ok else NAN)
    t_log = first[0] - first[2] * dt if log_ok else NAN
    h_log = min(first[1], t_log + first[3] * dt) if log_ok else NAN
    if deck_ok:
        t_deck = keys.get('StartSimTime', 0.0)
        h_deck = min(keys.get('MaxSimTime', math.inf), t_deck + keys['TimeStep'] * keys['MaxNumStep'])
    if deck_ok and (not log_ok or (agree(t_deck, t_log) or abs(t_deck - t_log) <= TOL_D4) and agree(h_deck, h_log)):
        setup['t_start'], setup['horizon'] = t_deck, h_deck
        setup['horizon_src'] = 'deck: min(MaxSimTime, StartSimTime + TimeStep x MaxNumStep)' + (', the log agrees' if log_ok else '')
    elif log_ok:
        setup['t_start'], setup['horizon'] = t_log, h_log
        setup['horizon_src'] = 'log: min(t_max, start + nsteps x dt) of the first time line' + (' (the deck disagrees)' if deck_ok else '')
    if last is not None and finite(setup['t_start']) and finite(dt):
        setup['t_run'] = setup['t_start'] + last[2] * dt
    for k in ('zmin', 'zmax'):
        if setup[k] is not None and not finite(setup[k]):
            setup[k] = None
    return setup


# ---------------------------------------------------------------------------------------------------
# protocol
# ---------------------------------------------------------------------------------------------------
def u_of_t(table, t):
    if t <= table[0][0]:
        return table[0][1]
    for (t0, u0), (t1, u1) in zip(table, table[1:]):
        if t <= t1:
            return u0 + (u1 - u0) * (t - t0) / (t1 - t0)
    return table[-1][1]


def inlet_history(table, t_start=NAN, horizon=NAN):
    """the rows of U(t) as the solver applies them between the start of the run and the horizon:
    the table, its first value held back to t_start, its last value held up to the horizon (cut at
    the horizon when the table is longer). Second value: True when the last value is held beyond the
    table with no known end (no horizon)."""
    if not table:
        return [], False
    knots = list(table)
    if finite(t_start) and t_start < knots[0][0] - tol_t(knots[0][0]):
        knots.insert(0, (t_start, knots[0][1]))
    if not finite(horizon):
        return knots, True
    if horizon > knots[-1][0] + tol_t(horizon):
        knots.append((horizon, knots[-1][1]))
    elif horizon < knots[-1][0] - tol_t(horizon):
        u_h = u_of_t(knots, horizon)
        knots = [k for k in knots if k[0] < horizon] + [(horizon, u_h)]
    return knots, False


def plateaus(table):
    """steps from the rows of U(t): intervals of constant U longer than 0.1 s. A step starts where the
    transition into it starts (the convention of the steps file of make_inlet_table.py), so the
    sub-windows 0.2-0.35 s and 0.35-0.5 s do not depend on the source of the step list. The
    preparation is merged with the first plateau; protocol_steps() evaluates the last 0.5 s of it."""
    out = []
    for i, ((t0, u0), (t1, u1)) in enumerate(zip(table, table[1:])):
        if u0 == u1 and t1 - t0 > 0.1:
            ramp = i > 0 and table[i - 1][1] != u0
            out.append((table[i - 1][0] if ramp else t0, t1, u0))
    return out


def protocol_steps(knots):
    steps = plateaus(knots)
    if steps and steps[0][1] - steps[0][0] > 0.6:
        # the first plateau of the protocol holds preparation + step 0: evaluate its last 0.5 s
        steps[0] = (steps[0][1] - 0.5, steps[0][1], steps[0][2])
    return steps


def lowest_plateau(knots):
    """(u, [(t_a, t_b)], [(t_transition_start, t_a)]): the lowest velocity held constant, its constant
    intervals and the transitions into them; None without a plateau"""
    held = []
    for i, ((t0, u0), (t1, u1)) in enumerate(zip(knots, knots[1:])):
        if u0 == u1 and t1 - t0 > 0.1:
            held.append((u0, t0, t1, knots[i - 1][0] if i > 0 and knots[i - 1][1] != u0 else t0))
    if not held:
        return None
    umin = min(h[0] for h in held)
    return (umin, [(h[1], h[2]) for h in held if h[0] == umin], [(h[3], h[1]) for h in held if h[0] == umin])


def read_steps(path):
    """rows 't_start t_end U' of a steps file; (rows, errors)"""
    rows, errors = [], []
    if not os.path.isfile(path):
        return rows, ['steps file %s not found' % path]
    for k, line in enumerate(open(path), 1):
        s = line.strip()
        if not s or s.startswith('#'):
            continue
        c = s.split()
        try:
            r = tuple(float(x) for x in c[:3])
        except ValueError:
            r = ()
        if len(r) != 3 or not all(finite(x) for x in r) or not r[1] > r[0]:
            errors.append('line %d %r: not three finite numbers with t_start < t_end' % (k, s))
            continue
        rows.append(r)
    return rows, errors


def compare_steps(listed, protocol):
    """a steps file against the protocol, all three columns (time within the print precision %.6g
    of the steps file, U relative 1e-6) and multiplicity: (exact one-to-one match, findings)"""
    def same(a, b):
        return (abs(a[0] - b[0]) <= 1e-5 * max(1.0, abs(b[0])) and abs(a[1] - b[1]) <= 1e-5 * max(1.0, abs(b[1]))
                and abs(a[2] - b[2]) <= 1e-6 * max(1.0, abs(b[2])))
    hits = [[i for i, p in enumerate(protocol) if same(s, p)] for s in listed]
    count = [sum(1 for h in hits if i in h) for i in range(len(protocol))]
    findings = []
    for s, h in zip(listed, hits):
        if not h:
            near = [p for p in protocol if abs(s[1] - p[1]) <= 1e-5 * max(1.0, abs(p[1])) and abs(s[2] - p[2]) <= 1e-6 * max(1.0, abs(p[2]))]
            findings.append('listed %s %s' % (_fs(s), 'differs from the plateau %s in the start time' % _fs(near[0])
                                                if near else 'is no plateau of the protocol'))
    for p, c in zip(protocol, count):
        if c == 0:
            findings.append('plateau %s missing from the list' % _fs(p))
        elif c > 1:
            findings.append('plateau %s listed %d times' % (_fs(p), c))
    return not findings, findings


def _fs(s):
    return '(%.4f, %.4f, %.4f)' % s


# ---------------------------------------------------------------------------------------------------
# sampling
# ---------------------------------------------------------------------------------------------------
def cadence(rows):
    """the audit cadence: median spacing of the distinct audit times; nan below two"""
    ts = sorted({r['t'] for r in rows if finite(r['t'])})
    d = sorted(b - a for a, b in zip(ts, ts[1:]))
    return d[len(d) // 2] if d else NAN


audit_interval = cadence              # name of round 2


def reached(rows, t_end, t_run=NAN):
    """the run reached t_end (print precision of the time): an audit line at t_end or later, or the time
    line of the last step of the log (t_run) at t_end or later. Whether the audit SAMPLED the interval
    up to there is a separate question (sampled())."""
    if not finite(t_end):
        return False
    tol = tol_t(t_end)
    return (bool(rows) and max(r['t'] for r in rows) >= t_end - tol) or (finite(t_run) and t_run >= t_end - tol)


def sampled(rows, t0, t1, cad):
    """(ok, reason): the audit samples [t0, t1] at its cadence: cadence at most a tenth of the window,
    the first sample within one cadence of t0 and the last within one cadence of t1, no gap wider than
    one cadence, n - 1 distinct samples for n cadence intervals (tolerance: print precision of the
    time)"""
    if not (finite(cad) and cad > 0):
        return False, 'no audit cadence'
    tol = 2 * tol_t(max(abs(t0), abs(t1)))
    if t1 <= t0 or cad > (t1 - t0) / 10 + tol:
        return False, 'audit cadence %.6g s wider than a tenth of the window' % cad
    ts = sorted({r['t'] for r in rows if t0 - tol <= r['t'] <= t1 + tol})
    if not ts:
        return False, 'no audit line in the window'
    if ts[0] - t0 > cad + tol or t1 - ts[-1] > cad + tol:
        return False, 'window %.4f-%.4f s sampled %.4f-%.4f s only' % (t0, t1, ts[0], ts[-1])
    gap = max((b - a for a, b in zip(ts, ts[1:])), default=0.0)
    if gap > cad + tol:
        return False, 'gap of %.6g s in the window (cadence %.6g s)' % (gap, cad)
    # n cadence intervals in the window hold at least n - 1 distinct samples (both ends may fall between
    # two audit lines); repeated time stamps count once
    need = max(1, int(math.floor((t1 - t0) / cad + 1e-6)) - 1)
    if len(ts) < need:
        return False, '%d distinct samples, %d needed' % (len(ts), need)
    return True, ''


def window(rows, t0, t1):
    tol = tol_t(max(abs(t0), abs(t1)))
    return [r for r in rows if t0 - tol <= r['t'] <= t1 + tol]


def stats(rows, key, t0, t1, open_end=False):
    """mean, std, count of rows[key] over t0 <= t <= t1 (t < t1 with open_end: disjoint windows)"""
    eps = 1e-9
    v = [r[key] for r in rows if t0 - eps <= r['t'] and (r['t'] < t1 - eps if open_end else r['t'] <= t1 + eps)]
    if not v:
        return NAN, NAN, 0
    m = sum(v) / len(v)
    s = math.sqrt(sum((x - m) ** 2 for x in v) / len(v))
    return m, s, len(v)


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


# ---------------------------------------------------------------------------------------------------
# the evaluation shared with fl11_movie.py
# ---------------------------------------------------------------------------------------------------
def prepare(rundir, d=0.87e-3, steps_path=None):
    """everything the step table needs: audit rows with U, H (= z99 + r - zmin), dp_force; area,
    buoyant weight W, W/A (nan when not positive finite), cadence, setup, the inlet history, the
    official protocol and how it was obtained ('echo' or 'steps file')"""
    audit, budget, table = parse(rundir)
    setup = parse_setup(rundir)
    ev = dict(audit=audit, budget=budget, table=table, setup=setup, d=d, r=d / 2, steps_path=steps_path,
              listed=None, listed_errors=[], listed_ok=None, listed_findings=[])
    if steps_path:
        ev['listed'], ev['listed_errors'] = read_steps(steps_path)
    knots, open_hold = inlet_history(table, setup['t_start'], setup['horizon'])
    ev['knots'], ev['open_hold'] = knots, open_hold
    if table:
        ev['source'] = 'echo'
        ev['protocol'] = protocol_steps(knots)
        if ev['listed'] is not None:
            ev['listed_ok'], ev['listed_findings'] = compare_steps(ev['listed'], ev['protocol'])
            if ev['listed_errors']:
                ev['listed_ok'] = False
    elif ev['listed'] is not None and not ev['listed_errors']:
        ev['source'] = 'steps file'
        ev['protocol'] = sorted(ev['listed'])
    else:
        ev['source'] = None
        ev['protocol'] = []
    zmin = setup['zmin'] if setup['zmin'] is not None else NAN
    area = pos(audit[-1]['a_in']) if audit else NAN    # nan without an inlet plane (no element face at zmin)
    for row in audit:
        row['U'] = u_of_t(knots, row['t']) if knots else NAN
        fb = budget.get(row['step'])
        row['dp_force'] = div(fb[0] + fb[1], area) if fb else NAN
        row['dp_drag'] = div(fb[0], area) if fb else NAN
        row['H'] = row['z99'] + ev['r'] - zmin                  # height of the bed above the inlet plane
    ev['area'] = area
    ev['W'] = pos(-audit[-1]['wz']) if audit else NAN       # nan for an empty bed
    ev['wa'] = div(ev['W'], area)
    ev['cad'] = cadence(audit)
    return ev


def step_table(ev, steps, flag=0.02):
    """per step: dp_dyn statistics (dp, std, sub-windows), dp_force, H, eps_bed, floor balance and the
    two completeness criteria (reached, sampled); guarded ratios dps (dp*), stds, m1s, m2s, mfs, fls"""
    audit, wa, W, cad = ev['audit'], ev['wa'], ev['W'], ev['cad']
    res = []
    for k, (t0, t1, u) in enumerate(steps):
        lo = min(t1 - 0.3, t0 + 0.2)
        got = reached(audit, t1, ev['setup']['t_run'])
        ok, why = sampled(audit, lo, t1, cad) if got else (False, 'the run did not reach %.4f s' % t1)
        m, s, cnt = stats(audit, 'dp_dyn', t1 - 0.3, t1)
        m1, _, _ = stats(audit, 'dp_dyn', t0 + 0.2, t0 + 0.35, open_end=True)
        m2, _, _ = stats(audit, 'dp_dyn', t0 + 0.35, t0 + 0.5)
        mf, _, _ = stats(audit, 'dp_force', t1 - 0.3, t1)
        h, _, _ = stats(audit, 'H', t1 - 0.3, t1)
        e, _, _ = stats(audit, 'eps', t1 - 0.3, t1)
        fl, _, _ = stats(audit, 'fbal', t1 - 0.3, t1)
        flpe, _, _ = stats(audit, 'floor', t1 - 0.3, t1)
        res.append(dict(k=k, t0=t0, t1=t1, U=u, dp=m, std=s, cnt=cnt, eps=e, H=h, fl=div(fl, W), flpe=div(flpe, W),
                        dps=div(m, wa), stds=div(s, wa), m1s=div(m1, wa), m2s=div(m2, wa), mfs=div(mf, wa),
                        flag=finite(wa) and abs(m1 - m2) > flag * wa, reached=got,
                        complete=ok and cnt > 0, why=why))
    return res


def print_step_table(res, title):
    print(f'\n{title}')
    print(' step  U [m/s]  window      dp*      std     sub1    sub2   flag  F/(W/A)  H [mm]  eps_bed  floor/W  dp*+floor  (floor_pe/W)')
    for x in res:
        print(f" {x['k']:3d}  {x['U']:6.3f}  {x['t1']-0.3:5.2f}-{x['t1']:5.2f}  {x['dps']:7.4f}  {x['stds']:7.4f}  {x['m1s']:6.3f}  {x['m2s']:6.3f}  "
              f"{'DIFF' if x['flag'] else '    '}  {x['mfs']:7.4f}  {x['H']*1e3:6.2f}  {x['eps']:7.4f}  {x['fl']:7.4f}  {x['dps'] + x['fl']:8.4f}   ({x['flpe']:.4f})"
              f"{'' if x['complete'] else '  (incomplete: ' + x['why'] + ')'}")
    print('  dp* = dp_dyn A / W over 0.2-0.5 s of the step; floor from the particle momentum balance; floor_pe = the pe '
          'wall-contact impulse (unrelaxed PGS impulse, 1/0.9 of the applied one)')


# ---------------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir', nargs='?')
    ap.add_argument('--steps', default=None)
    ap.add_argument('--d', type=float, default=0.87e-3)
    ap.add_argument('--rhop', type=float, default=860.0)
    ap.add_argument('--rhof', type=float, default=1.0)
    ap.add_argument('--mu', type=float, default=1.8e-5)
    ap.add_argument('--g', type=float, default=9.81)
    ap.add_argument('--lz', type=float, default=None)
    ap.add_argument('--closure', default='gidaspow', choices=CLOSURES)
    ap.add_argument('--flag', type=float, default=0.02)
    ap.add_argument('--series', default=None)
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if not a.rundir:
        ap.error('RUNDIR is required')

    ev = prepare(a.rundir, a.d, a.steps)
    audit, setup = ev['audit'], ev['setup']
    if not audit:
        sys.exit('no EL_BED_AUDIT lines (SimPar@ELBedAuditFreq not set?)')
    if ev['source'] is None:
        sys.exit('no EL_INLET_TABLE echo (SimPar@ELInletVelocityTable not set?) and no valid --steps file'
                 + (': ' + '; '.join(ev['listed_errors']) if ev['listed_errors'] else ''))
    label = '' if ev['source'] == 'echo' else FROM_STEPS
    area, W, wa = ev['area'], ev['W'], ev['wa']
    vp = math.pi * a.d ** 3 / 6
    n0 = audit[0]['n']
    tend = audit[-1]['t']
    print(f"{os.path.basename(os.path.abspath(a.rundir))}: {len(audit)} audit lines, t {audit[0]['t']:.4f} .. {tend:.4f} s, "
          f"N {n0}, A {area:.6e} m2, buoyant weight {W:.6e} N, W/A {wa:.3f} Pa, "
          f"hydrostatic part in dp {audit[-1]['dp_h']:.4f} Pa (subtracted: dp_dyn)")
    if not finite(wa):
        print('  W/A unavailable (' + ', '.join(x for x, ok in (('no inlet plane: area_in not positive', finite(area)),
                                                                ('empty bed: no buoyant weight', finite(W))) if not ok)
              + '): dp*, F1, F2, F3 unavailable')
    if setup['zmin'] is None:
        print('  column extents not in the log (Global CFD domain extents): bed height, F3 and F5 unavailable')
    if a.series:
        with open(a.series, 'w') as fh:
            fh.write('t,step,U,dp,dp_hydro,dp_dyn,dp_force,dp_drag,floor_z,z50,z95,z99,zmax,eps_bed,n\n')
            for row in audit:
                fh.write(','.join(str(row[k]) for k in ('t', 'step', 'U', 'dp', 'dp_h', 'dp_dyn', 'dp_force',
                                                        'dp_drag', 'floor', 'z50', 'z95', 'z99', 'zmax', 'eps', 'n')) + '\n')

    protocol = ev['protocol']
    if ev['source'] == 'echo':
        print(f"\nprotocol: the {len(protocol)} plateaus of the echoed table (each once); horizon "
              + (f"{setup['horizon']:.4f} s [{setup['horizon_src']}]" if finite(setup['horizon']) else
                 'unknown (no time line, no deck): the hold of the last table value has no known end'))
    else:
        print(f"\nprotocol: {len(protocol)} steps of the steps file {a.steps}{FROM_STEPS}")
    proto_end = max((p[1] for p in protocol), default=NAN)
    print(f"protocol end {proto_end:.4f} s, audit to {tend:.4f} s, audit cadence {ev['cad']:.6f} s"
          + ('' if reached(audit, proto_end, setup['t_run']) else ': the run has not reached the protocol end'))
    res = step_table(ev, protocol, a.flag)
    print_step_table(res, 'steps of the protocol' + label)

    if ev['source'] == 'echo' and ev['listed'] is not None:
        if ev['listed_ok']:
            print(f"steps file {a.steps}: identical to the echoed protocol (all three columns, each plateau once)")
        else:
            print(f"steps file {a.steps}: NOT the echoed protocol: " + '; '.join(ev['listed_errors'] + ev['listed_findings']))
            print_step_table(step_table(ev, ev['listed'], a.flag),
                             'steps-file summary, NOT gated (the gates below use the echoed protocol)')

    gates_f1_f3(a, res, ev, n0, vp, label)
    gate_f5(a, ev, label)


def gates_f1_f3(a, res, ev, n0, vp, label=''):
    """F1, F3, F4 and F2 (each 'unavailable' with the reason when its inputs are not there)"""
    audit, setup, wa, area, r = ev['audit'], ev['setup'], ev['wa'], ev['area'], ev['r']
    use = [x for x in res if x['complete'] and finite(x['dp'])]
    if not finite(wa):
        print('\nF1, F2, F3 unavailable: W/A has no positive finite value' + label)
        return
    if len(use) < 3:
        print(f'\nF1, F2, F3 unavailable: {len(use)} complete step(s) with a measured dp, three needed' + label)
        return
    # ---- F1 (spec rule of 2026-11-20): mean dp* over all steps with dp* > 0.95 ------------------
    flu = [x for x in use if x['dp'] > 0.95 * wa]
    if flu:
        plateau = sum(x['dp'] for x in flu) / len(flu)
        print(f"\nF1 plateau: mean dp* over {len(flu)} steps with dp* > 0.95 (U {min(x['U'] for x in flu):.2f} .. "
              f"{max(x['U'] for x in flu):.2f}): {plateau / wa:.6f}   (gate 1 +- 0.02: {'PASS' if abs(plateau / wa - 1) <= 0.02 else 'FAIL'}){label}; "
              f"dp* + floor/W over these {sum(x['dp']/wa + x['fl'] for x in flu)/len(flu):.4f}")
    else:
        plateau = wa
        print('\nF1 not evaluated: no step with dp* > 0.95; plateau taken as W/A' + label)
    fixed = [x for x in use if x['dp'] < 0.8 * plateau]
    if len(fixed) < 2:
        print(f'F2, F3 unavailable: {len(fixed)} fixed-branch step(s) (dp* < 0.8 of the plateau), two needed' + label)
        return
    efix = sum(x['eps'] for x in fixed) / len(fixed)
    uc = closure_umf(efix, a.d, a.rhop, a.rhof, a.mu, a.g, a.closure)
    # whole-bed porosity 1 - N V_p / (A H), H = z99 + r - zmin: the reference porosity of F3 (spec, 2026-11-20)
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
          f"ratio {r3:.4f}   (gate 0.05: {('PASS' if abs(r3-1) <= 0.05 else 'FAIL') if finite(r3) else 'unavailable'}){label}")
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
    gate_f2(a, fixed, ev, n0, vp, label)


def valid_eps(e):
    """a mapped porosity the closure can be evaluated at: finite and 0 < eps < 1 (the -1 sentinel of
    old binaries and an unclipped 1 are not)"""
    return finite(e) and 0 < e < 1


def gate_f2(a, fixed, ev, n0, vp, label=''):
    """F2 in readings: (k) the dynamic window gradient of the audit against the closure at the mapped
    gas fraction (the gate), (w) the same gradient at the geometric window porosity, (e) dp_dyn / H_eq
    with H_eq the height of a uniform bed at eps_bed, (h) dp_dyn / H with H = z99 + r - zmin against the
    closure at the mean porosity of that height (end layers included)"""
    audit, setup, area, r = ev['audit'], ev['setup'], ev['area'], ev['r']
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
    # -- the window margins: a finite lower, upper and reported-or-legacy margin on every contributing line
    required = max(2 * a.d, setup['width'] * a.d)
    tol = max(1e-6 * required, TOL_Z)
    zmin = setup['zmin']
    lows = [row['zlo'] - zmin if zmin is not None else NAN for row in rows_all]
    highs = [row['z99'] + r - row['zhi'] for row in rows_all]
    fields = [row['margin'] if row['margin'] is not None else 2 * a.d for row in rows_all]   # before 627caf47: 2 d
    if zmin is None:
        margin_ok = False
        reasons.append('missing geometry: no inlet-plane elevation in the log (Global CFD domain extents), '
                       'the lower window margin cannot be measured')
    elif not rows_all or not all(finite(x) for x in lows + highs + fields):
        margin_ok = False
        reasons.append('window bounds not finite on every contributing audit line')
    else:
        margin_ok = min(lows + highs + fields) >= required - tol
        if not margin_ok:
            reasons.append('window inside the mapping support')
    fin = [x for x in lows + highs + fields if finite(x)]
    print(f" F2: gradient {gtext}")
    print(f"     audit window margin over {len(rows_all)} lines: min {div(min(fin), a.d) if fin else NAN:.3f} d "
          f"(lower end - inlet plane min {_fmin(lows)/a.d:.3f} d, z99 + r - upper end min {_fmin(highs)/a.d:.3f} d, "
          f"{'win_margin field' if all(row['margin'] is not None for row in rows_all) else 'margin 2 d (binary before 627caf47)'} "
          f"min {_fmin(fields)/a.d:.3f} d); required max(2 d, mapping support {setup['width']:g} d "
          f"[{setup['width_src']}]) = {required/a.d:.3f} d: {'met' if margin_ok else 'NOT met'}")
    # -- the mapped samples, each validated before the mean
    bad_e = [row for row in rows_all if not valid_eps(row['ekern'])]
    bad_g = [row for row in rows_all if not finite(row[gkey])]
    if bad_e:
        reasons.append(f'{len(bad_e)} mapped gas-fraction sample(s) of the contributing audit lines not usable '
                       '(non-finite, outside 0 < eps < 1, or the -1 sentinel)')
    if bad_g:
        reasons.append(f'{len(bad_g)} window-gradient sample(s) not finite')
    print(' F2: step  U      eps_bed  window-grad  closure(eps_bed)  ratio_w | eps_kern ratio_k | dp/H_eq  ratio_e | eps_H   dp/H   closure(eps_H)  ratio_h')
    for x in fixed:
        x['g'], _, _ = stats(audit, gkey, x['t1'] - 0.3, x['t1'])
        x['ek'], _, _ = stats(audit, 'ekern', x['t1'] - 0.3, x['t1'])
    rks, rws, undefined = [], [], []
    for x in fixed:
        heq = div(n0 * vp, area * (1 - x['eps']))
        clo = closure_grad(x['U'], x['eps'], a.d, a.rhof, a.mu, a.closure)
        eh = 1 - div(n0 * vp, area * x['H'])
        cloh = closure_grad(x['U'], eh, a.d, a.rhof, a.mu, a.closure)
        clok = closure_grad(x['U'], x['ek'], a.d, a.rhof, a.mu, a.closure) if valid_eps(x['ek']) else NAN
        rw = div(x['g'], clo)
        rk = div(x['g'], clok)
        re_ = div(div(x['dp'], heq), clo)
        rh = div(div(x['dp'], x['H']), cloh)
        if not (finite(pos(clok)) and finite(rk)):
            undefined.append(f"step {x['k']} (U {x['U']:.3f}): closure {clok:.4g}, ratio {rk:.4g}")
        rks.append(rk)
        rws.append(rw)
        print(f"     {x['k']:3d}  {x['U']:.3f}  {x['eps']:.4f}  {x['g']:10.1f}  {clo:10.1f}  {rw:12.4f} | {x['ek']:.4f} {rk:7.4f} | {div(x['dp'], heq):8.1f} {re_:7.4f} | "
              f"{eh:.4f} {div(x['dp'], x['H']):8.1f} {cloh:10.1f}  {rh:9.4f}")
    if undefined:
        reasons.append('no positive closure and finite ratio_k on every fixed step: ' + '; '.join(undefined))
    worst_k = max(abs(x - 1) for x in rks) if all(finite(x) for x in rks) else NAN
    worst_w = max(abs(x - 1) for x in rws) if all(finite(x) for x in rws) else NAN
    if not reasons:
        print(f"F2 worst |ratio_k - 1| (window gradient against the closure at the mapped gas fraction) {worst_k:.4f}   "
              f"(gate 0.05: {'PASS' if worst_k <= 0.05 else 'FAIL'}){label}")
    else:
        print(f"F2 measured, not gated: worst |ratio_k - 1| {worst_k:.4f}")
        print(f"F2 verdict unavailable ({'; '.join(reasons)}){label}")
    if finite(worst_w):
        print(f"   reported, not gated: worst |ratio_w - 1| at the geometric window porosity (rule before 2026-11-20) {worst_w:.4f}")


def _fmin(v):
    v = [x for x in v if finite(x)]
    return min(v) if v else NAN


def gate_f5(a, ev, label=''):
    """F5: (a) the particle count over the whole run, (b) the containment on the lowest plateau of the
    protocol; independent of the other gates"""
    audit, setup, cad = ev['audit'], ev['setup'], ev['cad']
    nset = sorted({row['n'] for row in audit})
    # (a) count coverage over the run: from its start to the horizon at the audit cadence
    t_a, t_b = setup['t_start'], setup['horizon']
    if not (finite(t_a) and finite(t_b)):
        count_ok, count_txt = False, 'start / end of the run unknown (no time line, no deck)'
    else:
        ok, why = sampled(audit, t_a, t_b, cad) if reached(audit, t_b, setup['t_run']) else (False, f'the run did not reach {t_b:.4f} s')
        count_ok = ok and len(nset) == 1 and nset[0] > 0
        count_txt = (f"counts seen {nset} on {len(audit)} audit lines, run {t_a:.4f}-{t_b:.4f} s "
                     + ('covered' if ok else f'NOT covered ({why})'))
    print(f"F5 (a) particle count over the run: {count_txt}")
    # (b) the lowest plateau
    if ev['source'] == 'echo':
        low = lowest_plateau(ev['knots'])
        u_last = ev['knots'][-1][1] if ev['knots'] else NAN
        if ev['open_hold'] and (low is None or u_last <= low[0]):
            # the last value is held from the last row on, for as long as the run lasted: a plateau of
            # unknown end that is (or may be) the lowest one
            print(f"F5 withheld: the lowest velocity {u_last:.3f} m/s is the hold of the last table value, whose end "
                  f"is unknown (no evaluation horizon in the log or the deck){label}")
            return
        trans_txt = 'transition'
    else:
        lo = min(ev['protocol'], key=lambda s: s[2]) if ev['protocol'] else None
        low = (lo[2], [(lo[0], lo[1])], []) if lo else None
        trans_txt = None
    if low is None:
        print(f"F5 not evaluated: the protocol holds no velocity plateau{label}")
        return
    u, held, trans = low
    span = ', '.join(f'{ta:.3f}-{tb:.3f} s' for ta, tb in held)
    cover = [(sampled(audit, ta, tb, cad) if reached(audit, tb, setup['t_run']) else (False, f'the run did not reach {tb:.4f} s'))
             for ta, tb in held]
    if not all(c[0] for c in cover):
        have = window(audit, held[0][0], held[-1][1])
        print(f"F5 not evaluated: the audit does not cover the plateau of the lowest protocol velocity {u:.3f} m/s "
              f"({span}; {len(have)} audit lines in it: {'; '.join(c[1] for c in cover if not c[0])}){label}")
        return
    if setup['zmin'] is None or setup['zmax'] is None:
        print(f"F5 unavailable: no column extents in the log (Global CFD domain extents): mid-column unknown{label}")
        return
    zmin, lz = setup['zmin'], setup['zmax'] - setup['zmin']
    if a.lz is not None and abs(a.lz - lz) > 1e-9 * max(1.0, lz):
        print(f"   note: --lz {a.lz:g} differs from the column height {lz:g} m of the echoed extents (the extents are used)")
    mid = zmin + lz / 2
    rows = [row for ta, tb in held for row in window(audit, ta, tb)]
    zs = [row['zmax'] for row in rows]
    if not all(finite(z) for z in zs) or nset == [0]:
        print(f"F5 unavailable: no particle centre height on every audit line of the lowest plateau (particle counts {nset}){label}")
        return
    zm = max(zs)
    ok = count_ok and zm < mid
    verdict = 'PASS' if ok else ('FAIL' if zm >= mid or (len(nset) > 1) else 'unavailable (count over the run not established)')
    print(f"F5 (b) highest centre during the plateau of the lowest velocity {u:.3f} m/s ({span}, {len(rows)} audit lines): "
          f"{(zm - zmin)*1e3:.2f} mm above the inlet plane (mid-column {lz*500:.1f} mm)")
    print(f"F5 particle counts seen {nset}; containment {'met' if zm < mid else 'NOT met'}: {verdict}{label}")
    if trans_txt and trans:
        rows_t = [row for ta, tb in trans for row in window(audit, ta, tb)]
        zt = max((row['zmax'] for row in rows_t if finite(row['zmax'])), default=NAN)
        print(f"   stricter, reported: including the {1e3*sum(tb - ta for ta, tb in trans):.0f} ms transition into it "
              f"({', '.join(f'{ta:.3f}-{tb:.3f} s' for ta, tb in trans)}): highest centre {(max(zm, zt) - zmin)*1e3:.2f} mm")


# ---------------------------------------------------------------------------------------------------
# self-test: synthetic logs
# ---------------------------------------------------------------------------------------------------
def _e(x, w=14, p=6):
    return ('%*s' % (w, 'NaN')) if x != x else '%*.*E' % (w, p, x)


PROTO = [(0.5, 0.30), (1.0, 0.20), (1.5, 0.10), (2.0, 0.05)]          # (start of step, U); end 2.5 s
FULL_STEPS = [(0.5, 1.0, 0.30), (1.0, 1.5, 0.20), (1.5, 2.0, 0.10), (2.0, 2.5, 0.05)]


def _dps(u):
    """dp* of the synthetic bed: fixed below 0.18 m/s, a plateau that falls slightly with U above"""
    return min(1.0, u / 0.18) * ((1.0 - 0.1 * (u - 0.2)) if u > 0.18 else 1.0)


def _write_run(path, protocol=PROTO, cad=0.01, times=None, t_from=None, t_end=None, area=2.25e-4, n=9789,
               width=2.5, margin_d=2.5, zlo_d=None, dyn=True, dp_hydro=0.0, d=0.87e-3, steps_file=None,
               ekern=0.37, ekern_fn=None, extents=True, z0=0.0, drop_last_knot=False, horizon='log',
               zmax_fn=None, legacy=False, tstep=0.5, t_run=None, planned=None):
    """a synthetic run directory: start-up echo, table echo, the time line of the first step (or the
    deck), EL_BED_AUDIT and EL_FORCE_BUDGET lines. Heights above the inlet plane z0: bed top 24 mm
    (zmax_fn(t) for the highest centre), mapped porosity ekern (or ekern_fn(t)); legacy = a binary before
    c8541cfc (no floor_bal_z, eps_kern, win_margin, gradz_win_dyn)."""
    os.makedirs(os.path.join(path, '_data'), exist_ok=True)
    rhop, rhof, g, dt_sim = 860.0, 1.0, 9.81, 2.5e-4
    vp = math.pi * d ** 3 / 6
    W = n * (rhop - rhof) * vp * g
    tab = [(0.0, protocol[0][1])]
    for (t, u) in protocol[1:]:
        tab += [(t, tab[-1][1]), (t + 0.02, u)]
    t_proto = protocol[-1][0] + tstep
    tab.append((t_proto, tab[-1][1]))
    if drop_last_knot:
        tab = tab[:-1]
    t_end = t_proto if t_end is None else t_end
    t_run = t_end if t_run is None else t_run
    lines = ['EL kernel                 = deen_poly', 'EL kernel width factor    = ' + _e(width),
             'EL fluid gravity          = ' + ('T' if dp_hydro else 'F')]
    lines += ['  EL_INLET_TABLE %s %s' % (_e(t, 16, 8), _e(u, 16, 8)) for t, u in tab]
    if extents:
        lines += ['Global CFD domain extents: ' + ' '.join(_e(x) for x in (0.0, 0.015, 0.0, 0.015, z0, z0 + 0.1))]
    nplan = int(round((t_proto if planned is None else planned) / dt_sim))
    if horizon == 'log':
        lines.append('time:  %s/  %s | itns:       1/  %d  | dt:  %s' % (
            _e(dt_sim, 10, 4).replace('E', 'D'), _e(t_proto, 10, 4).replace('E', 'D'), nplan,
            _e(dt_sim, 10, 4).replace('E', 'D')))
    if horizon == 'deck':
        with open(os.path.join(path, '_data', 'q2p1_param.dat'), 'w') as fh:
            fh.write('SimPar@TimeStep = %gd0\nSimPar@MaxNumStep = %d\nSimPar@StartSimTime = 0d0\n'
                     % (dt_sim, nplan))
    if times is None:
        times, k = [], 1
        while k * cad <= t_end + 1e-12:
            times.append(k * cad)
            k += 1
    times = [t for t in times if (t_from is None or t >= t_from - 1e-12) and t <= t_end + 1e-12]
    for t in times:
        u = u_of_t(tab, t)
        dps = _dps(u)
        dp = dps * div(W, area) + dp_hydro if (area > 0 and n > 0) else NAN
        z99 = z0 + 0.024 if n > 0 else NAN
        m = margin_d * d
        zlo = z0 + (zlo_d if zlo_d is not None else margin_d) * d if n > 0 else NAN
        zhi = z99 + d / 2 - m if n > 0 else NAN
        zmx = (z0 + (zmax_fn(t) if zmax_fn else 0.0241)) if n > 0 else NAN
        grad = dp / 0.1 if dp == dp else NAN
        gdyn = grad - dp_hydro / 0.1 if dp == dp else NAN
        ek = (ekern_fn(t) if ekern_fn else ekern) if n > 0 else NAN
        line = ('EL_BED_AUDIT t= %s step= %d p_in_out= %s %s area_in_out= %s %s dp_hydro_dyn= %s %s %s n= %d '
                'weight_z_floor_z= %s %s z50_95_99_max= %s %s %s %s eps_bed= %s window= %s %s gradz_win= %s' % (
                    _e(t), int(round(t / dt_sim)), _e(dp, 16, 8), _e(0.0, 16, 8), _e(area), _e(area), _e(dp, 16, 8),
                    _e(dp_hydro, 16, 8), _e(dp - dp_hydro if dp == dp else NAN, 16, 8), n,
                    _e(-W if n > 0 else 0.0), _e(0.0), _e(z99 - 0.012), _e(z99 - 0.001), _e(z99), _e(zmx),
                    _e(0.37 if n > 0 else NAN), _e(zlo), _e(zhi), _e(grad, 16, 8)))
        if not legacy:
            line += ' floor_bal_z= %s eps_kern_slip_win= %s %s' % (_e((1 - dps) * W, 16, 8), _e(ek, 16, 8),
                                                                  _e(0.1 if n > 0 else NAN, 16, 8))
            if dyn:
                line += ' win_margin= %s gradz_win_dyn= %s' % (_e(m if n > 0 else NAN), _e(gdyn, 16, 8))
        lines.append(line)
        lines.append('EL_FORCE_BUDGET step= %d drag= %s %s %s pressure= %s %s %s lift= 0 0 0 grav_buoy= %s %s %s' % (
            int(round(t / dt_sim)), _e(0.0), _e(0.0), _e(dps * W * 0.9), _e(0.0), _e(0.0), _e(dps * W * 0.1),
            _e(0.0), _e(0.0), _e(-W)))
    if horizon == 'log':
        lines.append('time:  %s/  %s | itns:  %d/  %d  | dt:  %s' % (
            _e(t_run, 10, 4).replace('E', 'D'), _e(t_proto, 10, 4).replace('E', 'D'), int(round(t_run / dt_sim)),
            nplan, _e(dt_sim, 10, 4).replace('E', 'D')))
    with open(os.path.join(path, 'run_slurm.log'), 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    if steps_file is not None:
        with open(os.path.join(path, '_data', 'steps.txt'), 'w') as fh:
            fh.write('# t_start t_end U\n' + ''.join('%.6g %.6g %.6g\n' % s for s in steps_file))
    return path


GATE_LINE = re.compile(r'^(F1 plateau|F2 worst|F2 verdict|   U_mf at dp\* = 1|F5 particle counts)')


def _gates(text):
    return [l for l in text.splitlines() if GATE_LINE.match(l)]


def self_test_cases():
    """(name, generator keywords, extra arguments, strings that must appear, strings that must not,
    True when the gate lines must equal those of the reference case). The cases marked r3 are the
    counterexamples of the third review."""
    S = '{run}/_data/steps.txt'
    hi_after = lambda t: 0.060 if t > 2.3 else 0.0241       # noqa: E731  a centre at 60 mm late in the lowest plateau
    return [
        ('reference (all gates measurable)', dict(), [], ['F1 plateau', 'F3 fixed branch', 'F2 worst |ratio_k - 1|',
                                                          'F5 particle counts seen [9789]', 'PASS', 'stricter, reported',
                                                          '(a) particle count over the run: counts seen [9789] on 250 audit lines, run 0.0000-2.5000 s covered'],
         ['unavailable', 'Traceback', 'incomplete'], False),
        ('missing inlet plane (area_in = 0)', dict(area=0.0), [], ['W/A unavailable (no inlet plane', 'F1, F2, F3 unavailable',
                                                                   'F5 particle counts seen'], ['Traceback'], False),
        ('empty bed (n = 0)', dict(n=0), [], ['empty bed', 'F1, F2, F3 unavailable', 'F5 unavailable'], ['Traceback', 'PASS'], False),
        ('r2: steps file without the lowest step', dict(steps_file=FULL_STEPS[:3]), ['--steps', S],
         ['plateau (2.0000, 2.5000, 0.0500) missing from the list', 'steps-file summary, NOT gated',
          'lowest velocity 0.050 m/s'], ['Traceback'], True),
        ('r2: foreign step in the steps file', dict(steps_file=FULL_STEPS + [(2.5, 3.0, 0.01)]), ['--steps', S],
         ['(2.5000, 3.0000, 0.0100) is no plateau of the protocol', 'lowest velocity 0.050 m/s'], ['Traceback', '0.010 m/s'], True),
        ('r3: steps file with a duplicated row', dict(steps_file=FULL_STEPS + [FULL_STEPS[0]]), ['--steps', S],
         ['listed 2 times', 'steps-file summary, NOT gated'], ['Traceback'], True),
        ('r3: steps file with an omitted row', dict(steps_file=FULL_STEPS[1:]), ['--steps', S],
         ['plateau (0.5000, 1.0000, 0.3000) missing from the list', 'steps-file summary, NOT gated'], ['Traceback'], True),
        ('r3: steps file with a shifted start', dict(steps_file=[FULL_STEPS[0], (1.2, 1.5, 0.20)] + FULL_STEPS[2:]),
         ['--steps', S], ['differs from the plateau (1.0000, 1.5000, 0.2000) in the start time', 'steps-file summary, NOT gated'],
         ['Traceback'], True),
        ('steps file identical to the protocol', dict(steps_file=FULL_STEPS), ['--steps', S],
         ['identical to the echoed protocol'], ['NOT gated', 'Traceback'], True),
        ('short window (run ends inside the lowest plateau)', dict(t_end=2.3), [],
         ['has not reached the protocol end', 'incomplete: the run did not reach 2.5000 s', 'F5 not evaluated: the audit does not cover'],
         ['Traceback', 'PASS (', 'containment met: PASS'], False),
        ('sparse audit (cadence wider than a tenth of the window)', dict(cad=0.4), [], ['F1, F2, F3 unavailable', 'F5 not evaluated'],
         ['Traceback'], False),
        ('window margin 2 d with a 2.5 d support', dict(margin_d=2.0), [], ['NOT met', 'F2 verdict unavailable (window inside the mapping support)',
                                                                            'F5 particle counts seen'], ['Traceback'], False),
        ('lower window end at 2 d, field says 2.5 d', dict(zlo_d=2.0), [], ['F2 verdict unavailable (window inside the mapping support)'],
         ['Traceback'], False),
        ('old log: no win_margin / gradz_win_dyn, zero hydrostatic part', dict(margin_d=2.0, dyn=False), [],
         ['hydrostatic part zero', 'margin 2 d (binary before 627caf47)', 'F2 verdict unavailable (window inside the mapping support)'],
         ['Traceback'], False),
        ('old log with a 2 d support: margin met, total gradient admitted', dict(margin_d=2.0, dyn=False, width=2.0), [],
         ['hydrostatic part zero', 'F2 worst |ratio_k - 1|'], ['F2 verdict unavailable', 'Traceback'], False),
        ('old log with a hydrostatic part', dict(margin_d=2.0, dyn=False, width=2.0, dp_hydro=0.981), [],
         ['nonzero hydrostatic part', 'F2 verdict unavailable'], ['Traceback'], False),
        ('r3: mapped porosity 1.0 on the fixed steps (closure zero)', dict(ekern=1.0), [],
         ['mapped gas-fraction sample(s)', 'no positive closure and finite ratio_k', 'F2 verdict unavailable'],
         ['F2 worst', 'Traceback'], False),
        ('r3: the -1 sentinel of an old binary in the mapped samples', dict(ekern_fn=lambda t: -1.0 if 2.3 < t < 2.4 else 0.37), [],
         ['mapped gas-fraction sample(s)', 'F2 verdict unavailable'], ['F2 worst', 'Traceback'], False),
        ('r3: zero-velocity fixed plateau', dict(protocol=PROTO + [(2.5, 0.0)]), [],
         ['no positive closure and finite ratio_k on every fixed step: step 4 (U 0.000)', 'F2 verdict unavailable'],
         ['F2 worst', 'Traceback'], False),
        ('r3: terminal hold without its last knot, centre at 60 mm late in the hold', dict(drop_last_knot=True, zmax_fn=hi_after), [],
         ['lowest velocity 0.050 m/s (2.020-2.500 s', 'containment NOT met: FAIL'], ['0.100 m/s (', 'Traceback'], False),
        ('r3: the same table, horizon from the deck', dict(drop_last_knot=True, zmax_fn=hi_after, horizon='deck'), [],
         ['[deck: min(MaxSimTime', 'lowest velocity 0.050 m/s (2.020-2.500 s', 'containment NOT met: FAIL'], ['Traceback'], False),
        ('r3: the same table without a horizon', dict(drop_last_knot=True, zmax_fn=hi_after, horizon=None), [],
         ['F5 withheld: the lowest velocity 0.050 m/s is the hold of the last table value'], ['containment', 'Traceback'], False),
        ('MaxNumStep ends the run before the table (horizon from the step count)', dict(planned=2.3, t_end=2.3), [],
         ['horizon 2.3000 s', 'lowest velocity 0.050 m/s (2.020-2.300 s', 'containment met: PASS'], ['Traceback', 'incomplete'], False),
        ('r3: one audit line of a binary without floor_bal_z', dict(legacy=True, times=[0.75]), [],
         ['1 audit lines', 'F1, F2, F3 unavailable', 'F5'], ['Traceback'], False),
        ('r3: no column extents, lower window end at 1 d', dict(zlo_d=1.0, extents=False), [],
         ['missing geometry: no inlet-plane elevation', 'F2 verdict unavailable', 'F5 unavailable: no column extents'],
         ['F2 worst', 'Traceback', 'met: PASS'], False),
        ('r3: the last audit line dropped (the run reached the end)', dict(t_end=2.49, t_run=2.5), [],
         ['containment met: PASS', 'F1 plateau'], ['has not reached', 'not evaluated', 'incomplete', 'Traceback'], False),
        ('r3: 40 ms samples in one window of a 30 ms audit', dict(cad=0.03, times=sorted(
            [0.03 * k for k in range(1, 84) if not 1.68 < 0.03 * k < 2.0 + 1e-9] + [1.72 + 0.04 * k for k in range(8)]), t_run=2.5), [],
         ['incomplete: gap of'], ['Traceback'], False),
        ('r3: uniform 30 ms audit (every window qualifies)', dict(cad=0.03, t_run=2.5), [],
         ['F1 plateau', 'F2 worst'], ['incomplete', 'Traceback'], False),
        ('r3: only the last 0.48 s of the run in the log', dict(t_from=2.02), [],
         ['NOT covered', 'unavailable (count over the run not established)'], ['containment met: PASS', 'Traceback'], False),
        ('r3: column translated down by 20 mm, centre 60 mm above the floor late', dict(z0=-0.020, zmax_fn=hi_after), [],
         ['60.00 mm above the inlet plane', 'containment NOT met: FAIL', 'whole-bed porosity 0.3'], ['Traceback'], False),
    ]


def self_test(keep=None):
    """runs every case; keep = a directory that receives the synthetic runs (for a comparison with
    another version of the tool), else a temporary directory"""
    cases = self_test_cases()
    nfail = 0
    ref_gates = None
    ctx = tempfile.TemporaryDirectory(prefix='fl11_post_selftest_') if keep is None else contextlib.nullcontext(keep)
    with ctx as tmp:
        for i, (name, kw, extra, must, mustnot, same_gates) in enumerate(cases):
            run = _write_run(os.path.join(tmp, 'run%02d' % i), **kw)
            out = io.StringIO()
            try:
                with contextlib.redirect_stdout(out):
                    main([run] + [x.format(run=run) for x in extra])
            except BaseException as exc:                 # noqa: BLE001  (a crash is a failed check)
                import traceback
                out.write('Traceback\n' + traceback.format_exc())
            text = out.getvalue()
            if i == 0:
                ref_gates = _gates(text)
            bad = [m for m in must if m not in text] + [m for m in mustnot if m in text]
            if same_gates and _gates(text) != ref_gates:
                bad.append('gate lines differ from the reference case: %s' % _gates(text))
            print(f"self-test {i:2d} {'ok  ' if not bad else 'FAIL'} {name}" + ('' if not bad else f": {bad}"))
            if bad:
                nfail += 1
                print('    ' + text.replace('\n', '\n    '))
            with open(os.path.join(run, 'case.txt'), 'w') as fh:
                fh.write('%s\n%s\n' % (name, ' '.join(x.format(run=run) for x in extra)))
    print(f"fl11_post self-test: {len(cases) - nfail}/{len(cases)} passed")
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
