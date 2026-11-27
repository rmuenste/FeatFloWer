#!/usr/bin/env python3
"""Spec-11 (minimum fluidisation) gates of a q2p1_el_pipeflow run.

The tool gates spec 11 and nothing else: it RECOGNISES a spec-11 run (design note 7.4). The spec is one
`Spec` object (inlet table, physical values, gate constants); the spec-11 instance takes its table from
make_inlet_table.py (the function that writes the deck's table) and its physical values from the spec.
Official gates are issued only when recognise() returns no reason:
  - the echoed inlet table equals the spec table (row count, every time and velocity as decimals to the
    print precision of the echo);
  - the effective physical inputs (the last start-up assignment in echo order; the viscosity cross-checked
    against the conversion echoes) equal the spec: d, rho_p, rho_f, mu, g, drag model, kernel width factor
    (echo, else a readable deck), column extents;
  - a fresh start (StartingProc = 0, no dump load), exactly one invocation (one start-up banner, each
    start-up echo once, nothing the tool reads after the finish marker), no reload line of any kind;
  - one constant-dt clock t0 + k dt from the first table row (deck TimeStep / StartSimTime as decimals,
    every time line and audit row on it to its print precision);
  - finished (finish marker) at or after the last table row.
All time comparisons are exact (fractions of the printed decimals and of the deck strings). Outside the
domain there is no official verdict; a witnessed FAIL (F5 count change, F5 height exceedance inside the
spec's lowest hold) is kept only when the run fails nothing but "finished" (an unfinished or still running
spec-11 run), otherwise it is printed as a measurement.

Gates (spec section 5, amended 2026-11-20; F1 amended 2026-11-24, after the 80-row result had been seen):
F1 plateau = mean dp* over the measured steps with U above the run's own measured U_mf (the F3 crossing of
the fixed-branch fit with dp* = 1), 1 +- 0.02; the sets come from an iteration (branch_sets: start from the
steps with dp* > 0.95, fixed = dp* < 0.8 of the plateau and not in the fluidised set, fit, U_mf, fluidised
= U > U_mf, until both sets are unchanged; no F1-F3 verdict when they do not settle or contradict each
other, or when a step velocity coincides with U_mf within the resolution of the log); the mean over the
steps with dp* > 0.95 is a reported cross-check; a plateau within the standard error of its mean (sample
standard deviation of the step means over sqrt(n), n >= 3) of the nearer limit is annotated "at the gate"
next to the PASS / FAIL (f1_gate). F2 fixed-branch window gradient against the closure at the mapped
gas fraction +- 5 %, window margin max(2 d, support); F3 U_mf of the fit dp* = a U + b U^2 (two distinct
nonzero fixed velocities, full rank) against the closure U_mf at the whole-bed porosity +- 5 %; F4
reported; F5 count constant over the run, height: no centre above zmin + Lz/2 in the lowest hold.

usage: fl11_post.py RUNDIR [--series OUT.csv]
       fl11_post.py --self-test                 the machinery on synthetic logs with fixture specs
       fl11_post.py --self-test-real RUNDIR [--variants N] [--seed S]   metamorphic tests on a real log
fl11_movie.py uses evaluate() of this module."""
import argparse
import contextlib
import io
import math
import os
import random
import re
import sys
import tempfile
from decimal import Decimal, localcontext, ROUND_HALF_EVEN
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fl11_closures import CLOSURES, closure_grad, closure_umf   # noqa: E402
import make_inlet_table                                           # noqa: E402

EVAL_VERSION = 'fl11_post r9 2026-11-27'   # part of the movie cache signature
NAN = float('nan')
INF = float('inf')
GATES = ('F1', 'F2', 'F3', 'F5 height', 'F5 count')
NOT_OFFICIAL = 'not a verdict (outside the supported domain)'
OUTSIDE = ' (outside the supported domain)'
AT_GATE_NOTE = ' (at the gate)'            # F1 annotation of the verdict line (the verdict stays PASS / FAIL)
FINISH = 'PP3D_LES has successfully finished'
BANNER = 'Parellel Q2/P1 FEM Fluid Dynamics code'


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def pos(x):
    return x if finite(x) and x > 0 else NAN


def div(a, b):
    b = pos(b)
    return a / b if finite(b) and finite(a) else NAN


def es_prec(x, digits=6):
    """half the spacing of the ES(w).digits print of x (0 for an exact zero)"""
    if not finite(x) or x == 0:
        return 0.0
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(x))) - digits)


def mean(v):
    return sum(v) / len(v) if v else NAN


def exact_mean(v):
    """the correctly rounded mean of finite floats (exact rational sum): equal values give that value"""
    return float(sum(map(Fraction, v)) / len(v)) if v and all(finite(x) for x in v) else NAN


def _fmin(v):
    v = [x for x in v if finite(x)]
    return min(v) if v else NAN


def _fmax(v):
    v = [x for x in v if finite(x)]
    return max(v) if v else NAN


def _f(x, spec):
    return format(x, spec) if finite(x) else format('n/a', '>' + re.match(r'\d*', spec).group(0))


# ---- printed decimals ------------------------------------------------------------------------------
NUMG = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?'
NUM_RE = re.compile('^' + NUMG + '$')


def num(tok):
    u = tok.upper()
    if u in ('NAN', '+NAN', '-NAN'):
        return NAN
    if u in ('INF', '+INF', 'INFINITY', '+INFINITY'):
        return INF
    if u in ('-INF', '-INFINITY'):
        return -INF
    if re.fullmatch(r'\*+', tok):
        return NAN
    if not NUM_RE.match(tok):
        raise ValueError(tok)
    return float(tok.replace('D', 'E').replace('d', 'e'))


def integer(tok):
    if not re.fullmatch(r'[-+]?\d+', tok):
        raise ValueError(tok)
    return int(tok)


def dec(s):
    """the decimal a printed number stands for (Fortran D exponent accepted); ValueError otherwise"""
    s = s.strip().strip('"\'')
    if not NUM_RE.match(s):
        raise ValueError(s)
    return Decimal(s.replace('D', 'E').replace('d', 'e'))


def frac(s):
    return Fraction(dec(s))


def half_ulp(s):
    """half the spacing of the printed digits of s (exact)"""
    return Fraction(1, 2) * Fraction(10) ** dec(s).as_tuple().exponent


MIN_SIG = 4      # a printed token shorter than this many significant digits must be the reference exactly


def print_equal(s, ref):
    """the printed number s stands for the decimal ref: s is ref exactly as a decimal, or s carries at least MIN_SIG
    significant digits and equals ref rounded to its digits (a coarse token such as '0' or '1' never stands for
    0.3 or 0.9)"""
    try:
        p = dec(s)
    except ValueError:
        return False
    if not p.is_finite():
        return False
    with localcontext() as c:
        c.prec = 80
        r = Decimal(ref)
        if p == r:
            return True
        digits = p.as_tuple().digits
        sig = len(digits) - next((i for i, x in enumerate(digits) if x), len(digits))
        return sig >= MIN_SIG and r.quantize(Decimal(1).scaleb(p.as_tuple().exponent), rounding=ROUND_HALF_EVEN) == p


# ===================================================================================================
# the spec
# ===================================================================================================
class Spec:
    """what a run must reproduce to be gated, and the gate constants; defined once"""
    HOLD_MIN = Fraction(1, 10)      # s: a hold of this length or less belongs to a transition
    EVAL_LEN = Fraction(3, 10)      # s: evaluation window at the end of a step
    FIRST_LEN, FIRST_LONG = Fraction(1, 2), Fraction(3, 5)   # the first hold: its last 0.5 s when longer than 0.6 s
    SUB = (Fraction(1, 5), Fraction(7, 20), Fraction(1, 2))  # sub-windows 0.2-0.35-0.5 s after the step start
    MIN_SAMPLES, MIN_SUB = 10, 5    # expected audit rows in a window / in a sub-window
    FLU, FIX = 0.95, 0.8            # start set of the iteration: dp* > FLU; fixed: dp* < FIX x plateau (and not fluidised)
    ROUNDS = 10                     # branch_sets: the sets must settle within this many rounds
    F1_LIMITS = (0.98, 1.02)        # F1: the gate 1 +- 0.02 as its two limits (a plateau on a limit passes)
    TOL_F2, TOL_F3 = 0.05, 0.05
    # branch_sets: a step velocity within this relative distance of U_mf is neither above nor not above it. dp* is
    # the logged pressure drop over the logged weight per area; the audit line prints the weight with seven
    # significant digits (relative up to 5e-7, common to all steps, and so in U_mf) and the pressures with nine.
    # The band is twice that bound. Round-off of the fit (1e-15) is far inside it.
    UMF_TIE = 1e-6
    UMF_MEASURED = 0.182            # m/s, NETL measurement (F4, reported)
    FLAG = 0.02                     # sub-window difference flagged in the step table (display)
    REL_DET = 1e-6                  # F3: smallest relative determinant of the normal equations
    # set-up checks (they qualify the verdicts, they do not withhold them): every contact pair of the
    # EL_CONTACT_PARAMS echo, restitution switched on; the largest overlap of EL_CONTACT_STATS
    CONTACT = dict(mu_static=Decimal('0.3'), mu_dynamic=Decimal('0.3'), e=Decimal('0.9'))
    RESTITUTION = 'on'
    OVERLAP_MAX = Decimal('0.01')   # x d

    def __init__(self, name, table, d, rhop, rhof, mu, g, drag, width, extents, N=9789):
        self.name = name
        self.N = int(N)                 # the bed-R population at start-up (the sphere-creation echo)
        self.gravity = (Decimal(0), Decimal(0), -Decimal(g))
        self.table = [(Decimal(t), Decimal(u)) for t, u in table]
        self.knots = [(Fraction(t), float(u)) for t, u in self.table]
        self.decimal = dict(d=Decimal(d), rhop=Decimal(rhop), rhof=Decimal(rhof), mu=Decimal(mu), g=Decimal(g),
                            width=Decimal(width))
        self.d, self.rhop, self.rhof, self.mu, self.g, self.width = (float(self.decimal[k]) for k in
                                                                     ('d', 'rhop', 'rhof', 'mu', 'g', 'width'))
        self.drag = drag
        self.extents = tuple(Decimal(x) for x in extents)
        self.zmin, self.zmax = float(self.extents[4]), float(self.extents[5])
        if drag not in CLOSURES or len(self.table) < 2 or any(not q[0] > p[0] for p, q in zip(self.knots, self.knots[1:])):
            raise ValueError('spec %s: drag model or table not usable' % name)
        self.holds = holds_of(self.table, self.HOLD_MIN)
        u = min(h['u'] for h in self.holds) if self.holds else None
        low = [h for h in self.holds if h['u'] == u]
        if len(low) != 1 or low[0] is self.holds[0]:
            raise ValueError('spec %s: the lowest hold must be unique and follow the first hold' % name)
        self.lowest = low[0]
        self.steps = steps_of(self.holds, self)
        self.t_first, self.t_last = self.knots[0][0], self.knots[-1][0]


def holds_of(table, hold_min):
    """maximal constant-velocity intervals between the table rows (equal-valued segments merged, then
    holds of hold_min or less dropped); exact: a, b Fractions, u the Decimal velocity, uf its float"""
    merged = []
    for (ta, ua), (tb, ub) in zip(table, table[1:]):
        if ua != ub:
            merged.append(None)
        elif merged and merged[-1] is not None and merged[-1][2] == ua:
            merged[-1] = (merged[-1][0], Fraction(tb), ua)
        else:
            merged.append((Fraction(ta), Fraction(tb), ua))
    return [dict(a=h[0], b=h[1], u=h[2], uf=float(h[2])) for h in merged if h is not None and h[1] - h[0] > hold_min]


def steps_of(holds, spec):
    steps = []
    for i, h in enumerate(holds):
        t1 = h['b']
        t0 = holds[i - 1]['b'] if i > 0 else h['a']
        if i == 0 and t1 - t0 > spec.FIRST_LONG:
            t0 = t1 - spec.FIRST_LEN
        steps.append(dict(k=i, t0=t0, t1=t1, U=h['uf'], hold=h))
    return steps


def spec11():
    pts, _, _ = make_inlet_table.protocol(**make_inlet_table.DEFAULTS)
    return Spec('spec 11', [r.split() for r in make_inlet_table.deck_rows(pts)], d='0.87e-3', rhop='860',
                rhof='1.0', mu='1.8e-5', g='9.81', drag='gidaspow', width='2.5',
                extents=('0', '0.015', '0', '0.015', '0', '0.1'))


SPEC11 = spec11()


def u_of_t(knots, t):
    """the solver's rule: linear between the rows, constant beyond the ends (display, the U column)"""
    if t <= knots[0][0]:
        return knots[0][1]
    for (t0, u0), (t1, u1) in zip(knots, knots[1:]):
        if t <= t1:
            return u0 + (u1 - u0) * float((t - t0) / (t1 - t0))
    return knots[-1][1]


# ===================================================================================================
# record reader: one pass over the log
# ===================================================================================================
SPECS = {
    'EL_BED_AUDIT': ([('t', 1, num), ('step', 1, integer), ('p_in_out', 2, num), ('area_in_out', 2, num),
                      ('dp_hydro_dyn', 3, num), ('n', 1, integer), ('weight_z_floor_z', 2, num),
                      ('z50_95_99_max', 4, num), ('eps_bed', 1, num), ('window', 2, num)],
                     [('gradz_win', 1, num), ('floor_bal_z', 1, num), ('eps_kern_slip_win', 2, num),
                      ('win_margin', 1, num), ('gradz_win_dyn', 1, num)]),
    'EL_FORCE_BUDGET': ([('step', 1, integer), ('drag', 3, num), ('pressure', 3, num), ('lift', 3, num),
                         ('grav_buoy', 3, num)], [('total', 3, num)]),
    'EL_NEWTON_PAIR': ([('t', 1, num), ('step', 1, integer), ('dp', 3, num), ('expected', 3, num)],
                       [('mismatch', 3, num)]),
    'EL_CONTACT_STATS': ([('t', 1, num), ('ncontacts', 1, integer), ('max_overlap', 1, num)], [('Tgran', 1, num)]),
}
KEY_TOK = re.compile(r'^[A-Za-z_0-9]+=$')


def parse_keyed(tokens, spec):
    """(record, '') for a complete record, (None, reason) otherwise; record['t_s'] keeps the printed time"""
    req, opt = spec
    groups, cur = {}, None
    for tok in tokens:
        if KEY_TOK.match(tok):
            cur = tok[:-1]
            if cur in groups:
                return None, 'key %s= twice' % cur
            groups[cur] = []
        elif cur is None:
            return None, 'value before the first key'
        else:
            groups[cur].append(tok)
    rec = {}
    for i, (key, cnt, conv) in enumerate(req + opt):
        if key not in groups:
            if i < len(req):
                return None, 'key %s= missing' % key
            continue
        vals = groups[key]
        if len(vals) != cnt:
            return None, '%s= has %d value(s), %d expected' % (key, len(vals), cnt)
        try:
            rec[key] = [conv(v) for v in vals]
        except ValueError:
            return None, '%s= value not readable' % key
    if 't' in groups:
        rec['t_s'] = groups['t'][0]
    return rec, ''


PAT_TIME = re.compile(r'^\s*time:\s*(%s)\s*/\s*(%s)\s*\|\s*itns:\s*(\d+)\s*/\s*(\d+)\s*\|\s*dt:\s*(%s)'
                      % (NUMG, NUMG, NUMG))
PAT_EL = re.compile(r'^EL ([A-Za-z][^=:]*?)\s+=\s?(.*)$')            # start-up configuration echo
PAT_PROP = re.compile(r'^Prop\s+-\s+(.+?)\s+=\s+(.*)$')
PAT_SIMPAR = re.compile(r'^(StartingProc|TimeStep|StartSimTime|MaxSimTime|MaxNumStep)\s+=\s*(\S+)')
PAT_SPHERES = re.compile(r'^EL terminal-velocity setup:\s*(\d+) free sphere\(s\) created \(r=([^,\s]+), rho=([^,\s]+)')
PAT_NU = re.compile(r'^EL fluid viscosity: diffusion operator rebuilt with mu =\s*(\S+)\s*\(nu =\s*(\S+) for the closures\)')
PAT_EXT = re.compile(r'^Global CFD domain extents:\s*(.*)$')


class Log:
    """everything the tool reads from run_slurm.log, in order, from one pass"""

    def __init__(self, lines, last_complete=True):
        self.el = {}                     # start-up echo key -> [(line, value)]
        self.prop = []                   # (line, key, [values]) in echo order
        self.simpar = {}                 # key -> [(line, value)]
        self.spheres, self.nu_el, self.extents = [], [], []
        self.banners, self.reloads, self.dumps = [], [], []
        self.table, self.audit, self.times = [], [], []    # table rows (t, U strings); audit (line, rec, raw); time lines
        self.budget, self.pairs = {}, {}
        self.contact_params, self.contact_stats = [], []   # start-up echo lines; (line, rec, raw)
        self.finish = None               # line of the first finish marker
        self.after_finish = None         # first line the tool reads after it: (line, text)
        self.first_time = None           # line of the first time line (end of the start-up)
        bad = []
        good_last = -1
        n = len(lines)
        for i, line in enumerate(lines):
            unfinished = (i == n - 1 and not last_complete)
            s = line.rstrip('\r\n')
            st = s.strip()
            kind = None
            if st.startswith('EL_CONTACT_PARAMS'):
                self.contact_params.append((i, st))
                kind = 'start-up echo'
            elif st.startswith('EL_INLET_TABLE'):
                tok = st.split()
                if unfinished or len(tok) != 3 or not (NUM_RE.match(tok[1]) and NUM_RE.match(tok[2])):
                    bad.append((i, 'EL_INLET_TABLE', 'no line end' if unfinished else 'not two numbers t U'))
                else:
                    self.table.append((tok[1], tok[2]))
                    good_last, kind = i, 'record'
            elif s.split(' ', 1)[0] in SPECS:
                tag = s.split(' ', 1)[0]
                rec, why = (None, 'no line end') if unfinished else parse_keyed(s.split()[1:], SPECS[tag])
                if rec is None:
                    bad.append((i, tag, why))
                else:
                    good_last, kind = i, 'record'
                    if tag == 'EL_BED_AUDIT':
                        self.audit.append((i, rec, ' '.join(s.split()[1:])))
                    elif tag == 'EL_CONTACT_STATS':
                        self.contact_stats.append((i, rec, ' '.join(s.split()[1:])))
                    else:
                        store = self.budget if tag == 'EL_FORCE_BUDGET' else self.pairs
                        k = rec['step'][0]
                        if k not in store:
                            store[k] = rec
                        elif store[k] is not None and store[k] != rec:
                            store[k] = None
            elif st.startswith('time:'):
                m = PAT_TIME.match(s)
                if m and not unfinished:
                    self.times.append((m.group(1), m.group(2), int(m.group(3)), int(m.group(4)), m.group(5)))
                    if self.first_time is None:
                        self.first_time = i
                    good_last, kind = i, 'time line'
                else:
                    bad.append((i, 'time line', 'no line end' if unfinished else 'not readable'))
            elif FINISH in s and not unfinished:
                if self.finish is None:
                    self.finish = i
                    kind = None
                else:
                    kind = 'finish marker'
                good_last = i
            elif BANNER in s:
                self.banners.append(i)
                kind = 'start-up banner'
            elif st.startswith('Reloading '):
                self.reloads.append((i, st))
                kind = 'reload'
            elif 'Loading particle data from' in st:
                self.dumps.append((i, st))
                kind = 'dump load'
            else:
                m = PAT_EL.match(s) or PAT_PROP.match(s) or PAT_SIMPAR.match(st)
                if m and m.re is PAT_EL:
                    self.el.setdefault(m.group(1).strip(), []).append((i, m.group(2).strip()))
                    kind = 'start-up echo'
                elif m and m.re is PAT_PROP:
                    self.prop.append((i, m.group(1).strip(), m.group(2).split()))
                    kind = 'property echo'
                elif m:
                    self.simpar.setdefault(m.group(1), []).append((i, m.group(2)))
                    kind = 'start-up echo'
                else:
                    for pat, store in ((PAT_SPHERES, self.spheres), (PAT_NU, self.nu_el), (PAT_EXT, self.extents)):
                        m = pat.match(st)
                        if m:
                            store.append((i, m.groups()))
                            kind = 'start-up echo'
                            break
            if kind and self.finish is not None and i > self.finish and self.after_finish is None:
                self.after_finish = (i, '%s: %s' % (kind, st[:60]))
        self.trailing = [b for b in bad if b[0] > good_last]
        self.malformed = [b for b in bad if b[0] <= good_last]

    def el_last(self, key):
        v = self.el.get(key)
        return v[-1][1] if v else None


def read_deck(text):
    """SimPar / Prop keys of a deck with the application's syntax (last '@', last '=', first token; a later
    line overrides an earlier one); values: lists of text"""
    keys = {}
    for line in text.splitlines():
        if '@' not in line or '=' not in line:
            continue
        iat, ieq = line.rfind('@'), line.rfind('=')
        if ieq < iat:
            continue
        sec = line[:iat].split()
        key = line[iat + 1:ieq].split()
        if len(sec) != 1 or sec[0] not in ('SimPar', 'Prop') or not key:
            continue
        val = [x for x in re.split(r'[\s,/]+', line[ieq + 1:].strip()) if x]
        keys[key[0] if sec[0] == 'SimPar' else 'Prop@' + key[0]] = val[:1] if sec[0] == 'SimPar' else val
    return keys


def deck_str(keys, key):
    v = keys.get(key)
    return v[0].strip('"\'') if v else None


# ===================================================================================================
# clock and audit rows (exact)
# ===================================================================================================
class Clock:
    """t0 + k dt from the deck strings, checked against every time line; the validated audit rows"""

    def __init__(self, log, deck, deck_ok, spec, lz):
        self.segment = []
        f = log.el_last('bed audit freq')
        self.freq = int(f) if f is not None and re.fullmatch(r'\d+', f) and int(f) > 0 else None
        self.dt = self.t0 = None
        try:
            self.dt, self.t0 = frac(deck_str(deck, 'TimeStep') or 'x'), frac(deck_str(deck, 'StartSimTime') or 'x')
        except ValueError:
            pass
        self.exact = self.dt is not None and self.dt > 0
        if not deck_ok:
            self.segment.append('deck _data/q2p1_param.dat not readable (no exact clock)')
        elif not self.exact:
            self.segment.append('deck SimPar@TimeStep / SimPar@StartSimTime missing or not a positive decimal')
        times = log.times
        if not times:
            self.segment.append('no time lines in the log')
        dts = sorted({dec(x[4]) for x in times})
        if len(dts) > 1:
            self.segment.append('time step not constant (printed dt %s)' % ', '.join('%.4g' % x for x in dts))
        elif dts and self.exact and not print_equal(times[0][4], self.dt.numerator / Decimal(self.dt.denominator)):
            self.segment.append('printed dt %s differs from the deck TimeStep %s' % (times[0][4], deck_str(deck, 'TimeStep')))
        for key, val in (('TimeStep', self.dt), ('StartSimTime', self.t0)):
            for _, s in log.simpar.get(key, [])[:1]:
                if self.exact and not print_equal(s, val.numerator / Decimal(val.denominator)):
                    self.segment.append('start-up echo %s = %s differs from the deck %s' % (key, s, deck_str(deck, key)))
        ks = [x[2] for x in times]
        consecutive = ks == list(range(1, len(ks) + 1))
        if times and not consecutive:
            j = next((i for i, k in enumerate(ks) if k != i + 1), 0)
            self.segment.append('the time lines are not the consecutive steps 1, 2, ... (line %d shows step %d)' % (j + 1, ks[j]))
        if len({(x[1], x[3]) for x in times}) > 1:
            self.segment.append('MaxSimTime or MaxNumStep changed during the run')
        if self.exact:
            for t, _, k, _, _ in times:
                tk = self.t0 + k * self.dt
                if abs(frac(t) - tk) > half_ulp(t):
                    self.segment.append('the time line of step %d shows t = %s s, the clock t0 + k dt gives %.6g s'
                                        % (k, t, float(tk)))
                    break
            if self.t0 != spec.t_first:
                self.segment.append('the clock starts at %.4f s, the protocol (first table row) at %.4f s'
                                    % (float(self.t0), float(spec.t_first)))
        self.K_time = ks[-1] if times and consecutive else max(ks, default=0)
        self.finished = log.finish is not None
        self.validate(log, lz)

    def tk(self, k):
        return self.t0 + k * self.dt if self.exact else self.row_t.get(k)

    def validate(self, log, lz):
        groups = {}
        for ln, rec, raw in log.audit:
            groups.setdefault(rec['step'][0], []).append((ln, rec, raw))
        self.rejected, cand = [], {}
        for k in sorted(groups):
            g = groups[k]
            rec = g[0][1]
            if len({raw for _, _, raw in g}) > 1:
                why = 'conflicting copies of the audit row'
            elif self.freq and k % self.freq:
                why = 'step not a multiple of the audit frequency'
            elif k < 1:
                why = 'step below 1'
            else:
                try:
                    t = frac(rec['t_s'])
                except ValueError:
                    t = None
                why = ''
                if t is None:
                    why = 'time not a number'
                elif self.exact and abs(t - (self.t0 + k * self.dt)) > half_ulp(rec['t_s']):
                    why = 'time %s does not match the step clock %.6f s' % (rec['t_s'], float(self.t0 + k * self.dt))
            if not why:
                r = row_of(rec)
                r['tk'] = self.t0 + k * self.dt if self.exact else t
                why = row_defect(r, lz)
            if why:
                self.rejected.append((k, why))
            else:
                cand[k] = r
        rows, prev = {}, None
        keys = sorted(cand)
        for j, k in enumerate(keys):
            r = cand[k]
            if not self.exact and prev is not None and not r['tk'] > rows[prev]['tk']:
                self.rejected.append((k, 'time not increasing with the step'))
                continue
            ref = rows[prev] if prev is not None else (cand[keys[j + 1]] if j + 1 < len(keys) else None)
            if ref is not None and (r['n'] == ref['n']) != (r['wz'] == ref['wz']):
                self.rejected.append((k, 'particle count and buoyant weight change inconsistently'))
                continue
            rows[k] = r
            prev = k
        self.rejected.sort()
        self.rows = rows
        self.row_t = {k: r['tk'] for k, r in rows.items()}
        for r in rows.values():
            r['tkf'] = float(r['tk'])
        self.sorted = [rows[k] for k in sorted(rows)]
        self.K = max([self.K_time] + list(rows))
        self.t_end = self.tk(self.K) if self.K > 0 and (self.exact or self.K in rows) else None
        self.t_last = float(self.t_end) if self.t_end is not None else NAN
        self.rej = dict(self.rejected)

    def coverage(self, a, b, min_samples):
        """(ok, reason): every expected audit row with its step time in [a, b] validated (exact, arithmetic)"""
        if not self.freq:
            return False, 'audit frequency unknown'
        if not self.exact:
            return False, 'no exact clock: the expected audit steps are unknown'
        if a < self.t0:
            return False, 'the window %.4f-%.4f s starts before the clock start %.4f s' % (float(a), float(b), float(self.t0))
        f = self.freq
        k0, k1 = max(1, math.ceil((a - self.t0) / self.dt)), math.floor((b - self.t0) / self.dt)
        j0, j1 = -(-k0 // f), k1 // f
        if j1 - j0 + 1 < max(1, min_samples):
            return False, 'audit too sparse for %.4f-%.4f s (%d expected rows)' % (float(a), float(b), max(0, j1 - j0 + 1))
        if j1 * f > self.K:
            return False, 'the run has not reached %.4f s' % float(self.tk(max(j0, self.K // f + 1) * f))
        miss = [j * f for j in range(j0, j1 + 1) if j * f not in self.rows]
        if miss:
            return False, '%d of %d expected audit rows missing or rejected (first at %.4f s: %s)' % (
                len(miss), j1 - j0 + 1, float(self.tk(miss[0])), self.rej.get(miss[0], 'not in the log'))
        return True, ''

    def window(self, a, b, open_end=False):
        return [r for r in self.sorted if r['tk'] >= a and (r['tk'] < b if open_end else r['tk'] <= b)]


def row_of(rec):
    def g(key, i=0):
        v = rec.get(key)
        return v[i] if v else NAN
    return dict(t=g('t'), step=g('step'), p_in=g('p_in_out'), p_out=g('p_in_out', 1), a_in=g('area_in_out'),
                a_out=g('area_in_out', 1), dp=g('dp_hydro_dyn'), dp_h=g('dp_hydro_dyn', 1), dp_dyn=g('dp_hydro_dyn', 2),
                n=g('n'), wz=g('weight_z_floor_z'), floor=g('weight_z_floor_z', 1), z50=g('z50_95_99_max'),
                z95=g('z50_95_99_max', 1), z99=g('z50_95_99_max', 2), zmax=g('z50_95_99_max', 3), eps=g('eps_bed'),
                zlo=g('window'), zhi=g('window', 1), gw=g('gradz_win'), fbal=g('floor_bal_z'),
                has_fbal='floor_bal_z' in rec, ekern=g('eps_kern_slip_win'), skern=g('eps_kern_slip_win', 1),
                margin=g('win_margin') if 'win_margin' in rec else None, gw_dyn=g('gradz_win_dyn'),
                has_gw_dyn='gradz_win_dyn' in rec)


def row_defect(r, lz):
    """the writer's identities for one audit row (EL_BED_AUDIT), to print precision, and the count domain"""
    if not (isinstance(r['n'], int) and r['n'] >= 0):
        return 'particle count %s outside its domain (non-negative integer)' % r['n']
    p8 = lambda x: es_prec(x, 8)                       # noqa: E731  (ES16.8)
    if finite(r['dp']) != (finite(r['p_in']) and finite(r['p_out'])):
        return 'pressure identity dp = p_in - p_out violated (availability)'
    if finite(r['dp_dyn']) != (finite(r['dp']) and finite(r['dp_h'])):
        return 'pressure identity dp_dyn = dp - dp_hydro violated (availability)'
    if finite(r['dp']) and abs(r['dp'] - (r['p_in'] - r['p_out'])) > 2 * (p8(r['p_in']) + p8(r['p_out']) + p8(r['dp'])):
        return 'pressure identity dp = p_in - p_out violated (%.8g vs %.8g)' % (r['dp'], r['p_in'] - r['p_out'])
    if finite(r['dp_dyn']) and abs(r['dp_dyn'] - (r['dp'] - r['dp_h'])) > 2 * (p8(r['dp']) + p8(r['dp_h']) + p8(r['dp_dyn'])):
        return 'pressure identity dp_dyn = dp - dp_hydro violated (%.8g vs %.8g)' % (r['dp_dyn'], r['dp'] - r['dp_h'])
    if r['has_gw_dyn'] and finite(lz):
        if finite(r['gw_dyn']) != (finite(r['gw']) and finite(r['dp_h'])):
            return 'gradient identity gradz_win_dyn = gradz_win - dp_hydro / Lz violated (availability)'
        if finite(r['gw_dyn']) and abs(r['gw_dyn'] - (r['gw'] - r['dp_h'] / lz)) > 2 * (p8(r['gw']) + p8(r['gw_dyn']) + p8(r['dp_h']) / lz):
            return 'gradient identity gradz_win_dyn = gradz_win - dp_hydro / Lz violated (%.8g vs %.8g)' % (
                r['gw_dyn'], r['gw'] - r['dp_h'] / lz)
    return ''


# ===================================================================================================
# recognition
# ===================================================================================================
FINISH_REASONS = ('no finish marker', 'the run ended at')     # the reasons of an unfinished spec-11 run


def recognise_inputs(log, deck, deck_ok, spec):
    """(reasons, provenance): the effective physical inputs (last start-up assignment, cross-checked) against
    the spec; the start-up state (fresh start, one invocation, no reload, no dump load)"""
    out, prov = [], []
    startup = log.first_time if log.first_time is not None else float('inf')

    def props(key):
        return [(ln, v) for ln, k, v in log.prop if k == key and ln < startup]

    def effective(key, what):
        a = props(key)
        if not a:
            out.append('%s not echoed (Prop - %s)' % (what, key))
            return None
        if len(a) > 1:
            prov.append('%s assigned %d times at start-up (the last is effective)' % (key, len(a)))
        return a[-1][1]

    def against(label, s, ref):
        if not print_equal(s, ref):
            out.append('%s of the run %s differs from the spec %s' % (label, s, ref))
    if len(log.banners) != 1:
        out.append('no start-up banner in the log' if not log.banners else
                   'the start-up banner appears %d times (concatenated invocations)' % len(log.banners))
    for key, v in sorted(log.el.items()):
        if len(v) > 1:
            out.append("start-up echo 'EL %s' appears %d times (one invocation echoes it once)" % (key, len(v)))
            break
    sp = log.simpar.get('StartingProc')
    if not sp:
        out.append('no start-up mode echo (StartingProc)')
    elif sp[-1][1] != '0':
        out.append('start-up mode StartingProc = %s (a fresh start is StartingProc = 0)' % sp[-1][1])
    if log.dumps:
        out.append('a dump was loaded (%s, line %d)' % (log.dumps[0][1][:60], log.dumps[0][0] + 1))
    if log.reloads:
        out.append('reload during the run: %s (line %d)' % (log.reloads[0][1][:70], log.reloads[0][0] + 1))
    late = [(ln, k) for ln, k, _ in log.prop if ln > startup]
    if late:
        out.append('property echo during the run (Prop - %s, line %d)' % (late[0][1], late[0][0] + 1))
    if log.after_finish:
        out.append('solver output after the finish marker (line %d, %s)' % (log.after_finish[0] + 1, log.after_finish[1]))
    for what, store in (('sphere creation (EL terminal-velocity setup)', log.spheres), ('EL fluid viscosity', log.nu_el),
                        ('Global CFD domain extents', log.extents)):
        if len(store) > 1:
            out.append('start-up echo %s appears %d times (one invocation echoes it once)' % (what, len(store)))
    # particles: the population of the spec and its radius and density, from the sphere-creation echo
    sph = log.spheres[-1][1] if log.spheres else None
    if sph is None:
        out.append('particle population, radius and density not echoed (EL terminal-velocity setup)')
    else:
        if int(sph[0]) != spec.N:
            out.append('particle count of the sphere setup %s differs from the spec %d' % (sph[0], spec.N))
        against('particle radius r', sph[1], spec.decimal['d'] / 2)
        against('particle density rho_p', sph[2], spec.decimal['rhop'])
        prov.append('r %s rho_p %s N %s [sphere setup]' % (sph[1], sph[2], sph[0]))
    # fluid
    rho = effective('Density', 'fluid density')
    g = effective('Gravity', 'gravity')
    mu = effective('DynVisc', 'dynamic viscosity')
    nu = effective('Calculated Kinematic Viscosity', 'kinematic viscosity conversion')
    if rho:
        against('fluid density rho_f', rho[0], spec.decimal['rhof'])
    if g and len(g) != 3:
        out.append('gravity echo has %d component(s), three expected' % len(g))
    elif g:
        for i, (gs, gr) in enumerate(zip(g, spec.gravity)):
            against('gravity component %s' % 'xyz'[i], gs, gr)
    if mu:
        against('dynamic viscosity mu', mu[0], spec.decimal['mu'])
    nu_spec = spec.decimal['mu'] / spec.decimal['rhof']
    if nu:
        against('kinematic viscosity (conversion echo)', nu[0], nu_spec)
    if rho and mu and nu:
        try:
            if not print_equal(nu[0], dec(mu[0]) / dec(rho[0])):
                out.append('viscosity echoes disagree: DynVisc %s / Density %s vs the conversion echo %s' % (mu[0], rho[0], nu[0]))
        except (ValueError, ArithmeticError):
            out.append('viscosity echoes not readable')
    el_nu = log.nu_el[-1][1] if log.nu_el else None
    if el_nu is None:
        out.append('closure viscosity not echoed (EL fluid viscosity: ... mu = ... (nu = ...))')
    else:
        against('operator viscosity mu (EL echo)', el_nu[0], spec.decimal['mu'])
        against('closure viscosity nu (EL echo)', el_nu[1], nu_spec)
    if rho and mu:
        prov.append('rho_f %s mu %s nu %s [Prop, last start-up assignment; conversion %s, EL echo %s]'
                    % (rho[0], mu[0], el_nu[1] if el_nu else '?', nu[0] if nu else '?', el_nu[0] if el_nu else '?'))
    if g:
        prov.append('g_z %s [Prop]' % g[-1])
    # model
    dm = log.el_last('drag model')
    if dm is None:
        dm = deck_str(deck, 'ELDragModel') if deck_ok else None
        src = 'deck'
    else:
        src = 'echo'
    if dm is None:
        out.append('drag model not echoed and not in a readable deck')
    elif dm.lower() != spec.drag:
        out.append('drag model of the run %s [%s] differs from the spec %s' % (dm, src, spec.drag))
    w = log.el_last('kernel width factor')
    if w is not None:
        wsrc = 'echo'
    elif not deck_ok:
        wsrc = None
        out.append('kernel width factor not echoed and the deck is not readable')
    elif deck_str(deck, 'ELKernelWidthFactor') is not None:
        w, wsrc = deck_str(deck, 'ELKernelWidthFactor'), 'deck'
    else:
        w, wsrc = '2.5', 'EL_CONFIG default (key absent from echo and deck)'
    if wsrc:
        try:
            ok = dec(w) == spec.decimal['width'] if wsrc != 'echo' else print_equal(w, spec.decimal['width'])
        except ValueError:
            ok = False
        if not ok:
            out.append('kernel width factor of the run %s [%s] differs from the spec %s' % (w, wsrc, spec.decimal['width']))
        prov.append('drag %s [%s] kernel width %s [%s]' % (dm, src, w, wsrc))
    ext = log.extents[-1][1][0].split() if log.extents else None
    if not ext or len(ext) != 6:
        out.append('column extents not in the log (Global CFD domain extents)')
    else:
        if not all(print_equal(s, r) for s, r in zip(ext, spec.extents)):
            out.append('column extents of the run %s differ from the spec %s' % (' '.join(ext), ' '.join(map(str, spec.extents))))
        prov.append('column %s [echo]' % ' '.join(ext))
    return out, prov


def recognise_table(log, spec):
    out = []
    declared = log.el_last('inlet velocity table')
    m = re.search(r'\((\d+) rows', declared or '')
    if not log.table:
        out.append('no inlet table echo (EL_INLET_TABLE)')
        return out
    if not m:
        out.append('table echo without its row count (EL inlet velocity table = ... rows)')
    elif int(m.group(1)) != len(log.table):
        out.append('table echo incomplete: %d of %s rows readable' % (len(log.table), m.group(1)))
    if len(log.table) != len(spec.table):
        out.append('echoed table has %d rows, the spec table %d' % (len(log.table), len(spec.table)))
    for i, ((ts, us), (t, u)) in enumerate(zip(log.table, spec.table)):
        if not (print_equal(ts, t) and print_equal(us, u)):
            out.append('echoed table row %d (t %s, U %s) differs from the spec row (t %s, U %s)' % (i + 1, ts, us, t, u))
            break
    return out


# ===================================================================================================
# the evaluation (shared with fl11_movie.py)
# ===================================================================================================
def evaluate(rundir=None, lines=None, deck_text_=None, last_complete=True, name=None, spec=None):
    spec = spec or SPEC11
    if lines is None:
        with open(os.path.join(rundir, 'run_slurm.log'), errors='replace', newline='') as fh:
            text = fh.read()
        last_complete = text.endswith('\n')
        lines = text.split('\n')
        if lines and lines[-1] == '':
            lines.pop()
    deck_ok = deck_text_ is not None
    if deck_text_ is None and rundir:
        p = os.path.join(rundir, '_data', 'q2p1_param.dat')
        try:
            deck_text_ = open(p, errors='replace').read()
            deck_ok = True
        except OSError:
            deck_text_ = ''
    log = Log(lines, last_complete)
    deck = read_deck(deck_text_ or '')
    ext = log.extents[-1][1][0].split() if log.extents else []
    try:
        zmin, zmax = (num(ext[4]), num(ext[5])) if len(ext) == 6 else (NAN, NAN)
        ext_prec = (es_prec(zmin), es_prec(zmax))
    except ValueError:
        zmin = zmax = NAN
        ext_prec = (0.0, 0.0)
    lz = zmax - zmin if finite(zmax - zmin) and zmax > zmin else NAN
    clock = Clock(log, deck, deck_ok, spec, lz)
    in_reasons, provenance = recognise_inputs(log, deck, deck_ok, spec)
    ev = dict(name=name or os.path.basename(os.path.abspath(rundir or '.')), log=log, deck=deck, clock=clock, spec=spec,
              zmin=zmin, zmax=zmax, lz=lz, ext_prec=ext_prec, provenance=provenance)
    r = spec.d / 2
    for row in clock.sorted:
        fb = log.budget.get(row['step'])
        row['dp_force'] = div(fb['drag'][2] + fb['pressure'][2], row['a_in']) if fb else NAN
        row['H'] = row['z99'] + r - zmin
        if not row['has_fbal']:
            pr = log.pairs.get(row['step'])
            dtf = float(clock.dt) if clock.exact else NAN
            row['fbal'] = (pr['dp'][2] - pr['expected'][2]) / dtf if pr and finite(pos(dtf)) else NAN
    tab = recognise_table(log, spec)
    echo_ok = bool(log.table) and not any(x.startswith(('table echo', 'no inlet')) for x in tab)
    ev['echo_table'] = [(num(t), num(u)) for t, u in log.table] if echo_ok else []
    knots = [(Fraction(dec(t)), num(u)) for t, u in log.table] if echo_ok else []
    for row in clock.sorted:
        row['U'] = u_of_t(knots, row['tk']) if knots else NAN
    pre = tab + clock.segment + in_reasons
    if not clock.freq:
        pre.append('audit frequency unknown (EL bed audit freq)')
    if not clock.finished:
        pre.append('no finish marker: the run has not finished')
    elif clock.t_end is None or clock.t_end < spec.t_last:
        pre.append('the run ended at %.4f s, before the last table row at %.4f s' % (clock.t_last, float(spec.t_last)))
    # a positive contradiction of the step clock disqualifies the run (an absent row is a coverage hole)
    off = [(k, w) for k, w in clock.rejected if 'does not match the step clock' in w or w.startswith('step not a multiple')]
    if off:
        pre.append('audit row at step %d contradicts the step clock (%s)' % off[0])
    sph = log.spheres[-1][1] if log.spheres else None
    ev['N0'] = int(sph[0]) if sph else None
    ev['precondition'] = pre
    ev['official'] = not pre
    ev['witness'] = all(p.startswith(FINISH_REASONS) for p in pre)
    ev['finished'] = clock.finished
    # the finish of ONE invocation (the movie's render predicate): a finish marker, nothing read after it, one banner
    ev['finished_single'] = clock.finished and log.after_finish is None and len(log.banners) <= 1
    ev['table_ok'] = not tab
    ev['t_end'] = clock.t_end
    ev['solver_end'] = clock.t_last
    if not tab:
        ev['horizon'] = float(spec.t_last)
    elif echo_ok:
        ev['horizon'] = max(clock.t_last, ev['echo_table'][-1][0]) if finite(clock.t_last) else ev['echo_table'][-1][0]
    else:
        ev['horizon'] = clock.t_last
    ev['steps'] = [eval_step(ev, s) for s in spec.steps]
    branches(ev)
    last = clock.sorted[-1] if clock.sorted else None
    ev['area'] = pos(last['a_in']) if last else NAN
    ev['W'] = pos(-last['wz']) if last else NAN
    return ev


def eval_step(ev, s):
    clock, spec, h = ev['clock'], ev['spec'], s['hold']
    x = dict(s)
    x.update(status='measured', why='')

    def fail(why):
        if x['status'] == 'measured':
            x['status'], x['why'] = 'unmeasurable', why
    t0, t1 = x['t0'], x['t1']
    a = t1 - spec.EVAL_LEN
    if a < h['a']:
        fail('hold %.4f-%.4f s shorter than the evaluation window' % (float(h['a']), float(h['b'])))
    ok, why = clock.coverage(a, t1, spec.MIN_SAMPLES)
    if not ok:
        fail(why)
    rows = clock.window(a, t1)
    x['rows'] = rows
    if rows and not all(finite(r['dp_dyn']) for r in rows):
        fail('non-finite dp_dyn in the window')
    W = {-r['wz'] for r in rows}
    A = {r['a_in'] for r in rows}
    N = {r['n'] for r in rows}
    if len(W) > 1 or len(A) > 1 or len(N) > 1:
        fail('buoyant weight, area or particle count not constant in the window')
    Wk = pos(next(iter(W))) if len(W) == 1 else NAN
    Ak = pos(next(iter(A))) if len(A) == 1 else NAN
    if rows and not (finite(Wk) and finite(Ak)):
        fail('no positive finite buoyant weight and inlet area in the window')
    wa = div(Wk, Ak)
    x['W'], x['A'], x['n'], x['wa'] = Wk, Ak, (next(iter(N)) if len(N) == 1 else NAN), wa

    def st(key, lo, hi, open_end=False):
        v = [r[key] for r in clock.window(lo, hi, open_end)]
        if not v:
            return NAN, NAN
        m = sum(v) / len(v)
        return m, math.sqrt(sum((y - m) ** 2 for y in v) / len(v))

    def sub(lo, hi, open_end):
        if lo < h['a'] or hi > h['b'] or not clock.coverage(lo, hi, spec.MIN_SUB)[0]:
            return NAN
        return st('dp_dyn', lo, hi, open_end)[0]
    x['dp'], x['std'] = st('dp_dyn', a, t1)
    s0, s1, s2 = (t0 + q for q in spec.SUB)
    m1, m2 = sub(s0, s1, True), sub(s1, s2, False)
    x['mf'], _ = st('dp_force', a, t1)
    x['H'], _ = st('H', a, t1)
    x['eps'], _ = st('eps', a, t1)
    fl, _ = st('fbal', a, t1)
    flpe, _ = st('floor', a, t1)
    x.update(dps=div(x['dp'], wa), stds=div(x['std'], wa), m1s=div(m1, wa), m2s=div(m2, wa), mfs=div(x['mf'], wa),
             fl=div(fl, Wk), flpe=div(flpe, Wk),
             flag=finite(wa) and finite(m1) and finite(m2) and abs(m1 - m2) > spec.FLAG * wa)
    x['measured'] = x['status'] == 'measured'
    return x


def branch_sets(steps, spec):
    """the fluidised set (F1) and the fixed set (F2, F3) of the MEASURED steps, by the rule of the spec amendment
    of 2026-11-24 (F1 = the steps with U above the run's own measured U_mf):
      start    S_1 = measured steps with dp* > FLU; P_r = mean dp* over S_r;
      round r  X_r = measured steps not in S_r with dp* < FIX x P_r (the fixed set); the F3 fit dp* = a U + b U^2
               over X_r (fixed_fit, the code path of F3); U_mf,r = its crossing with dp* = 1 (smallest positive
               root); S_(r+1) = measured steps with U > U_mf,r (strictly: a step at U_mf is not above);
      settled  S_(r+1) = S_r (X and U_mf follow from S, so both sets are unchanged): the sets are S_r, X_r, U_mf,r,
               P_r; every other measured step is intermediate.
    Returns dict(flu, fixed, plateau, umf, rounds, start, reason). reason '' means settled, consistent sets; any
    other reason means the sets are NOT determined (no F1-F3 verdict, the caller prints the reason; there is no
    fallback to the threshold rule): no start set; in a round an empty fixed set, a fit that is not identified,
    a fit without a positive crossing, a measured step whose velocity coincides with U_mf within UMF_TIE (the
    log does not resolve on which side it lies, and round-off must not decide it), or no measured step above
    U_mf; a cycle of sets, or no settlement within
    ROUNDS rounds; at the settled state a step above U_mf with dp* < FIX x P (fluidised by U, fixed by dp*: the
    branches contradict each other)."""
    meas = [x for x in steps if x['measured']]
    start = [x for x in meas if x['dps'] > spec.FLU]
    out = dict(flu=[], fixed=[], plateau=NAN, umf=NAN, rounds=0, start=start, reason='')
    if not start:
        out['reason'] = f'no step with dp* > {spec.FLU} (no measured plateau to start the U_mf iteration from)'
        return out
    flu, seen = start, [frozenset(x['k'] for x in start)]
    for r in range(1, spec.ROUNDS + 1):
        out['rounds'] = r
        plateau = exact_mean([x['dps'] for x in flu])
        fixed = [x for x in meas if x['k'] not in seen[-1] and x['dps'] < spec.FIX * plateau]
        if not fixed:
            out['reason'] = (f'round {r} of the U_mf iteration: no fixed-branch step (dp* < {spec.FIX} of the plateau '
                             f'{plateau:.4f}, outside the fluidised set): no fit, no U_mf')
            return out
        fit = fixed_fit(fixed, spec)
        umf = fit['cross'](1.0) if not fit['why'] else NAN
        if fit['why'] or not finite(umf):
            out['reason'] = (f"round {r} of the U_mf iteration: {fit['why']}" if fit['why'] else
                             f"round {r} of the U_mf iteration: the fit dp* = {fit['a']:.4f} U + {fit['b']:.4f} U^2 of "
                             f"{len(fixed)} fixed steps has no positive crossing with dp* = 1 (no U_mf)")
            return out
        tie = [x for x in meas if abs(x['U'] - umf) <= spec.UMF_TIE * umf]
        if tie:
            out['reason'] = (f"round {r} of the U_mf iteration: " + ', '.join(f"step {x['k']} (U {x['U']:.3f})" for x in tie)
                             + f" coincides with U_mf {umf:.6f} m/s within the resolution of the logged pressures (relative "
                             f"{spec.UMF_TIE:g}): whether it lies above U_mf is not determined")
            return out
        new = [x for x in meas if x['U'] > umf]
        if not new:
            out['reason'] = f'round {r} of the U_mf iteration: no measured step above U_mf {umf:.4f} m/s'
            return out
        keys = frozenset(x['k'] for x in new)
        if keys == seen[-1]:
            out.update(flu=flu, fixed=fixed, plateau=plateau, umf=umf)
            bad = [x for x in flu if x['dps'] < spec.FIX * plateau]
            if bad:
                out['reason'] = ('the branches contradict each other: ' + '; '.join(
                    f"step {x['k']} (U {x['U']:.3f}) lies above U_mf {umf:.4f} m/s with dp* {x['dps']:.4f} < {spec.FIX} x "
                    f"the plateau {plateau:.4f}" for x in bad))
            return out
        if keys in seen:
            j = seen.index(keys) + 1
            out['reason'] = (f'the U_mf iteration does not settle: round {r} returns to the fluidised set of round {j} '
                             f'(a cycle of {r + 1 - j} sets)')
            return out
        seen.append(keys)
        flu = new
    out['reason'] = f'the U_mf iteration did not settle within {spec.ROUNDS} rounds'
    return out


def branches(ev):
    """the sets of branch_sets and the branch of every step, once (gates and display): fluidised (U above the
    measured U_mf), fixed, intermediate, unclassified (sets not determined), unmeasurable"""
    S = ev['sets'] = branch_sets(ev['steps'], ev['spec'])
    ev['plateau'] = S['plateau'] if not S['reason'] else NAN
    flu, fixed = {x['k'] for x in S['flu']}, {x['k'] for x in S['fixed']}
    for x in ev['steps']:
        if not x['measured']:
            x['branch'] = 'unmeasurable'
        elif S['reason']:
            x['branch'] = 'unclassified'
        elif x['k'] in flu:
            x['branch'] = 'fluidised'
        elif x['k'] in fixed:
            x['branch'] = 'fixed'
        else:
            x['branch'] = 'intermediate'
    ev['classes'] = [x['branch'] for x in ev['steps']]


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


# ===================================================================================================
# report and verdicts
# ===================================================================================================
def print_step_table(res, title):
    print(f'\n{title}')
    print(' step  U [m/s]  window      dp*      std     sub1    sub2   flag  F/(W/A)  H [mm]  eps_bed  floor/W  dp*+floor  (floor_pe/W)')
    for x in res:
        t1 = float(x['t1'])
        print(f" {x['k']:3d}  {x['U']:6.3f}  {t1 - 0.3:5.2f}-{t1:5.2f}  {x['dps']:7.4f}  {x['stds']:7.4f}  {_f(x['m1s'], '6.3f')}  {_f(x['m2s'], '6.3f')}  "
              f"{'DIFF' if x['flag'] else '    '}  {x['mfs']:7.4f}  {x['H']*1e3:6.2f}  {x['eps']:7.4f}  {x['fl']:7.4f}  {x['dps'] + x['fl']:8.4f}   ({x['flpe']:.4f})"
              f"{'' if x['measured'] else '  (unmeasurable: ' + x['why'] + ')'}")
    print('  dp* = dp_dyn A / W over 0.2-0.5 s of the step (W, A of that step); floor from the particle momentum '
          'balance; floor_pe = the pe wall-contact impulse (unrelaxed PGS impulse, 1/0.9 of the applied one); '
          'n/a = a sub-window outside the hold or without its own coverage')


def word(v, official):
    return v if official else NOT_OFFICIAL


def report(ev, a):
    """prints the evaluation, returns the verdicts {gate: PASS | FAIL | unavailable}"""
    log, clock, spec = ev['log'], ev['clock'], ev['spec']
    rows, pre, official = clock.sorted, ev['precondition'], ev['official']
    for b in log.trailing:
        print(f"note: unfinished trailing {b[1]} record at line {b[0] + 1} ignored ({b[2]})", file=sys.stderr)
    if rows:
        last = rows[-1]
        print(f"{ev['name']}: {len(rows)} audit lines, t {rows[0]['t']:.4f} .. {last['t']:.4f} s, "
              f"N {rows[0]['n']}, A {ev['area']:.6e} m2, buoyant weight {ev['W']:.6e} N, W/A {div(ev['W'], ev['area']):.3f} Pa, "
              f"hydrostatic part in dp {last['dp_h']:.4f} Pa (subtracted: dp_dyn)")
    else:
        print(f"{ev['name']}: no validated EL_BED_AUDIT line")
    if log.malformed:
        ex = '; '.join(f"{t} line {i + 1}: {w}" for i, t, w in log.malformed[:3])
        print(f"  note: {len(log.malformed)} malformed record(s) ignored ({ex}{' ...' if len(log.malformed) > 3 else ''})")
    print(f"clock: dt {float(clock.dt) if clock.exact else NAN:.6g} s, start {float(clock.t0) if clock.exact else NAN:.4f} s "
          f"[deck], {clock.K} steps completed ({'run finished' if clock.finished else 'no finish marker'}), audit every "
          f"{clock.freq or '?'} steps: {len(rows)} validated rows, {len(clock.rejected)} rejected")
    if clock.rejected:
        print('  rejected audit rows: ' + '; '.join(f"step {k}: {w}" for k, w in clock.rejected[:3])
              + (' ...' if len(clock.rejected) > 3 else ''))
    print('run inputs: ' + '; '.join(ev['provenance']))
    if a.series:
        with open(a.series, 'w') as fh:
            fh.write('t,step,U,dp,dp_hydro,dp_dyn,dp_force,floor_z,z50,z95,z99,zmax,eps_bed,n\n')
            for r in rows:
                fh.write(','.join(str(r[k]) for k in ('t', 'step', 'U', 'dp', 'dp_h', 'dp_dyn', 'dp_force', 'floor',
                                                      'z50', 'z95', 'z99', 'zmax', 'eps', 'n')) + '\n')
    print(f"protocol: {spec.name}, {len(spec.table)} table rows, {len(spec.holds)} holds; lowest hold "
          f"{float(spec.lowest['a']):.3f}-{float(spec.lowest['b']):.3f} s at {spec.lowest['uf']:.3f} m/s")
    if pre:
        print('\n=== NO OFFICIAL VERDICT: not recognised as a complete, uninterrupted %s run ===' % spec.name)
        for p in pre:
            print('   - ' + p)
        print('=== everything below is a MEASUREMENT; '
              + ('this run fails only "finished": witnessed failures (F5 count change, F5 height exceedance in the '
                 'lowest hold) are verdicts ===' if ev['witness'] else 'no verdict at all ==='))
    print_step_table(ev['steps'], 'steps of the protocol' + ('' if official else ' (measurements)'))
    V = {g: 'unavailable' for g in GATES}
    V.update(gates_f1_f3(ev, official))
    V['F5 height'] = gate_f5_height(ev, official)
    V['F5 count'] = gate_f5_count(ev, official)
    qual = setup_checks(ev)
    print(f"verdicts [{'set-up met' if not qual else 'set-up: ' + '; '.join(qual)}]: "
          + ' | '.join(f"{g} {V[g]}" + (AT_GATE_NOTE if g == 'F1' and ev['f1_at_gate'] else '') for g in GATES)
          + ('' if official else OUTSIDE))
    return V


def setup_checks(ev):
    """the set-up checks of the spec that the gates do not cover: prints the block, returns the qualifiers of
    the verdict line (empty: every check met)"""
    log, clock, spec = ev['log'], ev['clock'], ev['spec']
    qual = []
    print(f"\nset-up checks ({spec.name}; they qualify the verdicts, they do not withhold them):")
    # contacts: the single start-up EL_CONTACT_PARAMS line, every pair kind, the effective pair values
    if not log.contact_params:
        print('   contacts: not echoed (binary before pe pull request #59): contact friction and restitution of this run '
              'are not known to the tool')
        qual.append('contacts not echoed')
    elif len(log.contact_params) > 1:
        print(f"   contacts: EL_CONTACT_PARAMS echoed {len(log.contact_params)} times (one start-up line expected): not determined")
        qual.append('contacts not determined')
    else:
        parts = [p.strip() for p in log.contact_params[0][1][len('EL_CONTACT_PARAMS'):].split('|')]
        dev, pairs, rest = [], 0, None
        for p in parts:
            kv = dict(x.split('=', 1) for x in p.split() if '=' in x)
            if 'pair' in kv:
                if kv.get('wall') == 'none':
                    dev.append('no particle-wall pair (wall=none)')
                    continue
                pairs += 1
                name = kv['pair'] + (' ' + kv['wall'] if 'wall' in kv else '')
                for key, ref in spec.CONTACT.items():
                    if key not in kv or not print_equal(kv[key], ref):
                        dev.append(f"{name} {key}={kv.get(key, 'absent')}")
            elif 'restitution' in kv:
                rest = kv['restitution']
        if rest != spec.RESTITUTION:
            dev.insert(0, f"restitution={rest or 'absent'}")
        if not pairs:
            dev.append('no contact pair echoed')
        text = ', '.join(f"{k} {v}" for k, v in spec.CONTACT.items()) + f", restitution {spec.RESTITUTION}"
        if dev:
            print(f"   contacts: deviates: {'; '.join(dev)} (spec: {text})")
            only_rest = all(d.startswith('restitution=') for d in dev)
            qual.append(f"contacts inelastic ({dev[0]})" if only_rest and rest == 'off' else f"contacts deviate ({'; '.join(dev)})")
        else:
            print(f"   contacts: met ({pairs} pair kind(s): {text})")
    # overlap: the validated EL_CONTACT_STATS rows (their own cadence on the step clock)
    rows, expected = contact_rows(log, clock)
    if not rows:
        print('   overlap: not available (no validated EL_CONTACT_STATS row)')
        qual.append('overlap not available')
        return qual
    complete = expected is not None and sum(1 for r in rows if r[0] <= clock.t_end) == expected
    lim = spec.OVERLAP_MAX * spec.decimal['d']
    prep_end = spec.steps[0]['t0']
    wins = [(x['t1'] - spec.EVAL_LEN, x['t1']) for x in ev['steps'] if x['measured']]

    def worst(sel):
        sel = list(sel)
        return max(sel, key=lambda r: r[1]) if sel else None
    allw = worst(rows)
    after = worst(r for r in rows if r[0] >= prep_end)
    inwin = worst(r for r in rows if any(a <= r[0] <= b for a, b in wins))
    pct = lambda r: 100 * float(r[1] / spec.decimal['d'])          # noqa: E731
    desc = (f"largest {pct(allw):.3f} % of d at {float(allw[0]):.4f} s over {len(rows)} rows"
            + ('' if complete else f" ({expected if expected is not None else '?'} expected)") + "; after the preparation "
            f"(t >= {float(prep_end):.1f} s) " + (f"{pct(after):.3f} % at {float(after[0]):.4f} s" if after else 'no row')
            + '; inside the evaluation windows of the measured steps '
            + (f"{pct(inwin):.3f} % at {float(inwin[0]):.4f} s" if inwin else 'no row') + f"; limit {100 * spec.OVERLAP_MAX:g} % of d")
    if allw[1] <= lim and not complete:
        state = 'not available (coverage incomplete)'
        qual.append('overlap not available (coverage incomplete)')
    elif allw[1] <= lim:
        state = 'met'
    elif inwin and inwin[1] > lim:
        state = 'exceeded inside an evaluation window'
        qual.append(f"overlap {pct(inwin):.2f} % of d at {float(inwin[0]):.3f} s (evaluation window)")
    elif after and after[1] > lim:
        state = 'exceeded after the preparation, outside the evaluation windows'
        qual.append(f"overlap {pct(after):.2f} % of d at {float(after[0]):.3f} s (outside the evaluation windows)")
    else:
        state = 'exceeded during the preparation only'
        qual.append(f"overlap {pct(allw):.2f} % of d at {float(allw[0]):.3f} s (preparation)")
    if allw[1] > lim and not complete:
        state += ' (observed; coverage incomplete)'
        qual[-1] = qual[-1][:-1] + ', coverage incomplete)'
    print(f"   overlap: {desc}: {state}")
    return qual


def contact_rows(log, clock):
    """validated EL_CONTACT_STATS rows: (step time, max_overlap as a decimal); on the step clock at the stream's own
    cadence (EL momentum audit freq), identical copies collapsed, conflicting copies and non-finite or negative
    overlaps dropped; and the number of rows that cadence expects over the completed steps (None: unknown)"""
    if not clock.exact:
        return [], None
    f = log.el_last('momentum audit freq')
    freq = int(f) if f is not None and re.fullmatch(r'\d+', f) and int(f) > 0 else None
    groups = {}
    for _, rec, raw in log.contact_stats:
        groups.setdefault(rec['t_s'], set()).add(raw)
    out = {}
    for ts, raws in groups.items():
        if len(raws) > 1:
            continue
        try:
            t, ov = frac(ts), dec(rec_tok(next(iter(raws)), 'max_overlap'))
        except ValueError:
            continue
        k = round((t - clock.t0) / clock.dt)
        if k < 1 or abs(t - (clock.t0 + k * clock.dt)) > half_ulp(ts) or (freq and k % freq) or not ov >= 0 \
                or k in out:
            continue
        out[k] = (clock.t0 + k * clock.dt, ov)
    return [out[k] for k in sorted(out)], (clock.K // freq if freq else None)


def rec_tok(raw, key):
    tok = raw.split()
    return tok[tok.index(key + '=') + 1]


def gates_f1_f3(ev, official):
    spec, res = ev['spec'], ev['steps']
    V = dict(F1='unavailable', F2='unavailable', F3='unavailable')
    ev['f1_at_gate'] = False
    bad = [x for x in res if not x['measured']]
    if len(res) < 3:
        print(f'\nF1, F2, F3 unavailable: {len(res)} protocol step(s), three needed')
        return V
    S = ev['sets']
    thr = S['start']
    xcheck = (f"mean over the {len(thr)} steps with dp* > {spec.FLU}: {mean([x['dps'] for x in thr]):.6f} (reported, not gated)"
              if thr else f'no measured step with dp* > {spec.FLU}')
    if bad:
        # the dp* of an unmeasurable step is unknown, so its membership is undetermined in every set: in the start
        # set (dp* > FLU), at or below U_mf in the fixed set (dp* < FIX x plateau: it moves the fit, U_mf and with it
        # the F1 set), above U_mf in the F1 set; each one blocks F1, F2 and F3
        why = '; '.join(f"step {x['k']} (U {x['U']:.3f}): {x['why']}" for x in bad[:3]) + (' ...' if len(bad) > 3 else '')
        sub = (f"mean dp* over the {len(S['flu'])} measured steps above the U_mf {S['umf']:.4f} m/s of the measured "
               f"steps: {S['plateau']:.6f}" if not S['reason'] else f"sets of the measured steps not determined ({S['reason']})")
        print(f"\nF1, F2, F3 unavailable ({len(bad)} protocol step(s) unmeasurable, each could belong to the set of a "
              f"gate: {why})")
        print(f"   measured, not a verdict: {sub}")
        print(f"   cross-check, measured: {xcheck}")
        return V
    if S['reason']:
        print(f"\nF1, F2, F3 unavailable: {S['reason']}")
        if thr:
            print(f"   cross-check, measured: {xcheck}")
        return V
    flu, fixed, plateau = S['flu'], S['fixed'], S['plateau']
    G = f1_gate([x['dps'] for x in flu], spec)                    # its plateau is S['plateau'] (the same exact mean)
    v1, n, dist, s, se = G['verdict'], G['n'], G['dist'], G['s'], G['se']
    if official:
        V['F1'] = v1
    note = ''
    if G['at_gate']:
        note = (f", at the gate ({'inside' if dist >= 0 else 'outside'} by {abs(dist):.4f}; standard error of the {n} step "
                f"means {se:.4f})")
        ev['f1_at_gate'] = official
    print(f"\nF1 plateau: mean dp* over the {n} steps above U_mf {S['umf']:.4f} m/s (U {min(x['U'] for x in flu):.2f} .. "
          f"{max(x['U'] for x in flu):.2f}): {plateau:.6f}   (gate 1 +- 0.02: {word(v1, official)}{note}); "
          f"dp* + floor/W over these {sum(x['dps'] + x['fl'] for x in flu)/len(flu):.4f}")
    if n >= 3:
        print(f"   scatter of the {n} step means: sample standard deviation {s:.4f}, standard error {se:.4f}; "
              f"distance to the nearer limit {dist:+.4f} ({'inside' if dist >= 0 else 'outside'}); at the gate when within "
              f"the standard error: {'yes' if note else 'no'}")
    else:
        print(f"   no 'at the gate' annotation: {n} step(s) in the set, a scatter needs three")
    rest = [x for x in res if x['branch'] == 'intermediate']
    print(f"   sets settled in round {S['rounds']} of the U_mf iteration (start: the {len(thr)} steps with dp* > {spec.FLU})"
          + ('; in neither set: ' + ', '.join(f"U {x['U']:.2f} (dp* {x['dps']:.4f})" for x in rest) if rest else ''))
    print(f"   cross-check: {xcheck}")
    v3 = f3_block(ev, fixed, plateau, official)
    v2 = gate_f2(ev, fixed, official)
    if official:
        V['F3'], V['F2'] = v3, v2
    return V


def f1_gate(dps, spec):
    """gate F1 of the step means dps of the F1 set: dict(plateau, verdict, dist, n, s, se, at_gate). dist is the
    distance of the plateau to the nearer of the two limits F1_LIMITS, >= 0 inside the gate (PASS; a plateau on a
    limit passes with dist 0), < 0 outside (FAIL). s is the sample standard deviation (n - 1) of the step means,
    se = s / sqrt(n) the standard error of the plateau; at_gate (the annotation of the spec amendment of
    2026-11-24, it never changes the verdict): n >= 3 and |dist| <= se. With fewer than three steps s and se are
    NaN and at_gate is False."""
    lo, hi = spec.F1_LIMITS
    n, plateau = len(dps), exact_mean(dps)
    dist = min(plateau - lo, hi - plateau)
    s = math.sqrt(sum((x - plateau) ** 2 for x in dps) / (n - 1)) if n >= 3 else NAN
    se = s / math.sqrt(n) if n >= 3 else NAN
    return dict(plateau=plateau, verdict='PASS' if dist >= 0 else 'FAIL', dist=dist, n=n, s=s, se=se,
                at_gate=n >= 3 and abs(dist) <= se)


def fixed_fit(fixed, P):
    """the F3 fit dp* = a U + b U^2 through the origin over the fixed steps (least squares, normal equations):
    dict(a, b, cross, why); why is the reason when the fit is not identified (fewer than two distinct nonzero
    velocities or a relative determinant not above REL_DET), else ''; cross(level) is the smallest positive U of
    the fit at dp* = level (NaN when there is none). Used by the U_mf iteration (branch_sets) and by F3."""
    su2 = sum(x['U']**2 for x in fixed); su3 = sum(x['U']**3 for x in fixed); su4 = sum(x['U']**4 for x in fixed)
    sy1 = sum(x['U']*x['dps'] for x in fixed); sy2 = sum(x['U']**2*x['dps'] for x in fixed)
    det = su2*su4 - su3*su3
    distinct = sorted({x['U'] for x in fixed if x['U'] != 0})
    rel = det / (su2 * su4) if su2 > 0 and su4 > 0 else 0.0
    if len(distinct) < 2 or not rel > P.REL_DET:
        return dict(a=NAN, b=NAN, cross=lambda level: NAN,
                    why=(f"the fit dp* = a U + b U^2 is not identified (rank-deficient: {len(distinct)} distinct "
                         f"nonzero fixed velocit{'y' if len(distinct) == 1 else 'ies'}, relative determinant {rel:.3g})"))
    ca = (sy1*su4 - sy2*su3)/det
    cb = (su2*sy2 - su3*sy1)/det

    def cross(level):
        # smallest positive root of b U^2 + a U = level (cancellation-free pair of roots)
        if not (finite(ca) and finite(cb)): return NAN
        if cb == 0: return level/ca if ca > 0 else NAN
        disc = ca*ca + 4*cb*level
        if disc < 0: return NAN
        q = -0.5 * (ca + math.copysign(math.sqrt(disc), ca))
        roots = [x for x in ((q / cb) if q != 0 else NAN, (-level / q) if q != 0 else NAN) if finite(x) and x > 0]
        return min(roots) if roots else NAN
    return dict(a=ca, b=cb, cross=cross, why='')


def f3_block(ev, fixed, plateau, official):
    P = ev['spec']
    vp = math.pi * P.d ** 3 / 6
    if not all(finite(x['H']) and finite(pos(x['A'])) and finite(x['n']) for x in fixed):
        print("F3 unavailable: bed height, area or particle count not finite on every fixed step")
        return 'unavailable'
    efix = mean([x['eps'] for x in fixed])
    uc = closure_umf(efix, P.d, P.rhop, P.rhof, P.mu, P.g, P.drag)
    ebed = sum(1 - div(x['n'] * vp, x['A'] * x['H']) for x in fixed) / len(fixed)
    ub = closure_umf(ebed, P.d, P.rhop, P.rhof, P.mu, P.g, P.drag)
    fit = fixed_fit(fixed, P)
    if fit['why']:
        print(f"F3 unavailable: {fit['why']}")
        return 'unavailable'
    ca, cb, cross = fit['a'], fit['b'], fit['cross']
    u1 = cross(1.0); up = cross(plateau)
    r3 = div(u1, ub)
    v = ('PASS' if abs(r3 - 1) <= P.TOL_F3 else 'FAIL') if finite(r3) else 'unavailable'
    print(f"F3 fixed branch {len(fixed)} steps (U {min(x['U'] for x in fixed):.3f} .. {max(x['U'] for x in fixed):.3f}): "
          f"dp* = {ca:.4f} U + {cb:.4f} U^2 (through the origin)")
    if not finite(ub):
        print(f"F3 unavailable: no crossing of the {P.drag} closure with the buoyant weight established at the whole-bed "
              f"porosity {ebed:.4f}")
        return 'unavailable'
    print(f"   U_mf at dp* = 1: {u1:.4f} m/s; {P.drag} U_mf at the whole-bed porosity {ebed:.4f}: {ub:.4f} m/s; "
          f"ratio {r3:.4f}   (gate 0.05: {word(v, official)})")
    print(f"   reported, not gated: U_mf at the measured plateau {plateau:.4f}: {up:.4f} m/s = {div(up, ub):.4f} of the closure; "
          f"at the core porosity {efix:.4f} the closure gives {uc:.4f} m/s: ratios {div(u1, uc):.4f} / {div(up, uc):.4f}")
    print(f"F4 U_mf / 0.182 (NETL measurement): {u1/P.UMF_MEASURED:.4f} (dp* = 1) / {up/P.UMF_MEASURED:.4f} (plateau) (reported, not gated)")
    s, i = linfit([x['U'] for x in fixed], [x['dp'] for x in fixed])
    print(f"   secondary (rule before 2026-11-20): linear fit of the fixed branch, slope {s:.1f} Pa s/m (intercept {i:.2f} Pa)")
    return v


def valid_eps(e):
    return finite(e) and 0 < e < 1


def gate_f2(ev, fixed, official):
    P = ev['spec']
    r, vp, zmin = P.d / 2, math.pi * P.d ** 3 / 6, ev['zmin']
    reasons = []
    rows_all = [row for x in fixed for row in x['rows']]
    if rows_all and all(row['has_gw_dyn'] for row in rows_all):
        gkey, gtext = 'gw_dyn', 'gradz_win_dyn (dynamic window gradient)'
    elif rows_all and all(row['dp_h'] == 0.0 for row in rows_all):
        fg = ev['log'].el_last('fluid gravity')
        gkey = 'gw'
        gtext = ('gradz_win (log before 627caf47 without gradz_win_dyn; hydrostatic part zero: dp_hydro = 0 on all '
                 f"{len(rows_all)} audit lines used" + (f", start-up echo 'EL fluid gravity = {fg}'" if fg in ('T', 'F') else '') + ')')
    else:
        gkey, gtext = 'gw', 'gradz_win (total gradient)'
        reasons.append('no dynamic window gradient (gradz_win_dyn absent) and a nonzero hydrostatic part (dp_hydro)')
    required = max(2 * P.d, P.width * P.d)
    lows, highs, fields, short, undecided = [], [], [], 0, 0
    for row in rows_all:
        lo = row['zlo'] - zmin
        hi = row['z99'] + r - row['zhi']
        fm = row['margin'] if row['margin'] is not None else 2 * P.d
        lows.append(lo); highs.append(hi); fields.append(fm)
        if not all(finite(x) for x in (lo, hi, fm)):
            continue
        for m, tol in ((lo, es_prec(row['zlo']) + ev['ext_prec'][0]),
                       (hi, es_prec(row['z99']) + es_prec(row['zhi'])), (fm, es_prec(fm))):
            if tol > 0.05 * required:
                undecided += 1
            elif m < required - tol:
                short += 1
    if not finite(zmin):
        margin_ok = False
        reasons.append('missing geometry: no inlet-plane elevation in the log (Global CFD domain extents)')
    elif not rows_all or not all(finite(x) for x in lows + highs + fields):
        margin_ok = False
        reasons.append('window bounds not finite on every contributing audit line')
    elif undecided:
        margin_ok = False
        reasons.append('window margin undecidable: the coordinates are printed too coarsely for the margin')
    else:
        margin_ok = not short
        if short:
            reasons.append('window inside the mapping support')
    fin = [x for x in lows + highs + fields if finite(x)]
    print(f" F2: gradient {gtext}")
    print(f"     audit window margin over {len(rows_all)} lines: min {div(min(fin), P.d) if fin else NAN:.3f} d "
          f"(lower end - inlet plane min {_fmin(lows)/P.d:.3f} d, z99 + r - upper end min {_fmin(highs)/P.d:.3f} d, "
          f"{'win_margin field' if all(row['margin'] is not None for row in rows_all) else 'margin 2 d (binary before 627caf47)'} "
          f"min {_fmin(fields)/P.d:.3f} d); required max(2 d, mapping support {P.width:g} d "
          f"[spec]) = {div(required, P.d):.3f} d: {'met' if margin_ok else 'NOT met'}")
    bad_e = [row for row in rows_all if not valid_eps(row['ekern'])]
    bad_g = [row for row in rows_all if not finite(row[gkey])]
    if bad_e:
        reasons.append(f'{len(bad_e)} mapped gas-fraction sample(s) of the contributing audit lines not usable '
                       '(non-finite, outside 0 < eps < 1, or the -1 sentinel)')
    if bad_g:
        reasons.append(f'{len(bad_g)} window-gradient sample(s) not finite')
    print(' F2: step  U      eps_bed  window-grad  closure(eps_bed)  ratio_w | eps_kern ratio_k | dp/H_eq  ratio_e | eps_H   dp/H   closure(eps_H)  ratio_h')
    rks, rws, undefined = [], [], []
    for x in fixed:
        g = mean([row[gkey] for row in x['rows']])
        ek = mean([row['ekern'] for row in x['rows']])
        heq = div(x['n'] * vp, x['A'] * (1 - x['eps']))
        clo_b = closure_grad(x['U'], x['eps'], P.d, P.rhof, P.mu, P.drag)
        eh = 1 - div(x['n'] * vp, x['A'] * x['H'])
        cloh = closure_grad(x['U'], eh, P.d, P.rhof, P.mu, P.drag)
        clok = closure_grad(x['U'], ek, P.d, P.rhof, P.mu, P.drag) if valid_eps(ek) else NAN
        rw, rk = div(g, clo_b), div(g, clok)
        re_, rh = div(div(x['dp'], heq), clo_b), div(div(x['dp'], x['H']), cloh)
        if not (finite(pos(clok)) and finite(rk)):
            undefined.append(f"step {x['k']} (U {x['U']:.3f}): closure {clok:.4g}, ratio {rk:.4g}")
        rks.append(rk); rws.append(rw)
        print(f"     {x['k']:3d}  {x['U']:.3f}  {x['eps']:.4f}  {g:10.1f}  {clo_b:10.1f}  {rw:12.4f} | {ek:.4f} {rk:7.4f} | {div(x['dp'], heq):8.1f} {re_:7.4f} | "
              f"{eh:.4f} {div(x['dp'], x['H']):8.1f} {cloh:10.1f}  {rh:9.4f}")
    if undefined:
        reasons.append('no positive closure and finite ratio_k on every fixed step: ' + '; '.join(undefined))
    worst_k = max(abs(x - 1) for x in rks) if all(finite(x) for x in rks) else NAN
    worst_w = max(abs(x - 1) for x in rws) if all(finite(x) for x in rws) else NAN
    if not reasons:
        v = 'PASS' if worst_k <= P.TOL_F2 else 'FAIL'
        print(f"F2 worst |ratio_k - 1| (window gradient against the closure at the mapped gas fraction) {worst_k:.4f}   "
              f"(gate 0.05: {word(v, official)})")
    else:
        v = 'unavailable'
        print(f"F2 measured, not gated: worst |ratio_k - 1| {worst_k:.4f}")
        print(f"F2 verdict unavailable ({'; '.join(reasons)})")
    if finite(worst_w):
        print(f"   reported, not gated: worst |ratio_w - 1| at the geometric window porosity (rule before 2026-11-20) {worst_w:.4f}")
    return v


def gate_f5_count(ev, official):
    clock = ev['clock']
    ns = sorted({r['n'] for r in clock.sorted} | ({ev['N0']} if ev['N0'] is not None else set()))
    if len(ns) > 1:
        if official or ev['witness']:
            print(f"F5 count: particle counts {ns} (start-up population and the validated audit rows, witnessed): FAIL")
            return 'FAIL'
        print(f"F5 count: particle counts {ns} (start-up population and the validated audit rows; measurement, not a verdict)")
        return 'unavailable'
    if not official:
        print(f"F5 count: counts seen {ns}; no official verdict for this run: unavailable")
        return 'unavailable'
    ok, why = clock.coverage(clock.t0, clock.t_end, 1)
    if ok:
        print(f"F5 count: counts seen {ns} on {len(clock.sorted)} audit lines, run {float(clock.t0):.4f}-{clock.t_last:.4f} s "
              f"covered: PASS")
        return 'PASS'
    print(f"F5 count: counts seen {ns}; count over the run not established ({why}): unavailable")
    return 'unavailable'


def gate_f5_height(ev, official):
    """the spec's lowest hold: an exceedance in a validated row is a FAIL (official, or a run that fails only
    'finished'); a PASS needs the official precondition, full coverage and every reading at or below the mid-
    column; a reading above it by less than the print uncertainty is undecidable"""
    clock, spec = ev['clock'], ev['spec']
    h = spec.lowest
    span = f"{float(h['a']):.3f}-{float(h['b']):.3f} s"
    if not finite(ev['lz']):
        print("F5 height unavailable: no column extents in the log (Global CFD domain extents): mid-column unknown")
        return 'unavailable'
    rows = clock.window(h['a'], h['b'])
    ok, why = clock.coverage(h['a'], h['b'], spec.MIN_SAMPLES)
    reasons = [] if ok else [why]
    zmin, lz = ev['zmin'], ev['lz']
    mid = zmin + lz / 2
    pmid = 0.5 * (ev['ext_prec'][0] + ev['ext_prec'][1])
    zs = [r['zmax'] for r in rows]
    fin = [z for z in zs if finite(z)]
    if len(fin) < len(zs):
        reasons.append(f'{len(zs) - len(fin)} audit line(s) without a centre height')
    above = [z for z in fin if z - mid > pmid + es_prec(z)]
    band = [z for z in fin if 0 < z - mid <= pmid + es_prec(z)]
    zm = max(fin) if fin else NAN
    if fin:
        print(f"F5 (b) highest centre during the plateau of the lowest velocity {h['uf']:.3f} m/s ({span}, {len(rows)} audit lines): "
              f"{(zm - zmin)*1e3:.2f} mm above the inlet plane (mid-column {lz*500:.1f} mm)")
    if above and (official or ev['witness']):
        print("F5 height: a centre above mid-column observed in a validated row (witnessed): FAIL")
        v = 'FAIL'
    elif above:
        print(f"F5 height: {len(above)} reading(s) above mid-column (measurement, not a verdict)")
        v = 'unavailable'
    elif band:
        print(f"F5 height unavailable (not decidable at print precision: {len(band)} reading(s) above mid-column by less "
              f"than the print uncertainty)")
        v = 'unavailable'
    elif not official:
        print(f"F5 height: no centre above mid-column in {len(fin)} validated row(s) of the lowest hold ({span}); "
              "no official verdict for this run: unavailable")
        v = 'unavailable'
    elif reasons or not fin:
        print(f"F5 height unavailable ({'; '.join(reasons) if reasons else 'no observation'})")
        v = 'unavailable'
    else:
        print("F5 height: no centre above mid-column (a centre at mid-height is not above): PASS")
        v = 'PASS'
    prev = [g for g in spec.holds if g['b'] < h['a']]
    if fin and prev:
        ta, tb = prev[-1]['b'], h['a']
        zt = _fmax([r['zmax'] for r in clock.window(ta, tb)])
        print(f"   stricter, reported: including the {1e3*float(tb - ta):.0f} ms transition into it "
              f"({float(ta):.3f}-{float(tb):.3f} s): highest centre {(_fmax([zm, zt]) - zmin)*1e3:.2f} mm")
    return v


# ===================================================================================================
# command line
# ===================================================================================================
def parser():
    ap = argparse.ArgumentParser(description='spec-11 gates of a q2p1_el_pipeflow run')
    ap.add_argument('rundir', nargs='?')
    ap.add_argument('--series', default=None)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--self-test-real', default=None, metavar='RUNDIR')
    ap.add_argument('--variants', type=int, default=300)
    ap.add_argument('--seed', type=int, default=20261122)
    return ap


def main(argv=None, spec=None):
    ap = parser()
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if a.self_test_real:
        return self_test_real(a.self_test_real, a.variants, a.seed)
    if not a.rundir:
        ap.error('RUNDIR is required')
    report(evaluate(a.rundir, spec=spec), a)
    return 0


def run_lines(lines, deck, last_complete=True, name='run', spec=None):
    """the report of an in-memory log (the metamorphic tests): (stdout, verdicts, evaluation)"""
    a = parser().parse_args(['x'])
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        ev = evaluate(None, lines=lines, deck_text_=deck, last_complete=last_complete, name=name, spec=spec)
        V = report(ev, a)
    return out.getvalue(), V, ev


# ===================================================================================================
# self-test: the machinery on synthetic logs, each with its fixture spec
# ===================================================================================================
def _e(x, w=14, p=6):
    return ('%*s' % (w, 'NaN')) if x != x else '%*.*E' % (w, p, x)


def _d(x):
    return ('%12.4E' % x).replace('E', 'D')


def _fe(x):
    """Fortran E16.8: 0.dddddddd E+ee"""
    if x == 0:
        return '  0.00000000E+00'
    e = math.floor(math.log10(abs(x))) + 1
    m = round(x / 10 ** e, 8)
    if abs(m) >= 1:
        m, e = m / 10, e + 1
    return '%16s' % ('%.8fE%+03d' % (m, e))


PROTO = [(0.5, 0.30), (1.0, 0.20), (1.5, 0.10), (2.0, 0.05)]          # (start of step, U); end 2.5 s
DT = 2.5e-4
DECK = 'SimPar@TimeStep = 0.00025d0\nSimPar@StartSimTime = 0d0\n'


def _dps(u):
    return min(1.0, u / 0.18) * ((1.0 - 0.1 * (u - 0.2)) if u > 0.18 else 1.0)


def _table(protocol, tstep=0.5):
    tab = [(0.0, protocol[0][1])]
    for (t, u) in protocol[1:]:
        tab += [(t, tab[-1][1]), (t + 0.02, u)]
    t_proto = protocol[-1][0] + tstep
    tab.append((t_proto, tab[-1][1]))
    return tab, t_proto


def _fixture_spec(kw, over=None):
    """the spec of a fixture: the table and physical values the writer uses (over: replacements)"""
    tab = kw.get('table') or _table(kw.get('protocol', PROTO))[0]
    z0 = kw.get('z0', 0.0)
    w = kw.get('width', 2.5)
    s = dict(table=[(repr(float(t)), repr(float(u))) for t, u in tab], d=repr(kw.get('d', 0.87e-3)),
             rhop=repr(kw.get('rhop', 860.0)), rhof=repr(kw.get('rhof', 1.0)), mu=repr(kw.get('mu', 1.8e-5)), g='9.81',
             drag=kw.get('drag', 'gidaspow') or 'gidaspow', width=repr(w) if isinstance(w, float) else '2.5',
             extents=('0', '0.015', '0', '0.015', repr(z0), repr(z0 + 0.1)), N=9789)
    s.update(over or {})
    return Spec('fixture spec', **s)


def _lines(protocol=PROTO, cad=0.01, t_end=None, area=2.25e-4, n=9789, width=2.5, margin_d=2.5, zlo_d=None,
           dyn=True, dp_hydro=0.0, d=0.87e-3, ekern=0.37, ekern_fn=None, extents=True, z0=0.0, table=None,
           echo_table=True, zmax_fn=None, legacy=False, finished=True, time_lines=True, drag='gidaspow',
           width_echo=True, tmax=None, n_fn=None, grad_fn=None, mu=1.8e-5, rhof=1.0, rhop=860.0, clock_fn=None,
           dps_fn=None, densities=None, starting_proc='0', contact='met', overlap_fn=None, contact_stats=True,
           gravity=(0.0, 0.0, -9.81), sph_r=None):
    """the lines of a synthetic run: start-up echo (EL configuration, banner, SimPar, sphere setup, Prop,
    viscosity), table echo, a time line per step, EL_FORCE_BUDGET and EL_BED_AUDIT lines, the finish marker"""
    g = 9.81
    vp = math.pi * d ** 3 / 6
    tab, t_proto = _table(protocol)
    if table is not None:
        tab = table
    t_end = t_proto if t_end is None else t_end
    freq = int(round(cad / DT))
    nsteps = int(round(t_end / DT))
    L = ['EL kernel                 = deen_poly']
    if width_echo:
        L.append('EL kernel width factor    = ' + (_e(width) if not isinstance(width, str) else width))
    if drag:
        L.append('EL drag model             = ' + drag)
    L.append('EL fluid gravity          = ' + ('T' if dp_hydro else 'F'))
    if echo_table:
        L.append('EL inlet velocity table   = _data/inlet_table.txt (%d rows t / U; linear, constant beyond the ends)' % len(tab))
        L += ['  EL_INLET_TABLE %s %s' % (_e(t, 16, 8), _e(u, 16, 8)) for t, u in tab]
    L.append('EL bed audit freq         = %d' % freq)
    L.append('EL momentum audit freq    = %d' % freq)
    L += ['  |          %s                                                  |' % BANNER,
          'StartingProc = %s' % starting_proc, 'TimeStep =   %s' % _d(DT).strip(), 'StartSimTime =   %s' % _d(0.0).strip()]
    if extents:
        L.append('Global CFD domain extents: ' + ' '.join(_e(x) for x in (0.0, 0.015, 0.0, 0.015, z0, z0 + 0.1)))
    if contact:
        pair = 'mu_static=0.3 mu_dynamic=0.3 e=0.9'
        L.append({'met': 'EL_CONTACT_PARAMS pair=particle-particle %s | pair=particle-wall wall=zwalls %s | restitution=on '
                         'threshold=0.01 scheme=post-pass' % (pair, pair),
                  'off': 'EL_CONTACT_PARAMS pair=particle-particle %s | pair=particle-wall wall=zwalls %s | restitution=off '
                         'threshold=0 scheme=post-pass' % (pair, pair)}.get(contact, contact))
    L.append('EL terminal-velocity setup: %d free sphere(s) created (r=%s, rho=%g, seedMode=file, seedDomain=box), '
             'PE decomposition 3x3x4.' % (n, sph_r or repr(d / 2), rhop))
    L.append('Prop    - Gravity =' + ''.join(_fe(x) for x in gravity))
    for rho in (densities or [rhof]):
        L.append('Prop    - Density =' + _fe(rho) + _fe(rho))
    L.append('Prop    - DynVisc =' + _fe(mu) + _fe(1.0))
    L.append('Prop    - Calculated Kinematic Viscosity =' + _fe(mu / (densities or [rhof])[-1]) + _fe(1.0))
    L.append(' EL fluid viscosity: diffusion operator rebuilt with mu = %s (nu = %s for the closures)'
             % (_e(mu, 16, 8), _e(mu / rhof, 16, 8)))
    tmx = t_end if tmax is None else tmax
    W0 = n * (rhop - rhof) * vp * g
    for k in range(1, nsteps + 1):
        t, dtk = clock_fn(k) if clock_fn else (k * DT, DT)
        if k % freq == 0:
            u = u_of_t([(Fraction(repr(float(a))), b) for a, b in tab], Fraction(repr(t)))
            dps = dps_fn(u) if dps_fn else _dps(u)
            nn = n_fn(t) if n_fn else n
            W = nn * (rhop - rhof) * vp * g
            dp = dps * div(W, area) + dp_hydro if (area > 0 and nn > 0) else NAN
            z99 = z0 + 0.024 if nn > 0 else NAN
            m = margin_d * d
            zlo = z0 + (zlo_d if zlo_d is not None else margin_d) * d if nn > 0 else NAN
            zhi = z99 + d / 2 - m if nn > 0 else NAN
            zmx = (z0 + (zmax_fn(t) if zmax_fn else 0.0241)) if nn > 0 else NAN
            grad = (grad_fn(t, dp, u) if grad_fn else dp / 0.1) if dp == dp else NAN
            gdyn = grad - dp_hydro / 0.1 if dp == dp else NAN
            ek = (ekern_fn(t) if ekern_fn else ekern) if nn > 0 else NAN
            line = ('EL_BED_AUDIT t= %s step= %d p_in_out= %s %s area_in_out= %s %s dp_hydro_dyn= %s %s %s n= %d '
                    'weight_z_floor_z= %s %s z50_95_99_max= %s %s %s %s eps_bed= %s window= %s %s gradz_win= %s' % (
                        _e(t), k, _e(dp, 16, 8), _e(0.0, 16, 8), _e(area), _e(area), _e(dp, 16, 8),
                        _e(dp_hydro, 16, 8), _e(dp - dp_hydro if dp == dp else NAN, 16, 8), nn,
                        _e(-W if nn > 0 else 0.0), _e(0.0), _e(z99 - 0.012), _e(z99 - 0.001), _e(z99), _e(zmx),
                        _e(0.37 if nn > 0 else NAN), _e(zlo), _e(zhi), _e(grad, 16, 8)))
            if not legacy:
                line += ' floor_bal_z= %s eps_kern_slip_win= %s %s' % (_e((1 - dps) * W, 16, 8), _e(ek, 16, 8),
                                                                      _e(0.1 if nn > 0 else NAN, 16, 8))
                if dyn:
                    line += ' win_margin= %s gradz_win_dyn= %s' % (_e(m if nn > 0 else NAN), _e(gdyn, 16, 8))
            L.append('EL_FORCE_BUDGET step= %d drag= %s %s %s pressure= %s %s %s lift= 0 0 0 grav_buoy= %s %s %s' % (
                k, _e(0.0), _e(0.0), _e(dps * W0 * 0.9), _e(0.0), _e(0.0), _e(dps * W0 * 0.1), _e(0.0), _e(0.0), _e(-W0)))
            L.append(line)
            if contact_stats:
                L.append('EL_CONTACT_STATS t= %s ncontacts= 4 max_overlap= %s Tgran= %s' % (
                    _e(t), _e(overlap_fn(t) if overlap_fn else 3.8e-7), _e(1e-6)))
        if time_lines:
            L.append('time:%s/%s | itns:%8d/%8d | dt:%s%s%s' % (_d(t), _d(tmx), k, nsteps, _d(dtk), _d(dtk), _d(dtk)))
    if finished:
        L.append(' PP3D_LES has successfully finished. ')
    return L


def _write_run(path, deck=DECK, mutate=None, last_complete=True, **kw):
    os.makedirs(os.path.join(path, '_data'), exist_ok=True)
    L = _lines(**kw)
    if mutate:
        L = mutate(L)
    with open(os.path.join(path, 'run_slurm.log'), 'w') as fh:
        fh.write('\n'.join(L) + ('\n' if last_complete else ''))
    if deck is not None:
        with open(os.path.join(path, '_data', 'q2p1_param.dat'), 'w') as fh:
            fh.write(deck)
    return path


def _audit_at(lines, t):
    for i, l in enumerate(lines):
        if l.startswith('EL_BED_AUDIT') and abs(float(l.split()[2]) - t) < 1e-9:
            return i
    raise KeyError(t)


def _drop(*ts):
    return lambda L: [l for i, l in enumerate(L) if i not in {_audit_at(L, t) for t in ts}]


def _set_tok(line, key, idx, text):
    tok = line.split()
    tok[tok.index(key + '=') + 1 + idx] = text
    return ' '.join(tok)


def _set_field(t, key, idx, text):
    def m(L):
        L = list(L)
        i = _audit_at(L, t)
        L[i] = _set_tok(L[i], key, idx, text)
        return L
    return m


def _add_field(t, key, idx, add):
    def m(L):
        L = list(L)
        i = _audit_at(L, t)
        tok = L[i].split()
        j = tok.index(key + '=') + 1 + idx
        tok[j] = _e(num(tok[j]) + add, 16, 8)
        L[i] = ' '.join(tok)
        return L
    return m


def _dup(t, copies):
    def m(L):
        i = _audit_at(L, t)
        return L[:i + 1] + [L[i]] * copies + L[i + 1:]
    return m


def _conflict(t):
    def m(L):
        i = _audit_at(L, t)
        tok = L[i].split()
        j = tok.index('eps_bed=') + 1
        tok[j] = _e(num(tok[j]) * 1.001)
        return L[:i + 1] + [' '.join(tok)] + L[i + 1:]
    return m


def _thin(keep_every):
    def m(L):
        out, j = [], 0
        for l in L:
            if l.startswith('EL_BED_AUDIT'):
                j += 1
                if j % keep_every:
                    continue
            out.append(l)
        return out
    return m


def _append(text):
    return lambda L: L + [text]


def _after_step(k, new):
    """insert lines after the time line of step k"""
    def m(L):
        i = next(j for j, l in enumerate(L) if PAT_TIME.match(l) and int(PAT_TIME.match(l).group(3)) == k)
        return L[:i + 1] + list(new) + L[i + 1:]
    return m


def _reload_extension(t_first_max, until):
    def m(L):
        out = []
        for l in L:
            mm = PAT_TIME.match(l)
            if mm and num(mm.group(1)) <= until + 1e-9:
                l = l.replace('/' + _d(num(mm.group(2))), '/' + _d(t_first_max), 1)
            out.append(l)
        return out
    return m


def _renumber(L, k0):
    """the part of a log after the step k0, steps renumbered from 1 (a restart segment), echo kept"""
    head, out, started = [], [], False
    for l in L:
        mm = PAT_TIME.match(l)
        if not started:
            if mm and int(mm.group(3)) == k0:
                started = True
            elif not (l.startswith(('EL_BED_AUDIT', 'EL_FORCE_BUDGET', 'EL_NEWTON_PAIR')) or mm):
                head.append(l)
            continue
        if mm:
            l = re.sub(r'itns:\s*%d/' % int(mm.group(3)), 'itns:%8d/' % (int(mm.group(3)) - k0), l, count=1)
        elif l.startswith(('EL_BED_AUDIT', 'EL_FORCE_BUDGET', 'EL_NEWTON_PAIR')):
            k = int(re.search(r'step=\s*(\d+)', l).group(1))
            l = re.sub(r'step=\s*%d\b' % k, 'step= %d' % (k - k0), l, count=1)
        out.append(l)
    return head + out


def _time_line(l, t_new, dt_new=None):
    mm = PAT_TIME.match(l)
    dt = num(mm.group(5)) if dt_new is None else dt_new
    return 'time:%s/%s | itns:%8d/%8d | dt:%s%s%s' % (_d(t_new), _d(num(mm.group(2))), int(mm.group(3)), int(mm.group(4)),
                                                     _d(dt), _d(dt), _d(dt))


P14 = [(0.5 + 0.5 * i, round(0.30 - 0.02 * i, 2)) for i in range(14)]           # 0.30 .. 0.04 m/s, the spec's velocities
P5 = [(0.5, 0.30), (1.0, 0.25), (1.5, 0.20), (2.0, 0.10), (2.5, 0.05)]
P5b = [(0.5, 0.30), (1.0, 0.20), (1.5, 0.16), (2.0, 0.10), (2.5, 0.05)]
P7 = [(0.5, 0.30), (1.0, 0.25), (1.5, 0.20), (2.0, 0.18), (2.5, 0.14), (3.0, 0.10), (3.5, 0.05)]
CYCLE = [1.0, 1.0, 1.0, 0.90, 0.79, 0.10 / 0.19, 0.05 / 0.19]   # with 0.14 fixed U_mf 0.1722, without 0.19 (exact)
def _tie(slope):
    """per-step dp* of the 14-step protocol: plateau 0.982 above 0.20 m/s, 0.95 / 0.90 / 0.80 at 0.20 / 0.18 / 0.16,
    fixed branch slope x U from 0.14 m/s down (U_mf = 1 / slope)"""
    return [0.982] * 5 + [0.95, 0.90, 0.80] + [slope * round(0.14 - 0.02 * i, 2) for i in range(6)]


TIE5 = _tie(5.0)                # review of r8: U_mf = 0.20 exactly, the velocity of step 5
Z80 = [0.9853, 0.9729, 0.9765, 0.9805, 0.9784, 0.9831, 0.9892, 0.9582, 0.8274, 0.6996, 0.5748, 0.4534, 0.3354, 0.2209]
R5 = [0.9716, 0.9677, 0.9732, 0.9746, 0.9774, 0.9861, 0.9936, 0.9147, 0.7905, 0.6688, 0.5499, 0.4338, 0.3210, 0.2115]
PEER = {  # per-step dp* of the peer codes, U 0.30 .. 0.04 m/s (spec 11 brief of 2026-11-24)
    'CFDEM 40': [0.9643, 0.9587, 0.9501, 0.9558, 0.9481, 0.9529, 0.9357, 0.8212, 0.7089, 0.5991, 0.4921, 0.3879, 0.2864, 0.1880],
    'CFDEM 80': [0.9752, 0.9786, 0.9805, 0.9696, 0.9675, 0.9819, 0.9804, 0.9506, 0.8193, 0.6922, 0.5709, 0.4478, 0.3307, 0.2169],
    'DPMFoam': [0.9441, 0.9452, 0.9412, 0.9421, 0.9551, 0.9434, 0.9000, 0.8090, 0.7059, 0.6011, 0.4970, 0.3940, 0.2932, 0.1934],
    'MFiX': [0.9888, 0.9991, 0.9928, 0.9883, 1.0053, 0.9970, 0.9926, 0.9173, 0.7904, 0.6669, 0.5469, 0.4304, 0.3177, 0.2081]}


def _by_u(values, protocol=None):
    """dps_fn: the given dp* on the holds of the protocol (P14 by default), _dps(u) in the ramps"""
    table = {round(u, 3): v for (_, u), v in zip(protocol or P14, values)}
    return lambda u: table.get(round(u, 3), _dps(u))


def self_test_cases():
    """(name, writer keywords, expected verdicts, strings that must appear, strings that must not, the answer of
    an earlier tool where it differed). Keyword '_spec': the fixture spec ('base' = the four-step fixture with
    the default physics, a dict = replacements of the derived one; default: derived from the writer keywords).
    Keyword '_set': gate constants of that fixture spec (e.g. ROUNDS). The fixtures test the machinery with a
    fixture spec, not spec 11."""
    hi_after = lambda t: 0.060 if t > 2.3 else 0.0241       # noqa: E731
    REF = {'F1': 'PASS', 'F2': 'FAIL', 'F3': 'PASS', 'F5 height': 'PASS', 'F5 count': 'PASS'}
    U = 'unavailable'
    ALLU = {g: U for g in GATES}
    NO = 'NO OFFICIAL'
    base_tab = _table(PROTO)[0]
    tab_sub = base_tab[:-1] + [(2.1, 0.05), (2.18, 0.05), (2.26, 0.05), (2.34, 0.05), (2.42, 0.05), (2.5, 0.05)]
    tab_knot = sorted(base_tab + [(1.3, 0.20)])
    rep_proto = [(0.5, 0.30), (1.0, 0.10), (1.5, 0.30), (2.0, 0.10), (2.5, 0.0)]
    nofl_proto = [(0.5, 0.14), (1.0, 0.12), (1.5, 0.10), (2.0, 0.08)]
    mu_e = 2.4e-5
    grad_mu = lambda t, dp, u: closure_grad(u, 0.37, 0.87e-3, 1.0, mu_e, 'gidaspow')  # noqa: E731
    sub_tab = [(0.0, 0.3), (1.0, 0.3), (1.02, 0.2), (1.4, 0.2), (1.42, 0.1), (1.9, 0.1), (1.92, 0.05), (2.4, 0.05)]
    reset = lambda k: ((k * DT - 0.5) if k > 7000 else k * DT, DT)                  # noqa: E731
    dtchg = lambda k: ((k * DT) if k <= 9000 else 9000 * DT + (k - 9000) * 2 * DT, DT if k <= 9000 else 2 * DT)  # noqa: E731
    shifted = base_tab[:-1] + [(2.5000002, 0.05)]
    excursion = sorted(base_tab + [(1.2, 0.2), (1.2001, 0.01), (1.2999, 0.01), (1.3, 0.2)])
    tied = _table([(0.5, 0.05), (1.0, 0.20), (1.5, 0.10), (2.0, 0.05)])[0]
    reload_lines = ['Reloading Physical parameters from file: _data/q2p1_param.dat',
                    'Prop    - DynVisc =' + _fe(3.6e-5) + _fe(1.0)]
    two_runs = lambda L: L + _lines(n=9788, cad=0.02)                     # noqa: E731
    overlap = {0.3: 1.80, 0.2: 0.99, 0.1: 0.40, 0.05: 0.20}
    return [
        ('reference (all gates measurable)', dict(), REF, ['stricter, reported', 'run inputs:', 'verdicts [set-up met]:'],
         ['Traceback', NO], ''),
        ('missing inlet plane (area_in = 0)', dict(area=0.0), {'F1': U, 'F2': U, 'F3': U, 'F5 height': 'PASS'}, [], ['Traceback'], ''),
        ('empty bed (n = 0)', dict(n=0), ALLU, ['particle count of the sphere setup 0 differs from the spec 9789'], ['Traceback'], ''),
        ('short run (2.3 s, no finish marker)', dict(t_end=2.3, finished=False), ALLU, ['no finish marker: the run has not finished'], [], ''),
        ('sparse audit (cadence 0.4 s)', dict(cad=0.4), {'F1': U, 'F5 height': U}, ['audit too sparse'], [], ''),
        ('window margin 2 d with a 2.5 d support', dict(margin_d=2.0), {'F2': U, 'F1': 'PASS'}, ['window inside the mapping support'], [], ''),
        ('lower window end at 2 d, field says 2.5 d', dict(zlo_d=2.0), {'F2': U}, ['window inside the mapping support'], [], ''),
        ('old log: no win_margin / gradz_win_dyn, zero hydrostatic part', dict(margin_d=2.0, dyn=False), {'F2': U},
         ['hydrostatic part zero', 'margin 2 d (binary before 627caf47)'], [], ''),
        ('old log with a 2 d support: margin met, total gradient admitted', dict(margin_d=2.0, dyn=False, width=2.0),
         {'F2': 'FAIL'}, ['hydrostatic part zero'], [], ''),
        ('old log with a hydrostatic part', dict(margin_d=2.0, dyn=False, width=2.0, dp_hydro=0.981), {'F2': U}, ['nonzero hydrostatic part'], [], ''),
        ('mapped porosity 1.0 (closure zero)', dict(ekern=1.0), {'F2': U}, ['no positive closure'], [], ''),
        ('the -1 sentinel in the mapped samples', dict(ekern_fn=lambda t: -1.0 if 2.3 < t < 2.4 else 0.37), {'F2': U},
         ['mapped gas-fraction sample(s)'], [], ''),
        ('zero-velocity fixed plateau', dict(protocol=PROTO + [(2.5, 0.0)]), {'F2': U}, ['step 4 (U 0.000)'], [], ''),
        ('r6: last knot removed, 60 mm late, no finish marker', dict(table=base_tab[:-1], zmax_fn=hi_after, finished=False, _spec='base'),
         ALLU, ['echoed table has 7 rows, the spec table 8'], ['(witnessed)'], 'f8fb7fd8 (r6 P1-5): F5 FAIL on the 0.10 m/s hold'),
        ('no time lines and no finish marker', dict(time_lines=False, finished=False), ALLU, ['no time lines in the log'], [], ''),
        ('one audit line of a binary without floor_bal_z', dict(legacy=True, cad=2.5), {'F1': U, 'F5 height': U}, [], ['Traceback'], ''),
        ('no column extents', dict(zlo_d=1.0, extents=False), ALLU, ['column extents not in the log'], [], ''),
        ('30 ms audit (every window qualifies)', dict(cad=0.03), {'F1': 'PASS', 'F2': 'FAIL'}, [], ['unmeasurable'], ''),
        ('only the final rows in the log', dict(mutate=lambda L: [l for l in L if not l.startswith('EL_BED_AUDIT')
                                                                 or float(l.split()[2]) > 2.0]),
         {'F1': U, 'F5 count': U, 'F5 height': 'PASS'}, [], [], ''),
        ('column translated down by 20 mm, 60 mm late', dict(z0=-0.020, zmax_fn=hi_after), {'F5 height': 'FAIL', 'F3': 'PASS'},
         ['60.00 mm above the inlet plane'], [], ''),
        ('redundant knot in the 0.20 m/s hold of its own fixture spec (machinery)', dict(table=tab_knot), REF, [], [], 'same'),
        ('lowest hold subdivided into 80 ms pieces, 60 mm in it', dict(table=tab_sub, zmax_fn=hi_after), {'F5 height': 'FAIL'},
         ['lowest velocity 0.050 m/s (2.020-2.500 s'], [], ''),
        ('one audit row of a fluidised window removed', dict(mutate=_drop(0.85)),
         {'F1': U, 'F2': U, 'F3': U, 'F5 count': U, 'F5 height': 'PASS'}, ['missing or rejected'], [], ''),
        ('dp_dyn NaN on one fluidised row (pressure identity)', dict(mutate=_set_field(0.85, 'dp_hydro_dyn', 2, 'NaN')),
         {'F1': U, 'F5 count': U, 'F5 height': 'PASS'}, ['pressure identity'], [], ''),
        ('one fixed step with a raised gradient, one of its rows removed',
         dict(grad_fn=lambda t, dp, u: dp / 0.1 * (1.5 if 1.5 < t <= 2.0 else 1.0), mutate=_drop(1.8)), {'F2': U}, [], [], ''),
        ('clean stop at 2.15 s (MaxNumStep)', dict(t_end=2.15), ALLU, ['the run ended at 2.1500 s, before the last table row'], [], ''),
        ('every other audit row removed, 60 mm on a removed row',
         dict(zmax_fn=lambda t: 0.060 if abs(t - 2.41) < 1e-9 else 0.0241, mutate=_thin(2)), {'F5 height': U, 'F1': U}, [], [], ''),
        ('seven identical copies of one audit row', dict(mutate=_dup(0.79, 7)), REF, [], [], ''),
        ('conflicting copy of one audit row', dict(mutate=_conflict(0.79)), {'F1': U}, ['conflicting copies'], [], ''),
        ('MaxSimTime 2.4999, the solver completes 2.5 s, 60 mm on the last row',
         dict(tmax=2.4999, zmax_fn=lambda t: 0.060 if t > 2.499 else 0.0241), {'F5 height': 'FAIL'}, [], [], ''),
        ('Reload_SimPar extends MaxSimTime, 60 mm at 2.4 s',
         dict(mutate=lambda L: _after_step(9200, ['Reloading Simulation parameters from file: _data/q2p1_param.dat'])(
             _reload_extension(2.3, 2.3)(L)), zmax_fn=lambda t: 0.060 if abs(t - 2.4) < 1e-9 else 0.0241),
         ALLU, ['reload during the run', 'MaxSimTime or MaxNumStep changed', 'measurement, not a verdict'], [],
         'f8fb7fd8: F5 height FAIL (witnessed)'),
        ('unfinished trailing EL_NEWTON_PAIR', dict(mutate=_append('EL_NEWTON_PAIR t= 2.5 step=')), REF, [], ['Traceback'], ''),
        ('last line without its line end (run in progress)', dict(mutate=_append('EL_BED_AUDIT t=   2.502500E+00 step= 10010 p_in_out='),
                                                                  last_complete=False, finished=False), ALLU, [], ['Traceback'], ''),
        ('malformed interior table row', dict(mutate=lambda L: [l.replace('3.00000000E-01', '3.0000000xE-01', 1)
                                                                if l.lstrip().startswith('EL_INLET_TABLE') and ' 1.00000000E+00' in l else l for l in L]),
         ALLU, ['table echo incomplete'], ['Traceback'], ''),
        ('particles lost on the last audit row (count and weight, witnessed)', dict(n_fn=lambda t: 9789 if t < 2.495 else 9691),
         {'F1': U, 'F5 count': 'FAIL'}, [], [], ''),
        ('wenyu fixture spec', dict(drag='wenyu'), {'F1': 'PASS'}, ['wenyu U_mf at the whole-bed'], [], ''),
        ('wenyu run against the gidaspow fixture spec', dict(drag='wenyu', _spec='base'), ALLU,
         ['drag model of the run wenyu [echo] differs from the spec gidaspow'], [], ''),
        ('no width echo, deck 3.0d0 with a comment', dict(width_echo=False, deck=DECK + 'SimPar@ELKernelWidthFactor = 3.0d0 ! mapping width\n',
                                                          _spec=dict(width='3.0')), {'F2': U}, ['kernel width 3.0d0 [deck]'], [], ''),
        ('no width echo, unreadable deck value', dict(width_echo=False, deck=DECK + 'SimPar@ELKernelWidthFactor = wide\n'),
         ALLU, ['kernel width factor of the run wide [deck] differs'], [], ''),
        ('NaN width echo', dict(width='NaN'), ALLU, ['kernel width factor of the run NaN [echo]'], [], ''),
        ('centre exactly at mid-height', dict(zmax_fn=lambda t: 0.050 if t > 2.3 else 0.0241), {'F5 height': 'PASS'},
         ['50.00 mm above the inlet plane', 'not above'], [], ''),
        ('column translated by +10 m (7-digit coordinates)', dict(z0=10.0), REF, [], [], ''),
        ('clean stop at 1.09 s, 59.5 mm in the preparation', dict(t_end=1.09, zmax_fn=lambda t: 0.0595 if t < 0.5 else 0.0241),
         ALLU, ['the run ended at 1.0900 s, before the last table row'], ['(witnessed)'], 'r5: F5 FAIL on the preparation'),
        ('segment from 1.0 s with steps renumbered from 1', dict(mutate=lambda L: _renumber(L, 4000)), ALLU,
         ['the time line of step 1 shows t = '], [], 'r5: F1 over the remaining steps'),
        ('Reload_SimPar resets the time from 1.75 to 1.25 s, unchanged dt, 60 mm late',
         dict(clock_fn=reset, t_end=3.0, zmax_fn=lambda t: 0.06 if t > 2.0 else 0.0241), ALLU, ['the time line of step'], [],
         'r5: F5 height PASS at 0.10 m/s'),
        ('repeated fixed velocity 0.30 -> 0.10 -> 0.30 -> 0.10 -> 0', dict(protocol=rep_proto), {'F3': U, 'F1': U, 'F2': U},
         ['not identified (rank-deficient: 1 distinct nonzero fixed velocity'], [],
         'r5: U_mf from one velocity; r7: F1 PASS (r8: no U_mf, so no F1 set)'),
        ('no fluidised step (0.14 .. 0.08 m/s)', dict(protocol=nofl_proto), {'F1': U, 'F2': U, 'F3': U},
         ['no step with dp* > 0.95 (no measured plateau'], [], 'r5: F2 FAIL, F3 PASS from the nominal plateau'),
        ('centre 0.05000002 m (inside the print band)', dict(zmax_fn=lambda t: 0.05000002 if t > 2.3 else 0.0241),
         {'F5 height': U}, ['not decidable at print precision'], [], 'r5: PASS'),
        ('centre 0.0501 m (above the band)', dict(zmax_fn=lambda t: 0.0501 if t > 2.3 else 0.0241), {'F5 height': 'FAIL'}, [], [], ''),
        ('dt changed late, 60 mm in the lowest hold', dict(clock_fn=dtchg, zmax_fn=lambda t: 0.06 if 2.2 < t < 2.26 else 0.0241),
         ALLU, ['time step not constant', 'measurement, not a verdict'], [], 'r5b: F5 height FAIL (witnessed)'),
        ('mu 2.4e-5 fixture spec, gradients from it', dict(mu=mu_e, grad_fn=grad_mu), {'F2': 'PASS'}, ['mu 0.24000000E-04'], [], ''),
        ('mu 2.4e-5 run against the 1.8e-5 fixture spec', dict(mu=mu_e, grad_fn=grad_mu, _spec='base'), ALLU,
         ['dynamic viscosity mu of the run'], [], 'r5: F2 31.91 % FAIL'),
        ('dp_dyn of one row +200 Pa, pressures unchanged', dict(mutate=_add_field(0.85, 'dp_hydro_dyn', 2, 200.0)),
         {'F1': U}, ['pressure identity dp_dyn = dp - dp_hydro violated'], [], 'r5: F1 PASS'),
        ('every count -1', dict(mutate=lambda L: [l.replace(' n= 9789 ', ' n= -1 ') for l in L]), ALLU,
         ['particle count -1 outside its domain'], [], 'r5: F5 count PASS'),
        ('0.20 m/s hold 1.02-1.40 s, sub-window 2 would cross into the next step', dict(table=sub_tab, t_end=2.4), {},
         ['   n/a  '], ['DIFF'], 'r5: sub2 0.740 DIFF'),
        ('no table echo, count change (ordinary report)', dict(echo_table=False, n_fn=lambda t: 9789 if t < 2.0 else 9788), ALLU,
         ['no inlet table echo', 'F5 count: particle counts [9788, 9789] (start-up population and the validated audit rows; measurement'], [],
         'f8fb7fd8 (r6 P1-7): the CLI exited before F5 count'),
        ('constant inlet 0.05 m/s, start bed collapsing from 60 mm', dict(table=[(0.0, 0.05)], t_end=1.0,
                                                                        zmax_fn=lambda t: 0.060 if t < 0.2 else 0.0241, _spec='base'),
         ALLU, ['echoed table has 1 rows'], ['(witnessed)'], '788ffbaf: F5 height FAIL (witnessed)'),
        # ---- review r6 ------------------------------------------------------------------------------------
        ('r6: superseded Prop Density 2 before the effective 1', dict(densities=[2.0, 1.0]), REF,
         ['Density assigned 2 times at start-up (the last is effective)'], [NO], 'f8fb7fd8: F2 51.74 % FAIL, F3 1.9679 FAIL'),
        ('r6: last knot at 2.5000002 s', dict(table=shifted, _spec='base'), ALLU, ['echoed table row 8 (t 2.50000020E+00'], [],
         'f8fb7fd8: four official PASSes'),
        ('r6: extra 0.01 m/s excursion 1.20-1.30 s, 60 mm in it', dict(table=excursion, _spec='base',
                                                                      zmax_fn=lambda t: 0.06 if 1.21 < t < 1.29 else 0.0241),
         ALLU, ['echoed table has 12 rows, the spec table 8'], ['(witnessed)'], 'f8fb7fd8: F5 FAIL'),
        ('r6: Reload_Prop doubles the viscosity at 1.5 s', dict(mutate=_after_step(6000, reload_lines),
                                                               grad_fn=lambda t, dp, u: dp / 0.1 * (2.0 if t > 1.5 else 1.0)),
         ALLU, ['reload during the run: Reloading Physical parameters', 'property echo during the run'], [],
         'f8fb7fd8: F1 PASS, F2 FAIL, F3 PASS, F5 PASS'),
        ('r6: StartingProc = 1', dict(starting_proc='1'), ALLU, ['start-up mode StartingProc = 1'], [], 'f8fb7fd8: four official PASSes'),
        ('r6: start-up lines of a new invocation appended', dict(mutate=lambda L: L + L[:20]), ALLU,
         ['solver output after the finish marker', "start-up echo 'EL"], [], 'f8fb7fd8 (on R_id400_new): F5 count PASS'),
        ('r6: two invocations concatenated, counts 9789 / 9788', dict(mutate=two_runs), ALLU,
         ['the start-up banner appears 2 times'], ['(witnessed)'], 'f8fb7fd8: witnessed F5 count FAIL'),
        ('r6: tied lowest holds 0.05 -> 0.20 -> 0.10 -> 0.05, 60 mm at 2.4 s', dict(table=tied, _spec='base', zmax_fn=hi_after),
         ALLU, ['echoed table row 1 (t'], [], 'f8fb7fd8: F5 height unavailable (preparation exception)'),
        ('r6: kernel width echo removed, no deck', dict(width_echo=False, deck=None), ALLU,
         ['kernel width factor not echoed and the deck is not readable', 'deck _data/q2p1_param.dat not readable'], [],
         'f8fb7fd8: official verdicts with the assumed width 2.5'),
        ('r6: last knot at 1e6 s, 2.5 s observed (arithmetic coverage)', dict(table=base_tab[:-1] + [(1.0e6, 0.05)]), ALLU,
         ['before the last table row'], ['Traceback'], 'f8fb7fd8: MemoryError under a 96 MiB limit'),
        ('r6: fluidised 1.80 and 0.99, plateau 1.395 (branch precedence)', dict(dps_fn=lambda u: overlap.get(round(u, 3), _dps(u))),
         {'F1': 'FAIL'}, ['F3 fixed branch 3 steps (U 0.050 .. 0.200)', 'over the 1 steps above U_mf 0.2018 m/s',
                          'sets settled in round 2'], [],
         'f8fb7fd8: the fit used 3 steps, the display 2; r7: fixed 2 steps (r8: 0.20 m/s leaves the F1 set, joins the fixed set)'),
        ('r6: gradz_win_dyn of one fixed row changed', dict(mutate=_add_field(1.85, 'gradz_win_dyn', 0, 500.0)),
         {'F2': U, 'F1': U}, ['gradient identity gradz_win_dyn = gradz_win - dp_hydro / Lz violated'], [],
         'f8fb7fd8: F2 from the edited row'),
        ('r6: dump load during the run', dict(mutate=_after_step(6000, [' Loading particle data from _dump/01'])), ALLU,
         ['a dump was loaded'], [], ''),
        # ---- set-up checks (they qualify the verdicts, they do not withhold them) ----------------------------
        ('set-up: restitution off', dict(contact='off'), REF,
         ['contacts: deviates: restitution=off', 'verdicts [set-up: contacts inelastic (restitution=off)]: F1 PASS'], [], ''),
        ('set-up: particle-particle friction 0.15', dict(contact='EL_CONTACT_PARAMS pair=particle-particle mu_static=0.15 '
                                                                'mu_dynamic=0.3 e=0.9 | pair=particle-wall wall=zwalls mu_static=0.3 '
                                                                'mu_dynamic=0.3 e=0.9 | restitution=on threshold=0.01'), REF,
         ['verdicts [set-up: contacts deviate (particle-particle mu_static=0.15)]'], [], ''),
        ('set-up: contacts not echoed', dict(contact=None), REF,
         ['not echoed (binary before pe pull request #59)', 'verdicts [set-up: contacts not echoed]'], [], ''),
        ('set-up: overlap 1.2 % of d at 0.07 s (preparation)', dict(overlap_fn=lambda t: 1.044e-5 if abs(t - 0.07) < 1e-9 else 3.8e-7),
         REF, ['exceeded during the preparation only', 'overlap 1.20 % of d at 0.070 s (preparation)]: F1 PASS'], [], ''),
        ('set-up: overlap 1.2 % of d at 2.4 s (evaluation window)', dict(overlap_fn=lambda t: 1.044e-5 if abs(t - 2.4) < 1e-9 else 3.8e-7),
         REF, ['exceeded inside an evaluation window', '(evaluation window)]'], [], ''),
        ('set-up: overlap exactly 1 % of d', dict(overlap_fn=lambda t: 8.7e-6 if abs(t - 2.4) < 1e-9 else 3.8e-7), REF,
         ['verdicts [set-up met]'], [], ''),
        ('set-up: no EL_CONTACT_STATS rows', dict(contact_stats=False), REF,
         ['overlap: not available', 'verdicts [set-up: overlap not available]'], [], ''),
        # ---- review r7 ------------------------------------------------------------------------------------
        ('r7: contacts mu 0 / e 1 for both pair kinds', dict(contact='EL_CONTACT_PARAMS pair=particle-particle mu_static=0 '
                                                                     'mu_dynamic=0 e=1 | pair=particle-wall wall=zwalls mu_static=0 '
                                                                     'mu_dynamic=0 e=1 | restitution=on threshold=0.01'), REF,
         ['contacts deviate (particle-particle mu_static=0; particle-particle mu_dynamic=0; particle-particle e=1'], [],
         '77914cfa: contacts met'),
        ('r7: sphere radius printed r=0.0004', dict(sph_r='0.0004'), ALLU, ['particle radius r of the run 0.0004 differs'], [],
         '77914cfa: recognised'),
        ('r7: 4895 spheres with a consistent weight', dict(n=4895), ALLU,
         ['particle count of the sphere setup 4895 differs from the spec 9789'], [], '77914cfa: official F1 and F5 PASS'),
        ('r7: 9789 at start-up, 9788 at the first audit, 9787 later', dict(n_fn=lambda t: 9788 if t < 1.5 else 9787),
         {'F5 count': 'FAIL'}, ['particle counts [9787, 9788, 9789] (start-up population'], [NO],
         '77914cfa: every gate unavailable'),
        ('r7: a positive count followed by zero', dict(n_fn=lambda t: 9789 if t < 2.0 else 0), {'F5 count': 'FAIL'},
         ['particle counts [0, 9789]'], [NO], "77914cfa: 'no validated audit row with N > 0'"),
        ('r7: gravity (1, 0, -9.81)', dict(gravity=(1.0, 0.0, -9.81)), ALLU, ['gravity component x of the run'], [],
         '77914cfa: official verdicts'),
        ('r7: one measured fixed step (F2 without the F3 fit)', dict(protocol=[(0.5, 0.30), (1.0, 0.25), (1.5, 0.20), (2.0, 0.10)]),
         {'F1': U, 'F2': U, 'F3': U}, ['round 1 of the U_mf iteration: the fit dp* = a U + b U^2 is not identified '
                                       '(rank-deficient: 1 distinct nonzero fixed velocity'], ['F2 worst |ratio_k - 1|'],
         '77914cfa: F2 and F3 unavailable (two fixed steps needed); r7: F2 FAIL (r8: no U_mf, so no sets)'),
        ('r7: contact rows before 0.2 s missing, 1.2 % at 0.07 s unseen',
         dict(overlap_fn=lambda t: 1.044e-5 if abs(t - 0.07) < 1e-9 else 3.8e-7,
              mutate=lambda L: [l for l in L if not (l.startswith('EL_CONTACT_STATS') and float(l.split()[2]) < 0.2)]),
         REF, ['overlap not available (coverage incomplete)]'], [], '77914cfa: overlap met'),
        ('r7: contact rows before 0.2 s missing, 1.2 % at 2.4 s observed',
         dict(overlap_fn=lambda t: 1.044e-5 if abs(t - 2.4) < 1e-9 else 3.8e-7,
              mutate=lambda L: [l for l in L if not (l.startswith('EL_CONTACT_STATS') and float(l.split()[2]) < 0.2)]),
         REF, ['exceeded inside an evaluation window (observed; coverage incomplete)', '(evaluation window, coverage incomplete)]'],
         [], ''),
        ('r7: one preparation audit row off the step clock', dict(mutate=_set_field(0.3, 't', 0, '3.001000E-01')), ALLU,
         ['audit row at step 1200 contradicts the step clock'], [], '77914cfa: the row a hole, F2 and F5 height PASS'),
        ('r7: sphere-creation echo twice', dict(mutate=lambda L: [l for x in L for l in ([x, x] if 'free sphere(s) created' in x else [x])]),
         ALLU, ['start-up echo sphere creation (EL terminal-velocity setup) appears 2 times'], [], '77914cfa: recognised'),
        ('r7: EL viscosity echo twice', dict(mutate=lambda L: [l for x in L for l in ([x, x] if 'diffusion operator rebuilt' in x else [x])]),
         ALLU, ['start-up echo EL fluid viscosity appears 2 times'], [], ''),
        ('r7: domain extents twice', dict(mutate=lambda L: [l for x in L for l in ([x, x] if 'Global CFD domain extents' in x else [x])]),
         ALLU, ['start-up echo Global CFD domain extents appears 2 times'], [], ''),
        # ---- round 8: F1 over the steps above the measured U_mf (spec amendment 2026-11-24) -----------------
        # per-step dp* of 14-step protocols 0.30 .. 0.04 m/s: the two runs (4 digits of post_r7.txt) and the peer series;
        # expected numbers from an independent iteration (numpy lstsq; unit/post_r8/indep.py) and the peer's
        ('r8: 80-row run values: 0.16 m/s at 0.9582 below U_mf, at the gate inside', dict(protocol=P14, dps_fn=_by_u(Z80)),
         {'F1': 'PASS'}, ['over the 7 steps above U_mf 0.1663 m/s (U 0.18 .. 0.30): 0.980843',
                          'PASS, at the gate (inside by 0.0008; standard error of the 7 step means 0.0021)',
                          'sample standard deviation 0.0055, standard error 0.0021', 'sets settled in round 2',
                          'in neither set: U 0.16 (dp* 0.9582)', 'mean over the 8 steps with dp* > 0.95: 0.97801',
                          '(reported, not gated)', ': F1 PASS (at the gate) | F2 '], [],
         'r7: F1 0.978013 over 8 steps FAIL'),
        ('r8: 40-row record values, at the gate outside', dict(protocol=P14, dps_fn=_by_u(R5)), {'F1': 'FAIL'},
         ['over the 7 steps above U_mf 0.1733 m/s (U 0.18 .. 0.30): 0.977743',
          'FAIL, at the gate (outside by 0.0023; standard error of the 7 step means 0.0034)', 'sets settled in round 1',
          ': F1 FAIL (at the gate) | F2 '], [], ''),
        ('r8: peer CFDEM 80 rows (stderr): plain FAIL', dict(protocol=P14, dps_fn=_by_u(PEER['CFDEM 80'])), {'F1': 'FAIL'},
         ['over the 7 steps above U_mf 0.1672 m/s (U 0.18 .. 0.30): 0.976243   (gate 1 +- 0.02: FAIL); ',
          'distance to the nearer limit -0.0038 (outside); at the gate when within the standard error: no',
          'mean over the 8 steps with dp* > 0.95: 0.97303'], [AT_GATE_NOTE, ', at the gate ('], 'peer: 7 steps, 0.9762, 0.1672 m/s'),
        ('r8: peer CFDEM 40 rows: plateau deficit, all 6 steps above U_mf', dict(protocol=P14, dps_fn=_by_u(PEER['CFDEM 40'])),
         {'F1': 'FAIL'}, ['over the 6 steps above U_mf 0.1906 m/s (U 0.20 .. 0.30): 0.954983   (gate 1 +- 0.02: FAIL); ',
                          'mean over the 5 steps with dp* > 0.95: 0.956360', 'in neither set: U 0.18 (dp* 0.9357)'],
         [AT_GATE_NOTE], 'peer: 6 steps, 0.9550, 0.1906 m/s; r7: 5 steps'),
        ('r8: peer DPMFoam: plateau 0.945, 1 of 6 steps above 0.95', dict(protocol=P14, dps_fn=_by_u(PEER['DPMFoam'])), {'F1': 'FAIL'},
         ['over the 6 steps above U_mf 0.1941 m/s (U 0.20 .. 0.30): 0.945183', 'mean over the 1 steps with dp* > 0.95: 0.955100'],
         [AT_GATE_NOTE], 'peer: 6 steps, 0.9452, 0.1941 m/s; r7: 1 step'),
        ('r8: peer MFiX: plain PASS', dict(protocol=P14, dps_fn=_by_u(PEER['MFiX'])), {'F1': 'PASS'},
         ['over the 7 steps above U_mf 0.1727 m/s (U 0.18 .. 0.30): 0.994843   (gate 1 +- 0.02: PASS); ', 'sets settled in round 1'],
         [AT_GATE_NOTE], 'peer: 7 steps, 0.9948, 0.1727 m/s'),
        ('r8: two steps above U_mf, 0.0025 outside: no scatter, plain FAIL', dict(dps_fn=_by_u([0.985, 0.97, 0.52, 0.26], PROTO)),
         {'F1': 'FAIL'}, ["no 'at the gate' annotation: 2 step(s) in the set, a scatter needs three",
                          'over the 2 steps above U_mf 0.1923 m/s (U 0.20 .. 0.30): 0.977500   (gate 1 +- 0.02: FAIL); '],
         [AT_GATE_NOTE], ''),
        ('r8: a step above U_mf with dp* < 0.8 x plateau (contradiction)',
         dict(protocol=P5, dps_fn=_by_u([1.5, 1.5, 0.99, 0.52, 0.26], P5)), {'F1': U, 'F2': U, 'F3': U},
         ['F1, F2, F3 unavailable: the branches contradict each other: step 2 (U 0.200) lies above U_mf 0.1923 m/s with '
          'dp* 0.9900 < 0.8 x the plateau 1.3300', 'cross-check, measured: mean over the 3 steps'], [], ''),
        ('r8: every step above 0.95 (empty fixed set)', dict(dps_fn=_by_u([1.0, 0.99, 0.97, 0.96], PROTO)), {'F1': U, 'F2': U, 'F3': U},
         ['F1, F2, F3 unavailable: round 1 of the U_mf iteration: no fixed-branch step'], [], 'r7: F1 PASS, F2 / F3 unavailable'),
        ('r8: fixed-branch fit without a crossing of dp* = 1', dict(dps_fn=_by_u([1.0, 1.0, 0.5, 0.45], PROTO)), {'F1': U, 'F2': U, 'F3': U},
         ['has no positive crossing with dp* = 1 (no U_mf)'], [], 'r7: F1 PASS, F3 unavailable'),
        ('r8: sets that cycle (0.14 m/s in and out of the fixed set)', dict(protocol=P7, dps_fn=_by_u(CYCLE, P7)), {'F1': U, 'F2': U, 'F3': U},
         ['the U_mf iteration does not settle: round 2 returns to the fluidised set of round 1 (a cycle of 2 sets)'], [], ''),
        ('r8: the iteration moves 0.16 m/s out of the F1 set, settles in round 2', dict(protocol=P5b, dps_fn=_by_u([0.99, 1.0, 0.96, 0.52, 0.26], P5b)),
         {'F1': 'PASS'}, ['over the 2 steps above U_mf 0.1923 m/s (U 0.20 .. 0.30): 0.995000', 'sets settled in round 2',
                          'in neither set: U 0.16 (dp* 0.9600)', 'mean over the 3 steps with dp* > 0.95: 0.983333'], [],
         'r7: F1 0.983333 over 3 steps'),
        ('r8: the same with ROUNDS = 1 (no settlement within the round limit)',
         dict(protocol=P5b, dps_fn=_by_u([0.99, 1.0, 0.96, 0.52, 0.26], P5b), _set=dict(ROUNDS=1)), {'F1': U, 'F2': U, 'F3': U},
         ['F1, F2, F3 unavailable: the U_mf iteration did not settle within 1 rounds'], [], ''),
        ('r8: U_mf above every step (empty F1 set)', dict(dps_fn=_by_u([1.0, 0.97, 0.3, 0.15], PROTO)), {'F1': U, 'F2': U, 'F3': U},
         ['F1, F2, F3 unavailable: round 1 of the U_mf iteration: no measured step above U_mf 0.3333 m/s'], [],
         'r7: F1 0.985000 PASS'),
        # ---- round 9 (review of r8): a step velocity at U_mf is not decided by round-off -----------------------
        ('r9: fixed branch dp* = 5 U, U_mf = 0.20 = a step velocity: not determined', dict(protocol=P14, dps_fn=_by_u(TIE5)),
         {'F1': U, 'F2': U, 'F3': U}, ['F1, F2, F3 unavailable: round 1 of the U_mf iteration: step 5 (U 0.200) coincides with U_mf 0.200000 m/s',
                                       'whether it lies above U_mf is not determined'], [],
         'r8: U_mf 0.19999999999999896, the 0.20 step in the set, 0.976667 FAIL (strict rule: 5 steps, 0.982 PASS)'),
        ('r9: the same with dp* = 5.001 U (U_mf 0.19996): the 0.20 m/s step is above', dict(protocol=P14, dps_fn=_by_u(_tie(5.001))),
         {'F1': 'FAIL'}, ['over the 6 steps above U_mf 0.2000 m/s (U 0.20 .. 0.30): 0.976667', 'sets settled in round 2'], [], ''),
        ('r9: the same with dp* = 4.999 U (U_mf 0.20004): the 0.20 m/s step is not above', dict(protocol=P14, dps_fn=_by_u(_tie(4.999))),
         {'F1': 'PASS'}, ['over the 5 steps above U_mf 0.2000 m/s (U 0.22 .. 0.30): 0.982000', 'sets settled in round 1'], [], ''),
    ]


def _run_case(run, spec):
    out, err = io.StringIO(), io.StringIO()
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            main([run], spec=spec)
    except SystemExit:
        pass
    except BaseException:                            # noqa: BLE001
        import traceback
        out.write('Traceback\n' + traceback.format_exc())
    return out.getvalue(), err.getvalue()


def verdicts_of(text):
    m = re.search(r'^verdicts(?: \[[^]]*\])?: (.*)$', text, re.M)
    if not m:
        return {}
    V = {}
    for part in m.group(1).replace(OUTSIDE, '').split(' | '):
        g, v = part[:-len(AT_GATE_NOTE)].rsplit(' ', 1) if part.endswith(AT_GATE_NOTE) else part.rsplit(' ', 1)
        V[g] = v
    return V


def _unit_checks():
    """(name, ok, detail): the closure crossing (bracket from zero) and the spec-11 table definition"""
    out = []
    eps = 0.4
    u = closure_umf(eps, 0.87e-3, 1001.0, 1000.0, 0.01, 9.81)
    target = (1 - eps) * 1.0 * 9.81
    g = closure_grad(u, eps, 0.87e-3, 1000.0, 0.01)
    out.append(('closure crossing of a liquid (rho_p 1001, rho_f 1000, mu 0.01) bracketed from zero',
                finite(u) and abs(g / target - 1) < 1e-9, f'U_mf {u:.6g} m/s, gradient {g:.6g} vs {target:.6g} Pa/m '
                '(f8fb7fd8: the bracket end 1e-6)'))
    out.append(('no crossing for rho_p < rho_f: NaN', not finite(closure_umf(eps, 0.87e-3, 900.0, 1000.0, 0.01, 9.81)), ''))
    pe = [('0', '0.3', False), ('1', '0.9', False), ('0.0004', '0.000435', False), ('0.3', '0.3', True),
          ('0.2500D-03', '0.00025', True), ('0.10000000E+01', '1.0', True), ('2.50000020E+00', '2.5', False),
          ('0.000000E+00', '0', True), ('1.0005', '1.00049', True), ('1.02000000E+00', '1.02', True)]
    bad = [(a, b) for a, b, want in pe if print_equal(a, Decimal(b)) != want]
    out.append(('print_equal: a coarse token never stands for a finer reference (%d pairs)' % len(pe), not bad,
                ('wrong: %s' % bad) if bad else "('0', 0.3), ('1', 0.9), ('0.0004', 0.000435) refused (77914cfa: accepted)"))
    # round 9: gate F1 on a limit, and a step velocity at U_mf (review of r8)
    on = [f1_gate([lim] * 7, SPEC11) for lim in SPEC11.F1_LIMITS]
    out.append(('F1 of seven equal step means on either limit: PASS, distance 0, at the gate',
                all(g['verdict'] == 'PASS' and g['dist'] == 0 and g['se'] == 0 and g['at_gate'] for g in on)
                and exact_mean([0.98] * 7) == 0.98, f"plain sum: {sum([0.98] * 7) / 7!r} (r8: PASS without the annotation)"))
    off = [f1_gate(v, SPEC11) for v in ([0.9799] * 7, [1.0201] * 7, [0.98, 0.98], [0.985, 0.975, 0.98])]
    out.append(('F1 off a limit without scatter: FAIL, no annotation; two steps: none; three steps on the limit: annotated',
                [(g['verdict'], g['at_gate']) for g in off] == [('FAIL', False), ('FAIL', False), ('PASS', False), ('PASS', True)], ''))
    rows = lambda dps: [dict(k=i, U=u, dps=d, measured=True) for i, ((_, u), d) in enumerate(zip(P14, dps))]   # noqa: E731
    quad = [0.982] * 5 + [0.95, 0.90, 0.80] + [4 * u + 5 * u * u for u in (0.14, 0.12, 0.10, 0.08, 0.06, 0.04)]
    why = [branch_sets(rows(v), SPEC11)['reason'] for v in (TIE5, quad)]
    out.append(('a step velocity at U_mf (fits 5 U and 4 U + 5 U^2, both crossing at 0.20): sets not determined',
                all('coincides with U_mf' in w and 'step 5 ' in w for w in why), '; '.join(w[:44] for w in why)))
    near = [branch_sets(rows(_tie(c)), SPEC11) for c in (5 * (1 + 3e-6), 5 * (1 - 3e-6))]
    out.append(('three band widths either side of the tie: decided (6 steps with 0.20 m/s, 5 without)',
                [(b['reason'], len(b['flu'])) for b in near] == [('', 6), ('', 5)], ''))
    s = SPEC11
    out.append(('spec 11: 28 rows from make_inlet_table, 14 holds, lowest 7.02-7.5 s at 0.04 m/s',
                len(s.table) == 28 and len(s.holds) == 14 and (s.lowest['a'], s.lowest['b']) == (Fraction('7.02'), Fraction('7.5'))
                and s.lowest['uf'] == 0.04, f'{len(s.table)} rows, {len(s.holds)} holds'))
    return out


def self_test(keep=None):
    cases = self_test_cases()
    nfail = 0
    ref_out = None
    ctx = tempfile.TemporaryDirectory(prefix='fl11_post_selftest_') if keep is None else contextlib.nullcontext(keep)
    with ctx as tmp:
        for i, (name, kw, expect, must, mustnot, note) in enumerate(cases):
            kw = dict(kw)
            sp = kw.pop('_spec', None)
            attrs = kw.pop('_set', {})
            spec = _fixture_spec({}) if sp == 'base' else _fixture_spec(kw, sp)
            for key, val in attrs.items():         # gate constants of this fixture spec (Spec class attributes)
                setattr(spec, key, val)
            run = _write_run(os.path.join(tmp, 'run%02d' % i), **kw)
            out, err = _run_case(run, spec)
            text = out + err
            V = verdicts_of(out)
            body = '\n'.join(l for l in out.split('\n')[1:] if not l.startswith('protocol:'))
            if i == 0:
                ref_out = body
            bad = [f'{g}: {V.get(g)} instead of {v}' for g, v in expect.items() if V.get(g) != v]
            bad += [f'missing {m!r}' for m in must if m not in text] + [f'unexpected {m!r}' for m in mustnot if m in text]
            if 'PASS' in V.values() and 'NO OFFICIAL' in out:
                bad.append('a PASS from a run outside the supported domain')
            if 'FAIL' in V.values() and 'NO OFFICIAL' in out and 'fails only "finished"' not in out:
                bad.append('a FAIL from a run that fails more than "finished"')
            if note == 'same' and body != ref_out:
                bad.append('output differs from the reference case')
            print(f"self-test {i:2d} {'ok  ' if not bad else 'FAIL'} {name}" + ('' if not bad else f": {bad}"))
            if bad:
                nfail += 1
                print('    ' + text.replace('\n', '\n    '))
            with open(os.path.join(run, 'case.txt'), 'w') as fh:
                fh.write('%s\nold: %s\n' % (name, note))
    units = _unit_checks()
    for name, ok, detail in units:
        print(f"unit      {'ok  ' if ok else 'FAIL'} {name}" + (f" ({detail})" if detail else ''))
        nfail += not ok
    print(f"fl11_post self-test: {len(cases) + len(units) - nfail}/{len(cases) + len(units)} passed "
          f"({len(cases)} fixture logs with fixture specs, {len(units)} unit checks)")
    return 1 if nfail else 0


# ===================================================================================================
# metamorphic tests on a real log
# ===================================================================================================
KEEP_PREFIX = ('EL_BED_AUDIT', 'EL_FORCE_BUDGET', 'EL_NEWTON_PAIR', 'EL_INLET_TABLE', 'EL_CONTACT_', 'EL ', 'Global CFD domain extents',
               'time:', 'Prop', 'StartingProc', 'TimeStep', 'StartSimTime', 'MaxSimTime', 'MaxNumStep', 'Reloading')


def _filtered(lines):
    """the lines the tool reads (the rest of the log does not change its output; checked once)"""
    return [l for l in lines if l.lstrip().startswith(KEEP_PREFIX) or FINISH in l or BANNER in l
            or 'Loading particle data from' in l]


def _table_edit(lines, new_rows):
    out, done = [], False
    for l in lines:
        if l.lstrip().startswith('EL_INLET_TABLE'):
            if not done:
                out += ['  EL_INLET_TABLE %s %s' % (_e(t, 16, 8), _e(u, 16, 8)) for t, u in new_rows]
                done = True
            continue
        m = re.search(r'\((\d+) rows', l) if l.startswith('EL inlet velocity table') else None
        if m:
            l = l.replace('(%s rows' % m.group(1), '(%d rows' % len(new_rows))
        out.append(l)
    return out


AUDIT_FIELDS = [('t', 0), ('step', 0), ('p_in_out', 0), ('p_in_out', 1), ('area_in_out', 0), ('area_in_out', 1),
                ('dp_hydro_dyn', 0), ('dp_hydro_dyn', 1), ('dp_hydro_dyn', 2), ('n', 0), ('weight_z_floor_z', 0),
                ('weight_z_floor_z', 1), ('z50_95_99_max', 0), ('z50_95_99_max', 1), ('z50_95_99_max', 2),
                ('z50_95_99_max', 3), ('eps_bed', 0), ('window', 0), ('window', 1), ('gradz_win', 0), ('floor_bal_z', 0),
                ('eps_kern_slip_win', 0), ('eps_kern_slip_win', 1), ('win_margin', 0), ('gradz_win_dyn', 0)]
HOLE = {('t', 0), ('step', 0), ('p_in_out', 0), ('p_in_out', 1), ('dp_hydro_dyn', 0), ('dp_hydro_dyn', 1),
        ('dp_hydro_dyn', 2), ('n', 0), ('weight_z_floor_z', 0)}
F2_FIELDS = {('z50_95_99_max', 2), ('window', 0), ('window', 1), ('eps_kern_slip_win', 0), ('win_margin', 0)}


def self_test_real(rundir, nvar=300, seed=20261122):
    """metamorphic tests on a real log: equivalence (byte-identical, or the verdict line identical), degrading
    (the expected verdict of every gate from its own dependencies), run history (every official gate
    unavailable with the expected reason)"""
    rnd = random.Random(seed)
    with open(os.path.join(rundir, 'run_slurm.log'), errors='replace') as fh:
        full = fh.read().split('\n')
    if full and full[-1] == '':
        full.pop()
    deck = open(os.path.join(rundir, '_data', 'q2p1_param.dat'), errors='replace').read()
    base = _filtered(full)
    ref_text, ref_V, ev0 = run_lines(base, deck)
    full_text, _, _ = run_lines(full, deck)
    vline = lambda t: re.search(r'^verdicts.*$', t, re.M).group(0)   # noqa: E731
    print(f"real log {rundir}: {len(full)} lines, {len(base)} read by the tool; intact verdicts: "
          + ' | '.join(f'{g} {ref_V[g]}' for g in GATES) + ('' if ev0['official'] else '  (not recognised: '
                                                             + '; '.join(ev0['precondition'][:2]) + ')'))
    nfail = 0
    if full_text != ref_text:
        nfail += 1
        print('FAIL: the full log and the lines the tool reads give different output')
    spec, clock0 = ev0['spec'], ev0['clock']
    windows = [(x['t1'] - spec.EVAL_LEN, x['t1']) for x in ev0['steps']]
    fixedw = [(x['t1'] - spec.EVAL_LEN, x['t1']) for x in ev0['steps'] if x['branch'] == 'fixed']
    lowest = [(spec.lowest['a'], spec.lowest['b'])]
    audit_idx = [i for i, l in enumerate(base) if l.startswith('EL_BED_AUDIT')]
    time_idx = {int(PAT_TIME.match(l).group(3)): i for i, l in enumerate(base) if PAT_TIME.match(l)}
    marker = max((i for i, l in enumerate(base) if FINISH in l), default=len(base) - 1)
    dtc = float(clock0.dt)
    tk_of = lambda i: clock0.tk(int(re.search(r'step=\s*(\d+)', base[i]).group(1)))   # noqa: E731
    inside = lambda t, ivs: any(a <= t <= b for a, b in ivs)                          # noqa: E731
    ALLU = {g: 'unavailable' for g in GATES}
    gw_dyn_rows = all(r['has_gw_dyn'] for r in clock0.sorted)

    def hole_expect(ts):
        E = dict(ref_V)
        if any(inside(t, windows) for t in ts):
            E.update(F1='unavailable', F2='unavailable', F3='unavailable')
        if ref_V['F5 height'] != 'FAIL' and any(inside(t, lowest) for t in ts):
            E['F5 height'] = 'unavailable'
        if ref_V['F5 count'] != 'FAIL':
            E['F5 count'] = 'unavailable'
        return E

    def field_expect(t, key, idx):
        """the gates that depend on this field of a row at step time t"""
        if (key, idx) in HOLE or (gw_dyn_rows and key in ('gradz_win', 'gradz_win_dyn')):
            return hole_expect([t])
        E = dict(ref_V)
        if (key, idx) == ('area_in_out', 0) and inside(t, windows):
            E.update(F1='unavailable', F2='unavailable', F3='unavailable')
        if inside(t, fixedw):
            if (key, idx) in F2_FIELDS or (not gw_dyn_rows and (key, idx) == ('gradz_win', 0)):
                E['F2'] = 'unavailable'
            if (key, idx) == ('z50_95_99_max', 2):
                E['F3'] = 'unavailable'
        if (key, idx) == ('z50_95_99_max', 3) and inside(t, lowest) and ref_V['F5 height'] != 'FAIL':
            E['F5 height'] = 'unavailable'
        return E

    first_prop = next(i for i, l in enumerate(base) if l.startswith('Prop'))
    dens = next(i for i, l in enumerate(base) if re.match(r'^Prop\s+-\s+Density', l))
    eq_kinds = ['trailing record', 'row duplicated', 'superseded Prop line']
    dg_kinds = ['one row deleted', 'rows deleted', 'conflicting copy', 'field NaN', 'log truncated', 'audit thinned',
                'count perturbed', 'weight perturbed', 'finish marker removed']
    if ev0['official'] and ref_V['F5 height'] == 'PASS':
        dg_kinds += ['witness kept, other rows deleted', 'witness in an unfinished run']
    hist_kinds = ['clean early stop', 'suffix renumbered', 'time reset (unchanged dt)', 'dt changed', 'property reload',
                  'StartingProc 1', 'second invocation appended', 'second banner', 'knot inserted', 'last knot removed']
    print('transforms: equivalence %s; degrading %s; run history %s; skipped (not applicable): %s'
          % (', '.join(eq_kinds), ', '.join(dg_kinds), ', '.join(hist_kinds),
             'none' if len(dg_kinds) == 11 else 'witness transforms (no official F5 height PASS)'))

    def equiv(kind):
        L = list(base)
        if kind == 'row duplicated':
            i = rnd.choice(audit_idx)
            L.insert(i + 1, L[i])
            return L, True, 'bytes'
        if kind == 'superseded Prop line':
            L.insert(dens, 'Prop    - Density =' + _fe(rnd.choice([2.0, 0.5, 1000.0])) + _fe(1.0))
            return L, True, 'verdicts'
        tail = rnd.choice(['EL_NEWTON_PAIR t= 7.5 step=', 'EL_BED_AUDIT t=   7.502500E+00 step= 30010 p_in_out=',
                           'EL_FORCE_BUDGET step= 30001 drag=', 'time:  0.7500D+01/', '  EL_INLET_TABLE   7.6'])
        return base + [tail], rnd.random() < 0.5, 'bytes'

    def degrade(kind):
        L = list(base)
        if kind == 'one row deleted':
            i = rnd.choice(audit_idx)
            t = tk_of(i)
            del L[i]
            return L, True, hole_expect([t])
        if kind == 'rows deleted':
            ii = rnd.sample(audit_idx, min(len(audit_idx), rnd.randint(2, 20)))
            ts = [tk_of(i) for i in ii]
            for i in sorted(ii, reverse=True):
                del L[i]
            return L, True, hole_expect(ts)
        if kind == 'conflicting copy':
            i = rnd.choice(audit_idx)
            tok = L[i].split()
            present = [f for f in AUDIT_FIELDS if f[0] + '=' in tok and f[0] not in ('t', 'step', 'n')]
            key, idx = rnd.choice(present)
            j = tok.index(key + '=') + 1 + idx
            v = num(tok[j])
            tok[j] = _e(v * 1.001 + 1e-9) if finite(v) else _e(1.0)
            L.insert(i + 1, ' '.join(tok))
            return L, True, hole_expect([tk_of(i)])
        if kind == 'field NaN':
            i = rnd.choice(audit_idx)
            tok = L[i].split()
            key, idx = rnd.choice([f for f in AUDIT_FIELDS if f[0] + '=' in tok])
            t = tk_of(i)
            L[i] = _set_tok(L[i], key, idx, 'NaN')
            return L, True, field_expect(t, key, idx)
        if kind == 'log truncated':
            L = L[:rnd.randrange(max(1, marker)) + 1]
            done = True
            if rnd.random() < 0.5 and L[-1]:
                L[-1] = L[-1][:rnd.randrange(len(L[-1]))]
                done = False
            return L, done, dict(ALLU)
        if kind == 'audit thinned':
            m = rnd.randint(2, 5)
            drop = set(audit_idx[rnd.randrange(m)::m])
            return [l for i, l in enumerate(L) if i not in drop], True, hole_expect([tk_of(i) for i in drop])
        if kind in ('count perturbed', 'weight perturbed'):
            i = rnd.choice(audit_idx)
            tok = L[i].split()
            if kind == 'count perturbed':
                j = tok.index('n=') + 1
                tok[j] = str(int(tok[j]) + rnd.randint(1, 100))
            else:
                j = tok.index('weight_z_floor_z=') + 1
                tok[j] = _e(num(tok[j]) * (1 + rnd.choice([-1, 1]) * rnd.uniform(1e-4, 1e-2)))
            L[i] = ' '.join(tok)
            return L, True, hole_expect([tk_of(i)])
        if kind == 'finish marker removed':
            return [l for l in L if FINISH not in l], True, dict(ALLU)
        w = rnd.choice([i for i in audit_idx if inside(tk_of(i), lowest)])
        L[w] = _set_tok(L[w], 'z50_95_99_max', 3, _e(ev0['zmin'] + 0.06))
        if kind == 'witness in an unfinished run':
            return [l for l in L if FINISH not in l], True, dict(ALLU, **{'F5 height': 'FAIL'})
        others = rnd.sample([i for i in audit_idx if i != w], rnd.randint(1, 20))
        ts = [tk_of(i) for i in others]
        for i in sorted(others, reverse=True):
            del L[i]
        E = hole_expect(ts)
        E['F5 height'] = 'FAIL'
        return L, True, E

    def history(kind):
        L = list(base)
        K = max(time_idx)
        if kind == 'clean early stop':
            k = rnd.randrange(max(1, K // 10), K - 1)
            return L[:time_idx[k] + 1] + [' PP3D_LES has successfully finished. '], 'before the last table row'
        if kind == 'suffix renumbered':
            f = clock0.freq or 1
            k0 = (rnd.randrange(int(0.1 * K), int(0.6 * K)) // f) * f
            return _renumber(L, k0), 'the time line of step'
        if kind in ('time reset (unchanged dt)', 'dt changed'):
            k1 = rnd.randrange(int(0.2 * K), int(0.8 * K))
            sh = rnd.choice([0.25, 0.5, 1.0]) * max(clock0.t_last, 0.1)
            t1 = float(clock0.tk(k1))
            new_t = (lambda k: float(clock0.tk(k)) - sh) if kind.startswith('time') else (lambda k: t1 + (k - k1) * 2 * dtc)
            out = []
            for l in L:
                mm = PAT_TIME.match(l)
                if mm and int(mm.group(3)) > k1:
                    l = _time_line(l, new_t(int(mm.group(3))), None if kind.startswith('time') else 2 * dtc)
                elif l.startswith(('EL_BED_AUDIT', 'EL_NEWTON_PAIR')) and int(re.search(r'step=\s*(\d+)', l).group(1)) > k1:
                    tok = l.split()
                    tok[2] = _e(new_t(int(re.search(r'step=\s*(\d+)', l).group(1))))
                    l = ' '.join(tok)
                out.append(l)
            return out, 'the time line of step' if kind.startswith('time') else 'time step not constant'
        if kind == 'property reload':
            k = rnd.randrange(int(0.1 * K), int(0.9 * K))
            i = time_idx[k] + 1
            return L[:i] + ['Reloading Physical parameters from file: _data/q2p1_param.dat', L[dens]] + L[i:], 'reload during the run'
        if kind == 'StartingProc 1':
            return [re.sub(r'^StartingProc = 0', 'StartingProc = 1', l) for l in L], 'StartingProc = 1'
        if kind == 'second invocation appended':
            return L + L[:rnd.randint(1, first_prop)], 'after the finish marker'
        if kind == 'second banner':
            i = rnd.choice(audit_idx)
            return L[:i] + ['  |          %s     |' % BANNER] + L[i:], 'banner appears 2 times'
        table = [(float(num(t)), float(num(u))) for t, u in ev0['log'].table]
        if kind == 'knot inserted' and len(table) < 2:          # a one-row table: a knot after its row
            return _table_edit(L, table + [(table[-1][0] + 1.0, table[-1][1])]), 'echoed table has'
        if kind == 'knot inserted':
            j = rnd.randrange(len(table) - 1)
            (ta, ua), (tb, ub) = table[j], table[j + 1]
            t = round(rnd.uniform(ta, tb), 6)
            return _table_edit(L, sorted(table + [(t, ua + (ub - ua) * (t - ta) / (tb - ta))])), 'echoed table has'
        return _table_edit(L, table[:-1]), 'echoed table has'

    stats, n_un = {}, 0
    for v in range(nvar):
        cls = v % 3
        key, E = 'variant %d' % v, None
        try:
            if cls == 0:
                kind = rnd.choice(eq_kinds)
                key = 'equivalence: ' + kind
                L, complete, how = equiv(kind)
                text, V, _ = run_lines(L, deck, last_complete=complete)
                ok = text == ref_text if how == 'bytes' else vline(text) == vline(ref_text)
            elif cls == 1:
                kind = rnd.choice(dg_kinds)
                key = 'degrading: ' + kind
                L, complete, E = degrade(kind)
                text, V, _ = run_lines(L, deck, last_complete=complete)
                ok = all(V.get(g) == E[g] for g in GATES)
            else:
                kind = rnd.choice(hist_kinds)
                key = 'run history: ' + kind
                L, reason = history(kind)
                text, V, _ = run_lines(L, deck)
                ok = all(V.get(g) == 'unavailable' for g in GATES) and (reason in text or not ev0['official'])
        except BaseException:                        # noqa: BLE001
            import traceback
            text, V, ok = 'Traceback ' + traceback.format_exc(), {}, False
        st = stats.setdefault(key, [0, 0, 0])
        st[0] += 1
        if V and cls and any(V.get(g) == 'unavailable' for g in GATES):
            st[2] += 1
            n_un += 1
        if not ok:
            st[1] += 1
            nfail += 1
            if st[1] <= 2:
                print(f"FAIL variant {v} ({key}): verdicts {V}" + (f", expected {E}" if E else ''))
                if text.startswith('Traceback'):
                    print(text[-1500:])
                elif cls == 0:
                    import difflib
                    print('\n'.join(list(difflib.unified_diff(ref_text.splitlines(), text.splitlines(), lineterm=''))[:20]))
                else:
                    print('   reasons: ' + '; '.join(re.findall(r'^   - (.*)$', text, re.M))[:400])
    neq = sum(s[0] for k, s in stats.items() if k.startswith('equivalence'))
    ndg = sum(s[0] for k, s in stats.items() if k.startswith('degrading'))
    print(f"{nvar} variants ({neq} equivalence, {ndg} degrading, {nvar - neq - ndg} run history); non-equivalence variants "
          f"with at least one verdict unavailable: {n_un}; failures: {nfail}")
    for k in sorted(stats):
        print(f"   {k:46s} {stats[k][0]:4d} run, {stats[k][1]:3d} failed, {stats[k][2]:3d} with a verdict unavailable")
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
