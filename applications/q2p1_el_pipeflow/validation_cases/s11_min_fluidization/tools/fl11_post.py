#!/usr/bin/env python3
"""Case 11 (minimum fluidisation) evaluation of a q2p1_el_pipeflow run.

The tool is built from five pieces, each defined once (design note section 7.4):
  A  record reader   one tokenizer for every EL_ record the tool reads (EL_BED_AUDIT, EL_FORCE_BUDGET,
                     EL_NEWTON_PAIR, EL_INLET_TABLE): a record is complete or it is nothing; an unfinished
                     trailing record (after the last complete one, or a last line without its line end)
                     is ignored with a note on stderr, a malformed interior record is diagnosed. The
                     start-up echo (kernel width factor, drag model, fluid gravity, bed audit frequency,
                     table row count, domain extents), the time lines and the finish marker are read in
                     the same pass; the deck (_data/q2p1_param.dat) with the application's syntax
                     (SimPar@Key = value: last '@', last '=', first list-directed token of the value).
  B  protocol        the solver's rule U(t) (EL_INLET_TABLE_VALUE: linear between the rows, first and last
                     value held beyond the ends) on the echoed table, reduced to maximal constant-velocity
                     holds (adjacent equal-valued segments merged BEFORE the 0.1 s threshold), the terminal
                     hold up to the end of the run. Every gate and the movie use these holds; equivalent
                     tables give identical output.
  C  clock / audit   the expected audit steps (multiples of the echoed audit frequency up to the last
                     completed step) and their times t0 + k dt where a constant time step is established
                     (one printed dt on every time line; the exact deck value when it agrees). Observed
                     rows are validated against them (step id vs time, count vs weight), identical
                     duplicates collapsed, conflicting duplicates rejected; a missing expected row is a
                     coverage hole. The end of the run is the last step the solver completed (its loop
                     finishes the step that crosses MaxSimTime; a reload may move MaxSimTime); without the
                     finish marker the terminal hold is open.
  D  per step        the evaluation window (last 0.3 s of the step) must lie inside the hold; statistics,
                     coverage and extrema use the same validated rows; buoyant weight, area and particle
                     count must be constant in the window (they normalise that step). A step is
                     'measured' or 'unmeasurable' with the reason.
  E  verdicts        PASS / FAIL / unavailable (reason). A gate is unavailable whenever a protocol step that
                     could belong to its set is unmeasurable (the subset value is printed as a measurement,
                     not as a verdict). F5 height and F5 count are independent of F1-F3 and of each other.

Gates of the spec (section 5, as amended on 2026-11-20):
  F1  plateau dp A / buoyant weight, mean over the steps with dp* > 0.95           (1 +- 0.02)
  F2  fixed branch: dynamic window gradient -<dp/dz> (gradz_win_dyn; gradz_win for older logs when
      dp_hydro = 0 on every contributing line) against the closure at the MAPPED gas fraction (+-5 %),
      when every contributing window excludes max(2 d, mapping support) at both ends (lower end - inlet
      plane, z99 + r - upper end, win_margin or 2 d before 627caf47; tolerance from the print precision
      of the coordinates), every mapped sample is valid (finite, 0 < eps < 1), and the closure is positive
      and the ratio finite on every fixed step.
  F3  U_mf where the fit dp* = a U + b U^2 of the fixed branch (through the origin) meets dp* = 1, against
      the closure U_mf at the whole-bed porosity 1 - N V_p / (A (z99 + r - zmin))      (+-5 %)
  F4  U_mf against the measurement 0.182 m/s (reported)
  F5  count: particle count constant over the run (every expected audit row validated, run finished);
      height: no centre above the mid-column zmin + Lz/2 (echoed extents) during the lowest-velocity
      hold(s) of the protocol (all tied holds; a centre at mid-height is not above; tolerance: print
      precision of the coordinates). An observed exceedance or count change is a FAIL even when the
      coverage is incomplete.
The closure of F2 / F3 is the drag model of the start-up echo (deck otherwise); a --closure that differs
gives a labelled comparison, not a verdict.
A --steps file never feeds the official gates of a log with a table echo: it is validated (three
columns, multiplicity) and gives a separately labelled summary when it differs. A log without an echo
is evaluated from a valid steps file (no duplicate or overlapping rows), every verdict labelled.

usage: fl11_post.py RUNDIR [--steps FILE] [--d 0.87e-3] [--rhop 860] [--rhof 1.0] [--mu 1.8e-5]
                    [--g 9.81] [--lz LZ] [--closure NAME] [--flag 0.02] [--series OUT.csv]
       fl11_post.py --self-test                synthetic logs (the counterexamples of reviews r2-r4)
       fl11_post.py --self-test-real RUNDIR [--variants N] [--seed S]
                                               metamorphic tests on a real log: equivalence transforms
                                               give identical output, degrading transforms leave every
                                               verdict equal to the intact one or unavailable
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fl11_closures import CLOSURES, closure_grad, closure_umf   # noqa: E402  (re-exported)

EVAL_VERSION = 'fl11_post r4 2026-11-22'   # part of the movie cache signature
NAN = float('nan')
INF = float('inf')
FROM_STEPS = ' [protocol from the steps file, not from the echo]'
GATES = ('F1', 'F2', 'F3', 'F5 height', 'F5 count')
HOLD_MIN = 0.1           # s: a hold of this length or less is part of a transition, not a step
EVAL_LEN = 0.3           # s: evaluation window at the end of a step
FIRST_LEN = 0.5          # s: the first hold (preparation + step 0) is evaluated over its last 0.5 s
MIN_SAMPLES = 10         # expected audit rows in a window used for statistics or extrema
TOL_D4 = 5e-4            # relative print precision of the 4-digit D format of the time line


def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def pos(x):
    """x when it is a positive finite number (a usable denominator), else nan"""
    return x if finite(x) and x > 0 else NAN


def div(a, b):
    """a / b with b a positive finite number; nan otherwise (never an exception)"""
    b = pos(b)
    return a / b if finite(b) and finite(a) else NAN


def es_prec(x, digits=6):
    """half the spacing of the ES(w).digits print of x (0 for an exact zero)"""
    if not finite(x) or x == 0:
        return 0.0
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(x))) - digits)


def mean(v):
    return sum(v) / len(v) if v else NAN


def _fmin(v):
    v = [x for x in v if finite(x)]
    return min(v) if v else NAN


def _fmax(v):
    v = [x for x in v if finite(x)]
    return max(v) if v else NAN


# ===================================================================================================
# A. record reader
# ===================================================================================================
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
        return NAN                                    # Fortran overflow field: no value
    if not NUM_RE.match(tok):
        raise ValueError(tok)
    return float(tok.replace('D', 'E').replace('d', 'e'))


def integer(tok):
    if not re.fullmatch(r'[-+]?\d+', tok):
        raise ValueError(tok)
    return int(tok)


# tag: (required [(key, count, converter)], optional [...]); other keys are ignored
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
}
KEY_TOK = re.compile(r'^[A-Za-z_0-9]+=$')


def parse_keyed(tokens, spec):
    """(record, '') for a complete record, (None, reason) otherwise"""
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
    return rec, ''


PAT_TIME = re.compile(r'^\s*time:\s*(%s)\s*/\s*(%s)\s*\|\s*itns:\s*(\d+)\s*/\s*(\d+)\s*\|\s*dt:\s*(%s)'
                      % (NUMG, NUMG, NUMG))
ECHO = {
    'width': re.compile(r'^EL kernel width factor\s*=\s*(\S+)\s*$'),
    'drag': re.compile(r'^EL drag model\s*=\s*(\S+)\s*$'),
    'fluid_gravity': re.compile(r'^EL fluid gravity\s*=\s*(\S+)\s*$'),
    'audit_freq': re.compile(r'^EL bed audit freq\s*=\s*(\S+)\s*$'),
    'table_rows': re.compile(r'^EL inlet velocity table\s*=.*\((\d+) rows'),
    'extents': re.compile(r'^Global CFD domain extents:\s*(.*)$'),
}
FINISH = 'PP3D_LES has successfully finished'


class Log:
    """everything the tool reads from run_slurm.log, from one pass over its lines"""

    def __init__(self, lines, last_complete=True):
        self.echo = {}
        self.table = []                  # (t, u)
        self.audit = []                  # (line number, record, raw text)
        self.budget = {}                 # step -> record (None: conflicting copies)
        self.pairs = {}                  # step -> record (None: conflicting copies)
        self.times = []                  # (t, t_max, step, nsteps, dt)
        self.finished = False
        bad = []                         # (line number, tag, reason)
        good_last = -1
        n = len(lines)
        for i, line in enumerate(lines):
            unfinished = (i == n - 1 and not last_complete)
            s = line.rstrip('\r\n')
            st = s.lstrip()
            if st.startswith('EL_INLET_TABLE'):
                tok = st.split()
                try:
                    if unfinished or len(tok) != 3:
                        raise ValueError
                    self.table.append((num(tok[1]), num(tok[2])))
                    good_last = i
                except ValueError:
                    bad.append((i, 'EL_INLET_TABLE', 'no line end' if unfinished else 'not two numbers t U'))
                continue
            tag = s.split(' ', 1)[0]
            if tag in SPECS:
                rec, why = (None, 'no line end') if unfinished else parse_keyed(s.split()[1:], SPECS[tag])
                if rec is None:
                    bad.append((i, tag, why))
                    continue
                good_last = i
                if tag == 'EL_BED_AUDIT':
                    self.audit.append((i, rec, ' '.join(s.split()[1:])))
                else:
                    store = self.budget if tag == 'EL_FORCE_BUDGET' else self.pairs
                    k = rec['step'][0]
                    if k not in store:
                        store[k] = rec
                    elif store[k] is not None and store[k] != rec:
                        store[k] = None
                continue
            if st.startswith('time:'):
                m = PAT_TIME.match(s)
                if m and not unfinished:
                    self.times.append((num(m.group(1)), num(m.group(2)), int(m.group(3)), int(m.group(4)),
                                       num(m.group(5))))
                    good_last = i
                else:
                    bad.append((i, 'time line', 'no line end' if unfinished else 'not readable'))
                continue
            if FINISH in s and not unfinished:
                self.finished = True
                good_last = i
                continue
            if s.startswith(('EL ', 'Global CFD')):
                for key, pat in ECHO.items():
                    if key not in self.echo:
                        m = pat.match(s)
                        if m:
                            self.echo[key] = m.group(1).strip()
        self.trailing = [b for b in bad if b[0] > good_last]
        self.malformed = [b for b in bad if b[0] <= good_last]


def read_deck(text):
    """SimPar keys of a deck with the application's syntax: 'SimPar@Key = value', the last '@' and the
    last '=' of the line, the first list-directed token of the value (blank, comma or '/' ends it, so a
    trailing '! comment' is not part of it), a later line overrides an earlier one. Values stay text."""
    keys = {}
    for line in text.splitlines():
        if '@' not in line or '=' not in line:
            continue
        iat, ieq = line.rfind('@'), line.rfind('=')
        if ieq < iat:
            continue
        sec = line[:iat].split()
        key = line[iat + 1:ieq].split()
        if not sec or sec[-1] != 'SimPar' or len(sec) != 1 or not key:
            continue
        val = re.split(r'[\s,/]+', line[ieq + 1:].strip())
        keys[key[0]] = val[0] if val and val[0] else ''
    return keys


def read_deck_num(keys, key):
    """(value, status): status 'deck', 'absent' or 'unreadable'"""
    if key not in keys:
        return NAN, 'absent'
    try:
        return num(keys[key].strip('"\'')), 'deck'
    except ValueError:
        return NAN, 'unreadable'


# ===================================================================================================
# B. protocol
# ===================================================================================================
def u_of_t(table, t):
    """the solver's rule: linear between the rows, constant beyond the ends"""
    if t <= table[0][0]:
        return table[0][1]
    for (t0, u0), (t1, u1) in zip(table, table[1:]):
        if t <= t1:
            return u0 + (u1 - u0) * (t - t0) / (t1 - t0)
    return table[-1][1]


def holds_of(table, t_start, t_end):
    """maximal constant-velocity intervals of U(t) from t_start to t_end (None: the run has no known end,
    the last hold is open): equal-valued adjacent segments merged first, then holds of HOLD_MIN or less
    dropped. Each hold: dict(u, a, b, open)."""
    if not table:
        return []
    segs = [(-INF, table[0][0], table[0][1], table[0][1])]
    segs += [(p[0], q[0], p[1], q[1]) for p, q in zip(table, table[1:])]
    segs.append((table[-1][0], INF, table[-1][1], table[-1][1]))
    merged = []
    for a, b, ua, ub in segs:
        if ua != ub:
            merged.append(None)                         # a ramp separates two holds
        elif merged and merged[-1] is not None and merged[-1][2] == ua and merged[-1][1] == a:
            merged[-1] = (merged[-1][0], b, ua)
        else:
            merged.append((a, b, ua))
    lo = t_start if finite(t_start) else table[0][0]
    out = []
    for h in merged:
        if h is None:
            continue
        a, b, u = max(h[0], lo), h[1], h[2]
        is_open = t_end is None and b == INF
        if t_end is not None:
            b = min(b, t_end)
        if is_open or b - a > HOLD_MIN:
            out.append(dict(u=u, a=a, b=b, open=is_open))
    return out


def steps_of(holds):
    """protocol steps of the holds: t1 = end of the hold, t0 = end of the previous hold (the start of the
    transition into it); the first hold is evaluated over its last FIRST_LEN s when longer than 0.6 s"""
    steps = []
    for i, h in enumerate(holds):
        t1 = h['b']
        t0 = holds[i - 1]['b'] if i > 0 else h['a']
        if i == 0 and not h['open'] and t1 - t0 > 0.6:
            t0 = t1 - FIRST_LEN
        steps.append(dict(k=i, t0=t0, t1=t1, U=h['u'], hold=h))
    return steps


# ===================================================================================================
# C. clock and audit sequence
# ===================================================================================================
class Clock:
    """time step, start, completed steps and audit frequency of the run, and the validated audit rows"""

    def __init__(self, log, deck):
        self.notes = []
        self.freq, self.freq_src = None, 'unknown'
        f_echo = log.echo.get('audit_freq')
        if f_echo is not None and re.fullmatch(r'\d+', f_echo) and int(f_echo) > 0:
            self.freq, self.freq_src = int(f_echo), 'log echo'
        else:
            v, st = read_deck_num(deck, 'ELBedAuditFreq')
            if st == 'deck' and finite(v) and v > 0 and v == int(v):
                self.freq, self.freq_src = int(v), 'deck'
        times = log.times
        dts = sorted({x[4] for x in times})
        self.const = len(dts) <= 1
        dt_deck, _ = read_deck_num(deck, 'TimeStep')
        dt_print = dts[0] if len(dts) == 1 else NAN
        recs = sorted(((r['step'][0], r['t'][0]) for _, r, _ in log.audit if finite(r['t'][0])))
        dt_audit = NAN
        if len(recs) > 1 and recs[-1][0] > recs[0][0]:
            dt_audit = (recs[-1][1] - recs[0][1]) / (recs[-1][0] - recs[0][0])

        def agree(a, b):
            return finite(a) and finite(b) and abs(a - b) <= TOL_D4 * max(abs(a), abs(b))
        self.dt, self.dt_src = NAN, 'unknown'
        if not self.const:
            self.notes.append('time step not constant (printed dt %s): no step clock, step times from the audit '
                              'rows only' % ', '.join('%.4g' % x for x in dts))
        elif finite(pos(dt_deck)) and (not times or agree(dt_deck, dt_print)):
            self.dt, self.dt_src = dt_deck, 'deck TimeStep' + (', time lines agree' if times else '')
        elif finite(pos(dt_audit)) and (not times or agree(dt_audit, dt_print)):
            self.dt, self.dt_src = dt_audit, 'audit rows' + (', time lines agree' if times else '')
        elif finite(pos(dt_print)):
            self.dt, self.dt_src = dt_print, 'time lines (4 digits)'
        self.exact = finite(self.dt)
        t_line = times[0][0] - times[0][2] * self.dt if times and self.exact else NAN
        t_deck, _ = read_deck_num(deck, 'StartSimTime')
        tol0 = TOL_D4 * abs(times[0][0]) + 1e-12 if times else 0.0
        self.t0, self.t0_src = NAN, 'unknown'
        if self.exact:
            if finite(t_deck) and (not times or abs(t_deck - t_line) <= tol0):
                self.t0, self.t0_src = t_deck, 'deck StartSimTime'
            elif recs and (not times or abs(recs[0][1] - recs[0][0] * self.dt - t_line) <= tol0):
                self.t0, self.t0_src = recs[0][1] - recs[0][0] * self.dt, 'first audit row'
            elif times:
                self.t0, self.t0_src = t_line, 'first time line (4 digits)'
            else:
                self.exact = False
        self.K_time = max((x[2] for x in times), default=0)
        self.finished = log.finished
        self.validate(log)

    def tk(self, k):
        return self.t0 + k * self.dt if self.exact else self.row_t.get(k, NAN)

    def validate(self, log):
        """validated rows (step -> row), rejected (step, reason), the completed steps K"""
        groups = {}
        for ln, rec, raw in log.audit:
            groups.setdefault(rec['step'][0], []).append((ln, rec, raw))
        self.rejected = []
        cand = {}
        for k in sorted(groups):
            g = groups[k]
            if len({raw for _, _, raw in g}) > 1:
                self.rejected.append((k, 'conflicting copies of the audit row'))
                continue
            rec = g[0][1]
            if self.freq and k % self.freq:
                self.rejected.append((k, 'step not a multiple of the audit frequency'))
                continue
            if k < 1:
                self.rejected.append((k, 'step below 1'))
                continue
            t = rec['t'][0]
            if self.exact:
                tk = self.t0 + k * self.dt
                if not (finite(t) and abs(t - tk) <= 2 * es_prec(tk) + 1e-3 * self.dt):
                    self.rejected.append((k, 'time %s does not match the step clock %.6f s' % (t, tk)))
                    continue
            elif not finite(t):
                self.rejected.append((k, 'time not finite'))
                continue
            cand[k] = row_of(rec)
        rows, prev = {}, None
        keys = sorted(cand)
        for j, k in enumerate(keys):
            r = cand[k]
            if not self.exact and prev is not None and not r['t'] > rows[prev]['t']:
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
        self.row_t = {k: r['t'] for k, r in rows.items()}
        for k, r in rows.items():
            r['tk'] = self.tk(k)
        self.sorted = [rows[k] for k in sorted(rows)]
        self.K = max([self.K_time] + list(rows))
        self.end = self.tk(self.K) if self.finished and self.K > 0 else None
        if self.end is not None and not finite(self.end):
            self.end = None
        self.rej = dict(self.rejected)

    def expected(self, a, b):
        """expected audit steps with their time in [a, b] (None when not established)"""
        if not (self.exact and self.freq) or not (finite(a) and finite(b)):
            return None
        tol = 1e-3 * self.dt
        k0 = max(1, math.ceil((a - tol - self.t0) / self.dt))
        k1 = math.floor((b + tol - self.t0) / self.dt)
        j0 = -(-k0 // self.freq)
        return [j * self.freq for j in range(j0, k1 // self.freq + 1)]

    def coverage(self, a, b, min_samples=MIN_SAMPLES):
        """(ok, reason): every expected audit step of [a, b] validated"""
        if not self.freq:
            return False, 'audit frequency unknown'
        if not self.exact:
            return False, 'no constant time step: the expected audit steps are unknown'
        ks = self.expected(a, b)
        if ks is None:
            return False, 'window not finite'
        if len(ks) < max(1, min_samples):
            return False, 'audit too sparse for %.4f-%.4f s (%d expected rows)' % (a, b, len(ks))
        late = [k for k in ks if k > self.K]
        if late:
            return False, 'the run has not reached %.4f s' % self.tk(late[0])
        miss = [k for k in ks if k not in self.rows]
        if miss:
            return False, '%d of %d expected audit rows missing or rejected (first at %.4f s: %s)' % (
                len(miss), len(ks), self.tk(miss[0]), self.rej.get(miss[0], 'not in the log'))
        return True, ''

    def window(self, a, b, open_end=False):
        """validated rows with their step time in [a, b] (b excluded with open_end), in step order"""
        tol = 1e-3 * self.dt if self.exact else 1e-9
        return [r for r in self.sorted if r['tk'] >= a - tol and (r['tk'] < b - tol if open_end else r['tk'] <= b + tol)]


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


# ===================================================================================================
# the evaluation (shared with fl11_movie.py)
# ===================================================================================================
def setup_of(log, deck):
    """mapping support, drag model, fluid gravity, extents"""
    s = dict(fluid_gravity=None, zmin=None, zmax=None, extents_prec=(0.0, 0.0))
    w = log.echo.get('width')
    if w is not None:
        try:
            s['width'], s['width_src'] = num(w), 'log echo'
        except ValueError:
            s['width'], s['width_src'] = NAN, 'log echo (unreadable)'
    else:
        v, st = read_deck_num(deck, 'ELKernelWidthFactor')
        if st == 'deck':
            s['width'], s['width_src'] = v, 'deck'
        elif st == 'unreadable':
            s['width'], s['width_src'] = NAN, 'deck (value unreadable)'
        else:
            s['width'], s['width_src'] = 2.5, 'EL_CONFIG default (key absent from echo and deck)'
    dm = log.echo.get('drag')
    if dm is not None:
        s['drag'], s['drag_src'] = dm.lower(), 'log echo'
    elif 'ELDragModel' in deck:
        s['drag'], s['drag_src'] = deck['ELDragModel'].strip('"\'').lower(), 'deck'
    else:
        s['drag'], s['drag_src'] = None, 'unknown'
    fg = log.echo.get('fluid_gravity')
    if fg in ('T', 'F'):
        s['fluid_gravity'] = (fg == 'T')
    ex = log.echo.get('extents')
    if ex is not None:
        try:
            v = [num(x) for x in ex.split()]
            if len(v) == 6 and finite(v[4]) and finite(v[5]) and v[5] > v[4]:
                s['zmin'], s['zmax'] = v[4], v[5]
                s['extents_prec'] = (es_prec(v[4]), es_prec(v[5]))
        except ValueError:
            pass
    return s


def evaluate(rundir=None, d=0.87e-3, steps_path=None, lines=None, deck_text=None, last_complete=True, name=None,
             flag=0.02):
    """read (A), protocol (B), clock (C), per step (D); the verdicts (E) are report()"""
    if lines is None:
        with open(os.path.join(rundir, 'run_slurm.log'), errors='replace', newline='') as fh:
            text = fh.read()
        last_complete = text.endswith('\n')
        lines = text.split('\n')
        if lines and lines[-1] == '':
            lines.pop()
    if deck_text is None:
        p = os.path.join(rundir, '_data', 'q2p1_param.dat') if rundir else ''
        deck_text = open(p, errors='replace').read() if p and os.path.exists(p) else ''
    log = Log(lines, last_complete)
    deck = read_deck(deck_text)
    clock = Clock(log, deck)
    setup = setup_of(log, deck)
    ev = dict(name=name or os.path.basename(os.path.abspath(rundir)), log=log, deck=deck, clock=clock, setup=setup,
              d=d, r=d / 2, flag=flag, steps_path=steps_path, listed=None, listed_errors=[], listed_findings=[],
              listed_ok=None, protocol_errors=[])
    zmin = setup['zmin'] if setup['zmin'] is not None else NAN
    for r in clock.sorted:
        fb = log.budget.get(r['step'])
        r['dp_force'] = div(fb['drag'][2] + fb['pressure'][2], r['a_in']) if fb else NAN
        r['H'] = r['z99'] + ev['r'] - zmin
        if not r['has_fbal']:
            pr = log.pairs.get(r['step'])
            r['fbal'] = (pr['dp'][2] - pr['expected'][2]) / clock.dt if pr and finite(pos(clock.dt)) else NAN
    table = log.table
    declared = log.echo.get('table_rows')
    if table:
        if declared is not None and int(declared) != len(table):
            ev['protocol_errors'].append('table echo incomplete: %d of %s rows readable' % (len(table), declared))
        if any(not q[0] > p[0] for p, q in zip(table, table[1:])) or not all(finite(t) and finite(u) for t, u in table):
            ev['protocol_errors'].append('table echo not strictly increasing and finite')
        if declared is None and any(b[1] == 'EL_INLET_TABLE' for b in log.malformed):
            ev['protocol_errors'].append('a table echo row is malformed')
    if steps_path:
        ev['listed'], ev['listed_errors'] = read_steps(steps_path)
    if table:
        ev['source'] = 'echo'
        ev['holds'] = holds_of(table, clock.t0, clock.end)
        ev['protocol'] = steps_of(ev['holds'])
        if ev['listed'] is not None:
            ev['listed_ok'], ev['listed_findings'] = compare_steps(ev['listed'], ev['protocol'])
            if ev['listed_errors']:
                ev['listed_ok'] = False
    elif ev['listed'] is not None:
        ev['source'] = 'steps file'
        rows = sorted(ev['listed'])
        errs = list(ev['listed_errors'])
        for p, q in zip(rows, rows[1:]):
            if q[0] < p[1] - 1e-9:
                errs.append('rows %s and %s overlap or repeat' % (_fs(p), _fs(q)))
        ev['protocol_errors'] += errs
        ev['holds'] = [dict(u=u, a=t0, b=t1, open=False) for t0, t1, u in rows]
        ev['protocol'] = [dict(k=i, t0=t0, t1=t1, U=u, hold=ev['holds'][i]) for i, (t0, t1, u) in enumerate(rows)]
    else:
        ev['source'] = None
        ev['holds'], ev['protocol'] = [], []
    ev['knots'] = table
    for r in clock.sorted:
        r['U'] = u_of_t(table, r['t']) if table and finite(r['t']) else NAN
    ev['steps'] = [eval_step(ev, s) for s in ev['protocol']]
    last = clock.sorted[-1] if clock.sorted else None
    ev['area'] = pos(last['a_in']) if last else NAN
    ev['W'] = pos(-last['wz']) if last else NAN
    return ev


def eval_step(ev, s, hold_check=True):
    """D: one step; status 'measured' or the reason it is unmeasurable"""
    clock, flag = ev['clock'], ev['flag']
    h = s['hold']
    x = dict(s)
    x.update(status='measured', why='')

    def fail(why):
        if x['status'] == 'measured':
            x['status'], x['why'] = 'unmeasurable', why
    if h['open'] or not finite(s['t1']):
        fail('the run has not finished: the end of this hold is unknown')
        x['t1'] = clock.tk(clock.K) if clock.K > 0 and finite(clock.tk(clock.K)) else h['a']
    t0, t1 = x['t0'], x['t1']
    a = t1 - EVAL_LEN
    if hold_check and a < h['a'] - 1e-9:
        fail('hold %.4f-%.4f s shorter than the evaluation window %.4f-%.4f s' % (h['a'], h['b'], a, t1))
    ok, why = clock.coverage(a, t1)
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
    x['dp'], x['std'] = st('dp_dyn', a, t1)
    m1, _ = st('dp_dyn', t0 + 0.2, t0 + 0.35, True)
    m2, _ = st('dp_dyn', t0 + 0.35, t0 + 0.5)
    x['mf'], _ = st('dp_force', a, t1)
    x['H'], _ = st('H', a, t1)
    x['eps'], _ = st('eps', a, t1)
    fl, _ = st('fbal', a, t1)
    flpe, _ = st('floor', a, t1)
    x.update(dps=div(x['dp'], wa), stds=div(x['std'], wa), m1s=div(m1, wa), m2s=div(m2, wa), mfs=div(x['mf'], wa),
             fl=div(fl, Wk), flpe=div(flpe, Wk), flag=finite(wa) and abs(m1 - m2) > flag * wa)
    x['measured'] = x['status'] == 'measured'
    return x


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
    """a steps file against the protocol, all three columns (time within the %.6g print precision of the
    steps file, U relative 1e-6) and multiplicity: (exact one-to-one match, findings)"""
    prot = [(p['t0'], p['t1'], p['U']) for p in protocol]

    def close(a, b, rel):
        return abs(a - b) <= rel * max(1.0, abs(b))

    def same(a, b):
        return close(a[0], b[0], 1e-5) and close(a[1], b[1], 1e-5) and close(a[2], b[2], 1e-6)
    hits = [[i for i, p in enumerate(prot) if same(s, p)] for s in listed]
    count = [sum(1 for h in hits if i in h) for i in range(len(prot))]
    findings = []
    for s, h in zip(listed, hits):
        if not h:
            near = [p for p in prot if close(s[1], p[1], 1e-5) and close(s[2], p[2], 1e-6)]
            findings.append('listed %s %s' % (_fs(s), 'differs from the plateau %s in the start time' % _fs(near[0])
                                                if near else 'is no plateau of the protocol'))
    for p, c in zip(prot, count):
        if c == 0:
            findings.append('plateau %s missing from the list' % _fs(p))
        elif c > 1:
            findings.append('plateau %s listed %d times' % (_fs(p), c))
    return not findings, findings


def _fs(s):
    return '(%.4f, %.4f, %.4f)' % tuple(s)


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
# E. report and verdicts
# ===================================================================================================
def print_step_table(res, title):
    print(f'\n{title}')
    print(' step  U [m/s]  window      dp*      std     sub1    sub2   flag  F/(W/A)  H [mm]  eps_bed  floor/W  dp*+floor  (floor_pe/W)')
    for x in res:
        print(f" {x['k']:3d}  {x['U']:6.3f}  {x['t1']-EVAL_LEN:5.2f}-{x['t1']:5.2f}  {x['dps']:7.4f}  {x['stds']:7.4f}  {x['m1s']:6.3f}  {x['m2s']:6.3f}  "
              f"{'DIFF' if x['flag'] else '    '}  {x['mfs']:7.4f}  {x['H']*1e3:6.2f}  {x['eps']:7.4f}  {x['fl']:7.4f}  {x['dps'] + x['fl']:8.4f}   ({x['flpe']:.4f})"
              f"{'' if x['measured'] else '  (unmeasurable: ' + x['why'] + ')'}")
    print('  dp* = dp_dyn A / W over 0.2-0.5 s of the step (W, A of that step); floor from the particle momentum '
          'balance; floor_pe = the pe wall-contact impulse (unrelaxed PGS impulse, 1/0.9 of the applied one)')


def report(ev, a):
    """prints the evaluation, returns the verdicts {gate: PASS | FAIL | unavailable}"""
    log, clock, setup = ev['log'], ev['clock'], ev['setup']
    rows = clock.sorted
    label = '' if ev['source'] == 'echo' else FROM_STEPS
    V = {g: 'unavailable' for g in GATES}
    for b in log.trailing:
        print(f"note: unfinished trailing {b[1]} record at line {b[0] + 1} ignored ({b[2]})", file=sys.stderr)
    if rows:
        last = rows[-1]
        print(f"{ev['name']}: {len(rows)} audit lines, t {rows[0]['t']:.4f} .. {last['t']:.4f} s, "
              f"N {rows[0]['n']}, A {ev['area']:.6e} m2, buoyant weight {ev['W']:.6e} N, W/A {div(ev['W'], ev['area']):.3f} Pa, "
              f"hydrostatic part in dp {last['dp_h']:.4f} Pa (subtracted: dp_dyn)")
    else:
        print(f"{ev['name']}: no validated EL_BED_AUDIT line (SimPar@ELBedAuditFreq not set, or every row rejected)")
    if log.malformed:
        ex = '; '.join(f"{t} line {i + 1}: {w}" for i, t, w in log.malformed[:3])
        print(f"  note: {len(log.malformed)} malformed record(s) ignored ({ex}{' ...' if len(log.malformed) > 3 else ''})")
    if setup['zmin'] is None:
        print('  column extents not in the log (Global CFD domain extents): bed height, F3 and F5 height unavailable')
    ks = clock.expected(clock.t0, clock.tk(clock.K)) if clock.K > 0 else None
    print(f"clock: dt {clock.dt:.6g} s [{clock.dt_src}], start {clock.t0:.4f} s [{clock.t0_src}], {clock.K} steps completed "
          f"({'run finished' if clock.finished else 'no finish marker: the run has not finished'}), audit every "
          f"{clock.freq if clock.freq else '?'} steps [{clock.freq_src}]: {len(ks) if ks is not None else '?'} expected rows, "
          f"{len(rows)} validated, {len(clock.rejected)} rejected")
    for n_ in clock.notes:
        print('  note: ' + n_)
    if clock.rejected:
        print('  rejected audit rows: ' + '; '.join(f"step {k}: {w}" for k, w in clock.rejected[:3])
              + (' ...' if len(clock.rejected) > 3 else ''))
    if a.series:
        with open(a.series, 'w') as fh:
            fh.write('t,step,U,dp,dp_hydro,dp_dyn,dp_force,floor_z,z50,z95,z99,zmax,eps_bed,n\n')
            for r in rows:
                fh.write(','.join(str(r[k]) for k in ('t', 'step', 'U', 'dp', 'dp_h', 'dp_dyn', 'dp_force', 'floor',
                                                      'z50', 'z95', 'z99', 'zmax', 'eps', 'n')) + '\n')
    if ev['source'] == 'echo':
        tend = 'open (no finish marker)' if clock.end is None else f'{clock.end:.4f} s'
        print(f"protocol: {len(ev['protocol'])} holds of the echoed table (maximal constant-velocity intervals of the "
              f"solver's U(t)), terminal hold to the end of the run {tend}")
    else:
        print(f"protocol: {len(ev['protocol'])} steps of the steps file {a.steps}{FROM_STEPS}")
    for e in ev['protocol_errors']:
        print('  protocol NOT usable: ' + e)
    print_step_table(ev['steps'], 'steps of the protocol' + label)
    if ev['source'] == 'echo' and ev['listed'] is not None:
        if ev['listed_ok']:
            print(f"steps file {a.steps}: identical to the echoed protocol (all three columns, each plateau once)")
        else:
            print(f"steps file {a.steps}: NOT the echoed protocol: " + '; '.join(ev['listed_errors'] + ev['listed_findings']))
            listed = [dict(k=i, t0=t0, t1=t1, U=u, hold=dict(u=u, a=t0, b=t1, open=False))
                      for i, (t0, t1, u) in enumerate(ev['listed'])]
            print_step_table([eval_step(ev, s, hold_check=False) for s in listed],
                             'steps-file summary, NOT gated (the gates below use the echoed protocol)')
    if ev['protocol_errors']:
        print(f"\nF1, F2, F3, F5 height unavailable (protocol not usable: {'; '.join(ev['protocol_errors'])}){label}")
    else:
        V.update(gates_f1_f3(a, ev, label))
        V['F5 height'] = gate_f5_height(a, ev, label)
    V['F5 count'] = gate_f5_count(ev, label)
    print('verdicts: ' + ' | '.join(f"{g} {V[g]}" for g in GATES) + label)
    return V


def gates_f1_f3(a, ev, label):
    res = ev['steps']
    V = dict(F1='unavailable', F2='unavailable', F3='unavailable')
    bad = [x for x in res if not x['measured']]
    if len(res) < 3:
        print(f'\nF1, F2, F3 unavailable: {len(res)} protocol step(s), three needed' + label)
        return V
    if bad:
        why = '; '.join(f"step {x['k']} (U {x['U']:.3f}): {x['why']}" for x in bad[:3]) + (' ...' if len(bad) > 3 else '')
        flu = [x for x in res if x['measured'] and x['dps'] > 0.95]
        sub = (f"mean dp* over the {len(flu)} measured steps with dp* > 0.95: {mean([x['dps'] for x in flu]):.6f}"
               if flu else 'no measured step with dp* > 0.95')
        print(f"\nF1, F2, F3 unavailable ({len(bad)} protocol step(s) unmeasurable, each could belong to the set of a "
              f"gate: {why}){label}")
        print(f"   measured, not a verdict: {sub}")
        return V
    # ---- F1 (spec rule of 2026-11-20): mean dp* over all steps with dp* > 0.95 ------------------
    flu = [x for x in res if x['dps'] > 0.95]
    if flu:
        plateau = sum(x['dps'] for x in flu) / len(flu)
        V['F1'] = 'PASS' if abs(plateau - 1) <= 0.02 else 'FAIL'
        print(f"\nF1 plateau: mean dp* over {len(flu)} steps with dp* > 0.95 (U {min(x['U'] for x in flu):.2f} .. "
              f"{max(x['U'] for x in flu):.2f}): {plateau:.6f}   (gate 1 +- 0.02: {V['F1']}){label}; "
              f"dp* + floor/W over these {sum(x['dps'] + x['fl'] for x in flu)/len(flu):.4f}")
    else:
        plateau = 1.0
        print('\nF1 unavailable: no step with dp* > 0.95; plateau taken as W/A for the fixed branch' + label)
    fixed = [x for x in res if x['dps'] < 0.8 * plateau]
    if len(fixed) < 2:
        print(f'F2, F3 unavailable: {len(fixed)} fixed-branch step(s) (dp* < 0.8 of the plateau), two needed' + label)
        return V
    vp = math.pi * a.d ** 3 / 6
    clo_name = ev['setup']['drag']
    official = clo_name if clo_name in CLOSURES else None
    runs = [(official, True)] if official else []
    if not official:
        print(f"F2, F3 unavailable: the drag model of the run [{ev['setup']['drag_src']}] is "
              f"{clo_name if clo_name else 'unknown'}, not one of {', '.join(CLOSURES)}{label}")
    if a.closure and a.closure != official:
        runs.append((a.closure, False))
    for clo, gated in runs:
        tag = label if gated else f' [comparison with --closure {clo}, not gated; the run used {clo_name}]'
        if not gated:
            print(f"\n--- comparison with --closure {clo} (not gated) ---")
        v3 = f3_block(a, fixed, plateau, clo, vp, tag, gated)
        v2 = gate_f2(a, ev, fixed, clo, vp, tag, gated)
        if gated:
            V['F3'], V['F2'] = v3, v2
    return V


def f3_block(a, fixed, plateau, clo, vp, tag, gated):
    if not all(finite(x['H']) and finite(pos(x['A'])) and finite(x['n']) for x in fixed):
        print(f"F3 unavailable: bed height, area or particle count not finite on every fixed step{tag}")
        return 'unavailable'
    efix = mean([x['eps'] for x in fixed])
    uc = closure_umf(efix, a.d, a.rhop, a.rhof, a.mu, a.g, clo)
    ebed = sum(1 - div(x['n'] * vp, x['A'] * x['H']) for x in fixed) / len(fixed)
    ub = closure_umf(ebed, a.d, a.rhop, a.rhof, a.mu, a.g, clo)
    su2 = sum(x['U']**2 for x in fixed); su3 = sum(x['U']**3 for x in fixed); su4 = sum(x['U']**4 for x in fixed)
    sy1 = sum(x['U']*x['dps'] for x in fixed); sy2 = sum(x['U']**2*x['dps'] for x in fixed)
    det = su2*su4 - su3*su3
    ca = (sy1*su4 - sy2*su3)/det if det != 0 else NAN
    cb = (su2*sy2 - su3*sy1)/det if det != 0 else NAN

    def cross(level):
        # smallest positive root of b U^2 + a U = level (the cancellation-free pair of roots q / b and
        # -level / q); nan when the fit has no positive crossing
        if not (finite(ca) and finite(cb)): return NAN
        if cb == 0: return level/ca if ca > 0 else NAN
        disc = ca*ca + 4*cb*level
        if disc < 0: return NAN
        q = -0.5 * (ca + math.copysign(math.sqrt(disc), ca))
        roots = [x for x in ((q / cb) if q != 0 else NAN, (-level / q) if q != 0 else NAN) if finite(x) and x > 0]
        return min(roots) if roots else NAN
    u1 = cross(1.0); up = cross(plateau)
    r3 = div(u1, ub)
    v = ('PASS' if abs(r3 - 1) <= 0.05 else 'FAIL') if finite(r3) else 'unavailable'
    print(f"F3 fixed branch {len(fixed)} steps (U {min(x['U'] for x in fixed):.3f} .. {max(x['U'] for x in fixed):.3f}): "
          f"dp* = {ca:.4f} U + {cb:.4f} U^2 (through the origin)")
    print(f"   U_mf at dp* = 1: {u1:.4f} m/s; {clo} U_mf at the whole-bed porosity {ebed:.4f}: {ub:.4f} m/s; "
          f"ratio {r3:.4f}   (gate 0.05: {v if gated else 'comparison, not gated'}){tag}")
    print(f"   reported, not gated: U_mf at the measured plateau {plateau:.4f}: {up:.4f} m/s = {div(up, ub):.4f} of the closure; "
          f"at the core porosity {efix:.4f} the closure gives {uc:.4f} m/s: ratios {div(u1, uc):.4f} / {div(up, uc):.4f}")
    if gated:
        print(f"F4 U_mf / 0.182 (NETL measurement): {u1/0.182:.4f} (dp* = 1) / {up/0.182:.4f} (plateau) (reported, not gated)")
        s, i = linfit([x['U'] for x in fixed], [x['dp'] for x in fixed])
        print(f"   secondary (rule before 2026-11-20): linear fit of the fixed branch, slope {s:.1f} Pa s/m (intercept {i:.2f} Pa)")
    return v if gated else 'unavailable'


def valid_eps(e):
    """a mapped porosity the closure can be evaluated at: finite and 0 < eps < 1"""
    return finite(e) and 0 < e < 1


def gate_f2(a, ev, fixed, clo, vp, tag, gated):
    setup, r = ev['setup'], ev['r']
    reasons = []
    rows_all = [row for x in fixed for row in x['rows']]
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
    width, zmin = setup['width'], setup['zmin']
    if not finite(pos(width)):
        required = NAN
        reasons.append(f"mapping support not known as a positive finite number (kernel width factor {width} "
                       f"[{setup['width_src']}])")
    else:
        required = max(2 * a.d, width * a.d)
    lows, highs, fields, short, undecided = [], [], [], 0, 0
    for row in rows_all:
        lo = row['zlo'] - zmin if zmin is not None else NAN
        hi = row['z99'] + r - row['zhi']
        fm = row['margin'] if row['margin'] is not None else 2 * a.d
        lows.append(lo); highs.append(hi); fields.append(fm)
        if not (all(finite(x) for x in (lo, hi, fm)) and finite(required)):
            continue
        for m, tol in ((lo, es_prec(row['zlo']) + setup['extents_prec'][0]),
                       (hi, es_prec(row['z99']) + es_prec(row['zhi'])), (fm, es_prec(fm))):
            if tol > 0.05 * required:
                undecided += 1
            elif m < required - tol:
                short += 1
    if zmin is None:
        margin_ok = False
        reasons.append('missing geometry: no inlet-plane elevation in the log (Global CFD domain extents), '
                       'the lower window margin cannot be measured')
    elif not rows_all or not all(finite(x) for x in lows + highs + fields):
        margin_ok = False
        reasons.append('window bounds not finite on every contributing audit line')
    elif not finite(required):
        margin_ok = False
    elif undecided:
        margin_ok = False
        reasons.append('window margin undecidable: the coordinates are printed too coarsely for the margin')
    else:
        margin_ok = not short
        if short:
            reasons.append('window inside the mapping support')
    fin = [x for x in lows + highs + fields if finite(x)]
    print(f" F2: gradient {gtext}")
    print(f"     audit window margin over {len(rows_all)} lines: min {div(min(fin), a.d) if fin else NAN:.3f} d "
          f"(lower end - inlet plane min {_fmin(lows)/a.d:.3f} d, z99 + r - upper end min {_fmin(highs)/a.d:.3f} d, "
          f"{'win_margin field' if all(row['margin'] is not None for row in rows_all) else 'margin 2 d (binary before 627caf47)'} "
          f"min {_fmin(fields)/a.d:.3f} d); required max(2 d, mapping support {width:g} d "
          f"[{setup['width_src']}]) = {div(required, a.d):.3f} d: {'met' if margin_ok else 'NOT met'}")
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
        clo_b = closure_grad(x['U'], x['eps'], a.d, a.rhof, a.mu, clo)
        eh = 1 - div(x['n'] * vp, x['A'] * x['H'])
        cloh = closure_grad(x['U'], eh, a.d, a.rhof, a.mu, clo)
        clok = closure_grad(x['U'], ek, a.d, a.rhof, a.mu, clo) if valid_eps(ek) else NAN
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
    if not reasons and gated:
        v = 'PASS' if worst_k <= 0.05 else 'FAIL'
        print(f"F2 worst |ratio_k - 1| (window gradient against the closure at the mapped gas fraction) {worst_k:.4f}   "
              f"(gate 0.05: {v}){tag}")
    else:
        v = 'unavailable'
        print(f"F2 measured, not gated: worst |ratio_k - 1| {worst_k:.4f}")
        if reasons:
            print(f"F2 verdict unavailable ({'; '.join(reasons)}){tag}")
    if finite(worst_w):
        print(f"   reported, not gated: worst |ratio_w - 1| at the geometric window porosity (rule before 2026-11-20) {worst_w:.4f}")
    return v


def gate_f5_count(ev, label):
    clock = ev['clock']
    ns = sorted({r['n'] for r in clock.sorted})
    if len(ns) > 1:
        print(f"F5 count: particle counts {ns} among the validated audit rows: FAIL{label}")
        return 'FAIL'
    if not clock.finished:
        why = 'the run has not finished (no finish marker)'
    elif not clock.sorted:
        why = 'no validated audit row'
    elif ns == [0]:
        why = 'no particle'
    else:
        ok, why = clock.coverage(clock.t0, clock.tk(clock.K), min_samples=1)
        if ok:
            print(f"F5 count: counts seen {ns} on {len(clock.sorted)} audit lines, run {clock.t0:.4f}-{clock.tk(clock.K):.4f} s "
                  f"covered: PASS{label}")
            return 'PASS'
    print(f"F5 count: counts seen {ns}; count over the run not established ({why}): unavailable{label}")
    return 'unavailable'


def gate_f5_height(a, ev, label):
    clock, setup, holds = ev['clock'], ev['setup'], ev['holds']
    if not holds:
        print(f"F5 height unavailable: the protocol holds no velocity plateau{label}")
        return 'unavailable'
    u = min(h['u'] for h in holds)
    low = [h for h in holds if h['u'] == u]
    span = ', '.join(f"{h['a']:.3f}-{h['b']:.3f} s" if not h['open'] else f"{h['a']:.3f} s - open" for h in low)
    reasons, rows = [], []
    for h in low:
        b = h['b'] if not h['open'] else clock.tk(clock.K)
        if finite(b):
            rows += clock.window(h['a'], b)
        if h['open']:
            reasons.append('the lowest hold is the terminal hold of a run that has not finished')
        else:
            ok, why = clock.coverage(h['a'], h['b'])
            if not ok:
                reasons.append(why)
    if setup['zmin'] is None or setup['zmax'] is None:
        print(f"F5 height unavailable: no column extents in the log (Global CFD domain extents): mid-column unknown{label}")
        return 'unavailable'
    zmin, lz = setup['zmin'], setup['zmax'] - setup['zmin']
    if a.lz is not None and abs(a.lz - lz) > 1e-9 * max(1.0, lz):
        print(f"   note: --lz {a.lz:g} differs from the column height {lz:g} m of the echoed extents (the extents are used)")
    mid = zmin + lz / 2
    pmid = 0.5 * (setup['extents_prec'][0] + setup['extents_prec'][1])
    zs = [r['zmax'] for r in rows]
    fin = [z for z in zs if finite(z)]
    if len(fin) < len(zs):
        reasons.append(f'{len(zs) - len(fin)} audit line(s) without a centre height')
    above = [z for z in fin if z > mid + pmid + es_prec(z)]
    zm = max(fin) if fin else NAN
    if fin:
        print(f"F5 (b) highest centre during the plateau of the lowest velocity {u:.3f} m/s ({span}, {len(rows)} audit lines): "
              f"{(zm - zmin)*1e3:.2f} mm above the inlet plane (mid-column {lz*500:.1f} mm)")
    if above:
        v = 'FAIL'
        print(f"F5 height: a centre above mid-column observed: FAIL{label}")
    elif reasons or not fin:
        v = 'unavailable'
        print(f"F5 height unavailable ({'; '.join(reasons) if reasons else 'no observation'}){label}")
    else:
        v = 'PASS'
        print(f"F5 height: no centre above mid-column (a centre at mid-height is not above): PASS{label}")
    if ev['source'] == 'echo' and fin:
        trans = []
        for h in low:
            prev = [g for g in holds if g['b'] <= h['a'] + 1e-9]
            if prev and prev[-1]['b'] < h['a']:
                trans.append((prev[-1]['b'], h['a']))
        if trans:
            zt = _fmax([r['zmax'] for ta, tb in trans for r in clock.window(ta, tb)])
            print(f"   stricter, reported: including the {1e3*sum(tb - ta for ta, tb in trans):.0f} ms transition into it "
                  f"({', '.join(f'{ta:.3f}-{tb:.3f} s' for ta, tb in trans)}): highest centre "
                  f"{(_fmax([zm, zt]) - zmin)*1e3:.2f} mm")
    return v


# ===================================================================================================
# command line
# ===================================================================================================
def positive(name):
    def conv(s):
        try:
            v = float(s)
        except ValueError:
            raise argparse.ArgumentTypeError(f'{name} must be a number, got {s!r}')
        if not (math.isfinite(v) and v > 0):
            raise argparse.ArgumentTypeError(f'{name} must be finite and positive, got {s!r}')
        return v
    return conv


def nonneg(name):
    def conv(s):
        try:
            v = float(s)
        except ValueError:
            raise argparse.ArgumentTypeError(f'{name} must be a number, got {s!r}')
        if not (math.isfinite(v) and v >= 0):
            raise argparse.ArgumentTypeError(f'{name} must be finite and not negative, got {s!r}')
        return v
    return conv


def parser():
    ap = argparse.ArgumentParser(description='case 11 evaluation of a q2p1_el_pipeflow run')
    ap.add_argument('rundir', nargs='?')
    ap.add_argument('--steps', default=None)
    ap.add_argument('--d', type=positive('--d'), default=0.87e-3)
    ap.add_argument('--rhop', type=positive('--rhop'), default=860.0)
    ap.add_argument('--rhof', type=positive('--rhof'), default=1.0)
    ap.add_argument('--mu', type=positive('--mu'), default=1.8e-5)
    ap.add_argument('--g', type=positive('--g'), default=9.81)
    ap.add_argument('--lz', type=positive('--lz'), default=None)
    ap.add_argument('--closure', default=None, choices=CLOSURES)
    ap.add_argument('--flag', type=nonneg('--flag'), default=0.02)
    ap.add_argument('--series', default=None)
    ap.add_argument('--self-test', action='store_true')
    ap.add_argument('--self-test-real', default=None, metavar='RUNDIR')
    ap.add_argument('--variants', type=int, default=300)
    ap.add_argument('--seed', type=int, default=20261122)
    return ap


def main(argv=None):
    ap = parser()
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    if a.self_test_real:
        return self_test_real(a.self_test_real, a.variants, a.seed)
    if not a.rundir:
        ap.error('RUNDIR is required')
    ev = evaluate(a.rundir, a.d, a.steps, flag=a.flag)
    if ev['source'] is None:
        sys.exit('no EL_INLET_TABLE echo (SimPar@ELInletVelocityTable not set?) and no --steps file'
                 + (': ' + '; '.join(ev['listed_errors']) if ev['listed_errors'] else ''))
    report(ev, a)
    return 0


def run_lines(lines, deck_text, argv=(), last_complete=True, name='run'):
    """the report of an in-memory log (the metamorphic tests): (stdout, verdicts)"""
    a = parser().parse_args(list(argv) + ['x'])
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        ev = evaluate(None, a.d, a.steps, lines=lines, deck_text=deck_text, last_complete=last_complete, name=name,
                      flag=a.flag)
        V = report(ev, a) if ev['source'] else {g: 'unavailable' for g in GATES}
    return out.getvalue(), V


# ===================================================================================================
# self-test: synthetic logs
# ===================================================================================================
def _e(x, w=14, p=6):
    return ('%*s' % (w, 'NaN')) if x != x else '%*.*E' % (w, p, x)


def _d(x):
    return ('%12.4E' % x).replace('E', 'D')


PROTO = [(0.5, 0.30), (1.0, 0.20), (1.5, 0.10), (2.0, 0.05)]          # (start of step, U); end 2.5 s
FULL_STEPS = [(0.5, 1.0, 0.30), (1.0, 1.5, 0.20), (1.5, 2.0, 0.10), (2.0, 2.5, 0.05)]
DT = 2.5e-4


def _dps(u):
    """dp* of the synthetic bed: fixed below 0.18 m/s, a plateau that falls slightly with U above"""
    return min(1.0, u / 0.18) * ((1.0 - 0.1 * (u - 0.2)) if u > 0.18 else 1.0)


def _table(protocol, tstep=0.5):
    tab = [(0.0, protocol[0][1])]
    for (t, u) in protocol[1:]:
        tab += [(t, tab[-1][1]), (t + 0.02, u)]
    t_proto = protocol[-1][0] + tstep
    tab.append((t_proto, tab[-1][1]))
    return tab, t_proto


def _write_run(path, protocol=PROTO, cad=0.01, t_end=None, area=2.25e-4, n=9789, width=2.5, margin_d=2.5, zlo_d=None,
               dyn=True, dp_hydro=0.0, d=0.87e-3, steps_file=None, ekern=0.37, ekern_fn=None, extents=True, z0=0.0,
               table=None, echo_table=True, zmax_fn=None, legacy=False, finished=True, time_lines=True, deck=None,
               drag='gidaspow', width_echo=True, mutate=None, tmax=None, n_fn=None, grad_fn=None, last_complete=True):
    """a synthetic run directory: start-up echo, table echo, time lines, EL_FORCE_BUDGET and EL_BED_AUDIT lines.
    Heights above the inlet plane z0: bed top 24 mm (zmax_fn(t): highest centre), mapped porosity ekern;
    legacy = a binary before c8541cfc (no floor_bal_z, eps_kern, win_margin, gradz_win_dyn); mutate(lines)
    edits the lines before they are written"""
    os.makedirs(os.path.join(path, '_data'), exist_ok=True)
    rhop, rhof, g = 860.0, 1.0, 9.81
    vp = math.pi * d ** 3 / 6
    tab, t_proto = _table(protocol)
    if table is not None:
        tab = table
    t_end = t_proto if t_end is None else t_end
    freq = int(round(cad / DT))
    nsteps = int(round(t_end / DT))
    lines = ['EL kernel                 = deen_poly']
    if width_echo:
        lines.append('EL kernel width factor    = ' + (_e(width) if not isinstance(width, str) else width))
    if drag:
        lines.append('EL drag model             = ' + drag)
    lines.append('EL fluid gravity          = ' + ('T' if dp_hydro else 'F'))
    if echo_table:
        lines.append('EL inlet velocity table   = _data/inlet_table.txt (%d rows t / U; linear, constant beyond the ends)' % len(tab))
        lines += ['  EL_INLET_TABLE %s %s' % (_e(t, 16, 8), _e(u, 16, 8)) for t, u in tab]
    lines.append('EL bed audit freq         = %d' % freq)
    if extents:
        lines.append('Global CFD domain extents: ' + ' '.join(_e(x) for x in (0.0, 0.015, 0.0, 0.015, z0, z0 + 0.1)))
    tmx = t_end if tmax is None else tmax
    if time_lines:
        lines.append('time:%s/%s | itns:%8d/%8d | dt:%s%s%s' % (_d(DT), _d(tmx), 1, nsteps, _d(DT), _d(DT), _d(DT)))
    for k in range(freq, nsteps + 1, freq):
        t = k * DT
        u = u_of_t(tab, t)
        dps = _dps(u)
        nn = n_fn(t) if n_fn else n
        W = nn * (rhop - rhof) * vp * g
        dp = dps * div(W, area) + dp_hydro if (area > 0 and nn > 0) else NAN
        z99 = z0 + 0.024 if nn > 0 else NAN
        m = margin_d * d
        zlo = z0 + (zlo_d if zlo_d is not None else margin_d) * d if nn > 0 else NAN
        zhi = z99 + d / 2 - m if nn > 0 else NAN
        zmx = (z0 + (zmax_fn(t) if zmax_fn else 0.0241)) if nn > 0 else NAN
        grad = (grad_fn(t, dp) if grad_fn else dp / 0.1) if dp == dp else NAN
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
        lines.append('EL_FORCE_BUDGET step= %d drag= %s %s %s pressure= %s %s %s lift= 0 0 0 grav_buoy= %s %s %s' % (
            k, _e(0.0), _e(0.0), _e(dps * W * 0.9), _e(0.0), _e(0.0), _e(dps * W * 0.1), _e(0.0), _e(0.0), _e(-W)))
        lines.append(line)
        if time_lines:
            lines.append('time:%s/%s | itns:%8d/%8d | dt:%s%s%s' % (_d(t), _d(tmx), k, nsteps, _d(DT), _d(DT), _d(DT)))
    if finished:
        lines.append(' PP3D_LES has successfully finished. ')
    if mutate:
        lines = mutate(lines)
    with open(os.path.join(path, 'run_slurm.log'), 'w') as fh:
        fh.write('\n'.join(lines) + ('\n' if last_complete else ''))
    if deck is not None:
        with open(os.path.join(path, '_data', 'q2p1_param.dat'), 'w') as fh:
            fh.write(deck)
    if steps_file is not None:
        with open(os.path.join(path, '_data', 'steps.txt'), 'w') as fh:
            fh.write('# t_start t_end U\n' + ''.join('%.6g %.6g %.6g\n' % s for s in steps_file))
    return path


def _audit_at(lines, t):
    """index of the EL_BED_AUDIT line at time t"""
    for i, l in enumerate(lines):
        if l.startswith('EL_BED_AUDIT') and abs(float(l.split()[2]) - t) < 1e-9:
            return i
    raise KeyError(t)


def _drop(*ts):
    return lambda L: [l for i, l in enumerate(L) if i not in {_audit_at(L, t) for t in ts}]


def _nan_field(t, key, idx=0):
    def m(L):
        i = _audit_at(L, t)
        tok = L[i].split()
        tok[tok.index(key + '=') + 1 + idx] = 'NaN'
        L = list(L)
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
        j = tok.index('dp_hydro_dyn=') + 3
        tok[j] = _e(num(tok[j]) * 1.001, 16, 8)
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


def _reload_extension(t_first_max, until):
    """a Reload_SimPar extension: the time lines up to 'until' show the earlier MaxSimTime"""
    def m(L):
        out = []
        for l in L:
            mm = PAT_TIME.match(l)
            if mm and num(mm.group(1)) <= until + 1e-9:
                l = l.replace('/' + _d(num(mm.group(2))), '/' + _d(t_first_max), 1)
            out.append(l)
        return out
    return m


def self_test_cases():
    """(name, generator keywords, extra arguments, expected verdicts, strings that must appear, strings
    that must not, the answer of the tool before this round where it differed)"""
    S = '{run}/_data/steps.txt'
    hi_after = lambda t: 0.060 if t > 2.3 else 0.0241       # noqa: E731
    REF = {'F1': 'PASS', 'F2': 'FAIL', 'F3': 'PASS', 'F5 height': 'PASS', 'F5 count': 'PASS'}
    U = 'unavailable'
    base_tab = _table(PROTO)[0]
    tab_sub = base_tab[:-1] + [(2.1, 0.05), (2.18, 0.05), (2.26, 0.05), (2.34, 0.05), (2.42, 0.05), (2.5, 0.05)]
    tab_knot = sorted(base_tab + [(1.3, 0.20)])
    deck_full = 'SimPar@TimeStep = 0.00025d0\nSimPar@StartSimTime = 0d0\nSimPar@MaxNumStep = 10000\n'
    return [
        ('reference (all gates measurable)', dict(), [], REF, ['stricter, reported'], ['Traceback'], ''),
        ('missing inlet plane (area_in = 0)', dict(area=0.0), [], {'F1': U, 'F2': U, 'F3': U, 'F5 height': 'PASS'},
         [], ['Traceback'], ''),
        ('empty bed (n = 0)', dict(n=0), [], {'F1': U, 'F5 height': U, 'F5 count': U}, [], ['Traceback'], ''),
        ('r2: steps file without the lowest step', dict(steps_file=FULL_STEPS[:3]), ['--steps', S], REF,
         ['plateau (2.0000, 2.5000, 0.0500) missing from the list', 'steps-file summary, NOT gated'], [], ''),
        ('r2: foreign step in the steps file', dict(steps_file=FULL_STEPS + [(2.5, 3.0, 0.01)]), ['--steps', S], REF,
         ['is no plateau of the protocol'], [], ''),
        ('r3: steps file with a duplicated row', dict(steps_file=FULL_STEPS + [FULL_STEPS[0]]), ['--steps', S], REF,
         ['listed 2 times'], [], ''),
        ('r3: steps file with an omitted row', dict(steps_file=FULL_STEPS[1:]), ['--steps', S], REF, ['missing from the list'], [], ''),
        ('r3: steps file with a shifted start', dict(steps_file=[FULL_STEPS[0], (1.2, 1.5, 0.20)] + FULL_STEPS[2:]),
         ['--steps', S], REF, ['differs from the plateau (1.0000, 1.5000, 0.2000) in the start time'], [], ''),
        ('steps file identical to the protocol', dict(steps_file=FULL_STEPS), ['--steps', S], REF,
         ['identical to the echoed protocol'], ['NOT gated'], ''),
        ('short run (ends inside the lowest hold, no finish marker)', dict(t_end=2.3, finished=False), [],
         {'F1': U, 'F2': U, 'F3': U, 'F5 height': U, 'F5 count': U}, ['has not finished'], [], ''),
        ('sparse audit (cadence 0.4 s)', dict(cad=0.4), [], {'F1': U, 'F5 height': U}, ['audit too sparse'], [], ''),
        ('window margin 2 d with a 2.5 d support', dict(margin_d=2.0), [], {'F2': U, 'F1': 'PASS'},
         ['window inside the mapping support'], [], ''),
        ('lower window end at 2 d, field says 2.5 d', dict(zlo_d=2.0), [], {'F2': U}, ['window inside the mapping support'], [], ''),
        ('old log: no win_margin / gradz_win_dyn, zero hydrostatic part', dict(margin_d=2.0, dyn=False), [], {'F2': U},
         ['hydrostatic part zero', 'margin 2 d (binary before 627caf47)'], [], ''),
        ('old log with a 2 d support: margin met, total gradient admitted', dict(margin_d=2.0, dyn=False, width=2.0), [],
         {'F2': 'FAIL'}, ['hydrostatic part zero'], [], ''),
        ('old log with a hydrostatic part', dict(margin_d=2.0, dyn=False, width=2.0, dp_hydro=0.981), [], {'F2': U},
         ['nonzero hydrostatic part'], [], ''),
        ('r3: mapped porosity 1.0 (closure zero)', dict(ekern=1.0), [], {'F2': U}, ['no positive closure'], [], ''),
        ('r3: the -1 sentinel in the mapped samples', dict(ekern_fn=lambda t: -1.0 if 2.3 < t < 2.4 else 0.37), [], {'F2': U},
         ['mapped gas-fraction sample(s)'], [], ''),
        ('r3: zero-velocity fixed plateau', dict(protocol=PROTO + [(2.5, 0.0)]), [], {'F2': U}, ['step 4 (U 0.000)'], [], ''),
        ('r3: terminal hold without its last knot, 60 mm late in the hold', dict(table=base_tab[:-1], zmax_fn=hi_after), [],
         {'F5 height': 'FAIL', 'F1': 'PASS'}, ['lowest velocity 0.050 m/s (2.020-2.500 s'], [], ''),
        ('r3: the same without time lines and finish marker (end of the run unknown)',
         dict(table=base_tab[:-1], time_lines=False, finished=False), [],
         {'F5 height': U, 'F5 count': U}, ['terminal hold of a run that has not finished'], [], ''),
        ('r3: one audit line of a binary without floor_bal_z', dict(legacy=True, cad=2.5), [],
         {'F1': U, 'F5 height': U}, [], ['Traceback'], ''),
        ('r3: no column extents, lower window end at 1 d', dict(zlo_d=1.0, extents=False), [], {'F2': U, 'F5 height': U, 'F3': U},
         ['missing geometry'], [], ''),
        ('r3: 30 ms audit (every window qualifies)', dict(cad=0.03, t_end=2.49), [], {'F1': 'PASS', 'F2': 'FAIL'}, [], ['unmeasurable'], ''),
        ('r3: only the final rows in the log', dict(mutate=lambda L: [l for l in L if not l.startswith('EL_BED_AUDIT')
                                                                     or float(l.split()[2]) > 2.0]), [],
         {'F1': U, 'F5 count': U, 'F5 height': 'PASS'}, [], [], ''),
        ('r3: column translated down by 20 mm, 60 mm late', dict(z0=-0.020, zmax_fn=hi_after), [], {'F5 height': 'FAIL', 'F3': 'PASS'},
         ['60.00 mm above the inlet plane'], [], ''),
        # ---- review r4 -------------------------------------------------------------------------------
        ('r4: redundant knot inside the 0.20 m/s hold (equivalent table)', dict(table=tab_knot), [], REF, [], [], 'same'),
        ('r4: lowest hold subdivided into 80 ms pieces, 60 mm in it', dict(table=tab_sub, zmax_fn=hi_after), [],
         {'F5 height': 'FAIL'}, ['lowest velocity 0.050 m/s (2.020-2.500 s'], [], ''),
        ('r4: one audit row of a fluidised window removed', dict(mutate=_drop(0.85)), [], {'F1': U, 'F2': U, 'F3': U, 'F5 count': U},
         ['missing or rejected'], [], ''),
        ('r4: dp_dyn NaN on one fluidised row', dict(mutate=_nan_field(0.85, 'dp_hydro_dyn', 2)), [], {'F1': U, 'F5 count': 'PASS'},
         ['non-finite dp_dyn'], [], ''),
        ('r4: one fixed step with a raised gradient, one of its rows removed',
         dict(grad_fn=lambda t, dp: dp / 0.1 * (1.5 if 1.5 < t <= 2.0 else 1.0), mutate=_drop(1.8)), [], {'F2': U}, [], [], ''),
        ('r4: run stopped at 2.15 s inside the lowest hold (MaxNumStep)', dict(t_end=2.15), [],
         {'F1': U, 'F2': U, 'F3': U, 'F5 height': 'PASS'}, ['shorter than the evaluation window'], [], ''),
        ('r4: every other audit row removed, 60 mm on a removed row',
         dict(zmax_fn=lambda t: 0.060 if abs(t - 2.41) < 1e-9 else 0.0241, mutate=_thin(2)), [], {'F5 height': U, 'F1': U}, [], [], ''),
        ('r4: seven identical copies of one audit row', dict(mutate=_dup(0.79, 7)), [], REF, [], [], ''),
        ('r4: conflicting copy of one audit row', dict(mutate=_conflict(0.79)), [], {'F1': U}, ['conflicting copies'], [], ''),
        ('r4: MaxSimTime 2.4999, the solver completes 2.5 s, 60 mm on the last row',
         dict(tmax=2.4999, zmax_fn=lambda t: 0.060 if t > 2.499 else 0.0241), [], {'F5 height': 'FAIL'}, [], [], ''),
        ('r4: Reload_SimPar extends MaxSimTime from 2.3 to 2.5 s, 60 mm at 2.4 s',
         dict(mutate=_reload_extension(2.3, 2.3), zmax_fn=lambda t: 0.060 if abs(t - 2.4) < 1e-9 else 0.0241), [],
         {'F5 height': 'FAIL'}, [], [], ''),
        ('r4: unfinished trailing EL_NEWTON_PAIR', dict(mutate=_append('EL_NEWTON_PAIR t= 2.5 step=')), [], REF, [], ['Traceback'], ''),
        ('r4: last line without its line end', dict(mutate=_append('EL_BED_AUDIT t=   2.502500E+00 step= 10010 p_in_out='),
                                                    last_complete=False), [], REF, [], ['Traceback'], ''),
        ('r4: malformed interior table row', dict(mutate=lambda L: [l.replace('3.00000000E-01', '3.0000000xE-01', 1)
                                                                    if l.lstrip().startswith('EL_INLET_TABLE') and ' 1.00000000E+00' in l else l for l in L]), [],
         {'F1': U, 'F5 height': U}, ['table echo incomplete'], ['Traceback'], ''),
        ('r4: particles lost on the last audit row (count and weight)', dict(n_fn=lambda t: 9789 if t < 2.495 else 9691), [],
         {'F1': U, 'F5 count': 'FAIL'}, [], [], ''),
        ('r4: wenyu echo, default invocation', dict(drag='wenyu'), [], {'F1': 'PASS'}, ['wenyu U_mf at the whole-bed'], [], ''),
        ('r4: explicit --closure that differs from the echo', dict(), ['--closure', 'ergun'], REF,
         ['comparison with --closure ergun, not gated'], [], ''),
        ('r4: no width echo, deck 3.0d0 with a comment', dict(width_echo=False, deck=deck_full + 'SimPar@ELKernelWidthFactor = 3.0d0 ! mapping width in diameters\n'),
         [], {'F2': U}, ['mapping support 3 d [deck]'], [], ''),
        ('r4: no width echo, unreadable deck value', dict(width_echo=False, deck=deck_full + 'SimPar@ELKernelWidthFactor = wide\n'),
         [], {'F2': U}, ['value unreadable'], [], ''),
        ('r4: NaN width echo', dict(width='NaN'), [], {'F2': U}, ['mapping support not known'], [], ''),
        ('r4: centre exactly at mid-height', dict(zmax_fn=lambda t: 0.050 if t > 2.3 else 0.0241), [], {'F5 height': 'PASS'},
         ['50.00 mm above the inlet plane', 'not above'], [], ''),
        ('r4: no echo, steps file with two 0.05 m/s rows, 60 mm in the second',
         dict(echo_table=False, protocol=[(0.5, 0.3), (1.0, 0.2), (1.5, 0.05), (2.0, 0.1), (2.5, 0.05)],
              steps_file=[(0.5, 1.0, 0.3), (1.0, 1.5, 0.2), (1.5, 2.0, 0.05), (2.0, 2.5, 0.1), (2.5, 3.0, 0.05)],
              zmax_fn=lambda t: 0.060 if t > 2.8 else 0.0241), ['--steps', S], {'F5 height': 'FAIL'},
         ['[protocol from the steps file, not from the echo]'], [], ''),
        ('r4: no echo, overlapping steps-file rows', dict(echo_table=False, steps_file=FULL_STEPS + [(2.2, 2.4, 0.05)]),
         ['--steps', S], {'F1': U, 'F5 height': U}, ['overlap or repeat'], [], ''),
        ('r4: column translated by +10 m (7-digit coordinates)', dict(z0=10.0), [], REF, [], [], ''),
        ('r4: --d 0 refused at parsing', dict(), ['--d', '0'], {}, ['--d must be finite and positive'], ['ZeroDivision'], ''),
        ('r4: --mu 0 refused at parsing', dict(), ['--closure', 'wenyu', '--mu', '0'], {}, ['--mu must be finite and positive'],
         ['ZeroDivision'], ''),
    ]


def _run_case(run, extra):
    out, err = io.StringIO(), io.StringIO()
    try:
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            main([run] + [x.format(run=run) for x in extra])
    except SystemExit:
        pass
    except BaseException:                            # noqa: BLE001  (a crash is a failed check)
        import traceback
        out.write('Traceback\n' + traceback.format_exc())
    return out.getvalue(), err.getvalue()


def verdicts_of(text):
    m = re.search(r'^verdicts: (.*)$', text, re.M)
    if not m:
        return {}
    V = {}
    for part in m.group(1).replace(FROM_STEPS, '').split(' | '):
        g, v = part.rsplit(' ', 1)
        V[g] = v
    return V


def self_test(keep=None):
    cases = self_test_cases()
    nfail = 0
    ref_out = None
    ctx = tempfile.TemporaryDirectory(prefix='fl11_post_selftest_') if keep is None else contextlib.nullcontext(keep)
    with ctx as tmp:
        for i, (name, kw, extra, expect, must, mustnot, note) in enumerate(cases):
            run = _write_run(os.path.join(tmp, 'run%02d' % i), **kw)
            out, err = _run_case(run, extra)
            text = out + err
            V = verdicts_of(out)
            body = out.split('\n', 1)[1] if '\n' in out else out
            if i == 0:
                ref_out = body
            bad = [f'{g}: {V.get(g)} instead of {v}' for g, v in expect.items() if V.get(g) != v]
            bad += [f'missing {m!r}' for m in must if m not in text] + [f'unexpected {m!r}' for m in mustnot if m in text]
            if note == 'same' and body != ref_out:
                bad.append('output differs from the reference case')
            print(f"self-test {i:2d} {'ok  ' if not bad else 'FAIL'} {name}" + ('' if not bad else f": {bad}"))
            if bad:
                nfail += 1
                print('    ' + text.replace('\n', '\n    '))
            with open(os.path.join(run, 'case.txt'), 'w') as fh:
                fh.write('%s\n%s\n' % (name, ' '.join(x.format(run=run) for x in extra)))
    print(f"fl11_post self-test: {len(cases) - nfail}/{len(cases)} passed")
    return 1 if nfail else 0


# ===================================================================================================
# metamorphic tests on a real log
# ===================================================================================================
KEEP_PREFIX = ('EL_BED_AUDIT', 'EL_FORCE_BUDGET', 'EL_NEWTON_PAIR', 'EL_INLET_TABLE', 'EL ', 'Global CFD domain extents',
               'time:')


def _filtered(lines):
    """the lines the tool reads (the rest of the log does not change its output; checked once)"""
    return [l for l in lines if l.lstrip().startswith(KEEP_PREFIX) or FINISH in l]


def _table_edit(lines, new_rows):
    """replace the echoed table rows (and the row count of the echo)"""
    out, done = [], False
    for l in lines:
        if l.lstrip().startswith('EL_INLET_TABLE'):
            if not done:
                out += ['  EL_INLET_TABLE %s %s' % (_e(t, 16, 8), _e(u, 16, 8)) for t, u in new_rows]
                done = True
            continue
        m = ECHO['table_rows'].match(l)
        if m:
            l = l.replace('(%s rows' % m.group(1), '(%d rows' % len(new_rows))
        out.append(l)
    return out


def self_test_real(rundir, nvar=300, seed=20261122):
    rnd = random.Random(seed)
    with open(os.path.join(rundir, 'run_slurm.log'), errors='replace') as fh:
        full = fh.read().split('\n')
    if full and full[-1] == '':
        full.pop()
    deck = open(os.path.join(rundir, '_data', 'q2p1_param.dat'), errors='replace').read()
    base = _filtered(full)
    ref_text, ref_V = run_lines(base, deck)
    full_text, _ = run_lines(full, deck)
    print(f"real log {rundir}: {len(full)} lines, {len(base)} read by the tool; intact verdicts: "
          + ' | '.join(f'{g} {ref_V[g]}' for g in GATES))
    nfail = 0
    if full_text != ref_text:
        nfail += 1
        print('FAIL: the full log and the lines the tool reads give different output')
    table = Log(base).table
    inner = [h for h in holds_of(table, 0.0, None) if not h['open']]
    audit_idx = [i for i, l in enumerate(base) if l.startswith('EL_BED_AUDIT')]

    def equiv():
        kind = rnd.choice(['knot in a hold', 'subdivided hold', 'last knot removed', 'knot in a ramp', 'trailing record'])
        if kind == 'knot in a hold':
            h = rnd.choice(inner)
            t = round(rnd.uniform(h['a'] + 0.001, h['b'] - 0.001), 6)
            return kind, _table_edit(base, sorted(table + [(t, h['u'])])), True
        if kind == 'subdivided hold':
            h = rnd.choice(inner)
            k = rnd.randint(2, 6)
            extra = [(round(h['a'] + (h['b'] - h['a']) * j / k, 6), h['u']) for j in range(1, k)]
            return kind, _table_edit(base, sorted(table + extra)), True
        if kind == 'last knot removed':
            return kind, (_table_edit(base, table[:-1]) if len(table) > 1 and table[-1][1] == table[-2][1] else base), True
        if kind == 'knot in a ramp':
            (ta, ua), (tb, ub) = rnd.choice([(p, q) for p, q in zip(table, table[1:]) if p[1] != q[1]])
            t = 0.5 * (ta + tb)
            return kind, _table_edit(base, sorted(table + [(t, ua + (ub - ua) * (t - ta) / (tb - ta))])), True
        tail = rnd.choice(['EL_NEWTON_PAIR t= 7.5 step=', 'EL_BED_AUDIT t=   7.502500E+00 step= 30010 p_in_out=',
                           'EL_FORCE_BUDGET step= 30001 drag=', 'time:  0.7500D+01/', '  EL_INLET_TABLE   7.6'])
        return kind, base + [tail], rnd.random() < 0.5

    def degrade():
        kind = rnd.choice(['one row deleted', 'rows deleted', 'row duplicated', 'conflicting copy', 'field NaN',
                           'log truncated', 'audit thinned', 'count perturbed', 'weight perturbed', 'finish marker removed'])
        L = list(base)
        if kind == 'one row deleted':
            del L[rnd.choice(audit_idx)]
        elif kind == 'rows deleted':
            for i in sorted(rnd.sample(audit_idx, rnd.randint(2, 20)), reverse=True):
                del L[i]
        elif kind == 'row duplicated':
            i = rnd.choice(audit_idx)
            L.insert(i + 1, L[i])
        elif kind == 'conflicting copy':
            i = rnd.choice(audit_idx)
            tok = L[i].split()
            j = rnd.choice([k for k, x in enumerate(tok) if k > 4 and NUM_RE.match(x)])
            tok[j] = _e(num(tok[j]) * 1.001 + 1e-9)
            L.insert(i + 1, ' '.join(tok))
        elif kind == 'field NaN':
            i = rnd.choice(audit_idx)
            tok = L[i].split()
            tok[rnd.choice([k for k, x in enumerate(tok) if k > 0 and not KEY_TOK.match(x)])] = 'NaN'
            L[i] = ' '.join(tok)
        elif kind == 'log truncated':
            L = L[:rnd.randrange(len(L)) + 1]
            if rnd.random() < 0.5 and L[-1]:
                L[-1] = L[-1][:rnd.randrange(len(L[-1]))]
                return kind, L, False
        elif kind == 'audit thinned':
            m = rnd.randint(2, 5)
            drop = set(audit_idx[rnd.randrange(m)::m])
            L = [l for i, l in enumerate(L) if i not in drop]
        elif kind in ('count perturbed', 'weight perturbed'):
            i = rnd.choice(audit_idx)
            tok = L[i].split()
            if kind == 'count perturbed':
                j = tok.index('n=') + 1
                tok[j] = str(int(tok[j]) + rnd.choice([-1, 1]) * rnd.randint(1, 100))
            else:
                j = tok.index('weight_z_floor_z=') + 1
                tok[j] = _e(num(tok[j]) * (1 + rnd.choice([-1, 1]) * rnd.uniform(1e-4, 1e-2)))
            L[i] = ' '.join(tok)
        else:
            L = [l for l in L if FINISH not in l]
        return kind, L, True

    stats = {}
    n_any = 0
    for v in range(nvar):
        is_eq = v % 3 == 0
        kind, L, complete = equiv() if is_eq else degrade()
        try:
            text, V = run_lines(L, deck, last_complete=complete)
        except BaseException as exc:                 # noqa: BLE001
            text, V = 'Traceback ' + repr(exc), {}
        key = ('equivalence: ' if is_eq else 'degrading: ') + kind
        st = stats.setdefault(key, [0, 0, 0])
        st[0] += 1
        if is_eq:
            ok = text == ref_text
        else:
            ok = bool(V) and all(V.get(g) in (ref_V[g], 'unavailable') for g in GATES)
            if any(V.get(g) == 'unavailable' for g in GATES):
                st[2] += 1
                n_any += 1
        if not ok:
            st[1] += 1
            nfail += 1
            if st[1] <= 2:
                print(f"FAIL variant {v} ({key}): verdicts {V}")
                if is_eq:
                    import difflib
                    print('\n'.join(list(difflib.unified_diff(ref_text.splitlines(), text.splitlines(), lineterm=''))[:20]))
    neq = sum(s[0] for k, s in stats.items() if k.startswith('equivalence'))
    print(f"{nvar} variants ({neq} equivalence, {nvar - neq} degrading); degrading variants with at least one verdict "
          f"unavailable: {n_any}; failures: {nfail}")
    for k in sorted(stats):
        print(f"   {k:42s} {stats[k][0]:4d} run, {stats[k][1]:3d} failed"
              + (f", {stats[k][2]:3d} with a verdict unavailable" if k.startswith('degrading') else ''))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
