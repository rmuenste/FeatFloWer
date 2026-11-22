#!/usr/bin/env python3
"""Case 11 (minimum fluidisation) evaluation of a q2p1_el_pipeflow run.

Supported domain (design note section 7.4). The tool certifies ONE thing: a complete, uninterrupted
run of an echoed inlet protocol. The official gates F1-F5 are only issued when a precondition, checked
once before any gate (official_preconditions), holds: (a) the inlet table echo is present and complete,
(b) one clock segment (every ordered time line is t0 + k dt of one constant dt, steps 1, 2, ... without
a gap, t0 = the first table row, no reload-changed clock), (c) the run finished and completed the time
of the last table row, (d) every physical input of the closures (d, rho_p, rho_f, mu, g, drag model,
kernel support, column extents, N) is resolved from the run, with its provenance printed. Otherwise
EVERY official gate is unavailable with the reasons (fail closed) and the numbers are printed as
measurements of a partial or non-standard run; the only verdicts that remain are witnessed failures
(F5 count change, F5 height exceedance in a validated row inside the intended lowest hold). A PASS
never comes out of a run that fails the precondition.

The tool is built from five pieces, each defined once:
  A  record reader   one tokenizer for every EL_ record the tool reads (EL_BED_AUDIT, EL_FORCE_BUDGET,
                     EL_NEWTON_PAIR, EL_INLET_TABLE): complete or nothing; an unfinished trailing record
                     (after the last complete one, or a last line without its line end) is ignored with
                     a note on stderr, a malformed interior record is diagnosed. The start-up echo (kernel
                     width, drag model, fluid gravity, audit frequency, table row count, domain extents,
                     particle radius / density / count, fluid density, gravity, closure viscosity), the
                     time lines and the finish marker in the same pass; the deck with the application's
                     syntax; example.json for the particle properties when the log lacks them.
  B  protocol        the intended holds from the echoed table alone: the solver's rule U(t), maximal
                     constant-velocity holds (adjacent equal-valued segments merged before the 0.1 s
                     threshold) from the first row on, the last value held up to the end of a finished
                     run. Never clipped to the observed interval: what the protocol requires is fixed
                     before what the run observed.
  C  clock / audit   one clock segment from the time lines; the expected audit rows (multiples of the
                     echoed frequency up to the last completed step) at t0 + k dt; observed rows validated
                     (step vs time, pressure identities to print precision, count domain, count and weight
                     changing together), identical duplicates collapsed, conflicting ones rejected.
  D  per step        the evaluation window (last 0.3 s) inside the intended hold and covered by validated
                     rows; weight, area and count constant in the window; 'measured' or 'unmeasurable'
                     (reason). Each printed sub-window mean lies inside the hold with its own coverage.
  E  verdicts        PASS / FAIL / unavailable (reason). F1-F3 unavailable when any protocol step is
                     unmeasurable; F2 and F3 need the MEASURED plateau of F1; F3 needs two distinct nonzero
                     fixed velocities and a full-rank fit; F5 height: a centre exactly at mid-height passes,
                     one above it by less than the print uncertainty is undecidable, above that band FAIL.
                     Physical command-line arguments or a --closure that differ from the run's values give
                     a labelled comparison, never an official verdict.

Gates of the spec (section 5, as amended on 2026-11-20):
  F1  plateau dp A / buoyant weight, mean over the steps with dp* > 0.95           (1 +- 0.02)
  F2  fixed branch (dp* < 0.8 of the measured plateau): window gradient (gradz_win_dyn; gradz_win for
      older logs when dp_hydro = 0 on every contributing line) against the closure at the mapped gas
      fraction (+-5 %), with the window excluding max(2 d, mapping support) at both ends.
  F3  U_mf where the fit dp* = a U + b U^2 of the fixed branch (through the origin) meets dp* = 1, against
      the closure U_mf at the whole-bed porosity 1 - N V_p / (A (z99 + r - zmin))      (+-5 %)
  F4  U_mf against the measurement 0.182 m/s (reported)
  F5  count: particle count constant over the run; height: no centre above zmin + Lz/2 during the
      lowest-velocity hold(s) of the protocol.
A --steps file is a cross-check only (three columns, multiplicity); without a table echo there is no
official gate, the steps file then gives a labelled measurement summary.

usage: fl11_post.py RUNDIR [--steps FILE] [--d D] [--rhop R] [--rhof R] [--mu MU] [--g G] [--lz LZ]
                    [--closure NAME] [--flag 0.02] [--series OUT.csv]
       fl11_post.py --self-test                synthetic logs (the counterexamples of reviews r2-r5)
       fl11_post.py --self-test-real RUNDIR [--variants N] [--seed S]
                                               metamorphic tests on a real log (the transforms applicable
                                               to it; equivalence: identical output; degrading: the
                                               expected outcome per gate; coherent changes of run history:
                                               every official gate unavailable with its reason)
fl11_movie.py uses evaluate() and classify() of this module."""
import argparse
import contextlib
import io
import json
import math
import os
import random
import re
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fl11_closures import CLOSURES, closure_grad, closure_umf   # noqa: E402  (re-exported)

EVAL_VERSION = 'fl11_post r5 2026-11-22'   # part of the movie cache signature
NAN = float('nan')
INF = float('inf')
GATES = ('F1', 'F2', 'F3', 'F5 height', 'F5 count')
HOLD_MIN = 0.1           # s: a hold of this length or less is part of a transition, not a step
EVAL_LEN = 0.3           # s: evaluation window at the end of a step
FIRST_LEN = 0.5          # s: the first hold (preparation + step 0) is evaluated over its last 0.5 s
MIN_SAMPLES = 10         # expected audit rows in a window used for statistics or extrema
MIN_SUB = 5              # ... in a printed sub-window
TOL_D4 = 5e-4            # relative print precision of the 4-digit D format of the time line
REL_DET = 1e-6           # F3: smallest relative determinant of the normal equations (full rank)
PHYS = ('d', 'rhop', 'rhof', 'mu', 'g')
NOT_OFFICIAL = 'not a verdict (outside the supported domain)'
OUTSIDE = ' (outside the supported domain)'


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


def d4_prec(x):
    """half the spacing of the D12.4 print of the time line (0.dddd x 10^E)"""
    if not finite(x) or x == 0:
        return 0.0
    return 0.5 * 10.0 ** (math.floor(math.log10(abs(x))) - 3)


def mean(v):
    return sum(v) / len(v) if v else NAN


def _fmin(v):
    v = [x for x in v if finite(x)]
    return min(v) if v else NAN


def _fmax(v):
    v = [x for x in v if finite(x)]
    return max(v) if v else NAN


def _f(x, spec):
    """a number, or n/a for a value that is not printed"""
    return format(x, spec) if finite(x) else format('n/a', '>' + re.match(r'\d*', spec).group(0))


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
    'spheres': re.compile(r'^EL terminal-velocity setup:\s*(\d+) free sphere\(s\) created \(r=([^,\s]+), rho=([^,\s]+)'),
    'rhof': re.compile(r'^Prop\s+-\s+Density\s*=\s*(\S+)'),
    'gravity': re.compile(r'^Prop\s+-\s+Gravity\s*=\s*(\S+)\s+(\S+)\s+(\S+)'),
    'nu_closure': re.compile(r'^EL fluid viscosity: diffusion operator rebuilt with mu =\s*(\S+)\s*\(nu =\s*(\S+) for the closures\)'),
    'nu_kin': re.compile(r'^Prop\s+-\s+Calculated Kinematic Viscosity\s*=\s*(\S+)'),
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
            if st.startswith(('EL ', 'Global CFD', 'Prop')):
                for key, pat in ECHO.items():
                    if key not in self.echo:
                        m = pat.match(st)
                        if m:
                            self.echo[key] = m.groups()
        self.trailing = [b for b in bad if b[0] > good_last]
        self.malformed = [b for b in bad if b[0] <= good_last]

    def e(self, key, i=0):
        v = self.echo.get(key)
        return v[i].strip() if v else None


def read_deck(text):
    """SimPar / Prop keys of a deck with the application's syntax: 'Section@Key = value', the last '@'
    and the last '=' of the line, the value split at blanks / commas / '/' (a trailing '! comment' is
    not part of the first token), a later line overrides an earlier one. Values: lists of text."""
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
        name = key[0] if sec[0] == 'SimPar' else 'Prop@' + key[0]
        keys[name] = val[:1] if sec[0] == 'SimPar' else val
    return keys


def read_deck_num(keys, key, i=0):
    """(value, status): status 'deck', 'absent' or 'unreadable'"""
    if key not in keys:
        return NAN, 'absent'
    try:
        return num(keys[key][i].strip('"\'')), 'deck'
    except (ValueError, IndexError):
        return NAN, 'unreadable'


def deck_text(keys, key):
    v = keys.get(key)
    return v[0].strip('"\'') if v else None


# ===================================================================================================
# B. protocol: the intended holds
# ===================================================================================================
def u_of_t(table, t):
    """the solver's rule: linear between the rows, constant beyond the ends"""
    if t <= table[0][0]:
        return table[0][1]
    for (t0, u0), (t1, u1) in zip(table, table[1:]):
        if t <= t1:
            return u0 + (u1 - u0) * (t - t0) / (t1 - t0)
    return table[-1][1]


def holds_of(table, t_end=None):
    """the intended holds of the protocol: maximal constant-velocity intervals of U(t) from the first
    table row on, the last value held on to t_end (the end of a finished run) when that is later;
    equal-valued adjacent segments merged first, then holds of HOLD_MIN or less dropped.
    Each hold: dict(u, a, b)."""
    if not table:
        return []
    segs = [(p[0], q[0], p[1], q[1]) for p, q in zip(table, table[1:])]
    if t_end is not None and finite(t_end) and t_end > table[-1][0]:
        segs.append((table[-1][0], t_end, table[-1][1], table[-1][1]))
    merged = []
    for a, b, ua, ub in segs:
        if ua != ub:
            merged.append(None)                         # a ramp separates two holds
        elif merged and merged[-1] is not None and merged[-1][2] == ua and merged[-1][1] == a:
            merged[-1] = (merged[-1][0], b, ua)
        else:
            merged.append((a, b, ua))
    return [dict(u=h[2], a=h[0], b=h[1]) for h in merged if h is not None and h[1] - h[0] > HOLD_MIN]


def steps_of(holds):
    """protocol steps of the holds: t1 = end of the hold, t0 = end of the previous hold (the start of the
    transition into it); the first hold is evaluated over its last FIRST_LEN s when longer than 0.6 s"""
    steps = []
    for i, h in enumerate(holds):
        t1 = h['b']
        t0 = holds[i - 1]['b'] if i > 0 else h['a']
        if i == 0 and t1 - t0 > 0.6:
            t0 = t1 - FIRST_LEN
        steps.append(dict(k=i, t0=t0, t1=t1, U=h['u'], hold=h))
    return steps


# ===================================================================================================
# C. clock and audit sequence
# ===================================================================================================
class Clock:
    """time step, start, completed steps, audit frequency, the clock-segment findings of the time lines,
    and the validated audit rows"""

    def __init__(self, log, deck):
        self.notes = []
        self.segment = []                # reasons the time lines are not one clock segment
        self.freq, self.freq_src = None, 'unknown'
        f_echo = log.e('audit_freq')
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
        if not times:
            self.segment.append('no time lines in the log')
        if not self.const:
            self.segment.append('time step not constant (printed dt %s)' % ', '.join('%.4g' % x for x in dts))
            self.notes.append('time step not constant: no step clock, step times from the audit rows only')
        elif finite(pos(dt_deck)) and (not times or agree(dt_deck, dt_print)):
            self.dt, self.dt_src = dt_deck, 'deck TimeStep' + (', time lines agree' if times else '')
        elif finite(pos(dt_audit)) and (not times or agree(dt_audit, dt_print)):
            self.dt, self.dt_src = dt_audit, 'audit rows' + (', time lines agree' if times else '')
        elif finite(pos(dt_print)):
            self.dt, self.dt_src = dt_print, 'time lines (4 digits)'
        self.exact = finite(self.dt)
        t_line = times[0][0] - times[0][2] * self.dt if times and self.exact else NAN
        t_deck, _ = read_deck_num(deck, 'StartSimTime')
        tol0 = d4_prec(times[0][0]) + 1e-3 * self.dt + 1e-12 if times and self.exact else 0.0
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
        if not self.exact and self.const and times:
            self.segment.append('no constant time step established')
        # one clock segment: steps 1, 2, ... in order, one t_max / nsteps, t = t0 + k dt on every line
        if times:
            ks = [x[2] for x in times]
            if ks != list(range(1, len(ks) + 1)):
                j = next((i for i, k in enumerate(ks) if k != i + 1), 0)
                self.segment.append('the time lines are not the consecutive steps 1, 2, ... (line %d shows step %d)'
                                    % (j + 1, ks[j]))
            if len({(x[1], x[3]) for x in times}) > 1:
                self.segment.append('MaxSimTime or MaxNumStep changed during the run (Reload_SimPar)')
            if self.exact:
                for t, _, k, _, _ in times:
                    tk = self.t0 + k * self.dt
                    if abs(t - tk) > d4_prec(t) + d4_prec(tk) + 1e-3 * self.dt:
                        self.segment.append('the time line of step %d shows t = %.4g s, the clock t0 + k dt gives %.6g s'
                                            % (k, t, tk))
                        break
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
            r = row_of(rec)
            why = row_defect(r)
            if why:
                self.rejected.append((k, why))
                continue
            cand[k] = r
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
        self.t_last = self.tk(self.K) if self.K > 0 else NAN
        self.end = self.t_last if self.finished and finite(self.t_last) else None
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
        if a < self.t0 - 1e-3 * self.dt:
            return False, 'the window %.4f-%.4f s starts before the clock start %.4f s' % (a, b, self.t0)
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


def row_defect(r):
    """the identities the writer guarantees for one audit row (EL_BED_AUDIT: dp = p_in - p_out, dp_dyn =
    dp - dp_hydro, to print precision; a value is unavailable exactly where its inputs are) and the count
    domain; '' for a consistent row"""
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
    return ''


# ===================================================================================================
# the evaluation (shared with fl11_movie.py)
# ===================================================================================================
def setup_of(log, deck, rundir):
    """mapping support, drag model, fluid gravity, extents and the physical inputs with their provenance"""
    s = dict(fluid_gravity=None, zmin=None, zmax=None, extents_prec=(0.0, 0.0), phys={}, N=None, N_src=None)
    w = log.e('width')
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
    dm = log.e('drag')
    if dm is not None:
        s['drag'], s['drag_src'] = dm.lower(), 'log echo'
    elif deck_text(deck, 'ELDragModel'):
        s['drag'], s['drag_src'] = deck_text(deck, 'ELDragModel').lower(), 'deck'
    else:
        s['drag'], s['drag_src'] = None, 'unknown'
    fg = log.e('fluid_gravity')
    if fg in ('T', 'F'):
        s['fluid_gravity'] = (fg == 'T')
    ex = log.e('extents')
    if ex is not None:
        try:
            v = [num(x) for x in ex.split()]
            if len(v) == 6 and finite(v[4]) and finite(v[5]) and v[5] > v[4]:
                s['zmin'], s['zmax'] = v[4], v[5]
                s['extents_prec'] = (es_prec(v[4]), es_prec(v[5]))
        except ValueError:
            pass
    js = {}
    p = os.path.join(rundir, 'example.json') if rundir else ''
    if p and os.path.exists(p):
        try:
            js = json.load(open(p))
        except (ValueError, OSError):
            js = {}

    def first(*cands):
        for val, src in cands:
            try:
                v = num(val) if isinstance(val, str) else val
            except ValueError:
                continue
            if finite(v):
                return v, src
        return NAN, 'not resolved'

    def jnum(key):
        v = js.get(key)
        return v if isinstance(v, (int, float)) and not isinstance(v, bool) else NAN
    ph = s['phys']
    sph = log.echo.get('spheres')
    r_echo = sph[1] if sph and NUM_RE.match(sph[1]) else None
    ph['d'] = first((2 * num(r_echo) if r_echo else NAN, 'log echo (2 r of the sphere setup)'),
                    (2 * jnum('benchRadius_'), 'example.json benchRadius_'))
    ph['rhop'] = first((sph[2] if sph else 'x', 'log echo (sphere setup)'), (jnum('particleDensity_'), 'example.json particleDensity_'))
    ph['rhof'] = first((log.e('rhof') or 'x', 'log echo (Prop Density)'),
                       (deck['Prop@Density'][0] if deck.get('Prop@Density') else 'x', 'deck Prop@Density'))
    gdeck = deck.get('Prop@Gravity')
    gv = first((log.e('gravity', 2) or 'x', 'log echo (Prop Gravity)'), (gdeck[2] if gdeck and len(gdeck) > 2 else 'x', 'deck Prop@Gravity'))
    ph['g'] = (abs(gv[0]), gv[1]) if finite(gv[0]) else gv
    nu = first((log.e('nu_closure', 1) or 'x', 'log echo (nu of the closures)'),
               (log.e('nu_kin') or 'x', 'log echo (Prop kinematic viscosity)'))
    ph['mu'] = ((nu[0] * ph['rhof'][0], nu[1] + ' x rho_f') if finite(nu[0]) and finite(ph['rhof'][0]) else (NAN, 'not resolved'))
    if sph and re.fullmatch(r'\d+', sph[0]):
        s['N'], s['N_src'] = int(sph[0]), 'log echo (sphere setup)'
    return s


def evaluate(rundir=None, d=None, steps_path=None, lines=None, deck_text_=None, last_complete=True, name=None,
             flag=0.02):
    """read (A), clock (C), protocol (B), per step (D) and the precondition; the verdicts (E) are report()"""
    if lines is None:
        with open(os.path.join(rundir, 'run_slurm.log'), errors='replace', newline='') as fh:
            text = fh.read()
        last_complete = text.endswith('\n')
        lines = text.split('\n')
        if lines and lines[-1] == '':
            lines.pop()
    if deck_text_ is None:
        p = os.path.join(rundir, '_data', 'q2p1_param.dat') if rundir else ''
        deck_text_ = open(p, errors='replace').read() if p and os.path.exists(p) else ''
    log = Log(lines, last_complete)
    deck = read_deck(deck_text_)
    clock = Clock(log, deck)
    setup = setup_of(log, deck, rundir)
    d_run = setup['phys']['d'][0]
    d_eval = d_run if finite(pos(d_run)) else (d if d else 0.87e-3)
    ev = dict(name=name or os.path.basename(os.path.abspath(rundir)), log=log, deck=deck, clock=clock, setup=setup,
              d=d_eval, r=d_eval / 2, flag=flag, steps_path=steps_path, listed=None, listed_errors=[],
              listed_findings=[], listed_ok=None, protocol_errors=[])
    zmin = setup['zmin'] if setup['zmin'] is not None else NAN
    for r in clock.sorted:
        fb = log.budget.get(r['step'])
        r['dp_force'] = div(fb['drag'][2] + fb['pressure'][2], r['a_in']) if fb else NAN
        r['H'] = r['z99'] + ev['r'] - zmin
        if not r['has_fbal']:
            pr = log.pairs.get(r['step'])
            r['fbal'] = (pr['dp'][2] - pr['expected'][2]) / clock.dt if pr and finite(pos(clock.dt)) else NAN
    table = log.table
    declared = log.e('table_rows')
    if not table:
        ev['protocol_errors'].append('no inlet table echo (EL_INLET_TABLE)')
    else:
        if declared is None:
            ev['protocol_errors'].append('table echo without its row count (EL inlet velocity table = ... rows)')
        elif int(declared) != len(table):
            ev['protocol_errors'].append('table echo incomplete: %d of %s rows readable' % (len(table), declared))
        if any(not q[0] > p[0] for p, q in zip(table, table[1:])) or not all(finite(t) and finite(u) for t, u in table):
            ev['protocol_errors'].append('table echo not strictly increasing and finite')
    if steps_path:
        ev['listed'], ev['listed_errors'] = read_steps(steps_path)
    ev['source'] = 'echo' if table else None
    ev['knots'] = table
    ev['horizon'] = (max(clock.end, table[-1][0]) if clock.end is not None else table[-1][0]) if table else NAN
    ev['holds'] = holds_of(table, clock.end) if table else []
    ev['protocol'] = steps_of(ev['holds'])
    if ev['listed'] is not None and table:
        ev['listed_ok'], ev['listed_findings'] = compare_steps(ev['listed'], ev['protocol'])
        if ev['listed_errors']:
            ev['listed_ok'] = False
    for r in clock.sorted:
        r['U'] = u_of_t(table, r['t']) if table and finite(r['t']) else NAN
    ev['steps'] = [eval_step(ev, s) for s in ev['protocol']]
    last = clock.sorted[-1] if clock.sorted else None
    ev['area'] = pos(last['a_in']) if last else NAN
    ev['W'] = pos(-last['wz']) if last else NAN
    ev['precondition'] = official_preconditions(ev)
    ev['classes'] = classify(ev['steps'])
    return ev


PRECONDITION_REASONS = """the reasons official_preconditions() can emit (design note 7.4):
  no inlet table echo (EL_INLET_TABLE)
  table echo without its row count / table echo incomplete: N of M rows readable / not strictly increasing and finite
  no time lines in the log / time step not constant (printed dt ...) / no constant time step established
  the time lines are not the consecutive steps 1, 2, ... (line j shows step k)
  MaxSimTime or MaxNumStep changed during the run (Reload_SimPar)
  the time line of step k shows t = ..., the clock t0 + k dt gives ...
  the clock starts at ... s, the protocol (first table row) at ... s
  no finish marker: the run has not finished / the run ended at ... s, before the last table row at ... s
  audit frequency unknown (EL bed audit freq)
  particle diameter d / particle density rho_p / fluid density rho_f / fluid viscosity mu / gravity g not resolved
  drag model of the run ... not one of gidaspow, ergun, wenyu
  mapping support not a positive finite number (kernel width factor ...)
  column extents not in the log (Global CFD domain extents)
  particle count N not resolved (no validated audit row with N > 0) / particle count of the sphere setup differs"""


def official_preconditions(ev):
    """the reasons this run is outside the supported domain (empty: official gates may be issued)"""
    clock, setup, table = ev['clock'], ev['setup'], ev['knots']
    out = list(ev['protocol_errors'])
    out += clock.segment
    if table and clock.exact and finite(clock.t0):
        if abs(clock.t0 - table[0][0]) > 1e-3 * clock.dt + es_prec(table[0][0], 8):
            out.append('the clock starts at %.4f s, the protocol (first table row) at %.4f s' % (clock.t0, table[0][0]))
    if not clock.finished:
        out.append('no finish marker: the run has not finished')
    elif table and finite(clock.t_last) and clock.t_last < table[-1][0] - 1e-3 * (clock.dt if clock.exact else 1.0):
        out.append('the run ended at %.4f s, before the last table row at %.4f s' % (clock.t_last, table[-1][0]))
    if not clock.freq:
        out.append('audit frequency unknown (EL bed audit freq)')
    for key, what in (('d', 'particle diameter d'), ('rhop', 'particle density rho_p'), ('rhof', 'fluid density rho_f'),
                      ('mu', 'fluid viscosity mu'), ('g', 'gravity g')):
        if not finite(pos(setup['phys'][key][0])):
            out.append('%s not resolved from the run (log echo, deck, example.json)' % what)
    if setup['drag'] not in CLOSURES:
        out.append('drag model of the run %s [%s] not one of %s' % (setup['drag'] or 'unknown', setup['drag_src'], ', '.join(CLOSURES)))
    if not finite(pos(setup['width'])):
        out.append('mapping support not a positive finite number (kernel width factor %s [%s])' % (setup['width'], setup['width_src']))
    if setup['zmin'] is None:
        out.append('column extents not in the log (Global CFD domain extents)')
    ns = {r['n'] for r in clock.sorted}
    if not clock.sorted or min(ns) <= 0:
        out.append('particle count N not resolved (no validated audit row with N > 0)')
    elif setup['N'] is not None and setup['N'] != clock.sorted[0]['n']:
        out.append('particle count of the sphere setup %d differs from the first audit row %d' % (setup['N'], clock.sorted[0]['n']))
    return out


def eval_step(ev, s, hold_check=True):
    """D: one step; status 'measured' or the reason it is unmeasurable"""
    clock, flag = ev['clock'], ev['flag']
    h = s['hold']
    x = dict(s)
    x.update(status='measured', why='')

    def fail(why):
        if x['status'] == 'measured':
            x['status'], x['why'] = 'unmeasurable', why
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

    def sub(lo, hi, open_end):
        """a sub-window mean only inside the hold and with its own coverage"""
        if lo < h['a'] - 1e-9 or hi > h['b'] + 1e-9 or not clock.coverage(lo, hi, MIN_SUB)[0]:
            return NAN
        return st('dp_dyn', lo, hi, open_end)[0]
    x['dp'], x['std'] = st('dp_dyn', a, t1)
    m1 = sub(t0 + 0.2, t0 + 0.35, True)
    m2 = sub(t0 + 0.35, t0 + 0.5, False)
    x['mf'], _ = st('dp_force', a, t1)
    x['H'], _ = st('H', a, t1)
    x['eps'], _ = st('eps', a, t1)
    fl, _ = st('fbal', a, t1)
    flpe, _ = st('floor', a, t1)
    x.update(dps=div(x['dp'], wa), stds=div(x['std'], wa), m1s=div(m1, wa), m2s=div(m2, wa), mfs=div(x['mf'], wa),
             fl=div(fl, Wk), flpe=div(flpe, Wk), flag=finite(wa) and finite(m1) and finite(m2) and abs(m1 - m2) > flag * wa)
    x['measured'] = x['status'] == 'measured'
    return x


def classify(steps):
    """branch / display state of every step: 'fluidised' (dp* > 0.95), 'fixed' (dp* < 0.8 of the measured
    plateau), 'intermediate' (neither), 'unclassified' (measured, no measured plateau), 'unmeasurable'"""
    flu = [x for x in steps if x['measured'] and x['dps'] > 0.95]
    plateau = mean([x['dps'] for x in flu]) if flu else NAN
    out = []
    for x in steps:
        if not x['measured']:
            out.append('unmeasurable')
        elif x['dps'] > 0.95:
            out.append('fluidised')
        elif not finite(plateau):
            out.append('unclassified')
        elif x['dps'] < 0.8 * plateau:
            out.append('fixed')
        else:
            out.append('intermediate')
    return out


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
        print(f" {x['k']:3d}  {x['U']:6.3f}  {x['t1']-EVAL_LEN:5.2f}-{x['t1']:5.2f}  {x['dps']:7.4f}  {x['stds']:7.4f}  {_f(x['m1s'], '6.3f')}  {_f(x['m2s'], '6.3f')}  "
              f"{'DIFF' if x['flag'] else '    '}  {x['mfs']:7.4f}  {x['H']*1e3:6.2f}  {x['eps']:7.4f}  {x['fl']:7.4f}  {x['dps'] + x['fl']:8.4f}   ({x['flpe']:.4f})"
              f"{'' if x['measured'] else '  (unmeasurable: ' + x['why'] + ')'}")
    print('  dp* = dp_dyn A / W over 0.2-0.5 s of the step (W, A of that step); floor from the particle momentum '
          'balance; floor_pe = the pe wall-contact impulse (unrelaxed PGS impulse, 1/0.9 of the applied one); '
          'n/a = a sub-window outside the hold or without its own coverage')


class Phys:
    """the physical inputs of a closure evaluation"""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def run_phys(ev):
    ph = ev['setup']['phys']
    return Phys(**{k: ph[k][0] for k in PHYS})


def word(v, official):
    return v if official else NOT_OFFICIAL


def report(ev, a):
    """prints the evaluation, returns the verdicts {gate: PASS | FAIL | unavailable}"""
    log, clock, setup = ev['log'], ev['clock'], ev['setup']
    rows = clock.sorted
    pre = ev['precondition']
    official = not pre
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
    ks = clock.expected(clock.t0, clock.t_last) if clock.K > 0 else None
    print(f"clock: dt {clock.dt:.6g} s [{clock.dt_src}], start {clock.t0:.4f} s [{clock.t0_src}], {clock.K} steps completed "
          f"({'run finished' if clock.finished else 'no finish marker'}), audit every "
          f"{clock.freq if clock.freq else '?'} steps [{clock.freq_src}]: {len(ks) if ks is not None else '?'} expected rows, "
          f"{len(rows)} validated, {len(clock.rejected)} rejected")
    for n_ in clock.notes:
        print('  note: ' + n_)
    if clock.rejected:
        print('  rejected audit rows: ' + '; '.join(f"step {k}: {w}" for k, w in clock.rejected[:3])
              + (' ...' if len(clock.rejected) > 3 else ''))
    ph = setup['phys']
    print('physical inputs of the run: ' + '; '.join(f"{k} {ph[k][0]:.6g} [{ph[k][1]}]" for k in PHYS)
          + f"; drag model {setup['drag']} [{setup['drag_src']}]; kernel width factor {setup['width']:g} [{setup['width_src']}]"
          + (f"; column {setup['zmin']:g}-{setup['zmax']:g} m [log echo]" if setup['zmin'] is not None else '; column extents unknown')
          + (f"; N {setup['N']} [{setup['N_src']}]" if setup['N'] is not None else ''))
    if a.series:
        with open(a.series, 'w') as fh:
            fh.write('t,step,U,dp,dp_hydro,dp_dyn,dp_force,floor_z,z50,z95,z99,zmax,eps_bed,n\n')
            for r in rows:
                fh.write(','.join(str(r[k]) for k in ('t', 'step', 'U', 'dp', 'dp_h', 'dp_dyn', 'dp_force', 'floor',
                                                      'z50', 'z95', 'z99', 'zmax', 'eps', 'n')) + '\n')
    if ev['source'] == 'echo':
        tend = 'the last table row (the run has not finished)' if clock.end is None else f"{ev['horizon']:.4f} s"
        print(f"protocol: {len(ev['protocol'])} intended holds of the echoed table (maximal constant-velocity intervals of "
              f"the solver's U(t) from the first table row), terminal hold to {tend}")
    if pre:
        print('\n=== NO OFFICIAL VERDICT: this run is outside the supported domain (a complete, uninterrupted run of an '
              'echoed inlet protocol) ===')
        for p in pre:
            print('   - ' + p)
        print('=== everything below is a MEASUREMENT of a partial or non-standard run; only witnessed failures '
              '(F5 count change, F5 height exceedance) are verdicts ===')
    if ev['protocol']:
        print_step_table(ev['steps'], 'steps of the protocol' + ('' if official else ' (measurements)'))
    if ev['listed'] is not None and (not ev['listed_ok'] or not ev['source']):
        if ev['source']:
            print(f"steps file {a.steps}: NOT the echoed protocol: " + '; '.join(ev['listed_errors'] + ev['listed_findings']))
        listed = [dict(k=i, t0=t0, t1=t1, U=u, hold=dict(u=u, a=t0, b=t1)) for i, (t0, t1, u) in enumerate(ev['listed'])]
        print_step_table([eval_step(ev, s, hold_check=False) for s in listed],
                         'steps-file summary, NOT gated (measurements over the listed intervals)')
    elif ev['listed'] is not None:
        print(f"steps file {a.steps}: identical to the echoed protocol (all three columns, each plateau once)")
    V = {g: 'unavailable' for g in GATES}
    if ev['source'] and not ev['protocol_errors']:
        V.update(gates_f1_f3(a, ev, official))
        V['F5 height'] = gate_f5_height(a, ev, official)
    elif ev['source']:
        print(f"\nF1, F2, F3, F5 height unavailable (protocol not usable: {'; '.join(ev['protocol_errors'])})")
    else:
        print('\nF1, F2, F3, F5 height unavailable: no official gate without the table echo')
    V['F5 count'] = gate_f5_count(ev, official)
    print('verdicts: ' + ' | '.join(f"{g} {V[g]}" for g in GATES) + ('' if official else OUTSIDE))
    return V


def gates_f1_f3(a, ev, official):
    res = ev['steps']
    V = dict(F1='unavailable', F2='unavailable', F3='unavailable')
    bad = [x for x in res if not x['measured']]
    if len(res) < 3:
        print(f'\nF1, F2, F3 unavailable: {len(res)} protocol step(s), three needed')
        return V
    if bad:
        why = '; '.join(f"step {x['k']} (U {x['U']:.3f}): {x['why']}" for x in bad[:3]) + (' ...' if len(bad) > 3 else '')
        flu = [x for x in res if x['measured'] and x['dps'] > 0.95]
        sub = (f"mean dp* over the {len(flu)} measured steps with dp* > 0.95: {mean([x['dps'] for x in flu]):.6f}"
               if flu else 'no measured step with dp* > 0.95')
        print(f"\nF1, F2, F3 unavailable ({len(bad)} protocol step(s) unmeasurable, each could belong to the set of a "
              f"gate: {why})")
        print(f"   measured, not a verdict: {sub}")
        return V
    flu = [x for x in res if x['dps'] > 0.95]
    if not flu:
        nom = [x for x in res if x['dps'] < 0.8]
        print('\nF1, F2, F3 unavailable: no step with dp* > 0.95 (no measured plateau; the fixed branch is dp* < 0.8 '
              'of the measured plateau)')
        print(f"   measurement only, not a verdict: {len(nom)} step(s) with dp* < 0.8 of the nominal W/A "
              f"(U {', '.join('%.3f' % x['U'] for x in nom)})")
        return V
    plateau = sum(x['dps'] for x in flu) / len(flu)
    v1 = 'PASS' if abs(plateau - 1) <= 0.02 else 'FAIL'
    V['F1'] = v1 if official else 'unavailable'
    print(f"\nF1 plateau: mean dp* over {len(flu)} steps with dp* > 0.95 (U {min(x['U'] for x in flu):.2f} .. "
          f"{max(x['U'] for x in flu):.2f}): {plateau:.6f}   (gate 1 +- 0.02: {word(v1, official)}); "
          f"dp* + floor/W over these {sum(x['dps'] + x['fl'] for x in flu)/len(flu):.4f}")
    fixed = [x for x in res if x['dps'] < 0.8 * plateau]
    if len(fixed) < 2:
        print(f'F2, F3 unavailable: {len(fixed)} fixed-branch step(s) (dp* < 0.8 of the plateau), two needed')
        return V
    P = run_phys(ev)
    clo = ev['setup']['drag']
    if clo in CLOSURES and all(finite(pos(getattr(P, k))) for k in PHYS):
        v3 = f3_block(a, fixed, plateau, clo, P, official, '')
        v2 = gate_f2(a, ev, fixed, clo, P, official, '')
        if official:
            V['F3'], V['F2'] = v3, v2
    else:
        print(f"F2, F3 unavailable: drag model ({clo or 'unknown'}) or physical inputs of the run not resolved")
    # labelled comparison with physical arguments or a closure of the command line that differ from the run
    diff = {k: getattr(a, k) for k in PHYS if getattr(a, k) is not None and getattr(a, k) != getattr(P, k)}
    if diff or (a.closure and a.closure != clo):
        Q = Phys(**{k: diff.get(k, getattr(P, k)) for k in PHYS})
        cclo = a.closure if a.closure else clo
        what = ', '.join([f'{k} {v:g} (run {getattr(P, k):.6g})' for k, v in diff.items()]
                         + ([f'closure {a.closure} (run {clo})'] if a.closure and a.closure != clo else []))
        tag = f' [comparison with the command line {what}, not gated]'
        print(f"\n--- comparison with the command line: {what} (not gated) ---")
        if cclo in CLOSURES and all(finite(pos(getattr(Q, k))) for k in PHYS):
            f3_block(a, fixed, plateau, cclo, Q, False, tag)
            gate_f2(a, ev, fixed, cclo, Q, False, tag)
        else:
            print('   comparison not possible: a physical input is not known')
    return V


def f3_block(a, fixed, plateau, clo, P, official, tag):
    vp = math.pi * P.d ** 3 / 6
    if not all(finite(x['H']) and finite(pos(x['A'])) and finite(x['n']) for x in fixed):
        print(f"F3 unavailable: bed height, area or particle count not finite on every fixed step{tag}")
        return 'unavailable'
    efix = mean([x['eps'] for x in fixed])
    uc = closure_umf(efix, P.d, P.rhop, P.rhof, P.mu, P.g, clo)
    ebed = sum(1 - div(x['n'] * vp, x['A'] * x['H']) for x in fixed) / len(fixed)
    ub = closure_umf(ebed, P.d, P.rhop, P.rhof, P.mu, P.g, clo)
    su2 = sum(x['U']**2 for x in fixed); su3 = sum(x['U']**3 for x in fixed); su4 = sum(x['U']**4 for x in fixed)
    sy1 = sum(x['U']*x['dps'] for x in fixed); sy2 = sum(x['U']**2*x['dps'] for x in fixed)
    det = su2*su4 - su3*su3
    distinct = sorted({x['U'] for x in fixed if x['U'] != 0})
    rel = det / (su2 * su4) if su2 > 0 and su4 > 0 else 0.0
    if len(distinct) < 2 or not rel > REL_DET:
        print(f"F3 unavailable: the fit dp* = a U + b U^2 is not identified (rank-deficient: {len(distinct)} distinct "
              f"nonzero fixed velocit{'y' if len(distinct) == 1 else 'ies'}, relative determinant {rel:.3g}){tag}")
        return 'unavailable'
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
    u1 = cross(1.0); up = cross(plateau)
    r3 = div(u1, ub)
    v = ('PASS' if abs(r3 - 1) <= 0.05 else 'FAIL') if finite(r3) else 'unavailable'
    gated = official and not tag
    print(f"F3 fixed branch {len(fixed)} steps (U {min(x['U'] for x in fixed):.3f} .. {max(x['U'] for x in fixed):.3f}): "
          f"dp* = {ca:.4f} U + {cb:.4f} U^2 (through the origin)")
    print(f"   U_mf at dp* = 1: {u1:.4f} m/s; {clo} U_mf at the whole-bed porosity {ebed:.4f}: {ub:.4f} m/s; "
          f"ratio {r3:.4f}   (gate 0.05: {v if gated else ('comparison, not gated' if tag else NOT_OFFICIAL)}){tag}")
    print(f"   reported, not gated: U_mf at the measured plateau {plateau:.4f}: {up:.4f} m/s = {div(up, ub):.4f} of the closure; "
          f"at the core porosity {efix:.4f} the closure gives {uc:.4f} m/s: ratios {div(u1, uc):.4f} / {div(up, uc):.4f}")
    if not tag:
        print(f"F4 U_mf / 0.182 (NETL measurement): {u1/0.182:.4f} (dp* = 1) / {up/0.182:.4f} (plateau) (reported, not gated)")
        s, i = linfit([x['U'] for x in fixed], [x['dp'] for x in fixed])
        print(f"   secondary (rule before 2026-11-20): linear fit of the fixed branch, slope {s:.1f} Pa s/m (intercept {i:.2f} Pa)")
    return v if gated else 'unavailable'


def valid_eps(e):
    """a mapped porosity the closure can be evaluated at: finite and 0 < eps < 1"""
    return finite(e) and 0 < e < 1


def gate_f2(a, ev, fixed, clo, P, official, tag):
    setup = ev['setup']
    r = P.d / 2
    vp = math.pi * P.d ** 3 / 6
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
        required = max(2 * P.d, width * P.d)
    lows, highs, fields, short, undecided = [], [], [], 0, 0
    for row in rows_all:
        lo = row['zlo'] - zmin if zmin is not None else NAN
        hi = row['z99'] + r - row['zhi']
        fm = row['margin'] if row['margin'] is not None else 2 * P.d
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
        reasons.append('missing geometry: no inlet-plane elevation in the log (Global CFD domain extents)')
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
    print(f"     audit window margin over {len(rows_all)} lines: min {div(min(fin), P.d) if fin else NAN:.3f} d "
          f"(lower end - inlet plane min {_fmin(lows)/P.d:.3f} d, z99 + r - upper end min {_fmin(highs)/P.d:.3f} d, "
          f"{'win_margin field' if all(row['margin'] is not None for row in rows_all) else 'margin 2 d (binary before 627caf47)'} "
          f"min {_fmin(fields)/P.d:.3f} d); required max(2 d, mapping support {width:g} d "
          f"[{setup['width_src']}]) = {div(required, P.d):.3f} d: {'met' if margin_ok else 'NOT met'}")
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
        clo_b = closure_grad(x['U'], x['eps'], P.d, P.rhof, P.mu, clo)
        eh = 1 - div(x['n'] * vp, x['A'] * x['H'])
        cloh = closure_grad(x['U'], eh, P.d, P.rhof, P.mu, clo)
        clok = closure_grad(x['U'], ek, P.d, P.rhof, P.mu, clo) if valid_eps(ek) else NAN
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
    gated = official and not tag
    if not reasons:
        v = 'PASS' if worst_k <= 0.05 else 'FAIL'
        print(f"F2 worst |ratio_k - 1| (window gradient against the closure at the mapped gas fraction) {worst_k:.4f}   "
              f"(gate 0.05: {v if gated else ('comparison, not gated' if tag else NOT_OFFICIAL)}){tag}")
    else:
        v = 'unavailable'
        print(f"F2 measured, not gated: worst |ratio_k - 1| {worst_k:.4f}")
        print(f"F2 verdict unavailable ({'; '.join(reasons)}){tag}")
    if finite(worst_w):
        print(f"   reported, not gated: worst |ratio_w - 1| at the geometric window porosity (rule before 2026-11-20) {worst_w:.4f}")
    return v if gated else 'unavailable'


def gate_f5_count(ev, official):
    clock = ev['clock']
    ns = sorted({r['n'] for r in clock.sorted})
    if len(ns) > 1:
        print(f"F5 count: particle counts {ns} among the validated audit rows (witnessed): FAIL")
        return 'FAIL'
    if not official:
        print(f"F5 count: counts seen {ns}; no official verdict for this run: unavailable")
        return 'unavailable'
    if not clock.sorted:
        why = 'no validated audit row'
    elif ns == [0]:
        why = 'no particle'
    else:
        ok, why = clock.coverage(clock.t0, clock.t_last, min_samples=1)
        if ok:
            print(f"F5 count: counts seen {ns} on {len(clock.sorted)} audit lines, run {clock.t0:.4f}-{clock.t_last:.4f} s "
                  f"covered: PASS")
            return 'PASS'
    print(f"F5 count: counts seen {ns}; count over the run not established ({why}): unavailable")
    return 'unavailable'


def gate_f5_height(a, ev, official):
    """an exceedance in a validated row inside the intended lowest hold(s) is a FAIL whatever else is
    missing; a PASS needs the official precondition, full coverage and every reading at or below the mid-
    column; a reading above it by less than the print uncertainty is undecidable"""
    clock, setup, holds = ev['clock'], ev['setup'], ev['holds']
    if not holds:
        print("F5 height unavailable: the protocol holds no velocity plateau")
        return 'unavailable'
    if setup['zmin'] is None or setup['zmax'] is None:
        print("F5 height unavailable: no column extents in the log (Global CFD domain extents): mid-column unknown")
        return 'unavailable'
    u = min(h['u'] for h in holds)
    low = [h for h in holds if h['u'] == u]
    span = ', '.join(f"{h['a']:.3f}-{h['b']:.3f} s" for h in low)
    reasons, rows = [], []
    for h in low:
        rows += clock.window(h['a'], h['b'])
        ok, why = clock.coverage(h['a'], h['b'])
        if not ok:
            reasons.append(why)
    zmin, lz = setup['zmin'], setup['zmax'] - setup['zmin']
    if a.lz is not None and abs(a.lz - lz) > 1e-9 * max(1.0, lz):
        print(f"   note: --lz {a.lz:g} differs from the column height {lz:g} m of the echoed extents (the extents are used)")
    mid = zmin + lz / 2
    pmid = 0.5 * (setup['extents_prec'][0] + setup['extents_prec'][1])
    zs = [r['zmax'] for r in rows]
    fin = [z for z in zs if finite(z)]
    if len(fin) < len(zs):
        reasons.append(f'{len(zs) - len(fin)} audit line(s) without a centre height')
    above = [z for z in fin if z - mid > pmid + es_prec(z)]
    band = [z for z in fin if 0 < z - mid <= pmid + es_prec(z)]
    zm = max(fin) if fin else NAN
    if fin:
        print(f"F5 (b) highest centre during the plateau of the lowest velocity {u:.3f} m/s ({span}, {len(rows)} audit lines): "
              f"{(zm - zmin)*1e3:.2f} mm above the inlet plane (mid-column {lz*500:.1f} mm)")
    if above:
        print("F5 height: a centre above mid-column observed in a validated row (witnessed): FAIL")
        v = 'FAIL'
    elif band:
        print(f"F5 height unavailable (not decidable at print precision: {len(band)} reading(s) above mid-column by less "
              f"than the print uncertainty)")
        v = 'unavailable'
    elif not official and not fin:
        print(f"F5 height: no validated row in the intended lowest hold ({span}); no official verdict for this run: unavailable")
        v = 'unavailable'
    elif not official:
        print("F5 height: no centre above mid-column in the observed part; no official verdict for this run: unavailable")
        v = 'unavailable'
    elif reasons or not fin:
        print(f"F5 height unavailable ({'; '.join(reasons) if reasons else 'no observation'})")
        v = 'unavailable'
    else:
        print("F5 height: no centre above mid-column (a centre at mid-height is not above): PASS")
        v = 'PASS'
    if fin:
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
    for k, h in (('d', 'particle diameter'), ('rhop', 'particle density'), ('rhof', 'fluid density'),
                 ('mu', 'fluid viscosity'), ('g', 'gravity')):
        ap.add_argument('--' + k, type=positive('--' + k), default=None,
                        help=h + ' for a labelled comparison (the official gates use the value of the run)')
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
    if ev['source'] is None and ev['listed'] is None:
        sys.exit('no EL_INLET_TABLE echo (SimPar@ELInletVelocityTable not set?): no official gate')
    report(ev, a)
    return 0


def run_lines(lines, deck, argv=(), last_complete=True, name='run', rundir=None):
    """the report of an in-memory log (the metamorphic tests): (stdout, verdicts, evaluation)"""
    a = parser().parse_args(list(argv) + ['x'])
    out = io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
        ev = evaluate(rundir, a.d, a.steps, lines=lines, deck_text_=deck, last_complete=last_complete, name=name,
                      flag=a.flag)
        V = report(ev, a)
    return out.getvalue(), V, ev


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
               drag='gidaspow', width_echo=True, mutate=None, tmax=None, n_fn=None, grad_fn=None, last_complete=True,
               mu=1.8e-5, clock_fn=None):
    """a synthetic run directory with the lines the tool reads: start-up echo (with the sphere setup and the
    fluid properties), table echo, a time line for every step, EL_FORCE_BUDGET and EL_BED_AUDIT lines, the finish
    marker. Heights above the inlet plane z0: bed top 24 mm (zmax_fn(t): highest centre), mapped porosity ekern;
    legacy = a binary before c8541cfc; clock_fn(k) -> (t, dt) replaces t = k DT; mutate(lines) edits the lines"""
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
    lines.append('EL terminal-velocity setup: %d free sphere(s) created (r=%r, rho=%g, seedMode=file, seedDomain=box), '
                 'PE decomposition 3x3x4.' % (n, d / 2, rhop))
    lines += ['Prop    - Gravity =   0.00000000E+00  0.00000000E+00 -0.98100000E+01',
              'Prop    - Density =   0.10000000E+01  0.10000000E+01',
              ' EL fluid viscosity: diffusion operator rebuilt with mu = %s (nu = %s for the closures)'
              % (_e(mu * rhof, 16, 8), _e(mu, 16, 8))]
    tmx = t_end if tmax is None else tmax
    for k in range(1, nsteps + 1):
        t, dtk = clock_fn(k) if clock_fn else (k * DT, DT)
        if k % freq == 0:
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
            lines.append('EL_FORCE_BUDGET step= %d drag= %s %s %s pressure= %s %s %s lift= 0 0 0 grav_buoy= %s %s %s' % (
                k, _e(0.0), _e(0.0), _e(dps * W * 0.9), _e(0.0), _e(0.0), _e(dps * W * 0.1), _e(0.0), _e(0.0), _e(-W)))
            lines.append(line)
        if time_lines:
            lines.append('time:%s/%s | itns:%8d/%8d | dt:%s%s%s' % (_d(t), _d(tmx), k, nsteps, _d(dtk), _d(dtk), _d(dtk)))
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


def self_test_cases():
    """(name, generator keywords, extra arguments, expected verdicts, strings that must appear, strings that
    must not, the answer of the tool before this round where it differed)"""
    S = '{run}/_data/steps.txt'
    hi_after = lambda t: 0.060 if t > 2.3 else 0.0241       # noqa: E731
    REF = {'F1': 'PASS', 'F2': 'FAIL', 'F3': 'PASS', 'F5 height': 'PASS', 'F5 count': 'PASS'}
    U = 'unavailable'
    ALLU = {g: U for g in GATES}
    base_tab = _table(PROTO)[0]
    tab_sub = base_tab[:-1] + [(2.1, 0.05), (2.18, 0.05), (2.26, 0.05), (2.34, 0.05), (2.42, 0.05), (2.5, 0.05)]
    tab_knot = sorted(base_tab + [(1.3, 0.20)])
    deck_full = 'SimPar@TimeStep = 0.00025d0\nSimPar@StartSimTime = 0d0\nSimPar@MaxNumStep = 10000\n'
    rep_proto = [(0.5, 0.30), (1.0, 0.04), (1.5, 0.30), (2.0, 0.04)]
    nofl_proto = [(0.5, 0.14), (1.0, 0.12), (1.5, 0.10), (2.0, 0.08)]
    mu_e = 2.4e-5
    grad_mu = lambda t, dp, u: closure_grad(u, 0.37, 0.87e-3, 1.0, mu_e, 'gidaspow')  # noqa: E731
    sub_tab = [(0.0, 0.3), (1.0, 0.3), (1.02, 0.2), (1.4, 0.2), (1.42, 0.1), (1.9, 0.1), (1.92, 0.05), (2.4, 0.05)]
    reset = lambda k: ((k * DT - 0.5) if k > 7000 else k * DT, DT)                  # noqa: E731  time reset at 1.75 s
    dtchg = lambda k: ((k * DT) if k <= 9000 else 9000 * DT + (k - 9000) * 2 * DT, DT if k <= 9000 else 2 * DT)  # noqa: E731
    return [
        ('reference (all gates measurable)', dict(), [], REF, ['stricter, reported', 'physical inputs of the run'], ['Traceback', 'NO OFFICIAL'], ''),
        ('missing inlet plane (area_in = 0)', dict(area=0.0), [], {'F1': U, 'F2': U, 'F3': U, 'F5 height': 'PASS'}, [], ['Traceback'], ''),
        ('empty bed (n = 0)', dict(n=0), [], ALLU, ['particle count N not resolved'], ['Traceback'], ''),
        ('r2: steps file without the lowest step', dict(steps_file=FULL_STEPS[:3]), ['--steps', S], REF,
         ['plateau (2.0000, 2.5000, 0.0500) missing from the list', 'steps-file summary, NOT gated'], [], ''),
        ('r2: foreign step in the steps file', dict(steps_file=FULL_STEPS + [(2.5, 3.0, 0.01)]), ['--steps', S], REF,
         ['is no plateau of the protocol'], [], ''),
        ('r3: steps file with a duplicated row', dict(steps_file=FULL_STEPS + [FULL_STEPS[0]]), ['--steps', S], REF, ['listed 2 times'], [], ''),
        ('r3: steps file with an omitted row', dict(steps_file=FULL_STEPS[1:]), ['--steps', S], REF, ['missing from the list'], [], ''),
        ('r3: steps file with a shifted start', dict(steps_file=[FULL_STEPS[0], (1.2, 1.5, 0.20)] + FULL_STEPS[2:]),
         ['--steps', S], REF, ['differs from the plateau (1.0000, 1.5000, 0.2000) in the start time'], [], ''),
        ('steps file identical to the protocol', dict(steps_file=FULL_STEPS), ['--steps', S], REF, ['identical to the echoed protocol'], ['NOT gated'], ''),
        ('short run (ends inside the lowest hold, no finish marker)', dict(t_end=2.3, finished=False), [], ALLU,
         ['no finish marker: the run has not finished'], [], ''),
        ('sparse audit (cadence 0.4 s)', dict(cad=0.4), [], {'F1': U, 'F5 height': U}, ['audit too sparse'], [], ''),
        ('window margin 2 d with a 2.5 d support', dict(margin_d=2.0), [], {'F2': U, 'F1': 'PASS'}, ['window inside the mapping support'], [], ''),
        ('lower window end at 2 d, field says 2.5 d', dict(zlo_d=2.0), [], {'F2': U}, ['window inside the mapping support'], [], ''),
        ('old log: no win_margin / gradz_win_dyn, zero hydrostatic part', dict(margin_d=2.0, dyn=False), [], {'F2': U},
         ['hydrostatic part zero', 'margin 2 d (binary before 627caf47)'], [], ''),
        ('old log with a 2 d support: margin met, total gradient admitted', dict(margin_d=2.0, dyn=False, width=2.0), [],
         {'F2': 'FAIL'}, ['hydrostatic part zero'], [], ''),
        ('old log with a hydrostatic part', dict(margin_d=2.0, dyn=False, width=2.0, dp_hydro=0.981), [], {'F2': U}, ['nonzero hydrostatic part'], [], ''),
        ('r3: mapped porosity 1.0 (closure zero)', dict(ekern=1.0), [], {'F2': U}, ['no positive closure'], [], ''),
        ('r3: the -1 sentinel in the mapped samples', dict(ekern_fn=lambda t: -1.0 if 2.3 < t < 2.4 else 0.37), [], {'F2': U},
         ['mapped gas-fraction sample(s)'], [], ''),
        ('r3: zero-velocity fixed plateau', dict(protocol=PROTO + [(2.5, 0.0)]), [], {'F2': U}, ['step 4 (U 0.000)'], [], ''),
        ('r3: terminal hold without its last knot, 60 mm late in the hold', dict(table=base_tab[:-1], zmax_fn=hi_after), [],
         {'F5 height': 'FAIL', 'F1': 'PASS'}, ['lowest velocity 0.050 m/s (2.020-2.500 s'], [], ''),
        ('r3: the same without time lines and finish marker', dict(table=base_tab[:-1], time_lines=False, finished=False), [],
         ALLU, ['no time lines in the log'], [], ''),
        ('r3: one audit line of a binary without floor_bal_z', dict(legacy=True, cad=2.5), [], {'F1': U, 'F5 height': U}, [], ['Traceback'], ''),
        ('r3: no column extents, lower window end at 1 d', dict(zlo_d=1.0, extents=False), [], ALLU, ['column extents not in the log'], [], ''),
        ('r3: 30 ms audit (every window qualifies)', dict(cad=0.03), [], {'F1': 'PASS', 'F2': 'FAIL'}, [], ['unmeasurable'], ''),
        ('r3: only the final rows in the log', dict(mutate=lambda L: [l for l in L if not l.startswith('EL_BED_AUDIT')
                                                                     or float(l.split()[2]) > 2.0]), [],
         {'F1': U, 'F5 count': U, 'F5 height': 'PASS'}, [], [], ''),
        ('r3: column translated down by 20 mm, 60 mm late', dict(z0=-0.020, zmax_fn=hi_after), [], {'F5 height': 'FAIL', 'F3': 'PASS'},
         ['60.00 mm above the inlet plane'], [], ''),
        ('r4: redundant knot inside the 0.20 m/s hold (equivalent table)', dict(table=tab_knot), [], REF, [], [], 'same'),
        ('r4: lowest hold subdivided into 80 ms pieces, 60 mm in it', dict(table=tab_sub, zmax_fn=hi_after), [],
         {'F5 height': 'FAIL'}, ['lowest velocity 0.050 m/s (2.020-2.500 s'], [], ''),
        ('r4: one audit row of a fluidised window removed', dict(mutate=_drop(0.85)), [],
         {'F1': U, 'F2': U, 'F3': U, 'F5 count': U, 'F5 height': 'PASS'}, ['missing or rejected'], [], ''),
        ('r4: dp_dyn NaN on one fluidised row (pressure identity)', dict(mutate=_set_field(0.85, 'dp_hydro_dyn', 2, 'NaN')), [],
         {'F1': U, 'F5 count': U, 'F5 height': 'PASS'}, ['pressure identity'], [], ''),
        ('r4: one fixed step with a raised gradient, one of its rows removed',
         dict(grad_fn=lambda t, dp, u: dp / 0.1 * (1.5 if 1.5 < t <= 2.0 else 1.0), mutate=_drop(1.8)), [], {'F2': U}, [], [], ''),
        ('r4: run stopped at 2.15 s inside the lowest hold (MaxNumStep)', dict(t_end=2.15), [], ALLU,
         ['the run ended at 2.1500 s, before the last table row'], [], ''),
        ('r4: every other audit row removed, 60 mm on a removed row',
         dict(zmax_fn=lambda t: 0.060 if abs(t - 2.41) < 1e-9 else 0.0241, mutate=_thin(2)), [], {'F5 height': U, 'F1': U}, [], [], ''),
        ('r4: seven identical copies of one audit row', dict(mutate=_dup(0.79, 7)), [], REF, [], [], ''),
        ('r4: conflicting copy of one audit row', dict(mutate=_conflict(0.79)), [], {'F1': U}, ['conflicting copies'], [], ''),
        ('r4: MaxSimTime 2.4999, the solver completes 2.5 s, 60 mm on the last row',
         dict(tmax=2.4999, zmax_fn=lambda t: 0.060 if t > 2.499 else 0.0241), [], {'F5 height': 'FAIL'}, [], [], ''),
        ('r4: Reload_SimPar extends MaxSimTime from 2.3 to 2.5 s, 60 mm at 2.4 s (witnessed)',
         dict(mutate=_reload_extension(2.3, 2.3), zmax_fn=lambda t: 0.060 if abs(t - 2.4) < 1e-9 else 0.0241), [],
         dict(ALLU, **{'F5 height': 'FAIL'}), ['MaxSimTime or MaxNumStep changed'], [], ''),
        ('r4: unfinished trailing EL_NEWTON_PAIR', dict(mutate=_append('EL_NEWTON_PAIR t= 2.5 step=')), [], REF, [], ['Traceback'], ''),
        ('r4: last line without its line end (run in progress)', dict(mutate=_append('EL_BED_AUDIT t=   2.502500E+00 step= 10010 p_in_out='),
                                                                      last_complete=False, finished=False), [], ALLU, [], ['Traceback'], ''),
        ('r4: malformed interior table row', dict(mutate=lambda L: [l.replace('3.00000000E-01', '3.0000000xE-01', 1)
                                                                    if l.lstrip().startswith('EL_INLET_TABLE') and ' 1.00000000E+00' in l else l for l in L]), [],
         ALLU, ['table echo incomplete'], ['Traceback'], ''),
        ('r4: particles lost on the last audit row (count and weight, witnessed)', dict(n_fn=lambda t: 9789 if t < 2.495 else 9691), [],
         {'F1': U, 'F5 count': 'FAIL'}, [], [], ''),
        ('r4: wenyu echo, default invocation', dict(drag='wenyu'), [], {'F1': 'PASS'}, ['wenyu U_mf at the whole-bed'], [], ''),
        ('r4: explicit --closure that differs from the echo', dict(), ['--closure', 'ergun'], REF,
         ['comparison with the command line: closure ergun (run gidaspow)'], [], ''),
        ('r4: no width echo, deck 3.0d0 with a comment', dict(width_echo=False, deck=deck_full + 'SimPar@ELKernelWidthFactor = 3.0d0 ! mapping width in diameters\n'),
         [], {'F2': U}, ['mapping support 3 d [deck]'], [], ''),
        ('r4: no width echo, unreadable deck value', dict(width_echo=False, deck=deck_full + 'SimPar@ELKernelWidthFactor = wide\n'),
         [], ALLU, ['value unreadable'], [], ''),
        ('r4: NaN width echo', dict(width='NaN'), [], ALLU, ['mapping support not a positive finite number'], [], ''),
        ('r4: centre exactly at mid-height', dict(zmax_fn=lambda t: 0.050 if t > 2.3 else 0.0241), [], {'F5 height': 'PASS'},
         ['50.00 mm above the inlet plane', 'not above'], [], ''),
        ('r4: no echo, steps file with two 0.05 m/s rows, 60 mm in the second',
         dict(echo_table=False, protocol=[(0.5, 0.3), (1.0, 0.2), (1.5, 0.05), (2.0, 0.1), (2.5, 0.05)],
              steps_file=[(0.5, 1.0, 0.3), (1.0, 1.5, 0.2), (1.5, 2.0, 0.05), (2.0, 2.5, 0.1), (2.5, 3.0, 0.05)],
              zmax_fn=lambda t: 0.060 if t > 2.8 else 0.0241), ['--steps', S], ALLU,
         ['no official gate without the table echo', 'steps-file summary, NOT gated'], [], 'r4 (fallback): labelled FAIL'),
        ('r4: column translated by +10 m (7-digit coordinates)', dict(z0=10.0), [], REF, [], [], ''),
        ('r4: --d 0 refused at parsing', dict(), ['--d', '0'], {}, ['--d must be finite and positive'], ['ZeroDivision'], ''),
        ('r4: --mu 0 refused at parsing', dict(), ['--closure', 'wenyu', '--mu', '0'], {}, ['--mu must be finite and positive'], ['ZeroDivision'], ''),
        # ---- review r5 -----------------------------------------------------------------------------------
        ('r5: clean stop at 1.09 s (MaxNumStep), 59.5 mm in the preparation', dict(t_end=1.09, zmax_fn=lambda t: 0.0595 if t < 0.5 else 0.0241), [],
         ALLU, ['the run ended at 1.0900 s, before the last table row'], [], 'r5: F5 FAIL on the 0.30 m/s preparation'),
        ('r5: clean stop at 2.11 s, inside the lowest hold', dict(t_end=2.11), [], ALLU, ['before the last table row'], [],
         'r5: official F2 / F3 / F5 over the remaining steps'),
        ('r5: segment from 1.0 s with steps renumbered from 1', dict(mutate=lambda L: _renumber(L, 4000)), [], ALLU,
         ['the clock starts at 1.0000 s, the protocol (first table row) at 0.0000 s'], [], 'r5: F1 over the remaining steps'),
        ('r5: Reload_SimPar resets the time from 1.75 to 1.25 s, unchanged dt, 60 mm late',
         dict(clock_fn=reset, t_end=3.0, zmax_fn=lambda t: 0.06 if t > 2.0 else 0.0241), [], ALLU, ['the time line of step'], [],
         'r5: F5 height PASS at 0.10 m/s'),
        ('r5: repeated fixed velocity 0.30 -> 0.04 -> 0.30 -> 0.04', dict(protocol=rep_proto), [], {'F3': U, 'F1': 'PASS'},
         ['not identified (rank-deficient: 1 distinct nonzero fixed velocity'], [], 'r5: U_mf 0.5000, ratio 2.6490 FAIL'),
        ('r5: no fluidised step (0.14 .. 0.08 m/s)', dict(protocol=nofl_proto), [], {'F1': U, 'F2': U, 'F3': U},
         ['no step with dp* > 0.95 (no measured plateau'], [], 'r5: F2 FAIL, F3 PASS from the nominal plateau'),
        ('r5: centre 0.05000002 m (above mid-height by less than the print uncertainty)',
         dict(zmax_fn=lambda t: 0.05000002 if t > 2.3 else 0.0241), [], {'F5 height': U}, ['not decidable at print precision'], [], 'r5: PASS'),
        ('r5: centre 0.0501 m (above the band)', dict(zmax_fn=lambda t: 0.0501 if t > 2.3 else 0.0241), [], {'F5 height': 'FAIL'}, [], [], ''),
        ('r5: dt changed late, 60 mm in the lowest hold (witnessed)', dict(clock_fn=dtchg, zmax_fn=lambda t: 0.06 if 2.2 < t < 2.26 else 0.0241), [],
         dict(ALLU, **{'F5 height': 'FAIL'}), ['time step not constant'], [], 'r5: F5 height unavailable (observation lost)'),
        ('r5: echoed mu 2.4e-5, fixed gradients from it', dict(mu=mu_e, grad_fn=grad_mu), [], {'F2': 'PASS'}, ['mu 2.4e-05'], [],
         'r5: F2 31.91 % FAIL with the default mu'),
        ('r5: --mu 1.9e-5 that differs from the run', dict(mu=mu_e, grad_fn=grad_mu), ['--mu', '1.9e-5'], {'F2': 'PASS'},
         ['comparison with the command line: mu 1.9e-05 (run 2.4e-05)'], [], 'r5: official FAIL with the argument'),
        ('r5: dp_dyn of one row +200 Pa, pressures unchanged', dict(mutate=_add_field(0.85, 'dp_hydro_dyn', 2, 200.0)), [],
         {'F1': U}, ['pressure identity dp_dyn = dp - dp_hydro violated'], [], 'r5: F1 PASS from the edited row'),
        ('r5: every count -1', dict(mutate=lambda L: [l.replace(' n= 9789 ', ' n= -1 ') for l in L]), [], ALLU,
         ['particle count N not resolved'], [], 'r5: F5 count PASS'),
        ('r5: 0.20 m/s hold 1.02-1.40 s, sub-window 2 would cross into the next step', dict(table=sub_tab, t_end=2.4), [], {},
         ['   n/a  '], ['DIFF'], 'r5: sub2 0.740 DIFF'),
        ('r5: no echo, steps file only, 60 mm in the 20 ms transition', dict(echo_table=False, steps_file=FULL_STEPS,
                                                                            zmax_fn=lambda t: 0.06 if abs(t - 2.01) < 1e-9 else 0.0241),
         ['--steps', S], ALLU, ['no official gate without the table echo'], [], 'r5: labelled 60 mm FAIL'),
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
    for part in m.group(1).replace(OUTSIDE, '').split(' | '):
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
            if 'PASS' in V.values() and 'NO OFFICIAL VERDICT' in out:
                bad.append('a PASS from a run outside the supported domain')
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
               'time:', 'Prop')


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


AUDIT_FIELDS = [('t', 0), ('step', 0), ('p_in_out', 0), ('p_in_out', 1), ('area_in_out', 0), ('area_in_out', 1),
                ('dp_hydro_dyn', 0), ('dp_hydro_dyn', 1), ('dp_hydro_dyn', 2), ('n', 0), ('weight_z_floor_z', 0),
                ('weight_z_floor_z', 1), ('z50_95_99_max', 0), ('z50_95_99_max', 1), ('z50_95_99_max', 2),
                ('z50_95_99_max', 3), ('eps_bed', 0), ('window', 0), ('window', 1), ('gradz_win', 0), ('floor_bal_z', 0),
                ('eps_kern_slip_win', 0), ('eps_kern_slip_win', 1), ('win_margin', 0), ('gradz_win_dyn', 0)]
# NaN in one of these makes the row a hole (malformed record, failed identity, count / weight change)
HOLE_FIELDS = {('t', 0), ('step', 0), ('p_in_out', 0), ('p_in_out', 1), ('dp_hydro_dyn', 0), ('dp_hydro_dyn', 1),
               ('dp_hydro_dyn', 2), ('n', 0), ('weight_z_floor_z', 0)}


def self_test_real(rundir, nvar=300, seed=20261122):
    """metamorphic tests on a real log with the expected outcome per gate (design note 7.4)"""
    rnd = random.Random(seed)
    with open(os.path.join(rundir, 'run_slurm.log'), errors='replace') as fh:
        full = fh.read().split('\n')
    if full and full[-1] == '':
        full.pop()
    deck = open(os.path.join(rundir, '_data', 'q2p1_param.dat'), errors='replace').read()
    base = _filtered(full)
    ref_text, ref_V, ev0 = run_lines(base, deck, rundir=rundir)
    full_text, _, _ = run_lines(full, deck, rundir=rundir)
    print(f"real log {rundir}: {len(full)} lines, {len(base)} read by the tool; intact verdicts: "
          + ' | '.join(f'{g} {ref_V[g]}' for g in GATES))
    nfail = 0
    if full_text != ref_text:
        nfail += 1
        print('FAIL: the full log and the lines the tool reads give different output')
    table, holds, clock0 = ev0['knots'], ev0['holds'], ev0['clock']
    windows = [(x['t1'] - EVAL_LEN, x['t1']) for x in ev0['steps']]
    fixed = [(x['t1'] - EVAL_LEN, x['t1']) for x, c in zip(ev0['steps'], ev0['classes']) if c == 'fixed']
    lowest = [(h['a'], h['b']) for h in holds if h['u'] == min(g['u'] for g in holds)] if holds else []
    old_gw = not all(r['has_gw_dyn'] for r in clock0.sorted)
    audit_idx = [i for i, l in enumerate(base) if l.startswith('EL_BED_AUDIT')]
    time_idx = {int(PAT_TIME.match(l).group(3)): i for i, l in enumerate(base) if PAT_TIME.match(l)}
    marker = max((i for i, l in enumerate(base) if FINISH in l), default=len(base) - 1)
    dtc = clock0.dt
    tol = 1e-3 * dtc
    inside = lambda t, ivs: any(a - tol <= t <= b + tol for a, b in ivs)   # noqa: E731
    t_of = lambda i: num(base[i].split()[2])                              # noqa: E731
    ALLU = {g: 'unavailable' for g in GATES}

    def hole_expect(ts):
        """expected verdicts when a row of each time in ts becomes a hole"""
        E = {'F5 count': 'FAIL' if ref_V['F5 count'] == 'FAIL' else 'unavailable'}
        any_win = any(inside(t, windows) for t in ts)
        for g in ('F1', 'F2', 'F3'):
            E[g] = 'unavailable' if any_win else ref_V[g]
        E['F5 height'] = ('unavailable' if any(inside(t, lowest) for t in ts) and ref_V['F5 height'] != 'FAIL'
                          else ref_V['F5 height'])
        return E

    def field_expect(t, key, idx):
        if (key, idx) in HOLE_FIELDS:
            return hole_expect([t])
        E = dict(ref_V)
        if (key, idx) == ('area_in_out', 0) and inside(t, windows):
            E.update(F1='unavailable', F2='unavailable', F3='unavailable')
        if inside(t, fixed) and E['F2'] != 'unavailable':
            f2 = {('z50_95_99_max', 2), ('window', 0), ('window', 1), ('eps_kern_slip_win', 0), ('win_margin', 0),
                  ('gradz_win', 0) if old_gw else ('gradz_win_dyn', 0)}
            if (key, idx) in f2:
                E['F2'] = 'unavailable'
            if (key, idx) == ('z50_95_99_max', 2):
                E['F3'] = 'unavailable'
        if (key, idx) == ('z50_95_99_max', 3) and inside(t, lowest) and ref_V['F5 height'] != 'FAIL':
            E['F5 height'] = 'unavailable'
        return E

    # ---- the transforms applicable to this input -----------------------------------------------------
    inner = [h for h in holds if h['b'] - h['a'] > 0.01]
    ramps = [(p, q) for p, q in zip(table, table[1:]) if p[1] != q[1]]
    eq_kinds = ['trailing record', 'row duplicated']
    if inner:
        eq_kinds += ['knot in a hold', 'subdivided hold']
    if len(table) > 1 and table[-1][1] == table[-2][1]:
        eq_kinds.append('last knot removed')
    if ramps:
        eq_kinds.append('knot in a ramp')
    dg_kinds = ['one row deleted', 'rows deleted', 'conflicting copy', 'field NaN', 'log truncated', 'audit thinned',
                'count perturbed', 'weight perturbed', 'finish marker removed']
    if lowest and ref_V['F5 height'] == 'PASS':
        dg_kinds.append('witness kept, other rows deleted')
    hist_kinds = ['suffix renumbered', 'time reset (unchanged dt)', 'dt changed']
    if table and any(t > 0.25 * clock0.t_last for t, _ in table):
        hist_kinds.append('clean early stop')
    if sum(1 for c in ev0['classes'] if c == 'fixed') >= 2:
        hist_kinds.append('repeated fixed velocity')
    if ramps:
        hist_kinds.append('table without ramps')
    skipped = (sorted({'knot in a hold', 'subdivided hold', 'last knot removed', 'knot in a ramp'} - set(eq_kinds))
               + sorted({'witness kept, other rows deleted'} - set(dg_kinds))
               + sorted({'clean early stop', 'repeated fixed velocity', 'table without ramps'} - set(hist_kinds)))
    print('transforms applicable: equivalence %s; degrading %s; run history %s; skipped (not applicable to this input): %s'
          % (', '.join(eq_kinds), ', '.join(dg_kinds), ', '.join(hist_kinds), ', '.join(skipped) or 'none'))

    def equiv(kind):
        if kind == 'knot in a hold':
            h = rnd.choice(inner)
            t = round(rnd.uniform(h['a'] + 0.001, h['b'] - 0.001), 6)
            return _table_edit(base, sorted(set(table + [(t, h['u'])]))), True
        if kind == 'subdivided hold':
            h = rnd.choice(inner)
            k = rnd.randint(2, 6)
            extra = [(round(h['a'] + (h['b'] - h['a']) * j / k, 6), h['u']) for j in range(1, k)]
            return _table_edit(base, sorted(set(table + extra))), True
        if kind == 'last knot removed':
            return _table_edit(base, table[:-1]), True
        if kind == 'knot in a ramp':
            (ta, ua), (tb, ub) = rnd.choice(ramps)
            t = 0.5 * (ta + tb)
            return _table_edit(base, sorted(table + [(t, ua + (ub - ua) * (t - ta) / (tb - ta))])), True
        if kind == 'row duplicated':
            L = list(base)
            i = rnd.choice(audit_idx)
            L.insert(i + 1, L[i])
            return L, True
        tail = rnd.choice(['EL_NEWTON_PAIR t= 7.5 step=', 'EL_BED_AUDIT t=   7.502500E+00 step= 30010 p_in_out=',
                           'EL_FORCE_BUDGET step= 30001 drag=', 'time:  0.7500D+01/', '  EL_INLET_TABLE   7.6'])
        return base + [tail], rnd.random() < 0.5

    def degrade(kind):
        L = list(base)
        if kind == 'one row deleted':
            i = rnd.choice(audit_idx)
            del L[i]
            return L, True, hole_expect([t_of(i)])
        if kind == 'rows deleted':
            ii = rnd.sample(audit_idx, min(len(audit_idx), rnd.randint(2, 20)))
            ts = [t_of(i) for i in ii]
            for i in sorted(ii, reverse=True):
                del L[i]
            return L, True, hole_expect(ts)
        if kind == 'conflicting copy':
            i = rnd.choice(audit_idx)
            key, idx = rnd.choice([f for f in AUDIT_FIELDS if f[0] not in ('t', 'step', 'n')])
            tok = L[i].split()
            j = tok.index(key + '=') + 1 + idx
            v = num(tok[j])
            tok[j] = _e(v * 1.001 + 1e-9) if finite(v) else _e(1.0)
            L.insert(i + 1, ' '.join(tok))
            return L, True, hole_expect([t_of(i)])
        if kind == 'field NaN':
            i = rnd.choice(audit_idx)
            key, idx = rnd.choice(AUDIT_FIELDS)
            t = t_of(i)
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
            return [l for i, l in enumerate(L) if i not in drop], True, hole_expect([t_of(i) for i in drop])
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
            return L, True, hole_expect([t_of(i)])
        if kind == 'finish marker removed':
            return [l for l in L if FINISH not in l], True, dict(ALLU)
        # witness kept: a centre above mid-column on one row of the lowest hold survives deleting other rows
        w = rnd.choice([i for i in audit_idx if inside(t_of(i), lowest)])
        L[w] = _set_tok(L[w], 'z50_95_99_max', 3, _e(ev0['setup']['zmin'] + 0.06))
        others = rnd.sample([i for i in audit_idx if i != w], rnd.randint(1, 20))
        ts = [t_of(i) for i in others]
        for i in sorted(others, reverse=True):
            del L[i]
        E = hole_expect(ts)
        E['F5 height'] = 'FAIL'
        return L, True, E

    def history(kind):
        L = list(base)
        K = max(time_idx)
        if kind == 'clean early stop':
            last_t = max(t for t, _ in table)
            kmax = max(1, int((last_t - clock0.t0) / dtc) - 1)
            k = rnd.randrange(max(1, kmax // 10), kmax)
            return L[:time_idx[k] + 1] + [' PP3D_LES has successfully finished. '], 'before the last table row'
        if kind == 'suffix renumbered':
            f = clock0.freq or 1
            k0 = (rnd.randrange(int(0.1 * K), int(0.6 * K)) // f) * f
            return _renumber(L, k0), 'the clock starts at'
        if kind == 'time reset (unchanged dt)':
            k1 = rnd.randrange(int(0.2 * K), int(0.8 * K))
            sh = rnd.choice([0.25, 0.5, 1.0]) * max(clock0.t_last, 0.1)
            out = []
            for l in L:
                mm = PAT_TIME.match(l)
                if mm and int(mm.group(3)) > k1:
                    l = _time_line(l, num(mm.group(1)) - sh)
                elif l.startswith(('EL_BED_AUDIT', 'EL_NEWTON_PAIR')) and int(re.search(r'step=\s*(\d+)', l).group(1)) > k1:
                    tok = l.split()
                    tok[2] = _e(num(tok[2]) - sh)
                    l = ' '.join(tok)
                out.append(l)
            return out, 'the time line of step'
        if kind == 'dt changed':
            # a TimeStep reload to dt / 2 after step k1: the run stays inside the time span of the intact run, so
            # the transform creates no observation (a longer terminal hold) that the intact run does not have
            k1 = rnd.randrange(int(0.2 * K), int(0.8 * K))
            out = []
            for l in L:
                mm = PAT_TIME.match(l)
                if mm and int(mm.group(3)) > k1:
                    k = int(mm.group(3))
                    l = _time_line(l, clock0.t0 + k1 * dtc + (k - k1) * 0.5 * dtc, 0.5 * dtc)
                elif l.startswith(('EL_BED_AUDIT', 'EL_NEWTON_PAIR')):
                    k = int(re.search(r'step=\s*(\d+)', l).group(1))
                    if k > k1:
                        tok = l.split()
                        tok[2] = _e(clock0.t0 + k1 * dtc + (k - k1) * 0.5 * dtc)
                        l = ' '.join(tok)
                out.append(l)
            return out, 'time step not constant'
        if kind == 'repeated fixed velocity':
            fx = [x for x, c in zip(ev0['steps'], ev0['classes']) if c == 'fixed']
            uf = fx[0]['U']
            new = [(t, uf) if any(x['hold']['a'] - 1e-9 <= t <= x['hold']['b'] + 1e-9 for x in fx) else (t, u)
                   for t, u in table]
            out = []
            for j, p in enumerate(new):
                out.append(p)
                if j + 1 < len(new) and p[1] == uf and new[j + 1][1] == uf and new[j + 1][0] - p[0] < 0.05:
                    out.append((0.5 * (p[0] + new[j + 1][0]), uf * 0.99))     # a 20 ms dip: a transition, not a step
            return _table_edit(L, out), 'rank-deficient'
        rows = list(table)
        for j in range(1, len(rows)):
            if rows[j][1] != rows[j - 1][1] and rows[j][0] - rows[j - 1][0] < 0.05:
                rows[j] = (rows[j - 1][0], rows[j][1])                       # the ramp collapsed onto its start
        return _table_edit(L, rows), 'not strictly increasing'

    stats = {}
    n_un = 0
    for v in range(nvar):
        cls = v % 3
        key = 'variant %d' % v
        E = None
        try:
            if cls == 0:
                kind = rnd.choice(eq_kinds)
                key = 'equivalence: ' + kind
                L, complete = equiv(kind)
                text, V, _ = run_lines(L, deck, last_complete=complete, rundir=rundir)
                ok = text == ref_text
            elif cls == 1:
                kind = rnd.choice(dg_kinds)
                key = 'degrading: ' + kind
                L, complete, E = degrade(kind)
                text, V, _ = run_lines(L, deck, last_complete=complete, rundir=rundir)
                ok = all(V.get(g) == E[g] for g in GATES)
            else:
                kind = rnd.choice(hist_kinds)
                key = 'run history: ' + kind
                L, reason = history(kind)
                text, V, _ = run_lines(L, deck, rundir=rundir)
                if kind == 'repeated fixed velocity':
                    ok = V.get('F3') == 'unavailable' and reason in text
                else:
                    ok = all(V.get(g) == 'unavailable' for g in GATES) and reason in text
        except BaseException:                        # noqa: BLE001
            import traceback
            text, V, ok = 'Traceback ' + traceback.format_exc(), {}, False
        st = stats.setdefault(key, [0, 0, 0])
        st[0] += 1
        if V and any(V.get(g) == 'unavailable' for g in GATES) and cls:
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
    neq = sum(s[0] for k, s in stats.items() if k.startswith('equivalence'))
    ndg = sum(s[0] for k, s in stats.items() if k.startswith('degrading'))
    print(f"{nvar} variants ({neq} equivalence, {ndg} degrading, {nvar - neq - ndg} run history); non-equivalence variants "
          f"with at least one verdict unavailable: {n_un}; failures: {nfail}")
    for k in sorted(stats):
        print(f"   {k:46s} {stats[k][0]:4d} run, {stats[k][1]:3d} failed, {stats[k][2]:3d} with a verdict unavailable")
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
