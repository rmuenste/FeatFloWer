#!/usr/bin/env python3
"""Momentum balance of a q2p1_el_pipeflow run with SimPar@ELMomentumBalanceAudit = Yes (case 11).

Reads the EL_MOM_BALANCE_{P,GA,GD,B,C} lines (one set per macro step; meaning and sign of every
field in source/src_el/el_momentum_balance.f90) from RUNDIR/run_slurm.log (or a log file) and prints,
for a time window and for the whole log, the step means of the z components as mean force [N] and in
per cent of the buoyant bed weight W (W/A = 126.408 Pa on 15 mm x 15 mm unless given):

  particles   drag, pressure and viscous force armed in pe, hydro = H_p, weight, dp/dt, floor force
              from the particle balance, pe wall impulse (raw and relaxed), particle closure residual
  gas         every term of the discrete gas momentum equation over all rows (A), the Dirichlet rows
              (D) and the free rows (F = A - D), the derived hrz increments before / after the projection
  boundaries  plane pressures, (p_in - p_out) A, momentum flux and viscous traction of inlet and outlet,
              the Dirichlet-row constraint impulse J_D
  closures    gas balance on the free rows (solver tolerance), its round-off remainder, projection
              identity, J_D identity, particle balance
  F1 residual hydro - dpA (= H_p/W - D, the residual of RUNBOOK section 15) and its reconstruction by
              the ten terms rec_*; the pressure terms split once more into their all-row and Dirichlet
              parts

Only numpy and the standard library. usage:
  fl11_mombal.py RUNDIR|LOG [--window T0 T1] [--wa PA] [--area M2] [--csv OUT]
  fl11_mombal.py --self-test
"""
import argparse
import os
import re
import sys
import tempfile

import numpy as np

TAGS = ('EL_MOM_BALANCE_P', 'EL_MOM_BALANCE_GA', 'EL_MOM_BALANCE_GD', 'EL_MOM_BALANCE_B',
        'EL_MOM_BALANCE_C')
GAS = ('mom_start', 'mom_pred', 'mom_end', 'inertia_pred', 'inertia_proj', 'pgrad_old', 'mmp_old',
       'mmv_old', 'feedback', 'partner', 'rhs_other', 'dragimp_pred', 'visc_pred', 'conv_pred',
       'defect', 'rhs_full', 'pgrad_q', 'mmp_q', 'dragimp_end', 'pgrad_end', 'mmp_end', 'hrz_pred',
       'hrz_end')
REC = ('rec_pB', 'rec_pfp', 'rec_visc', 'rec_conv', 'rec_fbD', 'rec_hrz', 'rec_other', 'rec_solver',
       'rec_inert', 'rec_dep')
REC_TEXT = {
    'rec_pB': 'unweighted pressure term on free rows - (p_in - p_out) A',
    'rec_pfp': 'particle pressure force + model-A solid-pressure term (free rows)',
    'rec_visc': 'particle viscous force + gas viscous terms (free rows)',
    'rec_conv': 'gas convection on the free rows (momentum flux)',
    'rec_fbD': 'particle reaction deposited on Dirichlet (inlet) rows, lost',
    'rec_hrz': 'implicit hrz drag increment -B (u~ - u^n), free rows',
    'rec_other': 'other right-hand-side terms (body forces)',
    'rec_solver': 'residual of the momentum solve (free rows)',
    'rec_inert': 'gas inertia (free rows)',
    'rec_dep': 'armed drag + lift + deposited source (deposit closure)',
}
WA_DEFAULT = 126.408
AREA_DEFAULT = 0.015*0.015
NUMTOK = re.compile(r'^[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?$|^[-+]?nan$|^[-+]?inf', re.I)


def parse_line(line):
    """'TAG t= .. step= .. key= v v v ...' -> (tag, dict key -> list of floats or str)."""
    toks = line.split()
    tag, out, key = toks[0], {}, None
    for tok in toks[1:]:
        if tok.endswith('='):
            key = tok[:-1]
            out[key] = []
        elif key is not None:
            out[key].append(float(tok.replace('D', 'E').replace('d', 'e')) if NUMTOK.match(tok) else tok)
    return tag, out


def load(path):
    """step -> {tag: fields}; only complete steps (all five lines) are kept."""
    steps = {}
    with open(path, errors='replace') as f:
        for line in f:
            s = line.lstrip()
            if not s.startswith('EL_MOM_BALANCE_'):
                continue
            tag, d = parse_line(s)
            if tag not in TAGS or 'step' not in d:
                continue
            steps.setdefault(int(d['step'][0]), {})[tag] = d
    return {k: v for k, v in steps.items() if len(v) == len(TAGS)}


def series(steps):
    """Arrays (z components unless named otherwise) over the sorted steps."""
    ks = sorted(steps)
    s = {'step': np.array(ks, dtype=float), 't': np.array([steps[k]['EL_MOM_BALANCE_P']['t'][0] for k in ks])}

    def vec(tag, key, comp=2):
        return np.array([steps[k][tag][key][comp] for k in ks], dtype=float)
    for key in ('drag_expl', 'drag_impl', 'pressure', 'viscous', 'lift', 'weight', 'hydro', 'fb_used',
                'deposit', 'dp_dt', 'expected_dt', 'floor_bal', 'wall_contact_raw', 'wall_lub_raw',
                'res_p_raw', 'res_p'):
        s['P.' + key] = vec('EL_MOM_BALANCE_P', key)
    s['P.relax'] = np.array([steps[k]['EL_MOM_BALANCE_P']['relax'][0] for k in ks])
    for key in GAS:
        s['A.' + key] = vec('EL_MOM_BALANCE_GA', key)
        s['D.' + key] = vec('EL_MOM_BALANCE_GD', key)
        s['F.' + key] = s['A.' + key] - s['D.' + key]
    b = 'EL_MOM_BALANCE_B'
    s['B.p_in'] = vec(b, 'p_in_out', 0)
    s['B.p_out'] = vec(b, 'p_in_out', 1)
    s['B.area_in'] = vec(b, 'area_in_out', 0)
    for key in ('dpA', 'flux_in', 'flux_out', 'trac_in', 'trac_out'):
        s['B.' + key] = vec(b, key)
    # optional fields (absent in logs of the first audit binary): inlet / outlet gas fraction, inlet-row
    # reaction by dof class
    if all('eps_in_out' in steps[k][b] for k in ks):
        s['B.eps_in'] = vec(b, 'eps_in_out', 0)
        s['B.eps_out'] = vec(b, 'eps_in_out', 1)
        for i, c in enumerate(('vertex', 'edge', 'face')):
            s['B.fbD_' + c] = vec(b, 'fbD_class', i)
    for key in ('res_gas_F', 'defect_F', 'gas_roundoff', 'proj_F', 'proj_D', 'jd', 'jd_check',
                'hydro_dpA', 'rec_sum') + REC:
        s['C.' + key] = vec('EL_MOM_BALANCE_C', key)
    if 'B.eps_in' in s:
        # split of the pressure pair rec_pB + rec_pfp with the inlet-plane gas fraction (continuous
        # identity -sum_all fp(p) = -P_p + p_in (1 - eps_in) A):
        #   press_samp  = P_p + mmp_end(A) - p_in (1 - eps_in) A: particle pressure force against the
        #                 solid-fraction pressure term of the gas (model-A sampling / consistency)
        #   press_inlet = the rest = -[pgrad_end + mmp_end](D) + pgrad_end(A) - (eps_in p_in - p_out) A:
        #                 inlet-row representation of the eps-weighted pressure force against its plane value
        s['C.press_samp'] = s['P.pressure'] + s['A.mmp_end'] - s['B.p_in']*(1.0 - s['B.eps_in'])*s['B.area_in']
        s['C.press_inlet'] = s['C.rec_pB'] + s['C.rec_pfp'] - s['C.press_samp']
    return s


def window_mask(s, t0, t1):
    eps = 1e-9
    return (s['t'] >= t0 - eps) & (s['t'] <= t1 + eps)


def report(s, mask, W, label, out=sys.stdout):
    n = int(mask.sum())
    if n == 0:
        print(f'\n== {label}: no step', file=out)
        return {}
    m = {k: float(np.mean(v[mask])) for k, v in s.items() if k not in ('step', 't')}
    mx = {k: float(np.max(np.abs(v[mask]))) for k, v in s.items() if k not in ('step', 't')}
    pc = lambda v: 100.0*v/W

    def row(name, key, text=''):
        print(f'  {name:<34s} {m[key]: .6e} N {pc(m[key]): 9.4f} %W  {text}', file=out)
    print(f'\n== {label}: {n} steps, t {s["t"][mask][0]:.6f} .. {s["t"][mask][-1]:.6f} s, '
          f'W = {W:.6e} N (measured buoyant weight {-m["P.weight"]:.6e} N)', file=out)
    print('-- particles (z, mean force over the step)', file=out)
    for key, text in (('P.drag_expl', 'closure drag at the start of the step'),
                      ('P.drag_impl', 'drag armed (sub-cycled semi-implicit impulse / dt)'),
                      ('P.pressure', '-V grad p'), ('P.viscous', 'V div tau'), ('P.lift', ''),
                      ('P.hydro', 'H_p = drag_impl + pressure + viscous + lift'),
                      ('P.weight', 'buoyant weight'), ('P.dp_dt', 'particle momentum change / dt'),
                      ('P.floor_bal', 'wall force from the particle balance'),
                      ('P.wall_contact_raw', 'pe wall contact force, unrelaxed'),
                      ('P.res_p', 'particle closure with the relaxed wall impulse'),
                      ('P.res_p_raw', 'particle closure with the raw wall impulse'),
                      ('P.fb_used', 'feedback the particles deposit'), ('P.deposit', 'deposited source')):
        row(key, key, text)
    print('-- gas (z): all rows A / Dirichlet rows D / free rows F, mean force [N] and %W of F', file=out)
    for key in GAS:
        if key.startswith('mom_'):
            print(f'  {key:<14s} A {m["A."+key]: .6e} D {m["D."+key]: .6e}   [kg m/s]', file=out)
            continue
        print(f'  {key:<14s} A {m["A."+key]: .6e} D {m["D."+key]: .6e} F {m["F."+key]: .6e} '
              f'({pc(m["F."+key]): 9.4f} %W)', file=out)
    print('-- boundaries', file=out)
    print(f'  p_in {m["B.p_in"]:.6e} Pa  p_out {m["B.p_out"]:.6e} Pa  area_in {m["B.area_in"]:.6e} m2', file=out)
    for key, text in (('B.dpA', '(p_in - p_out) A'),
                      ('B.flux_in', 'momentum flux out of the domain through the inlet (gas gets -flux)'),
                      ('B.flux_out', 'momentum flux out of the domain through the outlet'),
                      ('B.trac_in', 'viscous traction on the gas at the inlet'),
                      ('B.trac_out', 'viscous traction on the gas at the outlet'),
                      ('C.jd', 'constraint impulse of the inlet rows J_D / dt')):
        row(key, key, text)
    print('-- closures (mean / max |.| over the steps, N)', file=out)
    for key, text in (('C.res_gas_F', 'free-row gas balance (= -defect_F: solver tolerance)'),
                      ('C.defect_F', 'momentum-solve residual, free rows'),
                      ('C.gas_roundoff', 'res_gas_F + defect_F (round-off)'),
                      ('C.proj_F', 'projection identity, free rows (round-off)'),
                      ('C.proj_D', 'projection on Dirichlet rows (exactly 0)'),
                      ('C.jd_check', 'J_D identity (round-off)'),
                      ('P.res_p', 'particle balance (pe bookkeeping)')):
        print(f'  {key:<16s} mean {m[key]: .3e}  max {mx[key]: .3e}  ({pc(m[key]): .4f} %W)  {text}', file=out)
    print('-- F1 residual hydro - dpA and its reconstruction (z, mean force, %W)', file=out)
    D = m['B.dpA']/W
    print(f'  D = dpA/W {D:.6f}   H_p/W {m["P.hydro"]/W:.6f}   floor/W {m["P.floor_bal"]/W:.6f}   '
          f'-inertia_p/W {-m["P.dp_dt"]/W:.6f}', file=out)
    row('hydro - dpA (= residual R)', 'C.hydro_dpA')
    for key in REC:
        row(key, 'C.' + key, REC_TEXT[key])
    row('sum of the ten', 'C.rec_sum')
    print(f'  {"hydro - dpA - sum":<34s} {m["C.hydro_dpA"] - m["C.rec_sum"]: .3e} N', file=out)
    print('-- the pressure terms split (z, %W)', file=out)
    pB_all = m['A.pgrad_end']
    pB_D = -m['D.pgrad_end'] - m['B.dpA']
    print(f'  rec_pB = pgrad_end(A) {pc(pB_all): .4f} + [-pgrad_end(D) - dpA] {pc(pB_D): .4f}  '
          f'(pgrad_end(D) {m["D.pgrad_end"]:.6e} N = -(inlet-row pressure force))', file=out)
    print(f'  rec_pfp = [pressure + mmp_end(A)] {pc(m["P.pressure"] + m["A.mmp_end"]): .4f} '
          f'+ [-mmp_end(D)] {pc(-m["D.mmp_end"]): .4f}', file=out)
    print(f'  rec_pB + rec_pfp = gas pressure term on free rows (pgrad_end + mmp_end, F) + particle pressure '
          f'force - dpA: {pc(m["C.rec_pB"] + m["C.rec_pfp"]): .4f}', file=out)
    if 'C.press_samp' in m:
        print(f'  pressure pair = press_samp {pc(m["C.press_samp"]): .4f} (particle pressure force vs the gas solid-fraction '
              f'pressure term) + press_inlet {pc(m["C.press_inlet"]): .4f} (inlet rows vs eps_in p_in A; eps_in '
              f'{m["B.eps_in"]:.5f})', file=out)
        print(f'  fbD by inlet dof class (z, %W): vertex {pc(-m["B.fbD_vertex"]): .4f} edge {pc(-m["B.fbD_edge"]): .4f} '
              f'face {pc(-m["B.fbD_face"]): .4f}', file=out)
    print(f'  coupling terms together rec_fbD + rec_hrz + rec_dep {pc(m["C.rec_fbD"] + m["C.rec_hrz"] + m["C.rec_dep"]): .4f} '
          f'(the split between them depends on how the explicit drag partner is formed)', file=out)
    print(f'  hrz increment: before projection A {pc(m["A.hrz_pred"]): .4f} F {pc(m["F.hrz_pred"]): .4f}; '
          f'end of step A {pc(m["A.hrz_end"]): .4f} F {pc(m["F.hrz_end"]): .4f}', file=out)
    return m


def write_csv(s, path):
    keys = ['step', 't'] + [k for k in s if k not in ('step', 't')]
    with open(path, 'w') as f:
        f.write(','.join(keys) + '\n')
        for i in range(len(s['step'])):
            f.write(','.join(repr(float(s[k][i])) for k in keys) + '\n')


def resolve(path):
    if os.path.isdir(path):
        return os.path.join(path, 'run_slurm.log')
    return path


def analyse(path, t0, t1, wa, area, csv=None, out=sys.stdout):
    steps = load(resolve(path))
    if not steps:
        print('no complete EL_MOM_BALANCE step in', path, file=out)
        return None
    s = series(steps)
    W = wa*area
    print(f'{resolve(path)}: {len(steps)} complete steps', file=out)
    mw = report(s, window_mask(s, t0, t1), W, f'window {t0:g} .. {t1:g} s', out)
    mf = report(s, np.ones(len(s['t']), dtype=bool), W, 'whole log', out)
    if csv:
        write_csv(s, csv)
    return s, mw, mf


def synthetic_log(n, dt, rng):
    """A log with known step means: z values = base + noise with zero mean over the steps."""
    lines = []
    base = {'drag_impl': 0.011, 'pressure': 0.0145, 'viscous': -1e-7, 'weight': -0.0284418}
    for i in range(1, n + 1):
        t = i*dt
        head = f' t= {t:24.16E} step= {i}'
        sgn = 1.0 if i % 2 else -1.0
        vz = {k: v + sgn*1e-4 for k, v in base.items()}
        hydro = vz['drag_impl'] + vz['pressure'] + vz['viscous']

        def v3(z):
            return f'{0.0:24.15E}{0.0:24.15E}{z:24.15E}'
        p = ['EL_MOM_BALANCE_P' + head + ' unit= N']
        for key in ('drag_expl', 'drag_impl', 'pressure', 'viscous', 'lift', 'weight', 'hydro', 'fb_used',
                    'deposit', 'dp_dt', 'expected_dt', 'floor_bal', 'wall_contact_raw', 'wall_lub_raw',
                    'res_p_raw', 'res_p'):
            z = {'drag_impl': vz['drag_impl'], 'drag_expl': vz['drag_impl'], 'pressure': vz['pressure'],
                 'viscous': vz['viscous'], 'weight': vz['weight'], 'hydro': hydro,
                 'floor_bal': -(hydro + vz['weight'])}.get(key, 0.0)
            p.append(f' {key}= ' + v3(z))
        p.append(f' relax= {0.9:12.4E}')
        lines.append(''.join(p))
        for tag, scale in (('EL_MOM_BALANCE_GA', 1.0), ('EL_MOM_BALANCE_GD', 0.25)):
            g = [tag + head + ' unit= N(mom_*:kg_m/s)']
            for j, key in enumerate(GAS):
                g.append(f' {key}= ' + v3(scale*(j + 1)*1e-3 + sgn*1e-5))
            lines.append(''.join(g))
        b = ['EL_MOM_BALANCE_B' + head, f' p_in_out= {112.0 + sgn:24.15E}{0.0:24.15E}',
             f' area_in_out= {2.25e-4:24.15E}{2.25e-4:24.15E}', ' dpA= ' + v3(112.0*2.25e-4 + sgn*2.25e-4)]
        for key in ('flux_in', 'flux_out', 'trac_in', 'trac_out'):
            b.append(f' {key}= ' + v3(1e-6))
        lines.append(''.join(b))
        c = ['EL_MOM_BALANCE_C' + head]
        recv = {k: (j + 1)*1e-5 for j, k in enumerate(REC)}
        for key in ('res_gas_F', 'defect_F', 'gas_roundoff', 'proj_F', 'proj_D', 'jd', 'jd_check',
                    'hydro_dpA', 'rec_sum') + REC:
            z = recv.get(key, 0.0)
            if key in ('hydro_dpA', 'rec_sum'):
                z = sum(recv.values())
            c.append(f' {key}= ' + v3(z))
        lines.append(''.join(c))
        lines.append(f'XxX unrelated solver line {i}')
    return '\n'.join(lines) + '\n'


def self_test():
    rng = np.random.default_rng(1)
    ok = True
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, 'run_slurm.log')
        with open(path, 'w') as f:
            f.write(synthetic_log(40, 0.025, rng))
            f.write('EL_MOM_BALANCE_P t= 1.1 step= 41 unit= N drag_impl= 1 2 3\n')   # incomplete step
        steps = load(path)
        ok &= len(steps) == 40
        s = series(steps)
        out = open(os.devnull, 'w')
        mw = report(s, window_mask(s, 0.5, 1.0), WA_DEFAULT*AREA_DEFAULT, 'w', out)
        mf = report(s, np.ones(40, dtype=bool), WA_DEFAULT*AREA_DEFAULT, 'f', out)
        # 40 steps of alternating +-1e-4 (even count): exact means of the base values
        ok &= abs(mf['P.drag_impl'] - 0.011) < 1e-15 and abs(mf['P.weight'] + 0.0284418) < 1e-15
        # window 0.5 .. 1.0 s inclusive: steps 20 .. 40, 21 steps (11 odd, 10 even)
        ok &= int(window_mask(s, 0.5, 1.0).sum()) == 21
        ok &= abs(mw["P.drag_impl"] - (0.011 - 1e-4/21)) < 1e-14
        # free rows = all - Dirichlet
        ok &= abs(mf['F.feedback'] - 0.75*9e-3) < 1e-15
        ok &= abs(mf['C.rec_sum'] - mf['C.hydro_dpA']) < 1e-18
        ok &= abs(mf['B.dpA'] - 112.0*2.25e-4) < 1e-15
        ok &= abs(mf['B.p_in'] - 112.0) < 1e-12
        csv = os.path.join(d, 'x.csv')
        write_csv(s, csv)
        ok &= sum(1 for _ in open(csv)) == 41
    print('fl11_mombal self-test', 'PASSED' if ok else 'FAILED')
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('path', nargs='?', help='run directory (run_slurm.log) or log file')
    ap.add_argument('--window', nargs=2, type=float, default=(0.70, 1.00), metavar=('T0', 'T1'),
                    help='evaluation window in s, both ends inclusive (default 0.70 1.00)')
    ap.add_argument('--wa', type=float, default=WA_DEFAULT, help='buoyant bed weight per area W/A [Pa]')
    ap.add_argument('--area', type=float, default=AREA_DEFAULT, help='cross-section A [m2]')
    ap.add_argument('--csv', help='write the per-step series (z components) to this file')
    ap.add_argument('--self-test', action='store_true', help='run the self-test on a synthetic log')
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.path:
        ap.error('a run directory or log is needed')
    return 0 if analyse(a.path, a.window[0], a.window[1], a.wa, a.area, a.csv) else 1


if __name__ == '__main__':
    sys.exit(main())
