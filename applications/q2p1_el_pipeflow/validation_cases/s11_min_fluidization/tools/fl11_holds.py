#!/usr/bin/env python3
"""Fixed-bed holds of case 11 (experiment E2): window audit, momentum balance, packing drift.

For one or more run directories with SimPar@ELBedWindowAudit = Yes and SimPar@ELMomentumBalanceAudit
= Yes (lines EL_BED_WINDOW, EL_BED_VERTEX_MODE, EL_BED_AUDIT, EL_MOM_BALANCE_*; meanings in
source/src_el/el_bed_audit.f90 and el_momentum_balance.f90) prints

  time series  per 10 ms: dp* = dpA/W, floor/W, window gradients, the free-row gas balance residual
               and R = H_p - dpA (to judge the start-up transient)
  table 1      per window (m2d, m2.5d, interior) over the evaluation window: gradient G, closure
               gradient at the mapped gas fraction G_cl (fl11_closures.closure_grad at the inlet
               velocity, mean particle diameter), G/G_cl, the two drag means (per particle, N), mean
               slip s_z and gas fraction, pressure and viscous force per particle; the vertex-mode
               spread of the interior; window-to-window differences within a run and run-to-run
               differences at a common window (when two runs are given)
  table 2      the momentum balance of the evaluation window in %W and in % of the run's own dpA
               (H_p, dpA, floor, R and the reconstruction terms of fl11_mombal)
  drift        particle displacement between the seed file (RUNDIR/particles.xyz, line k = ip k)
               and the last logged positions, the bed top z99 + r and the floor force over time
  record       (--record DIR --record-window T0 T1) the record's EL_BED_AUDIT means over its hold:
               dp*, window gradient, kernel gas fraction, bed porosity, floor share

usage: fl11_holds.py RUNDIR [RUNDIR2] [--window T0 T1] [--u 0.12] [--record DIR --record-window 5.2 5.5]
       fl11_holds.py --self-test
"""
import argparse
import glob
import math
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fl11_closures  # noqa: E402
import fl11_mombal  # noqa: E402

WINDOWS = ('m2d', 'm2.5d', 'interior')
CLASSES = ('vertex', 'edge', 'face', 'centre')
RHOF, MU = 1.0, 1.8e-5


def keyed(tokens):
    out, key = {}, None
    for tok in tokens:
        if tok.endswith('='):
            key = tok[:-1]
            out[key] = []
        elif key is not None:
            try:
                out[key].append(float(tok))
            except ValueError:
                out[key].append(tok)
    return out


def load_windows(log):
    """{window: {key: array over the audit steps}}, vertex mode {class: (n, mean, std, min, max)
    arrays}, bed audit {key: arrays}"""
    win = {w: {} for w in WINDOWS}
    vm = {'t': []}
    bed = {k: [] for k in ('t', 'dp', 'area', 'weight', 'floor', 'z99', 'eps_bed', 'gradz_win', 'eps_kern', 'slip')}
    with open(log, errors='replace') as f:
        for line in f:
            if line.startswith('EL_BED_WINDOW '):
                d = keyed(line.split()[1:])
                w = win[d['win'][0]]
                for k in ('t', 'gradz', 'n', 'eps_geo', 'eps_kern', 'drag_mean', 'drag_of_means', 'pres_mean',
                          'visc_mean'):
                    w.setdefault(k, []).append(d[k][0])
                w.setdefault('slip_z', []).append(d['slip'][2])
                w.setdefault('z_lo', []).append(d['z_lo_hi'][0])
                w.setdefault('z_hi', []).append(d['z_lo_hi'][1])
            elif line.startswith('EL_BED_VERTEX_MODE '):
                d = keyed(line.split()[1:])
                vm['t'].append(d['t'][0])
                for c in CLASSES:
                    vm.setdefault(c, []).append(d[c][:5])
            elif line.startswith('EL_BED_AUDIT '):
                d = keyed(line.split()[1:])
                bed['t'].append(d['t'][0])
                bed['dp'].append(d['dp_hydro_dyn'][2])
                bed['area'].append(d['area_in_out'][0])
                bed['weight'].append(-d['weight_z_floor_z'][0])
                bed['floor'].append(d['floor_bal_z'][0])
                bed['z99'].append(d['z50_95_99_max'][2])
                bed['eps_bed'].append(d['eps_bed'][0])
                bed['gradz_win'].append(d['gradz_win'][0])
                bed['eps_kern'].append(d['eps_kern_slip_win'][0])
                bed['slip'].append(d['eps_kern_slip_win'][1])
    win = {w: {k: np.array(v, dtype=float) for k, v in d.items()} for w, d in win.items()}
    vm = {k: np.array(v, dtype=float) for k, v in vm.items()}
    bed = {k: np.array(v, dtype=float) for k, v in bed.items()}
    return win, vm, bed


def tmask(t, t0, t1):
    return (t >= t0 - 1e-6) & (t <= t1 + 1e-6)


def table1(name, win, vm, t0, t1, u, d, out):
    print(f'\n== {name}: windows, mean over the audits in {t0:g} .. {t1:g} s (U = {u} m/s, d = {d:.4e} m)', file=out)
    print(f'  {"window":9s} {"z_lo":>9s} {"z_hi":>9s} {"n":>6s} {"G [Pa/m]":>11s} {"eps_kern":>9s} {"G_cl":>11s} '
          f'{"G/G_cl":>8s} {"eps_geo":>8s} {"s_z":>10s} {"<B s_z> [N]":>12s} {"B(<>)<s_z>":>12s} {"ratio":>8s} '
          f'{"F_p [N]":>11s} {"F_v [N]":>10s}', file=out)
    res = {}
    for w in WINDOWS:
        a = win[w]
        if 't' not in a or len(a['t']) == 0:
            continue
        m = tmask(a['t'], t0, t1)
        if not m.any():
            continue
        mean = {k: float(np.mean(v[m])) for k, v in a.items()}
        gcl = fl11_closures.closure_grad(u, mean['eps_kern'], d, RHOF, MU)
        mean['G_cl'] = gcl
        mean['ratio'] = mean['gradz']/gcl
        mean['drag_ratio'] = mean['drag_mean']/mean['drag_of_means']
        res[w] = mean
        print(f'  {w:9s} {mean["z_lo"]:9.3e} {mean["z_hi"]:9.3e} {mean["n"]:6.0f} {mean["gradz"]:11.5e} '
              f'{mean["eps_kern"]:9.6f} {gcl:11.5e} {mean["ratio"]:8.5f} {mean["eps_geo"]:8.5f} '
              f'{mean["slip_z"]:10.4e} {mean["drag_mean"]:12.5e} {mean["drag_of_means"]:12.5e} '
              f'{mean["drag_ratio"]:8.5f} {mean["pres_mean"]:11.4e} {mean["visc_mean"]:10.3e}', file=out)
    if len(vm.get('t', [])):
        m = tmask(vm['t'], t0, t1)
        print('  vertex mode, interior window, w over the Q2 dofs per class (means over the audits):', file=out)
        for c in CLASSES:
            v = vm[c][m]
            print(f'    {c:7s} n {v[:, 0].mean():8.1f}  mean {v[:, 1].mean(): .5e}  std {v[:, 2].mean(): .5e}  '
                  f'min {v[:, 3].mean(): .5e}  max {v[:, 4].mean(): .5e} m/s', file=out)
            res['vm_' + c] = [float(x) for x in v.mean(axis=0)]
    for a, b in (('m2d', 'm2.5d'), ('m2.5d', 'interior'), ('m2d', 'interior')):
        if a in res and b in res:
            print(f'  {a} - {b}: G/G_cl {res[a]["ratio"] - res[b]["ratio"]: .5f}, eps_kern '
                  f'{res[a]["eps_kern"] - res[b]["eps_kern"]: .6f}, G {res[a]["gradz"] - res[b]["gradz"]: .4e} Pa/m',
                  file=out)
    return res


def mombal_table(rundir, t0, t1, wa, area, out):
    steps = fl11_mombal.load(os.path.join(rundir, 'run_slurm.log'))
    if not steps:
        print('  (no EL_MOM_BALANCE lines)', file=out)
        return None
    s = fl11_mombal.series(steps)
    m = fl11_mombal.window_mask(s, t0, t1)
    W = wa*area
    mean = {k: float(np.mean(v[m])) for k, v in s.items() if k not in ('step', 't')}
    dpa = mean['B.dpA']
    print(f'  momentum balance, {int(m.sum())} steps in {t0:g} .. {t1:g} s; z, %W and % of dpA:', file=out)
    rows = [('H_p', mean['P.hydro']), ('dpA', dpa), ('floor (balance)', mean['P.floor_bal']),
            ('R = H_p - dpA', mean['C.hydro_dpA'])]
    rows += [(k, mean['C.' + k]) for k in fl11_mombal.REC]
    rows += [('pressure pair (pB + pfp)', mean['C.rec_pB'] + mean['C.rec_pfp'])]
    rows += [('coupling fbD + hrz + dep', mean['C.rec_fbD'] + mean['C.rec_hrz'] + mean['C.rec_dep'])]
    if 'C.press_samp' in mean:
        rows += [('  of which press_samp', mean['C.press_samp']), ('  of which press_inlet', mean['C.press_inlet'])]
    for c in ('vertex', 'edge', 'face'):
        if 'B.fbD_' + c in mean:
            rows.append(('fbD on inlet ' + c + ' dofs', -mean['B.fbD_' + c]))
    for name, v in rows:
        print(f'    {name:26s} {v: .6e} N  {100*v/W: 9.4f} %W  {100*v/dpa: 9.4f} %dpA', file=out)
    if 'B.eps_in' in mean:
        print(f'    gas fraction on the inlet plane (transport porosity, area mean) {mean["B.eps_in"]:.6f}', file=out)
    print(f'    closures: max |gas_roundoff| {np.max(np.abs(s["C.gas_roundoff"][m])):.2e} N, mean defect_F '
          f'{mean["C.defect_F"]:.2e} N, max |res_p| {np.max(np.abs(s["P.res_p"][m])):.2e} N', file=out)
    return s, mean


def series_report(name, s, bed, win, wa, area, out, dt_bin=0.01):
    W = wa*area
    print(f'\n== {name}: time series (means per {dt_bin*1000:.0f} ms)', file=out)
    print(f'  {"t_end":>6s} {"dp*":>8s} {"floor/W":>8s} {"G m2.5d":>10s} {"G int":>10s} {"res_gas_F/W":>12s} '
          f'{"R/W %":>8s} {"z99 mm":>8s}', file=out)
    tend = float(s['t'].max()) if s is not None else float(bed['t'].max())
    nb = int(round(tend/dt_bin))
    for k in range(nb):
        a, b = k*dt_bin, (k + 1)*dt_bin
        mb = (bed['t'] > a + 1e-9) & (bed['t'] <= b + 1e-9)
        ms = (s['t'] > a + 1e-9) & (s['t'] <= b + 1e-9) if s is not None else None
        w25 = win['m2.5d']
        wi = win['interior']
        g25 = np.mean(w25['gradz'][(w25['t'] > a + 1e-9) & (w25['t'] <= b + 1e-9)]) if 't' in w25 else float('nan')
        gi = np.mean(wi['gradz'][(wi['t'] > a + 1e-9) & (wi['t'] <= b + 1e-9)]) if 't' in wi else float('nan')
        dps = np.mean(bed['dp'][mb]*bed['area'][mb]/W) if mb.any() else float('nan')
        fl = np.mean(bed['floor'][mb]/W) if mb.any() else float('nan')
        z99 = np.mean(bed['z99'][mb])*1e3 if mb.any() else float('nan')
        rg = np.mean(s['C.res_gas_F'][ms])/W if s is not None and ms.any() else float('nan')
        rr = 100*np.mean(s['C.hydro_dpA'][ms])/W if s is not None and ms.any() else float('nan')
        print(f'  {b:6.3f} {dps:8.5f} {fl:8.5f} {g25:10.4e} {gi:10.4e} {rg:12.3e} {rr:8.4f} {z99:8.4f}', file=out)


def drift(rundir, r, out):
    seed = os.path.join(rundir, 'particles.xyz')
    files = sorted(glob.glob(os.path.join(rundir, 'el_particles_rank*.log')))
    if not os.path.exists(seed) or not files:
        print('  drift: no seed file or particle logs', file=out)
        return None
    x0 = np.loadtxt(seed)
    # the last logged time over all files, then every particle at that time
    tl = max(float(line.split()[0]) for f in files for line in open(f) if not line.startswith('#'))
    pos = {}
    for f in files:
        for line in open(f):
            if line.startswith('#'):
                continue
            t = line.split()
            if abs(float(t[0]) - tl) <= 1e-9:
                pos[int(t[1])] = (float(t[8]), float(t[9]), float(t[10]))
    if len(pos) != len(x0):
        print(f'  drift: {len(pos)} particles at t = {tl}, {len(x0)} seeded', file=out)
        return None
    x1 = np.array([pos[i + 1] for i in range(len(x0))])
    dx = x1 - x0
    dx[:, 0] -= 0.015*np.round(dx[:, 0]/0.015)
    dx[:, 1] -= 0.015*np.round(dx[:, 1]/0.015)
    dd = np.sqrt((dx*dx).sum(axis=1))
    d = 2*r
    print(f'  drift seed -> t = {tl:.6f} s: |dx| mean {dd.mean():.3e} m ({dd.mean()/d:.4f} d), rms {np.sqrt((dd**2).mean()):.3e}, '
          f'p99 {np.percentile(dd, 99):.3e}, max {dd.max():.3e} m ({dd.max()/d:.3f} d); mean dz {dx[:, 2].mean():.3e} m; '
          f'particles moved > 0.1 d: {(dd > 0.1*d).sum()}', file=out)
    return dd


def record_means(rundir, t0, t1, wa, area, out):
    _, _, bed = load_windows(os.path.join(rundir, 'run_slurm.log'))
    m = tmask(bed['t'], t0, t1)
    W = wa*area
    print(f'  record {os.path.basename(rundir.rstrip("/"))} EL_BED_AUDIT {t0:g} .. {t1:g} s ({int(m.sum())} audits): '
          f'dp* {np.mean(bed["dp"][m]*bed["area"][m])/W:.5f}, gradz_win {np.mean(bed["gradz_win"][m]):.5e} Pa/m, '
          f'eps_kern {np.mean(bed["eps_kern"][m]):.6f}, eps_bed {np.mean(bed["eps_bed"][m]):.6f}, '
          f'slip {np.mean(bed["slip"][m]):.5e}, floor/W {np.mean(bed["floor"][m])/W:.5f}, z99 '
          f'{1e3*np.mean(bed["z99"][m]):.4f} mm', file=out)
    return bed


def analyse(runs, t0, t1, u, wa, area, record=None, rwin=(5.2, 5.5), out=sys.stdout, d=8.7e-4):
    results = {}
    for rd in runs:
        log = os.path.join(rd, 'run_slurm.log')
        win, vm, bed = load_windows(log)
        name = os.path.basename(rd.rstrip('/'))
        mb = mombal_table_quiet(rd, t0, t1)
        series_report(name, mb[0] if mb else None, bed, win, wa, area, out)
        res = table1(name, win, vm, t0, t1, u, d, out)
        mombal_table(rd, t0, t1, wa, area, out)
        m = tmask(bed['t'], t0, t1)
        if m.any():
            print(f'  bed audit window means: dp* {np.mean(bed["dp"][m]*bed["area"][m])/(wa*area):.5f}, floor/W '
                  f'{np.mean(bed["floor"][m])/(wa*area):.5f}, z99 {1e3*np.mean(bed["z99"][m]):.4f} mm '
                  f'(first audit {1e3*bed["z99"][0]:.4f}, last {1e3*bed["z99"][-1]:.4f}), eps_bed '
                  f'{np.mean(bed["eps_bed"][m]):.6f}', file=out)
        drift(rd, d/2, out)
        results[name] = res
    names = list(results)
    if len(names) == 2:
        a, b = names
        print(f'\n== {a} - {b} at a common window:', file=out)
        for w in WINDOWS:
            if w in results[a] and w in results[b]:
                ra, rb = results[a][w], results[b][w]
                print(f'  {w:9s} G/G_cl {ra["ratio"] - rb["ratio"]: .5f}  G {ra["gradz"] - rb["gradz"]: .4e} Pa/m  '
                      f'eps_kern {ra["eps_kern"] - rb["eps_kern"]: .6f}  s_z {ra["slip_z"] - rb["slip_z"]: .4e}  '
                      f'<B s_z> {ra["drag_mean"] - rb["drag_mean"]: .4e} N  F_p {ra["pres_mean"] - rb["pres_mean"]: .4e} N',
                      file=out)
    if record:
        print('', file=out)
        record_means(record, rwin[0], rwin[1], wa, area, out)
    return results


def mombal_table_quiet(rd, t0, t1):
    steps = fl11_mombal.load(os.path.join(rd, 'run_slurm.log'))
    if not steps:
        return None
    s = fl11_mombal.series(steps)
    return s, None


def self_test():
    ok = True
    with tempfile.TemporaryDirectory() as d:
        log = os.path.join(d, 'run_slurm.log')
        lines = []
        for k in range(1, 21):
            t = 0.01*k
            for w, g, e in (('m2d', 3000.0, 0.37), ('m2.5d', 3100.0, 0.365), ('interior', 3200.0, 0.36)):
                lines.append(f'EL_BED_WINDOW t= {t:.8E} step= {k} win= {w} z_lo_hi= 1.0E-03 2.0E-02 gradz= {g + k:.8E} '
                             f'n= 5000 eps_geo= 3.7E-01 eps_kern= {e:.8E} slip= 0.0 0.0 3.0E-01 drag_mean= 2.0E-06 '
                             f'drag_of_means= 1.9E-06 pres_mean= 3.0E-06 visc_mean= 0.0')
            lines.append(f'EL_BED_VERTEX_MODE t= {t:.8E} step= {k} z_lo_hi= 5.0E-03 2.0E-02 vertex= 100 0.3 0.05 0.1 0.5 '
                         f'edge= 300 0.3 0.01 0.25 0.35 face= 300 0.3 0.01 0.25 0.35 centre= 100 0.3 0.01 0.25 0.35')
            lines.append(f'EL_BED_AUDIT t= {t:.6E} step= {k} p_in_out= 1.0E+02 0.0 area_in_out= 2.25E-04 2.25E-04 '
                         f'dp_hydro_dyn= 1.0E+02 0.0 1.0E+02 n= 9789 weight_z_floor_z= -2.844E-02 0.0 '
                         f'z50_95_99_max= 1.0E-02 2.2E-02 2.3E-02 2.4E-02 eps_bed= 3.7E-01 window= 2.0E-03 2.1E-02 '
                         f'gradz_win= 3.1E+03 floor_bal_z= 5.0E-03 eps_kern_slip_win= 3.65E-01 3.0E-01 win_margin= 2.0E-03 '
                         f'gradz_win_dyn= 3.1E+03')
        with open(log, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        win, vm, bed = load_windows(log)
        ok &= len(win['m2.5d']['t']) == 20 and vm['vertex'].shape == (20, 5) and len(bed['t']) == 20
        sink = open(os.devnull, 'w')
        res = table1('x', win, vm, 0.105, 0.2, 0.12, 8.7e-4, sink)
        # audits 0.11 .. 0.20: k = 11..20, gradz mean 3100 + 15.5
        ok &= abs(res['m2.5d']['gradz'] - 3115.5) < 1e-9
        gcl = fl11_closures.closure_grad(0.12, 0.365, 8.7e-4, RHOF, MU)
        ok &= abs(res['m2.5d']['ratio'] - 3115.5/gcl) < 1e-12
        ok &= abs(res['vm_vertex'][2] - 0.05) < 1e-12
        rb = record_means(d, 0.105, 0.2, 126.408, 2.25e-4, sink)
        ok &= len(rb['t']) == 20
        # drift: seed and one log line per particle
        np.savetxt(os.path.join(d, 'particles.xyz'), np.array([[0.001, 0.002, 0.003], [0.0149, 0.002, 0.004]]))
        with open(os.path.join(d, 'el_particles_rank0001.log'), 'w') as f:
            f.write('# time ip fx fy fz tx ty tz px py pz vx vy vz\n')
            f.write('      0.200000       1' + ' 0.0'*6 + ' 1.0E-03 2.0E-03 3.1E-03 0.0 0.0 0.0\n')
            f.write('      0.200000       2' + ' 0.0'*6 + ' 1.0E-04 2.0E-03 4.0E-03 0.0 0.0 0.0\n')
        dd = drift(d, 4.35e-4, sink)
        # particle 2 crossed the periodic x seam: 0.0149 -> 0.0001 is a move of +0.0002
        ok &= dd is not None and abs(dd[0] - 1e-4) < 1e-12 and abs(dd[1] - 2e-4) < 1e-12
    print('fl11_holds self-test', 'PASSED' if ok else 'FAILED')
    return 0 if ok else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('runs', nargs='*')
    ap.add_argument('--window', nargs=2, type=float, default=(0.10, 0.20), metavar=('T0', 'T1'))
    ap.add_argument('--u', type=float, default=0.12, help='superficial inlet velocity of the hold [m/s]')
    ap.add_argument('--wa', type=float, default=fl11_mombal.WA_DEFAULT)
    ap.add_argument('--area', type=float, default=fl11_mombal.AREA_DEFAULT)
    ap.add_argument('--record', help='record run directory for the comparison')
    ap.add_argument('--record-window', nargs=2, type=float, default=(5.2, 5.5))
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args()
    if a.self_test:
        return self_test()
    if not a.runs:
        ap.error('one or two run directories')
    analyse(a.runs, a.window[0], a.window[1], a.u, a.wa, a.area, a.record, a.record_window)
    return 0


if __name__ == '__main__':
    sys.exit(main())
