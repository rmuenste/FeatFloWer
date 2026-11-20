#!/usr/bin/env python3
"""Case 11 (minimum fluidisation) evaluation of a q2p1_el_pipeflow run.

Reads run_slurm.log of RUNDIR:
  EL_BED_AUDIT       (SimPar@ELBedAuditFreq > 0; EL_BED_AUDIT module): t, step, plane pressures,
                     dp, its hydrostatic part, buoyant weight, floor contact, bed heights,
                     window porosity
  EL_FORCE_BUDGET    (every step): drag, pressure force, grav_buoy summed over the particles
  EL_INLET_TABLE     (echo of SimPar@ELInletVelocityTable): U(t)
Per velocity step of the protocol (steps file of make_inlet_table.py --steps, or the plateaus
of the echoed table): mean and std of dp_dyn over the last 0.3 s and over the two sub-windows
0.2-0.35 s and 0.35-0.5 s of the step (flagged when they differ by more than --flag of the
plateau), the force-based dp (drag + pressure force)/A, bed height H = z99 + r, porosity of the
bed window.

Gates of the spec (section 5):
  F1  plateau dp A / buoyant weight, mean over the steps with U > 1.2 U_mf       (1 +- 0.02)
  F2  fixed branch: measured dp/dz against the closure gradient at the measured porosity (+-5 %);
      gate value: the window gradient -<dp/dz> of the audit against the closure at eps_bed (both
      over the bed window); also printed: dp_dyn / H_eq (H_eq = n V_p / (A (1 - eps_bed))) and
      dp_dyn / H (H = z99 + r) against the closure at 1 - n V_p / (A H)
  F3  U_mf (intersection of the fixed-branch fit, steps below 0.8 of the plateau, with the
      plateau) against the closure U_mf at the mean fixed-branch porosity               (+-5 %)
  F4  U_mf against the measurement 0.182 m/s (reported)
  F5  particle count constant; highest centre below mid-column at the lowest velocity
The closure is selected with --closure (gidaspow default; ergun, wenyu): dp/dz(U, eps) =
n F_d(U/eps) / eps with the per-particle forms of EL_FORCES (model A; identical to the spec's
umf_reference.py for ergun and wenyu).

usage: fl11_post.py RUNDIR [--steps FILE] [--d 0.87e-3] [--rhop 860] [--rhof 1.0] [--mu 1.8e-5]
                    [--g 9.81] [--lz 0.10] [--closure gidaspow] [--flag 0.02] [--series OUT.csv]
A run that did not reach the protocol end is evaluated over the steps it covers (partial steps
are listed but not used for the gates)."""
import argparse
import math
import os
import re
import sys

NUM = r'[-+]?[0-9.]+(?:[EeDd][-+]?[0-9]+)?'


def f(x):
    return float(x.replace('D', 'E').replace('d', 'e'))


def parse(rundir):
    audit, budget, table, pairs = [], {}, [], {}
    pat_a = re.compile(r'^EL_BED_AUDIT t=\s*(%s) step=\s*(\d+) p_in_out=\s*(%s)\s+(%s) area_in_out=\s*(%s)\s+(%s)'
                       r' dp_hydro_dyn=\s*(%s)\s+(%s)\s+(%s) n=\s*(\d+) weight_z_floor_z=\s*(%s)\s+(%s)'
                       r' z50_95_99_max=\s*(%s)\s+(%s)\s+(%s)\s+(%s) eps_bed=\s*(%s) window=\s*(%s)\s+(%s)'
                       r'(?: gradz_win=\s*(%s))?(?: floor_bal_z=\s*(%s))?(?: eps_kern_slip_win=\s*(%s)\s+(%s))?'
                       % ((NUM,) + (NUM,) * 7 + (NUM,) * 2 + (NUM,) * 4 + (NUM,) * 3 + (NUM,) * 4))
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
                              gw=f(g[19]) if g[19] else float('nan'),
                              fbal=f(g[20]) if g[20] else None,
                              ekern=f(g[21]) if g[21] else float('nan'),
                              skern=f(g[22]) if g[22] else float('nan')))
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
                row['fbal'] = (pr[0] - pr[1])/dt if pr else float('nan')
    return audit, budget, table


def u_of_t(table, t):
    if t <= table[0][0]:
        return table[0][1]
    for (t0, u0), (t1, u1) in zip(table, table[1:]):
        if t <= t1:
            return u0 + (u1 - u0) * (t - t0) / (t1 - t0)
    return table[-1][1]


def plateaus(table):
    """steps from the echoed table: intervals of constant U (the preparation is merged with the
    first plateau; the evaluation windows are taken from its end)"""
    out = []
    for (t0, u0), (t1, u1) in zip(table, table[1:]):
        if u0 == u1 and t1 - t0 > 0.1:
            out.append((t0, t1, u0))
    return out


def closure_grad(u, eps, d, rhof, mu, closure):
    """|dp/dz| of a fixed bed: n F_d(U/eps)/eps with drag_B of EL_FORCES (model-A form)"""
    vp = math.pi * d ** 3 / 6
    ur = u / eps

    def ergun():
        return vp * (150 * (1 - eps) * mu / (eps * d * d) + 1.75 * rhof * ur / d)

    def wenyu():
        re_ = rhof * eps * ur * d / mu
        cds = 24 * mu / (rhof * eps * d) * (1 + 0.15 * re_ ** 0.687) if re_ < 1000 else 0.44 * ur
        return vp * 0.75 * cds * rhof * eps * eps ** -2.65 / d
    if closure == 'ergun':
        b = ergun()
    elif closure == 'wenyu':
        b = wenyu()
    else:
        b = ergun() if eps <= 0.8 else wenyu()
    return (1 - eps) / vp * b * ur / eps


def closure_umf(eps, d, rhop, rhof, mu, g, closure):
    target = (1 - eps) * (rhop - rhof) * g
    a, b = 1e-6, 10.0
    for _ in range(200):
        m = 0.5 * (a + b)
        if closure_grad(m, eps, d, rhof, mu, closure) < target:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def stats(rows, key, t0, t1):
    v = [r[key] for r in rows if t0 <= r['t'] <= t1]
    if not v:
        return float('nan'), float('nan'), 0
    m = sum(v) / len(v)
    s = math.sqrt(sum((x - m) ** 2 for x in v) / len(v))
    return m, s, len(v)


def linfit(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if n < 2 or sxx == 0:
        return float('nan'), float('nan')
    s = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx
    return s, my - s * mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('rundir')
    ap.add_argument('--steps', default=None)
    ap.add_argument('--d', type=float, default=0.87e-3)
    ap.add_argument('--rhop', type=float, default=860.0)
    ap.add_argument('--rhof', type=float, default=1.0)
    ap.add_argument('--mu', type=float, default=1.8e-5)
    ap.add_argument('--g', type=float, default=9.81)
    ap.add_argument('--lz', type=float, default=0.10)
    ap.add_argument('--closure', default='gidaspow', choices=('gidaspow', 'ergun', 'wenyu'))
    ap.add_argument('--flag', type=float, default=0.02)
    ap.add_argument('--series', default=None)
    a = ap.parse_args()

    audit, budget, table = parse(a.rundir)
    if not audit:
        sys.exit('no EL_BED_AUDIT lines (SimPar@ELBedAuditFreq not set?)')
    if not table:
        sys.exit('no EL_INLET_TABLE echo (SimPar@ELInletVelocityTable not set?)')
    r = a.d / 2
    vp = math.pi * a.d ** 3 / 6
    area = audit[-1]['a_in']
    for row in audit:
        row['U'] = u_of_t(table, row['t'])
        fb = budget.get(row['step'])
        row['dp_force'] = (fb[0] + fb[1]) / area if fb else float('nan')
        row['dp_drag'] = fb[0] / area if fb else float('nan')
        row['H'] = row['z99'] + r
    W = -audit[-1]['wz']
    wa = W / area
    n0 = audit[0]['n']
    tend = audit[-1]['t']
    print(f"{os.path.basename(os.path.abspath(a.rundir))}: {len(audit)} audit lines, t {audit[0]['t']:.4f} .. {tend:.4f} s, "
          f"N {n0}, A {area:.6e} m2, buoyant weight {W:.6e} N, W/A {wa:.3f} Pa, "
          f"hydrostatic part in dp {audit[-1]['dp_h']:.4f} Pa (subtracted: dp_dyn)")
    if a.series:
        with open(a.series, 'w') as fh:
            fh.write('t,step,U,dp,dp_hydro,dp_dyn,dp_force,dp_drag,floor_z,z50,z95,z99,zmax,eps_bed,n\n')
            for row in audit:
                fh.write(','.join(str(row[k]) for k in ('t', 'step', 'U', 'dp', 'dp_h', 'dp_dyn', 'dp_force',
                                                        'dp_drag', 'floor', 'z50', 'z95', 'z99', 'zmax', 'eps', 'n')) + '\n')

    if a.steps:
        steps = [tuple(map(float, l.split()[:3])) for l in open(a.steps) if l.strip() and not l.startswith('#')]
    else:
        steps = plateaus(table)
        # the first plateau of the protocol holds preparation + step 0: evaluate its last 0.5 s
        if steps and steps[0][1] - steps[0][0] > 0.6:
            steps[0] = (steps[0][1] - 0.5, steps[0][1], steps[0][2])
    print('\n step  U [m/s]  window      dp*      std     sub1    sub2   flag  F/(W/A)  H [mm]  eps_bed  floor/W  dp*+floor  (floor_pe/W)')
    res = []
    for k, (t0, t1, u) in enumerate(steps):
        complete = tend >= t1 - 1e-9
        m, s, cnt = stats(audit, 'dp_dyn', t1 - 0.3, t1)
        m1, _, _ = stats(audit, 'dp_dyn', t0 + 0.2, t0 + 0.35)
        m2, _, _ = stats(audit, 'dp_dyn', t0 + 0.35, t0 + 0.5)
        mf, _, _ = stats(audit, 'dp_force', t1 - 0.3, t1)
        h, _, _ = stats(audit, 'H', t1 - 0.3, t1)
        e, _, _ = stats(audit, 'eps', t1 - 0.3, t1)
        fl, _, _ = stats(audit, 'fbal', t1 - 0.3, t1)
        flpe, _, _ = stats(audit, 'floor', t1 - 0.3, t1)
        gw, _, _ = stats(audit, 'gw', t1 - 0.3, t1)
        ek, _, _ = stats(audit, 'ekern', t1 - 0.3, t1)
        flag = abs(m1 - m2) > a.flag * wa
        res.append(dict(k=k, U=u, dp=m, std=s, cnt=cnt, eps=e, H=h, gw=gw, ek=ek, fl=fl/W,
                        complete=complete and cnt > 0))
        print(f" {k:3d}  {u:6.3f}  {t1-0.3:5.2f}-{t1:5.2f}  {m/wa:7.4f}  {s/wa:7.4f}  {m1/wa:6.3f}  {m2/wa:6.3f}  "
              f"{'DIFF' if flag else '    '}  {mf/wa:7.4f}  {h*1e3:6.2f}  {e:7.4f}  {fl/W:7.4f}  {m/wa + fl/W:8.4f}   ({flpe/W:.4f})"
              f"{'' if complete else '  (incomplete)'}")
    print('  dp* = dp_dyn A / W over 0.2-0.5 s of the step; floor from the particle momentum balance; floor_pe = the pe '
          'wall-contact impulse (unrelaxed PGS impulse, 1/0.9 of the applied one)')
    use = [x for x in res if x['complete']]
    if len(use) < 3:
        print('\nfewer than three complete steps: no gates')
        return
    # ---- F1 (spec rule of 2026-11-20): mean dp* over all steps with dp* > 0.95 ------------------
    flu = [x for x in use if x['dp'] > 0.95 * wa]
    if flu:
        plateau = sum(x['dp'] for x in flu) / len(flu)
        print(f"\nF1 plateau: mean dp* over {len(flu)} steps with dp* > 0.95 (U {min(x['U'] for x in flu):.2f} .. "
              f"{max(x['U'] for x in flu):.2f}): {plateau / wa:.4f}   (gate 1 +- 0.02: {'PASS' if abs(plateau / wa - 1) <= 0.02 else 'FAIL'}); "
              f"dp* + floor/W over these {sum(x['dp']/wa + x['fl'] for x in flu)/len(flu):.4f}")
    else:
        plateau = wa
        print('\nF1 not evaluated: no step with dp* > 0.95; plateau taken as W/A')
    fixed = [x for x in use if x['dp'] < 0.8 * plateau]
    if len(fixed) < 2:
        print('fewer than two fixed-branch steps (dp* < 0.8 of the plateau): no F2 / F3')
        return
    efix = sum(x['eps'] for x in fixed) / len(fixed)
    uc = closure_umf(efix, a.d, a.rhop, a.rhof, a.mu, a.g, a.closure)
    # ---- F3 (agreed with the peer): Ergun form dp* = a U + b U^2 through the origin over the fixed steps
    su2 = sum(x['U']**2 for x in fixed); su3 = sum(x['U']**3 for x in fixed); su4 = sum(x['U']**4 for x in fixed)
    sy1 = sum(x['U']*x['dp']/wa for x in fixed); sy2 = sum(x['U']**2*x['dp']/wa for x in fixed)
    det = su2*su4 - su3*su3
    ca = (sy1*su4 - sy2*su3)/det; cb = (su2*sy2 - su3*sy1)/det
    def cross(level):
        # positive root of b U^2 + a U = level; nan when the fit has no positive crossing
        if abs(cb) < 1e-14: return level/ca if ca > 0 else float('nan')
        disc = ca*ca + 4*cb*level
        if disc < 0: return float('nan')
        roots = [x for x in ((-ca + math.sqrt(disc))/(2*cb), (-ca - math.sqrt(disc))/(2*cb)) if x > 0]
        return min(roots) if roots else float('nan')
    u1 = cross(1.0); up = cross(plateau/wa)
    print(f"F3 fixed branch {len(fixed)} steps (U {min(x['U'] for x in fixed):.3f} .. {max(x['U'] for x in fixed):.3f}): "
          f"dp* = {ca:.4f} U + {cb:.4f} U^2 (through the origin)")
    print(f"   U_mf at dp* = 1: {u1:.4f} m/s, at the measured plateau {plateau/wa:.4f}: {up:.4f} m/s; {a.closure} U_mf at the "
          f"mean fixed-branch core porosity {efix:.4f}: {uc:.4f} m/s; ratios {u1/uc:.4f} / {up/uc:.4f}   (gate 0.05: "
          f"{'PASS' if abs(u1/uc-1) <= 0.05 else 'FAIL'} / {'PASS' if abs(up/uc-1) <= 0.05 else 'FAIL'})")
    print(f"F4 U_mf / 0.182 (NETL measurement): {u1/0.182:.4f} (dp* = 1) / {up/0.182:.4f} (plateau) (reported, not gated)")
    # ---- secondary (rules before 2026-11-20): F1 over U > 1.2 U_mf of the linear fit, linear-fit U_mf
    pl2 = wa
    umf2 = float('nan')
    for _ in range(3):
        fx2 = [x for x in use if x['dp'] < 0.8 * pl2]
        if len(fx2) >= 2:
            s2, i2 = linfit([x['U'] for x in fx2], [x['dp'] for x in fx2])
            umf2 = (pl2 - i2) / s2 if s2 else float('nan')
        fl2 = [x for x in use if (not math.isnan(umf2)) and x['U'] > 1.2 * umf2] or \
              [x for x in use if x['dp'] >= 0.8 * pl2]
        if fl2:
            pl2 = sum(x['dp'] for x in fl2) / len(fl2)
    s, i = linfit([x['U'] for x in fixed], [x['dp'] for x in fixed])
    print(f"   secondary: F1 over the {len(fl2)} steps with U > 1.2 U_mf(linear) {pl2/wa:.4f}; linear fit slope {s:.1f} Pa s/m "
          f"(intercept {i:.2f} Pa), U_mf(linear) {umf2:.4f} m/s = {umf2/uc:.4f} of the closure")
    # F2 in three readings: (w) the window gradient -<dp/dz> of the audit (bed core, the definition
    # that matches eps_bed; gate value when present), (e) dp_dyn / H_eq with H_eq the height of a
    # uniform bed at eps_bed, (h) dp_dyn / H with H = z99 + r against the closure at the mean
    # porosity of that height (end layers included)
    print(' F2: step  U      eps_bed  window-grad  closure(eps_bed)  ratio_w | eps_kern ratio_k | dp/H_eq  ratio_e | eps_H   dp/H   closure(eps_H)  ratio_h')
    worst = 0.0
    use_w = all(not math.isnan(x['gw']) for x in fixed)
    for x in fixed:
        heq = n0 * vp / (area * (1 - x['eps']))
        clo = closure_grad(x['U'], x['eps'], a.d, a.rhof, a.mu, a.closure)
        eh = 1 - n0 * vp / (area * x['H'])
        cloh = closure_grad(x['U'], eh, a.d, a.rhof, a.mu, a.closure)
        rw = x['gw'] / clo
        rk = x['gw'] / closure_grad(x['U'], x['ek'], a.d, a.rhof, a.mu, a.closure) if x['ek'] == x['ek'] and x['ek'] > 0 else float('nan')
        re_ = x['dp'] / heq / clo
        rh = x['dp'] / x['H'] / cloh
        worst = max(worst, abs((rw if use_w else re_) - 1))
        print(f"     {x['k']:3d}  {x['U']:.3f}  {x['eps']:.4f}  {x['gw']:10.1f}  {clo:10.1f}  {rw:12.4f} | {x['ek']:.4f} {rk:7.4f} | {x['dp']/heq:8.1f} {re_:7.4f} | "
              f"{eh:.4f} {x['dp']/x['H']:8.1f} {cloh:10.1f}  {rh:9.4f}")
    print(f"F2 worst |ratio - 1| ({'window gradient' if use_w else 'dp/H_eq'}) {worst:.4f}   (gate 0.05: {'PASS' if worst <= 0.05 else 'FAIL'})")
    print('   ratio_k: the same gradient against the closure at the mean kernel-sampled porosity of the window '
          'particles (the porosity the drag uses; like-for-like reference of an unresolved method)')
    nset = {row['n'] for row in audit}
    lowest = min(use, key=lambda x: x['U'])
    zlast = [row['zmax'] for row in audit if row['t'] >= tend - 1e-9][0]
    print(f"F5 particle counts seen {sorted(nset)}; highest centre at the end {zlast*1e3:.2f} mm "
          f"(mid-column {a.lz*500:.1f} mm; lowest velocity evaluated {lowest['U']:.3f} m/s): "
          f"{'PASS' if len(nset) == 1 and zlast < a.lz / 2 else 'FAIL'}")


if __name__ == '__main__':
    main()
