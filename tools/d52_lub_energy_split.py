#!/usr/bin/env python3
"""D5.2 closure twins: split the lubrication torque delta by energy balance.

The viscometer bob torque T (VISC_TORQUE_DNA, momentum-consistent grad(alpha)
volume form) times Omega is the power the bob feeds into the cell.  In a
statistically steady state that power is dissipated in the fluid and in the
sub-grid lubrication elements of the pe closure:

    Omega * T_L      = Phi_fluid,L + P_lub
    Omega * T_noL    = Phi_fluid,noL
    Omega * (T_L - T_noL) = dPhi_fluid + P_lub

The DNS_LUB record prints, per step, the summed work of the lubrication
impulses on the pairs (`dissipation` = sum J.g + L.w, <= 0, energy per step;
libs/pe/pe/core/lubrication/LubricationStage.h), so P_lub = -dissipation/dt.
The ratio P_lub / (Omega dT) is the DIRECT share of the closure delta (the
film dissipation itself); the rest is INDIRECT (the closure changes particle
kinematics/microstructure and with it the resolved fluid dissipation).

Hard-contact dissipation (inelastic HardContactAndFluid) is not logged and is
folded into the indirect part; negligible at phi <= 0.10 (gaps >= 0.05 d).

Usage:
  d52_lub_energy_split.py RUNDIR --eta-ref 1.1062 [--t0 270] [--T0 84.2296]
                          [--dt 0.005] [--omega 0.1] [--log run_slurm.log]

--T0 must be the empty-instrument torque of the SAME mesh level as RUNDIR
(L3 VISCO2_108: 84.2296 row d52_v20_baseline; L4 VISCO2_431: 83.77503 row
d52_v24f_baseline_hr) and --eta-ref the closure-off twin plateau on that level.
"""
import argparse, os, re, sys
import numpy as np

PAT_LUB = re.compile(r"DNS_LUB time=\s*(\S+)\s+F_total=\s*(\S+)\s+J_max=\s*(\S+)"
                     r"\s+dissipation=\s*(\S+)\s+n_pairs=\s*(\d+)\s+n_saturated=\s*(\d+)")
PAT_DNA = re.compile(r"VISC_TORQUE_DNA time=\s*(\S+)\s+Tz=\s*(\S+)")


def parse(path):
    lub, dna = {}, {}
    with open(path, errors="replace") as f:
        for line in f:
            m = PAT_LUB.search(line)
            if m:
                lub[round(float(m.group(1)), 4)] = tuple(float(x) for x in m.groups()[1:])
                continue
            m = PAT_DNA.search(line)
            if m:
                dna[round(float(m.group(1)), 4)] = float(m.group(2))
    return lub, dna


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundir")
    ap.add_argument("--eta-ref", type=float, required=True, help="closure-off twin plateau eta (same level)")
    ap.add_argument("--t0", type=float, default=270.0, help="window start (t >= t0)")
    ap.add_argument("--t1", type=float, default=1e30, help="window end (t < t1)")
    ap.add_argument("--T0", type=float, default=84.2296, help="empty-instrument torque of this level")
    ap.add_argument("--dt", type=float, default=0.005, help="pe step (json stepsize_ == deck TimeStep)")
    ap.add_argument("--omega", type=float, default=0.1)
    ap.add_argument("--log", default="run_slurm.log")
    a = ap.parse_args()

    lub, dna = parse(os.path.join(a.rundir, a.log))
    if not lub:
        sys.exit("no DNS_LUB records - closure off in this run?")
    tl = np.array(sorted(t for t in lub if a.t0 <= t < a.t1))
    td = np.array(sorted(t for t in dna if a.t0 <= t < a.t1))
    L = np.array([lub[t] for t in tl])           # F_total, J_max, dissipation, n_pairs, n_sat
    T = -np.array([dna[t] for t in td])          # |T_z|

    T_L = T.mean()
    eta_L = T_L / a.T0
    T_ref = a.eta_ref * a.T0
    dT = T_L - T_ref
    P_lub = -L[:, 2] / a.dt
    P = P_lub.mean()
    P_in = a.omega * T_L
    direct = P / (a.omega * dT)

    print(f"rundir            {a.rundir}   window t in [{tl[0]:.3f}, {tl[-1]:.3f}]  "
          f"n_lub={len(tl)} n_torque={len(td)}")
    print(f"T0 (this level)   {a.T0:.5f}   eta_ref (closure off) {a.eta_ref:.4f}  -> T_ref {T_ref:.3f}")
    print(f"closure-on        |T| = {T_L:.3f} (pstd {T.std()/T_L:.1e})  eta_L = {eta_L:.4f}  "
          f"delta = {100*(eta_L/a.eta_ref-1):+.2f} %")
    print(f"power in          Omega*T_L = {P_in:.4f}     Omega*dT = {a.omega*dT:.4f}")
    print(f"closure power     P_lub = {P:.5f} (pstd {P_lub.std()/P:.1e}), per active pair {P/L[:,3].mean():.2e}")
    print(f"                  pairs {L[:,3].mean():.1f} (saturated {L[:,4].mean():.1f}), F_total {L[:,0].mean():.3f}, "
          f"J_max {L[:,1].mean():.2e}")
    print(f"DIRECT film share of the delta   P_lub/(Omega dT) = {direct:.3f}   "
          f"(= {100*P/(a.omega*a.T0)/a.eta_ref:.3f} % of eta_ref)")
    print(f"INDIRECT (fluid dissipation change, incl. unlogged contacts) = {1-direct:.3f}   "
          f"(= {100*(dT-P/a.omega)/a.T0/a.eta_ref:.3f} % of eta_ref)")
    print(f"P_lub as share of total power in  {100*P/P_in:.3f} %")


if __name__ == "__main__":
    main()
