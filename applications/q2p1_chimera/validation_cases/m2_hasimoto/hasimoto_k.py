#!/usr/bin/env python3
"""Hasimoto drag ratio K from a q2p1_chimera periodic sphere-array run.

Reads the `ChimeraForce<k>:` and `ChimeraBulk:` lines (protocol file or
run log) and reports, for the LAST coupling update (or every --every
steps for the plateau history):

  F_s     = sum_k F_z(k)                 surface-traction force (Chimera)
  F_H     = F_s + f_z V_solid            Hasimoto convention: the sphere force
                                          balances the MEAN pressure gradient
                                          over the whole cell, i.e. it contains
                                          the mean-gradient "buoyancy" f V_solid
                                          that a body force acting on the fluid
                                          only does not exert on the surface
                                          (the FBM constraint force contains it
                                          automatically: the hole fluid is forced)
  U_sup   = <u>_z of the whole cell      (ChimeraBulk column 4, k = 0 mode)
  K_meas  = F_H / (6 pi mu R U_sup)      per sphere
  K_bal   = f_z V_cell / (6 pi mu R U_sup)   force from the exact momentum
                                          balance (F_s -> f_z V_fluid at steady
                                          state); K_meas - K_bal measures the
                                          discrete momentum leak of the coupling
  K_ref   = 1 / (1 - 1.7601 phi^(1/3) + phi - 1.5593 phi^2)   Hasimoto (1959)
  balance = F_s / (f_z V_fluid) - 1      (steady, conservative scheme: 0)

Conventions pinned in applications/q2p1_dns_drag/validation_cases/
d11_hasimoto/RUNBOOK.md: U is the k = 0 Fourier mode = cell volume average
(superficial velocity), F balances the mean pressure gradient over the
WHOLE cell (F = f V_cell for body-force driving of the whole cell).

Usage:
  hasimoto_k.py _data/prot.txt --phi 0.0193925 --mu 1 --fz 1e-2 [--every 50]
  (R defaults to the radius of phi in a unit cell with one sphere; pass --R
   and --nbody for lattices with several spheres per cell)
"""

import argparse
import math
import re
import sys


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("logfile")
    ap.add_argument("--phi", type=float, required=True, help="solid volume fraction")
    ap.add_argument("--mu", type=float, default=1.0, help="dynamic viscosity")
    ap.add_argument("--fz", type=float, default=1e-2, help="body force density f_z")
    ap.add_argument("--L", type=float, default=1.0, help="cell edge length")
    ap.add_argument("--nbody", type=int, default=1, help="spheres per cell")
    ap.add_argument("--R", type=float, default=None, help="sphere radius")
    ap.add_argument("--every", type=int, default=0, help="print every n-th update")
    a = ap.parse_args(argv)

    vcell = a.L ** 3
    R = a.R if a.R is not None else (a.phi * 3.0 * vcell / (4.0 * math.pi * a.nbody)) ** (1.0 / 3.0)
    phi = a.phi
    kref = 1.0 / (1.0 - 1.7601 * phi ** (1.0 / 3.0) + phi - 1.5593 * phi ** 2)

    force_re = re.compile(r"^ChimeraForce(\d+):\s+(.*)$")
    bulk_re = re.compile(r"^ChimeraBulk:\s+(.*)$")
    rows = []          # (time, Fz_sum, Uz, fluidfrac)
    fz_acc = {}
    with open(a.logfile) as f:
        for line in f:
            m = force_re.match(line)
            if m:
                vals = [float(x.replace("D", "E")) for x in m.group(2).split()]
                t = vals[0]
                fz_acc.setdefault(t, 0.0)
                fz_acc[t] += vals[5]
                continue
            m = bulk_re.match(line)
            if m:
                vals = [float(x.replace("D", "E")) for x in m.group(1).split()]
                t = vals[0]
                rows.append((t, fz_acc.get(t, vals[7]), vals[3], vals[4]))
    if not rows:
        sys.exit("no ChimeraBulk lines in %s" % a.logfile)

    vsolid = phi * vcell
    vfluid = vcell - vsolid
    print("# Hasimoto: phi=%.7g R=%.7g mu=%g fz=%g nbody=%d  K_ref=%.4f"
          % (phi, R, a.mu, a.fz, a.nbody, kref))
    print("# %10s %13s %13s %8s %8s %8s %8s %8s %9s" % (
        "time", "F_s(sum Fz)", "U_sup", "K_meas", "dK_meas", "K_bal", "dK_bal",
        "balance", "ffrac"))
    if a.every <= 0:
        sel = [rows[-1]]
    else:
        sel = rows[::a.every]
        if sel[-1] is not rows[-1]:
            sel.append(rows[-1])
    denom0 = 6.0 * math.pi * a.mu * R * a.nbody
    for t, fs, uz, ff in sel:
        if uz == 0.0:
            continue
        fh = fs + a.fz * vsolid
        kmeas = fh / (denom0 * uz)
        kbal = a.fz * vcell / (denom0 * uz)
        bal = fs / (a.fz * vfluid) - 1.0
        print("  %10.4f %13.6e %13.6e %8.4f %+8.4f %8.4f %+8.4f %+8.4f %9.6f"
              % (t, fs, uz, kmeas, kmeas / kref - 1.0, kbal, kbal / kref - 1.0, bal, ff))


if __name__ == "__main__":
    main(sys.argv[1:])
