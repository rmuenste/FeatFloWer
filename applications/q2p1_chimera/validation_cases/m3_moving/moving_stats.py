#!/usr/bin/env python3
"""Force/kinematics statistics of a q2p1_chimera moving-body run (Phase 6).

Reads `ChimeraForce<k>:`, `ChimeraBulk:` and `ChimeraBody<k>:` lines and
prints, over the time window [--tmin, tmax]:

  per body:  <F>, min/max, peak-to-peak/|<F>| of the force component --comp
             (force continuity across background-cell crossings), the
             final position/velocity, and the distance travelled in cells
             of size --h;
  cell:      <u>_comp of the composite bulk velocity (lab frame) and the
             relative superficial velocity U_rel = <u>_cell - U_body
             (single body), plus K = F_ref/(6 pi mu R U_rel) with F_ref =
             --fref (e.g. the external force of a sedimenting sphere) or,
             without --fref, the measured mean force.

Usage:
  moving_stats.py _data/prot.txt --tmin 2.0 --h 0.08333 --R 0.1666665 --mu 1
  moving_stats.py prot_SED_S.txt --tmin 2.0 --R 0.1666665 --fref 1e-2
"""

import argparse
import math
import re
import sys


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("logfile")
    ap.add_argument("--tmin", type=float, default=0.0)
    ap.add_argument("--comp", type=int, default=3, help="component 1..3 (default z)")
    ap.add_argument("--h", type=float, default=None, help="background cell size")
    ap.add_argument("--R", type=float, default=None)
    ap.add_argument("--mu", type=float, default=1.0)
    ap.add_argument("--fref", type=float, default=None)
    a = ap.parse_args(argv)
    c = a.comp - 1

    force_re = re.compile(r"^ChimeraForce(\d+):\s+(.*)$")
    body_re = re.compile(r"^ChimeraBody(\d+):\s+(.*)$")
    bulk_re = re.compile(r"^ChimeraBulk:\s+(.*)$")
    forces, bodies, bulk = {}, {}, []
    with open(a.logfile) as f:
        for line in f:
            m = force_re.match(line)
            if m:
                v = [float(x.replace("D", "E")) for x in m.group(2).split()]
                forces.setdefault(int(m.group(1)), []).append((v[0], v[3:6]))
                continue
            m = body_re.match(line)
            if m:
                v = [float(x.replace("D", "E")) for x in m.group(2).split()]
                bodies.setdefault(int(m.group(1)), []).append((v[0], v[1:4], v[4:7], v[7:10]))
                continue
            m = bulk_re.match(line)
            if m:
                v = [float(x.replace("D", "E")) for x in m.group(1).split()]
                bulk.append((v[0], v[1:4]))
    if not forces:
        sys.exit("no ChimeraForce lines")

    ubody = None
    for k in sorted(forces):
        sel = [fv for t, fv in forces[k] if t >= a.tmin]
        if not sel:
            continue
        vals = [fv[c] for fv in sel]
        mean = sum(vals) / len(vals)
        p2p = (max(vals) - min(vals)) / abs(mean) if mean != 0 else float("nan")
        line = "body %d: <F_%d> = %.7e  min %.7e  max %.7e  p2p/|mean| = %.2e  (n = %d, t >= %g)" % (
            k, a.comp, mean, min(vals), max(vals), p2p, len(vals), a.tmin)
        if k in bodies:
            t, X, U, W = bodies[k][-1]
            X0 = bodies[k][0][1]
            dist = math.sqrt(sum((X[i] - X0[i]) ** 2 for i in range(3)))
            line += "\n         final t = %.4f  X = (%.6f %.6f %.6f)  U = (%.6e %.6e %.6e)  |W| = %.3e" % (
                t, X[0], X[1], X[2], U[0], U[1], U[2], math.sqrt(sum(w * w for w in W)))
            if a.h:
                line += "  travelled %.2f cells" % (dist / a.h)
            ubody = U[c]
        print(line)
        if len(forces) == 1:
            fmean = mean

    selb = [u for t, u in bulk if t >= a.tmin]
    if selb:
        ub = sum(u[c] for u in selb) / len(selb)
        print("cell: <u_%d> (composite, lab frame) = %.7e over %d steps" % (a.comp, ub, len(selb)))
        if ubody is not None and a.R:
            urel = ub - ubody
            fref = a.fref if a.fref is not None else abs(fmean)
            K = fref / (6.0 * math.pi * a.mu * a.R * abs(urel)) if urel != 0 else float("nan")
            print("      U_body = %.7e  U_rel = <u>_cell - U_body = %.7e  K = F_ref/(6 pi mu R |U_rel|) = %.4f"
                  % (ubody, urel, K))


if __name__ == "__main__":
    main(sys.argv[1:])
