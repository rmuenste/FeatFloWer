#!/usr/bin/env python3
"""D6.1 Oberbeck spheroid drag - gate analysis (CASE_SPEC d61_oberbeck section 6).

Per run: steady wrench from particle_force.log (plateau check + torque null),
steady superficial velocity from bulk_flow.log, hydrodynamic radius R_h from the
self-consistent Hasimoto inversion (d11 convention: F = f*V_cell whole-cell
balance, U = superficial velocity), then the gates:

  G-ratio  (PRIMARY): R_h(V2)/R_h(V1) vs Y^A/X^A            (+-2%)
  G-abs             : R_h vs a*X^A / a*Y^A, raw and a_eff-corrected
                      (axis-wise a+0.14h, d11_rh_collapse)   (+-3% corrected)
  G-torque          : |T| / (a_major * |F|) < 1e-3 every run
  G-offdiag (V3)    : angle(F, axis) vs atan(Y^A/X^A * tan 45) = 48.88 deg;
                      drift of F from U = 3.88 deg           (+-0.5 deg)
  V0 anchor         : raw K = F/(6 pi mu r U_sup) vs the certified 1.7404 (L3)

Usage: d61_oberbeck_analysis.py <rundir_v0> <rundir_v1> <rundir_v2> [rundir_v3]
Each rundir needs particle_force.log and bulk_flow.log.
"""
import math
import sys

MU = 1.0           # Prop@Viscosity (campaign units, rho=1)
L_CELL = 1.0       # periodic cell edge
R_SPHERE = 1.0 / 6.0
RE_ASPECT = 2.0
B_MINOR = R_SPHERE / 2.0 ** (1.0 / 3.0)      # 0.13228342...
A_MAJOR = 2.0 * B_MINOR                      # 0.26456684...
H_L3 = 1.0 / 36.0
AEFF_C = 0.14                                # a_eff ~ a + 0.14 h (d11_rh_collapse)
K_L3_CERTIFIED = 1.7404                      # RUNBOOK L3 row (job 137540)
SOLID_NOMINAL = 0.019392547                  # (4/3) pi r^3 = (4/3) pi a b^2 (matched)
DLNK_DLNR = 1.87                             # d11 sensitivity dlnK/dlnr


def resistance_functions(re_aspect):
    e = math.sqrt(1.0 - 1.0 / re_aspect ** 2)
    L = math.log((1.0 + e) / (1.0 - e))
    XA = (8.0 / 3.0) * e ** 3 / (-2.0 * e + (1.0 + e * e) * L)
    YA = (16.0 / 3.0) * e ** 3 / (2.0 * e + (3.0 * e * e - 1.0) * L)
    return XA, YA


def hasimoto_K(phi):
    return 1.0 / (1.0 - 1.7601 * phi ** (1.0 / 3.0) + phi - 1.5593 * phi * phi)


def invert_rh(F, U):
    """Fixed point: F/(6 pi mu R U) = K_H(phi(R)), phi = (4/3) pi R^3 / L^3."""
    R = R_SPHERE
    for _ in range(200):
        phi = (4.0 / 3.0) * math.pi * R ** 3 / L_CELL ** 3
        R_new = F / (6.0 * math.pi * MU * U * hasimoto_K(phi))
        if abs(R_new - R) < 1e-14:
            R = R_new
            break
        R = 0.5 * (R + R_new)
    return R


def read_run(rundir):
    rows = []
    with open(f"{rundir}/particle_force.log") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append([float(x) for x in (v[0], *v[2:8])])
    t, fx, fy, fz, tx, ty, tz = rows[-1]
    prev_fz = rows[-2][3] if len(rows) > 1 else fz
    bulk = []
    with open(f"{rundir}/bulk_flow.log") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            v = line.split()
            bulk.append((float(v[0]), float(v[1]), float(v[3])))
    tb, usup, ffrac = bulk[-1]
    return dict(t=t, F=(fx, fy, fz), T=(tx, ty, tz),
                plateau_dF=abs(fz - prev_fz) / max(abs(fz), 1e-30),
                U=usup, fluid_frac=ffrac)


def main(argv):
    names = ["V0", "V1", "V2", "V3"][: len(argv)]
    runs = {n: read_run(d) for n, d in zip(names, argv)}
    XA, YA = resistance_functions(RE_ASPECT)
    print(f"X^A={XA:.6f}  Y^A={YA:.6f}  Y/X={YA/XA:.5f}")
    print(f"{'run':4} {'t_end':>6} {'|F|':>12} {'U_sup':>12} {'R_h':>10} "
          f"{'|T|/(a|F|)':>11} {'plateau dF/F':>13}")
    rh = {}
    rhc = {}
    for n in names:
        r = runs[n]
        Fmag = math.sqrt(sum(f * f for f in r["F"]))
        Tmag = math.sqrt(sum(x * x for x in r["T"]))
        rh[n] = invert_rh(Fmag, r["U"])
        # d11_rh_collapse volume correction: dr from the run's own fluid_frac,
        # K divided by (1+dr)^1.87 == invert with F/(1+dr)^1.87
        dr = ((1.0 - r["fluid_frac"]) / SOLID_NOMINAL) ** (1.0 / 3.0) - 1.0
        rhc[n] = invert_rh(Fmag / (1.0 + dr) ** DLNK_DLNR, r["U"])
        r["dr"] = dr
        print(f"{n:4} {r['t']:6.2f} {Fmag:12.6e} {r['U']:12.6e} {rh[n]:10.6f} "
              f"{Tmag/(A_MAJOR*Fmag):11.2e} {r['plateau_dF']:13.2e}  dr={dr:+.4%} "
              f"R_h_corr={rhc[n]:.6f}")

    print("\n--- gates ---")
    if "V0" in runs:
        K0 = (math.sqrt(sum(f * f for f in runs['V0']['F']))
              / (6.0 * math.pi * MU * R_SPHERE * runs['V0']['U']))
        print(f"V0 anchor: raw K = {K0:.4f} vs certified {K_L3_CERTIFIED} "
              f"({(K0/K_L3_CERTIFIED-1)*100:+.3f}%)")
    if "V1" in runs and "V2" in runs:
        ratio = rh["V2"] / rh["V1"]
        ratio_c = rhc["V2"] / rhc["V1"]
        tgt = YA / XA
        print(f"G-ratio : R_h(perp)/R_h(par) = {ratio:.5f} raw / {ratio_c:.5f} "
              f"vol-corr vs {tgt:.5f} ({(ratio/tgt-1)*100:+.2f}% / "
              f"{(ratio_c/tgt-1)*100:+.2f}%)  [band +-2%]")
        # absolutes, raw and a_eff-corrected
        for n, coef in (("V1", XA), ("V2", YA)):
            tgt_abs = A_MAJOR * coef
            a_e, b_e = A_MAJOR + AEFF_C * H_L3, B_MINOR + AEFF_C * H_L3
            XAe, YAe = resistance_functions(a_e / b_e)
            coef_e = XAe if n == "V1" else YAe
            tgt_eff = a_e * coef_e
            print(f"G-abs {n}: R_h raw {rh[n]:.6f} / vol-corr {rhc[n]:.6f} vs "
                  f"a*C = {tgt_abs:.6f} ({(rhc[n]/tgt_abs-1)*100:+.2f}% corr-vs-nominal) "
                  f"| a_eff target {tgt_eff:.6f} ({(rhc[n]/tgt_eff-1)*100:+.2f}%)  "
                  f"[band +-3%]")
    if "V3" in runs:
        fx, fy, fz = runs["V3"]["F"]
        axis = (1.0 / math.sqrt(2.0), 0.0, 1.0 / math.sqrt(2.0))
        Fmag = math.sqrt(fx * fx + fy * fy + fz * fz)
        cos_fa = (fx * axis[0] + fy * axis[1] + fz * axis[2]) / Fmag
        ang_axis = math.degrees(math.acos(max(-1.0, min(1.0, cos_fa))))
        ang_pred = math.degrees(math.atan(YA / XA))
        drift = ang_axis - 45.0
        print(f"G-offdiag [F-form, INVALID for body-force driving]: "
              f"angle(F,axis) = {ang_axis:.2f} deg (steady-state momentum "
              f"balance forces F -> f*V_cell, transverse F -> 0; the "
              f"off-diagonal mobility appears as a TRANSVERSE MEAN FLOW "
              f"U_x/U_z = (Y-X)/(Y+X) = {(YA-XA)/(YA+XA):+.4f} = "
              f"{math.degrees(math.atan((YA-XA)/(YA+XA))):+.2f} deg instead - "
              f"needs the 3-component bulk-flow diagnostic, gate redesigned)")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
