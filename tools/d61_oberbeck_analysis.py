#!/usr/bin/env python3
"""D6.1 Oberbeck spheroid drag - gate analysis (CASE_SPEC d61_oberbeck section 6).

Per run: steady wrench from particle_force.log (trailing-window mean, plateau
check + torque null), steady superficial velocity from bulk_flow.log (same
window), hydrodynamic radius R_h from the self-consistent Hasimoto inversion
(d11 convention: F = f*V_cell whole-cell balance, U = superficial velocity),
then the gates:

  G-ratio  (PRIMARY): R_h(V2)/R_h(V1) vs Y^A/X^A            (+-2%)
  G-abs             : R_h vs a*X^A / a*Y^A, raw and a_eff-corrected
                      (axis-wise a - 0.14h, d11_aeff_sign_erratum) (+-3% corrected)
  G-torque          : |T| / (a_major * |F|) < 1e-3 every run
  G-offdiag (V3)    : F-form is structurally null under body-force driving
                      (momentum balance) - reported for the record only
  V0 anchor         : raw K = F/(6 pi mu r U_sup) vs the certified 1.7404 (L3)

Usage: d61_oberbeck_analysis.py [--h H] [--scale S] [--window W] [--sensitivity]
                                <rundir_v0> <rundir_v1> <rundir_v2> [rundir_v3]
  --h            mesh spacing of the run level (1/36 at L3, 1/72 at L4);
                 REQUIRED knowledge for the a_eff correction - default L3
  --scale        body-size scale (V5 half-size rung: 0.5); env D61_SCALE also works
  --window       trailing time window (t.u.) averaged for F, T, U; 0 = last sample
                 (the pre-2026-09-10 behaviour). Default 0.5
  --sensitivity  print R_h and the ratio for windows 0 / 0.25 / 0.5 / 1.0 / 2.0
                 so the plateau claim is checked, not assumed (review finding 3)
Each rundir needs particle_force.log and bulk_flow.log.
"""
import argparse
import math
import os
import sys

MU = 1.0           # Prop@Viscosity (campaign units, rho=1)
L_CELL = 1.0       # periodic cell edge
R_SPHERE = 1.0 / 6.0
RE_ASPECT = 2.0
AEFF_C = -0.14                               # a_eff = a - 0.14 h (d11_aeff_sign_erratum)
K_L3_CERTIFIED = 1.7404                      # RUNBOOK L3 row (job 137540)
DLNK_DLNR = 1.87                             # d11 sensitivity dlnK/dlnr
SENS_WINDOWS = (0.0, 0.25, 0.5, 1.0, 2.0)


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


def load_run(rundir):
    rows = []
    with open(f"{rundir}/particle_force.log") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append([float(x) for x in (v[0], *v[2:8])])
    bulk = []
    with open(f"{rundir}/bulk_flow.log") as fh:
        for line in fh:
            if line.lstrip().startswith("#") or not line.strip():
                continue
            v = line.split()
            bulk.append((float(v[0]), float(v[1]), float(v[3])))
    return rows, bulk


def window_mean(samples, window, col):
    """Mean of column col over the trailing window (time in column 0); 0 = last."""
    t_end = samples[-1][0]
    if window <= 0:
        return samples[-1][col]
    sel = [s[col] for s in samples if s[0] > t_end - window - 1e-12]
    return sum(sel) / len(sel)


def reduce_run(rows, bulk, window):
    t = rows[-1][0]
    F = tuple(window_mean(rows, window, c) for c in (1, 2, 3))
    T = tuple(window_mean(rows, window, c) for c in (4, 5, 6))
    U = window_mean(bulk, window, 1)
    ffrac = window_mean(bulk, window, 2)
    fz_last, fz_prev = rows[-1][3], rows[-2][3] if len(rows) > 1 else rows[-1][3]
    return dict(t=t, F=F, T=T, U=U, fluid_frac=ffrac,
                plateau_dF=abs(fz_last - fz_prev) / max(abs(fz_last), 1e-30))


def rh_pair(r, solid_nominal):
    Fmag = math.sqrt(sum(f * f for f in r["F"]))
    rh = invert_rh(Fmag, r["U"])
    # d11_rh_collapse volume correction: dr from the run's own fluid_frac,
    # K divided by (1+dr)^1.87 == invert with F/(1+dr)^1.87
    dr = ((1.0 - r["fluid_frac"]) / solid_nominal) ** (1.0 / 3.0) - 1.0
    rhc = invert_rh(Fmag / (1.0 + dr) ** DLNK_DLNR, r["U"])
    return Fmag, rh, rhc, dr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rundirs", nargs="+")
    ap.add_argument("--h", type=float, default=1.0 / 36.0,
                    help="mesh spacing at the run level (L3 1/36, L4 1/72)")
    ap.add_argument("--scale", type=float,
                    default=float(os.environ.get("D61_SCALE", "1.0")))
    ap.add_argument("--window", type=float, default=0.5)
    ap.add_argument("--sensitivity", action="store_true")
    a = ap.parse_args()

    scale, h = a.scale, a.h
    b_minor = scale * R_SPHERE / 2.0 ** (1.0 / 3.0)
    a_major = 2.0 * b_minor
    solid_nominal = 0.019392547 * scale ** 3     # (4/3) pi a b^2 at this scale

    names = ["V0", "V1", "V2", "V3"][: len(a.rundirs)]
    data = {n: load_run(d) for n, d in zip(names, a.rundirs)}
    runs = {n: reduce_run(*data[n], a.window) for n in names}
    XA, YA = resistance_functions(RE_ASPECT)
    print(f"X^A={XA:.6f}  Y^A={YA:.6f}  Y/X={YA/XA:.5f}   "
          f"(scale {scale:g}, h {h:.6f}, window {a.window:g} t.u.)")
    print(f"{'run':4} {'t_end':>6} {'|F|':>12} {'U_sup':>12} {'R_h':>10} "
          f"{'|T|/(a|F|)':>11} {'plateau dF/F':>13}")
    rh, rhc = {}, {}
    for n in names:
        r = runs[n]
        Fmag, rh[n], rhc[n], dr = rh_pair(r, solid_nominal)
        Tmag = math.sqrt(sum(x * x for x in r["T"]))
        print(f"{n:4} {r['t']:6.2f} {Fmag:12.6e} {r['U']:12.6e} {rh[n]:10.6f} "
              f"{Tmag/(a_major*Fmag):11.2e} {r['plateau_dF']:13.2e}  dr={dr:+.4%} "
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
        # absolutes, raw and a_eff-corrected (a_eff = a - 0.14 h, axis-wise)
        a_e, b_e = a_major + AEFF_C * h, b_minor + AEFF_C * h
        XAe, YAe = resistance_functions(a_e / b_e)
        for n, coef, coef_e in (("V1", XA, XAe), ("V2", YA, YAe)):
            tgt_abs = a_major * coef
            tgt_eff = a_e * coef_e
            print(f"G-abs {n}: R_h raw {rh[n]:.6f} / vol-corr {rhc[n]:.6f} vs "
                  f"a*C = {tgt_abs:.6f} (raw {(rh[n]/tgt_abs-1)*100:+.2f}%, "
                  f"vol-corr {(rhc[n]/tgt_abs-1)*100:+.2f}%) "
                  f"| a_eff target {tgt_eff:.6f} (vol-corr {(rhc[n]/tgt_eff-1)*100:+.2f}%)  "
                  f"[band +-3%]")
    if "V3" in runs:
        fx, fy, fz = runs["V3"]["F"]
        axis = (1.0 / math.sqrt(2.0), 0.0, 1.0 / math.sqrt(2.0))
        Fmag = math.sqrt(fx * fx + fy * fy + fz * fz)
        cos_fa = (fx * axis[0] + fy * axis[1] + fz * axis[2]) / Fmag
        ang_axis = math.degrees(math.acos(max(-1.0, min(1.0, cos_fa))))
        print(f"G-offdiag [F-form, INVALID for body-force driving]: "
              f"angle(F,axis) = {ang_axis:.2f} deg (steady-state momentum "
              f"balance forces F -> f*V_cell, transverse F -> 0; the "
              f"off-diagonal mobility appears as a TRANSVERSE MEAN FLOW "
              f"U_x/U_z = (Y-X)/(Y+X) = {(YA-XA)/(YA+XA):+.4f} = "
              f"{math.degrees(math.atan((YA-XA)/(YA+XA))):+.2f} deg instead - "
              f"needs the 3-component bulk-flow diagnostic, gate redesigned)")

    if a.sensitivity and "V1" in runs and "V2" in runs:
        print("\n--- window sensitivity (trailing mean over W t.u.; 0 = last sample) ---")
        print(f"{'W':>5} {'R_h V1 raw':>11} {'R_h V2 raw':>11} {'ratio raw':>10} "
              f"{'ratio vol-corr':>15} {'dev vs Y/X':>11}")
        for w in SENS_WINDOWS:
            rr = {n: reduce_run(*data[n], w) for n in ("V1", "V2")}
            _, r1, r1c, _ = rh_pair(rr["V1"], solid_nominal)
            _, r2, r2c, _ = rh_pair(rr["V2"], solid_nominal)
            print(f"{w:5.2f} {r1:11.6f} {r2:11.6f} {r2/r1:10.5f} {r2c/r1c:15.5f} "
                  f"{(r2c/r1c/(YA/XA)-1)*100:+10.3f}%")
        # first-half vs second-half of the final 2 t.u.: a plateau has no trend
        t_end = data["V1"][0][-1][0]
        for n in ("V1", "V2"):
            rows = data[n][0]
            def mean_between(lo, hi):
                sel = [r[3] for r in rows if lo < r[0] <= hi]
                return sum(sel) / len(sel)
            f1, f2 = mean_between(t_end - 2.0, t_end - 1.0), mean_between(t_end - 1.0, t_end)
            print(f"{n}: F_z mean (t_end-2..t_end-1) = {f1:.6e}, (t_end-1..t_end) = {f2:.6e}, "
                  f"trend {(f2/f1-1)*100:+.3f}% per t.u.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main()
