#!/usr/bin/env python3
"""D3.1 steady-state extrapolation pipeline (CHERD supplementary numbers).

Finite-Re array runs at nu=0.01 approach steady state exponentially with
relaxation time tau ~ rho*U_ss/(f*(1-phi)); the dilute bands (phi<=0.10)
are up to ~13% below steady velocity at the t=4 endpoint. This tool
extrapolates every run's superficial velocity to steady state via the
geometric tail

    U_ss = U(te) + dU * r/(1-r),  r = [U(te)-U(te-s)]/[U(te-s)-U(te-2s)]

(s = 0.5 t.u.; cross-checked against s = 1.0), then recomputes:
  - Table S5: F* (steady, band mean +- seed sd) vs Beetstra eq17 at the
    band-mean steady Re  [deviation convention: mean-F* at mean-Re]
  - corner refinement shifts L4->L5 (steady, same seed/config)
  - Table S6: production deficits L2 vs L4 (steady), raw fixed-forcing
    and matched-Re (L2 point moved to the L4 steady Re along eq17)

Force convention: momentum balance F_d = (1-phi) f V / N throughout.
"""
import json, math, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def deck_val(rd, key):
    for line in open(f"{rd}/_data/q2p1_param.dat"):
        s = line.split("=")
        if len(s) == 2 and s[0].strip().endswith(key):
            return s[1].strip()

def bulk(rd):
    rows = []
    for line in open(f"{rd}/bulk_flow.log"):
        if line.startswith("#"):
            continue
        p = line.split()
        rows.append((float(p[0]), float(p[1])))
    return rows

def u_at(rows, t):
    return min(rows, key=lambda x: abs(x[0] - t))[1]

def u_steady(rd, spacing=0.5):
    rows = bulk(rd)
    te = rows[-1][0]
    u1, u2, u3 = (u_at(rows, te - 2 * spacing), u_at(rows, te - spacing),
                  u_at(rows, te))
    r = (u3 - u2) / (u2 - u1) if u2 != u1 else 0.0
    uss = u3 + (u3 - u2) * r / (1 - r) if 0 < r < 1 else u3
    v1, v2, v3 = u_at(rows, te - 4 * spacing), u_at(rows, te - 2 * spacing), u3
    r2 = (v3 - v2) / (v2 - v1) if v2 != v1 else 0.0
    uss2 = v3 + (v3 - v2) * r2 / (1 - r2) if 0 < r2 < 1 else v3
    return uss, abs(uss - uss2) / uss, (uss - u3) / uss

def eq6(p):
    return 10 * p / (1 - p) ** 2 + (1 - p) ** 2 * (1 + 1.5 * math.sqrt(p))

def eq17(p, Re):
    if Re < 1e-3:
        return eq6(p)
    add = (0.413 * Re / (24 * (1 - p) ** 2)
           * (((1 - p) ** -1 + 3 * p * (1 - p) + 8.4 * Re ** -0.343)
              / (1 + 10 ** (3 * p) * Re ** (-(1 + 4 * p) / 2))))
    return eq6(p) + add

def run_data(name, phi, N=16):
    rd = os.path.join(ROOT, f"q2p1_dns_rundir_{name}")
    f = float(deck_val(rd, "SimPar@ConstantForcing").split(",")[2].replace("d", "e"))
    mu = float(deck_val(rd, "Prop@Viscosity").split(",")[0].replace("d", "e"))
    d = 2 * json.load(open(f"{rd}/example.json"))["benchRadius_"]
    uss, unc, gap = u_steady(rd)
    Fs = (1 - phi) * f / N / (3 * math.pi * mu * d * uss)
    return {"F": Fs, "Re": uss * d / mu, "unc": unc, "gap": gap, "d": d,
            "mu": mu, "f": f, "uss": uss}

def l4_name(pp, s, tier):
    if pp == "p005" and s == "s1" and tier == "":
        return "d31_l4_img"
    return f"d31_{pp}_{s}{tier}"

def main():
    phis = [("p005", 0.05), ("p010", 0.10), ("p020", 0.20), ("p030", 0.30)]
    tiers = [("", "Stokes"), ("_re10", "Re~9"), ("_re30", "Re~27")]

    print("=== Table S5 (steady): band mean F* +- sd, dev vs eq17(mean Re) ===")
    worst_gap = worst_unc = 0.0
    for pp, phi in phis:
        cells = []
        for tier, lab in tiers:
            R = [run_data(l4_name(pp, s, tier), phi) for s in ("s1", "s2", "s3")]
            mF = sum(r["F"] for r in R) / 3
            sd = math.sqrt(sum((r["F"] - mF) ** 2 for r in R) / 3)
            mRe = sum(r["Re"] for r in R) / 3
            ref = eq17(phi, mRe if tier else 0)
            for r in R:
                worst_gap = max(worst_gap, abs(r["gap"]))
                worst_unc = max(worst_unc, r["unc"])
            cells.append(f"{lab}: {mF:.2f}+-{sd:.2f} ({mF/ref-1:+.1%}) "
                         f"Re {min(r['Re'] for r in R):.3g}-{max(r['Re'] for r in R):.3g}")
        print(f"phi={phi}: " + " | ".join(cells))
    print(f"worst endpoint-to-steady gap {worst_gap:.1%}, worst extrap disagreement {worst_unc:.1%}")

    print("\n=== corner refinement shifts L4->L5 (steady, seed 1) ===")
    for pp, phi in [("p005", 0.05), ("p010", 0.10), ("p020", 0.20)]:
        r4 = run_data(f"d31_{pp}_s1_re30", phi)
        r5 = run_data(f"d31_l5_{pp}_re30", phi)
        print(f"phi={phi}: F*4={r4['F']:.2f} -> F*5={r5['F']:.2f}  shift {r5['F']/r4['F']-1:+.1%}"
              f"  dev5 vs eq17 {r5['F']/eq17(phi, r5['Re'])-1:+.1%}")

    print("\n=== Table S6 (steady): L2 vs L4, raw and matched-Re ===")
    for pp, phi in [("p005", 0.05), ("p010", 0.10)]:
        for tier, lab in tiers:
            raw, mat, re4r, re2r = [], [], [], []
            for s in ("s1", "s2", "s3"):
                r4 = run_data(l4_name(pp, s, tier), phi)
                r2 = run_data(f"d31_{pp}_{s}{tier}_l2", phi)
                dr = r2["F"] / r4["F"] - 1
                # matched-Re transport: multiplicative, along the correlation's
                # relative Re-dependence (consistent with the deficit acting as a
                # factor on drag); see CHERD SI eq. si_matched_re
                ratio = (eq17(phi, r4["Re"]) / eq17(phi, r2["Re"])) if tier else 1.0
                dm = (r2["F"] * ratio) / r4["F"] - 1
                raw.append(dr); mat.append(dm)
                re4r.append(r4["Re"]); re2r.append(r2["Re"])
            print(f"{pp} {lab:6s}: raw {'/'.join(f'{x:+.1%}' for x in raw)} (mean {sum(raw)/3:+.1%})"
                  f" | matched-Re {'/'.join(f'{x:+.1%}' for x in mat)} (mean {sum(mat)/3:+.1%})"
                  f" | Re4 {min(re4r):.3g}-{max(re4r):.3g}, Re2 {min(re2r):.3g}-{max(re2r):.3g}")

    print("\n=== N=54 (steady = endpoint, Stokes) ===")
    for name, lvl in [("d31_n54_l4", "L4"), ("d31_n54_l5", "L5")]:
        r = run_data(name, 0.10, N=54)
        print(f"{lvl}: F*={r['F']:.3f}  dev vs eq6 {r['F']/eq6(0.10)-1:+.1%}")

if __name__ == "__main__":
    main()
