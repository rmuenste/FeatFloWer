#!/usr/bin/env python3
"""D6.2 Jeffery-orbit gate analysis (CASE_SPEC d62_jeffery section 1).

Reads the DNS_PART_AXIS trace from a run log, reconstructs the in-plane
orientation angle phi(t) = unwrap(atan2(axis_z, axis_x)), and gates:

  G-period  : T*gammadot vs 2 pi (r_e + 1/r_e)          [15.70796 at r_e=2]
  G-waveform: |dphi/dt| extrema vs Jeffery gammadot/(r_e^2+1) .. gammadot r_e^2/(r_e^2+1);
              fast/slow modulation ratio vs r_e^2 (4.0 at r_e=2)
  G-orient  : the waveform AGAINST orientation - |dphi/dt| binned by phi mod pi
              (phi from the flow axis) vs gammadot (cos^2 + r_e^2 sin^2)/(r_e^2+1):
              rate(phi~0)/rate(phi~pi/2) vs ~1/r_e^2 (slow at flow alignment) and
              the rms residual over the bins. Added after the 2026-09-10 review:
              extrema alone pass a phase-shifted orbit (fast at alignment).
  sense     : dphi/dt sign vs the expected omega_y = +gammadot/2 (phi decreasing)
  G-plane   : max |axis_y| over the trace (in-plane stability)

r_e = 1 degenerates to the uniform spin control (V0): dphi/dt = gammadot/2,
T*gammadot = 4 pi.

Usage: d62_jeffery_analysis.py <run_slurm.log> --gammadot 0.2 [--re 2.0]
                               [--plot fig.png]

--plot draws the classic Jeffery validation figure: Jeffery's closed-form
solution as continuous curves with the DNS samples overlaid as markers -
axis components a_x(t), a_z(t); rotation rate dphi/dt / gammadot vs t; and
the rate waveform vs phi. The analytic curve uses the THEORETICAL period
(phase anchored at the first sample), so a period error shows as a growing
phase lag, exactly what the figure is meant to expose.

Convention (flow u = gammadot z along x, vorticity +gammadot along +y): phi measured
from the flow axis in the x-z plane, phi(t) = atan2(sin(psi)/r_e, cos(psi))
unwrapped, psi = -gammadot r_e t/(r_e^2+1) + psi0, i.e. tan(phi) =
tan(psi)/r_e; rate gammadot (cos^2 phi + r_e^2 sin^2 phi)/(r_e^2+1) - slow
at flow alignment, fast through the gradient direction.
"""
import argparse
import math
import re
import sys


def jeffery_phi(tt, t0, phi0, sign, re_, g):
    """Closed-form Jeffery angle (continuous unwrap), anchored at phi(t0)=phi0."""
    r = re_
    # psi0 from phi0: tan(psi) = r tan(phi), same half-turn as phi0
    psi0 = math.atan2(r * math.sin(phi0), math.cos(phi0))
    psi0 += math.pi * round((phi0 - psi0) / math.pi)
    out = []
    for x in tt:
        psi = psi0 + sign * g * r * (x - t0) / (r * r + 1.0)
        p = math.atan2(math.sin(psi), r * math.cos(psi))
        p += math.pi * round((psi - p) / math.pi)
        out.append(p)
    return out


def jeffery_rate(phi_, re_, g):
    r2 = re_ * re_
    return [g * (math.cos(p) ** 2 + r2 * math.sin(p) ** 2) / (r2 + 1.0) for p in phi_]


def make_plot(path, t, ax, az, phi, re_, g, Tg_meas, Tg_jeff):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = np.asarray(t); ax = np.asarray(ax); az = np.asarray(az); phi = np.asarray(phi)
    sign = 1.0 if phi[-1] >= phi[0] else -1.0
    tf = np.linspace(t[0], t[-1], 4000)
    phf = np.asarray(jeffery_phi(tf, t[0], phi[0], sign, re_, g))
    rf = np.asarray(jeffery_rate(phf, re_, g)) / g
    # measured centred rate
    rm = np.abs(np.gradient(phi, t)) / g
    # marker subsampling: ~200 markers per panel
    step = max(1, len(t) // 200)
    sl = slice(0, None, step)
    # analytic components must follow the DNS axis sign convention (the
    # recorded axis is continuous, phi is unwrapped from it, so cos/sin match)
    fig, axs = plt.subplots(3, 1, figsize=(8.5, 10.5))
    a0, a1, a2 = axs
    a0.plot(tf, np.cos(phf), "-", color="tab:blue", lw=1.2, label="Jeffery $a_x$")
    a0.plot(tf, np.sin(phf), "-", color="tab:red", lw=1.2, label="Jeffery $a_z$")
    a0.plot(t[sl], ax[sl], "o", ms=3.5, mfc="none", color="tab:blue", label="DNS $a_x$")
    a0.plot(t[sl], az[sl], "s", ms=3.5, mfc="none", color="tab:red", label="DNS $a_z$")
    a0.set_ylabel("axis component"); a0.set_ylim(-1.15, 1.15); a0.grid(alpha=.3)
    a0.legend(ncol=4, fontsize=8, loc="upper right")
    a1.plot(tf, rf, "-", color="k", lw=1.2, label="Jeffery")
    a1.plot(t[sl], rm[sl], "o", ms=3.5, mfc="none", color="tab:green", label="DNS")
    a1.set_ylabel(r"$\dot\varphi/\dot\gamma$"); a1.grid(alpha=.3); a1.legend(fontsize=8)
    r2 = re_ * re_
    a1.set_ylim(0, max(1.05 * r2 / (r2 + 1), 0.6))
    a1.set_xlabel("t")
    a0.set_xlabel("t")
    # waveform vs phi (mod pi, measured from the flow axis)
    pw = np.linspace(0, math.pi, 400)
    a2.plot(pw, np.asarray(jeffery_rate(pw, re_, g)) / g, "-", color="k", lw=1.2, label="Jeffery")
    a2.plot(np.mod(phi[sl], math.pi), rm[sl], "o", ms=3.5, mfc="none", color="tab:green", label="DNS")
    a2.set_xlabel(r"$\varphi$ mod $\pi$ (from the flow axis)"); a2.set_ylabel(r"$\dot\varphi/\dot\gamma$")
    a2.set_xlim(0, math.pi); a2.set_ylim(0, max(1.05 * r2 / (r2 + 1), 0.6)); a2.grid(alpha=.3)
    a2.set_xticks([0, math.pi / 4, math.pi / 2, 3 * math.pi / 4, math.pi])
    a2.set_xticklabels(["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$", r"$\pi$"]); a2.legend(fontsize=8)
    ttl = (f"D6.2 Jeffery orbit, $r_e$={re_:g}, $\\dot\\gamma$={g:g}: "
           f"$T\\dot\\gamma$ = {Tg_jeff:.4f} (Jeffery)")
    if Tg_meas is not None:
        ttl += f", {Tg_meas:.4f} measured ({(Tg_meas / Tg_jeff - 1) * 100:+.2f}%)"
    fig.suptitle(ttl, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=150)
    print(f"plot written: {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--gammadot", type=float, required=True)
    ap.add_argument("--re", type=float, default=2.0)
    ap.add_argument("--tmin", type=float, default=0.0,
                    help="discard the startup transient before this time")
    ap.add_argument("--plot", default=None,
                    help="write the Jeffery-vs-DNS overlay figure (png/pdf)")
    a = ap.parse_args()

    pat = re.compile(r"DNS_PART_AXIS time=\s*(\S+)\s+ip=\s*\d+\s+axis=\s*(\S+)\s+(\S+)\s+(\S+)")
    t, ax, ay, az = [], [], [], []
    with open(a.log) as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                tt = float(m.group(1))
                if tt < a.tmin:
                    continue
                t.append(tt)
                ax.append(float(m.group(2)))
                ay.append(float(m.group(3)))
                az.append(float(m.group(4)))
    if len(t) < 10:
        sys.exit(f"only {len(t)} DNS_PART_AXIS samples found (tmin={a.tmin})")

    # unwrap phi; the axis is headless (a and -a identical), so unwrap modulo pi
    phi = [math.atan2(az[0], ax[0])]
    for i in range(1, len(t)):
        p = math.atan2(az[i], ax[i])
        prev = phi[-1]
        while p - prev > math.pi / 2:
            p -= math.pi
        while p - prev < -math.pi / 2:
            p += math.pi
        phi.append(p)

    total = abs(phi[-1] - phi[0])
    span = t[-1] - t[0]
    n_half = total / math.pi
    # Period from pi-crossing times: by orbit symmetry each advance of phi by
    # pi takes exactly T/2, so T = 2 * mean spacing of successive crossings.
    # (A mean-rate estimate is biased whenever the trace covers a non-integer
    # number of half-turns - validated on a synthetic trace.)
    # direction-agnostic: fold a decreasing phi (negative rotation sense, the
    # D6.2 convention: omega_y = +gammadot/2 makes phi decrease) onto an increasing one
    sgn = 1.0 if phi[-1] >= phi[0] else -1.0
    ph = [sgn * p for p in phi]
    base = math.floor(ph[0] / math.pi)
    crossings = []
    for i in range(1, len(t)):
        while ph[i] - (base + 1) * math.pi >= 0:
            base += 1
            f = ((base) * math.pi - ph[i - 1]) / (ph[i] - ph[i - 1])
            crossings.append(t[i - 1] + f * (t[i] - t[i - 1]))
    if len(crossings) >= 2:
        gaps = [crossings[i + 1] - crossings[i] for i in range(len(crossings) - 1)]
        T_mean = 2.0 * sum(gaps) / len(gaps)
        print(f"period from {len(crossings)} pi-crossings "
              f"(half-period spread {min(gaps):.3f}..{max(gaps):.3f})")
    else:
        print(f"WARNING: {len(crossings)} pi-crossings - period from mean rate (biased)")
        T_mean = span / (total / (2 * math.pi))
    Tg = T_mean * a.gammadot
    Tg_jeff = 2 * math.pi * (a.re + 1.0 / a.re)

    # waveform: centered dphi/dt vs phi
    rates, phis = [], []
    for i in range(1, len(t) - 1):
        dt = t[i + 1] - t[i - 1]
        if dt <= 0:
            continue
        rates.append(abs((phi[i + 1] - phi[i - 1]) / dt))
        phis.append(phi[i])
    rmax, rmin = max(rates), min(rates)
    mod_meas = rmax / rmin if rmin > 0 else float("inf")
    mod_jeff = a.re ** 2

    # predicted rate extremes
    g, r2 = a.gammadot, a.re ** 2
    fast, slow = g * r2 / (r2 + 1), g / (r2 + 1)

    # G-orient: the waveform AGAINST ORIENTATION (review 2026-09-10 finding 1:
    # extrema and their ratio alone cannot tell a phase-shifted orbit - one
    # that is fastest at flow alignment - from the real thing). Bin the
    # measured |dphi/dt| by phi mod pi (phi from the flow axis) and compare
    # bin means with Jeffery's rate at the bin centre; report the rms
    # residual and the slow/fast placement explicitly.
    nbin = 12
    sums = [0.0] * nbin
    cnts = [0] * nbin
    for p, r in zip(phis, rates):
        k = int((p % math.pi) / math.pi * nbin) % nbin
        sums[k] += r
        cnts[k] += 1
    resid2, nres = 0.0, 0
    bin_meas = {}
    for k in range(nbin):
        if cnts[k] == 0:
            continue
        pc = (k + 0.5) * math.pi / nbin
        r_j = g * (math.cos(pc) ** 2 + r2 * math.sin(pc) ** 2) / (r2 + 1.0)
        bin_meas[k] = sums[k] / cnts[k]
        resid2 += (bin_meas[k] - r_j) ** 2
        nres += 1
    rms_resid = math.sqrt(resid2 / nres) / g if nres else float("nan")
    # slow-at-alignment placement: mean rate in the two bins around phi=0
    # (mod pi) over the two bins around phi=pi/2; Jeffery expects ~1/r_e^2
    def bins_mean(keys):
        vals = [bin_meas[k] for k in keys if k in bin_meas]
        return sum(vals) / len(vals) if vals else float("nan")
    r_align = bins_mean((0, nbin - 1))
    r_grad = bins_mean((nbin // 2 - 1, nbin // 2))
    pc0, pc1 = 0.5 * math.pi / nbin, (nbin // 2 - 0.5) * math.pi / nbin
    place_jeff = ((math.cos(pc0) ** 2 + r2 * math.sin(pc0) ** 2)
                  / (math.cos(pc1) ** 2 + r2 * math.sin(pc1) ** 2))
    place_meas = r_align / r_grad if r_grad else float("nan")
    # signed rotation sense: with u = gammadot*z along x the vorticity is
    # +gammadot*y and phi = atan2(a_z, a_x) DECREASES (omega_y = +gammadot/2
    # for a sphere <-> dphi/dt = -gammadot/2)
    sense = -1.0 if (phi[-1] - phi[0]) < 0 else 1.0
    sense_ok = (sense < 0) == (a.gammadot > 0)

    ay_max = max(abs(v) for v in ay)

    print(f"samples {len(t)}  t=[{t[0]:.3f},{t[-1]:.3f}]  half-turns {n_half:.2f}  "
          f"rotation sense dphi/dt {'<' if sense < 0 else '>'} 0 "
          f"({'as expected' if sense_ok else 'WRONG SIGN'} for gammadot={a.gammadot:g}: "
          f"omega_y = +gammadot/2 -> phi decreasing)")
    print(f"G-period  : T*gammadot = {Tg:.4f} vs {Tg_jeff:.4f} "
          f"({(Tg/Tg_jeff-1)*100:+.2f}%)  [band +-3%]")
    print(f"G-waveform: |dphi/dt| in [{rmin:.5f},{rmax:.5f}] vs Jeffery "
          f"[{slow:.5f},{fast:.5f}]; modulation {mod_meas:.3f} vs {mod_jeff:.3f} "
          f"({(mod_meas/mod_jeff-1)*100:+.2f}%)  [band +-5%]")
    if math.isnan(place_meas):
        place_txt = ("rate(phi~0)/rate(phi~pi/2) = n/a (trace does not cover both "
                     "the aligned and the gradient orientation)")
    else:
        place_txt = (f"rate(phi~0)/rate(phi~pi/2) = {place_meas:.4f} vs Jeffery "
                     f"{place_jeff:.4f} ({(place_meas/place_jeff-1)*100:+.2f}%; slow at "
                     f"flow alignment) [band +-10%]")
    print(f"G-orient  : {place_txt}; rms waveform residual over {nres} phi-bins = "
          f"{rms_resid*100:.2f}% of gammadot [band < 3%]")
    print(f"G-plane   : max|axis_y| = {ay_max:.3e}  [band < 0.02]")

    if a.plot:
        make_plot(a.plot, t, ax, az, phi, a.re, a.gammadot,
                  Tg if len(crossings) >= 2 else None, Tg_jeff)


if __name__ == "__main__":
    main()
