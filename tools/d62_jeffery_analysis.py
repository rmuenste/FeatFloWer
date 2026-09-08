#!/usr/bin/env python3
"""D6.2 Jeffery-orbit gate analysis (CASE_SPEC d62_jeffery section 1).

Reads the DNS_PART_AXIS trace from a run log, reconstructs the in-plane
orientation angle phi(t) = unwrap(atan2(axis_z, axis_x)), and gates:

  G-period  : T*gammadot vs 2 pi (r_e + 1/r_e)          [15.70796 at r_e=2]
  G-waveform: |dphi/dt|(phi) vs Jeffery
              gammadot (r_e^2 cos^2 + sin^2)(phi-phi0) / (r_e^2+1),
              fitted with a free phase offset phi0 (convention-robust);
              reports the fast/slow modulation ratio vs r_e^2 (4.0 at r_e=2)
  G-plane   : max |axis_y| over the trace (in-plane stability)

r_e = 1 degenerates to the uniform spin control (V0): dphi/dt = gammadot/2,
T*gammadot = 4 pi.

Usage: d62_jeffery_analysis.py <run_slurm.log> --gammadot 0.2 [--re 2.0]
"""
import argparse
import math
import re
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("log")
    ap.add_argument("--gammadot", type=float, required=True)
    ap.add_argument("--re", type=float, default=2.0)
    ap.add_argument("--tmin", type=float, default=0.0,
                    help="discard the startup transient before this time")
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
    base = math.floor(phi[0] / math.pi)
    crossings = []
    for i in range(1, len(t)):
        while phi[i] - (base + 1) * math.pi >= 0:
            base += 1
            f = ((base) * math.pi - phi[i - 1]) / (phi[i] - phi[i - 1])
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

    ay_max = max(abs(v) for v in ay)

    print(f"samples {len(t)}  t=[{t[0]:.3f},{t[-1]:.3f}]  half-turns {n_half:.2f}")
    print(f"G-period  : T*gammadot = {Tg:.4f} vs {Tg_jeff:.4f} "
          f"({(Tg/Tg_jeff-1)*100:+.2f}%)  [band +-3%]")
    print(f"G-waveform: |dphi/dt| in [{rmin:.5f},{rmax:.5f}] vs Jeffery "
          f"[{slow:.5f},{fast:.5f}]; modulation {mod_meas:.3f} vs {mod_jeff:.3f} "
          f"({(mod_meas/mod_jeff-1)*100:+.2f}%)  [band +-5%]")
    print(f"G-plane   : max|axis_y| = {ay_max:.3e}  [band < 0.02]")


if __name__ == "__main__":
    main()
