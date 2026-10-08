#!/usr/bin/env python3
"""Happel & Brenner (1973) Stokes resistance of spheroids: table + closed form.

Source, read from the scanned page images (NOT from the OCR text layer, which
garbles the formulae) of `literature/happel_brenner_1973.pdf`:

  * Table 4-26.1 (PDF page 162 = book page 149) "Resistance of oblate and
    prolate spheroids expressed in terms of the Stokes' law correction factor
    for a sphere having the SAME EQUATORIAL RADIUS" - axisymmetric flow only,
    i.e. motion PARALLEL to the symmetry axis (broadside for an oblate body).
  * Table 5-11.1 (PDF page 236 = book page 223) "Values of equivalent radius
    for an ellipsoid of revolution", both orientations, as R/(equatorial
    radius), where R is the radius of the sphere of equal Stokes resistance:
    F = 6 pi mu R U.
  * Equations 5-11.18 / 5-11.20 (parallel, prolate / oblate) and
    5-11.22 / 5-11.24 (perpendicular, prolate / oblate), PDF pages 235 and 237.

Conventions used here, matching the book:

  parallel      flow along the symmetry axis; semiaxes b = c (equatorial
                radius c), polar semiaxis a, aspect ratio phi = a/c.
                phi < 1 oblate (BROADSIDE settling), phi > 1 prolate.
  perpendicular flow across the symmetry axis; semiaxes a = b (equatorial
                radius a), polar semiaxis c, aspect ratio phi = c/a.
                phi < 1 oblate (EDGEWISE settling), phi > 1 prolate.

The campaign writes the Stokes balance as u_t = (rho_p - rho_f) g V /
(3 pi mu D_eq K), i.e. K is referred to the EQUAL-VOLUME sphere of diameter
D_eq, while the book's R and K are referred to the equatorial radius. The
conversion is K_eq_volume = R / a_eq with a_eq = D_eq/2; `k_equal_volume()`
does it.

Run with no arguments for the self-test against both printed tables and the
D6.4 tier-S evaluation.
"""

import math

PDF = "literature/happel_brenner_1973.pdf"

# --- Table 4-26.1, PDF p.162 / book p.149 -----------------------------------
# K = R/a for flow PARALLEL to the symmetry axis, per equatorial radius.
# b/a is polar/equatorial. The prolate column is K for the prolate spheroid of
# the same b/a read as equatorial/polar (see the book's wording); it is kept
# here only as printed, the self-test uses the oblate column.
TABLE_4_26_1 = {
    # b/a : (K oblate, K prolate)
    0.0: (0.84882639, math.inf),
    0.1: (0.85245060, 2.6471358),
    0.2: (0.86145221, 1.7848095),
    0.3: (0.87394886, 1.4697413),
    0.4: (0.88880656, 1.3050489),
    0.5: (0.90530533, 1.2039411),
    0.6: (0.92296815, 1.1358194),
    0.7: (0.94146887, 1.0870324),
    0.8: (0.96057733, 1.0505422),
    0.9: (0.98012819, 1.0223468),
    1.0: (1.0, 1.0),
}

# --- Table 5-11.1, PDF p.236 / book p.223 -----------------------------------
# Perpendicular to the axis, R/a, discoid (oblate) branch, phi = c/a < 1:
T511_PERP_OBLATE = {
    1e-6: 0.5659, 1e-5: 0.5659, 1e-4: 0.5659, 1e-3: 0.5664, 1e-2: 0.5707,
    1e-1: 0.6133, 1.5e-1: 0.6366, 2e-1: 0.6596, 3e-1: 0.7049, 4e-1: 0.7492,
    5e-1: 0.7927, 6e-1: 0.8355, 7e-1: 0.8775, 8e-1: 0.9189, 9e-1: 0.9597,
    9.9e-1: 0.9960,
}
# Perpendicular to the axis, R/a, needlelike (prolate) branch, phi = c/a > 1:
T511_PERP_PROLATE = {
    1.01: 1.004, 1.05: 1.020, 1.1: 1.040, 1.5: 1.194, 2: 1.379, 5: 2.371,
    10: 3.812, 20: 6.365, 50: 13.06, 1e2: 23.00, 2e2: 41.08, 5e2: 90.00,
    1e3: 164.6, 1e4: 1281.0, 1e5: 10493.0, 1e6: 88838.0,
}
# Parallel to the axis, R/c, discoid (oblate) branch, phi = a/c < 1:
T511_PAR_OBLATE = {
    1e-6: 0.8488, 1e-5: 0.8488, 1e-4: 0.8488, 1e-3: 0.8488, 1e-2: 0.8489,
    1e-1: 0.8525, 1.5e-1: 0.8564, 2e-1: 0.8615, 3e-1: 0.8739, 4e-1: 0.8888,
    5e-1: 0.9053, 6e-1: 0.9230, 7e-1: 0.9415, 8e-1: 0.9606, 9e-1: 0.9801,
    9.9e-1: 0.9980,
}
# Parallel to the axis, R/c, needlelike (prolate) branch, phi = a/c > 1:
T511_PAR_PROLATE = {
    1.01: 1.002, 1.05: 1.010, 1.10: 1.020, 1.50: 1.102, 2.0: 1.204, 5.0: 1.785,
    10: 2.647, 20: 4.172, 50: 8.117, 1e2: 13.895, 2e2: 24.280, 5e2: 52.022,
    1e3: 93.881, 1e4: 708.92, 1e5: 5695.2, 1e6: 47590.0,
}


def r_over_equatorial_parallel(phi):
    """R/c for flow parallel to the symmetry axis; phi = a/c (polar/equatorial).

    Eq. (5-11.20) for phi < 1, Eq. (5-11.18) for phi > 1, Eq. (5-11.21) at 0.
    """
    if phi <= 0.0:
        return 8.0 / (3.0 * math.pi)                      # (5-11.21) disk
    if abs(phi - 1.0) < 1e-12:
        return 1.0
    if phi < 1.0:
        s = math.sqrt(1.0 - phi * phi)
        denom = (2.0 * phi / (1.0 - phi * phi)
                 + 2.0 * (1.0 - 2.0 * phi * phi) / (1.0 - phi * phi) ** 1.5
                 * math.atan2(s, phi))
    else:
        s = math.sqrt(phi * phi - 1.0)
        denom = (-2.0 * phi / (phi * phi - 1.0)
                 + (2.0 * phi * phi - 1.0) / (phi * phi - 1.0) ** 1.5
                 * math.log((phi + s) / (phi - s)))
    return (8.0 / 3.0) / denom


def r_over_equatorial_perpendicular(phi):
    """R/a for flow perpendicular to the symmetry axis; phi = c/a (polar/equatorial).

    Eq. (5-11.24) for phi < 1, Eq. (5-11.22) for phi > 1, Eq. (5-11.25) at 0.
    """
    if phi <= 0.0:
        return 16.0 / (9.0 * math.pi)                     # (5-11.25) disk
    if abs(phi - 1.0) < 1e-12:
        return 1.0
    if phi < 1.0:
        s = math.sqrt(1.0 - phi * phi)
        denom = (-phi / (1.0 - phi * phi)
                 - (2.0 * phi * phi - 3.0) / (1.0 - phi * phi) ** 1.5
                 * math.asin(s))
    else:
        s = math.sqrt(phi * phi - 1.0)
        denom = (phi / (phi * phi - 1.0)
                 + (2.0 * phi * phi - 3.0) / (phi * phi - 1.0) ** 1.5
                 * math.log(phi + s))
    return (8.0 / 3.0) / denom


def k_equal_volume(polar, equatorial, orientation):
    """Stokes correction K referred to the EQUAL-VOLUME sphere.

    F = 6 pi mu a_eq K U = 3 pi mu D_eq K U, with a_eq the equal-volume radius
    of the spheroid (semiaxes polar x equatorial x equatorial).
    """
    phi = polar / equatorial
    if orientation == "parallel":          # broadside for an oblate body
        r = r_over_equatorial_parallel(phi) * equatorial
    elif orientation == "perpendicular":   # edgewise for an oblate body
        r = r_over_equatorial_perpendicular(phi) * equatorial
    else:
        raise ValueError(orientation)
    a_eq = (polar * equatorial * equatorial) ** (1.0 / 3.0)
    return r / a_eq, r, a_eq


def _selftest():
    print(f"self-test of the closed forms against the printed tables ({PDF})")
    worst = 0.0
    for name, table, fn in (
        ("Table 5-11.1 parallel, oblate  (R/c)", T511_PAR_OBLATE, r_over_equatorial_parallel),
        ("Table 5-11.1 parallel, prolate (R/c)", T511_PAR_PROLATE, r_over_equatorial_parallel),
        ("Table 5-11.1 perp,     oblate  (R/a)", T511_PERP_OBLATE, r_over_equatorial_perpendicular),
        ("Table 5-11.1 perp,     prolate (R/a)", T511_PERP_PROLATE, r_over_equatorial_perpendicular),
    ):
        dev = max(abs(fn(phi) - val) / val for phi, val in table.items())
        worst = max(worst, dev)
        print(f"  {name}: {len(table):2d} entries, max relative deviation {dev:.2e}")
    dev = max(abs(r_over_equatorial_parallel(ba) - k) / k
              for ba, (k, _) in TABLE_4_26_1.items() if ba > 0.0)
    worst = max(worst, dev)
    print(f"  Table 4-26.1 oblate K (8 digits): max relative deviation {dev:.2e}")
    print(f"  worst overall: {worst:.2e}")
    return worst


def _d64_tier_s():
    # Geometry and fluid of the D6.4 tier-S oblate pair (CASE_SPEC.md 2,
    # datasheet row d64_s_oblate_r3_result).
    polar, equatorial = 0.2404, 0.7211
    mu, D_eq = 0.374, 1.0                       # rho_f = 1, nu = 0.374
    weight = 0.0733                             # (rho_p - rho_f) g V, rho_r 3, g 0.07
    meas = {"parallel": 0.011582, "perpendicular": 0.015151}   # 8 D_eq box

    print(f"\nD6.4 tier-S oblate: semiaxes {polar} x {equatorial} x {equatorial}, "
          f"c/a = {polar / equatorial:.5f}, mu = {mu}, net weight = {weight}")
    out = {}
    for orientation, label in (("parallel", "broadside"), ("perpendicular", "edgewise")):
        k, r, a_eq = k_equal_volume(polar, equatorial, orientation)
        u = weight / (6.0 * math.pi * mu * r)
        out[orientation] = (k, r, u)
        print(f"  {label:9s} ({orientation:13s}): R = {r:.6f}  a_eq = {a_eq:.6f}  "
              f"K = {k:.6f}  u_t(unbounded) = {u:.6f}  measured(8 D_eq) = {meas[orientation]:.6f}  "
              f"{(meas[orientation] / u - 1.0) * 100:+.1f} %")
    ratio_theory = out['parallel'][2] / out['perpendicular'][2]
    ratio_meas = meas['parallel'] / meas['perpendicular']
    print(f"  orientation ratio broadside/edgewise: theory {ratio_theory:.4f}  "
          f"measured {ratio_meas:.4f}  {(ratio_meas / ratio_theory - 1.0) * 100:+.2f} % "
          f"(gate +-2 %)")
    return out


if __name__ == "__main__":
    _selftest()
    _d64_tier_s()
