#!/usr/bin/env python3
"""Generate a coarse cubed-sphere SHELL mesh (.tri) -- the atmosphere of a sphere.

Design: chimera-integration-design.md v3, section 7 (tools/chimera_meshgen),
Phase 5 (sphere arrays).

The shell r_i <= r <= r_o is meshed by projecting the six faces of the
cube [-1,1]^3 radially onto the sphere (equiangular cubed-sphere map:
the face coordinates are tan(theta) with theta uniform in [-pi/4, pi/4],
which gives near-uniform surface cells) and stacking `nr` radial layers.
All vertices are shared between neighbouring cube faces (conforming,
closed shell: the only boundary faces are the inner and the outer
sphere, which is what chi_submesh's geometric classification expects).

Hex vertex ordering follows the FEAT convention (1..4 bottom, 5..8 top,
positive Jacobian).  Inner/outer vertices lie EXACTLY on the spheres
(direction vectors are normalised, radii applied exactly); the level-1
classification in chi_submesh.f90 relies on that.  At load time
CHI_SUBMESH_FIT_COARSE rescales [r_i, r_o] affinely to the per-body
[R, R+H] and translates to the body centre, so one shell file serves an
entire array with per-particle atmosphere widths.

Element count: 6 * n^2 * nr.  Radial grading: the layer thickness grows
geometrically from the inner surface by the factor `--grading` per layer
(1.0 = uniform).

Usage:
  sphere_shell_tri.py --ri 1.0 --ro 2.0 --n 4 --nr 2 --out sphere_shell_coarse.tri
"""

import argparse
import math
import sys

from annulus_tri import write_tri


def generate_sphere_shell(ri, ro, n, nr, grading=1.0):
    """Return (coords, kvert, knpr) for the cubed-sphere shell."""
    if not (ro > ri > 0.0):
        raise ValueError("need 0 < ri < ro")
    if n < 1 or nr < 1:
        raise ValueError("need n >= 1, nr >= 1")
    if grading <= 0.0:
        raise ValueError("need grading > 0")

    # radial layer radii (geometric grading from the inner surface)
    if abs(grading - 1.0) < 1e-14:
        radii = [ri + (ro - ri) * l / nr for l in range(nr + 1)]
    else:
        total = sum(grading ** l for l in range(nr))
        h0 = (ro - ri) / total
        radii = [ri]
        for l in range(nr):
            radii.append(radii[-1] + h0 * grading ** l)
    radii[0] = ri
    radii[-1] = ro

    # equiangular cube coordinate for lattice index i in 0..n
    def tcoord(i):
        if i == 0:
            return -1.0
        if i == n:
            return 1.0
        return math.tan(-0.25 * math.pi + 0.5 * math.pi * i / n)

    # unique surface lattice points: (i,j,k) with at least one index in {0,n}
    surf_id = {}
    directions = []
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                if not (i in (0, n) or j in (0, n) or k in (0, n)):
                    continue
                c = (tcoord(i), tcoord(j), tcoord(k))
                nrm = math.sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2])
                surf_id[(i, j, k)] = len(directions)
                directions.append((c[0] / nrm, c[1] / nrm, c[2] / nrm))
    nsurf = len(directions)

    coords = []
    knpr = []
    for l, r in enumerate(radii):
        for d in directions:
            coords.append((r * d[0], r * d[1], r * d[2]))
            knpr.append(1 if (l == 0 or l == nr) else 0)

    def vid(key, l):
        return l * nsurf + surf_id[key] + 1  # 1-based

    # the six cube faces: (fixed axis, fixed index, free axes)
    faces = []
    for axis in range(3):
        for side in (0, n):
            free = [a for a in range(3) if a != axis]
            faces.append((axis, side, free[0], free[1]))

    def key(axis, side, a, b, ia, ib):
        idx = [0, 0, 0]
        idx[axis] = side
        idx[a] = ia
        idx[b] = ib
        return tuple(idx)

    kvert = []
    for axis, side, a, b in faces:
        for ia in range(n):
            for ib in range(n):
                quad = [key(axis, side, a, b, ia, ib),
                        key(axis, side, a, b, ia + 1, ib),
                        key(axis, side, a, b, ia + 1, ib + 1),
                        key(axis, side, a, b, ia, ib + 1)]
                for l in range(nr):
                    bot = [vid(q, l) for q in quad]
                    top = [vid(q, l + 1) for q in quad]
                    hexa = bot + top
                    if jacobian_sign(coords, hexa) < 0.0:
                        hexa = [bot[0], bot[3], bot[2], bot[1],
                                top[0], top[3], top[2], top[1]]
                    if jacobian_sign(coords, hexa) <= 0.0:
                        raise RuntimeError("could not orient a hex positively")
                    kvert.append(tuple(hexa))
    return coords, kvert, knpr


def jacobian_sign(coords, hexa):
    """Sign of (v2-v1) x (v4-v1) . (v5-v1) for FEAT vertex order (1-based ids)."""
    p1 = coords[hexa[0] - 1]
    p2 = coords[hexa[1] - 1]
    p4 = coords[hexa[3] - 1]
    p5 = coords[hexa[4] - 1]
    e1 = [p2[i] - p1[i] for i in range(3)]
    e2 = [p4[i] - p1[i] for i in range(3)]
    e3 = [p5[i] - p1[i] for i in range(3)]
    cx = e1[1] * e2[2] - e1[2] * e2[1]
    cy = e1[2] * e2[0] - e1[0] * e2[2]
    cz = e1[0] * e2[1] - e1[1] * e2[0]
    return cx * e3[0] + cy * e3[1] + cz * e3[2]


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ri", type=float, required=True, help="inner radius")
    ap.add_argument("--ro", type=float, required=True, help="outer radius")
    ap.add_argument("--n", type=int, required=True, help="cells per cube edge")
    ap.add_argument("--nr", type=int, required=True, help="radial layers")
    ap.add_argument("--grading", type=float, default=1.0,
                    help="radial growth factor per layer (1 = uniform)")
    ap.add_argument("--out", required=True, help="output .tri path")
    a = ap.parse_args(argv)

    coords, kvert, knpr = generate_sphere_shell(a.ri, a.ro, a.n, a.nr, a.grading)
    comment = (f"chimera_meshgen sphere shell ri={a.ri} ro={a.ro} n={a.n} "
               f"nr={a.nr} grading={a.grading}")
    write_tri(a.out, coords, kvert, knpr, comment)
    print(f"wrote {a.out}: nel={len(kvert)} nvt={len(coords)} "
          f"boundary vertices={sum(knpr)}")


if __name__ == "__main__":
    main(sys.argv[1:])
