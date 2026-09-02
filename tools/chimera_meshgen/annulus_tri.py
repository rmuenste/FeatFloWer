#!/usr/bin/env python3
"""Generate coarse body-fitted shell meshes (.tri) for the Chimera component.

Design: chimera-integration-design.md v3, section 7 (tools/chimera_meshgen).

Currently supported shape: a z-aligned annulus (O-grid ring, extruded) --
the atmosphere of a cylinder along the z-axis. Node/element ordering uses
the FeatFloWer/FEAT conventions (hex vertex ordering 1..4 bottom CCW,
5..8 top; positive Jacobian). The theta direction closes periodically.

The .tri format matches readTriCoarse (source/src_mesh/mesh_refine.f90):
  line 1: comment
  line 2: comment
  line 3: NEL NVT NBCT NVE NEE NAE
  'DCORVG' + NVT coordinate lines
  'KVERT'  + NEL 8-vertex lines
  'KNPR'   + NVT nodal-property lines (1 = boundary node, 0 = interior)

The generator places all boundary vertices EXACTLY on the analytic
surfaces; the classification/projection in chi_submesh.f90 relies on
that (level-1 classification by exact radius, hierarchical propagation
during refinement).

Usage:
  annulus_tri.py --ri 1.0 --ro 2.0 --lz 1.0 --nr 2 --nt 12 --nz 2 \
                 --out annulus_coarse.tri
"""

import argparse
import math
import sys


def generate_annulus(ri, ro, lz, nr, nt, nz):
    """Return (coords, kvert, knpr) for the annulus O-grid."""
    if not (ro > ri > 0.0 and lz > 0.0):
        raise ValueError("need 0 < ri < ro and lz > 0")
    if nr < 1 or nt < 3 or nz < 1:
        raise ValueError("need nr >= 1, nt >= 3, nz >= 1")

    nvt = (nr + 1) * nt * (nz + 1)

    def vid(ir, it, iz):
        # 1-based FeatFloWer vertex id; theta wraps periodically
        return 1 + ir + (nr + 1) * ((it % nt) + nt * iz)

    coords = [None] * nvt
    knpr = [0] * nvt
    for iz in range(nz + 1):
        z = lz * iz / nz
        for it in range(nt):
            th = 2.0 * math.pi * it / nt
            for ir in range(nr + 1):
                r = ri + (ro - ri) * ir / nr
                i = vid(ir, it, iz) - 1
                coords[i] = (r * math.cos(th), r * math.sin(th), z)
                if ir == 0 or ir == nr or iz == 0 or iz == nz:
                    knpr[i] = 1

    kvert = []
    for iz in range(nz):
        for it in range(nt):
            for ir in range(nr):
                # FEAT ordering: bottom (r,t),(r+1,t),(r+1,t+1),(r,t+1),
                # then the same on the top face. local xi1 ~ radial,
                # xi2 ~ theta, xi3 ~ z: e_r x e_theta = +e_z => detJ > 0.
                kvert.append((
                    vid(ir,     it,     iz),
                    vid(ir + 1, it,     iz),
                    vid(ir + 1, it + 1, iz),
                    vid(ir,     it + 1, iz),
                    vid(ir,     it,     iz + 1),
                    vid(ir + 1, it,     iz + 1),
                    vid(ir + 1, it + 1, iz + 1),
                    vid(ir,     it + 1, iz + 1),
                ))
    return coords, kvert, knpr


def write_tri(path, coords, kvert, knpr, comment):
    nvt = len(coords)
    nel = len(kvert)
    with open(path, "w") as f:
        f.write(f"{comment}\n")
        f.write("Parametrisierung PARXC, PARYC, TMAXC\n")
        f.write(f"{nel} {nvt} 1 8 12 6     NEL NVT NBCT NVE NEE NAE\n")
        f.write("DCORVG\n")
        for x, y, z in coords:
            f.write(f"{x:.16e} {y:.16e} {z:.16e}\n")
        f.write("KVERT\n")
        for vs in kvert:
            f.write(" ".join(str(v) for v in vs) + "\n")
        f.write("KNPR\n")
        for k in knpr:
            f.write(f"{k}\n")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ri", type=float, required=True, help="inner radius")
    ap.add_argument("--ro", type=float, required=True, help="outer radius")
    ap.add_argument("--lz", type=float, required=True, help="extrusion length")
    ap.add_argument("--nr", type=int, default=2, help="radial cells")
    ap.add_argument("--nt", type=int, default=12, help="azimuthal cells")
    ap.add_argument("--nz", type=int, default=2, help="axial cells")
    ap.add_argument("--out", required=True, help="output .tri path")
    a = ap.parse_args(argv)

    coords, kvert, knpr = generate_annulus(a.ri, a.ro, a.lz, a.nr, a.nt, a.nz)
    comment = (f"chimera_meshgen annulus ri={a.ri} ro={a.ro} lz={a.lz} "
               f"nr={a.nr} nt={a.nt} nz={a.nz}")
    write_tri(a.out, coords, kvert, knpr, comment)
    print(f"wrote {a.out}: nel={len(kvert)} nvt={len(coords)}")


if __name__ == "__main__":
    main(sys.argv[1:])
