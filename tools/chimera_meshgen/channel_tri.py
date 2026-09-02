#!/usr/bin/env python3
"""Generate a structured box (channel) background mesh for Chimera cases.

Writes a FeatFloWer coarse .tri (see annulus_tri.py for the format), the
boundary .par files and a .prj project file, ready for tools/PyPartitioner.py.
The box [0,Lx] x [0,Ly] x [0,Lz] is meshed with nx x ny x nz hexes in FEAT
vertex ordering (positive Jacobian).  Boundary types follow the DFG 2D
flow-around-cylinder deck (_adc/2D_FAC): x=0 Inflow2 (parabolic profile of
GetVeloBCVal), x=Lx Symmetry011 (do-nothing in u), y=0 / y=Ly Wall,
z=0 / z=Lz Symmetry001 (the 2D-in-3D slab symmetry).  There is NO body in
the mesh - the cylinder is handled by the Chimera component.

Usage (the vendored applications/q2p1_chimera/_data/CHIMERA_FAC case):
  channel_tri.py --lx 2.2 --ly 0.41 --lz 0.05 --nx 44 --ny 8 --nz 1 \
                 --name chimera_fac --outdir <dir>
"""

import argparse
import os
import sys

from annulus_tri import write_tri


def generate_box(lx, ly, lz, nx, ny, nz):
    nvt = (nx + 1) * (ny + 1) * (nz + 1)

    def vid(ix, iy, iz):
        return 1 + ix + (nx + 1) * (iy + (ny + 1) * iz)

    coords = [None] * nvt
    knpr = [0] * nvt
    for iz in range(nz + 1):
        for iy in range(ny + 1):
            for ix in range(nx + 1):
                i = vid(ix, iy, iz) - 1
                coords[i] = (lx * ix / nx, ly * iy / ny, lz * iz / nz)
                if (ix in (0, nx)) or (iy in (0, ny)) or (iz in (0, nz)):
                    knpr[i] = 1
    kvert = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                kvert.append((
                    vid(ix, iy, iz), vid(ix + 1, iy, iz),
                    vid(ix + 1, iy + 1, iz), vid(ix, iy + 1, iz),
                    vid(ix, iy, iz + 1), vid(ix + 1, iy, iz + 1),
                    vid(ix + 1, iy + 1, iz + 1), vid(ix, iy + 1, iz + 1),
                ))
    return coords, kvert, knpr


def write_par(path, btype, verts):
    with open(path, "w") as f:
        f.write(f"{len(verts)} {btype}\n")
        f.write("' '\n")
        for v in verts:
            f.write(f"{v}\n")


def main(argv):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--lx", type=float, required=True)
    ap.add_argument("--ly", type=float, required=True)
    ap.add_argument("--lz", type=float, required=True)
    ap.add_argument("--nx", type=int, required=True)
    ap.add_argument("--ny", type=int, required=True)
    ap.add_argument("--nz", type=int, required=True)
    ap.add_argument("--name", required=True, help="base name of the .tri/.prj")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args(argv)

    coords, kvert, knpr = generate_box(a.lx, a.ly, a.lz, a.nx, a.ny, a.nz)
    os.makedirs(a.outdir, exist_ok=True)
    comment = (f"chimera_meshgen channel lx={a.lx} ly={a.ly} lz={a.lz} "
               f"nx={a.nx} ny={a.ny} nz={a.nz}")
    write_tri(os.path.join(a.outdir, a.name + ".tri"), coords, kvert, knpr, comment)

    tol = 1e-12
    sel = lambda pred: [i + 1 for i, c in enumerate(coords) if pred(c)]
    pars = [
        ("in",     "Inflow2",     sel(lambda c: abs(c[0]) < tol)),
        ("out",    "Symmetry011", sel(lambda c: abs(c[0] - a.lx) < tol)),
        ("bottom", "Wall",        sel(lambda c: abs(c[1]) < tol)),
        ("top",    "Wall",        sel(lambda c: abs(c[1] - a.ly) < tol)),
        ("wall1",  "Symmetry001", sel(lambda c: abs(c[2]) < tol)),
        ("wall2",  "Symmetry001", sel(lambda c: abs(c[2] - a.lz) < tol)),
    ]
    for name, btype, verts in pars:
        write_par(os.path.join(a.outdir, name + ".par"), btype, verts)
    with open(os.path.join(a.outdir, a.name + ".prj"), "w") as f:
        f.write(a.name + ".tri\n")
        for name, _, _ in pars:
            f.write(name + ".par\n")
    print(f"wrote {a.outdir}/{a.name}.tri: nel={len(kvert)} nvt={len(coords)}")


if __name__ == "__main__":
    main(sys.argv[1:])
