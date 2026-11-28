#!/usr/bin/env python3
"""Coarse hexahedral box mesh of case 11 (minimum fluidisation): laterally periodic column,
uniform inflow through the floor, open outlet on top.

Layout of the peer's make_box_mesh.py (case_specs/04_voidfraction/featflower; hex_ex.py export:
x-fastest node numbering, KVERT bottom face then top face, KNPR 0, one .par per face with the
1-based node list), boundary types of the case:
  xmin xmax ymin ymax : Periodic   (SimPar@PeriodicAxis = xy, PeriodicLength = LX = LY)
  zmin                : Inflow975  (uniform inlet, EL_CONFIG / EL_INLET_TABLE)
  zmax                : Outflow    (do-nothing outlet, SimPar@NoOutflow = No)
The parameter line of the periodic, inflow and outflow faces is '' (as in the s10 periodic box
and the s05 pipe).

usage: make_s11_mesh.py OUTDIR [--size LX LY LZ] [--cells NX NY NZ]
default: bed R, 0.015 x 0.015 x 0.10 m with 3 x 3 x 20 coarse cells (5 mm); MaxMeshLevel 2
gives 6 x 6 x 40 cells of 2.5 mm (2.9 d for d 0.87 mm), the mesh of the CFDEM run of the
peer (run_fluid.sh FL_NXY 6, FL_NZ 40)."""
import argparse
import os

TYPES = {'xmin': 'Periodic', 'xmax': 'Periodic', 'ymin': 'Periodic', 'ymax': 'Periodic',
         'zmin': 'Inflow975', 'zmax': 'Outflow'}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('outdir')
    ap.add_argument('--size', nargs=3, type=float, default=[0.015, 0.015, 0.10])
    ap.add_argument('--cells', nargs=3, type=int, default=[3, 3, 20])
    ap.add_argument('--name', default='box.tri')
    a = ap.parse_args()
    lx, ly, lz = a.size
    nx, ny, nz = a.cells
    os.makedirs(a.outdir, exist_ok=True)

    def nid(i, j, k):
        return 1 + i + (nx + 1) * (j + (ny + 1) * k)

    nvt = (nx + 1) * (ny + 1) * (nz + 1)
    nel = nx * ny * nz
    with open(os.path.join(a.outdir, a.name), 'w') as f:
        f.write('Coarse mesh exported by make_s11_mesh.py (hex_ex.py layout) \n')
        f.write('Parametrisierung PARXC, PARYC, TMAXC \n')
        f.write(f'{nel} {nvt} 0 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE \n')
        f.write('DCORVG\n')
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    f.write(f'{float(lx * i / nx)!r} {float(ly * j / ny)!r} {float(lz * k / nz)!r}\n')
        f.write('KVERT\n')
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    v = [nid(i, j, k), nid(i + 1, j, k), nid(i + 1, j + 1, k), nid(i, j + 1, k),
                         nid(i, j, k + 1), nid(i + 1, j, k + 1), nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1)]
                    f.write(' '.join(map(str, v)) + '\n')
        f.write('KNPR\n')
        for _ in range(nvt):
            f.write('0\n')

    sel = {'xmin': lambda i, j, k: i == 0, 'xmax': lambda i, j, k: i == nx,
           'ymin': lambda i, j, k: j == 0, 'ymax': lambda i, j, k: j == ny,
           'zmin': lambda i, j, k: k == 0, 'zmax': lambda i, j, k: k == nz}
    for name, s in sel.items():
        ids = [nid(i, j, k) for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1) if s(i, j, k)]
        with open(os.path.join(a.outdir, name + '.par'), 'w') as f:
            f.write(f'{len(ids)} {TYPES[name]}\n')
            f.write("''\n")
            for v in ids:
                f.write(f'{v}\n')
    with open(os.path.join(a.outdir, 'file.prj'), 'w') as f:
        f.write(a.name + '\n' + '\n'.join(n + '.par' for n in sel) + '\n')
    print(f'{a.outdir}: NEL={nel} NVT={nvt} coarse h=({lx/nx:g},{ly/ny:g},{lz/nz:g}) m; '
          f'faces ' + ', '.join(f'{n} {t}' for n, t in TYPES.items()))


if __name__ == '__main__':
    main()
