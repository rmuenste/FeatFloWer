#!/usr/bin/env python3
"""Convert the per-rank EL trajectory logs of a run into ParaView particle files.

usage: el_particles_to_vtp.py <rundir> [--radius 1.5e-3] [--out _vtk/particles]

Reads every el_particles_rank<NNNN>.log (written with SimPar@ELParticleLogFreq > 0; columns
time ip fx fy fz tx ty tz px py pz vx vy vz, locally owned particles only), merges the ranks per
time, and writes one PolyData file <out>_<k>.vtp per logged time (points = particle centres,
point data: id, radius, velocity, force) plus <out>.pvd, a time series ParaView opens directly.
Glyph the points with spheres scaled by 'radius' (scale factor 1, scale mode off or 'radius')
to see the particles at true size. ASCII XML, no third-party modules."""
import argparse, glob, os, sys
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("rundir")
ap.add_argument("--radius", type=float, default=1.5e-3, help="sphere radius written per point (m)")
ap.add_argument("--out", default="_vtk/particles", help="output prefix relative to rundir")
a = ap.parse_args()

logs = sorted(glob.glob(os.path.join(a.rundir, "el_particles_rank*.log")))
if not logs:
    sys.exit("no el_particles_rank*.log in " + a.rundir)
frames = defaultdict(dict)   # time -> ip -> (p, v, f)
for p in logs:
    with open(p) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            c = line.split()
            t = float(c[0]); ip = int(c[1])
            vals = [float(x) for x in c[2:14]]
            frames[t][ip] = (vals[6:9], vals[9:12], vals[0:3])
times = sorted(frames)
out = os.path.join(a.rundir, a.out)
os.makedirs(os.path.dirname(out), exist_ok=True)
pvd = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">', '  <Collection>']
for k, t in enumerate(times):
    ids = sorted(frames[t])
    n = len(ids)
    fn = f"{out}_{k:04d}.vtp"
    with open(fn, "w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <PolyData>\n    <Piece NumberOfPoints="{n}" NumberOfVerts="{n}">\n')
        f.write('      <Points>\n        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        for ip in ids:
            f.write("          %.9e %.9e %.9e\n" % tuple(frames[t][ip][0]))
        f.write('        </DataArray>\n      </Points>\n')
        f.write('      <Verts>\n        <DataArray type="Int64" Name="connectivity" format="ascii">\n          ' + " ".join(str(i) for i in range(n)) + '\n        </DataArray>\n')
        f.write('        <DataArray type="Int64" Name="offsets" format="ascii">\n          ' + " ".join(str(i + 1) for i in range(n)) + '\n        </DataArray>\n      </Verts>\n')
        f.write('      <PointData Scalars="radius" Vectors="velocity">\n')
        f.write('        <DataArray type="Int64" Name="id" format="ascii">\n          ' + " ".join(str(ip) for ip in ids) + '\n        </DataArray>\n')
        f.write('        <DataArray type="Float64" Name="radius" format="ascii">\n          ' + " ".join("%.6e" % a.radius for _ in ids) + '\n        </DataArray>\n')
        for name, j in (("velocity", 1), ("force", 2)):
            f.write(f'        <DataArray type="Float64" Name="{name}" NumberOfComponents="3" format="ascii">\n')
            for ip in ids:
                f.write("          %.9e %.9e %.9e\n" % tuple(frames[t][ip][j]))
            f.write('        </DataArray>\n')
        f.write('      </PointData>\n    </Piece>\n  </PolyData>\n</VTKFile>\n')
    pvd.append(f'    <DataSet timestep="{t:.6f}" group="" part="0" file="{os.path.basename(fn)}"/>')
pvd += ['  </Collection>', '</VTKFile>']
open(out + ".pvd", "w").write("\n".join(pvd) + "\n")
counts = sorted({len(frames[t]) for t in times})
print(f"{len(times)} frames, t = {times[0]:.4f} .. {times[-1]:.4f}, particles per frame {counts[0]}..{counts[-1]}, written {out}_NNNN.vtp + {out}.pvd")
