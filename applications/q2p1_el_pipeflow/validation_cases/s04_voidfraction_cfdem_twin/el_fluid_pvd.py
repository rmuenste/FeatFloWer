#!/usr/bin/env python3
"""Write a ParaView time series (_vtk/main.pvd) for the fluid frames of a run.

usage: el_fluid_pvd.py <rundir> --dt 0.001 [--first-time 0.001] [--out _vtk/main.pvd]

FeatFloWer writes _vtk/main.<k>.pvtu per output frame without a time stamp; the index is the
step counter minus one (postprocessing_app: iXgmv = istep_ns - 1), so frame k is the state after
step k+1, t = first_time + k dt, and a skipped output step leaves a gap in the numbering. With
OutputFreq slightly below TimeStep every step writes one frame (s04 RUNBOOK section 13: the
framework's output clock is single precision and skips a step when OutputFreq equals TimeStep).
The script checks that the numbering is gap-free and writes main.pvd with those times, so
ParaView's animation clock matches the particle series written by el_particles_to_vtp.py."""
import argparse, glob, os, re, sys

ap = argparse.ArgumentParser()
ap.add_argument("rundir")
ap.add_argument("--dt", type=float, required=True, help="time step (s), one frame per step")
ap.add_argument("--first-time", type=float, default=None, help="time of frame 0 (default dt)")
ap.add_argument("--out", default="_vtk/main.pvd")
a = ap.parse_args()
t0 = a.dt if a.first_time is None else a.first_time

files = glob.glob(os.path.join(a.rundir, "_vtk", "main.*.pvtu"))
idx = sorted(int(re.search(r"main\.(\d+)\.pvtu$", f).group(1)) for f in files)
if not idx:
    sys.exit("no _vtk/main.*.pvtu in " + a.rundir)
gaps = [k for k in range(idx[0], idx[-1] + 1) if k not in set(idx)]
if idx[0] != 0 or gaps:
    sys.exit(f"frame numbering is not gap-free (first {idx[0]}, missing {gaps[:5]}): the times cannot be assigned by index")
width = len(re.search(r"main\.(\d+)\.pvtu$", files[0]).group(1))
lines = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">', '  <Collection>']
for k in idx:
    lines.append(f'    <DataSet timestep="{t0 + k * a.dt:.6f}" group="" part="0" file="main.{k:0{width}d}.pvtu"/>')
lines += ['  </Collection>', '</VTKFile>']
out = os.path.join(a.rundir, a.out)
open(out, "w").write("\n".join(lines) + "\n")
print(f"{len(idx)} frames, t = {t0:.4f} .. {t0 + idx[-1] * a.dt:.4f}, written {out}")
