#!/usr/bin/env pvbatch
"""Render a DNS drag-setup scene: velocity-magnitude mid-plane slice + spheres.

Headless (osmesa) ParaView script. The solver writes no indicator field, so
spheres are drawn geometrically from the run config (single centered sphere
for the Hasimoto cell, particles.xyz incl. periodic images for the arrays).

Usage:
  pvbatch tools/render_setup_scene.py <pvtu> <out.png> \
      --sphere cx,cy,cz,r            (repeatable)  OR
      --xyz particles.xyz --radius R
      [--slice-y 0.5] [--vmax auto] [--title "..."]
"""
import argparse
import sys

from paraview.simple import (
    XMLPartitionedUnstructuredGridReader, Slice, Sphere, Outline, Show,
    Render, SaveScreenshot, GetActiveViewOrCreate, ColorBy,
    GetColorTransferFunction, GetScalarBar, GetDisplayProperties,
)

p = argparse.ArgumentParser()
p.add_argument("pvtu")
p.add_argument("out")
p.add_argument("--sphere", action="append", default=[])
p.add_argument("--xyz")
p.add_argument("--radius", type=float)
p.add_argument("--slice-y", type=float, default=0.5)
p.add_argument("--title", default="velocity magnitude")
p.add_argument("--front-cut", action="store_true",
               help="omit spheres entirely in front of the slice plane "
                    "(center y > slice_y + r) so the slice stays visible")
args = p.parse_args()

spheres = [tuple(float(x) for x in s.split(",")) for s in args.sphere]
if args.xyz:
    with open(args.xyz) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 3:
                spheres.append((float(parts[0]), float(parts[1]),
                                float(parts[2]), args.radius))

reader = XMLPartitionedUnstructuredGridReader(FileName=[args.pvtu])
reader.PointArrayStatus = ["Velocity"]
reader.UpdatePipeline()

view = GetActiveViewOrCreate("RenderView")
view.ViewSize = [1800, 1500]
view.Background = [1.0, 1.0, 1.0]
view.OrientationAxesVisibility = 0
view.UseColorPaletteForBackground = 0

sl = Slice(Input=reader)
sl.SliceType.Origin = [0.5, args.slice_y, 0.5]
sl.SliceType.Normal = [0.0, 1.0, 0.0]
disp = Show(sl, view)
ColorBy(disp, ("POINTS", "Velocity", "Magnitude"))
disp.SetRepresentationType("Surface")
lut = GetColorTransferFunction("Velocity")
lut.ApplyPreset("Viridis (matplotlib)", True)
disp.RescaleTransferFunctionToDataRange(True, False)
bar = GetScalarBar(lut, view)
bar.Title = args.title
bar.ComponentTitle = ""
bar.TitleColor = [0.1, 0.1, 0.1]
bar.LabelColor = [0.1, 0.1, 0.1]
bar.ScalarBarLength = 0.30
disp.SetScalarBarVisibility(view, True)

out = Outline(Input=reader)
odisp = Show(out, view)
odisp.AmbientColor = [0.25, 0.25, 0.25]
odisp.DiffuseColor = [0.25, 0.25, 0.25]
odisp.LineWidth = 1.5

if args.front_cut:
    spheres = [s for s in spheres if s[1] <= args.slice_y + s[3]]

for (cx, cy, cz, r) in spheres:
    sp = Sphere()
    sp.Center = [cx, cy, cz]
    sp.Radius = r
    sp.ThetaResolution = 48
    sp.PhiResolution = 48
    sdisp = Show(sp, view)
    sdisp.DiffuseColor = [0.62, 0.64, 0.68]
    sdisp.Specular = 0.35
    sdisp.SpecularPower = 40.0

view.CameraPosition = [2.35, -1.75, 1.85]
view.CameraFocalPoint = [0.5, 0.5, 0.5]
view.CameraViewUp = [0.0, 0.0, 1.0]
view.CameraViewAngle = 28.0

Render()
SaveScreenshot(args.out, view, ImageResolution=[1800, 1500],
               TransparentBackground=0)
print("wrote", args.out, "with", len(spheres), "spheres")
