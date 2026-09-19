"""ten Cate Fig. 9 style still: half-plane view of the sedimenting sphere at h/d_p ~ 0.5.

Usage: pvpython --mesa --force-offscreen-rendering tencate_fig9.py RUNDIR FRAME U_INF PZ OUT.png "label"
  RUNDIR : FeatFloWer rundir with _vtk/main.FRAME.pvtu
  U_INF  : normalisation velocity [m/s]
  PZ     : sphere centre height [m] from particle_force.log at the frame time
"""
import sys, os
from paraview.simple import *

rundir, frame, u_inf, pz, out, label = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5], sys.argv[6]
DP = 0.015
R = DP / 2
YS = -2.0e-5                # slice plane, a hair inside the domain (y <= 0)
FRAME = 2.7 * DP            # frame edge length [m]
pvtu = os.path.join(rundir, "_vtk", "main.%s.pvtu" % frame)

# ---- data ----------------------------------------------------------------
rd = XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
rd.PointArrayStatus = ["Velocity", "Mixer"]
rd.UpdatePipeline()
di = rd.GetDataInformation()
print("points:", di.GetNumberOfPoints(), "cells:", di.GetNumberOfCells())

calc = Calculator(Input=rd)
calc.ResultArrayName = "unorm"
calc.Function = "mag(Velocity)/%g" % u_inf

sl = Slice(Input=calc)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, YS, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]
sl.UpdatePipeline()
print("slice unorm range", sl.PointData["unorm"].GetRange())
# keep only the frame region (x in [-FRAME, 0], z in [0, FRAME])
sl = Clip(Input=sl, ClipType="Plane")
sl.ClipType.Origin = [-FRAME, 0.0, 0.0]
sl.ClipType.Normal = [1.0, 0.0, 0.0]
sl.Invert = 0
sl = Clip(Input=sl, ClipType="Plane")
sl.ClipType.Origin = [0.0, 0.0, FRAME]
sl.ClipType.Normal = [0.0, 0.0, 1.0]
sl.Invert = 1

iso = Contour(Input=sl)
iso.ContourBy = ["POINTS", "unorm"]
iso.Isosurfaces = [0.04, 0.08, 0.12, 0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

# regular PIV-like grid of direction arrows inside the frame
n = 28
grid = ResampleToImage(Input=calc)
grid.UseInputBounds = 0
grid.SamplingDimensions = [n, 1, n]
grid.SamplingBounds = [-FRAME + FRAME / (2 * n), -FRAME / (2 * n), YS, YS, FRAME / (2 * n), FRAME - FRAME / (2 * n)]
gmask = Threshold(Input=grid)
gmask.Scalars = ["POINTS", "Mixer"]
gmask.LowerThreshold = -1.0
gmask.UpperThreshold = 0.5
gmask.ThresholdMethod = "Between"
gmask2 = Threshold(Input=gmask)
gmask2.Scalars = ["POINTS", "unorm"]
gmask2.LowerThreshold = 0.01
gmask2.UpperThreshold = 1e9
gmask2.ThresholdMethod = "Between"
gl = Glyph(Input=gmask2, GlyphType="Arrow")
gl.GlyphType.TipResolution = 12
gl.GlyphType.TipRadius = 0.16
gl.GlyphType.TipLength = 0.4
gl.GlyphType.ShaftRadius = 0.05
gl.OrientationArray = ["POINTS", "Velocity"]
gl.ScaleArray = ["POINTS", "No scale array"]
gl.ScaleFactor = 0.62 * FRAME / n
gl.GlyphMode = "All Points"

# exact sphere at the logged centre; its y>0 half sits in front of the slice
sph = Sphere(Center=[0.0, 0.0, pz], Radius=R, ThetaResolution=96, PhiResolution=96)
wall = Line(Point1=[-0.05, YS + 1e-5, 0.0], Point2=[0.0, YS + 1e-5, 0.0])

# ---- view ----------------------------------------------------------------
v = GetActiveViewOrCreate("RenderView")
W, H, FPX = 2000, 1750, 1600   # canvas px, field px
v.ViewSize = [W, H]
v.OrientationAxesVisibility = 0
try:
    v.UseColorPaletteForBackground = 0
except Exception:
    pass
v.Background = [1.0, 1.0, 1.0]

lut = GetColorTransferFunction("unorm")
lut.AutomaticRescaleRangeMode = "Never"
try:
    ImportPresets(filename="/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tu-colormaps.json")
    lut.ApplyPreset("TU Petrol Sequential", True)
except Exception as e:
    print("preset fallback:", e)
    lut.ApplyPreset("Viridis (matplotlib)", True)
lut.RescaleTransferFunction(0.0, 1.0)

dfl = Show(sl, v)
ColorBy(dfl, ("POINTS", "unorm"))
lut.RescaleTransferFunction(0.0, 1.0)
GetOpacityTransferFunction("unorm").RescaleTransferFunction(0.0, 1.0)
dfl.SetScalarBarVisibility(v, True)
bar = GetScalarBar(lut, v)
bar.Title = "|u| / u_inf"
bar.ComponentTitle = ""
bar.TitleColor = [0, 0, 0]
bar.LabelColor = [0, 0, 0]
bar.TitleFontSize = 30
bar.LabelFontSize = 26
bar.WindowLocation = "Any Location"
bar.Orientation = "Vertical"
bar.Position = [0.84, 0.12]
bar.ScalarBarLength = 0.7
bar.ScalarBarThickness = 28
bar.UseCustomLabels = 1
bar.CustomLabels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
bar.AddRangeLabels = 0

diso = Show(iso, v)
ColorBy(diso, None)
diso.AmbientColor = [0.12, 0.12, 0.12]
diso.DiffuseColor = [0.12, 0.12, 0.12]
diso.LineWidth = 2.0

dgl = Show(gl, v)
ColorBy(dgl, None)
dgl.AmbientColor = [0, 0, 0]
dgl.DiffuseColor = [0, 0, 0]
dgl.Ambient = 1.0
dgl.Diffuse = 0.0

dsp = Show(sph, v)
dsp.ColorArrayName = ["POINTS", ""]
dsp.DiffuseColor = [0.97, 0.97, 0.97]
dsp.Specular = 0.25
dsp.SpecularPower = 40

dw = Show(wall, v)
dw.AmbientColor = [0, 0, 0]
dw.DiffuseColor = [0, 0, 0]
dw.LineWidth = 6.0

txt = Text(Text=label)
dt = Show(txt, v)
dt.Color = [0, 0, 0]
dt.FontSize = 32
dt.WindowLocation = "Upper Left Corner"

Render(v)  # first render resets the camera; set ours afterwards
wpp = FRAME / FPX                      # world units per pixel
scale = 0.5 * H * wpp                  # vertical half extent
fx, fz = -0.5 * W * wpp, scale - 20 * wpp   # x=0 at left edge, wall 20 px above bottom
v.InteractionMode = "2D"
v.CameraParallelProjection = 1
v.CameraPosition = [fx, 1.0, fz]
v.CameraFocalPoint = [fx, 0.0, fz]
v.CameraViewUp = [0.0, 0.0, 1.0]
v.CameraParallelScale = scale
Render(v)
SaveScreenshot(out, v, ImageResolution=[W, H])
print("wrote", out)
