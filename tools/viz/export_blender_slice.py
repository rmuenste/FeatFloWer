"""Export Blender resources for the ten Cate sedimentation still.

Usage:
  pvpython --mesa --force-offscreen-rendering export_blender_slice.py RUNDIR FRAME U_INF PZ OUTDIR

Produces, in OUTDIR:
  slice_y0_color.png     colour field |u|/u_inf on the plane y = 0, full column width, transparent outside
  slice_y0_overlay.png   isolines + direction arrows only, transparent background (same pixel frame)
  slice_y0_combined.png  colour + isolines + arrows
  colorbar.png           the |u|/u_inf legend, 0 to 1, transparent background
  isolines.obj           isolines as 3D polylines (metres, y = 0)
  arrows.ply             direction arrows as triangle meshes (metres, y = 0)

Pixel frame of the three slice images: x from -0.05 to +0.05 m left to right, z from 0 (bottom) to 0.16 m (top),
looking from -y toward +y. 2000 x 3200 px, i.e. 50 micrometres per pixel.
"""
import os
import sys
from paraview.simple import *

rundir, frame, u_inf, pz, outdir = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
# optional: colour map preset name, number of filled bands (0 = continuous), file suffix
CMAP = sys.argv[6] if len(sys.argv) > 6 else "TU Petrol Sequential"
BANDS = int(sys.argv[7]) if len(sys.argv) > 7 else 0
SUFFIX = sys.argv[8] if len(sys.argv) > 8 else ""
os.makedirs(outdir, exist_ok=True)

DP, R = 0.015, 0.0075
YS = -2.0e-5                       # slice plane, a hair inside the quarter domain (y <= 0)
XW, ZH = 0.10, 0.16                # full column width and height [m]
W, H = 2000, 3200                  # pixels; W/H == XW/ZH
LEVELS = [0.04, 0.08, 0.12, 0.16, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
if BANDS:
    LEVELS = [i / BANDS for i in range(1, BANDS)]   # isolines on the band edges
NX, NZ = 34, 54                    # arrow grid
PRESET = "/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tu-colormaps.json"

# ---- data: quarter domain -> full width by reflection across x = 0 ---------
rd = XMLPartitionedUnstructuredGridReader(FileName=[os.path.join(rundir, "_vtk", "main.%s.pvtu" % frame)])
rd.PointArrayStatus = ["Velocity", "Mixer"]
calc = Calculator(Input=rd, ResultArrayName="unorm", Function="mag(Velocity)/%g" % u_inf)
full = Reflect(Input=calc)
full.PlaneMode = "X Max"           # the quarter domain ends at x = 0
full.CopyInput = 1
full.ReflectAllInputArrays = 1     # flips the x component of Velocity in the mirrored half

sl = Slice(Input=full)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, YS, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]

# exact circular hole where the sphere sits (the Blender sphere fills it)
hole = Clip(Input=sl, ClipType="Sphere")
hole.ClipType.Center = [0.0, 0.0, pz]
hole.ClipType.Radius = R
hole.Invert = 0                    # keep outside

iso = Contour(Input=hole, ContourBy=["POINTS", "unorm"], Isosurfaces=LEVELS)

grid = ResampleToImage(Input=full)
grid.UseInputBounds = 0
grid.SamplingDimensions = [NX, 1, NZ]
dx, dz = XW / NX, ZH / NZ
grid.SamplingBounds = [-XW / 2 + dx / 2, XW / 2 - dx / 2, YS, YS, dz / 2, ZH - dz / 2]
gm = Threshold(Input=grid, Scalars=["POINTS", "Mixer"], LowerThreshold=-1.0, UpperThreshold=0.5, ThresholdMethod="Between")
gm2 = Threshold(Input=gm, Scalars=["POINTS", "unorm"], LowerThreshold=0.01, UpperThreshold=1e9, ThresholdMethod="Between")
gl = Glyph(Input=gm2, GlyphType="Arrow")
gl.GlyphType.TipResolution = 12
gl.GlyphType.TipRadius = 0.16
gl.GlyphType.TipLength = 0.4
gl.GlyphType.ShaftRadius = 0.05
gl.OrientationArray = ["POINTS", "Velocity"]
gl.ScaleArray = ["POINTS", "No scale array"]
gl.ScaleFactor = 0.62 * dx
gl.GlyphMode = "All Points"

# ---- geometry exports -----------------------------------------------------
isoMB = MergeBlocks(Input=iso)
isoM = ExtractSurface(Input=isoMB)      # OBJ/PLY writers need polydata
glMB = MergeBlocks(Input=gl)
glM = ExtractSurface(Input=glMB)
SaveData(os.path.join(outdir, "isolines%s.obj" % SUFFIX), proxy=isoM)
SaveData(os.path.join(outdir, "arrows.ply"), proxy=glM)

# ---- view ----------------------------------------------------------------
v = GetActiveViewOrCreate("RenderView")
v.ViewSize = [W, H]
v.OrientationAxesVisibility = 0
try:
    v.UseColorPaletteForBackground = 0
except Exception:
    pass
v.Background = [1.0, 1.0, 1.0]

lut = GetColorTransferFunction("unorm")
lut.AutomaticRescaleRangeMode = "Never"
ImportPresets(filename=PRESET)
lut.ApplyPreset(CMAP, True)
if BANDS:
    lut.Discretize = 1
    lut.NumberOfTableValues = BANDS
else:
    lut.Discretize = 0

dfl = Show(hole, v)
ColorBy(dfl, ("POINTS", "unorm"))
lut.RescaleTransferFunction(0.0, 1.0)
GetOpacityTransferFunction("unorm").RescaleTransferFunction(0.0, 1.0)
dfl.Ambient, dfl.Diffuse = 1.0, 0.0          # exact LUT colours, no shading

diso = Show(iso, v)
ColorBy(diso, None)
diso.AmbientColor = diso.DiffuseColor = [0.12, 0.12, 0.12]
diso.Ambient, diso.Diffuse = 1.0, 0.0
diso.LineWidth = 2.5

dgl = Show(gl, v)
ColorBy(dgl, None)
dgl.AmbientColor = dgl.DiffuseColor = [0.0, 0.0, 0.0]
dgl.Ambient, dgl.Diffuse = 1.0, 0.0

Render(v)
v.InteractionMode = "2D"
v.CameraParallelProjection = 1
v.CameraPosition = [0.0, -1.0, ZH / 2]       # look from -y toward +y: +x to the right, +z up
v.CameraFocalPoint = [0.0, 0.0, ZH / 2]
v.CameraViewUp = [0.0, 0.0, 1.0]
v.CameraParallelScale = ZH / 2
Render(v)


def shot(name, fl, il, gll):
    dfl.Visibility, diso.Visibility, dgl.Visibility = fl, il, gll
    Render(v)
    SaveScreenshot(os.path.join(outdir, name), v, ImageResolution=[W, H], TransparentBackground=1)
    print("wrote", name)


shot("slice_y0_combined%s.png" % SUFFIX, 1, 1, 1)
shot("slice_y0_color%s.png" % SUFFIX, 1, 0, 0)
shot("slice_y0_overlay%s.png" % SUFFIX, 0, 1, 1)

# ---- colour bar alone -------------------------------------------------------
dfl.Visibility = 1
dfl.SetScalarBarVisibility(v, True)
bar = GetScalarBar(lut, v)
bar.Title = "|u| / u_inf"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]
bar.TitleFontSize, bar.LabelFontSize = 36, 30
bar.WindowLocation = "Any Location"
bar.Orientation = "Vertical"
bar.Position = [0.25, 0.08]
bar.ScalarBarLength = 0.84
bar.ScalarBarThickness = 40
bar.UseCustomLabels = 1
bar.CustomLabels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
bar.AddRangeLabels = 0
v.ViewSize = [500, 1600]
v.CameraFocalPoint = [10.0, 0.0, 10.0]     # look at empty space: only the legend renders
v.CameraPosition = [10.0, -1.0, 10.0]
Render(v)
SaveScreenshot(os.path.join(outdir, "colorbar%s.png" % SUFFIX), v, ImageResolution=[500, 1600], TransparentBackground=1)
print("wrote colorbar.png")
