"""Export Blender resources for the D6.1 Oberbeck spheroid (fixed prolate spheroid, body-force-driven flow in a periodic unit cell).

Usage: pvpython --mesa --force-offscreen-rendering export_blender_oberbeck.py PVTU AX AY AZ OUTDIR TAG
Cell: [0,1]^3 periodic, body at (0.5, 0.5, 0.5), semi-axes a = 0.264567 (along the axis), b = c = 0.132283;
the driving force and the mean flow are along +z. Slice plane y = 0.5 contains the axis for V1 (+z) and V2 (+x).

Produces in OUTDIR:
  TAG_plane_{color,overlay,combined}.png   x,z in [0,1], 4000 x 4000 px; colour |u| / u_max of the run, streamlines + arrows of u
  TAG_plane_streamlines.png                streamlines only (transparent)
  colorbar_oberbeck.png
"""
import os
import sys
from paraview.simple import *

pvtu = sys.argv[1]
ex, ey, ez = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
outdir, tag = sys.argv[5], sys.argv[6]
os.makedirs(outdir, exist_ok=True)
A, B = 0.26456684, 0.13228342
C = (0.5, 0.5, 0.5)
PRESET = "/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tu-colormaps.json"
PX = 4000

rd = XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
rd.PointArrayStatus = ["Velocity"]
rd.UpdatePipeline()
umax = rd.PointData["Velocity"].GetRange(-1)[1]
print("cells", rd.GetDataInformation().GetNumberOfCells(), "umax", umax)

sl = Slice(Input=rd)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.5, 0.5, 0.5]
sl.SliceType.Normal = [0.0, 1.0, 0.0]
c1 = Calculator(Input=sl, ResultArrayName="unorm", Function="mag(Velocity)/%g" % umax)
dot = "((coordsX-%g)*%g+(coordsY-%g)*%g+(coordsZ-%g)*%g)" % (C[0], ex, C[1], ey, C[2], ez)
rr = "((coordsX-%g)^2+(coordsY-%g)^2+(coordsZ-%g)^2)" % C
c2 = Calculator(Input=c1, ResultArrayName="q", Function="%s*%s/%g + (%s-%s*%s)/%g" % (dot, dot, A * A, rr, dot, dot, B * B))
body = Clip(Input=c2, ClipType="Scalar")
body.Scalars = ["POINTS", "q"]
body.Value = 0.85
body.Invert = 0
body.UpdatePipeline()
print("outside-body unorm range", body.PointData["unorm"].GetRange())

iso = Contour(Input=body, ContourBy=["POINTS", "unorm"], Isosurfaces=[0.2, 0.4, 0.6, 0.8, 0.9, 0.95])

# streamlines of the in-plane velocity, seeded on a line across the bottom of the cell
st = StreamTracer(Input=c2, SeedType="Line")
st.Vectors = ["POINTS", "Velocity"]
st.SeedType.Point1 = [0.0125, 0.5, 0.01]
st.SeedType.Point2 = [0.9875, 0.5, 0.01]
st.SeedType.Resolution = 39   # seeds at 0.0125 + k*0.025: none on the body axis x = 0.5
st.MaximumStreamlineLength = 3.0
st.IntegrationDirection = "BOTH"
st.SurfaceStreamlines = 1

# direction arrows on a grid
n = 36
g = ResampleToImage(Input=c2)
g.UseInputBounds = 0
g.SamplingDimensions = [n, 1, n]
dd = 1.0 / n
g.SamplingBounds = [dd / 2, 1 - dd / 2, 0.5, 0.5, dd / 2, 1 - dd / 2]
m1 = Threshold(Input=g, Scalars=["POINTS", "q"], LowerThreshold=1.2, UpperThreshold=1e9, ThresholdMethod="Between")
gl = Glyph(Input=m1, GlyphType="Arrow")
gl.GlyphType.TipResolution = 12; gl.GlyphType.TipRadius = 0.16; gl.GlyphType.TipLength = 0.4; gl.GlyphType.ShaftRadius = 0.05
gl.OrientationArray = ["POINTS", "Velocity"]; gl.ScaleArray = ["POINTS", "No scale array"]
gl.ScaleFactor = 0.6 * dd; gl.GlyphMode = "All Points"

v = GetActiveViewOrCreate("RenderView")
v.ViewSize = [PX, PX]
v.OrientationAxesVisibility = 0
try:
    v.UseColorPaletteForBackground = 0
except Exception:
    pass
v.Background = [1.0, 1.0, 1.0]
ImportPresets(filename=PRESET)
lut = GetColorTransferFunction("unorm")
lut.AutomaticRescaleRangeMode = "Never"
lut.ApplyPreset("Viridis (matplotlib)", True)
lut.Discretize = 1; lut.NumberOfTableValues = 20
dfl = Show(body, v); ColorBy(dfl, ("POINTS", "unorm")); dfl.Ambient, dfl.Diffuse = 1.0, 0.0
lut.RescaleTransferFunction(0.0, 1.0); GetOpacityTransferFunction("unorm").RescaleTransferFunction(0.0, 1.0)
diso = Show(iso, v); ColorBy(diso, None); diso.AmbientColor = diso.DiffuseColor = [0.05, 0.05, 0.05]; diso.Ambient, diso.Diffuse = 1.0, 0.0; diso.LineWidth = 2.5
dst = Show(st, v); ColorBy(dst, None); dst.AmbientColor = dst.DiffuseColor = [0.95, 0.95, 0.95]; dst.Ambient, dst.Diffuse = 1.0, 0.0; dst.LineWidth = 3.0
dgl = Show(gl, v); ColorBy(dgl, None); dgl.AmbientColor = dgl.DiffuseColor = [0.0, 0.0, 0.0]; dgl.Ambient, dgl.Diffuse = 1.0, 0.0
Render(v)
v.InteractionMode = "2D"; v.CameraParallelProjection = 1
v.CameraPosition = [0.5, -10.0, 0.5]; v.CameraFocalPoint = [0.5, 0.5, 0.5]; v.CameraViewUp = [0, 0, 1]; v.CameraParallelScale = 0.5
for name, fl, il, sll, gll in (("combined", 1, 1, 1, 1), ("color", 1, 0, 0, 0), ("overlay", 0, 1, 0, 1), ("streamlines", 0, 0, 1, 0)):
    dfl.Visibility, diso.Visibility, dst.Visibility, dgl.Visibility = fl, il, sll, gll
    Render(v)
    SaveScreenshot(os.path.join(outdir, "%s_plane_%s.png" % (tag, name)), v, ImageResolution=[PX, PX], TransparentBackground=1)
    print("wrote", "%s_plane_%s.png" % (tag, name))

dfl.Visibility = 1; dfl.SetScalarBarVisibility(v, True)
bar = GetScalarBar(lut, v); bar.Title = "|u| / u_max"; bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]; bar.TitleFontSize, bar.LabelFontSize = 36, 30
bar.WindowLocation = "Any Location"; bar.Orientation = "Vertical"; bar.Position = [0.25, 0.08]
bar.ScalarBarLength = 0.84; bar.ScalarBarThickness = 40; bar.UseCustomLabels = 1; bar.CustomLabels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]; bar.AddRangeLabels = 0
v.ViewSize = [500, 1600]; v.CameraFocalPoint = [50.0, 0.5, 50.0]; v.CameraPosition = [50.0, -10.0, 50.0]
Render(v)
SaveScreenshot(os.path.join(outdir, "colorbar_oberbeck.png"), v, ImageResolution=[500, 1600], TransparentBackground=1)
print("wrote colorbar; umax", umax)
