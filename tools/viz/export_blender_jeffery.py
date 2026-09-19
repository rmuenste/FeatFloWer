"""Export Blender resources for the D6.2 Jeffery orbit still (planar Couette box, prolate spheroid).

Usage:
  pvpython --mesa --force-offscreen-rendering export_blender_jeffery.py PVTU AX AY AZ OUTDIR TAG

Box units (rho = nu = 1): box x in [-4,4], y in [-3,3], z in [-4,4]; walls at z = +-4 move at +-U x^, U = 0.8,
gammadot = 0.2; spheroid a = 0.5 (along the axis), b = c = 0.25, centre (0,0,0), axis (AX,AY,AZ).

Produces in OUTDIR (TAG = e.g. t090):
  TAG_full_dist_{color,overlay,combined}.png   plane y = 0, x,z in [-4,4], disturbance field u' = u - gammadot z x^,
                                               colour |u'|/(gammadot a), isolines, direction arrows;  4000 x 4000 px
  TAG_zoom_dist_{color,overlay,combined}.png   same, window x,z in [-2,2];  4000 x 4000 px
  TAG_full_total_color.png                     colour |u|/U of the total velocity (the linear shear dominates)
  TAG_zoom_isolines.obj, TAG_zoom_arrows.ply   geometry of the zoom window (box units, y = 0)
  colorbar_dist.png                            legend |u'|/(gammadot a)
"""
import os
import sys
from paraview.simple import *

pvtu = sys.argv[1]
ex, ey, ez = float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
outdir, tag = sys.argv[5], sys.argv[6]
os.makedirs(outdir, exist_ok=True)
# optional: colour map preset, number of filled bands (0 = continuous), output suffix
CMAP = sys.argv[7] if len(sys.argv) > 7 else "TU Petrol Sequential"
BANDS = int(sys.argv[8]) if len(sys.argv) > 8 else 0
SUFFIX = sys.argv[9] if len(sys.argv) > 9 else ""
LINECOL = [float(x) for x in sys.argv[10].split(",")] if len(sys.argv) > 10 else None   # optional arrow/isoline colour r,g,b

A, B = 0.5, 0.25
GD, U = 0.2, 0.8
UREF = GD * A                      # velocity scale of the disturbance field: 0.1
CMAX = 0.4                         # colour range top for |u'|/(gammadot a); peaks 0.47 (t=90), 0.34 (t=100)
LEVELS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8]
PRESET = "/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tu-colormaps.json"
PX = 4000

rd = XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
rd.PointArrayStatus = ["Velocity"]
rd.UpdatePipeline()
print("cells", rd.GetDataInformation().GetNumberOfCells(), "bounds", rd.GetDataInformation().GetBounds())

sl = Slice(Input=rd)
sl.SliceType = "Plane"
sl.SliceType.Origin = [0.0, 0.0, 0.0]
sl.SliceType.Normal = [0.0, 1.0, 0.0]

# disturbance velocity and its magnitude; total magnitude; ellipsoid inside test q (< 1 inside the body)
c1 = Calculator(Input=sl, ResultArrayName="udist", Function="Velocity - %g*coordsZ*iHat" % GD)
c2 = Calculator(Input=c1, ResultArrayName="dnorm", Function="mag(udist)/%g" % UREF)
c3 = Calculator(Input=c2, ResultArrayName="tnorm", Function="mag(Velocity)/%g" % U)
dot = "(coordsX*%g+coordsY*%g+coordsZ*%g)" % (ex, ey, ez)
c4 = Calculator(Input=c3, ResultArrayName="q",
                Function="%s*%s/%g + (coordsX*coordsX+coordsY*coordsY+coordsZ*coordsZ-%s*%s)/%g" % (dot, dot, A * A, dot, dot, B * B))
c4.UpdatePipeline()
print("dnorm range", c4.PointData["dnorm"].GetRange(), "tnorm range", c4.PointData["tnorm"].GetRange())

body = Clip(Input=c4, ClipType="Scalar")
body.Scalars = ["POINTS", "q"]
body.Value = 0.85                  # hole 8% inside the true ellipse so the Blender body always covers it
body.Invert = 0                    # keep q > 1: outside the spheroid
body.UpdatePipeline()
print("outside-body dnorm range", body.PointData["dnorm"].GetRange())

iso = Contour(Input=body, ContourBy=["POINTS", "dnorm"], Isosurfaces=LEVELS)

# ---- view ----------------------------------------------------------------
v = GetActiveViewOrCreate("RenderView")
v.ViewSize = [PX, PX]
v.OrientationAxesVisibility = 0
try:
    v.UseColorPaletteForBackground = 0
except Exception:
    pass
v.Background = [1.0, 1.0, 1.0]

ImportPresets(filename=PRESET)
lutd = GetColorTransferFunction("dnorm")
lutd.AutomaticRescaleRangeMode = "Never"
lutd.ApplyPreset(CMAP, True)
if BANDS:
    lutd.Discretize = 1
    lutd.NumberOfTableValues = BANDS
    LEVELS = [CMAX * i / BANDS for i in range(1, BANDS)]   # isolines on the band edges
    iso.Isosurfaces = LEVELS
else:
    lutd.Discretize = 0
lutt = GetColorTransferFunction("tnorm")
lutt.AutomaticRescaleRangeMode = "Never"
lutt.ApplyPreset("TU Petrol-Orange Diverging", True)

dfl = Show(body, v)
ColorBy(dfl, ("POINTS", "dnorm"))
lutd.RescaleTransferFunction(0.0, CMAX)
GetOpacityTransferFunction("dnorm").RescaleTransferFunction(0.0, CMAX)
dfl.Ambient, dfl.Diffuse = 1.0, 0.0

diso = Show(iso, v)
ColorBy(diso, None)
diso.AmbientColor = diso.DiffuseColor = LINECOL or [0.12, 0.12, 0.12]
diso.Ambient, diso.Diffuse = 1.0, 0.0
diso.LineWidth = 3.0

Render(v)
v.InteractionMode = "2D"
v.CameraParallelProjection = 1
v.CameraViewUp = [0.0, 0.0, 1.0]


def arrows(half, n):
    """direction arrows of the disturbance field on an n x n grid over |x|,|z| <= half"""
    g = ResampleToImage(Input=c4)
    g.UseInputBounds = 0
    g.SamplingDimensions = [n, 1, n]
    d = 2 * half / n
    g.SamplingBounds = [-half + d / 2, half - d / 2, 0.0, 0.0, -half + d / 2, half - d / 2]
    m1 = Threshold(Input=g, Scalars=["POINTS", "q"], LowerThreshold=1.15, UpperThreshold=1e9, ThresholdMethod="Between")
    m2 = Threshold(Input=m1, Scalars=["POINTS", "dnorm"], LowerThreshold=0.01, UpperThreshold=1e9, ThresholdMethod="Between")
    gl = Glyph(Input=m2, GlyphType="Arrow")
    gl.GlyphType.TipResolution = 12
    gl.GlyphType.TipRadius = 0.16
    gl.GlyphType.TipLength = 0.4
    gl.GlyphType.ShaftRadius = 0.05
    gl.OrientationArray = ["POINTS", "udist"]
    gl.ScaleArray = ["POINTS", "No scale array"]
    gl.ScaleFactor = 0.62 * d
    gl.GlyphMode = "All Points"
    return gl


def shots(prefix, half, n):
    gl = arrows(half, n)
    dgl = Show(gl, v)
    ColorBy(dgl, None)
    dgl.AmbientColor = dgl.DiffuseColor = LINECOL or [0.0, 0.0, 0.0]
    dgl.Ambient, dgl.Diffuse = 1.0, 0.0
    v.CameraPosition = [0.0, -10.0, 0.0]
    v.CameraFocalPoint = [0.0, 0.0, 0.0]
    v.CameraParallelScale = half
    for name, fl, il, gll in (("combined", 1, 1, 1), ("color", 1, 0, 0), ("overlay", 0, 1, 1)):
        dfl.Visibility, diso.Visibility, dgl.Visibility = fl, il, gll
        Render(v)
        SaveScreenshot(os.path.join(outdir, "%s_%s_dist_%s%s.png" % (tag, prefix, name, SUFFIX)), v, ImageResolution=[PX, PX], TransparentBackground=1)
        print("wrote", "%s_%s_dist_%s%s.png" % (tag, prefix, name, SUFFIX))
    Hide(gl, v)
    return gl


shots("full", 4.0, 48)
glz = shots("zoom", 2.0, 40)

# geometry of the zoom window
isoZ = Clip(Input=iso, ClipType="Box")
isoZ.ClipType.Position = [-2.0, -0.01, -2.0]
isoZ.ClipType.Length = [4.0, 0.02, 4.0]
isoZ.Invert = 1
isoM = MergeBlocks(Input=isoZ)
isoS = ExtractSurface(Input=isoM)
SaveData(os.path.join(outdir, "%s_zoom_isolines%s.obj" % (tag, SUFFIX)), proxy=isoS)
glM = MergeBlocks(Input=glz)
glS = ExtractSurface(Input=glM)
SaveData(os.path.join(outdir, "%s_zoom_arrows%s.ply" % (tag, SUFFIX)), proxy=glS)

# total velocity colour (full plane)
Hide(iso, v)
ColorBy(dfl, ("POINTS", "tnorm"))
lutt.RescaleTransferFunction(-1.0, 1.0)
GetOpacityTransferFunction("tnorm").RescaleTransferFunction(-1.0, 1.0)
# tnorm is a magnitude (>= 0); use the sequential map on it instead for readability
lutt.ApplyPreset("TU Petrol Sequential", True)
lutt.RescaleTransferFunction(0.0, 1.0)
dfl.Visibility = 1
v.CameraParallelScale = 4.0
Render(v)
SaveScreenshot(os.path.join(outdir, "%s_full_total_color%s.png" % (tag, SUFFIX)), v, ImageResolution=[PX, PX], TransparentBackground=1)
print("wrote total")

# colour bar for the disturbance field
ColorBy(dfl, ("POINTS", "dnorm"))
dfl.SetScalarBarVisibility(v, True)
bar = GetScalarBar(lutd, v)
bar.Title = "|u'| / (gammadot a)"
bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]
bar.TitleFontSize, bar.LabelFontSize = 36, 30
bar.WindowLocation = "Any Location"
bar.Orientation = "Vertical"
bar.Position = [0.25, 0.08]
bar.ScalarBarLength = 0.84
bar.ScalarBarThickness = 40
bar.UseCustomLabels = 1
bar.CustomLabels = [0.0, 0.1, 0.2, 0.3, 0.4]
bar.AddRangeLabels = 0
v.ViewSize = [500, 1600]
v.CameraFocalPoint = [100.0, 0.0, 100.0]
v.CameraPosition = [100.0, -10.0, 100.0]
Render(v)
SaveScreenshot(os.path.join(outdir, "colorbar_dist%s.png" % SUFFIX), v, ImageResolution=[500, 1600], TransparentBackground=1)
print("wrote colorbar")
