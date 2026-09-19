"""Export Blender resources for the D5.2 numerical viscometer (annular Couette, neutrally buoyant spheres).

Usage: pvpython --mesa --force-offscreen-rendering export_blender_viscometer.py PVTU OUTDIR TAG
Geometry (box units, d = 1): inner cylinder r_i = 5 rotating at Omega = 0.1 (wall speed 0.5), outer r_a = 10 fixed,
height z in [0, 10]; spheres d = 1.

Produces in OUTDIR:
  TAG_mid_{color,overlay,combined}.png   horizontal slice z = 5, x,y in [-10,10], 4000 x 4000 px (200 px per unit)
  TAG_vert_{color,overlay,combined}.png  vertical slice y = 0, x in [-10,10], z in [0,10], 4000 x 2000 px
  colorbar_visco.png                     legend |u| / U_wall, 0 to 1
Colour: |u| / U_wall with U_wall = 0.5. Sphere footprints are cut from the indicator field (Mixer) and are transparent.
Isolines at 0.1 .. 0.9. No arrows (the flow is azimuthal; the isolines carry the profile).
"""
import os
import sys
from paraview.simple import *

pvtu, outdir, tag = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)
# exact sphere holes from the logged centres (CSV written beside the textures): q = min_i |x - x_i|^2 / R^2
import numpy as np
cen = np.loadtxt(os.path.join(outdir, "%s_spheres_t250.csv" % tag), delimiter=",", skiprows=2, usecols=(1, 2, 3))
np.save(os.path.join("/tmp/rmuenste/claude-3086/-data-warehouse17-rmuenste-code-FF-EL-FeatFloWer/8f6ed8f2-2f2b-4009-9047-7ed8c2680fb3/scratchpad/viz", "_centres_%s.npy" % tag), cen)
CEN_NPY = os.path.join("/tmp/rmuenste/claude-3086/-data-warehouse17-rmuenste-code-FF-EL-FeatFloWer/8f6ed8f2-2f2b-4009-9047-7ed8c2680fb3/scratchpad/viz", "_centres_%s.npy" % tag)
RSPH = 0.5
UW = 0.5
PRESET = "/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/tu-colormaps.json"
LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

rd = XMLPartitionedUnstructuredGridReader(FileName=[pvtu])
rd.PointArrayStatus = ["Velocity", "Mixer"]
calc = Calculator(Input=rd, ResultArrayName="unorm", Function="mag(Velocity)/%g" % UW)

v = GetActiveViewOrCreate("RenderView")
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
lut.Discretize = 1
lut.NumberOfTableValues = 20
lut.RescaleTransferFunction(0.0, 1.0)
GetOpacityTransferFunction("unorm").RescaleTransferFunction(0.0, 1.0)


def plane(origin, normal, name, W, H, campos, focal, up, scale):
    sl = Slice(Input=calc)
    sl.SliceType = "Plane"
    sl.SliceType.Origin = origin
    sl.SliceType.Normal = normal
    pf = ProgrammableFilter(Input=sl)
    pf.Script = """
import numpy as np
from vtk.numpy_interface import dataset_adapter as dsa
inp = dsa.WrapDataObject(inputs[0].VTKObject)
pts = np.asarray(inp.Points)
cen = np.load(%r)
q = np.full(len(pts), 1e9)
for c in cen:
    d2 = ((pts - c) ** 2).sum(axis=1) / (%g ** 2)
    q = np.minimum(q, d2)
out = dsa.WrapDataObject(output.VTKObject)
out.ShallowCopy(inputs[0].VTKObject)
out.PointData.append(q, "q")
""" % (CEN_NPY, RSPH)
    fluid = Clip(Input=pf, ClipType="Scalar")
    fluid.Scalars = ["POINTS", "q"]
    fluid.Value = 0.85           # keep q > 0.85: fluid; the hole sits 8 percent inside each sphere
    fluid.Invert = 0
    fluid.UpdatePipeline()
    print(name, "unorm range", fluid.PointData["unorm"].GetRange())
    iso = Contour(Input=fluid, ContourBy=["POINTS", "unorm"], Isosurfaces=LEVELS)
    v.ViewSize = [W, H]
    dfl = Show(fluid, v); ColorBy(dfl, ("POINTS", "unorm")); dfl.Ambient, dfl.Diffuse = 1.0, 0.0
    lut.RescaleTransferFunction(0.0, 1.0)
    diso = Show(iso, v); ColorBy(diso, None); diso.AmbientColor = diso.DiffuseColor = [0.05, 0.05, 0.05]
    diso.Ambient, diso.Diffuse = 1.0, 0.0; diso.LineWidth = 2.5
    Render(v)
    v.InteractionMode = "2D"; v.CameraParallelProjection = 1
    v.CameraPosition = campos; v.CameraFocalPoint = focal; v.CameraViewUp = up; v.CameraParallelScale = scale
    for suffix, fl, il in (("combined", 1, 1), ("color", 1, 0), ("overlay", 0, 1)):
        dfl.Visibility, diso.Visibility = fl, il
        Render(v)
        SaveScreenshot(os.path.join(outdir, "%s_%s_%s.png" % (tag, name, suffix)), v, ImageResolution=[W, H], TransparentBackground=1)
        print("wrote", "%s_%s_%s.png" % (tag, name, suffix))
    Hide(fluid, v); Hide(iso, v)
    return dfl


# horizontal mid-height slice, seen from above (+z), +x right, +y up
plane([0, 0, 5.0], [0, 0, 1], "mid", 4000, 4000, [0, 0, 50.0], [0, 0, 5.0], [0, 1, 0], 10.0)
# vertical slice through the axis (y = 0), seen from -y, +x right, +z up
dfl = plane([0, 0, 0], [0, 1, 0], "vert", 4000, 2000, [0, -50.0, 5.0], [0, 0, 5.0], [0, 0, 1], 5.0)

# colour bar
dfl.Visibility = 1
dfl.SetScalarBarVisibility(v, True)
bar = GetScalarBar(lut, v)
bar.Title = "|u| / U_wall"; bar.ComponentTitle = ""
bar.TitleColor = bar.LabelColor = [0, 0, 0]; bar.TitleFontSize, bar.LabelFontSize = 36, 30
bar.WindowLocation = "Any Location"; bar.Orientation = "Vertical"; bar.Position = [0.25, 0.08]
bar.ScalarBarLength = 0.84; bar.ScalarBarThickness = 40
bar.UseCustomLabels = 1; bar.CustomLabels = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]; bar.AddRangeLabels = 0
v.ViewSize = [500, 1600]
v.CameraFocalPoint = [200.0, 0.0, 200.0]; v.CameraPosition = [200.0, -50.0, 200.0]
Render(v)
SaveScreenshot(os.path.join(outdir, "colorbar_visco.png"), v, ImageResolution=[500, 1600], TransparentBackground=1)
print("wrote colorbar")
