from paraview.simple import *
s = Sphere(Radius=1.0, ThetaResolution=64, PhiResolution=64)
v = GetActiveViewOrCreate('RenderView')
v.ViewSize = [800, 600]
d = Show(s, v); d.DiffuseColor = [0.85, 0.35, 0.2]
v.Background = [1,1,1]
try: v.UseColorPaletteForBackground = 0
except: pass
ResetCamera(v)
SaveScreenshot('/tmp/rmuenste/claude-3086/-data-warehouse17-rmuenste-code-FF-EL-FeatFloWer/8f6ed8f2-2f2b-4009-9047-7ed8c2680fb3/scratchpad/viz/smoke.png', v, ImageResolution=[800,600])
print("rendered")
