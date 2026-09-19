from paraview.simple import *
import paraview.servermanager as sm
lut = GetColorTransferFunction("unorm")
r = lut.ApplyPreset("Viridis (matplotlib)", True)
print("apply Viridis ->", r, "RGBPoints head", list(lut.RGBPoints)[:8])
r = lut.ApplyPreset("Inferno (matplotlib)", True); print("apply Inferno ->", r, list(lut.RGBPoints)[:4])
r = lut.ApplyPreset("Turbo", True); print("apply Turbo ->", r, list(lut.RGBPoints)[:4])
print("props:", [p for p in lut.ListProperties() if "iscret" in p or "Table" in p])
lut.Discretize = 1; lut.NumberOfTableValues = 25; print("Discretize", lut.Discretize, lut.NumberOfTableValues)
# list preset names containing keywords
try:
    from paraview import servermanager
    presets = servermanager.vtkSMTransferFunctionPresets.GetInstance()
    names = [presets.GetPresetName(i) for i in range(presets.GetNumberOfPresets())]
    print("n presets", len(names)); print([n for n in names if any(k in n.lower() for k in ("viridis","inferno","turbo","magma","plasma","cividis"))])
except Exception as e: print("preset list err", e)
