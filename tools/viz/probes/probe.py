from paraview.simple import *
s=Sphere(); r=Reflect(Input=s); print("Reflect:", r.ListProperties()); print("Plane options:", r.GetProperty("Plane").GetAvailable())
w=PLYWriter(Input=s, FileName="/dev/null"); print("PLYWriter:", w.ListProperties())
import paraview.simple as ps; print("exporters:", [n for n in dir(ps) if "Export" in n])
