from paraview.simple import *
import sys
rd = XMLPartitionedUnstructuredGridReader(FileName=[sys.argv[1]]); rd.PointArrayStatus=["Velocity"]; rd.UpdatePipeline()
di=rd.GetDataInformation(); print("cells",di.GetNumberOfCells(),"bounds",di.GetBounds()); print("Vel range",rd.PointData["Velocity"].GetRange(-1), "vz range", rd.PointData["Velocity"].GetRange(2))
