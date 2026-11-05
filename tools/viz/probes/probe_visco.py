from paraview.simple import *
import sys
rd = XMLPartitionedUnstructuredGridReader(FileName=[sys.argv[1]]); rd.PointArrayStatus=["Velocity","Mixer"]; rd.UpdatePipeline()
di=rd.GetDataInformation(); print("cells",di.GetNumberOfCells(),"points",di.GetNumberOfPoints(),"bounds",di.GetBounds())
print("Mixer range",rd.PointData["Mixer"].GetRange(),"Vel range",rd.PointData["Velocity"].GetRange(-1))
c=Calculator(Input=rd,ResultArrayName="r",Function="sqrt(coordsX^2+coordsY^2)")
for lo,hi,name in ((5.0,5.15,"inner"),(9.85,10.0,"outer"),(7.4,7.6,"mid")):
    th=Threshold(Input=c,Scalars=["POINTS","r"],LowerThreshold=lo,UpperThreshold=hi,ThresholdMethod="Between"); th.UpdatePipeline()
    print(name,"r in",lo,hi,"|u| range",th.PointData["Velocity"].GetRange(-1))
