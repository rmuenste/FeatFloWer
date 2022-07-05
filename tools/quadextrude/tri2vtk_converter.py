#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver script to perform a hex

"""
import os
import sys
import re
import getopt

from mesh import *

from shutil import copyfile

def convertTri2Vtk(triFile, vtkFile):
    """
    Converts a tri file to a vtk legacy file

    Attributes:
        triFile: The input tri file
        vtkFile: The name of the output vtk file
    """
    hexMesh = readTriFile(triFile)
    hexMesh.generateMeshStructures()
    writeHexMeshVTK(hexMesh, vtkFile)

def convertTri2VtkXml(triFile, vtkFile):
    """
    Converts a tri file to a vtk legacy file

    Attributes:
        triFile: The input tri file
        vtkFile: The name of the output vtk file
    """
    hexMesh = readTriFile(triFile)
    hexMesh.generateMeshStructures()
    writeHexMeshVTKXml(hexMesh, vtkFile)

def convertVtk2Tri(vtkFile, triFile):
    """
    Converts a vtk legacy file to a tri file
    Attributes:
        vtkFile: The name of the input vtk file
        triFile: The name of the output tri file
    """
    hexMesh = readMeshFromVTK(vtkFile)
    hexMesh.generateMeshStructures()
    writeTriFile(hexMesh, triFile)

def handleDevisorDir(dirName):
    """
    Converts files from a devisor directory to vtk
    Parameters:
        dirName: the path to the devisor directory
    """
    fullName = os.path.join(dirName, "NEWFAC")
    nameList = os.listdir(fullName)
    pvdList = []
    print(nameList)
    for item in nameList:
        print(os.path.splitext(item)[1])
        fileName = os.path.join(fullName, item)
        if os.path.isfile(fileName) and os.path.splitext(item)[1] == '.tri':
            baseName = os.path.splitext(fileName)[0] 
            vtkName = baseName + ".vtk"
            print("Found tri file %s \n" %item)
            print("Writing %s \n" %vtkName)
            convertTri2Vtk(fileName, vtkName)
        elif os.path.isdir(fileName):
            print("Found directory %s \n" %fileName)
#            subDir = os.path.join(fileName, item)
            subList = os.listdir(fileName)
            for subItem in subList:
                if os.path.splitext(subItem)[1] == '.tri' and len(os.path.splitext(subItem)[0]) > 4:
                    #print(os.path.join(fileName, subItem))

                    baseName = os.path.splitext(subItem)[0] 
                    vtkName = os.path.join(fileName, baseName + ".vtu")
                    pvdList.append(vtkName)
                    print("Writing %s \n" %vtkName)
                    convertTri2VtkXml(os.path.join(fileName, subItem), vtkName)
    print(pvdList)
    with open("pvdfile.pvd", "w") as f:
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        f.write("<Collection>\n")

        for fidx, file in enumerate(pvdList):
            f.write("<DataSet timestep=\"0\" part=\"%i\" file=\"%s\"/>\n" %(fidx, file))

        f.write("</Collection>\n")
        f.write("</VTKFile>\n")

def main():
    print(sys.argv)
    if len(sys.argv) < 1:
        raise RuntimeError("The script needs at least one command line parameters to run.")

    firstFile = sys.argv[1]
    if len(sys.argv) > 2:
        secondFile = sys.argv[2]

    if os.path.isdir(sys.argv[1]):
        handleDevisorDir(sys.argv[1])
    elif os.path.splitext(firstFile)[1] == '.tri' and os.path.splitext(secondFile)[1] == '.vtk':
        convertTri2Vtk(firstFile, secondFile)
    elif os.path.splitext(firstFile)[1] == '.vtk' and os.path.splitext(secondFile)[1] == '.tri':
        convertVtk2Tri(firstFile, secondFile)
    else:
        raise RuntimeError("File extensions not suitable for conversion")



if __name__ == "__main__":
    main()
