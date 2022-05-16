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

def main():
    print(sys.argv)
    if len(sys.argv) < 2:
        raise RuntimeError("The script needs at least two command line parameters to run.")

    firstFile = sys.argv[1]
    secondFile = sys.argv[2]

    if os.path.splitext(firstFile)[1] == '.tri' and os.path.splitext(secondFile)[1] == '.vtk':
        convertTri2Vtk(firstFile, secondFile)
    elif os.path.splitext(firstFile)[1] == '.vtk' and os.path.splitext(secondFile)[1] == '.tri':
        convertVtk2Tri(firstFile, secondFile)
    else:
        raise RuntimeError("File extensions not suitable for conversion")



if __name__ == "__main__":
    main()
