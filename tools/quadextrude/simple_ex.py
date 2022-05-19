#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver script to perform a hex

"""
import os
import sys
import re
import getopt
from tracemalloc import start
import numpy as np

#from mesh import *
#from mesh_io import *
from mesh import *

from shutil import copyfile

#===============================================================================
#                      Function: mkdir
#===============================================================================
def mkdir(dir):
    """
    Erzeugt nur dann das Verzeichnis "dir", wenn es noch nicht vorhanden ist.
    Falls eine Datei dieses Namens existieren sollte, wird sie durch das Verzeichnis ersetzt.
    """
    if os.path.exists(dir):
        if os.path.isdir(dir):
            return
        else:
            os.remove(dir)
    os.mkdir(dir)

#===============================================================================
#                        usage function
#===============================================================================
def usage():
    """
    Print out usage information
    """
    print("Usage: python hex_ex.py [options]")
    print("Where options can be:")
    print("[-d', '--distance-levels]: Comma-separated list of the extrusion distance of each level")
    print("[-e', '--extrusion-layers]: Comma-separated list of the number of extrusion layers on each level")
    print("[-f', '--mesh-file']: The path to the gmsh .msh input file")
    print("[-l', '--levels]: The number of extrusion levels")
    print("[-i', '--ids-levels]: The element ids that should be present on each level, this string should be put into quotation marks")
    print("[-o', '--output-dir]: The output directory for the .tri mesh")
    print("[-h', '--help']: prints this message")
    print("Example: python hex_ex.py  -f testfall.msh --extrusion-layers=3,15,3,5 --output-dir=meshDir --levels=4 --distance-levels=7.5,40.0,7.5,5.0 --ids-level='1:;2:3,4,5,6,7,8,9,10;3:2-13;4:11'")
    # Commands for cases:
    #python hex_ex.py  -f .\case1\mesh.msh --extrusion-layers=2,4,2,4 --levels=4 --distance-levels=5.0,5.0,5.0,5.0 --ids-level='1:1,3,4,5,6,7;2:3,4,5,6,7;3:2-8;4:8'
    #python hex_ex.py  -f .\case2\mesh.msh --extrusion-layers=2,4,2,4 --levels=4 --distance-levels=5.0,5.0,5.0,5.0 --ids-level='1:1,3,4,5,6,7;2:3,4,5,6,7;3:2-8;4:8'
    #python hex_ex.py  -f .\case3\mesh126.msh --extrusion-layers=2,4,2,4 --levels=4 --distance-levels=5.0,5.0,5.0,5.0 --ids-level='1:1,3,4,5,6,7,8,9,10;2:3,4,5,6,7,8,9,10;3:2-12;4:11'

#===============================================================================
#                        rad2deg function
#===============================================================================
def rad2deg(val):
    return 180.0/np.pi

#===============================================================================
#                        deg2rad function
#===============================================================================
def deg2rad(val):
    return np.pi/180.0

#===============================================================================
#                        deg2rad function
#===============================================================================
def filterFace(face, normal, tol, bndryCompList):
    angle = np.arccos(np.dot(face.normal, normal))
    if angle < tol:
        bndryCompList.append(face.idx)
        return False 
    else:
        return True

#===============================================================================
#                      parametrizeFaces function
#===============================================================================
def parametrizeFaces(hexMesh, degrees):

    workList = hexMesh.facesAtBoundary[:]
    boundaryComponents = []
    compIdx = 0
    tol = deg2rad(degrees)

    while not len(workList) == 0:
        firstFace = workList[0]
        boundaryComponents.append([])
        workList = list(filter( lambda face: filterFace(face, firstFace.normal, tol, boundaryComponents[compIdx]), workList))
        compIdx = compIdx + 1

    for item in boundaryComponents:
        item.sort()


    print("==============================================")
    print("Number of boundary components: ", len(boundaryComponents))
    return boundaryComponents

#===============================================================================
#                     writeBoundaryComponents
#===============================================================================
def writeBoundaryComponentsFaces(hexMesh, boundaryComponents):
    for idx, item in enumerate(boundaryComponents):
        with open("bcf%d.par" %idx, "w") as parFile:
            parFile.write("%d Wall\n" % len(item))
            normal = hexMesh.facesAtBoundary[item[0]].normal
            parFile.write("'4 %f %f %f'\n" % (normal[0], normal[1], normal[2]))
            for val in item:
                parFile.write("%d\n" % val)

#===============================================================================
#                        Main Script Function
#===============================================================================
def main():
    """
    The main function that controls the extrusion process

    Options:
        mesh-file: The file name of the input gmsh .msh file
        levels: The number of main extrusion levels
        layers: The number of subdivision layers on the levels
    """
    # Default output file name
    outputDirName = "meshDir"

    # Default output file name
    outputFileName = outputDirName + "/mesh.tri"

    # The number of layer extrusions on each level
    extrusionLayers = [3, 15, 3, 5]

    # The name of the gmsh input file
    inputName = "testfall.msh"

    meshQualityOK = True

    print(os.getcwd())
    # read the mesh from file
    hm = readTriFile("./test_meshes/building_block.tri")

    # generate basic mesh structures
    hm.generateMeshStructures()
    boundaryComponents = parametrizeFaces(hm, 20.0)
    hm.parametrizeVertices(boundaryComponents)
    #writeBoundaryComponentsFaces(hm, boundaryComponents)

    mkdir("NEWFAC")

    dz = 0.01333
    for i in range(1,13):
        mkdir("NEWFAC/sub%03d" %i)
        meshName = "NEWFAC/sub%03d/GRID%03d.tri" %(i, i)
        writeTriFile(hm, meshName)
        writeBoundaryComponents(hm, "NEWFAC/sub%03d" %(i), meshName)
        writeHexMeshVTK(hm, "GRID%03d.vtk" %(i))
        hm.translateMeshZ(dz)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
