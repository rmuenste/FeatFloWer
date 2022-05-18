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
#                        parametrize function
#===============================================================================
def parametrize(hexMesh):

    faceIndices = [[0, 1, 2, 3], [0, 1, 4, 5],
                   [1, 2, 5, 6], [3, 2, 6, 7],
                   [0, 3, 7, 4], [4, 5, 6, 7]]

    faceIdx = 0
    hidx = 0
    boundaryComponents = []
    candidates = []
    hex = hexMesh.hexas[hidx]

    startFace = Face()

    for fidx, item in enumerate(hex.neighIdx):
        if item == -1:

            bndryVertices = [hex.nodeIds[faceIndices[fidx][0]],
                             hex.nodeIds[faceIndices[fidx][1]],
                             hex.nodeIds[faceIndices[fidx][2]],
                             hex.nodeIds[faceIndices[fidx][3]]]

            startFace = Face(bndryVertices, faceIdx, 'boundaryFace')

            faceVerts = [hexMesh.nodes[startFace.nodeIds[i]] for i in startFace.nodeIds]
            hexVerts = [hexMesh.nodes[i] for i in hex.nodeIds]

            center = np.zeros(3)
            for vec in hexVerts:
                center = center + vec
            center = center * 0.125

            faceCenter = np.zeros(3)
            for vec in faceVerts:
                faceCenter = faceCenter + vec
            faceCenter = faceCenter * 0.25

            p1 = faceVerts[1] - faceVerts[0]
            p2 = faceVerts[3] - faceVerts[0]
            n0 = np.cross(p1, p2)
            n0 = n0 / np.linalg.norm(n0)

            v1 = n0 - faceCenter
            v2 = center - faceCenter
            if np.dot(v1, v2) > 0.0:
                n0 = -1.0 * n0

            startFace.normal = n0
    
#    compIdx = 0
#    candidates.append(startFace)
#    boundaryComponents.append(set())
#    for hidx in hexMesh.elementsAtBoundary:
#        hex = hexMesh.hexas[hidx]
#        for fidx, item in enumerate(hex.neighIdx):
#            if item == -1:
#                bndryVertices = [hex.nodeIds[faceIndices[fidx][0]],
#                                 hex.nodeIds[faceIndices[fidx][1]],
#                                 hex.nodeIds[faceIndices[fidx][2]],
#                                 hex.nodeIds[faceIndices[fidx][3]]]
#
#                currFace = Face(bndryVertices, faceIdx, 'boundaryFace')
#
#                faceVerts = [hexMesh.nodes[currFace.nodeIds[i]] for i in currFace.nodeIds]
#                hexVerts = [hexMesh.nodes[i] for i in hex.nodeIds]
#
#                center = np.zeros(3)
#                for vec in hexVerts:
#                    center = center + vec
#                center = center * 0.125
#
#                faceCenter = np.zeros(3)
#                for vec in faceVerts:
#                    faceCenter = faceCenter + vec
#                faceCenter = faceCenter * 0.25
#
#                p1 = faceVerts[1] - faceVerts[0]
#                p2 = faceVerts[3] - faceVerts[0]
#                n0 = np.cross(p1, p2)
#                n0 = n0 / np.linalg.norm(n0)
#
#                v1 = n0 - faceCenter
#                v2 = center - faceCenter
#                if np.dot(v1, v2) > 0.0:
#                    n0 = -1.0 * n0
#
#                currFace.normal = n0
#                angle = np.arccos(np.dot(currFace.normal, startFace.normal))
#                if angle < 0.2 * np.pi:
#                    boundaryComponents[compIdx].add(currFace)

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

    # read the mesh from file
    hm = readTriFile("./test_meshes/building_block.tri")

    # generate basic mesh structures
    hm.generateMeshStructures()
    parametrize(hm)

#    for hexa in hm.hexas:
#        print(hexa.nodeIds)
#
#    for i in range(len(hm.nodes)):
#        ex = (hm.nodes[i][0],hm.nodes[i][1],hm.nodes[i][2] + 2.0) 
#        hm.nodes.append(ex)
#
#    for node in hm.nodes:
#        print(node)

    mkdir("NEWFAC")

    dz = 0.01333
    for i in range(1,13):
        mkdir("NEWFAC/sub%03d" %i)
        writeTriFile(hm, "NEWFAC/sub%03d/GRID%03d.tri" %(i, i))
        writeHexMeshVTK(hm, "GRID%03d.vtk" %(i))
        hm.translateMeshZ(dz)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
