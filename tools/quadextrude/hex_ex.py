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
from mesh_io import *

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
#                        Function computeDzLayers
#===============================================================================
def computeDzLayers(layerSizeZ, nExtrusions):
    """
    Compute the extrusion delta for each level

    Args:
        layerSizeZ: A list of the lengths of each level
        nExtrusions: The number of subdivisions (layers) of each level

    Returns:
        layerDz: A list that contains the DeltaZ of each layer
    """
    layerDz = []
    for idx, item in enumerate(layerSizeZ):
        layerDz.append(item/nExtrusions[idx])

    return layerDz


#===============================================================================
#                      Function:  parseExtrusionLayers
#===============================================================================
def parseExtrusionLayers(argin):
    """
    Parses the input argument <extrusion-layers> that tells us
    the number of extrusions on each level
    """
    ex = argin.split(",")
    ex = [int(x) for x in ex]
    return ex


#===============================================================================
#                      Function:  parseLevelDistance
#===============================================================================
def parseLevelDistance(argin):
    """
    Parses the input argument <distance-levels> which tells us
    the total extrusion length of a particular level
    """
    d = argin.split(",")
    d = [float(x) for x in d]
    return d


#===============================================================================
#                      Function:  calculateSliceIds
#===============================================================================
def calculateSliceIds(extrusionLayers):
    """
    The id for a particular extrusion slice is
    calculated here
    """
    levels = len(extrusionLayers)
    slicesOnLevel = []
    sliceCount = 0
    for value in extrusionLayers:
        slices = []
        for i in range(sliceCount, sliceCount + value):
            slices.append(i)
            sliceCount = sliceCount + 1

        slicesOnLevel.append(slices)

    slicesOnLevel[levels-1].append(sliceCount)

    return slicesOnLevel


#===============================================================================
#                      Function:  parseLevelIds
#===============================================================================
def parseLevelIds(argin):
    """
    Parses the input argument <ids-level> that tells us
    which zone ids should be present on a particular extrusion level
    """
    idsList = []
    levels = argin.split(";")
    for item in levels:
        ids = item.split(":")
        exprList = ids[1].split(",")
        idsLevel = []
        for expr in exprList:
            if expr == "":
                idsLevel.extend([])
                break
            elif "-" in expr:
                values = expr.split('-')
                values = [int(x) for x in values]
                idsLevel.extend(range(values[0], values[1]+1))
            elif not expr.isdigit():
                print("Malformed --ids-level expression: " + expr)
                sys.exit(2)
            else:
                idsLevel.append(int(expr))

        idsList.append(idsLevel)

    return idsList


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
    print("[-h', '--help']: prints this message")
    print("Example: python hex_ex.py  -f testfall.msh --extrusion-layers=3,15,3,5 --levels=4 --distance-levels=7.5,40.0,7.5,5.0 --ids-level='1:;2:3,4,5,6,7,8,9,10;3:2-13;4:11'")
    #python hex_ex.py  -f .\case3\mesh126.msh --extrusion-layers=2,4,2,4 --levels=4 --distance-levels=5.0,5.0,5.0,5.0 --ids-level='1:1,3,4,5,6,7;2:3,4,5,6,7;3:2-8;4:8'
    #python hex_ex.py  -f .\case3\mesh126.msh --extrusion-layers=2,4,2,4 --levels=4 --distance-levels=5.0,5.0,5.0,5.0 --ids-level='1:1,3,4,5,6,7,8,9,10;2:3,4,5,6,7,8,9,10;3:2-12;4:11'

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
    outputFileName = "meshDir/mesh.tri"

    # Default value for levels
    levels = 4

    # The number of layer extrusions on each level
    extrusionLayers = [3, 15, 3, 5]

    # The name of the gmsh input file
    inputName = "testfall.msh"

    # Length of the individual levels
    levelLengthZ = [7.5, 40.0, 7.5, 5.0]

    # The number of layer extrusions on each level
    idsLevel = [[], list(range(3, 11)), list(range(2, 13)), [11]]

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'e:f:l:d:i:h',
                                   ['extrusion-layers=', 'mesh-file=',
                                    'levels=', 'distance-levels=',
                                    'ids-level=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-f', '--mesh-file'):
            inputName = arg
        elif opt in ('-l', '--levels'):
            levels = int(arg)
        elif opt in ('-e', '--extrusion-layers'):
            extrusionLayers = parseExtrusionLayers(arg)
        elif opt in ('-d', '--distance-levels'):
            levelLengthZ = parseLevelDistance(arg)
        elif opt in ('-i', '--ids-level'):
            idsLevel = parseLevelIds(arg)
        else:
            usage()
            sys.exit(2)

    # The slices on each level
    slicesOnLevel = calculateSliceIds(extrusionLayers)

    # Read the input mesh file
    quadMesh = readMeshFile(inputName)

    nNodes = len(quadMesh.nodes)

    qm2 = renumberQuadMesh(quadMesh, 1)
    qm2.generateMeshStructures()
    qm2.checkMeshValidity1()
    writeQuadMeshVTK(quadMesh, 1, "quadmesh.00.vtk")
    writeQuadMeshVTK(qm2, 0, "quadmesh.01.vtk")
#    sys.exit(2)

    offsetNodes = int(nNodes)
    layerNodes = offsetNodes
    offsetHex = 0
    nHex = len(quadMesh.elements)

    layerDz = computeDzLayers(levelLengthZ, extrusionLayers)

    # we already did one quad-extrusion, so
    # adjust to the number of hex-extrusions steps
    extrusionLayers[0] = extrusionLayers[0] - 1

    hm2 = extrudeQuadMeshToHexMesh(quadMesh, layerDz[0])
    hm2.slice = 1
    writeHexMeshVTK(hm2, "caseB.01.vtk")

    hm2.nodesLayer = layerNodes
    hm2.hexasLayer = nHex

    for ilevel in range(levels):
        for i in range(extrusionLayers[ilevel]):
            extrudeHexMeshZ(hm2, offsetHex, offsetNodes,
                            layerNodes, ilevel+1, layerDz[ilevel])
            offsetNodes = offsetNodes + layerNodes
            offsetHex = offsetHex + nHex

    for ilevel in range(levels):
        if not idsLevel[ilevel]:
            continue

        removeHexasLayer(hm2, ilevel+1, idsLevel[ilevel])

    renumberNodes(hm2)

    hm2.generateMeshStructures()

    mkdir("meshDir")
    writeParFiles(hm2, slicesOnLevel)

    writeHexMeshVTK(hm2, "caseB.00.vtk")
    writeTriFile(hm2, outputFileName)


#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
