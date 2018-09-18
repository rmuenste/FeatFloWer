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

myQuads = []
myNodes = []

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
#                      Function:  Write Tri File
#===============================================================================
def writeTriFile(hexMesh, fileName):
    """
    Writes out a hexMesh in TRI format

    Args:
        hexMesh: A reference to a HexMesh class
        fileName: The file name of the TRI file

    """
    with open(fileName, "w") as f:
        f.write("Coarse mesh exported by hex_ex.py \n")
        f.write("Parametrisierung PARXC, PARYC, TMAXC \n")
        f.write("%i %i 0 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE \n" % (len(hexMesh.hexas), len(hexMesh.nodes)))
        f.write("DCORVG\n")

        for n in hexMesh.nodes:
            f.write("%f %f %f\n" % (0.1 * n[0], 0.1 * n[1], 0.1 * n[2]))

        f.write("KVERT\n")
        for h in hexMesh.hexas:
            indices = (int(h.nodeIds[0]), int(h.nodeIds[1]),
                       int(h.nodeIds[2]), int(h.nodeIds[3]),
                       int(h.nodeIds[4]), int(h.nodeIds[5]),
                       int(h.nodeIds[6]), int(h.nodeIds[7]))

            f.write('%i %i %i %i %i %i %i %i\n' % (indices[0]+1, indices[1]+1,
                                                   indices[2]+1, indices[3]+1,
                                                   indices[4]+1, indices[5]+1,
                                                   indices[6]+1, indices[7]+1))

        f.write("KNPR\n")
        for n in hexMesh.nodes:
            f.write("0\n")


#===============================================================================
#               A very simple VTK polygon writer
#===============================================================================
def writeQuadMeshVTK(quadMesh):
    """
    Writes out a quadMesh in a very simple VTK format

    Args:
        quadMesh: A reference to a QuadMesh class
        fileName: The file name of the VTK file

    """

    with open("myvtk.00.vtk", "w") as f:
        f.write("# vtk DataFile Version 4.2 \n")
        f.write("vtk output \n")
        f.write("ASCII \n")

        nVertices = len(quadMesh.nodes)
        f.write("DATASET POLYDATA\n")
        f.write("POINTS " + str(nVertices) + " float\n")
        for n in quadMesh.nodes:
            f.write('%s %s %s\n' % (n[1], n[2], n[3]))

        nElem = len(quadMesh.quads)
        f.write("CELLS " + str(nElem) + " " + str(nElem * 5) + " \n")

        for q in quadMesh.quads:
            indices = (int(q[5])-1, int(q[6])-1, int(q[7])-1, int(q[8])-1)
            f.write('4 %i %i %i %i\n' % (indices[0], indices[1],
                                         indices[2], indices[3]))


#===============================================================================
#                A very simple VTK Unstructured Grid writer
#===============================================================================
def writeHexMeshVTK(hexMesh, fileName):
    """
    Writes out a hexMesh in a very simple VTK format

    Args:
        hexMesh: A reference to a HexMesh class
        fileName: The file name of the VTK file

    """

    with open(fileName, "w") as f:
        f.write("# vtk DataFile Version 4.2 \n")
        f.write("vtk output \n")
        f.write("ASCII \n")

        nVertices = len(hexMesh.nodes)
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS " + str(nVertices) + " float\n")
        for n in hexMesh.nodes:
            f.write('%s %s %s\n' % (n[0], n[1], n[2]))

        nElem = len(hexMesh.hexas)
        f.write("CELLS " + str(nElem) + " " + str(nElem * 9) + " \n")

        for h in hexMesh.hexas:
            indices = (int(h.nodeIds[0]), int(h.nodeIds[1]),
                       int(h.nodeIds[2]), int(h.nodeIds[3]),
                       int(h.nodeIds[4]), int(h.nodeIds[5]),
                       int(h.nodeIds[6]), int(h.nodeIds[7]))

            f.write('8 %i %i %i %i %i %i %i %i\n' % (indices[0], indices[1],
                                                     indices[2], indices[3],
                                                     indices[4], indices[5],
                                                     indices[6], indices[7]))

        f.write("CELL_TYPES " + str(nElem) + " \n")
        for h in hexMesh.hexas:
            f.write('12\n')

        f.write("CELL_DATA " + str(nElem) + " \n")
        f.write("SCALARS ZoneId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for h in hexMesh.hexas:
            f.write('%i\n' % (h.type))

        f.write("SCALARS LayerId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for h in hexMesh.hexas:
            f.write('%i\n' % (h.layerIdx))

        f.write("POINT_DATA " + str(nVertices) + " \n")
        f.write("SCALARS KNPR integer\n")
        f.write("LOOKUP_TABLE default\n")
        for val in hexMesh.verticesAtBoundary:
            f.write('%i\n' % (val))

        f.write("SCALARS SliceId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for val in hexMesh.nodesAtSlice:
            f.write('%i\n' % (val))


#===============================================================================
#                        Function readMeshFile
#===============================================================================
def readMeshFile(fileName):
    with open(fileName, "r") as f:

        while True:
            line = f.readline()

            if not line:
                break

            if re.match(r"^\$Nodes", line):
                readNodes(f)

            if re.match(r"^\$Elements", line):
                readElements(f)


#===============================================================================
#                        Function readNodes
#===============================================================================
def readNodes(f):
    """
    Reader for the nodes section of a .msh file

    Args:
        f: the file handle to the msh file 

    """

    line = f.readline()
    if not line:
        return

    line = f.readline()

    while line and not re.match(r"^\$EndNodes", line):

        if not line: break

        words = line.strip().split(" ")
        myNodes.append((words[0], words[1], words[2], words[3]))

        line = f.readline()


#===============================================================================
#                        Function readElements
#===============================================================================
def readElements(f):
    """
    Reader for the elements section of a .msh file

    Args:
        f: the file handle to the msh file 

    """

    line = f.readline()
    if not line:
        return

    line = f.readline()

    while line and not re.match(r"^\$EndElements", line):

        if not line: break

        words = line.strip().split(" ")
        if words[1] == "3":
            myQuads.append((words[0], words[1], words[2],
                            words[3], words[4], words[5],
                            words[6], words[7], words[8]))

        line = f.readline()


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
    readMeshFile(inputName)

    nNodes = len(myNodes)

    quadMesh = QuadMesh(myQuads, myNodes)
    offsetNodes = int(nNodes)
    layerNodes = offsetNodes
    offsetHex = 0
    nHex = len(quadMesh.quads)

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
