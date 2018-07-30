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
            f.write("%f %f %f\n" % (n[0], n[1], n[2]))

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


#===============================================================================
#                        Function readMeshFile
#===============================================================================
#def readMeshFile(fileName):


#===============================================================================
#                        Function readNodes
#===============================================================================
def readNodes(f):

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
    ex = argin.split(",")
    ex = [int(x) for x in ex]
    return ex

#===============================================================================
#                      Function:  parseLevelDistance
#===============================================================================
def parseLevelDistance(argin):
    d = argin.split(",")
    d = [float(x) for x in d]
    return d


#===============================================================================
#                      Function:  parseLevelDistance
#===============================================================================
def parseLevelIds(argin):
    idsList = []
    levels = argin.split(";")
    for item in levels:
        ids = item.split(":")
        if ids[1] == "":
            idsList.append([])
        elif "-" in ids[1]:
            values = ids[1].split('-')
            values = [int(x) for x in values]
            idsList.append(list(range(values[0],values[1])))
        elif "," in ids[1]:
            values = ids[1].split(',')
            values = [int(x) for x in values]
            idsList.append(values)
        else:
            if not ids[1].isdigit(): 
                print("Malformed --ids-level expression: " + ids[1])
                sys.exit(2)
            else:
                idsList.append([int(ids[1])])

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
    outputFileName = "mesh.tri"

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

    #readMeshFile(inputName)
    with open(inputName, "r") as f:

        while True:
            line = f.readline()

            if not line: break

            if re.match(r"^\$Nodes", line):
                readNodes(f)

            if re.match(r"^\$Elements", line):
                readElements(f)

    nNodes = len(myNodes)

    qm = QuadMesh(myQuads, myNodes)
    offsetNodes = int(nNodes)
    layerNodes = offsetNodes
    offsetHex = 0
    nHex = len(qm.quads)

    layerDz = computeDzLayers(levelLengthZ, extrusionLayers)

    # we already did one quad-extrusion, so
    # adjust to the number of hex-extrusions steps
    extrusionLayers[0] = extrusionLayers[0] - 1

    hm2 = extrudeQuadMeshToHexMesh(qm, layerDz[0])
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

    #writeHexMeshVTK(hm2, "hex.00.vtk")
    writeHexMeshVTK(hm2, "caseB.00.vtk")
    writeTriFile(hm2, outputFileName)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
