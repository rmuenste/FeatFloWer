#/usr/bin/env python
# vim: set filetype=python
"""
A module for input/output of different mesh formats
"""

import re

from mesh import *

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

    with open("quadmesh.00.vtk", "w") as f:
        f.write("# vtk DataFile Version 4.2 \n")
        f.write("vtk output \n")
        f.write("ASCII \n")

        nVertices = len(quadMesh.nodes)
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write("POINTS " + str(nVertices) + " float\n")
        for n in quadMesh.nodes:
            f.write('%s %s %s\n' % (n[1], n[2], n[3]))

        nElem = len(quadMesh.elements)
        f.write("CELLS " + str(nElem) + " " + str(nElem * 5) + " \n")

        for q in quadMesh.elements:
            indices = (q.nodeIds[0]-1, q.nodeIds[1]-1, q.nodeIds[2]-1, q.nodeIds[3]-1)
            f.write('4 %i %i %i %i\n' % (indices[0], indices[1],
                                         indices[2], indices[3]))

        f.write("CELL_TYPES " + str(nElem) + " \n")
        for q in quadMesh.elements:
            f.write('9\n')

        f.write("CELL_DATA " + str(nElem) + " \n")
        f.write("SCALARS ZoneId integer\n")
        f.write("LOOKUP_TABLE default\n")
        for e in quadMesh.elements:
            f.write('%i\n' % (e.zoneId))


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

    quadList = []
    nodesList = []
    with open(fileName, "r") as f:

        while True:
            line = f.readline()

            if not line:
                break

            if re.match(r"^\$Nodes", line):
                nodesList = readNodes(f)

            if re.match(r"^\$Elements", line):
                quadList = readElements(f)

    return QuadMesh(nodesList, quadList)


#===============================================================================
#                        Function readNodes
#===============================================================================
def readNodes(f):
    """
    Reader for the nodes section of a .msh file

    Args:
        f: the file handle to the msh file 

    The expected format of an entry in the $Nodes section of the file is:
    <node-number x y z>
    """

    meshNodes = []

    line = f.readline()
    if not line:
        return

    line = f.readline()

    while line and not re.match(r"^\$EndNodes", line):

        if not line: break

        words = line.strip().split(" ")
        meshNodes.append((words[0], words[1], words[2], words[3]))

        line = f.readline()

    return meshNodes

#===============================================================================
#                        Function readElements
#===============================================================================
def readElements(f):
    """
    Reader for the elements section of a .msh file

    Args:
        f: the file handle to the msh file 

    The expected format of an entry in the $Elements section of the file is:
    <elem-number elem-type number-of-tags 'number-of-tags tags ...' node-number-list>
    """

    quads = []
    line = f.readline()
    if not line:
        return

    line = f.readline()

    while line and not re.match(r"^\$EndElements", line):

        if not line: break

        quadCnt = 0
        words = line.strip().split(" ")
        if words[1] == "3":
            nodeIds = [int(words[i]) for i in range(5, 9)]
            quadElem = Quad(nodeIds, quadCnt, int(words[4]))
            quads.append(quadElem)
            quadCnt = quadCnt + 1

        line = f.readline()

    return quads

