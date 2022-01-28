#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver script to perform a hex

"""
import os
import sys
import re
import getopt
import json
from random import shuffle, randrange

from mesh import *
from mesh_io import *

def main():
    #(cells, nodes) = readTetMeshFromVTK("arch.vtk")
    #print("Found %i cells" %(len(cells)))
    #print(cells[0:10:1])

    hexMesh = readMeshFromVTK("DistanceMap.01.vtk")
    #hexMesh = readMeshFromVTK("testmesh.vtk")
    print("Hexas %i" %len(hexMesh.hexas))
    myList = range(0,len(hexMesh.hexas))
    #myList = range(0,11)
    shuffle(myList)
    #randrange(10)
    #print(myList)

    coords = []

    for idx in myList[0:1000]:
        idx1 = hexMesh.hexas[idx].nodeIds[0]
        idx2 = hexMesh.hexas[idx].nodeIds[1]
        idx3 = hexMesh.hexas[idx].nodeIds[2]
        idx4 = hexMesh.hexas[idx].nodeIds[3]
        idx5 = hexMesh.hexas[idx].nodeIds[4]
        idx6 = hexMesh.hexas[idx].nodeIds[5]
        idx7 = hexMesh.hexas[idx].nodeIds[6]
        idx8 = hexMesh.hexas[idx].nodeIds[7]

        centerX = (hexMesh.nodes[idx1][0] + hexMesh.nodes[idx2][0] + 
                   hexMesh.nodes[idx3][0] + hexMesh.nodes[idx4][0] +
                   hexMesh.nodes[idx5][0] + hexMesh.nodes[idx6][0] +
                   hexMesh.nodes[idx7][0] + hexMesh.nodes[idx8][0]) * 0.125
        centerY = (hexMesh.nodes[idx1][1] + hexMesh.nodes[idx2][1] + 
                   hexMesh.nodes[idx3][1] + hexMesh.nodes[idx4][1] +
                   hexMesh.nodes[idx5][1] + hexMesh.nodes[idx6][1] +
                   hexMesh.nodes[idx7][1] + hexMesh.nodes[idx8][1]) * 0.125
        centerZ = (hexMesh.nodes[idx1][2] + hexMesh.nodes[idx2][2] + 
                   hexMesh.nodes[idx3][2] + hexMesh.nodes[idx4][2] +
                   hexMesh.nodes[idx5][2] + hexMesh.nodes[idx6][2] +
                   hexMesh.nodes[idx7][2] + hexMesh.nodes[idx8][2]) * 0.125

        coord = [centerX, centerY, centerZ]
        coords.append(coord)

    myJson = []
    for coord in coords:
        entry = {
                "Type":"Sphere",
                "IsDynamic":"1",
                "Pos":[coord[0],coord[1],coord[2]],
                "Rot":[0.0,0.0,0.0],
                "Norm":[1.0,0.0,0.0],
                "Dim":[0.03,0.03,0.03]
                }
        myJson.append(entry)

    #print(json.dumps(myJson))
    with open("archimedes_test.json", "w") as myOut:
        json.dump(myJson, myOut, indent=4)

    writePointsVTK(coords, "points.vtk")

if __name__ == "__main__":
    main()