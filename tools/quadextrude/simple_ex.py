#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver script to perform a hex

"""
import os
import sys
import re
import getopt

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

#    writeHexMeshVTK(hm2, "caseB.01.vtk")
#    writeTriFile(hm2, outputFileName)
    hm = readTriFile("D:/code/FeatFloWer/Feat_FloWer/tools/quadextrude/meshDir/cube.tri")

    for hexa in hm.hexas:
        print(hexa.nodeIds)

    for i in range(len(hm.nodes)):
        ex = (hm.nodes[i][0],hm.nodes[i][1],hm.nodes[i][2] + 2.0) 
        hm.nodes.append(ex)

    for node in hm.nodes:
        print(node)

    writeHexMeshVTK(hm, "cube.vtk")


#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
