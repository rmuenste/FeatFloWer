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

def main():
    print(sys.argv)
    hexMesh = readTriFile(sys.argv[1])
    hexMesh.generateMeshStructures()
    writeHexMeshVTK(hexMesh, sys.argv[2])

if __name__ == "__main__":
    main()
