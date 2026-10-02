#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver script to perform a hex

"""
import os
import sys
import re
import argparse

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

def convertTri2VtkXml(triFile, vtkFile):
    """
    Converts a tri file to a vtk legacy file

    Attributes:
        triFile: The input tri file
        vtkFile: The name of the output vtk file
    """
    hexMesh = readTriFile(triFile)
    hexMesh.generateMeshStructures()
    writeHexMeshVTKXml(hexMesh, vtkFile)

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

def handleDevisorDir(dirName, outDir="./"):
    """
    Converts files from a devisor directory to vtk
    Parameters:
        dirName: the path to the devisor directory
    """
    fullName = os.path.join(dirName, "NEWFAC")
    nameList = os.listdir(fullName)
    pvdList = []
    print(nameList)
    for item in nameList:
#        print(os.path.splitext(item)[1])
        fileName = os.path.join(fullName, item)
        if os.path.isfile(fileName) and os.path.splitext(item)[1] == '.tri':
            baseName = os.path.splitext(fileName)[0] 
            vtkName = baseName + ".vtk"
            print("Found tri file %s \n" %item)
            print("Writing %s \n" %vtkName)
            convertTri2Vtk(fileName, vtkName)
        elif os.path.isdir(fileName):
            print("Found directory %s \n" %fileName)
#            subDir = os.path.join(fileName, item)
            subList = os.listdir(fileName)
            for subItem in subList:
                if re.search(r'GRID\d\d\d\.tri', subItem):
                #if os.path.splitext(subItem)[1] == '.tri' and len(os.path.splitext(subItem)[0]) > 4:
                #if os.path.splitext(subItem)[1] == '.tri':
                    #print(os.path.join(fileName, subItem))

                    baseName = os.path.splitext(subItem)[0] 
                    vtkName = os.path.join(fileName, baseName + ".vtu")
                    pvdList.append(vtkName)
                    print("Writing %s \n" %vtkName)
                    convertTri2VtkXml(os.path.join(fileName, subItem), vtkName)
    print(pvdList)
    with open(os.path.join(outDir, "pvdfile.pvd"), "w") as f:
        f.write("<?xml version=\"1.0\"?>\n")
        f.write("<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
        f.write("<Collection>\n")

        for fidx, file in enumerate(pvdList):
            f.write("<DataSet timestep=\"0\" part=\"%i\" file=\"%s\"/>\n" %(fidx, file))

        f.write("</Collection>\n")
        f.write("</VTKFile>\n")

def main():


    parser = argparse.ArgumentParser(description='Argument Parser Example')
    parser.add_argument('--input-file', help='Path to input file')
    parser.add_argument('--output-file', help='Path to output file')
    parser.add_argument('--input-dir', help='Path to input directory')
    parser.add_argument('--output-dir', help='Path to output directory')
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file
    input_dir = args.input_dir
    output_dir = args.output_dir

    # Use the arguments as needed in your code
    if input_file and output_file:
        print(f'Input file provided: {input_file}')
        print(f'Output file provided: {output_file}')
        if os.path.splitext(input_file)[1] == '.tri' and os.path.splitext(output_file)[1] == '.vtk':
            convertTri2Vtk(input_file, output_file)
        elif os.path.splitext(input_file)[1] == '.vtk' and os.path.splitext(output_file)[1] == '.tri':
            convertVtk2Tri(input_file, output_file)
    elif input_dir and output_dir:
        print(f'Input directory provided: {input_dir}')
        print(f'Output directory provided: {output_dir}')        
        handleDevisorDir(input_dir, output_dir)
    elif input_dir:
        handleDevisorDir(input_dir)
    else:
        print(f'Please provide input/output files or input/output directories')        
        sys.exit(2)

if __name__ == "__main__":
    main()
