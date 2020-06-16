#!/usr/bin/env python
# vim: set filetype=python
"""
A python launcher script for a FeatFloWer application
"""
from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
import os
import shutil
import sys
try:
    sys.path.append(os.environ['FF_PY_HOME'])
except:
    pass
import partitioner
import xml.etree.ElementTree as ET
import getopt
import platform
import subprocess
import re
import json
import pprint

if sys.version_info[0] < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

#===============================================================================
#               Remove off files from working dir
#===============================================================================
def cleanWorkingDir(workingDir):
    if not sys.platform == "win32":
        offList = list(workingDir.glob('*.off')) + list(workingDir.glob('*.OFF'))
    else: 
        offList = list(workingDir.glob('*.off'))

    for item in offList:
        os.remove(str(item))

#===============================================================================
#                          setup the case folder 
#===============================================================================
def folderSetup(workingDir, projectFolder):
    offList = []
    if not sys.platform == "win32":
        offList = list(Path(projectFolder).glob('*.off')) + list(Path(projectFolder).glob('*.OFF'))
    else: 
        offList = list(Path(projectFolder).glob('*.off'))

    for item in offList:
        shutil.copyfile(str(item), str(workingDir / item.name))

#===============================================================================
#                      Function: Usage
#===============================================================================
def usage():
    print("Usage: heat_start.py [options]")
    print("Where options can be:")
    print("[-h, --help]: prints this message")
    print("[-n, --num-processors]: defines the number of parallel jobs to be used")
    print("[-f, --in-folder]: defines the input project file (will be replaced in the q2p1-data-file)")

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                #new_file.write(line.replace(pattern, subst))
                if pattern in line:
                    new_file.write(line.replace(line, subst + "\n"))
                else:
                    new_file.write(line)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def main(argv):
    inputFile = '_data/meshDir/file.prj'
    inputCaseFolder = ''
    inputCaseFile = ''
    outputFile = ''
    useSrun = False

    numProcessors = -1

    try:
        opts, args = getopt.getopt(argv, "huf:n:", 
                                ["help", "use-srun", "in-folder=", "num-processors="])
    except getopt.GetoptError:
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-n", "--num-processors"):
            numProcessors = int(arg)
        elif opt in ("-f", "--in-folder"):
            inputCaseFolder = arg
        elif opt in ('-u', '--use-srun'):
            useSrun = True
        else:
            usage()
            sys.exit(2)

    print('Number Of Processors is: '+ str(numProcessors))
    
    inputCaseFile = inputCaseFolder+'/heat.s3d'
    print('Input case file is '+ inputCaseFile)

    workingDir = Path('.')
    folderSetup(workingDir, inputCaseFolder)
    
    shutil.copyfile(inputCaseFile, "_data/heat.s3d")

    #replace("_data/q2p1_param.dat",
            #"SimPar@ProjectFile = ",
            #"SimPar@ProjectFile = '" + inputFile + "'")

    if os.path.exists("_data/meshDir"):
        shutil.rmtree("_data/meshDir")

    subprocess.call(['./s3d_mesher -a %s' %('heat')], shell=True)

    # Call the partitioner
    partitioner.partition(numProcessors - 1, 1, 1, "NEWFAC", str(inputFile))

    # Configure the launch command and start a simulation
    launchCommand = ""
    if useSrun:
        launchCommand = "srun ./heat"
    else:
        launchCommand = 'mpirun -np %i ./%s' %(numProcessors, 'heat')

    # Start the simulation as a subprocess
    subprocess.call([launchCommand], shell=True)

    cleanWorkingDir(workingDir)

if __name__ == "__main__":
    main(sys.argv[1:])
