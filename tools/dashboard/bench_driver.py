#!/usr/bin/env python
# vim: set filetype=python
"""
A python launcher script for a FeatFloWer application
"""
import os

import sys
import getopt
import platform
import subprocess
import re
import json
sys.path.append(os.environ['FF_PY_HOME'])
import partitioner

#===============================================================================
#                      Function: Usage
#===============================================================================
def usage():
    print("Usage: configure [options]")
    print("Where options can be:")
    print("[-h, --help]: prints this message")
    print("[-h, --help]: prints this message")


#===============================================================================
#                      Function:  moveAndSetLevel
#===============================================================================
def moveAndSetLevel(fileIn, fileOut, level):
    maxLevelStr = "SimPar@MaxMeshLevel = " + str(level)
    with open(fileOut, "w") as newFile:
        with open(fileIn, "r") as file:
            for line in file:
                newLine = re.sub(r"^[\s]*SimPar@MaxMeshLevel[\s]*=(\s | \w)*", maxLevelStr, line)
                newFile.write(newLine)


#===============================================================================
#                      Function: Get log entry
#===============================================================================
def getLogEntry(fileName, varName):
    with open(fileName, "r") as sources:
        lines = sources.readlines()

    tFound = False

    for line in reversed(lines):
        m = re.match(varName, line)
        if m != None:
            splitLine = line.split(':')
            val = splitLine[1].strip()
            tFound = True
            break

    if tFound:
        return val
    else:
        return 0


#===============================================================================
#                              Main function
#===============================================================================
def main():
    """ Preliminary parameter list:
        ======================
        Simulation Parameters 
        ======================

        num-processors : Number of processors
        project-file : The path to the project file
        parameter-file : The path to the parameter file
        binary-name : The name of the executable 
        levels : The number of additional levels w.r.t the 
                 maxlevel given in the initial 
                 parameter-file that the test will be performed
                 on.
    
        ======================
         Dashboard Parameters 
        ======================
    
        benchName : The identifier of the benchmark 
        style : The visualization style of the test
        tableCaption : The caption if the table style is chosen 
        ouput-file : The path to the JSON output file 
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], 't:h', ['test-configuration=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
  
    testConfiguration = 16
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-t', '--test-configuration'):
            if arg.isdigit():
                testConfiguration = int(arg)
            else:
                print("--test-configuration= " + str(arg) + " is not a valid number.")
                sys.exit(2)
        else:
            usage()
            sys.exit(2)
    
    loadgccmpi81='module purge && module load cmake/3.13.2 gcc/8.1.0 boost/1.53 openmpi/2.1.5 && export CC=$(which mpicc) && export CXX=$(which mpicxx) && export FC=$(which mpif90) && module load slurm'

    print("Platform machine: " + platform.machine())
    print("Platform system: " + platform.system())
    print("System path: " + str(sys.path))

    dirPath = os.path.dirname(os.path.abspath(__file__)) + "/Feat_FloWer"
    dirPathBin = os.path.dirname(os.path.abspath(__file__)) + "/bin"
    #sys.exit(2)
    
    subprocess.call(['ctest -S ctest_driver.cmake -DBUILD_STRING=skylake-linux-gcc-release -DPROCS=5 -DSRC_DIR=%s -DBIN_DIR=%s -VV' %(dirPath, dirPathBin)], shell=True)

    # Remove the note json files from the binary dir
    potentialFiles = os.listdir(dirPathBin)

    for item in potentialFiles:
        m = re.search(r"note_single_\w+-bench.json$", item)
        if m:
            os.remove(dirPathBin + "/" + item)

#===============================================================================
#                           Main "Boiler Plate"
#===============================================================================
if __name__ == "__main__":
    main()
