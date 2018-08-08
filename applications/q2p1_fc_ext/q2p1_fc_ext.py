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
        opts, args = getopt.getopt(sys.argv[1:], 'n:p:h', ['num-processors=', 'params=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    params = ''
  
    numProcessors = 16
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-n', '--num-processors'):
            if arg.isdigit():
                numProcessors = int(arg)
            else:
                print("--num-processors= " + str(arg) + " is not a valid number.")
                sys.exit(2)
        elif opt in ('-p', '--params'):
            params = arg
        else:
            usage()
            sys.exit(2)
    
    if params != '':
        print("Parameter params = " + params)
    
    print("Platform machine: " + platform.machine())
    print("Platform system: " + platform.system())
    print("System path: " + str(sys.path))
  
  
    if numProcessors == 0:
        print("Number of processors is 0")
        usage()
        sys.exit(2)

    rowsArray = [] 
    timeEntry = [] 
    
    partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", "_adc/2D_FAC/2Dbench.prj")
    for ilevel in range(2, 4):
        moveAndSetLevel("_adc/2D_FAC/q2p1_param_2D.dat", "_data/q2p1_param.dat", ilevel)
        subprocess.call(['mpirun -np %i ./q2p1_fc_ext' %numProcessors], shell=True)
        force = getLogEntry("_data/prot.txt", "BenchForce:")
        force = force.split()
        timeEntry = getLogEntry("_data/Statistics.txt", " Overall time")
        timeEntry = timeEntry.split()
        rowsArray.append({"c": [{"v" : ilevel}, {"v" : force[1]}, {"v": force[2]}, {"v": timeEntry[0][:-3]}]})

    jsonData = {
        "benchName" : "NEWTFAC", 
        "tableCaption" : "Newtonian Flow Around A Cylinder", 
        "style" : "Table",
        "data" : {
            "cols": [
                {"label" : "Level", "type" : "number"},
                {"label" : "Drag", "type" : "number"},
                {"label" : "Lift", "type" : "number"},
                {"label" : "Time[s]", "type" : "number"}],
            "rows" : rowsArray}}
    
    #print(str(json.dumps(rows_array)))
    with open('../../note_single_fac2D-bench.json', 'w') as theFile:
        theFile.write(json.dumps(jsonData) + '\n')


#===============================================================================
#                           Main "Boiler Plate"
#===============================================================================
if __name__ == "__main__":
    main()
