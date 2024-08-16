#!/usr/bin/env python
# vim: set filetype=python
"""
A python launcher script for a FeatFloWer application
"""
import os
import xml.etree.ElementTree as ET
import sys
import getopt
import platform
import subprocess
import re
import json
sys.path.append(os.environ['FF_PY_HOME'])
import partitioner
import pprint

from tempfile import mkstemp
from shutil import move
from os import fdopen, remove

#===============================================================================
#                      Function: Usage
#===============================================================================
def usage():
    print("Usage: q2p1_fc_ext_start.py [options]")
    print("Where options can be:")
    print("[-h, --help]: prints this message")
    print("[-n, --num-processor]: defines the number of parallel jobs to be used")
    print("[-s, --subparts]: defines the number of subpartitions for recursive partitioning)")
    print("[-i, --in-file]: defines the input project file (will be replaced in the q2p1-data-file)")
    print("[-r, --rankfile]: defines the rank file for MPI")
    print("[-m, --machinefile]: defines the machine file for MPI")

def replace(file_path, pattern, subst):
    # Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh, 'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                if pattern in line:
                    new_file.write(line.replace(line, subst + "\n"))
                else:
                    new_file.write(line)
    # Remove original file
    remove(file_path)
    # Move new file
    move(abs_path, file_path)

def main(argv):
    inputfile = ''
    NumProcessor = 1
    subparts = 1
    rankfile = ''
    machinefile = ''
    
    try:
        opts, args = getopt.getopt(argv, "hi:n:s:r:m:", ["in-file=", "num-processor=", "subparts=", "rankfile=", "machinefile="])
    except getopt.GetoptError:
        sys.exit(2)
    
    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("-n", "--num-processor"):
            NumProcessor = int(arg)
        elif opt in ("-s", "--subparts"):
            subparts = int(arg)
        elif opt in ("-r", "--rankfile"):
            rankfile = "--rankfile " + arg
        elif opt in ("-m", "--machinefile"):
            machinefile = "--machinefile " + arg
        elif opt in ("-i", "--in-file"):
            inputfile = arg

    print('Number Of Processors is: ' + str(NumProcessor))
    print('Number Of Subpartitions is: ' + str(subparts))
    print('Input file is: ' + inputfile)
    if rankfile:
        print('Rankfile is: ' + rankfile)
    if machinefile:
        print('Machinefile is: ' + machinefile)
    
    replace("_data/q2p1_param.dat", "SimPar@ProjectFile = ", "SimPar@ProjectFile = '" + inputfile + "'")
    replace("_data/q2p1_param.dat", "SimPar@SubMeshNumber = ", "SimPar@SubMeshNumber = " + str(subparts))
    
    partitioner.partition(NumProcessor - 1, 1, subparts, "NEWFAC", str(inputfile))
    
    # Construct the mpirun command, include rankfile and machinefile if they are provided
    mpirun_command = f'mpirun -np {NumProcessor} {rankfile} {machinefile} ./q2p1_fc_ext'
    
    # Run the command
    print("Executing command:", mpirun_command)
    subprocess.call(mpirun_command, shell=True)

if __name__ == "__main__":
    main(sys.argv[1:])
