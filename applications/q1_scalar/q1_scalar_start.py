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
    print("Usage: q1_scalar_start.py [options]")
    print("Where options can be:")
    print("[-h, --help]: prints this message")
    print("[-n, --num-processor]: defines the number of parallel jobs to be used")
    print("[-i, --in-file]: defines the input project file (will be replaced in the q2p1-data-file)")

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
   inputfile = ''
   outputfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:n:",["in-file=","num-processor="])
   except getopt.GetoptError:
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
 #        print 'test.py -i <inputfile> -o <outputfile>'
         usage()
         sys.exit()
      elif opt in ("-n", "--num-processor"):
         NumProcessor = int(arg)
      elif opt in ("-i", "--in-file"):
         inputfile = arg

   print('Number Of Processors is: '+ str(NumProcessor))
   print('Input file is '+ inputfile)
   
   replace("_data/q2p1_param.dat","SimPar@ProjectFile = ","SimPar@ProjectFile = '" + inputfile + "'")
   
   partitioner.partition(NumProcessor-1, 1, 1, "NEWFAC", str(inputfile))
        
   subprocess.call(['mpirun -np %i ./%s' %(NumProcessor, 'q1_scalar')], shell=True)
   
if __name__ == "__main__":
   main(sys.argv[1:])