#!/usr/bin/env python
# vim: set filetype=python
import os
import shutil

import sys
import getopt
import platform
import subprocess
import re
import json
sys.path.append(os.environ['FF_PY_HOME'])
import partitioner

################
def usage():
  print("Usage: configure [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-n, --num-processors]: Total number of MPI processes to be used in the simulation")
################
def write_json_data(col, level):

  rows_array = [] 
  for i in range(len(col)-1):
    rows_array.append({"c": [{"v" : level}, {"v" : col[i][0]}, {"v": col[i][1]}] })

  d = {
   "ID" : "BENCHSED", 
   "Caption" : "Sedimentation Benchmark", 
   "data" : {
     "cols": [
     {"label" : "Level", "type" : "number"},
     {"label" : "Drag", "type" : "number"},
     {"label" : "Lift", "type" : "number"}
     ],
     "rows" : rows_array
   }
  }
###################################################################      
def moveAndSetLevel(file_in, file_out, level):
    maxLevelStr = "SimPar@MaxMeshLevel = " + str(level) 
    with open(file_out, "w") as n:
        with open(file_in, "r") as f:
            for line in f:
                new_line = re.sub(r"^[\s]*SimPar@MaxMeshLevel[\s]*=(\s | \w)*", maxLevelStr, line)
                n.write(new_line)

#===============================================================================
#                              Main function
#===============================================================================
def get_log_entry(file_name, var_name):
  with open(file_name, "r") as sources:
    lines = sources.readlines()
  
  it_found = False
  t_found = False

  for line in reversed(lines):
    m = re.match(var_name,line)
    if m != None:
      mysplit = line.split(':')
      val = mysplit[1].strip()
      t_found = True
      break

  if t_found:
    return val
  else :
    return 0

#===============================================================================
#                              Main function
#===============================================================================
def main():
  try:
      opts, args = getopt.getopt(sys.argv[1:], 'n:p:h', ['num-processors=', 'params=', 'help'])
  except getopt.GetoptError:
      usage()
      sys.exit(2)
  
  params=''

  numProcessors=16
  
  for opt, arg in opts:
      if opt in ('-h', '--help'):
          usage()
          sys.exit(2)
      elif opt in ('-n', '--num-processors'):
          numProcessors = int(arg)
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
  
  rows_array = [] 
  
  for l in range(2,4):
    moveAndSetLevel("_adc/FAC3Ds/q2p1_param_FAC3D_stat_CC.dat", "_data/q2p1_param.dat",l)
    partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", "_adc/FAC3Ds/fac3D_stat_Re20.prj")
    subprocess.call(['mpirun -np %i ./q2p1_cc' %numProcessors],shell=True)
    force = get_log_entry("_data/prot.txt", "BenchForce:")
    force = force.split()
    rows_array.append({"c": [{"v" : l}, {"v" : force[0]}, {"v": force[1]} ] })
  
  d = {
   "benchName" : "CC3DSTAT", 
   "tableCaption" : "Re20 Stationary FAC", 
   "style" : "Table",
   "data" : {
     "cols": [
     {"label" : "Level", "type" : "number"},
     {"label" : "Drag", "type" : "number"},
     {"label" : "Lift", "type" : "number"}
     ],
     "rows" : rows_array
   }
  }
  
  #print(str(json.dumps(rows_array)))
  with open('../../note_single_Re20FAC3DCC-bench.json','w') as f:
    f.write(json.dumps(d) + '\n')

#===============================================================================
#                           Main "Boiler Plate"
#===============================================================================
if __name__ == "__main__":
  main()

