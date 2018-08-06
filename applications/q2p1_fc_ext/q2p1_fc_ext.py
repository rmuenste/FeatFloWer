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
  print("[-h, --help]: prints this message")
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

###################################################################      

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
  
  params=''

  numProcessors=16
  
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
  
  rows_array = [] 
  timeEntry = [] 
  
  partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", "_adc/2D_FAC/2Dbench.prj")
  for l in range(2,4):
    moveAndSetLevel("_adc/2D_FAC/q2p1_param_2D.dat", "_data/q2p1_param.dat",l)
    subprocess.call(['mpirun -np %i ./q2p1_fc_ext' %numProcessors],shell=True)
    force = get_log_entry("_data/prot.txt", "BenchForce:")
    force = force.split()
    timeEntry = get_log_entry("_data/Statistics.txt", " Overall time")
    timeEntry = timeEntry.split()
    rows_array.append({"c": [{"v" : l}, {"v" : force[1]}, {"v": force[2]}, {"v": timeEntry[0][:-3]} ] })


  d = {
   "benchName" : "NEWTFAC", 
   "tableCaption" : "Newtonian Flow Around A Cylinder", 
   "style" : "Table",
   "data" : {
     "cols": [
     {"label" : "Level", "type" : "number"},
     {"label" : "Drag", "type" : "number"},
     {"label" : "Lift", "type" : "number"},
     {"label" : "Time[s]", "type" : "number"}
     ],
     "rows" : rows_array
   }
  }
  
  #print(str(json.dumps(rows_array)))
  with open('../../note_single_fac2D-bench.json','w') as f:
    f.write(json.dumps(d) + '\n')

#===============================================================================
#                           Main "Boiler Plate"
#===============================================================================
if __name__ == "__main__":
  main()
