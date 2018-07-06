#!/usr/bin/env python
# vim: set filetype=python
import os

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
      sim_time = mysplit[1].strip()
      #print("Force = " + str(sim_time))
      t_found = True
      break

  if t_found:
    return sim_time
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
  
  partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", "_adc/ViscoHex2/aaa.prj")
  for l in range(2,4):
    moveAndSetLevel("_adc/ViscoHex2/q2p1_param_visco_2D.dat", "_data/q2p1_param.dat",l)
    subprocess.call(['mpirun -np %i ./q2p1_fac_visco' %numProcessors],shell=True)
    force = get_log_entry("_data/prot.txt", "TimevsForce")
    print(str(force))
    force = force.split()
    timeEntry = get_log_entry("_data/Statistics.txt", " Overall time")
    timeEntry = timeEntry.split()
    rows_array.append({"c": [{"v" : l}, {"v" : force[1]}, {"v": force[2]}, {"v": timeEntry[0][:-3]} ] })
    
  
  d = {
   "benchName" : "VISCO-FAS", 
   "tableCaption" : "Visco-Elastic Flow Around A Cylinder", 
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
  
  with open('../../note_single_visco_fac2D-bench.json','w') as f:
    f.write(json.dumps(d) + '\n')

if __name__ == "__main__":
  main()
