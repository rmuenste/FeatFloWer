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

###################################################################      
def moveAndSetLevel(file_in, file_out, level):
    maxLevelStr = "SimPar@MaxMeshLevel = " + str(level) 
    with open(file_out, "w") as n:
        with open(file_in, "r") as f:
            for line in f:
                new_line = re.sub(r"^[\s]*SimPar@MaxMeshLevel[\s]*=(\s | \w)*", maxLevelStr, line)
                n.write(new_line)

###################################################################      

# Script begin

try:
    opts, args = getopt.getopt(sys.argv[1:], 'm:p:h', ['miner=', 'params=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

params=''
miner_name=''

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-m', '--miner'):
        miner_name = arg
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

rows_array = [] 

for l in range(2,4):
  moveAndSetLevel("_adc/2D_FAC/q2p1_param_nnewt_2D.dat", "_data/q2p1_param.dat",l)
  partitioner.partition(4, 1, 1, "NEWFAC", "_adc/2D_FAC/2Dbench.prj")
  subprocess.call(['mpirun -np 5 ./q2p1_fac_nnewt'],shell=True)
  force = get_log_entry("_data/prot.txt", "BenchForce:")
  force = force.split()
  rows_array.append({"c": [{"v" : l}, {"v" : force[1]}, {"v": force[2]} ] })
  

d = {
 "benchName" : "NON-NEWTFAC", 
 "tableCaption" : "Non-Newtonian Flow Around A Cylinder", 
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

print(str(json.dumps(d)))
with open('../../note_single_nnewt_fac2D-bench.json','w') as f:
  f.write(json.dumps(d) + '\n')

