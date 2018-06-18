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

shutil.copyfile("_adc/FAC3Ds/q2p1_param_nonstat.dat", "_data/q2p1_param.dat")
partitioner.partition(4, 1, 1, "NEWFAC", "_adc/FAC3Ds/fac3D_nonstat.prj")
subprocess.call(['mpirun -np 5 ./q2p1_fc_ext'],shell=True)
force = get_log_entry("_data/prot.txt", "BenchForce:")
force = force.split()
d = {'ID' : '3DFAC', 'Caption' : 'Full 3D Newtonian FAC', 
'Drag': force[1], 'Lift' : force[2]}

print(str(json.dumps(d)))
with open('../../note_single_fac3D-bench.json','w') as f:
  f.write(json.dumps(d) + '\n')

