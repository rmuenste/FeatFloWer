#!/usr/bin/env python
# vim: set filetype=python
import sys
import getopt
import platform
sys.path.append('/home/user/rmuenste/bin/partitioner')
import PyPartitioner
import subprocess
import re
import json

################
def usage():
  print("Usage: configure [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-h, --help]: prints this message")
################

def get_log_variable(file_name, var_name):
  with open(file_name, "r") as sources:
    lines = sources.readlines()
  
  it_found = False
  t_found = False

  for line in reversed(lines):
    m = re.match(" Force acting on the cylinder:",line)
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

module_string = "module purge && module load gcc/6.1.0 openmpi/gcc6.1.x/1.10.2/non-threaded/no-cuda/ethernet cmake && export CC=mpicc && export CXX=mpicxx && export FC=mpif90"

PyPartitioner.MainProcess(4, 1, 1, "NEWFAC", "_adc/ViscoHex2/aaa.prj")
subprocess.call(['mpirun -np 5 ./q2p1_fac_visco'],shell=True)
force = get_log_variable("_data/prot.txt", " Force acting on the cylinder:")
force = force.split()
d = {'ID' : 'VISCO-FAC', 'Caption' : 'Visco-Elastic Flow Around A Cylinder', 
'Drag': force[1], 'Lift' : force[2]}

print(str(json.dumps(d)))
with open('note.json','w') as f:
  f.write('[' + json.dumps(d) + ']\n')
with open('note_single.json','w') as f:
  f.write(json.dumps(d) + '\n')

