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

PyPartitioner.MainProcess(4, 1, 1, "NEWFAC", "_adc/2D_FAC/2Dbench.prj")
subprocess.call(['mpirun -np 5 ./q2p1_fc_ext'],shell=True)
force = get_log_entry("_data/prot.txt", " Force acting on the cylinder:")
force = force.split()
d = {'ID' : 'NEWTFAC', 'Caption' : 'Newtonian Flow Around A Cylinder', 
'Drag': force[1], 'Lift' : force[2]}

print(str(json.dumps(d)))
with open('note.json','w') as f:
  f.write('[' + json.dumps(d) + ']\n')
with open('note_single.json','w') as f:
  f.write(json.dumps(d) + '\n')

