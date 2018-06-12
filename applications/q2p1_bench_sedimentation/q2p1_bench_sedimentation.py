#!/usr/bin/env python
# vim: set filetype=python
import sys
import getopt
import platform
import subprocess
import re
import json
sys.path.append('/home/user/rmuenste/bin/partitioner')
import part_main

################
def usage():
  print("Usage: configure [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-h, --help]: prints this message")

###################################################################      

def get_col_data(file_name):
  with open(file_name, "r") as sources:
    lines = sources.readlines()
  
  it_found = False
  t_found = False
  
  col_data = []
  for line in lines:
    m = re.search("PartVel:",line)
    if m != None:
      val = line.split()
      col_data.append([val[4],val[3]])

  return col_data

###################################################################      

def write_json_data(col):

  rows_array = [] 
  for i in range(len(col)-1):
    rows_array.append({"c": [{"v" : col[i][0]}, {"v": col[i][1]}] })

  d = {
   "ID" : "BENCHSED", 
   "Caption" : "Sedimentation Benchmark", 
   "data" : {
     "cols": [
     {"label" : "Time", "type" : "number"},
     {"label" : "U_z", "type" : "number"}
     ],
     "rows" : rows_array
   }
  }

  with open("note_single.json","w") as f:
    json.dump(d,f)
    f.write("\n")

#  with open("note_single_alt.json","r") as f:
#    parsed = json.load(f)
#    print(json.dumps(parsed, indent=2, sort_keys=True))

#  with open("note_single.json","w") as f:
#    f.write("{\n")
#    f.write('"ID": "BENCHSED", "Caption": "Sedimentation Benchmark", "data": \n')
#    ###############
#    f.write("{\n")
#    f.write('"cols": [\n')
#    f.write('{"label": "Time", "type": "number"},\n')
#    f.write('{"label": "U_z", "type": "number"}\n')
#    f.write("],\n")
#    f.write('"rows": [\n')
#    for i in range(len(col)-1):
#        f.write('{"c":[{"v":' + col[i][0] + '},{"v":' + col[i][1] + '}]},\n')
#
#    f.write('{"c":[{"v":' + col[len(col)-1][0] + '},{"v":' + col[len(col)-1][1] + '}]}\n')
#    f.write("]\n")
#    f.write("}\n")
#    ###############
#    f.write("}\n")



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

part_main.MainProcess(4, 1, 1, "NEWFAC", "_adc/benchSym/bench.prj")
subprocess.call(['mpirun -np 5 ./q2p1_bench_sedimentation'],shell=True)
data = get_col_data("_data/prot.txt")
write_json_data(data)

