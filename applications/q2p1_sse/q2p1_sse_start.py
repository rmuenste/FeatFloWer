#!/usr/bin/env python
# vim: set filetype=python
import sys
import getopt
import platform
import os
import shutil
import subprocess
import re
import json
import part_main

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def usage():
  print("Usage: python q2p1_sse_start.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-f, --case-folder]: Configuration folder for the particular SSE case")
  print("[-n, --num-processors]: Total number of MPI processes to be used in the simulation")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

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


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def inspectCaseFolder(path):

  dir_list = os.listdir(path)
  print("Case folder: " + str(dir_list))
  prj_found = False

  project_file = ''

  for item in dir_list:
    m = re.search(r"Extrud3D.dat", item)
    if m: 
      prj_found = True
      project_file = item

  if not prj_found:
    print("Error: The case folder %s does not contain an Extrud3D.dat project file. " % path)
    sys.exit(2)
  else:
    return path + "/" + project_file
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
# Main script begin
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

try:
    opts, args = getopt.getopt(sys.argv[1:], 'f:n:h', ['case-folder=', 'num-processors=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

caseFolder=''
numProcessors=0

if len(opts) < 2:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-f', '--case-folder'):
        caseFolder = arg
    elif opt in ('-n', '--num-processors'):
        numProcessors = int(arg)
    else:
        usage()
        sys.exit(2)

print("Platform machine: " + platform.machine())
print("Platform system: " + platform.system())
print("Platform: " + sys.platform)

prj_file = inspectCaseFolder(caseFolder)
shutil.copyfile(prj_file, "_data/Extrud3D_0.dat")

delta = 40
nmin = 0
#nmax = nmin + 1
nmax = 2
start = 0.0

with open("_data/Extrud3D_0.dat", "a") as f:
  f.write("[E3DSimulationSettings]\n") 
  f.write("dAlpha=" + str(delta) + "\n") 

subprocess.call(["./s3d_mesher"])

if not os.path.exists("_data/meshDir"):
  shutil.copytree(caseFolder + "/meshDir","_data/meshDir")

#-------------------------------------------------------------------------------------
numProcessors = int(numProcessors)

part_main.mkdir("_mesh")

part_main.MainProcess(numProcessors-1, 1, 1, "NEWFAC", "_data/meshDir/file.prj")
#-------------------------------------------------------------------------------------

for i in range(nmin,nmax):
  angle = start + i * delta
  shutil.copyfile("_data/Extrud3D_0.dat","_data/Extrud3D.dat")
  with open("_data/Extrud3D.dat", "a") as f:
    f.write("Angle=" + str(angle) + "\n") 

  if sys.platform == "win32":
    shutil.copyfile("Release/q2p1_sse.exe", "./q2p1_sse.exe")
    print("Computing with angle: %f" % angle)
    print("Command:mpiexec -n %i ./q2p1_sse.exe" % numProcessors)
    subprocess.call(["mpiexec", "-n", str(numProcessors), "./q2p1_sse.exe"])
    iangle = int(angle)
    shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)
  else:
    subprocess.call(['mpirun -np %i ./q2p1_sse' % numProcessors],shell=True)
    iangle = int(angle)
    shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

