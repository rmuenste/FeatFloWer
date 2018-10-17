#!/usr/bin/env python
# vim: set filetype=python
"""
A python launcher script for the FeatFloWer SSE application
"""
import sys
import getopt
import platform
import os
import shutil
import subprocess
import re
import json
sys.path.append(os.environ['FF_PY_HOME'])
import partitioner

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def usage():
  print("Usage: python q2p1_sse_start.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-f, --case-folder]: Configuration folder for the particular SSE case")
  print("[-n, --num-processors]: Total number of MPI processes to be used in the simulation")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def inspectCaseFolder(path):
  """ This routine examines the case folder and tries to
      find a file named Extrud3D.dat that is needed for the simulation
      to execute properly. If this file is missing an error is returned

      Args:
          path: The string path to the particular case folder 
  """
  dir_list = os.listdir(path)
  print("Case folder: " + str(dir_list))
  prj_found = False

  project_file = ''

  for item in dir_list:
    m = re.search(r"^Extrud3D.dat$", item)
    if m: 
      prj_found = True
      project_file = item

  if not prj_found:
    print("Error: The case folder %s does not contain an Extrud3D.dat project file. " % path)
    sys.exit(2)
  else:
    return path + "/" + project_file
  
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#                              Main function
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
def main():
  """ This function is the entry point of the driver script. It first checks the
      inputs for general validity and then tries to start the parallel simulation.

      =====================
      Simulation
      =====================
      The SSE simulation is split up in several parallel simulations that run in a sequence.
      Each particular simulation in the sequence computes the configuration for a different
      angle. The sequence of computations is done in a loop in this routine.

      =====================
      Configuration Files
      =====================
      The main configuration file is the <Extrud3D.dat> file that has to be present in the
      case folder. In the initialization of the simulation the <Extrud3D.dat> is copied 
      to <_data/Extrud3D_0.dat>. This file is then copied to <_data/Extrud3D.dat> and a
      line that sets the current angle is added to it, i.e.:
      [E3DSimulationSettings]
      Angle=40 
      
      In each iteration of the angle loop this procedure is done and the <Angle> variable
      is set to the current value.

      =====================
      Meshing
      =====================
      The SSE application will use the mesh folder under the path <_data/meshDir>.
      In the initialization phase any <_data/meshDir> folder that may be present from a previous run
      is deleted. 
      This driver script will then call the external s3d_mesher program. The s3d_mesher will either
      generate a mesh if the user specified automatic mesh generation in the <Extrud3D.dat> file.
      If automatic mesh generation by the s3d_mesher was successful a mesh folder will be present
      under the path <_data/meshDir>. 
      If automatic mesh generation is not used this script will try to locate a <meshDir> folder
      inside of the case folder and copy it to <_data/meshDir>.
  """

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

  # Check the case folder for a Extrud3D.dat file
  # and returns it if found
  prj_file = inspectCaseFolder(caseFolder)

  shutil.copyfile(prj_file, "_data/Extrud3D_0.dat")

  delta = 40
  nmin = 0
  #nmax = nmin + 1
  nmax = 1
  start = 0.0

  with open("_data/Extrud3D_0.dat", "a") as f:
    f.write("\n[E3DSimulationSettings]\n") 
    f.write("dAlpha=" + str(delta) + "\n") 

  if os.path.exists("_data/meshDir"):
    shutil.rmtree("_data/meshDir")

  subprocess.call(["./s3d_mesher"])

  if not os.path.exists("_data/meshDir"):
    if os.path.exists(caseFolder + "/meshDir"):
      shutil.copytree(caseFolder + "/meshDir","_data/meshDir")
    else:
      print("Error: No mesh automatically generated and no <meshDir> " + 
            "folder present the case folder " + caseFolder)
      sys.exit(2)

  #-------------------------------------------------------------------------------------
  numProcessors = int(numProcessors)

  partitioner.mkdir("_mesh")

  partitioner.MainProcess(numProcessors-1, -3, 2, "NEWFAC", "_data/meshDir/file.prj")
  #-------------------------------------------------------------------------------------

  for i in range(nmin,nmax):  # nmax means the loop goes to nmax-1
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

#===============================================================================
#                           Main "Boiler Plate"
#===============================================================================
if __name__ == '__main__':
  main()
