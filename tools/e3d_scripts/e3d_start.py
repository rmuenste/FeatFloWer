#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver for the e3d simulator

"""

import os
import sys
import re
import getopt
import shutil
import subprocess
import math
try:
    sys.path.append(os.environ['FF_PY_HOME'])
except:
    print("Warning: Environment variable <FF_PY_HOME> is not set. "  +
          "Trying to use local install of FeatFloWer python modules.")
import partitioner

from pathlib2 import Path

#===============================================================================
#                        usage function
#===============================================================================
def usage():
    """
    Print out usage information
    """
    print("Usage: python e3d_start.py [options]")
    print("Where options can be:")
    print("[-f', '--project-folder]: Path to project folder containing a setup.e3d file")
    print("[-n', '--num-processors]: Number of processors to use")
    print("[-p', '--periodicity']: Periodicity of the solution (1, 2, 3, ... " +
          "usually the time flight number)")
    print("[-a', '--angle]: The angular step size between two simulations in " +
          "the sim loop (default 10)")
    print("[-c', '--host-conf]: A hostfile as input for the mpirun command")
    print("[-r', '--rank-file]: A rankfile as input for the mpirun command")
    print("[-t', '--time]: Number of time levels to complete a full 360 rotation")
    print("[-h', '--help']: prints this message")
    print("Example: python ./e3d_start.py -f myFolder -n 5 -t 0")


#===============================================================================
#                        Main Script Function
#===============================================================================
def main():
    """
    The main function that controls the extrusion process

    Options:
        project-folder: Path to project folder containing a setup.e3d file 
        num-processors: Number of processors to use
        periodicity: Periodicity of the solution (1, 2, 3, ... usually the time flight number)
        angle: The angular step size between two simulations in the sim loop (default 10)
        host-conf: A hostfile as input for the mpirun command
        rank-file: A rankfile as input for the mpirun command
        time: Number of time levels to complete a full 360 rotation 
    """
    # Angular step size
    deltaAngle = 10.0

    # Hostfile 
    hostFile = ""

    # Rankfile 
    rankFile = ""

    # timeLevels 
    timeLevels = 0

    # Periodicity 
    periodicity = 1

    # Number of processors 
    numProcessors = 5

    # The project folder
    projectFolder = ""

    # Number of simulation for a full rotation 
    nmax = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:p:a:c:r:t:h',
                                   ['num-processors=', 'project-folder=',
                                    'periodicity=', 'angle=',
                                    'host-conf=', 'rank-file=', 'time=', 'help'])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-f', '--project-folder'):
            projectFolder = arg
        elif opt in ('-n', '--num-processors'):
            numProcessors = int(arg)
        elif opt in ('-p', '--periodicity'):
            periodicity = int(arg)
        elif opt in ('-a', '--angle'):
            deltaAngle = float(arg)
        elif opt in ('-c', '--host-conf'):
            hostFile = arg
        elif opt in ('-r', '--rank-file'):
            rankFile = arg
        elif opt in ('-t', '--time'):
            timeLevels = int(arg)
        else:
            usage()
            sys.exit(2)

    if projectFolder == "":
        print("Error: no project folder specified.")
        usage()
        sys.exit(2)

    workingDir = Path('.')
    projectPath = Path(workingDir / projectFolder)
    projectFile = Path(projectPath / 'setup.e3d')

    if not projectFile.is_file():
        projectFile = Path(projectPath / 'Extrud3D.dat')
        if not projectFile.is_file():
            print("Could not find a valid parameter file in the project folder: " +
                  str(projectPath))
            sys.exit(2)

    backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU.dat")
    destDataFile = Path("_data") / Path("q2p1_param.dat")

    shutil.copyfile(str(backupDataFile), str(destDataFile))
    shutil.copyfile(str(projectFile), str(workingDir / Path("_data/Extrud3D_0.dat")))

    offList = []
    if not sys.platform == "win32":
        offList = list(Path(projectFolder).glob('*.off')) + list(Path(projectFolder).glob('*.OFF'))
    else: 
        offList = list(Path(projectFolder).glob('*.off'))

    print(offList)
    print(projectFolder)

    for item in offList:
        shutil.copyfile(str(item), str(workingDir / item.name))

    if Path("_data/meshDir").exists():
        print("meshDir exists")
        shutil.rmtree("_data/meshDir")

    subprocess.call(["./s3d_mesher"])

    partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", "_data/meshDir/file.prj")

    print("  ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("                                        Launching E3D...                                           ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("  ")
    print("  ")

    if not deltaAngle == 0.0 and timeLevels in (0, 1):
        nmax = 1
    elif not deltaAngle == 0.0:
        nmax = int(math.ceil(360.0 / periodicity / deltaAngle))
    else:
        deltaAngle = 360.0 / float(timeLevels)
        nmax = int(math.ceil(360.0 / periodicity / deltaAngle))

    mpiPath =Path("mpirun")
    if sys.platform == "win32":
        mpiPath = Path(os.environ['MSMPI_BIN']) / Path("mpiexec.exe")

    nmin = 0
    start = 0.0
    with open("_data/Extrud3D_0.dat", "a") as f:
        f.write("\n[E3DSimulationSettings]\n")
        f.write("dAlpha=" + str(deltaAngle) + "\n")
        f.write("Periodicity=" + str(periodicity) + "\n")
        f.write("nSolutions=" + str(timeLevels) + "\n")

    for i in range(nmin, nmax):  # nmax means the loop goes to nmax-1
        angle = start + i * deltaAngle
        shutil.copyfile("_data/Extrud3D_0.dat", "_data/Extrud3D.dat")

        with open("_data/Extrud3D.dat", "a") as f:
            f.write("Angle=" + str(angle) + "\n")

        if sys.platform == "win32":
            subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse.exe"])
        else:
            subprocess.call(['mpirun -np %i ./q2p1_sse' % numProcessors],shell=True)

        iangle = int(angle)
        shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

    if not sys.platform == "win32":
        offList = list(workingDir.glob('*.off')) + list(workingDir.glob('*.OFF'))
    else: 
        offList = list(workingDir.glob('*.off'))

    for item in offList:
        os.remove(str(item))


#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
