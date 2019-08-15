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
import partitioner

if sys.version_info[0] < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

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
    print("[-a', '--angle]: If this parameter is present a single simulation with " +
          "that specific angle will be done.")
    print("[-d', '--delta-angle]: The angular step size between two simulations " +
          "the sim loop (default 10)")
    print("[-c', '--host-conf]: A hostfile as input for the mpirun command")
    print("[-r', '--rank-file]: A rankfile as input for the mpirun command")
    print("[-t', '--time]: Number of time levels to complete a full 360 rotation")
    print("[-h', '--help']: prints this message")
    print("Example: python ./e3d_start.py -f myFolder -n 5 -t 0")

#===============================================================================
#                          setup the case folder 
#===============================================================================
def folderSetup(workingDir, projectFile, projectPath, projectFolder):
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

    for item in offList:
        shutil.copyfile(str(item), str(workingDir / item.name))

    if Path("_data/meshDir").exists():
        print("meshDir exists")
        shutil.rmtree("_data/meshDir")

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

    # Single angle 
    singleAngle = -10.0 

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

    skipSetup = False

    skipSimulation = False

    # Number of simulation for a full rotation 
    nmax = 0

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:p:d:a:c:r:t:smh',
                                   ['num-processors=', 'project-folder=',
                                    'periodicity=', 'delta-angle=', 'angle=',
                                    'host-conf=', 'rank-file=', 'time=', 'skip-setup',
                                    'skip-simulation', 'help'])

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
            singleAngle = float(arg)
        elif opt in ('-d', '--delta-angle'):
            deltaAngle = float(arg)
        elif opt in ('-c', '--host-conf'):
            hostFile = arg
        elif opt in ('-r', '--rank-file'):
            rankFile = arg
        elif opt in ('-t', '--time'):
            timeLevels = int(arg)
        elif opt in ('-s', '--skip-setup'):
            skipSetup = True
        elif opt in ('-m', '--skip-simulation'):
            skipSimulation = True
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

    if not skipSetup:
        folderSetup(workingDir, projectFile, projectPath, projectFolder)
        subprocess.call(["./s3d_mesher"])
        partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", "_data/meshDir/file.prj")

    if skipSimulation:
        sys.exit()

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

    if singleAngle >= 0.0:
        nmax = 1

    mpiPath = Path("mpirun")
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
        if singleAngle >= 0.0:
            angle = singleAngle
        else:
            angle = start + i * deltaAngle

        shutil.copyfile("_data/Extrud3D_0.dat", "_data/Extrud3D.dat")

        with open("_data/Extrud3D.dat", "a") as f:
            f.write("Angle=" + str(angle) + "\n")

        if sys.platform == "win32":
            subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse.exe"])
        else:
            #comm = subprocess.call(['mpirun', '-np', '%i' % numProcessors,  './q2p1_sse', '-a', '%d' % angle],shell=True)
            subprocess.call(['mpirun -np %i ./q2p1_sse -a %d' % (numProcessors, angle)],shell=True)

        iangle = int(angle)

#        if not singleAngle >= 0.0:
#            shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

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
