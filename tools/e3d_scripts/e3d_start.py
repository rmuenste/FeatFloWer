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

import datetime

if sys.version_info[0] < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

paramDict = {
    "deltaAngle": 10.0, # Angular step size
    "singleAngle": -10.0, # Single angle to compute 
    "hostFile" : "", # Hostfile
    "rankFile" : "" , # Rankfile 
    "timeLevels" : 36, # timeLevels
    "periodicity" : 1, # Periodicity 
    "numProcessors" : 5, # Number of processors 
    "projectFolder" : "", # The project folder
    "skipSetup" :  False,
    "skipSimulation" : False,
    "hasDeltaAngle": False,
    "hasTimeLevels": False,
    "temperature" : False
}
#===============================================================================
#                        version function
#===============================================================================
def version():
    """
    Print out version information
    """
    print("E3D + Reporter for SIGMA Version 19.10, Copyright 2019 IANUS Simulation")

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
    print("['-o', '--do-temperature']: The simulation loop will do a temperature simulation")
    print("[-h', '--help']: prints this message")
    print("[-v', '--version']: prints out version information")
    print("Example: python ./e3d_start.py -f myFolder -n 5 -t 0")

#===============================================================================
#                          custom mkdir
#===============================================================================
def mkdir(dir):
    if os.path.exists(dir):
        if os.path.isdir(dir):
            return
        else:
            os.remove(dir)
    os.mkdir(dir)

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
#                             Simulation Setup 
#===============================================================================
def simulationSetup(workingDir, projectFile, projectPath, projectFolder):
    folderSetup(workingDir, projectFile, projectPath, projectFolder)

    subprocess.call(["./s3d_mesher"], env={"LD_LIBRARY_PATH": os.environ['LD_LIBRARY_PATH'] + ':.', "PATH":os.environ['PATH']})
    #subprocess.call(["./s3d_mesher"])
    
    if not Path("_data/meshDir").exists():    
      meshDirPath = projectPath / Path("meshDir")
      if meshDirPath.exists():
         print('Copying meshDir from Project Folder!')
         shutil.copytree(str(meshDirPath), "_data/meshDir")
      else:
         print("Error: No mesh automatically generated and no <meshDir> " + 
               "folder present the case folder " + str(projectPath))
         sys.exit(2)    
    
    partitioner.partition(paramDict['numProcessors']-1, 1, 1, "NEWFAC", "_data/meshDir/file.prj")
#===============================================================================
#                Compute maximum number of simulation iterations
#===============================================================================
def calcMaxSimIterations():
    nmax = 0

    if paramDict['hasDeltaAngle']:
        paramDict['timeLevels'] = 360.0 / paramDict['deltaAngle'] 

    if paramDict['timeLevels'] == 1:
        nmax = 1
    else:
        paramDict['deltaAngle'] = 360.0 / float(paramDict['timeLevels'])
        nmax = int(math.ceil(360.0 / paramDict['periodicity'] / paramDict['deltaAngle']))

    if paramDict['singleAngle'] >= 0.0:
        nmax = 1
    
    print("nmax: ",nmax)

    return nmax

#===============================================================================
#                Compute maximum number of simulation iterations
#===============================================================================
def setupMPICommand():
    mpiPath = Path("mpirun")
    if sys.platform == "win32":
        mpiPath = Path(os.environ['MSMPI_BIN']) / Path("mpiexec.exe")

    paramDict['mpiCmd'] = mpiPath

#===============================================================================
#                The simulatio loop for velocity calculation
#===============================================================================
def simLoopVelocity(workingDir):
    nmax = calcMaxSimIterations()

    mpiPath = paramDict['mpiCmd']
    numProcessors = paramDict['numProcessors']

    nmin = 0
    start = 0.0
    with open("_data/Extrud3D_0.dat", "a") as f:
        f.write("\n[E3DSimulationSettings]\n")
        f.write("dAlpha=" + str(paramDict['deltaAngle']) + "\n")
        f.write("Periodicity=" + str(paramDict['periodicity']) + "\n")
        f.write("nSolutions=" + str(paramDict['timeLevels']) + "\n")

    for i in range(nmin, nmax):  # nmax means the loop goes to nmax-1
        if paramDict['singleAngle'] >= 0.0:
            angle = paramDict['singleAngle']
        else:
            angle = start + i * paramDict['deltaAngle']

        shutil.copyfile("_data/Extrud3D_0.dat", "_data/Extrud3D.dat")

        with open("_data/Extrud3D.dat", "a") as f:
            f.write("Angle=" + str(angle) + "\n")

        if sys.platform == "win32":
            subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse.exe"])
        else:
            #comm = subprocess.call(['mpirun', '-np', '%i' % numProcessors,  './q2p1_sse', '-a', '%d' % angle],shell=True)
            subprocess.call(['mpirun -np %i ./q2p1_sse' % (numProcessors)],shell=True)

        iangle = int(angle)
        if os.path.exists(Path("_data/prot.txt")):
            shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

#===============================================================================
#                The simulatio loop for velocity calculation
#===============================================================================
def cleanWorkingDir(workingDir):
    if not sys.platform == "win32":
        offList = list(workingDir.glob('*.off')) + list(workingDir.glob('*.OFF'))
    else: 
        offList = list(workingDir.glob('*.off'))

    for item in offList:
        os.remove(str(item))

#===============================================================================
#                The simulatio loop for velocity calculation
#===============================================================================
def simLoopTemperatureCombined(workingDir):
 
    print("Temperature simulation is activated!")
    
    numProcessors = paramDict['numProcessors']
    mpiPath = paramDict['mpiCmd']
    maxIterations = 2
    for iter in range(maxIterations):
        backupVeloFile = Path("_data_BU") / Path("q2p1_paramV_%01d.dat" % iter)
        backupTemperatureFile = Path("_data_BU") / Path("q2p1_paramT_%01d.dat" % iter)
        veloDestFile = Path("_data") / Path("q2p1_param.dat")
        temperatureDestFile = Path("_data") / Path("q2p1_paramT.dat")
        print("Copying: ", backupVeloFile, veloDestFile)
        print("Copying: ", backupTemperatureFile, temperatureDestFile)
        shutil.copyfile(str(backupVeloFile), str(veloDestFile))
        shutil.copyfile(str(backupTemperatureFile), str(temperatureDestFile))
        simLoopVelocity(workingDir)
        print("temperature simulation")

        if sys.platform == "win32":
            subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse_temp.exe"])
        else:
            #comm = subprocess.call(['mpirun', '-np', '%i' % numProcessors,  './q2p1_sse', '-a', '%d' % angle],shell=True)
            subprocess.call(['mpirun -np %i ./q2p1_sse_temp ' % (numProcessors)],shell=True)
        
        dirName = Path("_prot%01d" % iter)
        mkdir(dirName)
        protList = list(Path("_data").glob('prot*'))
        print(protList)
        for item in protList:
            shutil.copy(str(item), dirName)
            os.remove(item)

    backupVeloFile = Path("_data_BU") / Path("q2p1_paramV_%01d.dat" % maxIterations)
    veloDestFile = Path("_data") / Path("q2p1_param.dat")
    print("Copying: ", backupVeloFile, veloDestFile)
    shutil.copyfile(str(backupVeloFile), str(veloDestFile))
    simLoopVelocity(workingDir)

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

    try:
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:p:d:a:c:r:t:smhov',
                                   ['num-processors=', 'project-folder=',
                                    'periodicity=', 'delta-angle=', 'angle=',
                                    'host-conf=', 'rank-file=', 'time=', 'skip-setup',
                                    'skip-simulation', 'help','do-temperature','version'])

    except getopt.GetoptError:
        usage()
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-f', '--project-folder'):
            paramDict['projectFolder'] = arg
        elif opt in ('-n', '--num-processors'):
            paramDict['numProcessors'] = int(arg)
        elif opt in ('-p', '--periodicity'):
            paramDict['periodicity'] = int(arg)
        elif opt in ('-a', '--angle'):
            paramDict['singleAngle'] = float(arg)
        elif opt in ('-d', '--delta-angle'):
            paramDict['deltaAngle'] = float(arg)
            paramDict['hasDeltaAngle'] = True 
            if paramDict['deltaAngle'] <= 0.0:
                print("Parameter deltaAngle is set to a number <= 0 which is invalid. Please enter a number > 0")
                sys.exit(2)
        elif opt in ('-c', '--host-conf'):
            paramDict['hostFile'] = arg
        elif opt in ('-r', '--rank-file'):
            paramDict['rankFile'] = arg
        elif opt in ('-t', '--time'):
            paramDict['timeLevels'] = int(arg)
            paramDict['hasTimeLevels'] = True 
            if paramDict['timeLevels'] == 0:
                print("Parameter timeLevels is set to a number <= 0 which is invalid. Please enter a number > 0")
                sys.exit(2)
        elif opt in ('-s', '--skip-setup'):
            paramDict['skipSetup'] = True
        elif opt in ('-s', '--skip-simulation'):
            paramDict['skipSimulation'] = True
        elif opt in ('-o', '--do-temperature'):
            paramDict['temperature'] = True
        elif opt in ('-v', '--version'):
            version()
            sys.exit(2)
        else:
            usage()
            sys.exit(2)

    if paramDict['projectFolder'] == "":
        print("Error: no project folder specified.")
        usage()
        sys.exit(2)

    if paramDict['hasDeltaAngle'] and paramDict['hasTimeLevels']:
        print("Error: Specifying both deltaAngle and timeLevels at the same time is error-prone and therefore prohibited.")
        sys.exit(2)

    # Get the case/working dir paths
    projectFolder = paramDict['projectFolder'] 
    workingDir = Path('.')
    projectPath = Path(workingDir / projectFolder)
    projectFile = Path(projectPath / 'setup.e3d')

    setupMPICommand()

    if not paramDict['skipSetup']:
        simulationSetup(workingDir, projectFile, projectPath, projectFolder)

    if paramDict['skipSimulation']:
        sys.exit()

    print("  ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("                                        Launching E3D...                                           ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("  ")
    print("  ")

    if not paramDict['temperature']:
        simLoopVelocity(workingDir)
        cleanWorkingDir(workingDir)
    else:
        simLoopTemperatureCombined(workingDir)
        cleanWorkingDir(workingDir)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    with open("e3d.log", "w") as log:
        log.write("E3D started at: " + str(datetime.datetime.now()) + "\n")
        log.write("E3D running\n")
    main()
    with open("e3d.log", "a") as log:
        log.write("E3D stopped at: " + str(datetime.datetime.now()) + "\n")
