#!/usr/bin/env python
# vim: set filetype=python
"""
This module is the driver for the s3d simulator

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

class S3dLog:
    def __init__(self):
#        print("The S3dLog constructor")
        self.fileName = "s3d.log"
        self.fileHandle = ""
        self.currPos = 0
        self.statusLineLength = 0
        self.statusLinePos = 0
        self.fileContents = []

    def __def__(self):
        #print("The S3dLog deconstructor")
        self.fileHandle.close()

    def openFileHandle(self):
        self.fileHandle = open(self.fileName, "w")

    def closeFileHandle(self):
        self.fileHandle.close()

    def appendLine(self, msg):
        with open(self.fileName, "a") as f:
            f.write(msg + "\n")

    def writeNumAnglePositions(self, numAngPos):
        self.fileHandle.write("NumOfAnglePositions=%i\n" %numAngPos)        
        self.fileHandle.write("Periodicity=%i\n" %paramDict['periodicity'])        
        self.currPos = self.fileHandle.tell()
        self.fileHandle.flush()

    def writeGeneralInfo(self, parameterDict):

        projectFolder = parameterDict['projectFolder'] 
        workingDir = Path('.')
        projectPath = Path(workingDir / projectFolder)
        projectFile = Path(projectPath / 'setup.e3d')

        self.fileHandle.write("[Extrud3DFileInfo]\n")
        self.fileHandle.write("FileType=Logging\n")
        self.fileHandle.write("FileVersion=Extrud3D 2020\n")
        today = datetime.date.today()
        self.fileHandle.write("Date=%s \n" % today.strftime("%d/%m/%Y"))
        self.fileHandle.write("Extrud3DVersion=Extrud3D 2020.02\n")
        self.fileHandle.write("[SimulationStatus]\n")
        tempValue = ""
        if paramDict["temperature"]:
            tempValue = "true"
        else:
            tempValue = "false"            
        self.fileHandle.write("PathToS3DFile=%s\n" %str(projectFile))
        self.fileHandle.write("TemperatureCalculation=%s\n" %tempValue)        
        self.fileHandle.write("NumOfCpus=%i\n" %parameterDict['numProcessors'])
        self.fileHandle.write("StartingTime=" + str(datetime.datetime.now()) + "\n")
        self.fileHandle.flush()

        self.currPos = self.fileHandle.tell()

    def updateStatusLine(self, msg):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-1] 
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(msg)

    def writeExitMsg(self):
        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-1] 
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write("CurrentStatus=finished\n")
            f.write("FinishingTime=" + str(datetime.datetime.now()) + "\n")

    def writeStatusLine(self):
        self.openFileHandle()
        self.fileHandle.write("CurrentStatus=running")
        self.closeFileHandle()

    def writeStatusLine2(self):
        self.fileHandle.write("CurrentStatus=running")
        self.currPos = self.fileHandle.tell()

    def logErrorExit(self, message, errorCode):
        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-1] 
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(message + "\n")
            f.write("ErrorCode=%i\n" % errorCode)
            f.write("FinishingTime=" + str(datetime.datetime.now()) + "\n")
            sys.exit(2)        

myLog = S3dLog()

paramDict = {
    "deltaAngle": 10.0, # Angular step size
    "singleAngle": 0.0, # Single angle to compute 
    "hostFile" : "", # Hostfile
    "rankFile" : "" , # Rankfile 
    "timeLevels" : 1, # timeLevels
    "periodicity" : 1, # Periodicity 
    "numProcessors" : 5, # Number of processors 
    "projectFolder" : "", # The project folder
    "skipSetup" :  False,
    "shortTest" :  False, 
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
    print("S3D + Reporter for REX Version 19.10, Copyright 2019 IANUS Simulation")

#===============================================================================
#                        usage function
#===============================================================================
def usage():
    """
    Print out usage information
    """
    print("Usage: python s3d_start.py [options]")
    print("Where options can be:")
    print("[-f', '--project-folder]: Path to project folder containing a setup.e3d file")
    print("[-n', '--num-processors]: Number of processors to use")
    #print("[-p', '--periodicity']: Periodicity of the solution (1, 2, 3, ... " +
          #"usually the time flight number)")
    #print("[-a', '--angle]: If this parameter is present a single simulation with " +
          #"that specific angle will be done.")
    #print("[-d', '--delta-angle]: The angular step size between two simulations " +
          #"the sim loop (default 10)")
    print("[-c', '--host-conf]: A hostfile as input for the mpirun command")
    print("[-r', '--rank-file]: A rankfile as input for the mpirun command")
    #print("[-t', '--time]: Number of time levels to complete a full 360 rotation")
    #print("['-o', '--do-temperature']: The simulation loop will do a temperature simulation")
    print("[-h', '--help']: prints this message")
    print("[-v', '--version']: prints out version information")
    print("[-x', '--short-test']: configures the program for a short test")
    print("Example: python ./s3d_start.py -f myFolder -n 5")

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

    if paramDict['shortTest']:
        backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU_test.dat")
    else:
        backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU.dat")
    
#    backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU.dat")
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

    myLog.updateStatusLine("CurrentStatus=running Mesher")

    if sys.platform == "win32":
        exitCode = subprocess.call(["./s3d_mesher"])        
    else:
        #comm = subprocess.call(['mpirun', '-np', '%i' % numProcessors,  './q2p1_sse', '-a', '%d' % angle],shell=True)
        exitCode = subprocess.call(["./s3d_mesher"], shell=True)


    if exitCode != 0:
        myLog.logErrorExit("CurrentStatus=abnormal Termination Mesher", exitCode)

    if not Path("_data/meshDir").exists():
        meshDirPath = projectPath / Path("meshDir")
        if meshDirPath.exists():
            print('Copying meshDir from Project Folder!')
            shutil.copytree(str(meshDirPath), "_data/meshDir")
        else:
            print("Error: No mesh automatically generated and no <meshDir> " + 
                  "folder present the case folder " + str(projectPath))
            sys.exit(2)
    
#    input("Press key to continue to Partitioner")
    try:
        myLog.updateStatusLine("CurrentStatus=running Partitioner")
        partitioner.partition(paramDict['numProcessors']-1, 1, 1, "NEWFAC", "_data/meshDir/file.prj")
    except:
        myLog.logErrorExit("CurrentStatus=abnormal Termination Partitioner", 2)
    
    return exitCode
    
#===============================================================================
#                Compute maximum number of simulation iterations
#===============================================================================
def calcMaxSimIterations():
    nmax = 1

    #if paramDict['hasDeltaAngle']:
        #paramDict['timeLevels'] = 360.0 / paramDict['deltaAngle'] 

    #if paramDict['timeLevels'] == 1:
        #nmax = 1
    #else:
        #paramDict['deltaAngle'] = 360.0 / float(paramDict['timeLevels'])
        #nmax = int(math.ceil(360.0 / paramDict['periodicity'] / paramDict['deltaAngle']))

    #if paramDict['singleAngle'] >= 0.0:
        #nmax = 1
    
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
        f.write("\n[S3DSimulationSettings]\n")
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

#        input("Press key to continue to MomentumSolver")
        myLog.updateStatusLine("CurrentStatus=running Momentum Solver")
        if sys.platform == "win32":
            exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse.exe"])
        else:
            #comm = subprocess.call(['mpirun', '-np', '%i' % numProcessors,  './q2p1_sse', '-a', '%d' % angle],shell=True)
            #exitCode = subprocess.call(['mpirun -np %i ./q2p1_sse' % (numProcessors)], shell=True)
            if paramDict['singleAngle'] >= 0.0 :
                exitCode = subprocess.call(['mpirun -np %i ./q2p1_sse -a %d' % (numProcessors, angle)], shell=True)
            else:
                exitCode = subprocess.call(['mpirun -np %i ./q2p1_sse' % (numProcessors)], shell=True)


        if exitCode != 0:
            myLog.logErrorExit("CurrentStatus=abnormal Termination Momentum Solver", exitCode)

        iangle = int(angle)
        if os.path.exists(Path("_data/prot.txt")):
            shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

    return exitCode    

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
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:p:d:a:c:r:t:smxhov',
                                   ['num-processors=', 'project-folder=',
                                    'periodicity=', 'delta-angle=', 'angle=',
                                    'host-conf=', 'rank-file=', 'time=', 'skip-setup',
                                    'skip-simulation','short-test', 'help','do-temperature','version'])

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
        #elif opt in ('-p', '--periodicity'):
            #paramDict['periodicity'] = int(arg)
        #elif opt in ('-a', '--angle'):
            #paramDict['singleAngle'] = float(arg)
        #elif opt in ('-d', '--delta-angle'):
            #paramDict['deltaAngle'] = float(arg)
            #paramDict['hasDeltaAngle'] = True 
            #if paramDict['deltaAngle'] <= 0.0:
                #print("Parameter deltaAngle is set to a number <= 0 which is invalid. Please enter a number > 0")
                #sys.exit(2)
        elif opt in ('-c', '--host-conf'):
            paramDict['hostFile'] = arg
        elif opt in ('-r', '--rank-file'):
            paramDict['rankFile'] = arg
        #elif opt in ('-t', '--time'):
            #paramDict['timeLevels'] = int(arg)
            #paramDict['hasTimeLevels'] = True 
            #if paramDict['timeLevels'] == 0:
                #print("Parameter timeLevels is set to a number <= 0 which is invalid. Please enter a number > 0")
                #sys.exit(2)
        elif opt in ('-s', '--skip-setup'):
            paramDict['skipSetup'] = True
        elif opt in ('-m', '--skip-simulation'):
            paramDict['skipSimulation'] = True
        elif opt in ('-v', '--version'):
            version()
            sys.exit(2)
        elif opt in ('-x', '--short-test'):
            paramDict['shortTest'] = True
        else:
            usage()
            sys.exit(2)

    if paramDict['projectFolder'] == "":
        print("Error: no project folder specified.")
        usage()
        sys.exit(2)

    if (paramDict['hasDeltaAngle'] and paramDict['hasTimeLevels']):
        print("Error: Specifying both deltaAngle and timeLevels at the same time is error-prone and therefore prohibited.")
        sys.exit(2)
        
    if (paramDict['singleAngle'] >= 0.0 and  paramDict['temperature']) :
        print("Error: Specifying both singleAngle and Temperature Simulation at the same time is prohibited.")
        sys.exit(2)
        
    # Get the case/working dir paths
    projectFolder = paramDict['projectFolder'] 
    workingDir = Path('.')
    projectPath = Path(workingDir / projectFolder)
    projectFile = Path(projectPath / 'setup.e3d')

    myLog.openFileHandle()
    myLog.writeGeneralInfo(paramDict)
    myLog.writeNumAnglePositions(calcMaxSimIterations())
    myLog.writeStatusLine2()
    myLog.closeFileHandle()

    setupMPICommand()

    if not paramDict['skipSetup']:
        exitCode = simulationSetup(workingDir, projectFile, projectPath, projectFolder)

    if paramDict['skipSimulation']:
        sys.exit()

    print("  ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("                                        Launching S3D...                                           ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("  ")
    print("  ")

    simLoopVelocity(workingDir)
    cleanWorkingDir(workingDir)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
    myLog.writeExitMsg()
