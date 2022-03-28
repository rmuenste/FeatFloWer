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
import fileinput
import datetime
import configparser
from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver  
from watchdog.events import FileSystemEventHandler

if sys.version_info[0] < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

debugNoSim = False
debugOutput = False

class E3dLog:
    def __init__(self):
#        print("The E3dLog constructor")
        self.fileName = "e3d.log"
        self.fileHandle = ""
        self.currPos = 0
        self.statusLineLength = 0
        self.statusLinePos = 0
        self.fileContents = []

    def __def__(self):
        #print("The E3dLog deconstructor")
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
        self.fileHandle.write("FileVersion=Extrud3D 2022\n")
        today = datetime.date.today()
        self.fileHandle.write("Date=%s \n" % today.strftime("%d/%m/%Y"))
        self.fileHandle.write("Extrud3DVersion=Extrud3D 2022.01\n")
        self.fileHandle.write("[SimulationStatus]\n")
        tempValue = ""
        if paramDict["temperature"]:
            tempValue = "true"
        else:
            tempValue = "false"            
        self.fileHandle.write("PathToE3DFile=%s\n" %str(projectFile))
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

    def updateStatusLineInnerIteration(self, msg):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-3] 

        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(msg)
            

    def updateStatusLineIteration(self, msg, itCurr):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        if itCurr == 0:
            self.fileContents = self.fileContents[:-1] 
        else:
            self.fileContents = self.fileContents[:-5] 
        
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(msg)
            

    def updateStatusLineHeatIteration(self, msg, numLines=3):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-3] 
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(msg)

    def writeStatusHeat(self, msg):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-1] 
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(msg)

    def popLinesAndWrite(self, numLines, msg):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-numLines] 
        with open(self.fileName, "w") as f:
            for line in self.fileContents:
                f.write(line)

            f.write(msg)


    def popLinesBack(self, numLines):

        with open(self.fileName, "r") as f:
            self.fileContents = f.readlines()

        self.fileContents = self.fileContents[:-numLines] 
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

myLog = E3dLog()

paramDict = {
    "deltaAngle": 10.0, # Angular step size
    "singleAngle": -10.0, # Single angle to compute 
    "hostFile" : "", # Hostfile
    "rankFile" : "" , # Rankfile 
    "timeLevels" : 36, # timeLevels
    "periodicity" : 1, # Periodicity 
    "numProcessors" : 0, # Number of processors 
    "projectFolder" : "", # The project folder
    "skipSetup" :  False,
    "shortTest" :  False, 
    "skipSimulation" : False,
    "hasDeltaAngle": False,
    "hasTimeLevels": False,
    "useSrun": False,
    "dieSimulation": False,
    "temperature" : False,
    "partialFilling" : False,
    "onlyMeshCreation" : False,
    "retryDeformation" : False
}

class ProtocolObserver(FileSystemEventHandler):
    def on_created(self, event):
        fileBaseName = os.path.basename(event.src_path)

    def on_modified(self, event):
        fileBaseName = os.path.basename(event.src_path)
        patternFound = False
        if fileBaseName == "prot.txt":
          pattern = "itns:\s*([0-9]+)[/]\s*([0-9]+)"
          with open(event.src_path, "r") as f:
              fileContents = f.readlines()
          for line in reversed(fileContents):
              matchObj = re.search(pattern, line)
              if matchObj: 
                  patternFound = True
                  statusMsg = "CurrentInnerIteration=%i\n"\
                              "MaxInnerIteration=%i\n"\
                              "CurrentStatus=running Momentum Solver" %(int(matchObj.group(1)), int(matchObj.group(2)))
                  myLog.updateStatusLineInnerIteration(statusMsg)
                  break
#===============================================================================


#===============================================================================
#                        version function
#===============================================================================
def version():
    """
    Print out version information
    """
    print("E3D + Reporter for SIGMA Version 2020.11.5, Copyright 2019 IANUS Simulation")
#===============================================================================


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
    print("[-x', '--short-test']: configures the program for a short test")
    print("[-u', '--use-srun']: Uses the srun launch mechanism")
    print("['--die-simulation']: fires up a single angle DIE sim with the corresponding datafile")
    print("Example: python ./e3d_start.py -f myFolder -n 5 -t 0")
#===============================================================================


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


#===============================================================================
#                          simple file in-situ replacement method
#===============================================================================
def replace_in_file(file_path, search_text, new_text):
    with fileinput.input(file_path, inplace=True) as f:
        for line in f:
            new_line = line.replace(search_text, new_text)
            print(new_line, end='')
#===============================================================================


#===============================================================================
#                     Parse MaxNumStep from q2p1_param.dat
#===============================================================================
def parseMaxNumSteps(file_path):
    maxIters = 300
    pattern = re.compile(r"SimPar@MaxNumStep\s+[=]\s+([0-9]+)")
    for line in open(file_path):
        for match in re.finditer(pattern, line):
            maxIters = int(match.group(1))
            return maxIters

    return maxIters
#===============================================================================


#===============================================================================
#                           e3dToDict 
#===============================================================================
def e3dToDict(pathName):
    config = configparser.ConfigParser()
    config.read(pathName)

    if not pathName.exists():
        raise FileNotFoundError("File {0} was not found.".format(pathName))

    return config
#===============================================================================


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
        if (paramDict['dieSimulation']) :
            backupDataFile = Path("_data_BU") / Path("q2p1_paramV_DIE_test.dat")
        elif (paramDict['partialFilling']) :
            backupDataFile = Path("_data_PF") / Path("q2p1_paramV_0_test.dat")
        else:
            backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU_test.dat")
    else:
        if (paramDict['dieSimulation']) :
            backupDataFile = Path("_data_BU") / Path("q2p1_paramV_DIE.dat")
        elif (paramDict['partialFilling']) :
            backupDataFile = Path("_data_PF") / Path("q2p1_paramV_0.dat")
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

#===============================================================================
#                           Mesher Step 
#===============================================================================
def mesherStep(workingDir, projectFile, projectPath, projectFolder):

    myLog.updateStatusLine("CurrentStatus=running Mesher")

    if sys.platform == "win32":
        exitCode = subprocess.call(["./s3d_mesher"])        
    else:
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
    
    return exitCode
#===============================================================================


#===============================================================================
#                           Partitioner Step 
#===============================================================================
def partitionerStep(workingDir, projectFile, projectPath, projectFolder):
    partitionerParameters = [1, 1]
    if paramDict['partialFilling']:
        partitionerParameters = [-3, 2]

    print("Partitioner parameters: ",-1, partitionerParameters)
    try:
        myLog.updateStatusLine("CurrentStatus=running Partitioner")
        partitioner.partition(paramDict['numProcessors']-1, partitionerParameters[0], partitionerParameters[1], "NEWFAC", "_data/meshDir/file.prj")
    except:
        myLog.logErrorExit("CurrentStatus=abnormal Termination Partitioner", 2)
    
#===============================================================================


#===============================================================================
#                             Simulation Setup 
#===============================================================================
def simulationSetup(workingDir, projectFile, projectPath, projectFolder):

    folderSetup(workingDir, projectFile, projectPath, projectFolder)

    exitCode = mesherStep(workingDir, projectFile, projectPath, projectFolder)
    
    partitionerStep(workingDir, projectFile, projectPath, projectFolder)
    
    return exitCode
#===============================================================================


#===============================================================================
#                         Only Mesh Creation
#===============================================================================
def onlyMeshCreation(workingDir, projectFile, projectPath, projectFolder):

    folderSetup(workingDir, projectFile, projectPath, projectFolder)

    exitCode = mesherStep(workingDir, projectFile, projectPath, projectFolder)
    
    return exitCode
#===============================================================================

    
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
    
    #print("nmax: ",nmax)

    return nmax
#===============================================================================


#===============================================================================
#                Compute maximum number of simulation iterations
#===============================================================================
def setupMPICommand():
    mpiPath = Path("mpirun")
    if sys.platform == "win32":
        mpiPath = Path(os.environ['MSMPI_BIN']) / Path("mpiexec.exe")

    paramDict['mpiCmd'] = mpiPath
#===============================================================================


#===============================================================================
#                The simulation loop for velocity calculation
#===============================================================================
def simLoopVelocity(workingDir):
    nmax = calcMaxSimIterations()

    mpiPath = paramDict['mpiCmd']
    numProcessors = paramDict['numProcessors']

    nmin = 0
    start = 0.0
    maxInnerIters = parseMaxNumSteps("_data/q2p1_param.dat")
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

        statusMsg = "CurrentAngleIteration=%i\n"\
                    "MaxAngleIteration=%i\n"\
                    "CurrentInnerIteration=%i\n"\
                    "MaxInnerIteration=%i\n"\
                    "CurrentStatus=running Momentum Solver" %(i+1, nmax, 1, maxInnerIters)
        myLog.updateStatusLineIteration(statusMsg, i)

        workingDir = os.getcwd()
        protocolFilePath = os.path.join(workingDir, "_data")

        eventHandler = ProtocolObserver()
        observer = ""   
        if sys.platform == "win32":
            observer = PollingObserver()
        else:
            observer = Observer()

        observer.schedule(eventHandler, path=protocolFilePath, recursive=False)
        observer.start()

        if sys.platform == "win32":
            exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse.exe"])
        else:
            launchCommand = ""

            if paramDict['useSrun']:
                launchCommand = "srun " + os.getcwd() + "/q2p1_sse"
                if paramDict['singleAngle'] >= 0.0:
                    launchCommand = launchCommand + " -a %d" %(angle)
            else:
                launchCommand = "mpirun -np " + str(numProcessors) + " " + os.getcwd() + "/q2p1_sse"
                if paramDict['singleAngle'] >= 0.0 :
                    launchCommand = launchCommand + " -a %d" %(angle)

            exitCode = subprocess.call([launchCommand], shell=True)

            if paramDict['retryDeformation'] and exitCode == 55:
                with open("_data/q2p1_param.dat", "r") as f:
                  for l in f:
                    if "SimPar@UmbrellaStepM" in l:
                        orig_umbrella = int(l.split()[2])
                UmbrellaStepM = orig_umbrella
                while exitCode == 55 and UmbrellaStepM != 0:
                    replace_in_file("_data/q2p1_param.dat", "SimPar@UmbrellaStepM = "+str(UmbrellaStepM), "SimPar@UmbrellaStepM = "+str(int(UmbrellaStepM/2)))
                    UmbrellaStepM = int(UmbrellaStepM / 2)
                    exitCode = subprocess.call([launchCommand], shell=True)
                replace_in_file("_data/q2p1_param.dat", "SimPar@UmbrellaStepM = "+str(UmbrellaStepM), "SimPar@UmbrellaStepM = "+str(orig_umbrella))

        # Here the observer can be turned off
        observer.stop()

        if exitCode == 88:
          myLog.logErrorExit("CurrentStatus=the screw could not be created: wrong angle", exitCode)

        if exitCode != 0:
            myLog.logErrorExit("CurrentStatus=abnormal Termination Momentum Solver", exitCode)

        # Write final inner iteration state
        statusMsg = "CurrentInnerIteration=%i\n"\
                    "MaxInnerIteration=%i\n"\
                    "CurrentStatus=running Momentum Solver" %(maxInnerIters, maxInnerIters)
        myLog.updateStatusLineInnerIteration(statusMsg)

        iangle = int(angle)
        if os.path.exists(Path("_data/prot.txt")):
            shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

    return exitCode    
#===============================================================================


#===============================================================================
#                The simulatio loop for velocity calculation
#===============================================================================
def cleanWorkingDir(workingDir):
    if not sys.platform == "win32":
        offList = list(workingDir.glob('*.off')) + list(workingDir.glob('*.OFF'))
    else: 
        offList = list(workingDir.glob('*.off'))

    # temporarily blocked ==> should be released later
    #for item in offList:
        #os.remove(str(item))
#===============================================================================


#===============================================================================
#                The simulatio loop for velocity calculation
#===============================================================================
def simLoopTemperatureCombined(workingDir):

    numProcessors = paramDict['numProcessors']
    mpiPath = paramDict['mpiCmd']
    maxIterations = 2
    for iter in range(maxIterations):
     
        if paramDict['shortTest']:
            backupVeloFile = Path("_data_BU") / Path("q2p1_paramV_%01d_test.dat" % iter)
        else:
            backupVeloFile = Path("_data_BU") / Path("q2p1_paramV_%01d.dat" % iter)
            
        if paramDict['shortTest']:
            backupTemperatureFile = Path("_data_BU") / Path("q2p1_paramT_test.dat")
        else:
            backupTemperatureFile = Path("_data_BU") / Path("q2p1_paramT_%01d.dat" % iter)
            
        veloDestFile = Path("_data") / Path("q2p1_param.dat")
        temperatureDestFile = Path("_data") / Path("q2p1_paramT.dat")
        print("Copying: ", backupVeloFile, veloDestFile)
        print("Copying: ", backupTemperatureFile, temperatureDestFile)
        shutil.copyfile(str(backupVeloFile), str(veloDestFile))
        shutil.copyfile(str(backupTemperatureFile), str(temperatureDestFile))

        if iter > 0:
            myLog.updateStatusLineHeatIteration("CurrentHeatIteration=%i\nHeatMaxIteration=%i\nCurrentStatus=running Heat Solver" %(iter+1, maxIterations))

        else:
            myLog.writeStatusHeat("CurrentHeatIteration=%i\nMaxHeatIteration=%i\nCurrentStatus=running Heat Solver" %(iter+1, maxIterations))
            
        exitCode = simLoopVelocity(workingDir)

#        statusMsg = "CurrentIteration=%i\nMaxIteration=%i\nCurrentStatus=running Momentum Solver" %(i+1, nmax)

        if sys.platform == "win32":
            exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse_temp.exe"])
        else:

            launchCommand = ""

            if paramDict['useSrun']:
                launchCommand = "srun " + os.getcwd() + "/q2p1_sse_temp"
                if paramDict['singleAngle'] >= 0.0:
                    launchCommand = launchCommand + " -a %d" %(angle)
            else:
                launchCommand = "mpirun -np " + str(numProcessors) + " " + os.getcwd() + "/q2p1_sse_temp"
                if paramDict['singleAngle'] >= 0.0 :
                    launchCommand = launchCommand + " -a %d" %(angle)

            exitCode = subprocess.call([launchCommand], shell=True)

        if exitCode != 0:
            myLog.logErrorExit("CurrentStatus=abnormal Termination Heat Solver", exitCode)

        myLog.popLinesAndWrite(5, "CurrentStatus=running Heat Solver")
        
        dirName = Path("_prot%01d" % iter)
        mkdir(dirName)
        protList = list(Path("_data").glob('prot*'))
        
        for item in protList:
            shutil.copy(str(item), dirName)
            os.remove(item)

    if paramDict['shortTest']:
        backupVeloFile = Path("_data_BU") / Path("q2p1_paramV_%01d_test.dat" % maxIterations)
    else:
        backupVeloFile = Path("_data_BU") / Path("q2p1_paramV_%01d.dat" % maxIterations)
        
    veloDestFile = Path("_data") / Path("q2p1_param.dat")
    print("Copying: ", backupVeloFile, veloDestFile)
    shutil.copyfile(str(backupVeloFile), str(veloDestFile))

    #myLog.popLinesAndWrite(3, "CurrentStatus=running Heat Solver")
    exitCode = simLoopVelocity(workingDir)
#===============================================================================


#===============================================================================
#                The cfd simulation loop for partial filling
#===============================================================================
def simLoopMainPartialFilling(workingDir, it, loops):

    nmax = calcMaxSimIterations()

    mpiPath = paramDict['mpiCmd']
    numProcessors = paramDict['numProcessors']

    nmin = 0
    start = 0.0
    maxInnerIters = parseMaxNumSteps("_data/q2p1_param.dat")
    if debugNoSim:
        print("Line 676, Updating maxInnerIters: ", maxInnerIters)
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

        statusMsg = "CurrentAngleIteration=%i\n"\
                    "MaxAngleIteration=%i\n"\
                    "CurrentInnerIteration=%i\n"\
                    "MaxInnerIteration=%i\n"\
                    "CurrentStatus=running Momentum Solver" %(it+1, loops+1, 1, maxInnerIters)

        myLog.updateStatusLineIteration(statusMsg, it)

        workingDir = os.getcwd()
        protocolFilePath = os.path.join(workingDir, "_data")

        eventHandler = ProtocolObserver()
        observer = ""   
        if sys.platform == "win32":
            observer = PollingObserver()
        else:
            observer = Observer()

        observer.schedule(eventHandler, path=protocolFilePath, recursive=False)
        observer.start()

        if not debugNoSim:
            if sys.platform == "win32":
                exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse_partfil.exe"])
            else:
                launchCommand = ""

                if paramDict['useSrun']:
                    launchCommand = "srun " + os.getcwd() + "/q2p1_sse"
                    if paramDict['singleAngle'] >= 0.0:
                        launchCommand = launchCommand + " -a %d" %(angle)
                else:
                    launchCommand = "mpirun -np " + str(numProcessors) + " " + os.getcwd() + "/q2p1_sse_partfil"
                    if paramDict['singleAngle'] >= 0.0 :
                        launchCommand = launchCommand + " -a %d" %(angle)

                exitCode = subprocess.call([launchCommand], shell=True)

                if paramDict['retryDeformation'] and exitCode == 55:
                    with open("_data/q2p1_param.dat", "r") as f:
                        for l in f:
                            if "SimPar@UmbrellaStepM" in l:
                                orig_umbrella = int(l.split()[2])
                        UmbrellaStepM = orig_umbrella
                        while exitCode == 55 and UmbrellaStepM != 0:
                            replace_in_file("_data/q2p1_param.dat", "SimPar@UmbrellaStepM = "+str(UmbrellaStepM), "SimPar@UmbrellaStepM = "+str(int(UmbrellaStepM/2)))
                            UmbrellaStepM = int(UmbrellaStepM / 2)
                            exitCode = subprocess.call([launchCommand], shell=True)
                        replace_in_file("_data/q2p1_param.dat", "SimPar@UmbrellaStepM = "+str(UmbrellaStepM), "SimPar@UmbrellaStepM = "+str(orig_umbrella))
        else:
            input("Line 743, DebugMode, SimLoopMainPartialFilling: q2p1_sse_partfil")
            exitCode = 0   
        # Here the observer can be turned off
        observer.stop()
        
        if exitCode == 88:
          myLog.logErrorExit("CurrentStatus=the screw could not be created: wrong angle", exitCode)

        if exitCode != 0:
            myLog.logErrorExit("CurrentStatus=abnormal Termination Momentum Solver", exitCode)

        # Write final inner iteration state
        statusMsg = "CurrentInnerIteration=%i\n"\
                    "MaxInnerIteration=%i\n"\
                    "CurrentStatus=running Multiphase Momentum Solver" %(maxInnerIters, maxInnerIters)
        myLog.updateStatusLineInnerIteration(statusMsg)

        iangle = int(angle)
        if os.path.exists(Path("_data/prot.txt")):
            shutil.copyfile("_data/prot.txt", "_data/prot_%04d.txt" % iangle)

    return exitCode    
#===============================================================================


#===============================================================================
#                The simulation loop for partial filling
#===============================================================================
def simLoopPartialFilling(workingDir):

    mpiPath = paramDict['mpiCmd']
    numProcessors = paramDict['numProcessors']

    if paramDict['shortTest']:
        nLoops = 2 
    else:
        nLoops = 4 

    # Start the initial loop
    it = 0
    if debugNoSim:
        input("Line 784, DebugMode, SimLoopPartialFilling: Initial q2p1_sse_partfil iteration")
    simLoopMainPartialFilling(workingDir, it, nLoops)
    it = it + 1

    
    for i in range(nLoops):
        if paramDict['shortTest']:
            sourceParamFile = Path("_data_PF") / Path("q2p1_paramAlpha_test.dat")
        else:
            sourceParamFile = Path("_data_PF") / Path("q2p1_paramAlpha.dat")
        
        destBackupFile = Path("_data") / Path("q2p1_param.dat")
        shutil.copyfile(str(sourceParamFile), str(destBackupFile))
        msg = "Copied %s => %s " %(str(sourceParamFile), str(destBackupFile))
        if debugNoSim:
            print(msg)

        sourceParamFile = Path("_data_PF") / Path("mesh_names.offs")
        destBackupFile = Path("./mesh_names.offs")
        shutil.copyfile(str(sourceParamFile), str(destBackupFile))
        msg = "Copied %s => %s " %(str(sourceParamFile), str(destBackupFile))
        #print(msg)

        if not debugNoSim:
            if sys.platform == "win32":
                exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q1_scalar_partfil.exe"])
            else:
                launchCommand = ""

                if paramDict['useSrun']:
                    launchCommand = "srun " + os.getcwd() + "/q1_scalar_partfil"
                else:
                    launchCommand = "mpirun -np " + str(numProcessors) + " " + os.getcwd() + "/q1_scalar_partfil"

                exitCode = subprocess.call([launchCommand], shell=True)        
        else:
            input("Line 824, DebugMode, SimLoopPartialFilling: q1_scalar_partfil")
            exitCode = 0   

        if paramDict['shortTest']:
            sourceParamFile = Path("_data_PF") / Path("q2p1_paramV_1_test.dat")
        else:
            sourceParamFile = Path("_data_PF") / Path("q2p1_paramV_1.dat")

        destBackupFile = Path("_data") / Path("q2p1_param.dat")
        shutil.copyfile(str(sourceParamFile), str(destBackupFile))
        msg = "Copied %s => %s " %(str(sourceParamFile), str(destBackupFile))
        if debugNoSim:
            print(msg)

        simLoopMainPartialFilling(workingDir, it, nLoops)
        it = it + 1
#===============================================================================
    

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
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:p:d:a:c:r:t:smxhovu',
                                   ['num-processors=', 'project-folder=',
                                    'periodicity=', 'delta-angle=', 'angle=',
                                    'host-conf=', 'rank-file=', 'time=', 'skip-setup','die-simulation',
                                    'skip-simulation','short-test', 'help',
                                    'do-temperature','version', 'use-srun',
                                    'retry-deformation', 'partial-filling','only-mesh-creation'])

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
        elif opt in ('-m', '--skip-simulation'):
            paramDict['skipSimulation'] = True
        elif opt in ('-o', '--do-temperature'):
            paramDict['temperature'] = True
        elif opt in ('-v', '--version'):
            version()
            sys.exit(2)
        elif opt in ('-x', '--short-test'):
            paramDict['shortTest'] = True
        elif opt in ('-u', '--use-srun'):
            paramDict['useSrun'] = True
        elif opt in ('--die-simulation'):
            paramDict['dieSimulation'] = True
        elif opt in ('--retry-deformation'):
            paramDict['retryDeformation'] = True
        elif opt in ('--partial-filling'):
            paramDict['partialFilling'] = True
            paramDict['singleAngle'] = 0.0 
        elif opt in ('--only-mesh-creation'):
            paramDict['onlyMeshCreation'] = True
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
        
    if (paramDict['dieSimulation']) :
        print("Switching to 'DIE' simulation !")
        paramDict['singleAngle'] = 0

    if (paramDict['numProcessors'] < 3 and not paramDict['onlyMeshCreation']) :
        print("Number of processors should be > 3")
        sys.exit(2)

#    if (paramDict['onlyMeshCreation']) :
#        print("Only doing mesh creation")
#        sys.exit(2)
     
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

    e3dSetupDict = e3dToDict(projectFile)

    if not paramDict['hasTimeLevels']:
        if "timeLevels" in e3dSetupDict['SimodSettings']:
          if e3dSetupDict['SimodSettings']['timeLevels'].isnumeric():
              paramDict['timeLevels'] = int(e3dSetupDict['SimodSettings']['timeLevels'])
              paramDict['hasTimeLevels'] = True
          else:
              raise TypeError("e3d.setup ['SimodSettings']['timeLevels'] is not a numeric entry. Pls enter a number > 0.")


    if not paramDict['skipSetup']:
        if paramDict['onlyMeshCreation']:
            exitCode = onlyMeshCreation(workingDir, projectFile, projectPath, projectFolder)
            return
        else:
            exitCode = simulationSetup(workingDir, projectFile, projectPath, projectFolder)

    if paramDict['skipSimulation']:
        sys.exit()

    print("  ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("                                        Launching E3D...                                           ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("  ")
    print("  ")

    if paramDict['partialFilling']:
        simLoopPartialFilling(workingDir)
        cleanWorkingDir(workingDir)
    elif paramDict['temperature']:
        simLoopTemperatureCombined(workingDir)
        cleanWorkingDir(workingDir)
    else:
        simLoopVelocity(workingDir)
        cleanWorkingDir(workingDir)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
    myLog.writeExitMsg()
