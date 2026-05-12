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
import yaml
from watchdog.observers import Observer
from watchdog.observers.polling import PollingObserver  
from watchdog.events import FileSystemEventHandler

if sys.version_info[0] < 3:
    from pathlib2 import Path
else:
    from pathlib import Path

debugNoSim = False
debugOutput = False

# Defaults required for reliable MPI file I/O on OpenMPI based runs
MPI_ENV_DEFAULTS = {
    "OMPI_MCA_io": "romio321",
    "ROMIO_CB_BUFFER_SIZE": "16777216",
    "ROMIO_DS_WRITE": "enable",
}

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
            sys.exit(errorCode)        

myLog = E3dLog()

paramDict = {
    "deltaAngle": 10.0, # Angular step size
    "singleAngle": -10.0, # Single angle to compute 
    "startAngle": 0.0, # Stage-local start angle for periodic runs
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
    "hasPeriodicity": False,
    "hasTimeLevels": False,
    "useSrun": False,
    "dieSimulation": False,
    "temperature" : False,
    "onlyMeshCreation" : False,
    "retryDeformation" : False,
    "resolutionLevel" : -1,
    "partitionFormat": "legacy",
    "yamlConfig": ""
}

XSE_PARAM_DIR = Path("_data_BU") / Path("XSE")
DIE_PARAM_DIR = Path("_data_BU") / Path("DIE")
stage_counter = 0

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
#                    Stage output archiving helper
#===============================================================================
def archive_stage_outputs(stage_idx, include_heat):
    stage_dir = Path(f"_prot_{stage_idx}")
    if stage_dir.exists():
        if stage_dir.is_dir():
            shutil.rmtree(stage_dir)
        else:
            stage_dir.unlink()
    stage_dir.mkdir()
    any_copied = False
    for item in Path("_data").glob('prot*'):
        target_name = "prot_TEMP.out" if include_heat and item.name == "prot.txt" else item.name
        shutil.copy(str(item), stage_dir / target_name)
        os.remove(item)
        any_copied = True
    if not any_copied:
        print(f"[Stage {stage_idx}] Warning: no protocol files found to archive.")
#===============================================================================


#===============================================================================
#                    Solver execution helpers
#===============================================================================
def run_momentum_step(param_file, workingDir):
    source = Path(param_file)
    if not source.exists():
        raise FileNotFoundError(f"Momentum param file missing: {param_file}")
    shutil.copyfile(str(source), "_data/q2p1_param.dat")
    apply_resolution_level_to_param_file(Path("_data") / Path("q2p1_param.dat"))
    return simLoopVelocity(workingDir)


def run_material_step(param_file):
    source = Path(param_file)
    if not source.exists():
        raise FileNotFoundError(f"MaterialDistribution param file missing: {param_file}")
    shutil.copyfile(str(source), "_data/q2p1_paramAlpha.dat")
    apply_resolution_level_to_param_file(Path("_data") / Path("q2p1_paramAlpha.dat"))

    names_source = Path("_data_BU/DIE/mesh_names.offs")
    if names_source.exists():
        shutil.copyfile(str(names_source), "mesh_names.offs")
    else:
        print("Warning: missing _data_BU/DIE/mesh_names.offs; q1_scalar_multimat may fail.")

    numProcessors = paramDict['numProcessors']
    mpiPath = paramDict['mpiCmd']
    angle = paramDict['singleAngle'] if paramDict['singleAngle'] >= 0.0 else 0.0

    if sys.platform == "win32":
        exec_name = "./q1_scalar_multimat.exe"
        exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  exec_name])
    else:
        if paramDict['useSrun']:
            launchCommand = "srun " + os.getcwd() + "/q1_scalar_multimat"
        else:
            launchCommand = "mpirun -np " + str(numProcessors) + " " + os.getcwd() + "/q1_scalar_multimat"
        if paramDict['singleAngle'] >= 0.0:
            launchCommand = launchCommand + " -a %d" %(angle)
        exitCode = subprocess.call([launchCommand], shell=True)

    if exitCode != 0:
        myLog.logErrorExit("CurrentStatus=abnormal Termination Material Distribution Solver", exitCode)

    prot_src = Path("_data/prot.txt")
    if prot_src.exists():
        prot_alpha = Path("_data/prot_ALPHA.txt")
        shutil.copyfile(str(prot_src), str(prot_alpha))
        prot_src.unlink()

    return exitCode


def run_heat_step(param_file):
    source = Path(param_file)
    if not source.exists():
        raise FileNotFoundError(f"Heat param file missing: {param_file}")
    shutil.copyfile(str(source), "_data/q2p1_paramT.dat")
    apply_resolution_level_to_param_file(Path("_data") / Path("q2p1_paramT.dat"))

    numProcessors = paramDict['numProcessors']
    mpiPath = paramDict['mpiCmd']
    angle = paramDict['singleAngle'] if paramDict['singleAngle'] >= 0.0 else 0.0

    if sys.platform == "win32":
        exitCode = subprocess.call([r"%s" % str(mpiPath), "-n",  "%i" % numProcessors,  "./q2p1_sse_temp.exe"])
    else:
        if paramDict['useSrun']:
            launchCommand = "srun " + os.getcwd() + "/q2p1_sse_temp"
        else:
            launchCommand = "mpirun -np " + str(numProcessors) + " " + os.getcwd() + "/q2p1_sse_temp"
        if paramDict['singleAngle'] >= 0.0:
            launchCommand = launchCommand + " -a %d" %(angle)
        exitCode = subprocess.call([launchCommand], shell=True)

    if exitCode != 0:
        myLog.logErrorExit("CurrentStatus=abnormal Termination Heat Solver", exitCode)
    return exitCode


def execute_yaml_plan(stages, workingDir):
    global stage_counter
    stage_counter = 0
    for phase in ['init', 'main', 'final']:
        if phase not in stages:
            continue
        stage_def = stages.get(phase, {})
        iterations = stage_def.get('loop', 1)
        try:
            iterations = int(iterations)
        except (TypeError, ValueError):
            raise ValueError(f"Stage '{phase}' has non-integer loop count '{iterations}'")
        if iterations < 1:
            print(f"[Stage {phase}] loop count {iterations} <= 0, skipping")
            continue
        steps = stage_def.get('steps', [])
        if not isinstance(steps, list):
            raise ValueError(f"Stage '{phase}' must define a list of steps")
        if not steps:
            print(f"[Stage {phase}] Warning: no steps defined, skipping.")
            continue

        stage_start_angle = 0.0
        if phase == 'init':
            stage_start_angle = stage_def.get('startangle', stage_def.get('start-angle', 0.0))
            try:
                stage_start_angle = float(stage_start_angle)
            except (TypeError, ValueError):
                raise ValueError("Stage 'init' option 'startangle' must be a number")
            if stage_start_angle < 0.0:
                raise ValueError("Stage 'init' option 'startangle' must be >= 0")
            if paramDict['dieSimulation'] and stage_start_angle > 0.0:
                raise ValueError("Stage 'init' option 'startangle' is not supported for die simulations")
            if paramDict['singleAngle'] >= 0.0 and stage_start_angle > 0.0:
                raise ValueError("Stage 'init' option 'startangle' cannot be combined with a single-angle run")
        elif 'startangle' in stage_def or 'start-angle' in stage_def:
            raise ValueError(f"Stage '{phase}' does not support the 'startangle' option")

        for iter_idx in range(iterations):
            heat_present = any(step.get('solver') == 'HeatEquation' for step in steps)
            print(f"[Stage {stage_counter}] Phase '{phase}' iteration {iter_idx+1}/{iterations}")
            paramDict['startAngle'] = stage_start_angle
            for step in steps:
                solver = step.get('solver')
                param_file = step.get('param_file')
                print(f"  -> {solver} using {param_file}")
                if solver == 'MomentumEquation':
                    run_momentum_step(param_file, workingDir)
                elif solver == 'HeatEquation':
                    run_heat_step(param_file)
                elif solver == 'MaterialDistribution':
                    run_material_step(param_file)
                else:
                    raise ValueError(f"Unknown solver '{solver}' in YAML plan")
            archive_stage_outputs(stage_counter, include_heat=heat_present)
            print(f"[Stage {stage_counter}] archived outputs to _prot_{stage_counter}")
            stage_counter += 1
#===============================================================================


def apply_yaml_options(plan):
    if not isinstance(plan, dict):
        raise ValueError("YAML configuration must contain a mapping at the top level")

    yaml_options = {}
    if isinstance(plan.get('options'), dict):
        yaml_options.update(plan['options'])
    if isinstance(plan.get('settings'), dict):
        yaml_options.update(plan['settings'])

    # Allow both top-level keys and grouped options/settings mappings.
    for key in ('die-simulation', 'dieSimulation'):
        if key in plan:
            yaml_options[key] = plan[key]

    die_simulation = yaml_options.get('die-simulation', yaml_options.get('dieSimulation'))
    if die_simulation is not None:
        if not isinstance(die_simulation, bool):
            raise ValueError("YAML option 'die-simulation' must be true or false")
        if die_simulation:
            paramDict['dieSimulation'] = True

    periodicity = yaml_options.get('periodicity')
    if periodicity is not None and not paramDict['hasPeriodicity']:
        try:
            periodicity = int(periodicity)
        except (TypeError, ValueError):
            raise ValueError("YAML option 'periodicity' must be an integer")
        if periodicity < 1:
            raise ValueError("YAML option 'periodicity' must be >= 1")
        paramDict['periodicity'] = periodicity

    delta_angle = yaml_options.get('delta-angle', yaml_options.get('deltaAngle'))
    if delta_angle is not None and not paramDict['hasDeltaAngle']:
        try:
            delta_angle = float(delta_angle)
        except (TypeError, ValueError):
            raise ValueError("YAML option 'delta-angle' must be a number")
        if delta_angle <= 0.0:
            raise ValueError("YAML option 'delta-angle' must be > 0")
        paramDict['deltaAngle'] = delta_angle
        paramDict['hasDeltaAngle'] = True

    resolution_level = yaml_options.get(
        'ResolutionLevel',
        yaml_options.get('resolution-level', yaml_options.get('resolutionLevel'))
    )
    if resolution_level is not None:
        try:
            resolution_level = int(resolution_level)
        except (TypeError, ValueError):
            raise ValueError("YAML option 'ResolutionLevel' must be an integer")
        if resolution_level < 2:
            raise ValueError("YAML option 'ResolutionLevel' must be >= 2")
        paramDict['resolutionLevel'] = resolution_level
#===============================================================================


#===============================================================================
#                    YAML input preflight validation
#===============================================================================
def validate_yaml_stage_inputs(plan, yaml_path):
    stages = plan.get('stages')
    if stages is None:
        raise ValueError("YAML configuration missing 'stages' section")
    if not isinstance(stages, dict):
        raise ValueError("YAML option 'stages' must be a mapping")

    defined_stages = [phase for phase in ['init', 'main', 'final'] if phase in stages]
    if not defined_stages:
        raise ValueError("YAML stages must define at least one of 'init', 'main', or 'final'")

    missing_files = []
    seen_files = set()

    for phase in defined_stages:
        stage_def = stages.get(phase, {})
        if not isinstance(stage_def, dict):
            raise ValueError(f"Stage '{phase}' must be a mapping")

        steps = stage_def.get('steps', [])
        if not isinstance(steps, list):
            raise ValueError(f"Stage '{phase}' must define a list of steps")
        if not steps:
            print(f"[Stage {phase}] Warning: no steps defined, skipping.")
            continue

        for idx, step in enumerate(steps, start=1):
            if not isinstance(step, dict):
                raise ValueError(f"Stage '{phase}' step {idx} must be a mapping")

            param_file = step.get('param_file')
            if not isinstance(param_file, str) or not param_file.strip():
                raise ValueError(f"Stage '{phase}' step {idx} must define a non-empty 'param_file'")

            if param_file in seen_files:
                continue
            seen_files.add(param_file)

            source = Path(param_file)
            if not source.exists():
                missing_files.append(param_file)

    if missing_files:
        print(f"Error: missing YAML input files referenced by {yaml_path}:")
        for missing_file in missing_files:
            print(f"  - {missing_file}")
        sys.exit(2)

    print(f"YAML preflight check passed: {len(seen_files)} input parameter file(s) found.")
#===============================================================================


#===============================================================================
#                    Parse partition format from parameter file
#===============================================================================
def parsePartitionFormat(paramFile):
    """
    Returns 'legacy' or 'json' based on SimPar@PartitionFormat entry.
    Defaults to legacy on missing/invalid values.
    """
    fmt = "legacy"
    try:
        with open(paramFile, "r") as f:
            for raw_line in f:
                line = raw_line.split("!")[0].strip()
                if not line:
                    continue
                lower_line = line.lower()
                if lower_line.startswith("simpar@partitionformat"):
                    parts = line.split("=", 1)
                    if len(parts) == 2:
                        candidate = parts[1].strip().strip('"').strip("'").lower()
                        if candidate in ("legacy", "json"):
                            fmt = candidate
                    break
    except OSError:
        pass
    return fmt
#===============================================================================



#===============================================================================
#                    Detect number of allocated nodes
#===============================================================================
def detect_node_count():
    """
    Returns the number of compute nodes assigned to this run.
    Tries scheduler env vars first, then host/rank files.
    """
    for env_var in ("SLURM_STEP_NUM_NODES", "SLURM_JOB_NUM_NODES"):
        val = os.environ.get(env_var)
        if val:
            try:
                nodes = int(val)
                if nodes > 0:
                    return nodes
            except ValueError:
                pass

    nodelist = os.environ.get("SLURM_JOB_NODELIST")
    if nodelist:
        try:
            output = subprocess.check_output(
                ["scontrol", "show", "hostnames", nodelist],
                universal_newlines=True,
            )
            hosts = [line.strip() for line in output.splitlines() if line.strip()]
            if hosts:
                return len(set(hosts))
        except (OSError, subprocess.SubprocessError):
            pass

    hostfile = paramDict.get('hostFile', "")
    if hostfile:
        try:
            hosts = set()
            with open(hostfile, "r") as f:
                for raw in f:
                    line = raw.split("#", 1)[0].strip()
                    if not line:
                        continue
                    token = line.split()[0]
                    hosts.add(token.split(":", 1)[0])
            if hosts:
                return len(hosts)
        except OSError:
            pass

    rankfile = paramDict.get('rankFile', "")
    if rankfile:
        try:
            hosts = set()
            pattern = re.compile(r"rank\s+\d+\s*=\s*([^\s]+)")
            with open(rankfile, "r") as f:
                for raw in f:
                    match = pattern.search(raw)
                    if match:
                        hosts.add(match.group(1))
            if hosts:
                return len(hosts)
        except OSError:
            pass

    return 0
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
#                     Apply mesh resolution override to param files
#===============================================================================
def get_effective_resolution_level():
    if paramDict['resolutionLevel'] != -1:
        return paramDict['resolutionLevel']
    return -1


def apply_resolution_level_to_param_file(file_path):
    resolution_level = get_effective_resolution_level()
    if resolution_level == -1:
        return
    lower_resolution_level = resolution_level - 1

    filepath = Path(file_path)
    if not filepath.exists():
        raise FileNotFoundError(f"ResolutionLevel target file missing: {str(filepath)}")

    with open(filepath, 'r') as file:
        filedata = file.read()

    original_data = filedata
    replacements = [
        (re.compile(r"SimPar@MaxMeshLevel\s*=\s*([0-9]+)"),
         f"SimPar@MaxMeshLevel = {resolution_level}"),
        (re.compile(r"Temper@MGMaxLev\s*=\s*([0-9]+)"),
         f"Temper@MGMaxLev = {resolution_level}"),
        (re.compile(r"Temper@MGMinLev\s*=\s*([0-9]+)"),
         f"Temper@MGMinLev = {lower_resolution_level}"),
        (re.compile(r"Temper@MGMedLev\s*=\s*([0-9]+)"),
         f"Temper@MGMedLev = {lower_resolution_level}")
    ]

    replaced_keys = []
    for pattern, replacement in replacements:
        filedata, num_replaced = re.subn(pattern, replacement, filedata)
        if num_replaced > 0:
            replaced_keys.append(replacement.split("=", 1)[0].strip())

    if filedata != original_data:
        print(f"Setting resolution level {resolution_level} in {str(filepath)} ({', '.join(replaced_keys)})")
        with open(filepath, 'w') as file:
            file.write(filedata)
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
        else:
            backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU_test.dat")
    else:
        if (paramDict['dieSimulation']) :
            backupDataFile = Path("_data_BU") / Path("q2p1_paramV_DIE.dat")
        else:
            backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU.dat")
    
#    backupDataFile = Path("_data_BU") / Path("q2p1_paramV_BU.dat")
    destDataFile = Path("_data") / Path("q2p1_param.dat")

    shutil.copyfile(str(backupDataFile), str(destDataFile))
    paramDict['partitionFormat'] = parsePartitionFormat(str(destDataFile))
    paramDict['recursivePartitioning'] = True
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
    detected_nodes = detect_node_count()
    if detected_nodes <= 0:
        detected_nodes = 1
        print("Could not detect node count automatically; defaulting to a single group.")
    else:
        print("Detected %d compute nodes for partitioning." % detected_nodes)
    partitionerParameters[1] = detected_nodes

    print("Partitioner parameters: ",-1, partitionerParameters)
    try:
        myLog.updateStatusLine("CurrentStatus=running Partitioner")
        partitioner.partition(
            paramDict['numProcessors']-1,
            partitionerParameters[0],
            partitionerParameters[1],
            "NEWFAC",
            "_data/meshDir/file.prj",
            partition_format=paramDict.get('partitionFormat', 'legacy')
        )
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
#                Compute periodic simulation start and iteration count
#===============================================================================
def calcSimLoopWindow():
    start_angle = 0.0
    nmax = calcMaxSimIterations()

    if paramDict['singleAngle'] >= 0.0:
        return start_angle, nmax

    start_angle = float(paramDict.get('startAngle', 0.0))
    span = 360.0 / float(paramDict['periodicity'])
    tolerance = 1e-9

    if start_angle >= span - tolerance:
        return start_angle, 0

    remaining_span = span - start_angle
    nmax = int(math.ceil(remaining_span / paramDict['deltaAngle']))

    return start_angle, nmax
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
#                Configure MPI environment variables
#===============================================================================
def configureMPIEnvironment():
    """
    Ensure the MPI file I/O environment matches the settings required
    for the die simulations, replacing the shell logic in RunnerGenDIE.
    """
    print("Configuring MPI environment variables for ROMIO:")
    for var, desired in MPI_ENV_DEFAULTS.items():
        prev = os.environ.get(var)
        os.environ[var] = desired
        if prev is None:
            status = "set"
        elif prev == desired:
            status = "kept"
        else:
            status = f"overridden (was {prev})"
        print(f"  {var}={desired} [{status}]")
#===============================================================================
#===============================================================================

#===============================================================================
#                The simulation loop for velocity calculation
#===============================================================================
def simLoopVelocity(workingDir):
    start, nmax = calcSimLoopWindow()

    mpiPath = paramDict['mpiCmd']
    numProcessors = paramDict['numProcessors']

    nmin = 0
    maxInnerIters = parseMaxNumSteps("_data/q2p1_param.dat")
    with open("_data/Extrud3D_0.dat", "a") as f:
        f.write("\n[E3DSimulationSettings]\n")
        f.write("dAlpha=" + str(paramDict['deltaAngle']) + "\n")
        f.write("Periodicity=" + str(paramDict['periodicity']) + "\n")
        f.write("nSolutions=" + str(paramDict['timeLevels']) + "\n")
        if paramDict.get('startAngle', 0.0) > 0.0:
            f.write("StartAngle=" + str(paramDict['startAngle']) + "\n")

    if nmax <= 0:
        print("No angle positions left to compute for this stage.")
        return 0

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
            
            if paramDict['rankFile'] == "":
             rankfileCommand = " "
            else:
             rankfileCommand = " -r " + paramDict['rankFile'] + " "

            if paramDict['useSrun']:
                launchCommand = "srun " + os.getcwd() + "/q2p1_sse"
                if paramDict['singleAngle'] >= 0.0:
                    launchCommand = launchCommand + " -a %d" %(angle)
            else:
                launchCommand = "mpirun -np " + str(numProcessors) + rankfileCommand + os.getcwd() + "/q2p1_sse"
                if paramDict['singleAngle'] >= 0.0 :
                    launchCommand = launchCommand + " -a %d" %(angle)

            exitCode = subprocess.call([launchCommand], shell=True)

            if paramDict['retryDeformation'] and exitCode == 55:
                with open("_data/q2p1_param.dat", "r") as f:
                  for line in f:
                    if "SimPar@UmbrellaStepM" in line:
                        orig_umbrella = int(line.split()[2])
                UmbrellaStepM = orig_umbrella
                while exitCode == 55 and UmbrellaStepM != 0:
                    replace_in_file("_data/q2p1_param.dat", "SimPar@UmbrellaStepM = "+str(UmbrellaStepM), "SimPar@UmbrellaStepM = "+str(int(UmbrellaStepM/2)))
                    UmbrellaStepM = int(UmbrellaStepM / 2)
                    print("Retrying deformation with UmbrellaStepsM = %d" % UmbrellaStepM)
                    exitCode = subprocess.call([launchCommand], shell=True)
                replace_in_file("_data/q2p1_param.dat", "SimPar@UmbrellaStepM = "+str(UmbrellaStepM), "SimPar@UmbrellaStepM = "+str(orig_umbrella))
                if UmbrellaStepM == 0:
                    print("UmbrellaStepsM reduced to 0 during retry")

        # Here the observer can be turned off
        observer.stop()

        if exitCode == 88:
          myLog.logErrorExit("CurrentStatus=the screw could not be created: wrong angle", exitCode)

        if exitCode == 57:
          if paramDict['temperature'] or paramDict['yamlConfig']:
            print("Warning: momentum solver did not converge (exit 57); continuing")
            exitCode = 0
          else:
            myLog.logErrorExit("CurrentStatus=run out of iterations", exitCode)
          
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
            os.remove("_data/prot.txt")

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
#                The simulation loop for coupled XSE velocity/temperature runs
#===============================================================================
def simLoopTemperatureXSE(workingDir):

    global stage_counter
    numProcessors = paramDict['numProcessors']
    mpiPath = paramDict['mpiCmd']
    maxIterations = 2
    for iter in range(maxIterations):
        angle = paramDict['singleAngle'] if paramDict['singleAngle'] >= 0.0 else 0.0
     
        if paramDict['shortTest']:
            backupVeloFile = XSE_PARAM_DIR / Path("q2p1_paramV_%01d_test.dat" % iter)
        else:
            backupVeloFile = XSE_PARAM_DIR / Path("q2p1_paramV_%01d.dat" % iter)
            
        if paramDict['shortTest']:
            backupTemperatureFile = XSE_PARAM_DIR / Path("q2p1_paramT_test.dat")
        else:
            backupTemperatureFile = XSE_PARAM_DIR / Path("q2p1_paramT_%01d.dat" % iter)
            
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

        archive_stage_outputs(stage_counter, include_heat=True)
        stage_counter = stage_counter + 1

    if paramDict['shortTest']:
        backupVeloFile = XSE_PARAM_DIR / Path("q2p1_paramV_%01d_test.dat" % maxIterations)
    else:
        backupVeloFile = XSE_PARAM_DIR / Path("q2p1_paramV_%01d.dat" % maxIterations)
        
    veloDestFile = Path("_data") / Path("q2p1_param.dat")
    print("Copying: ", backupVeloFile, veloDestFile)
    shutil.copyfile(str(backupVeloFile), str(veloDestFile))

    exitCode = simLoopVelocity(workingDir)
    archive_stage_outputs(stage_counter, include_heat=False)
    stage_counter = stage_counter + 1
#===============================================================================


#===============================================================================
#                Coupled loop for DIE temperature simulations
#===============================================================================
def simLoopTemperatureDIE(workingDir):

    global stage_counter
    numProcessors = paramDict['numProcessors']
    mpiPath = paramDict['mpiCmd']
    maxIterations = 2

    for iter in range(maxIterations):
        angle = paramDict['singleAngle'] if paramDict['singleAngle'] >= 0.0 else 0.0
        backupVeloFile = DIE_PARAM_DIR / Path("q2p1_paramV_DIE_%01d.dat" % iter)
        backupTemperatureFile = DIE_PARAM_DIR / Path("q2p1_paramT_DIE_%01d.dat" % iter)

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

        archive_stage_outputs(stage_counter, include_heat=True)
        stage_counter = stage_counter + 1

    backupVeloFile = DIE_PARAM_DIR / Path("q2p1_paramV_DIE_2.dat")
    veloDestFile = Path("_data") / Path("q2p1_param.dat")
    print("Copying: ", backupVeloFile, veloDestFile)
    shutil.copyfile(str(backupVeloFile), str(veloDestFile))

    exitCode = simLoopVelocity(workingDir)
    archive_stage_outputs(stage_counter, include_heat=False)
    stage_counter = stage_counter + 1
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
        opts, args = getopt.getopt(sys.argv[1:], 'n:f:p:d:a:c:r:t:smxhovury:',
                                   ['num-processors=', 'project-folder=',
                                    'periodicity=', 'delta-angle=', 'angle=',
                                    'host-conf=', 'rank-file=', 'time=', 'skip-setup','die-simulation',
                                    'skip-simulation','short-test', 'help',
                                    'do-temperature','version', 'use-srun',
                                    'retry-deformation', 'only-mesh-creation',
                                    'mesh-reduction','yaml='])

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
            paramDict['hasPeriodicity'] = True
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
        elif opt in ('--only-mesh-creation'):
            paramDict['onlyMeshCreation'] = True
        elif opt in ('-y', '--yaml'):
            paramDict['yamlConfig'] = arg
        else:
            usage()
            sys.exit(2)

    if paramDict['projectFolder'] == "":
        print("Error: no project folder specified.")
        usage()
        sys.exit(2)

    if paramDict['yamlConfig'] == "":
        print("Error: no YAML configuration specified (use -y/--yaml).")
        sys.exit(2)

    with open(paramDict['yamlConfig'], 'r') as cfg:
        plan = yaml.safe_load(cfg)
    if plan is None:
        plan = {}
    apply_yaml_options(plan)
    validate_yaml_stage_inputs(plan, paramDict['yamlConfig'])

    if (paramDict['hasDeltaAngle'] and paramDict['hasTimeLevels']):
        print("Error: Specifying both deltaAngle and timeLevels at the same time is error-prone and therefore prohibited.")
        sys.exit(2)
        
    if (paramDict['singleAngle'] >= 0.0 and  paramDict['temperature']) :
        print("Error: Specifying both singleAngle and Temperature Simulation at the same time is prohibited.")
        sys.exit(2)
        
    if (paramDict['dieSimulation']) :
        print("Switching to 'DIE' simulation !")
        paramDict['singleAngle'] = 0
        paramDict['periodicity'] = 36

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
    configureMPIEnvironment()

    e3dSetupDict = e3dToDict(projectFile)

    if not paramDict['hasTimeLevels']:
        if "SimodSetting" in e3dSetupDict:
            if "time_levels" in e3dSetupDict['SimodSetting']:
                if e3dSetupDict['SimodSetting']['time_levels'].isnumeric():
                    paramDict['timeLevels'] = int(e3dSetupDict['SimodSetting']['time_levels'])
                    paramDict['hasTimeLevels'] = True
                else:
                    raise TypeError("e3d.setup ['SimodSettings']['time_levels'] is not a numeric entry. Pls enter a number > 0.")


    if not paramDict['skipSetup']:
        if paramDict['onlyMeshCreation']:
            exitCode = onlyMeshCreation(workingDir, projectFile, projectPath, projectFolder)
            return
        else:
            exitCode = simulationSetup(workingDir, projectFile, projectPath, projectFolder)

        apply_resolution_level_to_param_file(Path("_data") / Path("q2p1_param.dat"))


    if paramDict['skipSimulation']:
        sys.exit()

    print("  ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("                                        Launching E3D...                                           ")
    print("  ----------------------------------------------------------------------------------------------------- ")
    print("  ")
    print("  ")

    stages = plan.get('stages')
    if 'main' in stages:
        stages['main'].setdefault('loop', 1)

    execute_yaml_plan(stages, workingDir)
    cleanWorkingDir(workingDir)

#===============================================================================
#                             Main Boiler Plate
#===============================================================================
if __name__ == "__main__":
    main()
    myLog.writeExitMsg()
