#!/usr/bin/env python
# vim: set filetype=python
"""
A python launcher script for a FeatFloWer application
"""
import os

import xml.etree.ElementTree as ET
import sys
import getopt
import platform
import subprocess
import requests
import re
import json
import partitioner
import pprint
import time
from datetime import datetime

url = 'http://127.0.0.1:3000/api/user/1'

#===============================================================================
#                      Function: Usage
#===============================================================================
def usage():
    print("Usage: configure [options]")
    print("Where options can be:")
    print("[-h, --help]: prints this message")

#===============================================================================
#                      Function: generateSlurmScript
#===============================================================================
def generateSlurmScript():
    partition = "med"
    constraint = "[bttf]"
    nodes = "1"
    ntasksPerNode = "16"
    time = "08:00:00"
    mem = "5G"
    name = "FF-Bench"
    executable = "q2p1_fc_ext"
    slurmString = f"""#!/bin/bash
#SBATCH --partition={partition}
#SBATCH --constraint={constraint}
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={ntasksPerNode}
#SBATCH --time={time}
#SBATCH --mem-per-cpu={mem}
#SBATCH --job-name={name}
shopt -s expand_aliases
source ~/.bashrc
mpirun -np {ntasksPerNode} ./{executable} 
"""
    return slurmString
#===============================================================================
#                      Function:  submitAndObserve
#===============================================================================
def submitAndObserve():
    # Submit job using sbatch
    slurmString = generateSlurmScript()

    sbatch_command = 'sbatch myjob.sh'
    process = subprocess.Popen(sbatch_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f'Error submitting job: {stderr.decode("utf-8")}')
        exit(1)
    job_id = stdout.decode("utf-8").strip()  # get job ID from sbatch output

    # Query job status using sacct
    sacct_command = f'sacct --format=State --jobs={job_id}'
    interval_seconds = 5
    while True:
        process = subprocess.Popen(sacct_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            print(f'Error querying job status: {stderr.decode("utf-8")}')
            exit(1)
        job_status = stdout.decode("utf-8").strip().split('\n')[1].split()[0]  # get job status from sacct output
        print(f'Job {job_id} status: {job_status}')
        if job_status == 'COMPLETED':
            print('Job completed successfully')
            break
        elif job_status in ['FAILED', 'CANCELLED']:
            print('Job failed or was cancelled')
            break
        time.sleep(interval_seconds)  # wait for interval_seconds before checking job status again

#===============================================================================
#                      Function:  submitAndObserveSync
#===============================================================================
def submitAndObserveSync():

    slurmString = generateSlurmScript()
    with open("myjob.sh", "w") as f:
        f.write(slurmString)

    # Submit job using sbatch
    sbatch_command = 'sbatch myjob.sh'

    try:
        output = subprocess.check_output(sbatch_command.split(), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        print(f'Error submitting job: {e.output.decode("utf-8")}')
        exit(1)
    job_id = output.decode("utf-8").strip().split()[-1]  # get job ID from sbatch output

    # Query job status using sacct
    sacct_command = f'sacct --format=State --jobs={job_id}'
    interval_seconds = 10
    while True:
        try:
            output = subprocess.check_output(sacct_command.split(), stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print(f'Error querying job status: {e.output.decode("utf-8")}')
            exit(1)
        #job_status = output.decode("utf-8").strip().split('\n')[1].split()[0]  # get job status from sacct output
        job_status = output.decode("utf-8").strip().split('\n')  # get job status from sacct output
        job_status = [x.strip() for x in job_status]
        if len(job_status) < 3:
            job_status = "UNINITIALIZED"
        else:
            job_status = job_status[2]

        print(f'Job {job_id} status: {job_status}')
        if job_status == 'COMPLETED':
            print('Job completed successfully')
            break
        elif job_status in ['FAILED', 'CANCELLED']:
            print('Job failed or was cancelled')
            break
        else:
            print('Job not FAILED, CANCELLED or COMPLETED')
        time.sleep(interval_seconds)  # wait for interval_seconds before checking job status again        
    return job_status

#===============================================================================
#                      Function:  moveAndSetLevel
#===============================================================================
def moveAndSetLevel(fileIn, fileOut, level):
    maxLevelStr = "SimPar@MaxMeshLevel = " + str(level)
    with open(fileOut, "w") as newFile:
        with open(fileIn, "r") as file:
            for line in file:
                newLine = re.sub(r"^[\s]*SimPar@MaxMeshLevel[\s]*=(\s | \w)*", maxLevelStr, line)
                newFile.write(newLine)


#===============================================================================
#                      Function: readTestConfiguration
#===============================================================================
def readTestConfiguration(fileName):
    print("TestXML: " + fileName)

    tree = ET.parse(fileName)

    root = tree.getroot()

    myDict = {}

    for detailedConf in root.findall('DetailedConfiguration'):
        print(detailedConf.tag, detailedConf.attrib)
        for k, v in detailedConf.attrib.items():
            myDict[k] = v
#            print("Key : {0}, Value : {1}".format(k,v))

    for child in root.findall('DashBoardConfiguration'):
        print(child.tag, child.attrib, "bla")
        for c2 in child.iter('DashBoardAttributes'):
#            print(c2.tag, c2.attrib['benchName'])
            myDict['benchName'] = c2.attrib['benchName']
            for c3 in child.iter('DashBoardData'):
                print("DashBoardDataTag: ", c3.tag, c3.attrib)
                myDict['benchKeyword'] = c3.attrib['benchKeyword']
                indices = []
                for value in c3.iter('value'):
                    indices.append(int(value.text))
                myDict['lineIndices'] = indices
#                print(myDict['lineIndices'])
            for c3 in c2.iter('DashBoardVisualization'):
                print("C3 ", c3.tag, c3.attrib)
                for c4 in c3.iter('DashBoardTable'):
                    print("C4 ", c4.tag, c4.attrib)
                    for c5 in c4.iter('TableColumn'):
                        print("C5 ", c5.tag, c5.attrib)


    return myDict


#===============================================================================
#                      Function: get_col_data 
#===============================================================================
def get_col_data(file_name, dataGenerator):
  with open(file_name, "r") as sources:
    lines = sources.readlines()
  
  it_found = False
  t_found = False

  lineIndices = dataGenerator['DataValues']['value']
  
  col_data = []
  for line in lines:
    m = re.search(dataGenerator['benchKeyword'],line)
    if m != None:
      val = line.split()
      colArray = []
      colObject = {"c" : colArray}
      for idx in lineIndices:
          colArray.append({"v" : val[int(idx)]}) 
      col_data.append(colObject)

#      col_data.append({"c" : [ {"v" : val[int(x)] } ]} for x in lineIndices] })

  return col_data


#===============================================================================
#                      Function: Get time entry
#===============================================================================
def getTimeEntry(fileName, varName):
    with open(fileName, "r") as sources:
        lines = sources.readlines()

    tFound = False

    for line in reversed(lines):
        m = re.match(varName, line)
        if m != None:
            splitLine = line.split(':')
            val = splitLine[1].strip()
            tFound = True
            break

    if tFound:
        return val
    else:
        return 0


#===============================================================================
#                      Function: Get log entry
#===============================================================================
def getLogEntry(fileName, dataGenerator):
    with open(fileName, "r") as sources:
        lines = sources.readlines()

    tFound = False

    lineIndices = dataGenerator['DataValues']['value']

    for line in reversed(lines):
        m = re.search(dataGenerator['benchKeyword'], line)
        if m != None:
            splitLine = line.split(':')
            val = splitLine[1].strip()
            tFound = True
            break

    if tFound:
        val = line.split()
        print("Force entry ", val)
        result = [{"v" : val[int(x)]} for x in lineIndices]
        return result
    else:
        return 0

#===============================================================================
#                      Function: executeTest 
#===============================================================================
def executeTest(fileName):

    with open(fileName) as f:
        paramDict = json.loads(f.read())

    numProcessors = int(paramDict["TestConfiguration"]["DetailedConfiguration"]["processorCount"])
    testLevels = int(paramDict["TestConfiguration"]["DetailedConfiguration"]["testLevels"])
    dataFile = paramDict["TestConfiguration"]["DetailedConfiguration"]["testDataFile"]
    binaryName = paramDict["TestConfiguration"]["DetailedConfiguration"]["testBinaryName"]
    noteName = paramDict["TestConfiguration"]["DetailedConfiguration"]["testNoteName"]
    meshProjectFile = paramDict["TestConfiguration"]["DetailedConfiguration"]["meshProjectFile"]
    benchKeyword = paramDict["TestConfiguration"]["DashBoardConfiguration"]["DashBoardAttributes"]["DashBoardData"]["benchKeyword"]
    lineIndices = paramDict["TestConfiguration"]["DashBoardConfiguration"]["DashBoardAttributes"]["DashBoardData"]["DataValues"]["value"]

    lineIndices = [int(x) for x in lineIndices]

    DashBoardViz = paramDict["TestConfiguration"]["DashBoardConfiguration"]["DashBoardAttributes"]["DashBoardVisualization"]

    Diagramms = []

    allTables = []

    if "DashBoardDiagramm" in DashBoardViz:
        Diagramms = DashBoardViz["DashBoardDiagramm"]["Diagramms"]
        print(Diagramms)
    else:
        print("no diagramms found")

    if "DashBoardTable" in DashBoardViz: 
        allTables = DashBoardViz['DashBoardTable']['Tables']
    else:
        print("no tables found")

    rowsArray = [] 
    timeEntry = [] 
    jsonData = {}

    jsonData['DashBoardVisualization'] = paramDict["TestConfiguration"] \
                                                  ["DashBoardConfiguration"] \
                                                  ["DashBoardAttributes"] \
                                                  ["DashBoardVisualization"] 
    
    partitioner.partition(numProcessors-1, 1, 1, "NEWFAC", meshProjectFile)

    result = ""
    for ilevel in range(2, 2 + testLevels):

        # copy the data file to the default location and set the level
        moveAndSetLevel(dataFile, "_data/q2p1_param.dat", ilevel)

        # call the specified binary using numProcessors processes 
        #subprocess.call(['mpirun -np %i ./%s' %(numProcessors, binaryName)], shell=True)
        result = submitAndObserveSync()
        if result != "COMPLETED":
          break

        """ 
            append a row to the array of rows, the format of a row is in general:
         
            {
                "c": [                     # The "c" is a json array of the of all column entries in the row
                       {"v" : value},      # An element of the "c" array is an object that stores
                                           # the value of the particular column entry under "v"
                       ...
                ]
            }

        """ 

        """ TODO
            This actually needs to be not the allTables, but the table for the particular test
            meanging something like allTables['fac2dTest']

            TABLE AND DIAGRAMM DATA:

            What is already ok, is that the table data is generated in the level loop 
            so that the progression of the result with the level can be seen

            Probably the diagramm data should also be generated on each level
        """

        for idx, item in enumerate(allTables):

            rowList = []
            rowList.append({"v" : ilevel})

            generator = item['DataGenerator']
            logData = getLogEntry("_data/prot.txt", generator)

            for i in logData:
                rowList.append(i)

            # get the total time from the statistics file
            timeEntry = getTimeEntry("_data/Statistics.txt", " Overall time")
            timeEntry = timeEntry.split()

            rowList.append({"v": timeEntry[0][:-3]})

            # append the entry to the row array
            rowsArray.append({"c": rowList})


    if result == "COMPLETED":
        for idx, item in enumerate(allTables):
            jsonData['DashBoardVisualization']['DashBoardTable']['Tables'][idx]['data']['rows']=rowsArray

            for entry in jsonData['DashBoardVisualization']['DashBoardTable']['Tables'][idx]['data']['cols']:
                if entry['label'] == "Time[s]":
                    entry['label'] = entry['label'] + ' ' + str(numProcessors) + 'P' 

        for idx, item in enumerate(Diagramms):
            generator = item['DataGenerator']
            rows = get_col_data("_data/prot.txt", generator)
            jsonData['DashBoardVisualization']['DashBoardDiagramm']['Diagramms'][idx]['data']['rows']=rows



    timeStamp = round(time.time())

    now = datetime.now()

    if result == "COMPLETED":
        jsonData["DashBoardVisualization"]["timeStamp"] = timeStamp 
        jsonData["DashBoardVisualization"]["date"] = now.strftime("%d/%m/%Y %H:%M:%S") 

        jsonData["DashBoardVisualization"]["git"] = True 
        jsonData["DashBoardVisualization"]["cmake"] = True 
        jsonData["DashBoardVisualization"]["build"] = True 
        jsonData["DashBoardVisualization"]["runtime"] = True 
        print(jsonData)
    else:
        jsonData["DashBoardVisualization"]["timeStamp"] = timeStamp 
        jsonData["DashBoardVisualization"]["date"] = now.strftime("%d/%m/%Y %H:%M:%S") 

        jsonData["DashBoardVisualization"]["git"] = True 
        jsonData["DashBoardVisualization"]["cmake"] = True 
        jsonData["DashBoardVisualization"]["build"] = True 
        jsonData["DashBoardVisualization"]["runtime"] = False 
        print(jsonData)


    finalName = f"{noteName}-{timeStamp}"

    headers = {
        'Content-Type': 'application/json'
    }

    try:
        response = requests.post(url, data=json.dumps(jsonData), headers=headers)
        response.raise_for_status() # This line will raise an HTTPError if the status code is not in the 200 range
        print("Request successful")
        # Process the response here
    except requests.exceptions.HTTPError as err:
        print("Error: Request failed with status code", err.response.status_code)
        # Handle the error here
    except Exception as err:
        print("Error: An error occurred during the request:", err)
        # Handle any other errors here

    with open('%s.json' %finalName, 'w') as theFile:
        theFile.write(json.dumps(jsonData) + '\n')
    
#===============================================================================
#                              Main function
#===============================================================================
def main():
    """ Preliminary parameter list:
        ======================
        Simulation Parameters 
        ======================

        num-processors : Number of processors
        project-file : The path to the project file
        parameter-file : The path to the parameter file
        binary-name : The name of the executable 
        levels : The number of additional levels w.r.t the 
                 maxlevel given in the initial 
                 parameter-file that the test will be performed
                 on.
    
        ======================
         Dashboard Parameters 
        ======================
    
        benchName : The identifier of the benchmark 
        style : The visualization style of the test
        tableCaption : The caption if the table style is chosen 
        ouput-file : The path to the JSON output file 
    """
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'n:p:h', ['num-processors=', 'params=', 'help'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    params = ''
  
    numProcessors = 16
    
    for opt, arg in opts:
        if opt in ('-h', '--help'):
            usage()
            sys.exit(2)
        elif opt in ('-n', '--num-processors'):
            if arg.isdigit():
                numProcessors = int(arg)
            else:
                print("--num-processors= " + str(arg) + " is not a valid number.")
                sys.exit(2)
        elif opt in ('-p', '--params'):
            params = arg
        else:
            usage()
            sys.exit(2)

    script_path = os.path.abspath(__file__)
    
    # get the directory that the script is located in
    script_dir = os.path.dirname(script_path)
    
    os.chdir(script_dir)

    potentialTests = os.listdir("tests")
    
    for item in potentialTests:
        m = re.search(r"^test-\w+.json$", item) 
        if m:
            print('Executing test: ' + item)
            executeTest("tests/" + item)


#===============================================================================
#                           Main "Boiler Plate"
#===============================================================================
if __name__ == "__main__":
    main()
