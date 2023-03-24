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
    print("[-h, --help]: prints this message")


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

    for ilevel in range(2, 2 + testLevels):

        # copy the data file to the default location and set the level
        moveAndSetLevel(dataFile, "_data/q2p1_param.dat", ilevel)

        # call the specified binary using numProcessors processes 
        subprocess.call(['mpirun -np %i ./%s' %(numProcessors, binaryName)], shell=True)


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


    for idx, item in enumerate(allTables):
        jsonData['DashBoardVisualization']['DashBoardTable']['Tables'][idx]['data']['rows']=rowsArray

        for entry in jsonData['DashBoardVisualization']['DashBoardTable']['Tables'][idx]['data']['cols']:
            if entry['label'] == "Time[s]":
                entry['label'] = entry['label'] + ' ' + str(numProcessors) + 'P' 

    for idx, item in enumerate(Diagramms):
        generator = item['DataGenerator']
        rows = get_col_data("_data/prot.txt", generator)
        jsonData['DashBoardVisualization']['DashBoardDiagramm']['Diagramms'][idx]['data']['rows']=rows


    print(jsonData)

    timeStamp = round(time.time())

    now = datetime.now()
    jsonData["DashBoardVisualization"]["timeStamp"] = timeStamp 
    jsonData["DashBoardVisualization"]["date"] = now.strftime("%d/%m/%y %H:%M:%S") 


    finalName = f"{noteName}-{timeStamp}"

    headers = {
        'Content-Type': 'application/json'
    }

    response = requests.post(url, data=json.dumps(jsonData), headers=headers)

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
