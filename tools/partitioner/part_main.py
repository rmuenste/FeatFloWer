#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from shutil import copy
from .part import *

# Dokumentation
__doc__ = \
"""
Partition a grid for multiprocessing.
Calling convention: ./PyPartitioner.py NPart PartMethod NSubPart MeshName ProjektFile
Example: ./PyPartitioner.py 12 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
"""

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Hilfsroutine
def mkdir(dir):
    """
    Erzeugt nur dann das Verzeichnis "dir", wenn es noch nicht vorhanden ist.
    Falls eine Datei dieses Namens existieren sollte, wird sie durch das Verzeichnis ersetzt.
    """
    if os.path.exists(dir):
        if os.path.isdir(dir):
            return
        else:
            os.remove(dir)
    os.mkdir(dir)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def checkParameters(params):
    """
    Diese Funktion ueberprueft die Eingabeparameter auf Gueltigkeit
    """
    # Format-Check der Parameter
    if not(len(params)==6 and params[1].isdigit() and params[3].isdigit()):
        sys.exit(__doc__)
    if not os.path.exists(params[5]):
        sys.exit("Projekt file '%s' does not exist!" % params[5])

    # Sanity-Check der Parameter
    NPart=int(params[1])
    PartMethod=int(params[2])
    NSubPart=int(params[3])
    if NPart <1:
        sys.exit("Number of Partitions has to be >=1 !")
    if NSubPart<1:
        sys.exit("There has to be at least one subgrid!")
    if not (PartMethod in (1,2,3,11,12,13) or str(-PartMethod).strip("123") == ""):
        sys.exit("Only integer numbers 1,2,3 (+10) or negative numbers containing " +
                 "the digits 1,2,3 are valid partitioning methods!")
    if PartMethod<0 and NSubPart==1:
        sys.exit("Partitionig Method %d requires more than 1 subgrid!"%PartMethod)
    MeshName=params[4]
    ProjektFile=params[5]

    # Rueckgabe der Parameter
    return NPart, PartMethod, NSubPart, MeshName, ProjektFile

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Hauptroutine
def MainProcess(nnPart,pMethod,nSubMesh,MeshName,ProjektFile):
    # Lese Projektdatei und extrahiere Gitter und Parameterdateinamen
    (nParFiles,myGridFile,myParFiles,myParNames)=GetFileList(ProjektFile)

    # Erzeuge übergeordnetes Verzeichnis, in dem gearbeitet werden soll
    workPath=os.path.join("_mesh",MeshName)
    mkdir(workPath)

    # Kopiere alle benötigten Dateien in dieses Verzeichnis
    copy(myGridFile,os.path.join(workPath,"GRID.tri"))
    copy(ProjektFile,os.path.join(workPath,"GRID.prj"))
    for iPar in range(nParFiles):
        copy(myParFiles[iPar],os.path.join(workPath,myParNames[iPar]+".par"))

    # Erzeuge zusätzliche Unterverzeichnisse falls Untergitter erzeugt werden sollen
    for i in range(1,nSubMesh+1):
        mkdir(os.path.join(workPath,"sub%04d" % i))

    # Bestimme, ob die Untergitter in umgekehrter Reihenfolge abgespeichert werden sollen.
    if pMethod in (11,12,13):
        bReversed=True
        pMethod-=10
    else:
        bReversed=False

    # Speziller Marker für atomares Splitting, der nötig wird, wenn nSubMesh>1 ist
    # und nnPart == #Gitterzellen des Hauptgitters
    bAtomicSplitting=False
    # Falls mit Untergittern gearbeitet werden soll, so spalte das Hauptgitter auf,
    # ansonsten nehme das Hauptgitter als Untergitter Nr.1
    if nSubMesh > 1:
        # Lese Gitter ein
        myGrid=GetGrid(os.path.join(workPath,"GRID.tri"))
        # (Anzahl der Gitterzellen == nnPart) => aktiviere atomares Splitting
        if myGrid[0]==nnPart:
            bAtomicSplitting=True
        # Erzeuge Nachbarschaftsinformationen für das Gitter
        myNeigh=GetNeigh(myGrid)
        # Lese Parametrisierungen und Ränder ein
        myParTypes=[]
        myParameters=[]
        myBoundaries=[]

        for iPar in range(nParFiles):
            ParName=os.path.join(workPath,myParNames[iPar]+".par")
            (ParType,Parameter,Boundary)=GetPar(ParName,myGrid[1])
            myParTypes.append(ParType)
            myParameters.append(Parameter)
            myBoundaries.append(Boundary)

        # Aufspaltung in Untergitter
        if pMethod in (1,2,3):
            myPart=GetParts(myNeigh,nSubMesh,pMethod)
        else:
            try:
                myPart=PartitionAlongAxis(myGrid,nSubMesh,pMethod)
            except AssertionError as ErrorInst:
                sys.exit("Error at creating subgrids along axis: %s"%ErrorInst)
            pMethod=1

        # Schreibe die Gitter und Parametrisierungen der einzelnen Rechengebiete
        myParam=(myParNames,myParTypes,myParameters,myBoundaries)
        GetSubs(workPath,myGrid,nSubMesh,myPart,myNeigh,nParFiles,myParam,False)
    elif nSubMesh==1:
        copy(myGridFile,os.path.join(workPath,"sub0001","GRID.tri"))
        copy(ProjektFile,os.path.join(workPath,"sub0001","GRID.prj"))
        for iPar in range(nParFiles):
            copy(myParFiles[iPar],os.path.join(workPath,"sub0001",myParNames[iPar]+".par"))

    # Im Grunde "kSubPart=int(math.ceil(nnPart/float(nSubMesh)))"
    kSubPart= nnPart//nSubMesh if nnPart%nSubMesh==0 else nnPart//nSubMesh+1
    iPart=0
    rIter = range(nSubMesh,0,-1) if bReversed else range(1,nSubMesh+1)
    for i in rIter:
        subPath=os.path.join(workPath,"sub%04d"%i)
        myGrid=GetGrid(os.path.join(subPath,"GRID.tri"))
        myNeigh=GetNeigh(myGrid)
        myParTypes=[]
        myParameters=[]
        myBoundaries=[]

        for iPar in range(nParFiles):
            ParName = os.path.join(subPath,myParNames[iPar] + ".par")
            (ParType, Parameter, Boundary)=GetPar(ParName, myGrid[1])
            myParTypes.append(ParType)
            myParameters.append(Parameter)
            myBoundaries.append(Boundary)

        nPart = min(iPart+kSubPart, nnPart) - iPart
        # Partitionierung mittel verschiedener Methoden (WIP)
        if pMethod in (1,2,3):
            if bAtomicSplitting:
                myPart=GetAtomicSplitting(len(myNeigh))
                nPart=max(myPart)
            else:
                myPart=GetParts(myNeigh,nPart,pMethod)
        else:
            sys.exit("Partitioning method %d is not available for subgrids!"%pMethod)
        # Schreibe die Gitter und Parametrisierungen der einzelnen Rechengebiete
        myParam=(myParNames,myParTypes,myParameters,myBoundaries)
        GetSubs(subPath,myGrid,nPart,myPart,myNeigh,nParFiles,myParam,True)
        iPart+=nPart
        if iPart==nnPart:
            break
        pass

    print("The partitioning was successful!")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def partition(NPart, PartMethod, NSubPart, MeshName, ProjektFile):

    # Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
    mkdir("_mesh")

    # Aufruf der Hauptroutine
    MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
