#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from shutil import copy
from .part import *

SUPPORTED_PARTITION_FORMATS = ("legacy", "json")

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
    raw_args = params[1:]
    positional = []
    partition_format = "legacy"
    for arg in raw_args:
        if arg.startswith("--"):
            if arg.startswith("--partition-format="):
                partition_format = arg.split("=", 1)[1].lower()
            else:
                sys.exit("Unknown option '%s'\n%s" % (arg, __doc__))
        else:
            positional.append(arg)

    if len(positional) != 5:
        sys.exit(__doc__)

    if partition_format not in SUPPORTED_PARTITION_FORMATS:
        sys.exit("Unsupported partition-format '%s' (expected one of %s)" %
                 (partition_format, ", ".join(SUPPORTED_PARTITION_FORMATS)))

    if not (positional[0].isdigit() and positional[2].isdigit()):
        sys.exit(__doc__)
    if not os.path.exists(positional[4]):
        sys.exit("Projekt file '%s' does not exist!" % positional[4])

    # Sanity-Check der Parameter
    NPart=int(positional[0])
    PartMethod=int(positional[1])
    NSubPart=int(positional[2])
    if NPart <1:
        sys.exit("Number of Partitions has to be >=1 !")
    if NSubPart<1:
        sys.exit("There has to be at least one subgrid!")
    if not (PartMethod in (1,2,3,11,12,13) or str(-PartMethod).strip("123") == ""):
        sys.exit("Only integer numbers 1,2,3 (+10) or negative numbers containing " +
                 "the digits 1,2,3 are valid partitioning methods!")
    if PartMethod<0 and NSubPart==1:
        sys.exit("Partitionig Method %d requires more than 1 subgrid!"%PartMethod)
    MeshName=positional[3]
    ProjektFile=positional[4]

    # Rueckgabe der Parameter
    return NPart, PartMethod, NSubPart, MeshName, ProjektFile, partition_format

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Hauptroutine
def MainProcess(nnPart,pMethod,nSubMesh,MeshName,ProjektFile,partition_format="legacy"):
    # Lese Projektdatei und extrahiere Gitter und Parameterdateinamen
    (nParFiles,myGridFile,myParFiles,myParNames)=GetFileList(ProjektFile)

    # Erzeuge übergeordnetes Verzeichnis, in dem gearbeitet werden soll
    workPath=os.path.join("_mesh",MeshName)
    mkdir(workPath)
    init_partition_output(workPath, partition_format)

    # Kopiere alle benötigten Dateien in dieses Verzeichnis
    gridPath=os.path.join(workPath,"GRID.tri")
    if partition_format == "legacy":
        copy(myGridFile, gridPath)
    copy(ProjektFile,os.path.join(workPath,"GRID.prj"))
    for iPar in range(nParFiles):
        targetParPath = os.path.join(workPath,myParNames[iPar]+".par")
        if partition_format == "legacy":
            copy(myParFiles[iPar], targetParPath)

    base_grid_source = gridPath if partition_format == "legacy" else myGridFile
    baseGrid = GetGrid(base_grid_source)
    cache_grid_data(gridPath, baseGrid)
    if partition_format == "json":
        OutputGrid(gridPath, baseGrid, sourceName="GRID.tri")

    baseParCache = []
    for iPar in range(nParFiles):
        parPath = os.path.join(workPath,myParNames[iPar]+".par")
        par_read_path = parPath if partition_format == "legacy" else myParFiles[iPar]
        parData = GetPar(par_read_path, baseGrid[1])
        baseParCache.append(parData)
        cache_par_data(parPath, parData)
        if partition_format == "json":
            OutputParFile(
                parPath,
                parData[0],
                parData[1],
                sorted(list(parData[2])),
                source_name="%s.par" % myParNames[iPar],
            )

    # Erzeuge zusätzliche Unterverzeichnisse falls Untergitter erzeugt werden sollen
    if partition_format == "legacy":
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
        myGrid=baseGrid
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
        GetSubs(workPath,myGrid,nSubMesh,myPart,myNeigh,nParFiles,myParam,False,sourceName="GRID.tri")
    elif nSubMesh==1:
        subPathBase = os.path.join(workPath,"sub0001")
        if partition_format == "legacy":
            copy(myGridFile,os.path.join(subPathBase,"GRID.tri"))
            copy(ProjektFile,os.path.join(subPathBase,"GRID.prj"))
            for iPar in range(nParFiles):
                copy(myParFiles[iPar],os.path.join(subPathBase,myParNames[iPar]+".par"))
        else:
            subGridPath = os.path.join(subPathBase, "GRID.tri")
            OutputGrid(subGridPath, baseGrid, sourceName="GRID.tri")
            cache_grid_data(subGridPath, baseGrid)
            for iPar in range(nParFiles):
                parName = "%s.par" % myParNames[iPar]
                parPath = os.path.join(subPathBase, parName)
                parType, parParam, parBoundary = baseParCache[iPar]
                OutputParFile(parPath, parType, parParam, sorted(list(parBoundary)),
                              source_name=parName)
                cache_par_data(parPath, (parType, parParam, parBoundary))

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
        GetSubs(subPath,myGrid,nPart,myPart,myNeigh,nParFiles,myParam,True,sourceName="GRID.tri")
        iPart+=nPart
        if iPart==nnPart:
            break
        pass

    print("The partitioning was successful!")
    finalize_partition_output()

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
def partition(NPart, PartMethod, NSubPart, MeshName, ProjektFile, partition_format="legacy"):

    # Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
    mkdir("_mesh")

    # Aufruf der Hauptroutine
    MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile,partition_format=partition_format)
