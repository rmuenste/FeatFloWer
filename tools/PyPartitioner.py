#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from shutil import copy
from part import *

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
  for iPar in xrange(nParFiles):
    copy(myParFiles[iPar],os.path.join(workPath,myParNames[iPar]+".par"))
  # Erzeuge zusätzliche Unterverzeichnisse falls Untergitter erzeugt werden sollen
  for i in xrange(1,nSubMesh+1):
    mkdir(os.path.join(workPath,"sub%03d" % i))
  # Bestimme, ob die Untergitter in umgekehrter Reihenfolge abgespeichert werden sollen.
  if pMethod in (11,12,13):
    bReversed=True
    pMethod-=10
  else:
    bReversed=False
  if nSubMesh > 1:
    # Lese Gitter ein
    myGrid=GetGrid(os.path.join(workPath,"GRID.tri"))
    # Erzeuge Nachbarschaftsinformationen für das Gitter
    myNeigh=GetNeigh(myGrid)
    # Lese Parametrisierungen und Ränder ein
    myParTypes=[]
    myParameters=[]
    myBoundaries=[]
    for iPar in xrange(nParFiles):
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
    copy(myGridFile,os.path.join(workPath,"sub001","GRID.tri"))
    copy(ProjektFile,os.path.join(workPath,"sub001","GRID.prj"))
    for iPar in xrange(nParFiles):
      copy(myParFiles[iPar],os.path.join(workPath,"sub001",myParNames[iPar]+".par"))
  # Im Grunde "kSubPart=int(math.ceil(nnPart/float(nSubMesh)))"
  kSubPart= nnPart//nSubMesh if nnPart%nSubMesh==0 else nnPart//nSubMesh+1
  iPart=0
  rIter = xrange(nSubMesh,0,-1) if bReversed else xrange(1,nSubMesh+1)
  for i in rIter:
    subPath=os.path.join(workPath,"sub%03d"%i)
    myGrid=GetGrid(os.path.join(subPath,"GRID.tri"))
    myNeigh=GetNeigh(myGrid)
    myParTypes=[]
    myParameters=[]
    myBoundaries=[]
    for iPar in xrange(nParFiles):
      ParName=os.path.join(subPath,myParNames[iPar]+".par")
      (ParType,Parameter,Boundary)=GetPar(ParName,myGrid[1])
      myParTypes.append(ParType)
      myParameters.append(Parameter)
      myBoundaries.append(Boundary)
    nPart=min(iPart+kSubPart,nnPart)-iPart
    # Partitionierung mittel verschiedener Methoden (WIP)    
    if pMethod in (1,2,3):
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
  print "The partitioning was successful!"

# Dokumentation
__doc__ = \
"""
Partition a grid for multiprocessing.
Calling convention: ./PyPartitioner.py NPart PartMethod NSubPart MeshName ProjektFile
Exanple: ./PyPartitioner.py 12 1 1 NEWFAC _adc/2D_FAC/2Dbench.prj
"""
# Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
#mkdir("_adc")
mkdir("_mesh")

# Format-Check der Parameter
if not(len(sys.argv)==6 and sys.argv[1].isdigit() and sys.argv[3].isdigit()):
  sys.exit(__doc__)
if not os.path.exists(sys.argv[5]):
  sys.exit("Projekt file '%s' does not exist!" % sys.argv[5])

# Sanity-Check der Parameter
NPart=int(sys.argv[1])
PartMethod=int(sys.argv[2])
NSubPart=int(sys.argv[3])
if NPart <1:
  sys.exit("Number of Partitions has to be >=1 !")
if NSubPart<1:
  sys.exit("There has to be at least one subgrid!")
if not (PartMethod in (1,2,3,11,12,13) or str(-PartMethod).strip("123")==""):
  sys.exit("Only integer numbers 1,2,3 (+10) or negative numbers containing the digits 1,2,3 are valid partitioning methods!")
if PartMethod<0 and NSubPart==1:
  sys.exit("Partitionig Method %d requires more than 1 subgrid!"%PartMethod)
MeshName=sys.argv[4]
ProjektFile=sys.argv[5]

# Aufruf der Hauptroutine
MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
