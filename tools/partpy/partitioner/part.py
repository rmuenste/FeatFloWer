#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Alles was man zum Laden von Metis benötigt
from ctypes import CDLL, c_int, POINTER, byref

#from future.standard_library import install_aliases
#install_aliases()

from functools import reduce

from six.moves import zip

from itertools import repeat, count

from collections import Counter

from math import sqrt
import os
import sys
from shutil import copy

metis=None
metis_func=[]

#Kleine private Hilfsroutinen
def _readAfterKeyword(fh,keyword):
  s=fh.readline()
  while keyword not in s:
    s=fh.readline()

def _try_in_place_first(name):
  tmp=os.path.join(os.curdir,name)

  if not os.path.exists(tmp):
    tmp=os.path.join(os.curdir, "../lib64", name)

  if not os.path.exists(tmp):
    tmp=name

  try:
    return CDLL(tmp)
  except OSError:
    print("An error of type OSError occurred while trying to find the metis library:")
    print("The metis library was neither found in the current folder nor in the system library path.")
  except:
    print("An error occurred loading the metis library:")
    print("The metis library was neither found in the current folder nor in the system library path.")

# Ab hier kommen die öffentlichen Funktionen des Moduls

def GetFileList(cProjName):
  """
  Auslesen der Projektdatei. Laden des Gitternamens und der Namen der Parametrisierungsdateien.
  """
  nPar=0
  myGridFile=""
  myParFiles=[]
  myParNames=[]
  myParExts=[]
  ProjektDir=os.path.dirname(cProjName)
  fProjekt=open(cProjName,'r')
  for s in fProjekt:
    s=s.strip()
    if '.' in s:
      (prefix,sep,ext)=s.rpartition('.')
      ext=ext.lower()
      if ext=="tri":
        myGridFile=os.path.join(ProjektDir,s)
      elif ext in ["pls","par"]:
        nPar+=1
        myParFiles.append(os.path.join(ProjektDir,s))
        myParNames.append(prefix)
        myParExts.append(ext)
  fProjekt.close()
  print("The projekt folder consists of the following files:")
  print("- Grid File:", myGridFile)
  print("- Boundary Files:")
  print("\n".join(map(lambda x: "  * %s" % x,myParFiles)))
  return (nPar,myGridFile,myParFiles,myParNames,myParExts)

def GetGrid(GridFileName):
    """
    Liest ein Gitter aus der Datei "GridFileName".
    Der Rückgabewert hat die Struktur: (NEL,NVT,Coord,KVert,Knpr)
    """
    print("Grid input file: '%s'" % GridFileName)
    f=open(GridFileName,'r')
    f.readline()
    f.readline()
    # Lese Anzahl der Zellen und Knoten
    g=f.readline().split()
    NEL=int(g[0])
    NVT=int(g[1])
    # Lese die Koordinaten
    _readAfterKeyword(f,"DCORVG")
    Coord=[]
    for i in range(NVT):
        g=f.readline().split()
        Coord.append(tuple(map(float,g)))
    # Lese die Zelldaten
    _readAfterKeyword(f,"KVERT")
    Kvert=[]
    for i in range(NEL):
        g=f.readline().split()
        Kvert.append(tuple(map(int,g)))
    # Lese die Randdaten
    _readAfterKeyword(f,"KNPR")
    g=f.read().split()
    Knpr=tuple(map(int,g))
    # Fertig
    #f.close()
    return (NEL,NVT,tuple(Coord),tuple(Kvert),Knpr)


def _GetPar(ParFileName):
    """
    Lese Randbeschreibungen aus einer Parameterdatei. Maximale Knotenzahl NVT
    bestimmt zudem die Länge der Randliste.
    Rückgabe: (Name des Randes, Daten des Randes, Boolsche Liste für alle Knoten)
    """
    print("Parameter input file: '%s'" % ParFileName)
    with open(ParFileName,'r') as f:
        g=f.readline().split()
        pPar=int(g[0])
        Type=g[1]
        Parameter=f.readline().strip()
        #if len(Parameter)>=2:
        #  Parameter=Parameter[1:-1]
        if not Parameter:
            Parameter="0"
        Boundary=(True,set(map(int,f.read().split())))
    return (Type,Parameter,Boundary)

def _GetPls(ParFileName):
    """
    Lese Randbeschreibungen aus einer Parameterdatei. Maximale Knotenzahl NVT
    bestimmt zudem die Länge der Randliste.
    Rückgabe: (Name des Randes, Daten des Randes, Boolsche Liste für alle Knoten)
    """
    print("Parameter input file: '%s'" % ParFileName)
    with open(ParFileName,'r') as f:
        g=f.readline().split()
        pPar=int(g[0])
        Type=g[1]
        Parameter=f.readline().strip()
        #if len(Parameter)>=2:
        #  Parameter=Parameter[1:-1]
        if not Parameter:
            Parameter="0"
        SFaces=set()
        for i in range(pPar):
            s=f.readline()
            Face=tuple(map(int,s.replace(","," ").split()))
            SFaces.add(Face)
        Boundary=(False,SFaces)
    return (Type,Parameter,Boundary)

def GetParametrisation(ParFileName,Extension):
    """
    Lese eine Parametrisierung aus einer Datei.
    Format wird über "Extension" bestimmt.
    """
    switch = {
        "par": _GetPar,
        "pls": _GetPls
    }
    return (switch[Extension](ParFileName))

def _build_line_by_format_list(format,L,sep=" "):
  return sep.join(map(lambda x: format % (x,),L))+"\n"

def _OutputParFile(Name,Type,Parameters,Boundary):
  assert Boundary[0], "Fehler: PAR kann nur Randbedingungen auf Knoten abspeichern!"
  Boundary=Boundary[1]
  print("Output parameter file: " + Name)
  with open(Name,"w") as f:
    f.write("%d %s\n"%(len(Boundary),Type))
    f.write(Parameters+"\n")
    f.write(_build_line_by_format_list("%d",Boundary,"\n"))

def _OutputPlsFile(Name,Type,Parameters,Boundary):
  assert not Boundary[0], "Fehler: PLS kann nur Randbedingungen auf Seitenflaechen abspeichern!"
  Boundary=Boundary[1]
  print("Output parameter file: " + Name)
  with open(Name,"w") as f:
    f.write("%d %s\n"%(len(Boundary),Type))
    f.write(Parameters+"\n")
    for Face in Boundary:
      f.write(_build_line_by_format_list("%d",Face,","))
      
def OutputParametrisation(Name,Type,Parameters,Boundary,Extension):
    """
    Schreibe eine Parametrisierung aus einer Datei.
    Format wird über "Extension" bestimmt.
    """
    switch = {
        "par": _OutputParFile,
        "pls": _OutputPlsFile
    }
    switch[Extension](Name,Type,Parameters,Boundary)

def OutputGrid(Name,Grid):
  # Auspacken der Gitterstruktur in einzelne Variablen
  (nel,nvt,coord,kvert,knpr)=Grid
  print("Output grid file: " + Name)
  with open(Name,'w') as f:
    f.write("Coarse mesh exported by Partitioner\n")
    f.write("Parametrisierung PARXC, PARYC, TMAXC\n")
    f.write("%d %d" % (nel,nvt))
    f.write(" 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE\nDCORVG\n")
    for node in coord:
      f.write(_build_line_by_format_list("%.17f",node))
    f.write("KVERT\n")
    for elem in kvert:
      f.write(_build_line_by_format_list("%d",elem))
    f.write("KNPR\n")
    # Auf einer Zeile
    #f.write(_build_line_by_format_list("%d",knpr))
    # Jeder Eintrag auf einer eigenen Zeile
    f.write(_build_line_by_format_list("%d",knpr,"\n"))

def GetNeigh(Grid):
  """
  Bestimme für jedes Element eines Gitters eine Liste mit Nachbarelementen.
  """
  face=((0,1,2,3),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(4,5,6,7))
  (NEL,NVT,KVert)=Grid[:2]+Grid[3:4]
  # Erzeuge für jeden Knoten eine Liste mit Elementen,
  # welche diesen Knoten enthalten.
  AuxStruct=[set() for i in range(NVT)]
  for (Elem_Num,Elem) in enumerate(KVert,1):
    for Vert in Elem:
      AuxStruct[Vert-1].add(Elem_Num)
  # Baue Liste mit Nachbarschaften auf
  Neigh=[]
  for (Elem_Num,Elem) in enumerate(KVert,1):
    n=[0,]*6
    for j in range(6):
      # Die folgenden vier Zeilen ersetzen "NeighFinder"
      s=reduce(set.intersection,[AuxStruct[Elem[i]-1] for i in face[j]])
      s.discard(Elem_Num)
      if s:
        n[j]=s.pop()
    Neigh.append(tuple(n))     
  return tuple(Neigh)

def _print_c_array(A):
  print("("+", ".join(map(str,A))+")")

def GetAtomicSplitting(Num):
  return tuple(range(1,Num+1))

def GetParts(Neigh,nPart,Method):
  # Falls nPart==1 ist, erzeuge direkt eine Dummy-Partitionierung
  if nPart==1:
    return (1,)*len(Neigh)
  # Falls die Anzahl der gesuchten Unterteilungen größer oder gleich der Anzahl der Zellen ist,
  # führen eine atomare Aufteilung des Gitters in einzelne Zellen durch.
  # Dies behebt ein merkwürdiges Verhalten von Metis, dass in diesem Fall Unterteilungen
  # mit 0 Elementen und mehr als einem Element erzeugt.
  if len(Neigh)<=nPart:
    return GetAtomicSplitting(len(Neigh))
  # Ein paar Einstellungsparameter
  cOpts=(c_int * 5)(0,100,4,1,1)
  cNum=c_int(1) # Nummerierung beginnt mit 1
  cWeight=c_int(0) # Keine Gewichte

  # Zähle alle nichtnull Elemente der Nachbarschaftsliste
  iCount=sum(list(map(lambda x: list(map(lambda y: bool(y),x)).count(True),Neigh)))

  # Alloziere die Listen MetisA, MetisB und Part
  NEL=len(Neigh)
  MetisA=(c_int * (NEL+1))()
  MetisB=(c_int * iCount)()
  Part=(c_int * NEL)()
  # Baue die komprimierte Graphenstruktur auf
  iOffset=1
  for (Idx,Elem_Neigh) in enumerate(Neigh):
    MetisA[Idx]=iOffset
    for iNeigh in Elem_Neigh:
      if iNeigh:
        MetisB[iOffset-1]=iNeigh
        iOffset+=1
  MetisA[NEL]=iCount+1
  # Rufe Metis auf
  null_ptr = POINTER(c_int)()
  cNEL=c_int(NEL)
  cnPart=c_int(nPart)
  EdgeCut=c_int()
  print("Calling Metis...")
  metis_func[Method-1](byref(cNEL),MetisA,MetisB,null_ptr,null_ptr,\
                       byref(cWeight),byref(cNum),byref(cnPart),cOpts,\
                       byref(EdgeCut),Part)
  print("%d edges were cut by Metis." % EdgeCut.value)
  # Fertig
  return tuple(Part)

def Flatten3dArray(maxZ, maxY, i, j, k):  
  idx1D = (i - 1) * maxZ * maxY + (j - 1) * maxZ + k - 1 
  return idx1D

def GetSubs(BaseName,Grid,nPart,Part,Neigh,nParFiles,Param,bSub, nSubMesh): # to be changed!
  face=((0,1,2,3),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(4,5,6,7))
  
  if isinstance(nSubMesh, int):
    subMeshes = nSubMesh
    GetSubsClassic(BaseName,Grid,nPart,Part,Neigh,nParFiles,Param,bSub)
    return
  else:
    subMeshes = nSubMesh[0] * nSubMesh[1] * nSubMesh[2]
    partX = nSubMesh[0]
    partY = nSubMesh[1]
    partZ = nSubMesh[2]

  # Auspacken der Gitterstruktur in einzelne Variablen
  (nel,nvt,coord,kvert,knpr)=Grid
  # Auspacken der Parametrisierungen
  (ParNames,ParExts,ParTypes,Parameters,Boundaries,NodeFlags)=Param
  # Add new boundary nodes at partition borders
  new_knpr=list(knpr)
  for (iPart,iNeigh,iElem) in zip(Part,Neigh,kvert):
    for (Idx,f) in zip(iNeigh,face):
      if Idx>0 and Part[Idx-1]!=iPart:
        for k in range(4):
          new_knpr[iElem[f[k]]-1]=1

  if not isinstance(nSubMesh, int):
    print("Partitioning scheme: {}x, {}y, {}z\n".format(nSubMesh[0], nSubMesh[1], nSubMesh[2]))
  # Für alle Rechengebiete
  # loop from [0, 0, 0] to [n, n, n]
  for iPartX in range(1,partX+1):
    for iPartY in range(1,partY+1):
      for iPartZ in range(1,partZ+1):
        iPart = [iPartZ, iPartY, iPartX]
        # Bestimme, welche Zellen und Knoten in diesem Gebiet liegen 
        iElem=tuple(eNum for (eNum,p) in enumerate(Part) if p==iPart)
        iCoor=set(vert-1 for eNum in iElem for vert in kvert[eNum])
        # Erzeuge Lookup-Listen: Neue-Idx->Alte Idx
        iCoor=list(iCoor)
        iCoor.sort()
        iCoor=tuple(iCoor)
        # Mappe Knotenkoordinaten und Knoteneigenschaften
        dCoor=tuple(coord[Idx] for Idx in iCoor)
        dKnpr=tuple(new_knpr[Idx] for Idx in iCoor)
        # Erzeuge Lookup-Liste: Alte Knotennummern->Neue Knotennummern
        LookUp=dict((k+1,v) for (v,k) in enumerate(iCoor,1))
        # Mappe die Knoten der Elemente
        dKvert=tuple(tuple(map(lambda x:LookUp[x],kvert[Idx])) for Idx in iElem)
        # Gitterausgabe
        localGrid=(len(dKvert),len(dCoor),dCoor,dKvert,dKnpr)
        if bSub:
          idx1D2 = Flatten3dArray(partZ, partY, iPart[2], iPart[1], iPart[0])  
          idx1D2 = idx1D2 + 1
          localGridName=os.path.join(BaseName,"GRID%04d.tri"%idx1D2)
        else:
          # 3D->1D map
          # [ix * (yMax * zMax)] + (iy * zMax)  + iz
          idx1D2 = Flatten3dArray(partZ, partY, iPart[2], iPart[1], iPart[0])  
          idx1D2 = idx1D2 + 1
          localGridName=os.path.join(BaseName,"sub%04d"%idx1D2,"GRID.tri")
        OutputGrid(localGridName,localGrid)

        if not isinstance(nSubMesh, int):
          id = 1
          localGridName=os.path.join(BaseName,"sub%04d" %idx1D2, "GRID%04d.tri"%id)
          OutputGrid(localGridName,localGrid)

        ###

        localRestriktion=set(LookUp.keys())
        for iPar in range(nParFiles):
          if bSub:
            idx1D2 = Flatten3dArray(partZ, partY, iPart[2], iPart[1], iPart[0])  
            idx1D2 = idx1D2 + 1
            localParName=os.path.join(BaseName,"%s_%04d.%s"%(ParNames[iPar],idx1D2,ParExts[iPar]))
          else:
            idx1D2 = Flatten3dArray(partZ, partY, iPart[2], iPart[1], iPart[0])  
            idx1D2 = idx1D2 + 1
            localParName=os.path.join(BaseName,"sub%04d"%idx1D2,"%s.%s"%(ParNames[iPar],ParExts[iPar]))
          # Alle Knoten, welche innerhalb des neuen Teilgebietes liegen, werden mit der lokalen Knotennummer versehen.  
          # Wenn alle Knoten einer Seite Teil der alten Randparametrisierung ist und ganz neuen Teilgebiet
          # liegt dann gehoert diese dort auch zur Randparametrisierung
          (NodeFlag,Bd)=Boundaries[iPar]
          if NodeFlag:
            localBoundary=(True,[LookUp[i] for i in (Bd&localRestriktion)])
          else:
            localBoundary=(False,[tuple(LookUp[i] for i in Face) for Face in Bd if localRestriktion.issuperset(Face)])
          localBoundary[1].sort()
          OutputParametrisation(localParName,ParTypes[iPar],Parameters[iPar],localBoundary,ParExts[iPar])

          if not isinstance(nSubMesh, int):
            id = 1
            localParName=os.path.join(BaseName,"sub%04d" %idx1D2,"%s_%04d.%s"%(ParNames[iPar],id,ParExts[iPar]))
            OutputParametrisation(localParName,ParTypes[iPar],Parameters[iPar],localBoundary,ParExts[iPar])

def GetSubsClassic(BaseName,Grid,nPart,Part,Neigh,nParFiles,Param,bSub): # to be changed!
  face=((0,1,2,3),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(4,5,6,7))
  # Auspacken der Gitterstruktur in einzelne Variablen
  (nel,nvt,coord,kvert,knpr)=Grid
  # Auspacken der Parametrisierungen
  (ParNames,ParExts,ParTypes,Parameters,Boundaries)=Param
  # Add new boundary nodes at partition borders
  new_knpr=list(knpr)

  for (iPart,iNeigh,iElem) in zip(Part,Neigh,kvert):
    for (Idx,f) in zip(iNeigh,face):
      if Idx>0 and Part[Idx-1]!=iPart:
        for k in range(4):
          new_knpr[iElem[f[k]]-1]=1
  # Für alle Rechengebiete
  for iPart in range(1,nPart+1):
    # Bestimme, welche Zellen und Knoten in diesem Gebiet liegen 
    iElem=tuple(eNum for (eNum,p) in enumerate(Part) if p==iPart)
#    print(len(iElem))
    iCoor=set(vert-1 for eNum in iElem for vert in kvert[eNum])
    # Erzeuge Lookup-Listen: Neue-Idx->Alte Idx
    iCoor=list(iCoor)
    iCoor.sort()
    iCoor=tuple(iCoor)
    # Mappe Knotenkoordinaten und Knoteneigenschaften
    dCoor=tuple(coord[Idx] for Idx in iCoor)
    dKnpr=tuple(new_knpr[Idx] for Idx in iCoor)
    # Erzeuge Lookup-Liste: Alte Knotennummern->Neue Knotennummern
    LookUp=dict((k+1,v) for (v,k) in enumerate(iCoor,1))
    # Mappe die Knoten der Elemente
    dKvert=tuple(tuple(map(lambda x:LookUp[x],kvert[Idx])) for Idx in iElem)
    # Gitterausgabe
    localGrid=(len(dKvert),len(dCoor),dCoor,dKvert,dKnpr)
    if bSub:
      localGridName=os.path.join(BaseName,"GRID%04d.tri"%iPart)
    else:
      localGridName=os.path.join(BaseName,"sub%04d"%iPart,"GRID.tri")
    OutputGrid(localGridName,localGrid)

    ###

    localRestriktion=set(LookUp.keys())
    for iPar in range(nParFiles):
      if bSub:
        localParName=os.path.join(BaseName,"%s_%04d.%s"%(ParNames[iPar],iPart,ParExts[iPar]))
      else:
        localParName=os.path.join(BaseName,"sub%04d"%iPart,"%s.%s"%(ParNames[iPar],ParExts[iPar]))
      # Alle Knoten, welche innerhalb des neuen Teilgebietes liegen, werden mit der lokalen Knotennummer versehen.  
      # Wenn alle Knoten einer Seite Teil der alten Randparametrisierung ist und ganz neuen Teilgebiet
      # liegt dann gehoert diese dort auch zur Randparametrisierung
      (NodeFlag,Bd)=Boundaries[iPar]
      if NodeFlag:
        localBoundary=(True,[LookUp[i] for i in (Bd&localRestriktion)])
      else:
        localBoundary=(False,[tuple(LookUp[i] for i in Face) for Face in Bd if localRestriktion.issuperset(Face)])
      localBoundary[1].sort()
      OutputParametrisation(localParName,ParTypes[iPar],Parameters[iPar],localBoundary,ParExts[iPar])

def MultPartitionAlongAxis(Grid,nSubMesh,Method):
  (nel,nvt,coord,kvert,knpr)=Grid
  # An array that tells you in which partition the i-th is
  # Let Part be a list of tuples (x, y, z) where x, y, z are the
  # cartesian indices of the partition
  Part=[0,]*nel
  Dir = 2
  zCoords = [p[2] for p in coord]
  numCoords = len(zCoords)  
  zCoords.sort()
  zMin = zCoords[0]
  zMax = zCoords[numCoords-1]

  # The delta for the current subdivion
  dZ = (zMax - zMin) / nSubMesh
  theList = [i * dZ for i in range(1, nSubMesh + 1)]
  print(zMin)
  print(zMax)
  print(dZ)
  print(theList)
  PosFak=1
  for (ElemIdx,Elem) in enumerate(kvert):
    for idx, val in enumerate(theList):
      if all([(coord[Vert-1][Dir] -val <= 1e-5) for Vert in Elem]):
        Part[ElemIdx]=idx + 1
        break

  return tuple(Part)

def AxisBasedPartitioning(Grid,nSubMesh,Method):
  (nel,nvt,coord,kvert,knpr)=Grid
  # An array that tells you in which partition the i-th is
  # Let Part be a list of tuples (x, y, z) where x, y, z are the
  # cartesian indices of the partition
  #Part=[[0, 0, 0],]*nel
  Part = []
  for i in range(nel):
    Part.append([0,0,0])

#  Part = [ [0,0,0] for i in nel ]
  Dir = 2
  order = 0
  zCoords = [p[Dir] for p in coord]
  numCoords = len(zCoords)  
  zCoords.sort()
  zMin = zCoords[0]
  zMax = zCoords[numCoords-1]

  # The delta for the z-subdivision
  dZ = (zMax - zMin) / nSubMesh[Dir]
  theList = [zMin + i * dZ for i in range(1, nSubMesh[Dir] + 1)]
  print(zMin)
  print(zMax)
  print(dZ)
  print(theList)
  PosFak=1
  valCount = [0] * len(theList)
  for (ElemIdx,Elem) in enumerate(kvert):
    for idx, val in enumerate(theList):
        count = 0
        for Vert in Elem:
          dist = coord[Vert-1][Dir] - (val + zMin)  
          if dist <= 1e-5:
            count = count + 1
            valCount[idx] = valCount[idx] + 1
        if count == 8:
          Part[ElemIdx][order]=idx + 1
          break
  print(valCount)
  
  # y-subdivision
  Dir = 1
  order = 1
  yCoords = [p[Dir] for p in coord]
  numCoords = len(yCoords)  
  yCoords.sort()
  yMin = yCoords[0]
  yMax = yCoords[numCoords-1]

  # The delta for the y-subdivision
  dY = (yMax - yMin) / nSubMesh[Dir]
  theList = [i * dY for i in range(1, nSubMesh[Dir] + 1)]
  print(yMin)
  print(yMax)
  print(dY)
  print(theList)
  PosFak=1
  for (ElemIdx,Elem) in enumerate(kvert):
    for idx, val in enumerate(theList):
      count = 0
      for Vert in Elem:
        dist = coord[Vert-1][Dir] - (yMin + val)
        if dist <= 1e-5:
          count = count + 1
      if count == 8:
        Part[ElemIdx][order]=idx + 1
        break
  # x-subdivision
  Dir = 0
  order = 2
  xCoords = [p[Dir] for p in coord]
  numCoords = len(xCoords)  
  xCoords.sort()
  xMin = xCoords[0]
  xMax = xCoords[numCoords-1]

  # The delta for the x-subdivision
  dX = (xMax - xMin) / nSubMesh[Dir]
  theList = [i * dX for i in range(1, nSubMesh[Dir] + 1)]
  print(xMin)
  print(xMax)
  print(dX)
  print(theList)
  PosFak=1
  for (ElemIdx,Elem) in enumerate(kvert):
    for idx, val in enumerate(theList):
      count = 0
      for Vert in Elem:
        dist = coord[Vert-1][Dir] - (xMin + val)
        if dist <= 1e-5:
          count = count + 1
      if count == 8:
        Part[ElemIdx][order]=idx + 1
        break

  return tuple(Part)

def PartitionAlongAxis(Grid,nSubMesh,Method):
  # Berechne 1D Median einer Liste (die Liste wird dabei sortiert)
  def median(L):
    Length=len(L)
    assert Length>0, "Only for non-empty lists can a median be computed!"
    L.sort()
    Idx=(Length-1)//2
    return (L[Idx]+L[Idx+1])/2.0 if Length%2==0 else L[Idx]
  # Ab hier fängt die eigentliche Routine an.
  # Bestimme zuerst, ob die Parameter gültig sind
  assert Method<0, "Only Methods <0 are valid!"
  tmp=str(-Method)
  assert tmp.strip("1234")=="", "Only 1, 2, 3 or 4 are valid axis!"
  Axis=list(map(lambda char: char in tmp,"1234"))

  if -Method == 4:
    return AxisBasedPartitioning(Grid,nSubMesh,Method)

  NumAxis=sum(Axis)
  nSub=2**NumAxis

  assert nSub==nSubMesh, "Your subgrid splitting choice requires exactly %d subgrids!"%nSub  
  # Entpacke die Informationen in Parameter Grid
  (nel,nvt,coord,kvert,knpr)=Grid
  # Initialisiere Gebietsaufteilung der Elemente (Am Anfang nur ein Gebiet!)
  Part=[1,]*nel
  # Spaltungsprozedere für alle gewählten Richtungen
  PosFak=1
  for Dir in range(3):
    if Axis[Dir]:
      # Bestimme Median für die gewählte Richtung
      Mid=median([p[Dir] for p in coord])
      # Teile die Elemente dahingehend auf, ob alle Knoten <=Mid sind oder nicht
      for (ElemIdx,Elem) in enumerate(kvert):
        if all([(coord[Vert-1][Dir]>=Mid) for Vert in Elem]):
          Part[ElemIdx]+=PosFak
      # Bestimme nächste 2er Potenz im Stellenwertsystem
      PosFak*=2
  return tuple(Part)

# Dokumentation
__doc__ = \
"""
Dieses Modul führt die Partitionierung eines Gitters mittels der Metis-Bibliothek durch.
"""
# Startroutine des Moduls, die Metis lädt.
if os.name=="posix":
  metis=_try_in_place_first("libmetis.so")
elif os.name=="nt":
  metis=_try_in_place_first("metis.dll")
else:
  sys.exit("Loading of Metis not yet implemented for platform '%s'!"%os.name)

if metis==None:
  sys.exit("Could not load the Metis library!")

# Füge Aufrufparameter von den drei verwendeten Metis-Funktionen hinzu
_pidx=POINTER(c_int)
_pint=POINTER(c_int)
_PartArgs=(_pint,_pidx,_pidx,_pidx,_pidx,_pint,_pint,_pint,_pint,_pint,_pidx)
metis.METIS_PartGraphRecursive.argtypes=_PartArgs
metis.METIS_PartGraphVKway.argtypes=_PartArgs
metis.METIS_PartGraphKway.argtypes=_PartArgs
metis_func=(metis.METIS_PartGraphRecursive,metis.METIS_PartGraphVKway,metis.METIS_PartGraphKway)

if __name__=="__main__":
  if metis!=None:
    print("Metis has been loaded.")

