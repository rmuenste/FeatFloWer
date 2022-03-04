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
  ProjektDir=os.path.dirname(cProjName)
  fProjekt=open(cProjName,'r')
  for s in fProjekt:
    s=s.strip()
    if '.' in s:
      (prefix,sep,ext)=s.rpartition('.')
      if ext=="tri":
        myGridFile=os.path.join(ProjektDir,s)
      elif ext=="par":
        nPar+=1
        myParFiles.append(os.path.join(ProjektDir,s))
        myParNames.append(prefix)
  fProjekt.close()
  print("The projekt folder consists of the following files:")
  print("- Grid File:", myGridFile)
  print("- Boundary Files:")
  print("\n".join(map(lambda x: "  * %s" % x,myParFiles)))
  return (nPar,myGridFile,myParFiles,myParNames)

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

def GetPar(ParFileName,NVT):
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
        Boundary=set(map(int,f.read().split()))
    return (Type,Parameter,Boundary)

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

def GetSubs(BaseName,Grid,nPart,Part,Neigh,nParFiles,Param,bSub):
  face=((0,1,2,3),(0,1,5,4),(1,2,6,5),(2,3,7,6),(3,0,4,7),(4,5,6,7))
  # Auspacken der Gitterstruktur in einzelne Variablen
  (nel,nvt,coord,kvert,knpr)=Grid
  # Auspacken der Parametrisierungen
  (ParNames,ParTypes,Parameters,Boundaries)=Param
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
      localGridName=os.path.join(BaseName,"GRID%03d.tri"%iPart)
    else:
      localGridName=os.path.join(BaseName,"sub%03d"%iPart,"GRID.tri")
    OutputGrid(localGridName,localGrid)

    ###

    localRestriktion=set(LookUp.keys())
    for iPar in range(nParFiles):
      if bSub:
        localParName=os.path.join(BaseName,"%s_%03d.par"%(ParNames[iPar],iPart))
      else:
        localParName=os.path.join(BaseName,"sub%03d"%iPart,"%s.par"%ParNames[iPar])
      localBoundary=[LookUp[i] for i in (Boundaries[iPar]&localRestriktion)]
      localBoundary.sort()
      OutputParFile(localParName,ParTypes[iPar],Parameters[iPar],localBoundary)

def _build_line_by_format_list(format,L,sep=" "):
  return sep.join(map(lambda x: format % (x,),L))+"\n"

def OutputParFile(Name,Type,Parameters,Boundary):
  print("Output parameter file: " + Name)
  with open(Name,"w") as f:
    f.write("%d %s\n"%(len(Boundary),Type))
    f.write(Parameters+"\n")
    f.write(_build_line_by_format_list("%d",Boundary,"\n"))
  pass

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
  assert tmp.strip("123")=="", "Only 1, 2, or 3 are valid axis!"
  Axis=list(map(lambda char: char in tmp,"123"))
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

