#!/usr/bin/env python
# -*- coding: utf-8 -*-

# METIS library loaded lazily via prioritised search (Option C)
from ._libmetis import get_metis_func

from ctypes import c_int, c_float, POINTER, byref
from functools import reduce
from itertools import repeat, count
from collections import Counter
from math import sqrt
import os
import sys
from shutil import copy


# Private helpers

def _readAfterKeyword(fh, keyword):
    s = fh.readline()
    while keyword not in s:
        s = fh.readline()


# Public module functions

def GetFileList(cProjName):
    """
    Read project file; return grid filename and parameter file names.
    """
    nPar = 0
    myGridFile = ""
    myParFiles = []
    myParNames = []
    ProjektDir = os.path.dirname(cProjName)
    fProjekt = open(cProjName, 'r')
    for s in fProjekt:
        s = s.strip()
        if '.' in s:
            (prefix, sep, ext) = s.rpartition('.')
            if ext == "tri":
                myGridFile = os.path.join(ProjektDir, s)
            elif ext == "par":
                nPar += 1
                myParFiles.append(os.path.join(ProjektDir, s))
                myParNames.append(prefix)
    fProjekt.close()
    print("The projekt folder consists of the following files:")
    print("- Grid File:", myGridFile)
    print("- Boundary Files:")
    print("\n".join(map(lambda x: "  * %s" % x, myParFiles)))
    if not myGridFile:
        raise RuntimeError(
            "No .tri grid file entry found in project file '%s'.\n"
            "Check that the file is not empty and contains a line ending in '.tri'." % cProjName
        )
    if not os.path.exists(myGridFile):
        raise RuntimeError(
            "Grid file '%s' listed in '%s' does not exist." % (myGridFile, cProjName)
        )
    missing_par = [f for f in myParFiles if not os.path.exists(f)]
    if missing_par:
        raise RuntimeError(
            "The following boundary file(s) listed in '%s' do not exist:\n  %s"
            % (cProjName, "\n  ".join(missing_par))
        )
    return (nPar, myGridFile, myParFiles, myParNames)


def GetGrid(GridFileName):
    """
    Read a grid from file.  Returns (NEL, NVT, Coord, KVert, Knpr).
    """
    print("Grid input file: '%s'" % GridFileName)
    f = open(GridFileName, 'r')
    f.readline()
    f.readline()
    g = f.readline().split()
    NEL = int(g[0])
    NVT = int(g[1])
    _readAfterKeyword(f, "DCORVG")
    Coord = []
    for i in range(NVT):
        g = f.readline().split()
        Coord.append(tuple(map(float, g)))
    _readAfterKeyword(f, "KVERT")
    Kvert = []
    for i in range(NEL):
        g = f.readline().split()
        Kvert.append(tuple(map(int, g)))
    _readAfterKeyword(f, "KNPR")
    g = f.read().split()
    Knpr = tuple(map(int, g))
    return (NEL, NVT, tuple(Coord), tuple(Kvert), Knpr)


def GetPar(ParFileName, NVT):
    """
    Read boundary descriptions from a parameter file.
    """
    print("Parameter input file: '%s'" % ParFileName)
    with open(ParFileName, 'r') as f:
        g = f.readline().split()
        pPar = int(g[0])
        Type = g[1]
        Parameter = f.readline().strip()
        if not Parameter:
            Parameter = "0"
        Boundary = set(map(int, f.read().split()))
    return (Type, Parameter, Boundary)


def GetNeigh(Grid):
    """
    Build neighbour list for every element.
    """
    face = ((0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5),
            (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7))
    (NEL, NVT, KVert) = Grid[:2] + Grid[3:4]
    AuxStruct = [set() for i in range(NVT)]
    for (Elem_Num, Elem) in enumerate(KVert, 1):
        for Vert in Elem:
            AuxStruct[Vert - 1].add(Elem_Num)
    Neigh = []
    for (Elem_Num, Elem) in enumerate(KVert, 1):
        n = [0, ] * 6
        for j in range(6):
            s = reduce(set.intersection, [AuxStruct[Elem[i] - 1] for i in face[j]])
            s.discard(Elem_Num)
            if s:
                n[j] = s.pop()
        Neigh.append(tuple(n))
    return tuple(Neigh)


def _print_c_array(A):
    print("(" + ", ".join(map(str, A)) + ")")


def GetAtomicSplitting(Num):
    return tuple(range(1, Num + 1))


def GetParts(Neigh, nPart, Method):
    if nPart == 1:
        return (1,) * len(Neigh)
    if len(Neigh) <= nPart:
        return GetAtomicSplitting(len(Neigh))

    iCount = sum(sum(1 for y in x if y) for x in Neigh)

    NEL = len(Neigh)
    MetisA = (c_int * (NEL + 1))()
    MetisB = (c_int * iCount)()
    Part = (c_int * NEL)()

    iOffset = 0
    for (Idx, Elem_Neigh) in enumerate(Neigh):
        MetisA[Idx] = iOffset
        for iNeigh in Elem_Neigh:
            if iNeigh:
                MetisB[iOffset] = iNeigh - 1
                iOffset += 1
    MetisA[NEL] = iCount

    cNEL    = c_int(NEL)
    cNcon   = c_int(1)
    cnPart  = c_int(nPart)
    EdgeCut = c_int()

    func_idx = 0 if Method == 1 else 1
    metis_func = get_metis_func()
    print("Calling Metis...")
    metis_func[func_idx](
        byref(cNEL), byref(cNcon),
        MetisA, MetisB,
        None, None, None,
        byref(cnPart),
        None, None,
        None,
        byref(EdgeCut), Part)
    print("%d edges were cut by Metis." % EdgeCut.value)

    return tuple(p + 1 for p in Part)


def GetSubs(BaseName, Grid, nPart, Part, Neigh, nParFiles, Param, bSub):
    face = ((0, 1, 2, 3), (0, 1, 5, 4), (1, 2, 6, 5),
            (2, 3, 7, 6), (3, 0, 4, 7), (4, 5, 6, 7))
    (nel, nvt, coord, kvert, knpr) = Grid
    (ParNames, ParTypes, Parameters, Boundaries) = Param
    new_knpr = list(knpr)

    for (iPart, iNeigh, iElem) in zip(Part, Neigh, kvert):
        for (Idx, f) in zip(iNeigh, face):
            if Idx > 0 and Part[Idx - 1] != iPart:
                for k in range(4):
                    new_knpr[iElem[f[k]] - 1] = 1

    for iPart in range(1, nPart + 1):
        iElem = tuple(eNum for (eNum, p) in enumerate(Part) if p == iPart)
        print(len(iElem))
        iCoor = set(vert - 1 for eNum in iElem for vert in kvert[eNum])
        iCoor = list(iCoor)
        iCoor.sort()
        iCoor = tuple(iCoor)
        dCoor = tuple(coord[Idx] for Idx in iCoor)
        dKnpr = tuple(new_knpr[Idx] for Idx in iCoor)
        LookUp = dict((k + 1, v) for (v, k) in enumerate(iCoor, 1))
        dKvert = tuple(tuple(map(lambda x: LookUp[x], kvert[Idx])) for Idx in iElem)
        localGrid = (len(dKvert), len(dCoor), dCoor, dKvert, dKnpr)
        if bSub:
            localGridName = os.path.join(BaseName, "GRID%04d.tri" % iPart)
        else:
            localGridName = os.path.join(BaseName, "sub%04d" % iPart, "GRID.tri")
        OutputGrid(localGridName, localGrid)

        localRestriktion = set(LookUp.keys())
        for iPar in range(nParFiles):
            if bSub:
                localParName = os.path.join(BaseName, "%s_%04d.par" % (ParNames[iPar], iPart))
            else:
                localParName = os.path.join(BaseName, "sub%04d" % iPart, "%s.par" % ParNames[iPar])
            localBoundary = [LookUp[i] for i in (Boundaries[iPar] & localRestriktion)]
            localBoundary.sort()
            OutputParFile(localParName, ParTypes[iPar], Parameters[iPar], localBoundary)


def _build_line_by_format_list(format, L, sep=" "):
    return sep.join(map(lambda x: format % (x,), L)) + "\n"


def OutputParFile(Name, Type, Parameters, Boundary):
    with open(Name, "w") as f:
        f.write("%d %s\n" % (len(Boundary), Type))
        f.write(Parameters + "\n")
        f.write(_build_line_by_format_list("%d", Boundary, "\n"))


def OutputGrid(Name, Grid):
    (nel, nvt, coord, kvert, knpr) = Grid
    with open(Name, 'w') as f:
        f.write("Coarse mesh exported by Partitioner\n")
        f.write("Parametrisierung PARXC, PARYC, TMAXC\n")
        f.write("%d %d" % (nel, nvt))
        f.write(" 1 8 12 6     NEL,NVT,NBCT,NVE,NEE,NAE\nDCORVG\n")
        for node in coord:
            f.write(_build_line_by_format_list("%.17f", node))
        f.write("KVERT\n")
        for elem in kvert:
            f.write(_build_line_by_format_list("%d", elem))
        f.write("KNPR\n")
        f.write(_build_line_by_format_list("%d", knpr, "\n"))


def MultPartitionAlongAxis(Grid, nSubMesh, Method):
    (nel, nvt, coord, kvert, knpr) = Grid
    Part = [0, ] * nel
    Dir = 2
    zCoords = [p[2] for p in coord]
    numCoords = len(zCoords)
    zCoords.sort()
    zMin = zCoords[0]
    zMax = zCoords[numCoords - 1]
    dZ = (zMax - zMin) / nSubMesh
    theList = [i * dZ for i in range(1, nSubMesh + 1)]
    print(zMin)
    print(zMax)
    print(dZ)
    print(theList)
    for (ElemIdx, Elem) in enumerate(kvert):
        for idx, val in enumerate(theList):
            if all([(coord[Vert - 1][Dir] - val <= 1e-5) for Vert in Elem]):
                Part[ElemIdx] = idx + 1
                break
    return tuple(Part)


def PartitionAlongAxis(Grid, nSubMesh, Method):
    def median(L):
        Length = len(L)
        assert Length > 0, "Only for non-empty lists can a median be computed!"
        L.sort()
        Idx = (Length - 1) // 2
        return (L[Idx] + L[Idx + 1]) / 2.0 if Length % 2 == 0 else L[Idx]

    assert Method < 0, "Only Methods <0 are valid!"
    tmp = str(-Method)
    assert tmp.strip("123") == "", "Only 1, 2, or 3 are valid axis!"
    Axis = list(map(lambda char: char in tmp, "123"))
    NumAxis = sum(Axis)
    nSub = 2 ** NumAxis

    if nSub != nSubMesh:
        return MultPartitionAlongAxis(Grid, nSubMesh, Method)

    assert nSub == nSubMesh, "Your subgrid splitting choice requires exactly %d subgrids!" % nSub
    (nel, nvt, coord, kvert, knpr) = Grid
    Part = [1, ] * nel
    PosFak = 1
    for Dir in range(3):
        if Axis[Dir]:
            Mid = median([p[Dir] for p in coord])
            for (ElemIdx, Elem) in enumerate(kvert):
                if all([(coord[Vert - 1][Dir] >= Mid) for Vert in Elem]):
                    Part[ElemIdx] += PosFak
            PosFak *= 2
    return tuple(Part)


__doc__ = """
This module partitions a mesh using the METIS library.
"""
