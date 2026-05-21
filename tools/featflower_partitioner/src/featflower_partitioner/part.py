#!/usr/bin/env python
# -*- coding: utf-8 -*-

# METIS library loaded lazily via prioritised search (Option C)
from ._libmetis import get_metis_func

from ctypes import c_int, c_float, POINTER, byref
from functools import reduce
from itertools import repeat, count
from collections import Counter
import json
from math import sqrt
import os
import sys
from shutil import copy

PARTITION_FORMAT = "legacy"
JSON_WRITER = None
BASE_OUTPUT_PATH = None
GRID_CACHE = {}
PAR_CACHE = {}


# Private helpers

def _readAfterKeyword(fh, keyword):
    s = fh.readline()
    while keyword not in s:
        s = fh.readline()


def cache_grid_data(path, grid):
    GRID_CACHE[os.path.abspath(path)] = grid


def cache_par_data(path, par_data):
    PAR_CACHE[os.path.abspath(path)] = par_data


class JsonPartitionWriter:
    def __init__(self, base_path):
        self.base_path = base_path
        self.files = {}

    def _get_file_entry(self, source_name, file_type):
        entry = self.files.get(source_name)
        if entry is None:
            entry = {
                "type": file_type,
                "master": None,
                "subs": {},
            }
            self.files[source_name] = entry
        return entry

    def _classify_path(self, rel_path):
        normalized = rel_path.replace("\\", "/")
        parts = normalized.split("/")
        if len(parts) > 1 and parts[0].startswith("sub") and parts[0][3:].isdigit():
            return parts[0], "/".join(parts[1:])
        return None, normalized

    def _normalize_part_name(self, base_name, remainder):
        if remainder is None:
            return remainder
        if "/" in remainder or "\\" in remainder:
            return remainder
        base_root, base_ext = os.path.splitext(base_name)
        rem_root, rem_ext = os.path.splitext(remainder)
        if rem_ext.lower() != base_ext.lower():
            return remainder
        suffix_start = len(rem_root)
        while suffix_start > 0 and rem_root[suffix_start - 1].isdigit():
            suffix_start -= 1
        if suffix_start == len(rem_root):
            return remainder
        suffix = rem_root[suffix_start:]
        prefix = rem_root[:suffix_start]
        if prefix.lower()[:len(base_root)] != base_root.lower():
            return remainder
        extra = prefix[len(base_root):]
        if extra not in ("", "_"):
            return remainder
        if not suffix:
            return remainder
        return f"{base_root}.{suffix}{base_ext}"

    def add_grid(self, source_name, rel_path, grid):
        entry = self._get_file_entry(source_name, "tri")
        (nel, nvt, coord, kvert, knpr) = grid
        payload = {
            "nel": nel,
            "nvt": nvt,
            "dcorvg": [list(node) for node in coord],
            "kvert": [list(elem) for elem in kvert],
            "knpr": list(knpr),
        }
        sub_key, remainder = self._classify_path(rel_path)
        if sub_key is None:
            entry["master"] = payload
        else:
            sub_bucket = entry["subs"].setdefault(sub_key, {"coarse": None, "parts": {}})
            normalized = self._normalize_part_name(source_name, remainder)
            if normalized.lower() == source_name.lower():
                sub_bucket["coarse"] = payload
            else:
                sub_bucket["parts"][normalized] = payload

    def add_par(self, source_name, rel_path, par_type, parameter, boundary):
        entry = self._get_file_entry(source_name, "par")
        cleaned = parameter
        if len(cleaned) >= 2 and cleaned[0] == "'" and cleaned[-1] == "'":
            cleaned = cleaned[1:-1]
        payload = {
            "type": par_type,
            "parameter": cleaned,
            "nodes": list(boundary),
        }
        sub_key, remainder = self._classify_path(rel_path)
        base_name = os.path.basename(source_name)
        normalized = self._normalize_part_name(base_name, remainder)
        if sub_key is None:
            entry["master"] = payload
        else:
            sub_bucket = entry["subs"].setdefault(sub_key, {"coarse": None, "parts": {}})
            if normalized.lower() == base_name.lower():
                sub_bucket["coarse"] = payload
            else:
                sub_bucket["parts"][normalized] = payload

    def write_files(self):
        for source_name, payload in self.files.items():
            out_path = os.path.join(self.base_path, "%s.json" % source_name)
            tmp_path = out_path + ".tmp"
            data = {
                "format": "FeatFloWerPartition",
                "version": 1,
                "source": source_name,
                "kind": payload["type"],
                "master": payload["master"],
                "subs": payload["subs"],
            }
            with open(tmp_path, "w") as f:
                json.dump(data, f, indent=2)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp_path, out_path)

        try:
            dir_fd = os.open(self.base_path, os.O_RDONLY)
        except OSError:
            return
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)


def init_partition_output(base_path, partition_format):
    global PARTITION_FORMAT, JSON_WRITER, BASE_OUTPUT_PATH, GRID_CACHE, PAR_CACHE
    PARTITION_FORMAT = partition_format
    BASE_OUTPUT_PATH = base_path
    if PARTITION_FORMAT == "json":
        JSON_WRITER = JsonPartitionWriter(base_path)
    else:
        JSON_WRITER = None
    GRID_CACHE = {}
    PAR_CACHE = {}


def finalize_partition_output():
    if JSON_WRITER is not None:
        JSON_WRITER.write_files()


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
    abs_name = os.path.abspath(GridFileName)
    if PARTITION_FORMAT == "json" and abs_name in GRID_CACHE:
        return GRID_CACHE[abs_name]
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
    abs_name = os.path.abspath(ParFileName)
    if PARTITION_FORMAT == "json" and abs_name in PAR_CACHE:
        return PAR_CACHE[abs_name]
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


def GetSubs(BaseName, Grid, nPart, Part, Neigh, nParFiles, Param, bSub, sourceName="GRID.tri"):
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
        OutputGrid(localGridName, localGrid, sourceName=sourceName)
        if PARTITION_FORMAT == "json":
            cache_grid_data(localGridName, localGrid)

        localRestriktion = set(LookUp.keys())
        for iPar in range(nParFiles):
            if bSub:
                localParName = os.path.join(BaseName, "%s_%04d.par" % (ParNames[iPar], iPart))
            else:
                localParName = os.path.join(BaseName, "sub%04d" % iPart, "%s.par" % ParNames[iPar])
            localBoundary = [LookUp[i] for i in (Boundaries[iPar] & localRestriktion)]
            localBoundary.sort()
            baseParName = "%s.par" % ParNames[iPar]
            OutputParFile(localParName, ParTypes[iPar], Parameters[iPar], localBoundary,
                          source_name=baseParName)
            if PARTITION_FORMAT == "json":
                cache_par_data(localParName, (ParTypes[iPar], Parameters[iPar], set(localBoundary)))


def _build_line_by_format_list(format, L, sep=" "):
    return sep.join(map(lambda x: format % (x,), L)) + "\n"


def OutputParFile(Name, Type, Parameters, Boundary, source_name=None):
    if PARTITION_FORMAT == "json":
        if JSON_WRITER is None or BASE_OUTPUT_PATH is None:
            raise RuntimeError("JSON writer not initialised")
        rel_path = os.path.relpath(Name, BASE_OUTPUT_PATH)
        entry_name = source_name if source_name is not None else os.path.basename(Name)
        JSON_WRITER.add_par(entry_name, rel_path, Type, Parameters, Boundary)
        cache_par_data(Name, (Type, Parameters, set(Boundary)))
        return
    with open(Name, "w") as f:
        f.write("%d %s\n" % (len(Boundary), Type))
        f.write(Parameters + "\n")
        f.write(_build_line_by_format_list("%d", Boundary, "\n"))


def OutputGrid(Name, Grid, sourceName=None):
    (nel, nvt, coord, kvert, knpr) = Grid
    if PARTITION_FORMAT == "json":
        if JSON_WRITER is None or BASE_OUTPUT_PATH is None:
            raise RuntimeError("JSON writer not initialised")
        rel_path = os.path.relpath(Name, BASE_OUTPUT_PATH)
        entry_name = sourceName if sourceName is not None else os.path.basename(Name)
        JSON_WRITER.add_grid(entry_name, rel_path, Grid)
        cache_grid_data(Name, Grid)
        return
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
