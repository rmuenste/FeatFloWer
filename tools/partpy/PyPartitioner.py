#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Tested with:
#         NPart  Meth Subs                   Prj-Path
# "args": ["15", "1", "1", "NEWFAC", "./2D_FAC/2Dbench.prj"]
# "args": ["27", "-4", "x3-y3-z3", "NEWFAC", "./dev3x3x3/dev3x3x3.prj"]
# "args": ["2", "-3", "2", "NEWFAC2", "./2D_FAC/2Dbench.prj"]

import sys
import partitioner


# Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
partitioner.mkdir("_mesh")

NPart, PartMethod, NSubPart, MeshName, ProjektFile = partitioner.checkParameters(sys.argv) 

# Aufruf der Hauptroutine
partitioner.MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
