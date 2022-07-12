#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Tested with:
# "args": ["8", "-3", "3", "NEWFAC", "./dev_fine/dev_fine.prj"]

import sys
import partitioner


# Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
partitioner.mkdir("_mesh")

NPart, PartMethod, NSubPart, MeshName, ProjektFile = partitioner.checkParameters(sys.argv) 

# Aufruf der Hauptroutine
partitioner.MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
