#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import partitioner


# Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
partitioner.mkdir("_mesh")

NPart, PartMethod, NSubPart, MeshName, ProjektFile = partitioner.checkParameters(sys.argv) 

# Aufruf der Hauptroutine
partitioner.MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
