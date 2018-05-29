#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
from part_main import *


# Erzeuge benötigte Verzeichnisse, falls noch nicht vorhanden
mkdir("_mesh")

NPart, PartMethod, NSubPart, MeshName, ProjektFile = checkParameters(sys.argv) 

# Aufruf der Hauptroutine
MainProcess(NPart,PartMethod,NSubPart,MeshName,ProjektFile)
