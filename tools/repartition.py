#!/usr/bin/env python
# vim: set filetype=python

import os
import sys

dir_list = os.listdir("./_dump")

step = 1

num_processors = len(dir_list)

print("Dump directory: " + str(dir_list))

print("Number of processors: " + str(num_processors))

dump_path = "./_dump/processor_" + str(1) + "/" + str(step)

print("Dump load path: " + dump_path)

print("Backup files: " + str(os.listdir(dump_path)))

# The format of a .dmp file is as follows:
# - it is a list of coarse element entries
# - a coarse element entry consists of multiple lines
# - in the first line a single integer values denotes the global coarse element index
# - after the first line there is a line for each component of the output field
# - in each 'component line' the values of fine mesh dofs in the coarse element are given in the FEATFLOW multi-level ordering

