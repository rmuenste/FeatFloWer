#!/usr/bin/env python
# vim: set filetype=python

import os
import sys
import re

# Scan what is in the _dump/processor_X folder 
# for each .dmp field found there
# read and join the partitioned files to
# a single file per field 

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   

def cmp_entries(a, b):
# sort the tuples in the list element_entries 
# according to the global element index
  if a[0] > b[0]:
    return 1
  elif a[0] == b[0]:
    return 0
  else:
    return -1

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   

def readPartitionedField(elem_entries, fileName):
# reads a single field from the file <fileName>
# and stores it in elem_entries         

  header_info = {}
  header_line = ""

  with open(fileName, "r") as f:
    header = f.readline()
    header_line = header
    print(header)
    split_header = header.strip().split(",")
    print(split_header)
    header_info = {(pair[0]) : (pair[1]) for pair in (l.split(":") for l in split_header)}
    print(header_info['Components'])

    while True:
      line = f.readline()

      if not line: break

      element_idx = int(line)
      element_components = []
      entry = (element_idx, element_components)

      for i in range(int(header_info['Components'])):
        line = f.readline()
        if not line: break
        line = line.strip()
        values = line.split(" ")
        entry[1].append(values)

      elem_entries.append(entry)

    return header_line

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   

def writeMergedField(elem_entries, header, components, fileName):
  with open(fileName, "w") as f:
    f.write(header)
    for e in elem_entries:
      f.write(str(e[0]) + '\n')
      for i in range(components):
        f.write(" ".join(e[1][i]) + '\n')

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   



#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-   
        
# Parameters for reading and merging the field:
# Number of processors: 
# - equal to number of processor_* directories in the folder
# - the dump folder index 

dir_list = os.listdir("./_dump")

dump_folder_idx = 1

step = dump_folder_idx

num_processors = 0

for e in dir_list:
  m = re.search(r"processor_\d+", e)
  if m:
    num_processors = num_processors + 1

print("Found " + str(num_processors) + " processor_ directories")

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

header_line = ""
element_entries = []

#for i in range(1,5):
#  fileToRead = "./_dump/processor_" + str(i) + "/" + str(step) + "/" + "velocity.dmp"
#  header_line = readPartitionedField(element_entries, fileToRead)
#
#    #print("Using directories: " + "./_dump/processor_" + str(i) + "/" + str(step))
#
#element_entries.sort(cmp_entries)
#
#print("Number of elements entries: " + str(len(element_entries)))
#    
#writeMergedField(element_entries, header_line, header_line['Components'], "velocity.dmp"):

