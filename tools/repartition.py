#!/usr/bin/env python
# vim: set filetype=python

import os
import sys
import re
import getopt

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
    #print(header)
    split_header = header.strip().split(",")
    #print(split_header)
    header_info = {(pair[0]) : (pair[1]) for pair in (l.split(":") for l in split_header)}
    #print(header_info['Components'])

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

    return header_line, header_info

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def writeCombinedField(elem_entries, header, components, fileName):

  with open(fileName, "w") as f:
    f.write(header)
    for e in elem_entries:
      f.write(str(e[0]) + '\n')
      for i in range(components):
        f.write(" ".join(e[1][i]) + '\n')

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def combineField(nprocs,fieldName, path, out_idx):

  element_entries = []

  for i in range(1,nprocs+1):
    fileToRead = path + "/processor_" + str(i) + "/" + str(out_idx) + "/" + fieldName

    header_line, header_info = readPartitionedField(element_entries, fileToRead)

    element_entries.sort(cmp_entries)

    print("Number of elements entries: " + str(len(element_entries)))

    writeCombinedField(element_entries, header_line, int(header_info['Components']), fieldName)

  print(fieldName + "\n" + header_info['Components'] + "\n" + header_line)

  del header_info['DofsTotal']

  header_info['NEL'] = len(element_entries)
  header_string = str(header_info)

  print(header_string)

  header_bla=[]
  for k, v in header_info.iteritems():
    entry=[]
    entry.append(str(k))
    entry.append(str(v))
    str_entry = ":".join(entry)
    header_bla.append(str_entry)

  print(header_bla) 
  print(",".join(header_bla)) 


#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def getDumpFields(nprocs, path, out_idx):

  print("Found " + str(nprocs) + " processor_ directories")

  print("Dump directory: " + str(dir_list))

  print("Number of processors: " + str(num_processors))

  dump_path = path + "/processor_" + str(1) + "/" + str(out_idx)

  print("Dump load path: " + dump_path)

  potential_files = os.listdir(dump_path)

  files = []

  for item in potential_files:
    m = re.search(r"\w+.dmp", item)
    if m and item != "time.dmp":
      files.append(item)

  return files

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def processDumpDir(nprocs, path, out_idx):

  print("Found " + str(nprocs) + " processor_ directories")

  print("Dump directory: " + str(dir_list))

  print("Number of processors: " + str(num_processors))

  dump_path = dump_folder_path + "/processor_" + str(1) + "/" + str(dump_folder_idx)

  print("Dump load path: " + dump_path)

  print("Backup files: " + str(os.listdir(dump_path)))

  header_line = ""
  element_entries = []

  fields = getDumpFields(nprocs, path, out_idx)

  for f in fields:
    combineField(num_processors, f, dump_folder_path, dump_folder_idx)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def usage():
  print("Usage: python combine_fields.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-d, --dump-path]: path to the dump folder")
  print("[-i, --idx]: index of the dump folder")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

# Parameters for reading and merging the field:
# A)Number of processors:
# - equal to number of processor_* directories in the folder
# B)the dump folder index

# The format of a .dmp file is as follows:
# - it is a list of coarse element entries
# - a coarse element entry consists of multiple lines
# - in the first line a single integer values denotes the global coarse element index
# - after the first line there is a line for each component of the output field
# - in each 'component line' the values of fine mesh dofs in the coarse element are given in the FEATFLOW multi-level ordering

# Script begin
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

try:
    opts, args = getopt.getopt(sys.argv[1:], 'i:d:h', ['idx=','dump-path=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

dump_folder_idx = 1

dump_folder_path = "./_dump"

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-i', '--idx'):
        dump_folder_idx = arg
    elif opt in ('-d', '--dump-path'):
        dump_folder_idx = arg
    else:
        usage()
        sys.exit(2)

dir_list = os.listdir(dump_folder_path)

num_processors = 0

for e in dir_list:
  m = re.search(r"processor_\d+", e)
  if m:
    num_processors = num_processors + 1

processDumpDir(num_processors, dump_folder_path, dump_folder_idx)
