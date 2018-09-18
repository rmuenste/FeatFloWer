#!/usr/bin/env python
# vim: set filetype=python

import os
import sys
import re
import getopt
from shutil import copyfile

# Scan what is in the _dump/processor_X folder
# for each .dmp field found there
# read and join the partitioned files to
# a single file per field

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def cmp_entries(a, b):
# sort the tuples in the list element_entries
# according to the global element index
  if a[0] > b[0]:
    return 1
  elif a[0] == b[0]:
    return 0
  else:
    return -1

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def readPartitionedField(elem_entries, fileName):
  """ This routine reads a single field from the file <fileName>
      and stores it in elem_entries
      .dmp solution file into a combined single solution file.

      Args:
          fieldName: The name of the field that is processed
          elem_entries: The list of entries that is appended to
      Returns:
          header_line: The header line as a string
          header_info: The header as a list
          
  """

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

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def writeCombinedField(elem_entries, header, components, fileName):
  """ This routine handles the writing to a combined solution .dmp file
      from a list of element entries.

      Args:
          elem_entries: The list of entries that are written to the file
          header: The header info line
          components: The number of scalar components of the solution
          fieldName: The name of the combined output file
  """

  print("Writing combined files: " + str(fileName))

  with open(fileName, "w") as f:
    f.write(header)
    for e in elem_entries:
      f.write(str(e[0]) + '\n')
      for i in range(components):
        f.write(" ".join(e[1][i]) + '\n')

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def combineField(nprocs,fieldName, path, out_idx):
  """ This routine handles the reading and combining of a partitioned
      .dmp solution file into a combined single solution file.

      Args:
          nprocs: The number of processor directories for the solution
          fieldName: The name of the field that is processed
          path: The path to the dump files
          out_idx: The index of the solution directory
  """

  element_entries = []

  for i in range(1,nprocs+1):
    fileToRead = path + "/processor_" + str(i) + "/" + str(out_idx) + "/" + fieldName

    header_line, header_info = readPartitionedField(element_entries, fileToRead)

    element_entries.sort(cmp_entries)

    #print("Number of elements entries: " + str(len(element_entries)))

  del header_info['DofsTotal']

  header_info['NEL'] = len(element_entries)
  header_string = str(header_info)

  header_l=[]
  for k, v in header_info.iteritems():
    entry=[]
    entry.append(str(k))
    entry.append(str(v))
    str_entry = ":".join(entry)
    header_l.append(str_entry)

  header_mod = ",".join(header_l) 
  header_mod = header_mod + "\n"

  myPath = "./_dump/" + str(out_idx)

  if not os.path.exists(myPath):
    os.mkdir("./_dump/" + str(out_idx))

  outName = myPath + "/" + fieldName
  writeCombinedField(element_entries, header_mod, int(header_info['Components']), outName)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def getDumpFields(nprocs, path, out_idx):
  """ Scans the directory for .dmp dump files and adds
      them to a list.

      Args:
          nprocs: The number of processor directories for the solution
          path: The path to the dump files
          out_idx: The index of the solution directory
  """

  dump_path = path + "/processor_" + str(1) + "/" + str(out_idx)

  print("Dump load path: " + dump_path)

  potential_files = os.listdir(dump_path)

  files = []

  for item in potential_files:
    m = re.search(r"\w+.dmp", item)
    if m: 
      if item != "time.dmp":
        files.append(item)

  return files

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def processDumpDir(nprocs, path, out_idx):
  """ Entry point routine that directs the reading and combining of the
      dump folders and files.

      Args:
          nprocs: The number of processor directories for the solution
          path: The path to the dump files
          out_idx: The index of the solution directory
  """

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

  timeFile = dump_path + "/time.dmp"
  if os.path.exists(timeFile):
    print("Writing time file: " + dump_folder_path + "/" + str(dump_folder_idx) + "/time.dmp")
    copyfile(timeFile, dump_folder_path + "/" + str(dump_folder_idx) + "/time.dmp")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

def usage():
  print("Usage: python combine_fields.py [options]")
  print("Where options can be:")
  print("[-h, --help]: prints this message")
  print("[-d, --dump-path]: path to the dump folder, default is ./_dump")
  print("[-i, --idx]: index of the dump folder, default is 1")
  print("example: python ./combine_fields.py --dump-path=_dump --idx=2")

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

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
        dump_folder_path = arg
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
