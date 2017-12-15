import os
import sys

newt_ok = False
non_newt_ok = False
visco_ok = False
sedimentation_ok = False

if not(len(sys.argv) == 2):
  sys.exit("Too few arguments, you need to pass the bin dir.")
else:
  bin_folder = sys.argv[1]

try:
  with open(bin_folder + '/applications/q2p1_fc_ext/note_single.json','r') as f:
    newt = f.read()
    newt_ok = True
except:
  print("Exception while reading file newt")

try:
  with open(bin_folder + '/applications/q2p1_fac_nnewt/note_single.json','r') as f:
    non_newt = f.read()
    non_newt_ok = True
except:
  print("Exception while reading file non_newt")


try:
  with open(bin_folder + '/applications/q2p1_fac_visco/note_single.json','r') as f:
    visco = f.read()
    visco_ok = True
except:
  print("Exception while reading file visco_newt")


try:
  with open(bin_folder + '/applications/q2p1_bench_sedimentation/note_single.json','r') as f:
    sedimentation = f.read()
    sedimentation_ok = True
except:
  print("Exception while reading file sedimentation")

with open(bin_folder + '/note.json','w') as f:
  f.write('[\n')
  if newt_ok == True:
    f.write(newt + ',')
  if non_newt_ok == True:
    f.write(non_newt + ',')
  if visco_ok == True:
    f.write(visco + ',')
  if sedimentation_ok == True:
    f.write(sedimentation)
  f.write(']\n')
