import os
import fnmatch
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-n', '--name', type=str, default="output.pvd", help='The output name')
parser.add_argument('-p', '--path', type=str, default="./_vtk", help='The input directory path')
parser.add_argument('-f', '--frequency', type=int, default=20, help='The output frequency')

args = parser.parse_args()

# Look for files matching the pattern in the _vtk directory
files = []

out_name = args.name

dir_path = args.path 

if os.path.exists(dir_path) and os.path.isdir(dir_path):
    for file in os.listdir(dir_path):
        if fnmatch.fnmatch(file, "main.*.pvtu"):
            files.append(os.path.join(dir_path, file))

# Sort the files by name to ensure they are in the correct order
files.sort()
freq = args.frequency

# Generate the XML file
with open(out_name, "w") as f:
    f.write('<?xml version="1.0"?>\n')
    f.write('<!-- This pvd file references all the written time steps. -->\n')
    f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
    f.write('<Collection>\n')
    timestep = 0
    for file in files:
        f.write('<DataSet timestep="%d" part="0" file="%s"/>\n' % (timestep, file))
        timestep += freq
    f.write('</Collection>\n')
    f.write('</VTKFile>\n')
