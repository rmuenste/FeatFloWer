

import os
import glob
import argparse

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description='Create a PVD file from a sequence of VTK files.')
parser.add_argument('--step', type=int, default=20, help='Timestep increment between files (default: 20)')
args = parser.parse_args()

# Directory where the PVD file will be saved
output_directory = '.'
pvd_filename = os.path.join(output_directory, 'main_files.pvd')

# Directory containing the files to be included in the PVD file
vtk_directory = '_vtk'

# Find all main.* files
file_pattern = os.path.join(vtk_directory, 'main.*.pvtu')
files = sorted(glob.glob(file_pattern))

# Generate PVD content
pvd_content = '<?xml version="1.0"?>\n'
pvd_content += '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
pvd_content += '<Collection>\n'

timestep = 0
for f in files:
    # The path in the PVD file should be relative to the PVD file itself.
    relative_path = os.path.relpath(f, output_directory)
    pvd_content += f'<DataSet timestep="{timestep}" part="0" file="{relative_path}"/>\n'
    timestep += args.step

pvd_content += '</Collection>\n'
pvd_content += '</VTKFile>\n'

# Write the PVD file
with open(pvd_filename, 'w') as f_out:
    f_out.write(pvd_content)

print(f"Successfully created {pvd_filename}")

