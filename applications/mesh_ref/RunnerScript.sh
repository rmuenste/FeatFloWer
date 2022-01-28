#!/usr/bin/env bash

# NLMAX=3
SimFolder="VEKA_S"

# loop for the resolution

# create the area distribution into the same folder
python3 computeAreas.py ${SimFolder}

# we need to create an empty mesh_names.offs file that the meshrefiner can use it
echo "0" > mesh_names.offs

# run the mesh refinement binary
./meshref -f ${SimFolder}

# extract the number of generated elements
# NEL=`grep NEL ${SimFolder}/meshDir/Merged_Mesh.tri | awk -F " " '{print  $1}'`

# if NEL > 80,000 we have to go back to beginning and reset in the ${SimFolder}/info.json the "extrusion"@"minGap" to half value! ==> at the same time NLMAX = NLAMX + 1
# end loop for the resolution

# unfortunately we can practically afford only resolution 3 and 4, not higher! so NLMAX = MIN(4,NLMAX) ! but we have to save a warning signal that the resolution might be too low and pass it back to SR ...
# update SimPar@MaxMeshLevel = NLMAX in _data_BU/q2p1_paramV_DIE.dat and _data_BU/q2p1_paramV_DIE_test.dat

exit 0
