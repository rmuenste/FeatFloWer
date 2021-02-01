#!/usr/bin/env bash

# NLMAX=3
SimFolder="VEKA_S"

# direct the input and outp folder of the meshrefinemnt binary to the SimFolder
cp param_BU.cfg param.cfg
sed -i 's/SimFolder/'${SimFolder}'/g' param.cfg

# loop for the resolution

# create the area distribution into the same folder
python3 computeAreas.py ${SimFolder}

# we need to create an empty mesh_names.offs file that the meshrefiner can use it
echo "0" > mesh_names.offs

# run the mesh refinement binary
./meshref

# extract the number of generated elements
NEL=`grep NEL ${SimFolder}/meshDir/Merged_Mesh.tri | awk -F " " '{print  $1}'`

# if NEL > 80,000 we have to go back to beginning and reset in the ${SimFolder}/info.json the "extrusion"@"minGap" to half value! ==> at the same time NLMAX = NLAMX + 1
# end loop for the resolution

# unfortunately we can practically afford only resolution 3 and 4, not higher! so NLMAX = MIN(4,NLMAX) ! but we have to save a warning signal that the resolution might be too low and pass it back to SR ...
# update SimPar@MaxMeshLevel = NLMAX in _data_BU/q2p1_paramV_DIE.dat and _data_BU/q2p1_paramV_DIE_test.dat

#NP=16
# nFineElem = $NEL*8^(NLMAX-1) is needed to determine the number of processors for the job == > $NP

#python3 ./e3d_start.py -n $NP -f ${SimFolder} --die-simulation

exit 0
