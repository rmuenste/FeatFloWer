#!/usr/bin/env bash

usage()
{
cat << EOF
usage: $0 options
 [./RunnerScript.sh -f e3d_setup -n N]

OPTIONS:
   -f      simulation folder name
   -n      number of parallel processes
EOF
}

#SimFolder="GEALAN_6506"

if [ $# -eq 0 ]
then
  usage
  exit 1
fi

while [ $# != 0 ]; do
 flag="$1"
 case "$flag" in
  -f)  if [ $# -gt 1 ]; then
        SimFolder="$2"
        shift
       else
        echo "You did not provide an argument for the -f flag"
        usage
        exit 1
       fi
#       echo "You supplied an argument for the -nm flag: $basename"
   ;;
  -n)  if [ $# -gt 1 ]; then
        NProcs="$2"
        shift
       else
        echo "You did not provide an argument for the -n flag"
        usage
        exit 1
       fi
#       echo "You supplied an argument for the -nm flag: $basename"
   ;;
   *) echo "Unrecognized flag or argument: $flag"
      usage
      exit 1
   ;;
 esac
 shift
done

GendieProcs=$NProcs
MeshRefProcs=$NProcs

rm -fr ${SimFolder}/Coarse_meshDir
rm -fr ${SimFolder}/meshDir_BU
rm -fr ${SimFolder}/meshDir
rm -fr ${SimFolder}/_vtk

# for G-Elit we will have to activate these lines as well 

#cp ${SimFolder}"/surfaceX.off" ${SimFolder}"/surface.off"
./STLvsTRI -f ${SimFolder}
#cp ${SimFolder}"/surfaceS.off" ${SimFolder}"/surface.off"

echo "1" > mesh_names.offs
echo  ${SimFolder}"/surface.off" >> mesh_names.offs

./meshref -f  ${SimFolder}

python3 ./e3d_start.py -n ${MeshRefProcs} -f ${SimFolder} --mesh-reduction
rm -fr ReducedMeshDir
mv _vtk ${SimFolder}


# Running gendie as before:
cp _data_BU/q2p1_paramV_DIE_0.dat _data_BU/q2p1_paramV_DIE.dat
python3 ./e3d_start.py -n ${GendieProcs} -f ${SimFolder} --die-simulation

cp _data_BU/q2p1_paramAlpha.dat _data/q2p1_param.dat
cp _data_BU/mesh_names.offs .
mpirun -np ${GendieProcs} q1_scalar_multimat

cp _data_BU/q2p1_paramV_DIE_1.dat _data_BU/q2p1_paramV_DIE.dat
python3 ./e3d_start.py -n ${GendieProcs} -f ${SimFolder} --die-simulation


exit 0
