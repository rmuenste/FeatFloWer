#!/usr/bin/env bash

usage()
{
cat << EOF
usage: $0 options
 [./RunnerScript.sh -f e3d_setup ]

OPTIONS:
   -f      simulation folder name
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
   *) echo "Unrecognized flag or argument: $flag"
      usage
      exit 1
   ;;
 esac
 shift
done

rm -fr ${SimFolder}/Coarse_meshDir
rm -fr ${SimFolder}/meshDir_BU
rm -fr ${SimFolder}/meshDir
rm -fr ${SimFolder}/_vtk

./STLvsTRI -f ${SimFolder}

echo "1" > mesh_names.offs
echo  ${SimFolder}"/surface.off" >> mesh_names.offs

./meshref -f  ${SimFolder} -o meshDir

exit 0
