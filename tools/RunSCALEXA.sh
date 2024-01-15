#!/usr/bin/env bash

#configure the simulation 
##############################################################
#set the configuration folder to the 'SMALL' or 'LARGE' case
f='SMALL'

#set the number of MPI processes 
n=16

# set the required MG mesh resolution 
# for the SMALL case :: testing at 2,3; sharp sim at 4,5
# for the LARGE case :: testing at 2; sharp sim at 3,4,5
XMSHLVEL=3

# this parameter might be useful for the recursive partitioning
# in case of using a cluster with nodes having the same number of cores then set here the number of nodes
# the partitioning will first decompose the mesh to this number of submeshes and then decompose further 
# given number of sub-submeshes.

# for example there are 72 cores per node then in case of running the code on 5 nodes 
# the XSUBMESHNUMBER = 5 and n=5x72=360
# by default it is set to 1 and then it operates almost as before (almost, because the order of the 
# partitions is reversed)

XSUBMESHNUMBER=1

###############################################################

# copy the parameter file for the simulation
cp SCALEXA/q2p1_paramV_DIE_0.dat _data_BU/q2p1_paramV_DIE.dat

# set the resolution level to the required value
sed -i s/XMSHLVEL/${XMSHLVEL}/g _data_BU/q2p1_paramV_DIE.dat

# set the recursive partitioning for the data file and for the partitioner, as well
sed -i s/XSUBMESHNUMBER/${XSUBMESHNUMBER}/g _data_BU/q2p1_paramV_DIE.dat
sed -i '0,/partitionerParameters = \[[^]]*\]/{s/partitionerParameters = \[[^]]*\]/partitionerParameters = [1,'${XSUBMESHNUMBER}']/}' e3d_start.py

./RankFileGenerator.sh -n ${n} -s ${XCOREPERNODES}

# run the simulation
python3 ./e3d_start.py -n $n -f SCALEXA/$f --die-simulation


