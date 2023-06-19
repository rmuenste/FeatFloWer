#!/usr/bin/env bash

#configure the simulation 
##############################################################
#set the configuration folder to the 'SMALL' or 'LARGE' case
f='SMALL'

#set the number of MPI processes 
n=16

#set the required MG mesh resolution 
# for the SMALL case :: testing at 2,3; sharp sim at 4,5
# for the LARGE case :: testing at 2; sharp sim at 3,4,5
XMSHLVEL=3

###############################################################

# copy the parameter file for the simulation
cp SCALEXA/q2p1_paramV_DIE_0.dat _data_BU/q2p1_paramV_DIE.dat

# set the resolution level to the required value
sed -i s/XMSHLVEL/${XMSHLVEL}/g _data_BU/q2p1_paramV_DIE.dat

# run the simulation
python3 ./e3d_start.py -n $n -f SCALEXA/$f --die-simulation


