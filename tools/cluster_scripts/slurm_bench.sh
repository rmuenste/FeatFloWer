#!/bin/bash

#SBATCH --partition=med
#SBATCH --constraint=bttf
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --job-name=FFBench

shopt -s expand_aliases
source ~/.bashrc 

# User specific aliases and functions
load-gcc-mpi81;
USER_BUILD_NAME="xeon-linux-gcc81-release"

echo $USER_BUILD_NAME
cd $SLURM_SUBMIT_DIR

ctest -S ctest_driver.cmake -DBUILD_STRING=xeon-linux-gcc-release -DUSER_BUILD_NAME=$USER_BUILD_NAME -DPROCS=$SLURM_NTASKS -DSRC_DIR=$SLURM_SUBMIT_DIR/Feat_FloWer -DBIN_DIR=$SLURM_SUBMIT_DIR/bin-bttf-gcc81 -VV

load-intel-mpi18;

USER_BUILD_NAME="xeon-linux-intel18-release"
ctest -S ctest_driver.cmake -DBUILD_STRING=xeon-linux-intel-release -DUSER_BUILD_NAME=$USER_BUILD_NAME -DPROCS=$SLURM_NTASKS -DSRC_DIR=$SLURM_SUBMIT_DIR/Feat_FloWer -DBIN_DIR=$SLURM_SUBMIT_DIR/bin-bttf-intel18 -VV
