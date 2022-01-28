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
load-gcc-mpi9-py;

cd $SLURM_SUBMIT_DIR

#python bench_driver.py -p 16 --build-id=epyc16core-linux-gcc-release --module-conf=gcc81 --bin-dir=bin-epyx
#python ./bench_driver.py -p $SLURM_NTASKS -c $SLURM_CONSTRAINT
python ./bench_driver.py -p $SLURM_NTASKS -c bttf --module-conf=gcc9 -d bin
