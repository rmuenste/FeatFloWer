#!/bin/bash
#SBATCH --job-name=gpuP_cleanup
#SBATCH --partition=med
#SBATCH --time=00:20:00
#SBATCH --constraint=ampere
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/cleanup.log
# Runs cleanup_scratch.sh on tardis's node-local /scratch after all results were copied to results/.
bash /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/cleanup_scratch.sh
