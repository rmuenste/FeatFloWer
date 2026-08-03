#!/bin/bash
#SBATCH --partition=long
#SBATCH --constraint=nx
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --job-name=ATC_rpm40

# Modules loaded explicitly (no ~/.bashrc aliases; matches guide 04 / build env)
source /etc/profile.d/modules.sh
module purge
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6

cd "$SLURM_SUBMIT_DIR"

# np = 63 CFD subdomains (_mesh/NEWFAC/sub0001) + 1 master
mpirun -np 64 ./q2p1_ATC | tee run_slurm.log
