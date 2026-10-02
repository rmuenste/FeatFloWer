#!/bin/bash
#SBATCH --job-name=dns_e3_l4
#SBATCH --partition=long
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=60G
#SBATCH --prefer=nx
#SBATCH --time=06:00:00
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4/slurm-%j.out
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e3_l4
mpirun -np 32 /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/build-dns-pe-serial/applications/q2p1_bench_sedimentation/q2p1_bench_sedimentation > run_slurm.log 2>&1
