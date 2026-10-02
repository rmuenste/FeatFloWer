#!/bin/bash
#SBATCH --job-name=dns_e1_l3_g3def
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --prefer=nx
#SBATCH --time=02:00:00
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def/slurm-%j.out
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/q2p1_dns_rundir_e1_l3_g3def
mpirun -np 32 /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/build-twin-pe6e5ba9c/applications/q2p1_bench_sedimentation/q2p1_bench_sedimentation > run_slurm.log 2>&1
