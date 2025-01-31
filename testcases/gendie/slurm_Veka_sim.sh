#!/bin/sh -l

#SBATCH --time=08:00:00
#SBATCH --nodes=1 --ntasks-per-node=64 --cpus-per-task=1
#SBATCH --mem=1000GB
#SBATCH --partition=ext_math_prio
#SBATCH --mail-user=omierka@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=veka
#SBATCH --no-requeue

cd $SLURM_SUBMIT_DIR

module purge 
module load gcc/6.5.0 python/3.7.7 git intel/studio-xe/19.0.3.203 openmpi/mpi_thread_multiple/no_cuda/4.0.3 cmake/3.13.2

echo $SLURM_SUBMIT_DIR
echo `module list`


./RunnerGenDIE.sh


