#!/bin/sh
#SBATCH --time=01:00:00
# The following line asks for exclusive usage of 1 compute node using
# 20 cores per node for a non-threaded code. See 'man sbatch' for details.
#SBATCH --nodes=1 --ntasks-per-node=5 --cpus-per-task=1 --exclusive
#SBATCH --partition=short
# Possible 'constraint' values:
#               all, public
#               cgpu01  OR  gpu     OR  tesla_k40
#               cquad01 OR  xeon_e54640v4
#               cquad02 OR  xeon_e54640v4
#               cstd01  OR  xeon_e52640v4   OR  ib_1to3
#               cstd02  OR  xeon_e52640v4   OR  ib_1to1 OR  nonblocking_comm
#SBATCH --constraint="cstd01"
# Maximum 'mem' values depending on constraint (values in MB):
#               cstd01/xeon_e52640v4/ib_1to3 AND
#               cstd02/xeon_e52640v4/ib_1to1/nonblocking_comm: 62264
#               cquad01: 255800
#               cquad02: 1029944
#SBATCH --mem=62264
#SBATCH --mem_bind=verbose,local --hint=memory_bound
#SBATCH --mail-user=raphael.muenster@math.tu-dortmund.de
# Possible 'mail-type' values: NONE, BEGIN, END, FAIL, ALL (= BEGIN,END,FAIL)
#SBATCH --mail-type=BEGIN,END,FAIL

# User specific aliases and functions
module purge && module load cmake gcc/6.2.0 openmpi/no_mpi_thread_multiple/no_cuda/2.1.2 && export CC=mpicc && export CXX=mpicxx && export FC=mpif90	

cd /home/smramuen/nobackup/FeatFlower/bin/applications/q2p1_fc_ext

mpirun -np 5 ./q2p1_fc_ext
