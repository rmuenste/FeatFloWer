#!/bin/bash
#SBATCH --job-name=gpuP_stage3c
#SBATCH --partition=med
#SBATCH --time=02:00:00
#SBATCH --constraint=ampere   # = tardis
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/stage3c_%j.log
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.13.5 cuda/12.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100
export OMP_NUM_THREADS=32
INST=/scratch/rmuenste/gpu-poisson/p1_160x160x480.npz
echo "stage 3c (full storage at full size, geo_lin only, clean memory) on $(hostname), job $SLURM_JOB_ID, $(date -Is)"
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 5 > results/stage3c_nvsmi_${SLURM_JOB_ID}.csv &
NVPID=$!
set -x
python3 gpu_amg_cg.py --instance $INST --storage full --matvec kernel --no-jacobi --amg geo_lin --tag stage3-full-geo --out results/stage3_p1_160x160x480_full_kernel_geo_lin.json
set +x
kill $NVPID 2>/dev/null
echo "stage 3c done $(date -Is)"
