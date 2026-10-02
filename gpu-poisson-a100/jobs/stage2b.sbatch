#!/bin/bash
#SBATCH --job-name=gpuP_stage2b
#SBATCH --partition=med
#SBATCH --time=04:00:00
#SBATCH --constraint=ampere   # = tardis
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/stage2b_%j.log
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.13.5 cuda/12.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100
export OMP_NUM_THREADS=32
INST=/scratch/rmuenste/gpu-poisson/p1_128x128x384.npz
echo "stage 2b (full storage rerun, slabs <= 2^30 nnz) on $(hostname), job $SLURM_JOB_ID, $(date -Is)"; ls -la $INST
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 5 > results/stage2b_nvsmi_${SLURM_JOB_ID}.csv &
NVPID=$!
set -x
python3 gpu_amg_cg.py --instance $INST --storage full --matvec cusparse --block-jacobi --amg ${AMG_PRIMARY:-sa_bsr_lin},sa_csr_spec,geo_lin --tag stage2 --out results/gpu_p1_128x128x384.json
python3 gpu_amg_cg.py --instance $INST --storage full --matvec kernel --no-jacobi --amg ${AMG_PRIMARY:-sa_bsr_lin} --tag stage2-full-kernel --out results/stage2_p1_128x128x384_full_kernel.json
set +x
kill $NVPID 2>/dev/null
echo "stage 2b done $(date -Is)"
