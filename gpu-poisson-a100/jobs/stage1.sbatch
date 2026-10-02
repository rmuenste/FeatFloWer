#!/bin/bash
#SBATCH --job-name=gpuP_stage1
#SBATCH --partition=med
#SBATCH --time=03:00:00
#SBATCH --constraint=ampere   # = tardis (only ampere node); --nodelist=tardis is rejected by the job_submit plugin
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/stage1_%j.log
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.13.5 cuda/12.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100
export OMP_NUM_THREADS=16
echo "stage 1 on $(hostname), job $SLURM_JOB_ID, $(date -Is)"; nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 5 > results/stage1_nvsmi_${SLURM_JOB_ID}.csv &
NVPID=$!
ALL=sa_csr_spec,sa_csr_mean,sa_bsr_mean,sa_bsr_lin,geo_lin,geo_lin_jac
set -x
python3 gpu_amg_cg.py --instance instances/p1_20x20x60.npz  --storage full --validate-sym --rough-rhs --block-jacobi --amg $ALL --pyamg-cpu-solve --tag stage1 --out results/gpu_p1_20x20x60.json
python3 gpu_amg_cg.py --instance instances/p1_40x40x120.npz --storage full --validate-sym --rough-rhs --block-jacobi --amg $ALL --pyamg-cpu-solve --tag stage1 --out results/gpu_p1_40x40x120.json
python3 gpu_amg_cg.py --instance instances/p1_40x40x120.npz --storage full --matvec kernel --amg geo_lin --tag stage1-full-kernel --out results/stage1_p1_40x40x120_full_kernel.json
python3 gpu_amg_cg.py --instance instances/p1_40x40x120.npz --storage sym --amg geo_lin,sa_bsr_lin --tag stage1-sym --out results/stage1_p1_40x40x120_sym.json
python3 gpu_amg_cg.py --instance instances/p1_40x40x120.npz --storage full --no-jacobi --amg geo_lin,sa_bsr_lin --omega-factor 1.0 --tag stage1-omega1 --out results/stage1_p1_40x40x120_omega1.json
set +x
kill $NVPID 2>/dev/null
echo "stage 1 done $(date -Is)"
