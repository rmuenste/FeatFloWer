#!/bin/bash
#SBATCH --job-name=gpuP_stage3
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --constraint=ampere   # = tardis (only ampere node); --nodelist=tardis is rejected by the job_submit plugin
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/stage3_%j.log
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.13.5 cuda/12.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100
export OMP_NUM_THREADS=32
SCR=/scratch/rmuenste/gpu-poisson; mkdir -p $SCR/tmp
INST=$SCR/p1_160x160x480.npz
echo "stage 3 on $(hostname), job $SLURM_JOB_ID, $(date -Is)"; nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv; df -h /scratch | tail -1; free -g | head -2
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 5 > results/stage3_nvsmi_${SLURM_JOB_ID}.csv &
NVPID=$!
set -x
# 3a: full size (D6.4 fine level), SYMMETRIC storage (upper triangle + warp kernel with f64 atomics), build slab-wise on the node
python3 gpu_amg_cg.py --build 160 160 480 --save $INST --instance $INST --workers 8 --tmpdir $SCR/tmp \
        --storage sym --amg ${AMG_PRIMARY:-sa_bsr_lin},geo_lin,sa_csr_spec --tag stage3-sym --out results/gpu_p1_160x160x480.json
# 3b: full size in FULL storage with the warp kernel (int64 row pointers; 61 GB matrix) -- expected to fit only with the small
#     cell-block hierarchy (~18 GB sa_bsr_lin, ~10 GB geo_lin; sa_csr_spec cannot run: pyamg core is int32-only); recorded either way (OOM is reported, not hidden)
python3 gpu_amg_cg.py --instance $INST --storage full --matvec kernel --amg ${AMG_PRIMARY:-sa_bsr_lin},geo_lin --tag stage3-full --out results/stage3_p1_160x160x480_full_kernel.json
set +x
kill $NVPID 2>/dev/null
ls -la $SCR; df -h /scratch | tail -1
echo "stage 3 done $(date -Is)"
