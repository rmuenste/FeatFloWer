#!/bin/bash
#SBATCH --job-name=gpuP_stage2
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --constraint=ampere   # = tardis (only ampere node); --nodelist=tardis is rejected by the job_submit plugin
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/stage2_%j.log
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.13.5 cuda/12.6
cd /data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100
export OMP_NUM_THREADS=32
SCR=/scratch/rmuenste/gpu-poisson; mkdir -p $SCR/tmp
INST=$SCR/p1_128x128x384.npz
echo "stage 2 on $(hostname), job $SLURM_JOB_ID, $(date -Is)"; nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv; df -h /scratch | tail -1; free -g | head -2
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 5 > results/stage2_nvsmi_${SLURM_JOB_ID}.csv &
NVPID=$!
set -x
# half size, FULL storage (cuSPARSE on row slabs < 2^31 nnz); build slab-wise on the node, instance to node-local scratch
python3 gpu_amg_cg.py --build 128 128 384 --save $INST --instance $INST --workers 8 --tmpdir $SCR/tmp \
        --storage full --block-jacobi --amg ${AMG_PRIMARY:-sa_bsr_lin},sa_csr_spec,geo_lin --tag stage2 --out results/gpu_p1_128x128x384.json
# same instance, symmetric storage (rehearsal of the stage-3 path at half size)
python3 gpu_amg_cg.py --instance $INST --storage sym --amg ${AMG_PRIMARY:-sa_bsr_lin} --tag stage2-sym --out results/stage2_p1_128x128x384_sym.json
set +x
kill $NVPID 2>/dev/null
ls -la $SCR; df -h /scratch | tail -1
echo "stage 2 done $(date -Is)"
