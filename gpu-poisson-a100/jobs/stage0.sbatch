#!/bin/bash
#SBATCH --job-name=gpuP_stage0
#SBATCH --partition=med
#SBATCH --time=00:20:00
#SBATCH --constraint=ampere   # tardis is the only ampere node; --nodelist=tardis is rejected by the job_submit plugin
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/stage0_%j.log
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load python/3.13.5 cuda/12.6
OUT=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/environment.txt
{
echo "# stage 0 environment, $(date -Is), host $(hostname), job $SLURM_JOB_ID"
echo "## modules"; module list 2>&1
echo "## nvidia-smi"; nvidia-smi
nvidia-smi --query-gpu=name,memory.total,memory.used,driver_version,compute_cap --format=csv
echo "## python"; which python3; python3 --version
python3 - <<'PY'
import cupy, pyamg, numpy, scipy, sys
print("cupy", cupy.__version__, "pyamg", pyamg.__version__, "numpy", numpy.__version__, "scipy", scipy.__version__)
d = cupy.cuda.Device(0); print("compute_capability", d.compute_capability)
free, total = d.mem_info; print(f"device mem free {free/2**30:.1f} GiB total {total/2**30:.1f} GiB")
print("cuda runtime", cupy.cuda.runtime.runtimeGetVersion(), "driver", cupy.cuda.runtime.driverGetVersion())
import cupyx.scipy.sparse as csp
x = cupy.arange(10.0); A = csp.identity(10, format='csr'); print("spmv ok", float((A @ x).sum()))
PY
echo "## host"; nproc; free -g; cat /proc/cpuinfo | grep "model name" | head -1
echo "## filesystems"; df -h 2>/dev/null | grep -Ev "tmpfs|overlay|^none|/boot|/snap"; echo; ls -ld /scratch /tmp /local* /data/scratch* 2>&1
echo "## scratch write test"; for d in /scratch /tmp; do if [ -d $d ]; then mkdir -p $d/rmuenste 2>&1 && echo "$d/rmuenste writable: $(touch $d/rmuenste/.probe && echo yes || echo no) $(df -h $d | tail -1)"; fi; done
} > "$OUT" 2>&1
cat "$OUT"
