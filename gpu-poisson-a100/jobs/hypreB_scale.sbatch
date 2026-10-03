#!/bin/bash
#SBATCH --job-name=gpuP_hypreB
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --constraint=ampere   # = tardis
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/hypreB_%j.log
# Part C(2,3) of the BoomerAMG follow-up: half size (more variants), FULL size attempt (expected not to fit in device memory
# beside the 63 GB matrix -- the OOM is recorded precisely), then the 3/4 grid 150x150x450 (built with the harness slab builder)
# so BoomerAMG has a number at the largest size that fits. Ranks per instance keep the local nnz < 2^30 (mixed-int hypre).
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load gcc/13.2.0 cuda/12.6 openmpi/4.1.6 python/3.13.5
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader OMP_NUM_THREADS=1
GP=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100; RES=$GP/results; JOB=$SLURM_JOB_ID
SCR=/scratch/rmuenste/gpu-poisson; DRV=$SCR/hypre_boomer_pcg
cd $GP
echo "hypre B (scale) on $(hostname), job $JOB, $(date -Is)"; nvidia-smi --query-gpu=name,memory.total,memory.used --format=csv; df -h /scratch | tail -1; free -g | head -2; ls -la $SCR $SCR/bin
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 2 > $RES/hypreB_nvsmi_${JOB}.csv &
NVPID=$!
HY=$SCR/hypre-cuda
# recompile the driver (Chebyshev option added after job A)
mpicc -O2 -Wall -Wno-unused-result -o $DRV hypre_boomer_pcg.c -I$HY/include -I$CUDA_HOME/include -L$HY/lib -lHYPRE -L$CUDA_HOME/lib64 -lcudart -lcusparse -lcurand -lcublas -lm -Wl,-rpath,$HY/lib -Wl,-rpath,$CUDA_HOME/lib64 2>&1 | grep -v "^$" | head; ls -la $DRV || { kill $NVPID; exit 1; }
run() { # label nranks flags...
  local L=$1 N=$2; shift 2; echo; echo "===== [$INST] variant $L ($N rank(s)) : $* : $(date -Is)"; local T=$(date +%s)
  timeout 5400 mpirun -np $N --bind-to none $DRV --bin $B --out ${P}${L}.json --label $L "$@" 2>&1 | grep -v "^ *[0-9]\+ \+[0-9.]\+e[+-][0-9]\+"
  echo "variant $L rc=${PIPESTATUS[0]} ($(( $(date +%s) - T )) s)"; nvidia-smi --query-gpu=memory.used --format=csv,noheader
}
# The variant set at scale (chosen from the 40x40x120 sweep in job A): PRIMARY = the brief's GPU settings with num_functions 4;
# BEST = the best-converging device variant of the sweep; LEAN = the smallest-footprint variant for the full size.
PRIMARY="${PRIMARY:-}"                       # default flags of the driver (PMIS, ext+i, l1-Jacobi, agg 1, nf 4)
BEST="${BEST:---agg-levels 0 --theta 0.5 --relax 16 --cheby-order 2}"   # provisional; job-B step 0 refines the choice
LEAN="${LEAN:---keep-transpose 0}"

# ---- 0. small instance: remaining smoother/coarsening combinations (1 rank) and the multi-rank cost, with and without MPS ----
INST=p1_40x40x120; B=$SCR/bin/$INST; P=$RES/hypre_${INST}__
run nf4_agg0_cheby2          1 --agg-levels 0 --relax 16 --cheby-order 2
run nf4_agg0_cheby3          1 --agg-levels 0 --relax 16 --cheby-order 3
run nf4_agg0_theta05_cheby2  1 --agg-levels 0 --theta 0.5 --relax 16 --cheby-order 2
run nf4_agg0_theta05_sweeps2 1 --agg-levels 0 --theta 0.5 --sweeps 2
run nf4_agg0_jac7_sweeps2    1 --agg-levels 0 --relax 7 --sweeps 2
run nf4_agg0_theta05_jac7_sw2 1 --agg-levels 0 --theta 0.5 --relax 7 --sweeps 2
run nf4_agg0_theta05_pmax8   1 --agg-levels 0 --theta 0.5 --pmax 8
run nf4_agg0_theta05_cheby2_pmax8 1 --agg-levels 0 --theta 0.5 --relax 16 --cheby-order 2 --pmax 8
run nf4_agg0_interp15        1 --agg-levels 0 --interp 15   # adaptive weights (GPU only)
run nf4_agg0_theta05_interp15_cheby2 1 --agg-levels 0 --theta 0.5 --interp 15 --relax 16 --cheby-order 2
run nf1_agg0_theta05_cheby2  1 --agg-levels 0 --theta 0.5 --relax 16 --cheby-order 2 --num-functions 1
run nf4_3ranks               3 $BEST
# MPS (lets the ranks' kernels run concurrently instead of time-slicing); user-level daemon, if the binary exists on the node
MPSCTL=$(command -v nvidia-cuda-mps-control || ls /usr/bin/nvidia-cuda-mps-control 2>/dev/null | head -1)
if [ -n "$MPSCTL" ]; then
  export CUDA_MPS_PIPE_DIRECTORY=$SCR/mps/pipe CUDA_MPS_LOG_DIRECTORY=$SCR/mps/log; mkdir -p $CUDA_MPS_PIPE_DIRECTORY $CUDA_MPS_LOG_DIRECTORY
  $MPSCTL -d && echo "MPS daemon started ($MPSCTL)" || echo "MPS daemon failed to start"
  run nf4_3ranks_mps 3 $BEST
  run nf4_1rank_mps  1 $BEST
  if grep -q "variant nf4_3ranks_mps rc=0" $RES/hypreB_${JOB}.log 2>/dev/null; then echo "keeping MPS on for the scale runs"; else echo quit | $MPSCTL; unset CUDA_MPS_PIPE_DIRECTORY CUDA_MPS_LOG_DIRECTORY; echo "MPS off"; fi
else echo "no nvidia-cuda-mps-control on this node"; fi
python3 hypre_merge_json.py --out $RES/hypre_${INST}.json --primary nf4_pmis_exti_l1jac ${P}*.json

# ---- 1. half size (bins from job A), 3 ranks ----
INST=p1_128x128x384; B=$SCR/bin/$INST; P=$RES/hypre_${INST}__
if [ -f $B.meta.txt ]; then
  run best      3 $BEST
  run nf1_best  3 $BEST --num-functions 1
  run best_jac7 3 $BEST --relax 7
  run lean      3 $LEAN
  python3 hypre_merge_json.py --out $RES/hypre_${INST}.json --primary nf4_pmis_exti_l1jac ${P}*.json
fi

# ---- 2. FULL size 160x160x480 (49.15 M unknowns, 5.21 G nnz = 62.9 GB CSR), 6 ranks ----
INST=p1_160x160x480; B=$SCR/bin/$INST; P=$RES/hypre_${INST}__
if [ ! -f $B.meta.txt ]; then python3 npz_to_bin.py --npz $SCR/$INST.npz --out $B 2>&1 | tail -3; fi
df -h /scratch | tail -1
run jacobi   6 --jacobi --no-amg                      # matrix + vectors only: must fit (63 GB + ~3 GB)
run primary  6 $PRIMARY                               # BoomerAMG setup beside the 63 GB matrix: fits or OOM (recorded)
run best     6 $BEST
run lean     6 $LEAN
python3 hypre_merge_json.py --out $RES/hypre_${INST}.json --primary primary ${P}*.json

# ---- 3. 3/4 grid 150x150x450 (40.5 M unknowns, ~4.3 G nnz = ~52 GB CSR): build (harness slab builder), dump, 5 ranks ----
INST=p1_150x150x450
if [ ! -f $SCR/$INST.npz ]; then
  python3 gpu_amg_cg.py --build 150 150 450 --save $SCR/$INST.npz --workers 8 --tmpdir $SCR/tmp --storage full --matvec kernel --amg none --no-jacobi --tag build-3q --out $RES/gpu_${INST}_build.json 2>&1 | grep -v "^\[.*\[jacobi"
fi
B=$SCR/bin/$INST; P=$RES/hypre_${INST}__
if [ ! -f $B.meta.txt ]; then python3 npz_to_bin.py --npz $SCR/$INST.npz --out $B 2>&1 | tail -3; fi
df -h /scratch | tail -1
run jacobi   5 --jacobi --no-amg
run primary  5 $PRIMARY
run best     5 $BEST
run lean     5 $LEAN
python3 hypre_merge_json.py --out $RES/hypre_${INST}.json --primary primary ${P}*.json

# ---- 4. control experiment: a BIGINT build (HYPRE_Int 64-bit, 16 B/nnz) runs the half size on ONE rank -> isolates the cost of
#         time-slicing the GPU between ranks (mixed-int needs 3 ranks there). Guarded: if hypre refuses bigint+CUDA, it is recorded.
HYB=$SCR/hypre-cuda-bigint; DRVB=$SCR/hypre_boomer_pcg_bigint
if [ ! -f $HYB/lib/libHYPRE.so ]; then
  mkdir -p $SCR/hypre-src-bigint && cp -r $SCR/hypre-src/src $SCR/hypre-src-bigint/ && cd $SCR/hypre-src-bigint/src && make distclean > /dev/null 2>&1
  CONFB="./configure --prefix=$HYB --with-cuda --with-cuda-home=$CUDA_HOME --with-gpu-arch=80 --enable-unified-memory=no --enable-shared --enable-bigint --enable-cusparse --enable-curand --enable-cublas --without-umpire --with-MPI CC=mpicc CXX=mpicxx"
  echo "bigint configure: $CONFB"; $CONFB > $RES/hypre_configure_bigint_${JOB}.log 2>&1; echo "configure rc=$?"; grep -n "HYPRE_BIGINT\|HYPRE_MIXEDINT\|HYPRE_USING_CUDA " HYPRE_config.h
  T=$(date +%s); make -j32 > $RES/hypre_make_bigint_${JOB}.log 2>&1; echo "make rc=$? ($(( $(date +%s) - T )) s)"; grep -i "error:" $RES/hypre_make_bigint_${JOB}.log | head -5; make install >> $RES/hypre_make_bigint_${JOB}.log 2>&1; echo "install rc=$?"
  cd $GP
fi
if [ -f $HYB/lib/libHYPRE.so ]; then
  mpicc -O2 -Wall -Wno-unused-result -o $DRVB hypre_boomer_pcg.c -I$HYB/include -I$CUDA_HOME/include -L$HYB/lib -lHYPRE -L$CUDA_HOME/lib64 -lcudart -lcusparse -lcurand -lcublas -lm -Wl,-rpath,$HYB/lib -Wl,-rpath,$CUDA_HOME/lib64 2>&1 | grep -v "^$" | head -5
  { echo "## bigint control build (job $JOB): $CONFB"; grep "define HYPRE_BIGINT\|define HYPRE_USING_CUDA " $HYB/include/HYPRE_config.h; } >> $RES/hypre_build.txt
  DRV=$DRVB; INST=p1_128x128x384; B=$SCR/bin/$INST; P=$RES/hypre_${INST}__bigint_
  run jacobi_1rank 1 --jacobi --no-amg
  run best_1rank   1 $BEST
  run primary_1rank 1 $PRIMARY
  python3 hypre_merge_json.py --out $RES/hypre_${INST}_bigint_1rank.json --primary best_1rank ${P}*.json
else echo "bigint build not available"; fi

kill $NVPID 2>/dev/null
ls -la $SCR $SCR/bin; df -h /scratch | tail -1; du -sh $SCR
echo "hypre B done $(date -Is)"
