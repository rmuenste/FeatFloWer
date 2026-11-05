#!/bin/bash
#SBATCH --job-name=gpuP_hypreA
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --constraint=ampere   # = tardis (only ampere node); --nodelist=tardis is rejected by the job_submit plugin
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8            # MPI slots for the multi-rank runs (ranks time-share the one GPU)
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --output=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100/results/hypreA_%j.log
# Part A/B/C(1) of the BoomerAMG follow-up: build CUDA hypre 3.1.0 on node-local scratch, smoke test (test/ij),
# compile the driver, run the 40x40x120 correctness variants (device + host), dump the half-size instance and run
# its first BoomerAMG numbers.
set -u
source /etc/profile.d/modules.sh 2>/dev/null
module load gcc/13.2.0 cuda/12.6 openmpi/4.1.6 python/3.13.5
export OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,vader   # single node: shared memory only
export OMP_NUM_THREADS=1
GP=/data/warehouse17/rmuenste/code/FF-EL/FeatFloWer/gpu-poisson-a100; RES=$GP/results; JOB=$SLURM_JOB_ID
SCR=/scratch/rmuenste/gpu-poisson; mkdir -p $SCR/bin
HY=$SCR/hypre-cuda; DRV=$SCR/hypre_boomer_pcg
TARDIR=/data/warehouse17/rmuenste/gpu-poisson-hypre; TAR=$TARDIR/hypre-v3.1.0-src.tar.gz
cd $GP
echo "hypre A (build + small instance) on $(hostname), job $JOB, $(date -Is)"; module list 2>&1; which mpicc nvcc gcc python3
nvidia-smi --query-gpu=name,memory.total,memory.used,driver_version --format=csv; df -h /scratch | tail -1; free -g | head -2; ls -la $SCR
nvidia-smi --query-gpu=timestamp,memory.used,utilization.gpu --format=csv -l 2 > $RES/hypreA_nvsmi_${JOB}.csv &
NVPID=$!

# ---- 1. source on scratch (tarball staged from the login node; git clone as fallback) ----
if [ ! -f $SCR/hypre-src/src/configure ]; then
  if [ -f $TAR ]; then cp $TAR $SCR/ && tar xzf $SCR/hypre-v3.1.0-src.tar.gz -C $SCR && rm -f $SCR/hypre-v3.1.0-src.tar.gz && rm -rf $TARDIR && echo "source from tarball ($TAR), staging dir removed"
  else (cd $SCR && git clone -q --depth 1 --branch v3.1.0 https://github.com/hypre-space/hypre.git hypre-src) && echo "source cloned on the node"; fi
fi
ls $SCR/hypre-src/src/configure || { echo "NO SOURCE"; kill $NVPID; exit 1; }

# ---- 2. configure + build + install ----
CONFIGURE_LINE="./configure --prefix=$HY --with-cuda --with-cuda-home=$CUDA_HOME --with-gpu-arch=80 --enable-unified-memory=no --enable-shared --enable-mixedint --enable-bigint=no --enable-cusparse --enable-curand --enable-cublas --without-umpire --with-MPI CC=mpicc CXX=mpicxx"   # hypre 3.x auto-enables Umpire for GPU builds and errors without its paths; opted out (plain cudaMalloc)
if [ ! -f $HY/lib/libHYPRE.so ]; then
  cd $SCR/hypre-src/src
  echo "configure: $CONFIGURE_LINE"; T=$(date +%s)
  $CONFIGURE_LINE > $RES/hypre_configure_${JOB}.log 2>&1; echo "configure rc=$? ($(( $(date +%s) - T )) s)"; tail -5 $RES/hypre_configure_${JOB}.log
  grep -n "HYPRE_USING_CUDA\|HYPRE_MIXEDINT\|HYPRE_BIGINT\|UNIFIED\|HYPRE_USING_CUSPARSE\|HYPRE_USING_CURAND\|HYPRE_USING_CUBLAS\|HYPRE_HAVE_MPI" HYPRE_config.h
  T=$(date +%s); make -j32 > $RES/hypre_make_${JOB}.log 2>&1; echo "make rc=$? ($(( $(date +%s) - T )) s)"; grep -i " error\|error:" $RES/hypre_make_${JOB}.log | head -5
  make install >> $RES/hypre_make_${JOB}.log 2>&1; echo "install rc=$?"
  ls -la $HY/lib | head; du -sh $HY
  cd $GP
fi

# ---- 3. build record ----
{
  echo "# CUDA hypre build record, job $JOB, $(date -Is), host $(hostname)"
  echo "## source"; echo "hypre v3.1.0, tag v3.1.0 = commit 9dc9e18aed6a945a95f966e57daacfb1c269f6ec (Release 3.1.0, 2026-01-23), github.com/hypre-space/hypre; source tree $SCR/hypre-src (scratch), install $HY"
  echo "## configure line (run in $SCR/hypre-src/src)"; echo "$CONFIGURE_LINE"
  echo "## toolchain"; module list 2>&1; gcc --version | head -1; nvcc --version | tail -2; mpicc --version | head -1; mpirun --version | head -1; ompi_info 2>/dev/null | grep -i "configure command line" | head -1
  nvidia-smi --query-gpu=name,driver_version,compute_cap --format=csv
  echo "## HYPRE_config.h (defined)"; grep "^#define" $HY/include/HYPRE_config.h
  echo "## HYPRE_config.h (GPU/int related, all)"; grep -n "CUDA\|GPU\|CUSPARSE\|CURAND\|CUBLAS\|UNIFIED\|DEVICE\|MIXEDINT\|BIGINT\|MPI\|OPENMP" $HY/include/HYPRE_config.h
  echo "## installed"; ls -la $HY/lib; du -sh $HY $SCR/hypre-src
} > $RES/hypre_build.txt 2>&1
echo "build record -> $RES/hypre_build.txt"

# ---- 4. smoke test: hypre's own ij driver on the GPU ----
cd $SCR/hypre-src/src/test && make ij > $RES/hypre_make_ij_${JOB}.log 2>&1; echo "make test/ij rc=$?"; ls -la ./ij
echo "== ij -solver 1 -n 50 50 50 (AMG-PCG, 7-pt Laplacian, device memory/exec by default in a GPU build)"
timeout 600 mpirun -np 1 ./ij -solver 1 -n 50 50 50 -print_system 0 2>&1 | tee $RES/hypre_ij_smoke_${JOB}.log | grep -i "Iterations\|Final Relative\|Complexity\|wall clock\|memory\|error\|CUDA" | head -20
echo "== ij -solver 1 -n 100 100 100 -agg_nl 1 -pmis -interptype 6 -relax 18"
timeout 600 mpirun -np 1 ./ij -solver 1 -n 100 100 100 -agg_nl 1 -pmis -interptype 6 -rlx 18 2>&1 | tee -a $RES/hypre_ij_smoke_${JOB}.log | grep -i "Iterations\|Final Relative\|Complexity\|wall clock\|error\|CUDA" | head -12
cd $GP

# ---- 5. compile the driver ----
CL="mpicc -O2 -Wall -Wno-unused-result -o $DRV hypre_boomer_pcg.c -I$HY/include -I$CUDA_HOME/include -L$HY/lib -lHYPRE -L$CUDA_HOME/lib64 -lcudart -lcusparse -lcurand -lcublas -lm -Wl,-rpath,$HY/lib -Wl,-rpath,$CUDA_HOME/lib64"
echo "driver compile: $CL"; $CL 2>&1 | grep -v "^$" | head -20; ls -la $DRV || { echo "DRIVER BUILD FAILED"; kill $NVPID; exit 1; }
{ echo "## driver compile line (hypre_boomer_pcg.c -> $DRV)"; echo "$CL"; } >> $RES/hypre_build.txt

# ---- 6. small instance: dump + variants ----
python3 npz_to_bin.py --npz instances/p1_40x40x120.npz --out $SCR/bin/p1_40x40x120 2>&1 | tail -2
B=$SCR/bin/p1_40x40x120; P=$RES/hypre_p1_40x40x120__
run() { # label nranks flags...
  local L=$1 N=$2; shift 2; echo; echo "===== variant $L ($N rank(s)) : $*"; local T=$(date +%s)
  timeout 2400 mpirun -np $N --bind-to none $DRV --bin $B --out ${P}${L}.json --label $L "$@" 2>&1 | grep -v "^ *[0-9]\+ \+[0-9.]\+e[+-][0-9]\+"   # drop hypre's per-iteration residual lines
  echo "variant $L rc=${PIPESTATUS[0]} ($(( $(date +%s) - T )) s)"
}
# device, 1 rank
run jacobi                 1 --jacobi --rough-rhs --no-amg
run nf4_pmis_exti_l1jac    1                                   # the brief's GPU settings + num_functions 4 (unknown-based)
run nf1_pmis_exti_l1jac    1 --num-functions 1                 # scalar AMG, block structure ignored
run nf4_exte18             1 --interp 18                       # ext+e (matrix-matrix form)
run nf4_jacobi7            1 --relax 7
run nf4_agg0               1 --agg-levels 0                    # no aggressive coarsening
run nf4_agg0_theta05       1 --agg-levels 0 --theta 0.5
run nf4_agg0_theta07       1 --agg-levels 0 --theta 0.7
run nf4_agg0_ext14         1 --agg-levels 0 --interp 14          # extended interpolation
run nf4_agg0_sweeps2       1 --agg-levels 0 --sweeps 2
run nf4_agg0_pmax8         1 --agg-levels 0 --pmax 8
run nf1_agg0               1 --agg-levels 0 --num-functions 1
run nf4_agg0_theta05_jac7  1 --agg-levels 0 --theta 0.5 --relax 7
run nf4_aggint5            1 --agg-interp 5
run nf4_aggint7            1 --agg-interp 7
run nf4_aggint8            1 --agg-interp 8
run nf4_theta05            1 --theta 0.5
run nf4_sweeps2            1 --sweeps 2
run nf4_ij                 1 --ij                               # HYPRE_IJMatrixSetValues2 from device arrays (must match nf4_pmis_exti_l1jac)
run nf4_spmv_hypre         1 --spmv-vendor 0
run nf4_keepT0             1 --keep-transpose 0
run nf4_2ranks             2                                   # multi-rank on one GPU (offd path), must match the 1-rank count
run nf4_nodal1_DEVICE      1 --nodal 1                         # expected unsupported on the device in hypre 3.1.0 (host-only hypre_BoomerAMGCreateNodalA); recorded either way
# host (reference for the nodal / interpolation-vector variants that only exist on the host)
run host_nf4               1 --host
run host_nf4_nodal1        1 --host --nodal 1
run host_nf4_nodal1_agg0   1 --host --nodal 1 --agg-levels 0
run host_nf4_hmis_agg0     1 --host --coarsen 10 --agg-levels 0   # HMIS (host-only coarsening)
run host_nf4_nodal4        1 --host --nodal 4
run host_nf4_nodal1_gm3    1 --host --nodal 1 --interp-vectors 3 --interp-vec-variant 2      # GM2 with [x,y,z] in the P1-disc basis
run host_nf4_nodal1_gm4    1 --host --nodal 1 --interp-vectors 4 --interp-vec-variant 2      # GM2 with [e,x,y,z]
run host_nf4_nodal1_ln3    1 --host --nodal 1 --interp-vectors 3 --interp-vec-variant 3      # LN
run host_nf4_nodal1_gm3_q4 1 --host --nodal 1 --interp-vectors 3 --interp-vec-variant 2 --interp-vec-qmax 4
python3 hypre_merge_json.py --out $RES/hypre_p1_40x40x120.json --primary nf4_pmis_exti_l1jac ${P}*.json

# ---- 7. half size: dump + first numbers (3 ranks: 2.61 G nnz -> 0.87 G per rank < 2^30) ----
if [ -f $SCR/p1_128x128x384.npz ]; then
  python3 npz_to_bin.py --npz $SCR/p1_128x128x384.npz --out $SCR/bin/p1_128x128x384 2>&1 | tail -3
  B=$SCR/bin/p1_128x128x384; P=$RES/hypre_p1_128x128x384__
  run jacobi              3 --jacobi --no-amg
  run nf4_pmis_exti_l1jac 3
  run nf1_pmis_exti_l1jac 3 --num-functions 1
  python3 hypre_merge_json.py --out $RES/hypre_p1_128x128x384.json --primary nf4_pmis_exti_l1jac ${P}*.json
else echo "half-size npz missing on scratch"; fi

kill $NVPID 2>/dev/null
ls -la $SCR $SCR/bin; df -h /scratch | tail -1
echo "hypre A done $(date -Is)"
