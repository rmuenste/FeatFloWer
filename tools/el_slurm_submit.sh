#!/bin/sh
# Submit an already-staged EL validation run dir (see el_stage_rundir.sh) as a
# single-node np=28 Slurm job. Usage:
#   tools/el_slurm_submit.sh <rundir> <jobname> <walltime e.g. 1-12:00:00> [partition]
# The binary is taken from EL_BIN or the default build tree. Output log:
# <rundir>/run_slurm.log, Slurm log <rundir>/slurm-%j.out.
set -eu
RUNDIR=$1; NAME=$2; WALL=$3; PART=${4:-long}
BIN=${EL_BIN:-/home/user/rmuenste/nobackup/code/FF-EL/FeatFloWer/build-el-phase2-pe-gcc14/applications/q2p1_el_pipeflow/q2p1_el_pipeflow}
[ -d "$RUNDIR/_data" ] || { echo "not a staged rundir: $RUNDIR" >&2; exit 1; }
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 1; }

cat > "$RUNDIR/job.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=$NAME
#SBATCH --partition=$PART
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --cpus-per-task=1
#SBATCH --time=$WALL
#SBATCH --mem=25G
#SBATCH --output=$RUNDIR/slurm-%j.out
source /etc/profile >/dev/null 2>&1 || true
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
BIN=$BIN
MPIPREFIX=\$(ldd "\$BIN" | grep -m1 'libmpi\\.so' | awk '{print \$3}' | sed 's#/lib/libmpi.*##')
export PATH="\$MPIPREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="\$MPIPREFIX/lib:\$LD_LIBRARY_PATH"
cd $RUNDIR
mpirun -np 28 "\$BIN" > run_slurm.log 2>&1
EOF
sbatch "$RUNDIR/job.sbatch"
