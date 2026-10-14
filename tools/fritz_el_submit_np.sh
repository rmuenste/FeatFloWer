#!/bin/sh
# fritz_el_submit_np.sh -- submit a staged EL rundir on NHR@FAU Fritz.
# Usage: fritz_el_submit_np.sh <rundir> <jobname> <walltime> <ntasks> <nodes> [partition]
# EL_BIN must point at the q2p1_el_pipeflow binary (build-el-fritz). Fritz rules:
# 24 h limit everywhere, --export=NONE + unset SLURM_EXPORT_ENV, env from ~/ff-env-gcc14.sh.
set -eu
[ $# -ge 5 ] || { echo "usage: $0 <rundir> <jobname> <walltime> <ntasks> <nodes> [partition]" >&2; exit 2; }
RUNDIR=$(CDPATH= cd -- "$1" && pwd); NAME=$2; WALL=$3; NP=$4; NODES=$5; PART=${6:-multinode}
BIN=${EL_BIN:?set EL_BIN}
[ -d "$RUNDIR/_data" ] || { echo "not a staged rundir: $RUNDIR" >&2; exit 1; }
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 1; }
PERNODE=$(( (NP + NODES - 1) / NODES ))
cat > "$RUNDIR/job.sbatch" <<EOT
#!/bin/bash
#SBATCH --job-name=$NAME
#SBATCH --partition=$PART
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NP
#SBATCH --ntasks-per-node=$PERNODE
#SBATCH --cpus-per-task=1
#SBATCH --time=$WALL
#SBATCH --export=NONE
#SBATCH --output=$RUNDIR/slurm-%j.out
unset SLURM_EXPORT_ENV
source ~/ff-env-gcc14.sh
cd $RUNDIR
mpirun -np $NP $BIN > run_slurm.log 2>&1
echo "exit=\$?"
EOT
sbatch "$RUNDIR/job.sbatch"
