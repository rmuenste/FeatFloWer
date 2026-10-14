#!/bin/sh
# Submit an already-staged EL validation run dir with an arbitrary rank/node
# shape (the sibling el_slurm_submit.sh is fixed at np=28 on one node).
# Usage:
#   tools/el_slurm_submit_np.sh <rundir> <jobname> <walltime> <partition> <ntasks> <nodes> [mem]
# The binary comes from EL_BIN (required). Output: <rundir>/run_slurm.log and
# <rundir>/slurm-%j.out. Modules match the gcc-14.3 OpenMPI build.
set -eu
[ $# -ge 6 ] || { echo "usage: $0 <rundir> <jobname> <walltime> <partition> <ntasks> <nodes> [mem]" >&2; exit 2; }
RUNDIR=$(CDPATH= cd -- "$1" && pwd); NAME=$2; WALL=$3; PART=$4; NP=$5; NODES=$6; MEM=${7:-25G}
BIN=${EL_BIN:?set EL_BIN to the q2p1_el_pipeflow binary}
[ -d "$RUNDIR/_data" ] || { echo "not a staged rundir: $RUNDIR" >&2; exit 1; }
[ -x "$BIN" ] || { echo "binary not found: $BIN" >&2; exit 1; }
PERNODE=$(( (NP + NODES - 1) / NODES ))
cat > "$RUNDIR/job.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=$NAME
#SBATCH --partition=$PART
#SBATCH --nodes=$NODES
#SBATCH --ntasks=$NP
#SBATCH --ntasks-per-node=$PERNODE
#SBATCH --cpus-per-task=1
#SBATCH --time=$WALL
#SBATCH --mem=$MEM
#SBATCH --prefer=nx
#SBATCH --exclude=worldgames
#SBATCH --output=$RUNDIR/slurm-%j.out
source /etc/profile >/dev/null 2>&1 || true
module load gcc/14.3.0 openmpi/4.1.6
BIN=$BIN
MPIPREFIX=\$(ldd "\$BIN" | grep -m1 'libmpi\\.so' | awk '{print \$3}' | sed 's#/lib/libmpi.*##')
export PATH="\$MPIPREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="\$MPIPREFIX/lib:\$LD_LIBRARY_PATH"
cd $RUNDIR
mpirun -np $NP "\$BIN" > run_slurm.log 2>&1
echo "exit=\$?"
EOF
sbatch "$RUNDIR/job.sbatch"
