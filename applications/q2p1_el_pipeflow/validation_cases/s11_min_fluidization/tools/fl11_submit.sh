#!/bin/bash
# Write job.sbatch + versions.txt in a staged case-11 run directory and submit it (pattern of
# runs/pcoarse/tools/pc_submit.sh). Jobs point to a COPY of the binary (never to a build tree).
# usage: fl11_submit.sh RUNDIR BINARY_COPY [time hh:mm:ss] [partition]   (FL11_NOSUBMIT=1: write only)
#   ranks = SubMeshNumber + 1 of the deck; BINARY_COPY.commit (text) is recorded if present.
set -eu
R=$(cd "$1" && pwd); BIN=$2; TM=${3:-00:30:00}; PART=${4:-short}
W=$(cd "$(dirname "$0")/../../../../.." && pwd)
[ -x "$BIN" ] || { echo "no binary $BIN" >&2; exit 1; }
# one run per directory: a second run appends its particle frames to the rank logs of the first. The lock
# directory (mkdir is atomic) refuses a second submission even before the first job has written a log; the
# generated job checks again before it truncates run_slurm.log (a requeue or a job started by hand).
if compgen -G "$R/el_particles_rank*.log" > /dev/null; then
  echo "refused: $R already holds particle rank logs (el_particles_rank*.log) of an earlier run; a new run would" \
       "append to them and fl11_movie.py would refuse the mixed frames. Stage a fresh run directory, or move the" \
       "old rank logs (and run_slurm.log) out of $R, then submit again." >&2
  exit 6
fi
if [ "${FL11_NOSUBMIT:-0}" != 1 ] && ! mkdir "$R/.fl11_submit.lock" 2>/dev/null; then
  echo "refused: $R/.fl11_submit.lock exists (a job was already submitted into this directory: $(cat "$R/.fl11_submit.lock/job" 2>/dev/null || echo unknown)). Stage a fresh run directory." >&2
  exit 7
fi
NP=$(( $(grep '^SimPar@SubMeshNumber' "$R/_data/q2p1_param.dat" | sed 's/.*= *//') + 1 ))
{
  echo "FeatFloWer worktree $W"
  echo "binary source $(cat "$BIN.commit" 2>/dev/null || echo unknown)"
  echo "worktree HEAD at submission $(git -C "$W" rev-parse HEAD) ($(git -C "$W" rev-parse --abbrev-ref HEAD)), dirty $(git -C "$W" -c core.autocrlf=false status --short | grep -vc '^??' || true)"
  echo "libs/pe of the worktree at submission: $(cat "$W/libs/pe/.pin" 2>/dev/null || echo '43e8fb541edc6d0c0d533a49ce4246786b36ff18 (rsync of the FeatFloWer-perkernel submodule tree)'); the pe of the BINARY is in its source line above"
  echo "binary $BIN md5 $(md5sum "$BIN" | cut -c1-32) built on $(cat "$BIN.host" 2>/dev/null || echo ?)"
  echo "bed md5 $(md5sum "$R/particles.xyz" | cut -c1-12), $(wc -l < "$R/particles.xyz") particles"
  echo "deck keys: $(grep -hE '^SimPar@(TimeStep|MaxNumStep|ELDragModel|ELMomentumModel|ELMMPressureOperator|ELInletVelocityTable|ELBedAuditFreq|ELDragImplicitWeights|ELVelocityAuditFreq|ELParticleLogFreq|PeriodicAxis|PeriodicLength|NoOutflow)|^Pres@MGCrsSolverType' "$R/_data/q2p1_param.dat" | tr '\n' ';')"
  echo "inlet table md5 $(md5sum "$R/_data/inlet_table.txt" | cut -c1-12)"
  echo "ranks $NP, partition $PART, time $TM; written $(date '+%F %T') by fl11_submit.sh"
} > "$R/versions.txt"
cat > "$R/job.sbatch" <<EOT
#!/bin/bash
#SBATCH --job-name=$(basename "$R")
#SBATCH --partition=$PART
#SBATCH --nodes=1
#SBATCH --ntasks=$NP
#SBATCH --ntasks-per-node=$NP
#SBATCH --cpus-per-task=1
#SBATCH --time=$TM
#SBATCH --mem=40G
#SBATCH --prefer=nx
#SBATCH --exclude=worldgames,82ndalltheway,christmastruce,gottmituns,redbaron,soldierofheaven,steelcommanders,stormtroopers,wolfpack
#SBATCH --output=$R/slurm-%j.out
source /etc/profile >/dev/null 2>&1 || true
module load gcc/14.3.0 openmpi/4.1.6
BIN=$BIN
MPIPREFIX=\$(ldd "\$BIN" | grep -m1 'libmpi\.so' | awk '{print \$3}' | sed 's#/lib/libmpi.*##')
export PATH="\$MPIPREFIX/bin:\$PATH"
export LD_LIBRARY_PATH="\$MPIPREFIX/lib:\$LD_LIBRARY_PATH"
cd $R
if compgen -G "el_particles_rank*.log" > /dev/null; then
  echo "refused: $R holds particle rank logs of an earlier run (requeue or second start): not run, run_slurm.log kept"
  exit 6
fi
echo "host \$(hostname) job \$SLURM_JOB_ID start \$(date '+%F %T')"
mpirun -np $NP "\$BIN" > run_slurm.log 2>&1
status=\$?
echo "exit=\$status end \$(date '+%F %T')"
exit "\$status"
EOT
[ "${FL11_NOSUBMIT:-0}" = 1 ] && { echo "staged, not submitted: $R/job.sbatch"; exit 0; }
cd "$R" && out=$(sbatch --no-requeue job.sbatch) && echo "$out" | tee "$R/.fl11_submit.lock/job"
