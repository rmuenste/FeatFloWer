#!/bin/bash
# d32_stage_wide_rundir.sh -- stage a D3.2 wide-column (12d x 12d x 24d) DNS rundir
# by cloning the inputs of a finished narrow-column (6d x 6d x 24d) twin.
#
# The narrow rundir is treated as READ-ONLY: only the input files are copied
# (_data/q2p1_param.dat, example.json, cube.json, particles.xyz, start/).
# Results (_dump, _vtk, _gmv, solution, logs) are NOT copied.
#
# The point of the design is that particles.xyz is byte-identical to the twin,
# so wall distance is the only changed variable.
#
# Usage:
#   d32_stage_wide_rundir.sh <src_rundir> <dst_rundir> <job_name> [walltime]
#
# Env overrides:
#   MESHDIR   wide coarse mesh dir       (default <repo>/d32_wide_mesh_v1/coarseB)
#   MESHNAME  partitioned mesh folder    (default DKT_W_108)
#   BINARY    solver binary              (default build-dns-pe-serial/.../q2p1_dkt)
set -euo pipefail

SRC=${1:?usage: $0 <src_rundir> <dst_rundir> <job_name> [walltime]}
DST=${2:?usage: $0 <src_rundir> <dst_rundir> <job_name> [walltime]}
JOBNAME=${3:?usage: $0 <src_rundir> <dst_rundir> <job_name> [walltime]}
WALLTIME=${4:-32:00:00}

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
MESHDIR=${MESHDIR:-$REPO/d32_wide_mesh_v1/coarseB}
MESHNAME=${MESHNAME:-DKT_W_108}
BINARY=${BINARY:-$REPO/build-dns-pe-serial/applications/q2p1_dkt/q2p1_dkt}

[ -d "$SRC" ]                    || { echo "no such src rundir: $SRC"; exit 1; }
[ -d "$MESHDIR/_mesh/$MESHNAME" ] || { echo "no partitioned mesh: $MESHDIR/_mesh/$MESHNAME"; exit 1; }
[ -x "$BINARY" ]                 || { echo "no solver binary: $BINARY"; exit 1; }

mkdir -p "$DST"/{_data,_dump,_vtk,_gmv,_ns,solution,testresults,_mesh}

# --- inputs only, from the narrow twin -------------------------------------
cp -p "$SRC/_data/q2p1_param.dat" "$DST/_data/"
# MG.dat is an INPUT: CreateDumpStructures (source/OutputProfiles.f90:1802-1810) reads a
# refinement-topology lookup table from it. It depends only on NLMAX, not on the mesh, so
# the twin's copy is reused verbatim. Omitting it aborts with "End of file" during init.
cp -p "$SRC/_data/MG.dat" "$DST/_data/"
cp -p "$SRC/example.json" "$SRC/cube.json" "$SRC/particles.xyz" "$DST/"
[ -d "$SRC/start" ] && cp -rp "$SRC/start" "$DST/"

# --- wide mesh: coarse case folder + its 108-way partitioning ---------------
mkdir -p "$DST/coarseB"
for f in "$MESHDIR"/*; do [ -f "$f" ] && cp -p "$f" "$DST/coarseB/"; done
cp -rp "$MESHDIR/_mesh/$MESHNAME" "$DST/_mesh/"

# --- point the deck at the wide mesh; every other parameter untouched -------
sed -i "s|^SimPar@MeshFolder = .*|SimPar@MeshFolder = \"$MESHNAME\"|" "$DST/_data/q2p1_param.dat"
sed -i "s|^SimPar@ProjectFile = .*|SimPar@ProjectFile = \"coarseB/file.prj\"|" "$DST/_data/q2p1_param.dat"

# --- job script -------------------------------------------------------------
cat > "$DST/job.sbatch" <<EOF
#!/bin/bash
#SBATCH --job-name=$JOBNAME
#SBATCH --partition=long
#SBATCH --nodes=2
#SBATCH --ntasks=109
#SBATCH --ntasks-per-node=55
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --prefer=nx
#SBATCH --exclude=worldgames
#SBATCH --time=$WALLTIME
#SBATCH --output=$DST/slurm-%j.out
module load gcc/latest-v13 openmpi/options/interface/ethernet openmpi/4.1.6
cd $DST
FREE_GB=\$(df -BG --output=avail /data/warehouse17 | tail -1 | tr -dc 0-9); [ "\$FREE_GB" -lt "\${FF_MIN_FREE_GB:-100}" ] && { echo "REFUSING to run: only \${FREE_GB}G free on warehouse17"; exit 1; }
mpirun -np 109 $BINARY > run_slurm.log 2>&1
EOF
chmod +x "$DST/job.sbatch"

# --- verify -----------------------------------------------------------------
echo "staged $DST"
echo "  particles.xyz identical to twin : $(cmp -s "$SRC/particles.xyz" "$DST/particles.xyz" && echo YES || echo NO)"
echo "  cube.json     identical to twin : $(cmp -s "$SRC/cube.json"     "$DST/cube.json"     && echo YES || echo NO)"
echo "  example.json  identical to twin : $(cmp -s "$SRC/example.json"  "$DST/example.json"  && echo YES || echo NO)"
echo "  MeshFolder  : $(grep '^SimPar@MeshFolder'  "$DST/_data/q2p1_param.dat")"
echo "  ProjectFile : $(grep '^SimPar@ProjectFile' "$DST/_data/q2p1_param.dat")"
echo "  TimeStep    : $(grep '^SimPar@TimeStep'    "$DST/_data/q2p1_param.dat")"
echo "  MaxNumStep  : $(grep '^SimPar@MaxNumStep'  "$DST/_data/q2p1_param.dat")"
echo "  subdomains  : $(ls -d "$DST/_mesh/$MESHNAME"/sub*/ | wc -l)"
