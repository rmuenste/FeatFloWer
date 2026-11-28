#!/bin/bash
# Stage a run directory of case 11 (minimum fluidisation) for q2p1_el_pipeflow. Nothing is run
# or submitted (fl11_submit.sh does that).
#
# usage: stage_s11.sh RUNDIR [KEY=VALUE ...]
#   KEY=VALUE   deck keys (e.g. SimPar@MaxNumStep=400 SimPar@TimeStep=0.0001d0): an existing line
#               "KEY = ..." is replaced, a missing SimPar@ key is added after SimPar@MaxNumStep,
#               a missing Pres@ key after Pres@MGCrsSolverType. The json stepsize_ / timesteps_
#               follow SimPar@TimeStep / SimPar@MaxNumStep.
# env:
#   S11_BED=R|A          bed R 15 x 15 mm (default) or bed A 25 x 25 mm, column 100 mm
#   S11_BEDMODE=rsa|lattice   initial bed (make_bed_s11.py; default rsa = the peer's CFDEM bed)
#   S11_SEED=1           bed seed
#   S11_NZ=4             fluid / pe partitions along z (3 x 3 laterally: pe needs >= 3 per periodic axis)
#   S11_INLET=protocol | const:U   inlet table (make_inlet_table.py; default the spec protocol)
# Mesh: coarse cubic cells of about 5 mm, MaxMeshLevel 2 (2.5 mm = 2.9 d for bed R; 2.08 mm =
# 2.4 d for bed A, whose 25 mm do not split into 3 x 3 parts of 5 mm cells).
# Partition: PyPartitioner axis_uniform from the pe_partpy branch fix/axis-uniform-boundary-precision
# (copy in runs/batch10/tools/pe_partpy_axisfix); NPart = 9 S11_NZ = ranks - 1.
set -eu
HERE=$(cd "$(dirname "$0")" && pwd)
CASE=$(cd "$HERE/.." && pwd)
REPO=$(cd "$CASE/../../../.." && pwd)
PART=${S11_PARTITIONER:-/data/warehouse17/rmuenste/runs/batch10/tools/pe_partpy_axisfix/PyPartitioner.py}
R=$1; shift
[ -e "$R" ] && { echo "exists: $R" >&2; exit 1; }
BED=${S11_BED:-R}; MODE=${S11_BEDMODE:-rsa}; SEED=${S11_SEED:-1}; NZ=${S11_NZ:-4}; INLET=${S11_INLET:-protocol}
case $BED in
  R) L=0.015; NXY=3 ;;
  A) L=0.025; NXY=6 ;;
  *) echo "S11_BED must be R or A" >&2; exit 1 ;;
esac
LZ=0.10
NZC=$(python3 -c "print(round($LZ/($L/$NXY)))")
[ $((NZC % NZ)) -eq 0 ] || { echo "coarse z cells $NZC not divisible by S11_NZ=$NZ" >&2; exit 1; }
NP=$((9 * NZ)); MESH=FL${NP}P
mkdir -p "$R"/{_data,_dump,_sol,_vtk,_gmv,_ns,protocols,solution,testresults}
R=$(cd "$R" && pwd)
python3 -I "$HERE/make_s11_mesh.py" "$R/fl_box" --size $L $L $LZ --cells $NXY $NXY $NZC
(cd "$R" && python3 "$PART" $NP axis_uniform x3-y3-z$NZ $MESH fl_box/file.prj > log.partition 2>&1) || \
  { echo "partition failed, see $R/log.partition" >&2; exit 1; }
n=$(/bin/ls -d "$R/_mesh/$MESH"/sub* | wc -l)
[ "$n" -eq $NP ] || { echo "partition: $n sub directories, expected $NP" >&2; exit 1; }
python3 -I "$HERE/make_bed_s11.py" "$R/particles.xyz" --mode $MODE --lx $L --ly $L --seed $SEED
if [ "$INLET" = protocol ]; then
  python3 -I "$HERE/make_inlet_table.py" "$R/_data/inlet_table.txt" --steps "$R/_data/steps.txt"
else
  python3 -I "$HERE/make_inlet_table.py" "$R/_data/inlet_table.txt" --const "${INLET#const:}"
fi
cp "$REPO/_data/MG.dat" "$R/_data/MG.dat"
D=$R/_data/q2p1_param.dat
sed -e "s|^SimPar@MeshFolder = .*|SimPar@MeshFolder = \"$MESH\"|" \
    -e "s|^SimPar@SubMeshNumber = .*|SimPar@SubMeshNumber = $NP|" \
    -e "s|^SimPar@PeriodicLength = .*|SimPar@PeriodicLength = ${L}d0|" "$CASE/q2p1_param_s11.dat" > "$D"
for kv in "$@"; do
  k=${kv%%=*}; v=${kv#*=}
  if grep -q "^$k *=" "$D"; then sed -i "s|^$k *=.*|$k = $v|" "$D"
  else case $k in
    Pres@*) sed -i "/^Pres@MGCrsSolverType/a $k = $v" "$D" ;;
    *) sed -i "/^SimPar@MaxNumStep/a $k = $v" "$D" ;;
  esac; fi
done
DT=$(grep '^SimPar@TimeStep' "$D" | sed 's/.*= *//; s/[dD]/e/')
NS=$(grep '^SimPar@MaxNumStep' "$D" | sed 's/.*= *//')
python3 - "$CASE/example_s11.json" "$R/example.json" "$DT" "$NS" "$NZ" <<'EOF'
import json, sys
j = json.load(open(sys.argv[1]))
j['stepsize_'] = float(sys.argv[3]); j['timesteps_'] = int(sys.argv[4]); j['processesZ_'] = int(sys.argv[5])
json.dump(j, open(sys.argv[2], 'w'), indent=1)
EOF
{
  echo "staged $(date '+%F %T') by stage_s11.sh from $(git -C "$REPO" rev-parse --short HEAD) ($REPO)"
  echo "bed $BED L $L m, mode $MODE seed $SEED, N $(wc -l < "$R/particles.xyz"), md5 $(md5sum "$R/particles.xyz" | cut -c1-12)"
  echo "mesh coarse $NXY x $NXY x $NZC, MaxMeshLevel $(grep '^SimPar@MaxMeshLevel' "$D" | sed 's/.*= *//'), partition $MESH x3-y3-z$NZ ($NP + 1 ranks)"
  echo "inlet $INLET; dt $DT; steps $NS"
  echo "deck changes against q2p1_param_s11.dat:"; diff "$CASE/q2p1_param_s11.dat" "$D" | grep '^[<>]' || true
} > "$R/staged.txt"
cat "$R/staged.txt"
echo "ranks $((NP + 1))"
