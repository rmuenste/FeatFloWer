#!/bin/bash
# Stage the column rung set on Fritz: four variants x (u_t + 12 DNS clouds).
#   cont    : ELVolumeAveraged = continuity
#   std     : ELDragModel = difelice_std
#   stdcont : both
#   all     : both + Model-A pressure force (ELFluidGravity=No, ELPressureForce=Yes)
#   contpf  : campaign drag (difelice) + continuity + Model-A pressure force
# VARIANTS="..." overrides the variant list (default: cont std stdcont all).
set -eu
W=/anvme/workspace/k115ce12-featflower/FF/FeatFloWer-el
C=$W/applications/q2p1_el_pipeflow/validation_cases/v2w_column_settling
ROOT=$W/build-el-fritz/applications; mkdir -p $ROOT
STAGE=$W/tools/el_stage_rundir.sh
apply() { # variant deck
  case $1 in
    cont)    printf 'SimPar@ELVolumeAveraged = continuity\n' >> $2 ;;
    std)     sed -i 's/^SimPar@ELDragModel = .*/SimPar@ELDragModel = difelice_std/' $2 ;;
    stdcont) sed -i 's/^SimPar@ELDragModel = .*/SimPar@ELDragModel = difelice_std/' $2; printf 'SimPar@ELVolumeAveraged = continuity\n' >> $2 ;;
    contpf)  sed -i 's/^SimPar@ELFluidGravity = .*/SimPar@ELFluidGravity = No/; s/^SimPar@ELPressureForce = .*/SimPar@ELPressureForce = Yes/' $2; printf 'SimPar@ELVolumeAveraged = continuity\n' >> $2 ;;
    all)     sed -i 's/^SimPar@ELDragModel = .*/SimPar@ELDragModel = difelice_std/; s/^SimPar@ELFluidGravity = .*/SimPar@ELFluidGravity = No/; s/^SimPar@ELPressureForce = .*/SimPar@ELPressureForce = Yes/' $2; printf 'SimPar@ELVolumeAveraged = continuity\n' >> $2 ;;
  esac
}
n=0
for v in ${VARIANTS:-cont std stdcont all}; do
  for r in ut n20_s1 n20_s2 n20_s3 n40_s1 n40_s2 n40_s3 n80_s1 n80_s2 n80_s3 n120_s1 n120_s2 n120_s3; do
    d=$ROOT/q2p1_el_rundir_v2wF_${v}_$r
    EL_REF_DIR=/nonexistent $STAGE -f v2w_column_settling $d > /dev/null
    [ $r = ut ] || cp $C/clouds/$r.xyz $d/particles.xyz
    apply $v $d/_data/q2p1_param.dat
    n=$((n+1))
  done
done
echo "staged $n rundirs under $ROOT"
grep -E 'ELDragModel|ELVolumeAveraged|ELFluidGravity|ELPressureForce' $ROOT/q2p1_el_rundir_v2wF_${v}_n40_s1/_data/q2p1_param.dat
ls $ROOT/q2p1_el_rundir_v2wF_cont_n120_s3/_mesh/DKT_B_108 | wc -l; wc -l < $ROOT/q2p1_el_rundir_v2wF_cont_n120_s3/particles.xyz
du -sh $ROOT
