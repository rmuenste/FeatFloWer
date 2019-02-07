#!/usr/bin/env bash

nProc=15
inPath='/home/user/omierka/nobackup/Fortran/NewGen_Q2P1/Feat_FloWer_v3.0/bin-intel-epyx/applications/q2p1_devel/'
outPath='.'
inDump=4
outDump=11
fields=('velocity' 'pressure' 'coordinates')

for ((j=1; j<=${nProc};j++))
do
 cOut=$outPath'/_dump/processor_'$j'/'$outDump
 mkdir $outPath'/_dump/processor_'$j
 cIn=$inPath'_dump/processor_'$j'/'$inDump
#  echo $cIn, $cOut, ${#fields[@]}
 mkdir $cOut
 for ((k=0; k<${#fields[@]};k++))
 do
   cmd='cp '$cIn'/'${fields[$k]}'.dmp '$cOut'/.'
   echo $k $cmd
   `$cmd`
 done
done

exit 0
