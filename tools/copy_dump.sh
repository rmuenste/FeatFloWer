#!/usr/bin/env bash

nProc=15
inPath='/data/warehouse14/rjendrny/IANUS/IANUS_Feat/bin_arrow/applications/q2p1_devel/'
outPath='.'
inDump=6
outDump=11
fields=('velocity' 'pressure' 'coordinates')

for ((j=1; j<=${nProc};j++))
do
 cOut=$outPath'/_dump/processor_'$j'/'$outDump
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